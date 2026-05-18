#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import logging
import time
from pathlib import Path

import requests
from tqdm import tqdm

GDC_API = "https://api.gdc.cancer.gov"


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1_048_576), b""):
            h.update(chunk)
    return h.hexdigest()


def _query_files(project: str, data_type: str, workflow: str, max_files: int) -> list[dict]:
    payload = {
        "filters": {
            "op": "and",
            "content": [
                {"op": "=", "content": {"field": "cases.project.project_id",  "value": project}},
                {"op": "=", "content": {"field": "data_type",                 "value": data_type}},
                {"op": "=", "content": {"field": "analysis.workflow_type",    "value": workflow}},
                {"op": "=", "content": {"field": "cases.samples.sample_type", "value": "Primary Tumor"}},
            ],
        },
        "fields": "file_id,file_name,md5sum,file_size,"
                  "cases.case_id,cases.submitter_id",
        "format": "json",
        "size":   max_files,
    }

    response = requests.post(f"{GDC_API}/files", json=payload, timeout=60)
    response.raise_for_status()

    hits = response.json().get("data", {}).get("hits", [])
    logging.info("[TCGA] %d expression file(s) found for %s", len(hits), project)
    return hits


def _query_clinical(project: str) -> list[dict]:
    payload = {
        "filters": {
            "op": "=",
            "content": {"field": "project.project_id", "value": project},
        },
        "fields": (
            "case_id,"
            "submitter_id,"
            "diagnoses.vital_status,"
            "diagnoses.days_to_death,"
            "diagnoses.days_to_last_follow_up,"
            "diagnoses.age_at_diagnosis,"
            "diagnoses.primary_diagnosis,"
            "demographic.gender,"
            "demographic.race,"
            "demographic.ethnicity"
        ),
        "format": "json",
        "size":   10_000,
    }

    response = requests.post(f"{GDC_API}/cases", json=payload, timeout=60)
    response.raise_for_status()

    hits = response.json().get("data", {}).get("hits", [])
    logging.info("[TCGA] %d clinical record(s) found for %s", len(hits), project)
    return hits


def _flatten_clinical(hits: list[dict]) -> list[dict]:
    rows = []
    for case in hits:
        diagnoses = case.get("diagnoses") or [{}]
        diag = diagnoses[0]
        demo = case.get("demographic") or {}

        rows.append({
            "case_id":                case.get("case_id", ""),
            "submitter_id":           case.get("submitter_id", ""),
            "vital_status":           diag.get("vital_status", ""),
            "days_to_death":          diag.get("days_to_death"),
            "days_to_last_follow_up": diag.get("days_to_last_follow_up"),
            "age_at_diagnosis":       diag.get("age_at_diagnosis"),
            "primary_diagnosis":      diag.get("primary_diagnosis", ""),
            "gender":                 demo.get("gender", ""),
            "race":                   demo.get("race", ""),
            "ethnicity":              demo.get("ethnicity", ""),
        })
    return rows


def _write_chunks(response: requests.Response, dest: Path) -> None:
    total = int(response.headers.get("content-length", 0)) or None
    
    with open(dest, "wb") as f, tqdm(
        desc=dest.name, total=total,
        unit="B", unit_scale=True, unit_divisor=1024,
        ascii=False, miniters=1,
    ) as pbar:
        for chunk in response.iter_content(chunk_size=1_048_576):
            if chunk:
                f.write(chunk)
                pbar.update(len(chunk))


def _download_file(file_id: str, dest: Path, expected_md5: str, retries: int = 3) -> bool:
    if dest.exists() and expected_md5 and _md5(dest) == expected_md5:
        logging.info("[TCGA] Already verified — skipping: %s", dest.name)
        return True

    url = f"{GDC_API}/data/{file_id}"

    for attempt in range(1, retries + 1):
        try:
            response = requests.get(url, stream=True, timeout=(15, 120))
            response.raise_for_status()

            _write_chunks(response, dest)

            if expected_md5 and _md5(dest) != expected_md5:
                raise ValueError(f"MD5 mismatch for {dest.name}")

            return True

        except Exception:
            if attempt == retries:
                logging.exception("[TCGA] Failed: %s", dest.name)
                return False
            
            logging.warning("[TCGA] Retry %d/%d", attempt, retries, exc_info=True)
            time.sleep(10 * attempt)

    return False


def download_tcga(
    project: str   = "TCGA-BRCA",
    data_type: str = "Gene Expression Quantification",
    workflow: str  = "STAR - Counts",
    out_dir: str   = "data/tcga",
    max_files: int = 1_000,
    force: bool    = False,
) -> dict[str, Path]:
    out_path = Path(out_dir).resolve()
    expr_dir = out_path / "expression"
    expr_dir.mkdir(parents=True, exist_ok=True)

    logging.info("[TCGA] Querying %s — %s / %s", project, data_type, workflow)
    hits = _query_files(project, data_type, workflow, max_files)

    if not hits:
        raise RuntimeError(f"No files found for {project} / {data_type} / {workflow}")

    manifest: dict[str, str] = {}
    ok, failed = 0, 0

    for hit in hits:
        file_id   = hit["file_id"]
        file_name = hit["file_name"]
        md5sum    = hit.get("md5sum", "")
        dest      = expr_dir / file_name

        if not force and dest.exists() and md5sum and _md5(dest) == md5sum:
            logging.info("[TCGA] Already verified — skipping: %s", file_name)
            manifest[file_id] = str(dest)
            ok += 1
            continue

        success = _download_file(file_id, dest, md5sum)
        if success:
            manifest[file_id] = str(dest)
            ok += 1
        else:
            failed += 1

    logging.info("[TCGA] Expression download: %d ok | %d failed", ok, failed)

    manifest_path = out_path / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    logging.info("[TCGA] Manifest → %s", manifest_path)

    logging.info("[TCGA] Downloading clinical data for %s", project)
    clinical_hits = _query_clinical(project)

    if clinical_hits:
        flat = _flatten_clinical(clinical_hits)

        clinical_path = out_path / "clinical.json"
        with open(clinical_path, "w") as f:
            json.dump(flat, f, indent=2)
        logging.info("[TCGA] Clinical data (%d records) → %s", len(flat), clinical_path)

    return {fid: Path(p) for fid, p in manifest.items()}


if __name__ == "__main__":
    import argparse

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)-7s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    p = argparse.ArgumentParser()
    p.add_argument("--project",   default="TCGA-BRCA")
    p.add_argument("--data-type", default="Gene Expression Quantification")
    p.add_argument("--workflow",  default="STAR - Counts")
    p.add_argument("--out-dir",   default="data/tcga")
    p.add_argument("--max-files", type=int, default=1_000)
    p.add_argument("--force",     action="store_true")
    args = p.parse_args()

    download_tcga(
        project   = args.project,
        data_type = args.data_type,
        workflow  = args.workflow,
        out_dir   = args.out_dir,
        max_files = args.max_files,
        force     = args.force,
    )