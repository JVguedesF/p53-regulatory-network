from __future__ import annotations

import csv
import hashlib
import logging
import sys
import time
from pathlib import Path
from typing import Sequence

import requests
import yaml
from tqdm import tqdm

from src.resolver import infer_sample_type, resolve_sample_type

TSV_HEADER = [
    "Condition",
    "Accession",
    "Sample_Type",
    "Pair_ID",
    "Paired_End",
    "Replicate",
    "Input_Accession",
    "FASTQ_Path",
    "Download_Links",
    "File_Sizes",
    "MD5_Checksums",
    "QC_Status",
    "Trimmed_Path",
    "BAM_Path",
    "Alignment_Rate",
    "Peak_File",
    "N_Peaks",
    "Annotated",
    "DEG_Status",
]

def load_config(config_path: Path) -> dict:
    try:
        with open(config_path, encoding="utf-8") as f:
            return yaml.safe_load(f) or {}
    except FileNotFoundError:
        return {}

def setup_environment(project_dir: Path, config: dict) -> tuple[Path, Path]:
    paths = config.get("paths", {})
    log_dir = project_dir / paths.get("log_dir", "logs")
    tsv_dir = project_dir / paths.get("tsv_dir", "data/tsv")
    raw_dir = project_dir / paths.get("raw_dir", "data/raw")

    for d in (log_dir, tsv_dir, raw_dir):
        d.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)-7s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.FileHandler(log_dir / "downloader.log"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    return tsv_dir, raw_dir

def write_tsv(rows: list[dict], filepath: Path) -> None:
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=TSV_HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

class BaseDownloader:
    def __init__(self, raw_dir: Path, sample_type_overrides: dict[str, str] | None = None) -> None:
        self.raw_dir = raw_dir
        self.overrides: dict[str, str] = sample_type_overrides or {}

    def fetch(self, experiment_id: str) -> tuple[list[dict], str | None]:
        raise NotImplementedError

    def request_data(self, url: str, headers: dict | None = None) -> tuple:
        try:
            response = requests.get(url, headers=headers, timeout=(15, 60))
            response.raise_for_status()
            return response.json(), None
        except requests.exceptions.RequestException as exc:
            return None, str(exc)
        except ValueError as exc:
            return None, f"JSON decode error: {exc}"

    def calculate_md5(self, file_path: Path) -> str:
        h = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(1_048_576), b""):
                h.update(chunk)
        return h.hexdigest()

    def _file_is_valid(self, dest_path: Path, expected_size: int, expected_md5: str) -> bool:
        if not dest_path.exists():
            return False
        if expected_size > 0 and dest_path.stat().st_size != expected_size:
            dest_path.unlink()
            return False
        if expected_md5 and self.calculate_md5(dest_path) != expected_md5:
            dest_path.unlink()
            return False
        return True

    def _download(self, url: str, dest_path: Path, expected_size: int, expected_md5: str) -> None:
        initial_pos = dest_path.stat().st_size if dest_path.exists() else 0
        headers = {"Range": f"bytes={initial_pos}-"} if initial_pos > 0 else {}

        response = requests.get(url, headers=headers, stream=True, timeout=30)
        response.raise_for_status()

        if initial_pos > 0 and response.status_code == 200:
            initial_pos = 0
            open(dest_path, "wb").close()

        content_len = int(response.headers.get("content-length", 0))
        total = content_len + initial_pos if content_len else None

        mode = "ab" if initial_pos > 0 else "wb"
        with open(dest_path, mode) as f, tqdm(
            desc=dest_path.name,
            total=total,
            unit="B", unit_scale=True, unit_divisor=1024,
            initial=initial_pos, ascii=False, miniters=1,
        ) as pbar:
            for chunk in response.iter_content(chunk_size=1_048_576):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))

        if expected_size > 0 and dest_path.stat().st_size != expected_size:
            raise ValueError(f"Size mismatch for {dest_path.name}")
        if expected_md5 and self.calculate_md5(dest_path) != expected_md5:
            raise ValueError(f"MD5 mismatch for {dest_path.name}")

    def download_file(
        self,
        url: str,
        dest_path_str: str,
        expected_size: str,
        expected_md5: str,
        retries: int = 3,
    ) -> bool:
        dest_path = Path(dest_path_str)
        exp_size = int(expected_size) if str(expected_size).isdigit() else 0

        if self._file_is_valid(dest_path, exp_size, expected_md5):
            sys.stdout.write(f"{dest_path.name}: já verificado, pulando.\n")
            return True

        for attempt in range(retries):
            try:
                self._download(url, dest_path, exp_size, expected_md5)
                return True
            except Exception as exc:
                if attempt == retries - 1:
                    logging.error(f"Falha ao baixar {dest_path.name}: {exc}")
                    return False
                time.sleep(10)
        return False

    def _empty_row(self) -> dict[str, str]:
        return dict.fromkeys(TSV_HEADER, "")

class EncodeDownloader(BaseDownloader):
    def _resolve_sample_type(self, file_data: dict, condition: str) -> str:
        acc = file_data.get("accession", "")
        override = self.overrides.get(acc, "")
        if file_data.get("control_type"):
            return override or "Input"
        target = file_data.get("target") or {}
        target_label = target.get("label", "") if isinstance(target, dict) else ""
        return infer_sample_type(condition + " " + target_label, override)

    def _extract_fastq_rows(self, data: dict, condition: str) -> list[dict]:
        base_url = "https://www.encodeproject.org"
        rows: list[dict] = []

        for file in data.get("files", []):
            is_fastq = file.get("file_type") == "fastq" or file.get("file_format") == "fastq"
            if not (is_fastq and file.get("status") == "released"):
                continue

            rep = file.get("replicate") or {}
            pair_val = str(file.get("paired_end", ""))
            acc = file.get("accession", "N/A")
            is_pe = pair_val in ("1", "2")
            acc_key = f"{acc}_{pair_val}" if is_pe else acc

            row = self._empty_row()
            row.update({
                "Condition": condition,
                "Accession": acc_key,
                "Sample_Type": self._resolve_sample_type(file, condition),
                "Pair_ID": pair_val if is_pe else "",
                "Paired_End": str(is_pe),
                "Replicate": str(rep.get("biological_replicate_number", "N/A")),
                "FASTQ_Path": str(self.raw_dir / f"{acc_key}.fastq.gz"),
                "Download_Links": base_url + file.get("href", ""),
                "File_Sizes": str(file.get("file_size", "0")),
                "MD5_Checksums": file.get("md5sum", ""),
            })
            rows.append(row)
        return rows

    def fetch(self, experiment_id: str) -> tuple[list[dict], str | None]:
        url = f"https://www.encodeproject.org/experiments/{experiment_id}/?format=json"
        data, err = self.request_data(url, headers={"accept": "application/json"})
        if err or not isinstance(data, dict):
            return [], err or "Resposta inválida"

        target = data.get("target") or {}
        label = target.get("label", "") if isinstance(target, dict) else ""
        condition = f"{data.get('assay_term_name', 'Unknown')}_{label}".strip("_")

        rows = self._extract_fastq_rows(data, condition)
        return (rows, None) if rows else ([], "Nenhum FASTQ liberado encontrado")


class EnaDownloader(BaseDownloader):
    def _process_entry(self, entry: dict, condition: str) -> list[dict]:
        ftp_raw = entry.get("fastq_ftp", "") or entry.get("sra_ftp", "")
        if not ftp_raw:
            return []

        ftp_links = [f"http://{l}" if not l.startswith("http") else l for l in ftp_raw.split(";") if l]
        sizes = (entry.get("fastq_bytes", "") or "").split(";")
        md5s = (entry.get("fastq_md5", "") or "").split(";")

        run_acc = entry.get("run_accession", "N/A")
        replicate = entry.get("sample_accession", "N/A")
        sample_type = resolve_sample_type(run_acc, self.overrides.get(run_acc, ""))

        is_pe = entry.get("library_layout") == "PAIRED" and len(ftp_links) == 2
        if is_pe:
            return self._pe_rows(condition, sample_type, run_acc, replicate, ftp_links, sizes, md5s)
        return [self._se_row(condition, sample_type, run_acc, replicate, ftp_links, sizes, md5s)]

    def _pe_rows(
        self, condition: str, sample_type: str, run_acc: str, replicate: str,
        links: Sequence[str], sizes: Sequence[str], md5s: Sequence[str],
    ) -> list[dict]:
        rows = []
        for i, (link, pair_id) in enumerate(zip(links, ("1", "2"))):
            acc_key = f"{run_acc}_{pair_id}"
            row = self._empty_row()
            row.update({
                "Condition": condition,
                "Accession": acc_key,
                "Sample_Type": sample_type,
                "Pair_ID": pair_id,
                "Paired_End": "True",
                "Replicate": replicate,
                "FASTQ_Path": str(self.raw_dir / f"{acc_key}.fastq.gz"),
                "Download_Links": link,
                "File_Sizes": sizes[i] if i < len(sizes) else "",
                "MD5_Checksums": md5s[i] if i < len(md5s) else "",
            })
            rows.append(row)
        return rows

    def _se_row(
        self, condition: str, sample_type: str, run_acc: str, replicate: str,
        links: Sequence[str], sizes: Sequence[str], md5s: Sequence[str],
    ) -> dict:
        row = self._empty_row()
        row.update({
            "Condition": condition,
            "Accession": run_acc,
            "Sample_Type": sample_type,
            "Pair_ID": "",
            "Paired_End": "False",
            "Replicate": replicate,
            "FASTQ_Path": str(self.raw_dir / f"{run_acc}.fastq.gz"),
            "Download_Links": links[0] if links else "",
            "File_Sizes": sizes[0] if sizes else "",
            "MD5_Checksums": md5s[0] if md5s else "",
        })
        return row

    def fetch(self, experiment_id: str) -> tuple[list[dict], str | None]:
        url = (
            f"https://www.ebi.ac.uk/ena/portal/api/filereport"
            f"?accession={experiment_id}&result=read_run"
            f"&fields=run_accession,sample_accession,library_layout,"
            f"fastq_ftp,fastq_bytes,fastq_md5,sra_ftp,sra_bytes,sra_md5"
            f"&format=json"
        )
        data, err = self.request_data(url)
        if err:
            return [], err
        if not data:
            return [], "Resposta vazia"

        rows: list[dict] = []
        for entry in data:
            rows.extend(self._process_entry(entry, experiment_id))
        return (rows, None) if rows else ([], "Nenhuma entrada válida encontrada")

DOWNLOADER_REGISTRY: dict[str, type[BaseDownloader]] = {
    "ENC": EncodeDownloader,
    "GSE": EnaDownloader,
    "PRJ": EnaDownloader,
    "SRR": EnaDownloader,
    "ERR": EnaDownloader,
}

def main(target_ids: list[str] | None = None) -> None:
    project_dir = Path(__file__).parent.parent
    config_path = project_dir / "config" / "config.yaml"
    
    config = load_config(config_path)
    tsv_dir, raw_dir = setup_environment(project_dir, config)

    experiments: dict = config.get("downloads", {})
    ids = target_ids or list(experiments.keys())

    logging.info("Iniciando pipeline de download")

    for eid in ids:
        downloader_class = DOWNLOADER_REGISTRY.get(eid[:3].upper())
        if not downloader_class:
            logging.warning(f"[PULAR] {eid}: prefixo desconhecido.")
            continue

        exp_info = experiments.get(eid, {})
        overrides = exp_info.get("overrides", {})

        try:
            downloader = downloader_class(raw_dir=raw_dir, sample_type_overrides=overrides)
            rows, error = downloader.fetch(eid)

            if error:
                logging.error(f"[{eid}] Falha na busca: {error}")
                continue
            if not rows:
                logging.warning(f"[{eid}] Nenhum dado válido encontrado.")
                continue

            tsv_path = tsv_dir / f"{eid}_metadata.tsv"
            write_tsv(rows, tsv_path)
            logging.info(f"[{eid}] TSV gerado: {tsv_path.name} ({len(rows)} linhas)")

            for row in rows:
                downloader.download_file(
                    row["Download_Links"],
                    row["FASTQ_Path"],
                    row["File_Sizes"],
                    row["MD5_Checksums"],
                )

            time.sleep(5)

        except Exception as exc:
            logging.critical(f"[{eid}] Falha crítica: {exc}")

if __name__ == "__main__":
    main(sys.argv[1:] or None)