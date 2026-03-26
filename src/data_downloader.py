from __future__ import annotations

import csv
import hashlib
import logging
import os
import sys
import time

import requests
import yaml
from pathlib import Path
from tqdm import tqdm

from src.resolver import INPUT_KEYWORDS, infer_sample_type, EntrezResolver

# ── Config ────────────────────────────────────────────────────────────────────
_CONFIG_PATH = Path(__file__).parent.parent / "config" / "config.yaml"

def _load_config() -> dict:
    try:
        with open(_CONFIG_PATH, encoding="utf-8") as f:
            return yaml.safe_load(f) or {}
    except FileNotFoundError:
        logging.warning(f"config.yaml not found at {_CONFIG_PATH}. Using defaults.")
        return {}

_config  = _load_config()
_paths   = _config.get("paths", {})

PROJECT_DIR = Path(__file__).parent.parent
LOG_DIR     = PROJECT_DIR / _paths.get("log_dir",  "logs")
TSV_DIR     = PROJECT_DIR / _paths.get("tsv_dir",  "data/tsv")
RAW_DIR     = PROJECT_DIR / _paths.get("raw_dir",  "data/raw")

for _dir in [LOG_DIR, TSV_DIR, RAW_DIR]:
    _dir.mkdir(parents=True, exist_ok=True)

# ── Logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler(LOG_DIR / "downloader.log"),
        logging.StreamHandler(sys.stdout),
    ],
)

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
    "QC_Status",         # scripts/01_qc_trimming.sh
    "Trimmed_Path",      # scripts/01_qc_trimming.sh
    "BAM_Path",          # scripts/02_alignment.sh
    "Alignment_Rate",    # scripts/02_alignment.sh
    "Peak_File",         # scripts/03_peak_calling.sh (ChIP-seq only)
    "N_Peaks",           # scripts/03_peak_calling.sh (ChIP-seq only)
    "Annotated",         # src/chipseeker.R           (ChIP-seq only)
    "DEG_Status",        # src/deseq2.R               (RNA-seq only)
]


def write_tsv(rows: list[dict], filepath: Path | str) -> None:
    """Write a list of row dicts to a TSV file using TSV_HEADER column order."""
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=TSV_HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


# ── Base downloader ───────────────────────────────────────────────────────────
class BaseDownloader:
    """Shared network and I/O utilities for all downloaders."""

    def __init__(self, sample_type_overrides: dict[str, str] | None = None) -> None:
        self.overrides: dict[str, str] = sample_type_overrides or {}

    def fetch(self, experiment_id: str) -> tuple[list[dict], str | None]:
        """Fetch metadata for a given experiment ID."""
        raise NotImplementedError("Subclasses must implement the 'fetch' method.")

    def request_data(self, url: str, headers: dict | None = None) -> tuple:
        try:
            response = requests.get(url, headers=headers, timeout=(15, 60))
            response.raise_for_status()
            return response.json(), None
        except requests.exceptions.RequestException as e:
            return None, str(e)
        except ValueError as e:
            return None, f"JSON Decode Error: {e}"

    def calculate_md5(self, file_path: str) -> str:
        hash_md5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(1048576), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()

    def _check_existing_file(self, dest_path: str, expected_size: int, expected_md5: str) -> bool:
        if not os.path.exists(dest_path):
            return False
        if expected_size > 0 and os.path.getsize(dest_path) != expected_size:
            os.remove(dest_path)
            return False
        if expected_md5:
            if self.calculate_md5(dest_path) == expected_md5:
                return True
            os.remove(dest_path)
            return False
        return True

    def _perform_download_and_verify(
        self, url: str, dest_path: str, expected_size: int, expected_md5: str
    ) -> bool:
        file_name   = os.path.basename(dest_path)
        initial_pos = os.path.getsize(dest_path) if os.path.exists(dest_path) else 0
        headers     = {"Range": f"bytes={initial_pos}-"} if initial_pos > 0 else {}

        response = requests.get(url, headers=headers, stream=True, timeout=30)
        response.raise_for_status()

        if initial_pos > 0 and response.status_code == 200:
            initial_pos = 0
            open(dest_path, "wb").close()

        total_size = int(response.headers.get("content-length", 0)) + initial_pos

        with open(dest_path, "ab" if initial_pos > 0 else "wb") as f, tqdm(
            desc=file_name, total=total_size, unit="B",
            unit_scale=True, unit_divisor=1024,
            initial=initial_pos, ascii=False, miniters=1,
        ) as pbar:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))

        if expected_size > 0 and os.path.getsize(dest_path) != expected_size:
            raise ValueError(f"Size mismatch for {file_name}")
        if expected_md5 and self.calculate_md5(dest_path) != expected_md5:
            raise ValueError(f"MD5 mismatch for {file_name}")
        return True

    def download_file(
        self, url: str, dest_path: str, expected_size: str, expected_md5: str, retries: int = 3
    ) -> bool:
        exp_size  = int(expected_size) if str(expected_size).isdigit() else 0
        file_name = os.path.basename(dest_path)

        for attempt in range(retries):
            try:
                if self._check_existing_file(dest_path, exp_size, expected_md5):
                    sys.stdout.write(f"{file_name}: Already verified. Skipping.\n")
                    return True
                return self._perform_download_and_verify(url, dest_path, exp_size, expected_md5)
            except Exception as e:
                if attempt == retries - 1:
                    logging.error(f"Failed to download {file_name}: {e}")
                    return False
                time.sleep(10)
        return False

    def process_download_queue(self, tsv_path: str, output_dir: str) -> None:
        """Download all files listed in a metadata TSV using FASTQ_Path as destination."""
        os.makedirs(output_dir, exist_ok=True)
        with open(tsv_path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                self.download_file(
                    row["Download_Links"],
                    row["FASTQ_Path"],
                    row["File_Sizes"],
                    row["MD5_Checksums"],
                )

    def _empty_row(self) -> dict[str, str]:
        return dict.fromkeys(TSV_HEADER, "")


# ── ENCODE downloader ─────────────────────────────────────────────────────────
class EncodeDownloader(BaseDownloader):
    """Fetch metadata from the ENCODE REST API.

    Sample_Type is resolved from ENCODE's control_type field directly,
    so Entrez is not needed.
    """

    def _resolve_sample_type(self, file_data: dict, condition: str) -> str:
        acc      = file_data.get("accession", "")
        override = self.overrides.get(acc, "")
        if file_data.get("control_type"):
            return override or "Input"
        target       = file_data.get("target") or {}
        target_label = target.get("label", "") if isinstance(target, dict) else ""
        return infer_sample_type(condition + " " + target_label, override)

    def _extract_fastq_info(self, data: dict, condition: str) -> list[dict]:
        base_url = "https://www.encodeproject.org"
        rows: list[dict] = []

        for file in data.get("files", []):
            is_fastq = file.get("file_type") == "fastq" or file.get("file_format") == "fastq"
            if not (is_fastq and file.get("status") == "released"):
                continue

            rep       = file.get("replicate") or {}
            pair_val  = str(file.get("paired_end", ""))
            acc       = file.get("accession", "N/A")
            is_paired = pair_val in ("1", "2")
            acc_key   = f"{acc}_{pair_val}" if is_paired else acc

            row = self._empty_row()
            row.update({
                "Condition":      condition,
                "Accession":      acc_key,
                "Sample_Type":    self._resolve_sample_type(file, condition),
                "Pair_ID":        pair_val if is_paired else "",
                "Paired_End":     str(is_paired),
                "Replicate":      str(rep.get("biological_replicate_number", "N/A")),
                "FASTQ_Path":     str(RAW_DIR / f"{acc_key}.fastq.gz"),
                "Download_Links": base_url + file.get("href", ""),
                "File_Sizes":     str(file.get("file_size", "0")),
                "MD5_Checksums":  file.get("md5sum", ""),
            })
            rows.append(row)
        return rows

    def fetch(self, experiment_id: str) -> tuple[list[dict], str | None]:
        url       = f"https://www.encodeproject.org/experiments/{experiment_id}/?format=json"
        data, err = self.request_data(url, headers={"accept": "application/json"})
        if err or not isinstance(data, dict):
            return [], err or "Invalid response"

        target    = data.get("target") or {}
        label     = target.get("label", "") if isinstance(target, dict) else ""
        condition = f"{data.get('assay_term_name', 'Unknown')}_{label}".strip("_")

        rows = self._extract_fastq_info(data, condition)
        return (rows, None) if rows else ([], "No released FASTQs found")


# ── ENA downloader ─────────────────────────────────────────────────────────────
class EnaDownloader(BaseDownloader):
    """Fetch metadata from the ENA portal.

    Sample_Type is resolved via Entrez for each run accession.
    Each paired-end run produces two rows (one per FASTQ file).
    """

    def __init__(self, sample_type_overrides: dict[str, str] | None = None) -> None:
        super().__init__(sample_type_overrides)
        self._entrez = EntrezResolver()

    def _process_ena_entry(self, entry: dict, condition: str) -> list[dict]:
        """Processes a single entry from ENA API response and returns corresponding rows."""
        ftp_raw = entry.get("fastq_ftp", "") or entry.get("sra_ftp", "")
        if not ftp_raw:
            run_acc_log = entry.get("run_accession", "N/A")
            logging.warning(f"[{condition}] No download links for {run_acc_log}. Skipping.")
            return []

        bytes_raw = entry.get("fastq_bytes", "") or entry.get("sra_bytes", "")
        md5_raw = entry.get("fastq_md5", "") or entry.get("sra_md5", "")

        ftp_links = [f"http://{l}" if not l.startswith("http") else l for l in ftp_raw.split(";") if l]
        sizes = bytes_raw.split(";") if bytes_raw else []
        md5s = md5_raw.split(";") if md5_raw else []

        run_acc = entry.get("run_accession", "N/A")
        override = self.overrides.get(run_acc, "")
        sample_type = self._entrez.resolve(run_acc, override)
        replicate = entry.get("sample_accession", "N/A")
        logging.info(f"[Entrez] {run_acc} → {sample_type}")

        is_paired = entry.get("library_layout") == "PAIRED"
        if is_paired and len(ftp_links) == 2:
            return self._create_pe_rows(condition, sample_type, run_acc, replicate, ftp_links, sizes, md5s)

        return [self._create_se_row(condition, sample_type, run_acc, replicate, ftp_links, sizes, md5s)]

    def _create_pe_rows(
        self, condition: str, sample_type: str, run_acc: str, replicate: str,
        ftp_links: list[str], sizes: list[str], md5s: list[str]
    ) -> list[dict]:
        """Creates row data for paired-end reads."""
        rows = []
        for i, (link, pair_id) in enumerate(zip(ftp_links, ["1", "2"])):
            acc_key = f"{run_acc}_{pair_id}"
            row = self._empty_row()
            row.update({
                "Condition":      condition,
                "Accession":      acc_key,
                "Sample_Type":    sample_type,
                "Pair_ID":        pair_id,
                "Paired_End":     "True",
                "Replicate":      replicate,
                "FASTQ_Path":     str(RAW_DIR / f"{acc_key}.fastq.gz"),
                "Download_Links": link,
                "File_Sizes":     sizes[i] if i < len(sizes) else "",
                "MD5_Checksums":  md5s[i] if i < len(md5s) else "",
            })
            rows.append(row)
        return rows

    def _create_se_row(
        self, condition: str, sample_type: str, run_acc: str, replicate: str,
        ftp_links: list[str], sizes: list[str], md5s: list[str]
    ) -> dict:
        """Creates row data for a single-end read."""
        row = self._empty_row()
        row.update({
            "Condition":      condition,
            "Accession":      run_acc,
            "Sample_Type":    sample_type,
            "Pair_ID":        "",
            "Paired_End":     "False",
            "Replicate":      replicate,
            "FASTQ_Path":     str(RAW_DIR / f"{run_acc}.fastq.gz"),
            "Download_Links": ftp_links[0] if ftp_links else "",
            "File_Sizes":     sizes[0] if sizes else "",
            "MD5_Checksums":  md5s[0] if md5s else "",
        })
        return row

    def _extract_fastq_info(self, data: list[dict], condition: str) -> list[dict]:
        rows: list[dict] = []
        for entry in data:
            rows.extend(self._process_ena_entry(entry, condition))
        return rows

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
            return [], "Empty response"
        return self._extract_fastq_info(data, experiment_id), None


class GeoDownloader(EnaDownloader):
    """Download interface for GEO experiments (GSE/PRJ prefixes)."""


class SraDownloader(EnaDownloader):
    """Download interface for SRA runs (SRR/ERR prefixes)."""


class TcgaDownloader(BaseDownloader):
    """Placeholder for future TCGA GDC-Client integration."""

    def fetch(self, experiment_id: str) -> tuple[list[dict], str | None]:
        raise NotImplementedError("TCGA requires GDC-Client.") from None


DOWNLOADER_REGISTRY: dict[str, type[BaseDownloader]] = {
    "ENC": EncodeDownloader,
    "GSE": GeoDownloader,
    "PRJ": GeoDownloader,
    "SRR": SraDownloader,
    "ERR": SraDownloader,
    "TCG": TcgaDownloader,
}


# ── Entry point ───────────────────────────────────────────────────────────────
def main(target_ids: list[str] | None = None) -> None:
    """Run the download pipeline for the given experiment IDs.

    If target_ids is None, all experiments in config.yaml are processed.
    """
    experiments: dict = _config.get("downloads", {})
    ids = target_ids or list(experiments.keys())

    logging.info("Starting processing")

    for eid in ids:
        downloader_class = DOWNLOADER_REGISTRY.get(eid[:3].upper())
        if not downloader_class:
            logging.warning(f"[SKIP] {eid}: Unknown prefix.")
            continue

        exp_info  = experiments.get(eid, {})
        overrides = exp_info.get("overrides", {})

        try:
            downloader  = downloader_class(sample_type_overrides=overrides)
            rows, error = downloader.fetch(eid)

            if error:
                logging.error(f"[{eid}] {error}")
                continue
            if not rows:
                logging.warning(f"[{eid}] No valid data found.")
                continue

            tsv_path = TSV_DIR / f"{eid}_metadata.tsv"
            write_tsv(rows, tsv_path)
            logging.info(f"[{eid}] TSV → {tsv_path} ({len(rows)} rows)")

            logging.info(f"[{eid}] Downloading {len(rows)} file(s)...")
            downloader.process_download_queue(str(tsv_path), str(RAW_DIR))
            logging.info(f"[{eid}] Done")
            time.sleep(5)

        except NotImplementedError as e:
            logging.warning(f"[{eid}] {e}")
        except Exception as e:
            logging.critical(f"[{eid}] Critical failure: {e}")

    logging.info("Processing complete")


if __name__ == "__main__":
    main(sys.argv[1:] or None)