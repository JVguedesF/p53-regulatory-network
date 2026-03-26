from __future__ import annotations

import gzip
import logging
import os
import shutil

import requests
from tqdm import tqdm

_GENOME_REGISTRY: dict[str, tuple[str, str]] = {
    "hg38": ("homo_sapiens", "GRCh38"),
    "mm39": ("mus_musculus",  "GRCm39"),
}


def _ensembl_url(species: str, assembly: str, release: str) -> str:
    """Build the Ensembl FTP URL for a primary-assembly FASTA."""
    genus, epithet = species.split("_", 1)
    fname = f"{genus.capitalize()}_{epithet}.{assembly}.dna.primary_assembly.fa.gz"
    base = (
        "https://ftp.ensembl.org/pub/current_fasta"
        if release == "current"
        else f"https://ftp.ensembl.org/pub/release-{release}/fasta"
    )
    return f"{base}/{species}/dna/{fname}"


def _remote_size(url: str) -> int:
    """Return Content-Length for url, or 0 when the header is absent."""
    try:
        resp = requests.head(url, timeout=15, allow_redirects=True)
        resp.raise_for_status()
        return int(resp.headers.get("Content-Length", 0))
    except requests.RequestException:
        return 0


def download_genome(
    target: str = "hg38",
    release: str = "current",
    out_dir: str = "data/genome/fasta",
    force: bool = False,
) -> str:
    """Download and decompress a reference genome FASTA from Ensembl FTP.

    The genome is downloaded once and decompressed in place. If the final
    FASTA already exists the function returns early unless force=True.

    Unlike FASTQ downloads, no retry or resume logic is applied here —
    the genome is a one-time download and can be re-run manually if needed.

    Returns the path to the decompressed FASTA file.
    """
    if target not in _GENOME_REGISTRY:
        raise ValueError(
            f"Unknown target '{target}'. Supported: {sorted(_GENOME_REGISTRY.keys())}"
        )

    species, assembly = _GENOME_REGISTRY[target]
    url = _ensembl_url(species, assembly, release)

    out_dir    = os.path.abspath(out_dir)
    fasta_path = os.path.join(out_dir, f"{target}.fa")
    gz_path    = fasta_path + ".gz"

    os.makedirs(out_dir, exist_ok=True)

    if os.path.exists(fasta_path) and not force:
        logging.info(f"[GENOME] Already exists — skipping: {fasta_path}")
        return fasta_path

    if force and os.path.exists(fasta_path):
        logging.info(f"[GENOME] --force enabled → removing {fasta_path}")
        os.remove(fasta_path)

    total_size = _remote_size(url)
    if not total_size:
        logging.warning("[GENOME] Could not determine remote file size.")

    logging.info(f"[GENOME] Downloading {target} from Ensembl...")

    response = requests.get(url, stream=True, timeout=60)
    response.raise_for_status()

    with open(gz_path, "wb") as f, tqdm(
        desc=os.path.basename(gz_path),
        total=total_size or None,
        unit="B", unit_scale=True, unit_divisor=1024,
        ascii=False,
    ) as pbar:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            if chunk:
                f.write(chunk)
                pbar.update(len(chunk))

    logging.info(f"[GENOME] Extracting {os.path.basename(gz_path)} ...")
    with gzip.open(gz_path, "rb") as src, open(fasta_path, "wb") as dst:
        shutil.copyfileobj(src, dst)

    os.remove(gz_path)
    logging.info(f"[GENOME] FASTA ready: {fasta_path}")
    return fasta_path