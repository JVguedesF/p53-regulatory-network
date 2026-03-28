from __future__ import annotations

import gzip
import logging
import shutil
from pathlib import Path

import requests
from tqdm import tqdm

_GENOME_REGISTRY: dict[str, tuple[str, str]] = {
    "hg38": ("homo_sapiens", "GRCh38"),
    "mm39": ("mus_musculus",  "GRCm39"),
}

def _ensembl_url(species: str, assembly: str, release: str) -> str:
    genus, epithet = species.split("_", 1)
    fname = f"{genus.capitalize()}_{epithet}.{assembly}.dna.primary_assembly.fa.gz"
    base = (
        "https://ftp.ensembl.org/pub/current_fasta"
        if release == "current"
        else f"https://ftp.ensembl.org/pub/release-{release}/fasta"
    )
    return f"{base}/{species}/dna/{fname}"

def download_genome(
    target: str = "hg38",
    release: str = "current",
    out_dir: str = "data/genome/fasta",
    force: bool = False,
) -> str:
    if target not in _GENOME_REGISTRY:
        raise ValueError(
            f"Unknown target '{target}'. Supported: {sorted(_GENOME_REGISTRY.keys())}"
        )

    species, assembly = _GENOME_REGISTRY[target]
    url = _ensembl_url(species, assembly, release)

    out_path = Path(out_dir).resolve()
    fasta_path = out_path / f"{target}.fa"
    gz_path = out_path / f"{target}.fa.gz"

    out_path.mkdir(parents=True, exist_ok=True)

    if fasta_path.exists():
        if not force:
            logging.info(f"[GENOME] Already exists — skipping: {fasta_path}")
            return str(fasta_path)
        logging.info(f"[GENOME] --force enabled → removing {fasta_path}")
        fasta_path.unlink()

    logging.info(f"[GENOME] Downloading {target} from Ensembl...")

    response = requests.get(url, stream=True, timeout=60)
    response.raise_for_status()

    content_len = int(response.headers.get("content-length", 0)) or None

    with open(gz_path, "wb") as f, tqdm(
        desc=gz_path.name,
        total=content_len,
        unit="B", unit_scale=True, unit_divisor=1024,
        ascii=False,
    ) as pbar:
        for chunk in response.iter_content(chunk_size=1_048_576):
            if chunk:
                f.write(chunk)
                pbar.update(len(chunk))

    logging.info(f"[GENOME] Extracting {gz_path.name} ...")
    with gzip.open(gz_path, "rb") as src, open(fasta_path, "wb") as dst:
        shutil.copyfileobj(src, dst)

    gz_path.unlink()
    logging.info(f"[GENOME] FASTA ready: {fasta_path}")
    return str(fasta_path)