#!/usr/bin/env python3
from src.utils import download_genome

import argparse
import logging
import os
import sys
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

_CONFIG_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "config", "config.yaml",
)


def _load_genome_defaults() -> dict:
    try:
        with open(_CONFIG_PATH, encoding="utf-8") as f:
            config = yaml.safe_load(f) or {}
        return config.get("genome", {})
    except FileNotFoundError:
        return {}


def _build_parser(defaults: dict) -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Download a reference genome FASTA from Ensembl FTP.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "-t", "--target",
        default=defaults.get("target", "hg38"),
        choices=["hg38", "mm39"],
        metavar="GENOME",
        help="Genome shorthand.",
    )
    p.add_argument(
        "-r", "--release",
        default=defaults.get("release", "current"),
        metavar="REL",
        help="Ensembl release number or 'current'.",
    )
    p.add_argument(
        "-o", "--out-dir",
        default=defaults.get("out_dir", "data/genome/fasta"),
        metavar="DIR",
        dest="out_dir",
        help="Output directory.",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Force re-download even if the FASTA already exists.",
    )
    return p


def _configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)-7s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stdout,
    )


if __name__ == "__main__":
    _configure_logging()
    defaults = _load_genome_defaults()
    args     = _build_parser(defaults).parse_args()

    try:
        fasta = download_genome(
            target=args.target,
            release=args.release,
            out_dir=args.out_dir,
            force=args.force,
        )
        logging.info(f"[GENOME] Done → {fasta}")
    except (ValueError, RuntimeError) as e:
        logging.error(f"[GENOME] {e}")
        sys.exit(1)