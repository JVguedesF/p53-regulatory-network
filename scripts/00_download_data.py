#!/usr/bin/env python3
from src.data_downloader import main

import argparse
import logging
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Download FASTQ files from ENCODE / ENA / GEO / SRA.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "accessions",
        nargs="*",
        metavar="ID",
        help="Experiment IDs to download. If omitted, uses all from config.yaml.",
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
    args = _build_parser().parse_args()
    main(args.accessions or None)