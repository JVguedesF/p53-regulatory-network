#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import sys
import tempfile


def update_tsv(tsv_path: str, accession: str, **updates: str) -> None:
    """Update specific columns for a given accession in a pipeline TSV."""
    abs_path = os.path.abspath(tsv_path)
    tmp_fd, tmp_path = tempfile.mkstemp(dir=os.path.dirname(abs_path), suffix=".tmp")

    try:
        with (
            open(tmp_fd, "w", newline="", encoding="utf-8") as tmp_f,
            open(abs_path, newline="", encoding="utf-8") as src_f,
        ):
            reader = csv.DictReader(src_f, delimiter="\t")
            fieldnames = reader.fieldnames
            if not fieldnames:
                raise ValueError(f"TSV has no header: {tsv_path}")

            unknown = set(updates) - set(fieldnames)
            if unknown:
                raise ValueError(f"Unknown columns {unknown} — valid: {list(fieldnames)}")

            writer = csv.DictWriter(
                tmp_f,
                fieldnames=fieldnames,
                delimiter="\t",
                extrasaction="ignore",
                lineterminator="\n",
            )
            writer.writeheader()

            found = False
            for row in reader:
                if row["Accession"] == accession:
                    row.update(updates)
                    found = True
                writer.writerow(row)

        if not found:
            raise KeyError(f"Accession '{accession}' not found in {tsv_path}")

        os.replace(tmp_path, abs_path)

    except Exception:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise


def _cli() -> None:
    if len(sys.argv) < 4:
        print(
            "Usage: python src/tsv_updater.py <tsv_path> <accession> COL=VALUE ...",
            file=sys.stderr,
        )
        sys.exit(1)

    tsv_path  = sys.argv[1]
    accession = sys.argv[2]
    updates: dict[str, str] = {}

    for arg in sys.argv[3:]:
        if "=" not in arg:
            print(f"ERROR: expected COL=VALUE, got: {arg!r}", file=sys.stderr)
            sys.exit(1)
        col, _, val = arg.partition("=")
        updates[col.strip()] = val

    try:
        update_tsv(tsv_path, accession, **updates)
    except (KeyError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    _cli()