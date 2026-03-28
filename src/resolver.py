from __future__ import annotations

import logging
import os
import time
import xml.etree.ElementTree as ET

from Bio import Entrez

Entrez.email = os.environ.get("NCBI_EMAIL", "your-email@example.com")

INPUT_KEYWORDS = frozenset({"input", "control", "igg", "mock"})
RNA_KEYWORDS   = frozenset({"rna-seq", "rnaseq", "polya", "poly-a", "mrna"})

_cache: dict[str, str] = {}


def infer_sample_type(text: str, override: str = "") -> str:
    if override:
        return override
    lower = text.lower()
    if any(kw in lower for kw in RNA_KEYWORDS):
        return "RNA"
    if any(kw in lower for kw in INPUT_KEYWORDS):
        return "Input"
    return "IP"


def resolve_sample_type(run_accession: str, override: str = "") -> str:
    if override:
        return override
    if run_accession not in _cache:
        _cache[run_accession] = _query_entrez(run_accession)
    return _cache[run_accession]


def _query_entrez(run_accession: str) -> str:
    try:
        handle = Entrez.esearch(db="sra", term=run_accession, retmax=1)
        data   = Entrez.read(handle)
        handle.close()

        ids = list(data.get("IdList", [])) if isinstance(data, dict) else []
        if not ids:
            logging.warning(
                f"[Entrez] No record for {run_accession} "
                "— using fallback."
            )
            return infer_sample_type(run_accession)

        handle  = Entrez.efetch(db="sra", id=ids[0], rettype="full", retmode="xml")
        xml_raw = handle.read()
        handle.close()
        time.sleep(0.4)

        xml_str = xml_raw if isinstance(xml_raw, str) else xml_raw.decode("utf-8")
        return _parse_xml(xml_str, run_accession)

    except Exception as exc:
        logging.warning(
            f"[Entrez] Query failed for {run_accession}: {exc} "
            "— using fallback."
        )
        return infer_sample_type(run_accession)


def _is_input_chip(root: ET.Element) -> bool:
    for attr in root.iter("SAMPLE_ATTRIBUTE"):
        tag = (attr.findtext("TAG") or "").lower()
        val = (attr.findtext("VALUE") or "").lower()
        if any(kw in tag or kw in val for kw in INPUT_KEYWORDS):
            return True
    return False


def _parse_xml(xml_str: str, run_accession: str) -> str:
    try:
        root = ET.fromstring(xml_str)

        for strategy_el in root.iter("LIBRARY_STRATEGY"):
            strategy = (strategy_el.text or "").strip()

            if strategy == "RNA-Seq":
                return "RNA"

            if strategy == "ChIP-Seq":
                return "Input" if _is_input_chip(root) else "IP"

        logging.warning(
            f"[Entrez] Ambiguous XML for {run_accession} "
            "— using fallback."
        )
        return infer_sample_type(run_accession)

    except ET.ParseError as exc:
        logging.warning(
            f"[Entrez] XML parse error for {run_accession}: {exc}"
        )
        return infer_sample_type(run_accession)