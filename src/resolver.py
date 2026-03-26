from __future__ import annotations

import logging
import os
import time

from Bio import Entrez

Entrez.email = os.environ.get("NCBI_EMAIL", "your-email@example.com")

INPUT_KEYWORDS = {"input", "control", "igg", "mock"}
RNA_KEYWORDS   = {"rna-seq", "rna_seq", "rnaseq", "polya", "poly-a", "mrna"}


def infer_sample_type(text: str, override: str = "") -> str:
    """Infer Sample_Type from a text string using keyword matching.

    Returns 'IP', 'Input', or 'RNA'. Override takes precedence.
    Used as fallback when Entrez is unavailable.
    """
    if override:
        return override
    text_lower = text.lower()
    if any(kw in text_lower for kw in RNA_KEYWORDS):
        return "RNA"
    if any(kw in text_lower for kw in INPUT_KEYWORDS):
        return "Input"
    return "IP"


class EntrezResolver:
    """Resolve Sample_Type for SRR/ERR accessions via NCBI Entrez SRA XML.

    Falls back to keyword inference if Entrez is unavailable.
    Results are cached to avoid redundant API calls within a run.
    """

    def __init__(self) -> None:
        self._cache: dict[str, str] = {}

    def resolve(self, run_accession: str, override: str = "") -> str:
        """Return Sample_Type for a run accession. Override takes precedence."""
        if override:
            return override
        if run_accession not in self._cache:
            self._cache[run_accession] = self._query_entrez(run_accession)
        return self._cache[run_accession]

    def _query_entrez(self, run_accession: str) -> str:
        """Fetch SRA XML for a run accession and parse Sample_Type."""
        try:
            handle      = Entrez.esearch(db="sra", term=run_accession, retmax=1)
            parsed_data = Entrez.read(handle)
            handle.close()

            ids = list(parsed_data.get("IdList", [])) if isinstance(parsed_data, dict) else []
            if not ids:
                logging.warning(f"[Entrez] No record for {run_accession}. Using fallback.")
                return infer_sample_type(run_accession)

            handle  = Entrez.efetch(db="sra", id=ids[0], rettype="full", retmode="xml")
            xml_raw = handle.read()
            handle.close()
            time.sleep(0.4)  # respect NCBI rate limit (3 req/s without API key)

            xml_str = xml_raw if isinstance(xml_raw, str) else xml_raw.decode("utf-8")
            return self._parse_xml(xml_str, run_accession)

        except Exception as e:
            logging.warning(f"[Entrez] Query failed for {run_accession}: {e}. Using fallback.")
            return infer_sample_type(run_accession)

    def _parse_xml(self, xml: str, run_accession: str) -> str:
        """Parse SRA XML string to determine Sample_Type.

        Priority: RNA-Seq strategy → input keywords in characteristics/title
        → ChIP-Seq strategy → keyword fallback.
        """
        try:
            xl = xml.lower()

            if "<library_strategy>rna-seq</library_strategy>" in xl:
                return "RNA"
            if any(f"<value>{kw}" in xl for kw in INPUT_KEYWORDS):
                return "Input"
            if any(f">{kw}<" in xl or f">{kw} " in xl or f" {kw}<" in xl for kw in INPUT_KEYWORDS):
                return "Input"
            if "<library_strategy>chip-seq</library_strategy>" in xl:
                return "IP"

            logging.warning(f"[Entrez] Ambiguous XML for {run_accession}. Using fallback.")
            return infer_sample_type(xl)

        except Exception as e:
            logging.warning(f"[Entrez] XML parse failed for {run_accession}: {e}")
            return infer_sample_type(run_accession)