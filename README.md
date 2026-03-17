# p53-regulatory-network

> ⚠️ **This is a learning project, actively under development.** Methods, results, and documentation are subject to change as the work evolves.

Integration of ChIP-seq and RNA-seq to map the p53 regulatory network in cancer and its clinical relevance.

## Overview

This project integrates ChIP-seq and RNA-seq data to map the p53 regulatory network in tumor cell lines (e.g., MCF-7). Rather than analyzing DNA binding or transcription in isolation, the goal is to close the full regulatory loop:

**DNA binding (ChIP-seq) → transcriptional regulation (RNA-seq) → clinical outcome (TCGA survival analysis)**

This multi-omics approach distinguishes direct p53 targets — where the protein physically binds and alters expression — from indirect ones, providing a systemic view of p53 as a tumor suppressor and tumor microenvironment modulator.

## Repository Structure
```
p53-regulatory-network/
├── data/
│   ├── raw/          # Raw FASTQ files (not tracked by Git)
│   └── processed/    # BAM, BED, count matrices
├── scripts/          # Bash pipelines for primary processing
├── src/              # Python and R source code
├── notebooks/        # Jupyter Notebooks and RMarkdown for EDA
├── results/          # Figures, network files, survival statistics
└── docs/             # Methodology, parameter rationale, biological interpretation
```

## Tools

| Category | Tools |
|---|---|
| Genomic Processing | FastQC, Trim Galore, Bowtie2, STAR, Salmon, SAMtools, Picard |
| Peak Calling & Annotation | MACS2, ChIPseeker |
| Differential Expression | DESeq2 (R) |
| Motif Analysis | MEME Suite, custom probabilistic implementation (Python) |
| Integration & Networks | pandas, pybedtools, networkx |
| Survival & TCGA | TCGAbiolinks, survival, survminer (R) |
| Reproducibility | Nextflow, Docker |

## How to Run

### Option 1 — Automated (Docker + Nextflow)
Recommended for reproducibility on clusters and cloud servers.
*(Full container instructions and Nextflow commands coming soon.)*

### Option 2 — Manual (Conda)
Recommended for local development and step-by-step exploration.
```bash
# Clone the repository
git clone https://github.com/JVguedesF/p53-regulatory-network.git
cd p53-regulatory-network

# Create and activate the environment
conda env create -f environment.yml
conda activate p53-regulatory-network
```

Run the primary processing scripts in the order they appear in the `scripts/` directory.

## Key Findings

*This section will be updated as the project progresses with integrative findings, direct target validation, and clinical impact in the TCGA cohort.*