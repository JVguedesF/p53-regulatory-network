#!/usr/bin/env python3

import argparse
import json
import os
import sys
import urllib.request
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy.stats import hypergeom, fisher_exact

IMMUNE_GENE_SETS = {
    "checkpoint": [
        "PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "VSIR",
        "PDCD1LG2", "CD80", "CD86", "ICOS", "ICOSLG",
    ],
    "cytokines_receptors": [
        "IFNG", "TNF", "IL2", "IL6", "IL10", "IL12A", "IL12B",
        "IL15", "IL18", "TGFB1", "CXCL9", "CXCL10", "CXCL11",
        "CCL2", "CCL5", "IFNGR1", "IFNGR2", "IL2RA",
    ],
    "antigen_presentation": [
        "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1",
        "B2M", "TAP1", "TAP2", "TAPBP", "NLRC5",
    ],
    "innate_sensing": [
        "STING1", "CGAS", "IRF3", "IRF7", "MAVS", "DDX58",
        "IFIH1", "TLR3", "TLR4", "TLR9", "MYD88",
    ],
    "nk_cytotoxicity": [
        "KLRK1", "KLRD1", "NCR1", "NCR3", "MICA", "MICB",
        "ULBP1", "ULBP2", "ULBP3", "RAET1E",
    ],
    "t_cell_trafficking": [
        "CXCR3", "CCR5", "CXCR6", "CCR2", "SELL",
        "ITGAL", "ITGB2", "ICAM1", "VCAM1",
    ],
    "apoptosis_immune": [
        "FASLG", "FAS", "TRAIL", "TNFRSF10A", "TNFRSF10B",
        "CASP8", "CASP3", "CFLAR",
    ],
}

ALL_IMMUNE_GENES = sorted({g for genes in IMMUNE_GENE_SETS.values() for g in genes})

BACKGROUND_SIZE = 20_000


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--conserved",
                   default=os.path.join("pipeline_outputs", "multilinhagem", "conserved_targets_all_runs.tsv"))
    p.add_argument("--survival",
                   default=None,
                   help="Optional cox_results.tsv from script 11")
    p.add_argument("--cancer-type", default="TCGA-LUAD")
    p.add_argument("--fdr-survival", type=float, default=0.05)
    return p.parse_args()


def load_conserved(path):
    if not os.path.exists(path):
        sys.exit(f"File not found: {path}\nRun 10_multilinhagem.R first.")
    df = pd.read_csv(path, sep="\t")
    if "symbol" not in df.columns:
        sys.exit("Expected column 'symbol' not found.")
    return df["symbol"].dropna().unique().tolist()


def load_survival(path, fdr_threshold):
    if path is None or not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t")
    return df[df["padj"] < fdr_threshold]["gene"].tolist()


def classify_targets(target_genes, immune_sets):
    rows = []
    for gene in target_genes:
        matched = [cat for cat, members in immune_sets.items() if gene in members]
        rows.append({
            "gene": gene,
            "is_immune": len(matched) > 0,
            "immune_categories": ";".join(matched) if matched else "",
        })
    return pd.DataFrame(rows)


def enrichment_test(target_genes, immune_genes, background_size):
    k = len(set(target_genes) & set(immune_genes))
    n = len(target_genes)
    K = len(immune_genes)
    N = background_size

    table = [[k, n - k], [K - k, N - K - (n - k)]]
    table[1][1] = max(table[1][1], 0)
    pval   = float(fisher_exact(table, alternative="greater")[1])  # type: ignore[arg-type]
    pval_hg = float(hypergeom.sf(k - 1, N, K, n))

    return {
        "n_targets": n,
        "n_immune_background": K,
        "background_size": N,
        "overlap": k,
        "expected": round(n * K / N, 2),
        "fold_enrichment": round((k / n) / (K / N), 3) if K > 0 and n > 0 else 0,
        "pvalue_fisher": round(float(pval), 6),
        "pvalue_hypergeom": round(float(pval_hg), 6),
    }


def per_category_enrichment(target_genes, immune_sets, background_size):
    rows = []
    for cat, members in immune_sets.items():
        overlap = set(target_genes) & set(members)
        k = len(overlap)
        n = len(target_genes)
        K = len(members)
        N = background_size
        table = [[k, n - k], [K - k, max(N - K - (n - k), 0)]]
        pval = float(fisher_exact(table, alternative="greater")[1])  # type: ignore[arg-type]
        rows.append({
            "category": cat,
            "n_category": K,
            "overlap": k,
            "genes": ";".join(sorted(overlap)),
            "fold_enrichment": round((k / n) / (K / N), 3) if k > 0 and n > 0 else 0.0,
            "pvalue": round(float(pval), 6),
        })
    df = pd.DataFrame(rows).sort_values("pvalue")
    df["padj"] = np.minimum(df["pvalue"] * len(df) / (df.index + 1).values, 1.0)
    return df


def plot_category_enrichment(cat_df, fig_path):
    df = cat_df.sort_values("fold_enrichment", ascending=True).copy()
    colors = ["#E74C3C" if p < 0.05 else "#95A5A6" for p in df["pvalue"]]

    fig, ax = plt.subplots(figsize=(7, max(3, 0.5 * len(df))))
    ax.barh(df["category"], df["fold_enrichment"], color=colors, edgecolor="white")
    ax.axvline(1.0, color="grey", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Fold Enrichment")
    ax.set_title("Immune Category Enrichment of p53 Targets")
    sig_patch = mpatches.Patch(color="#E74C3C", label="p < 0.05")
    ns_patch  = mpatches.Patch(color="#95A5A6", label="n.s.")
    ax.legend(handles=[sig_patch, ns_patch], loc="lower right", fontsize=9)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)


def plot_overlap_heatmap(target_genes, immune_sets, fig_path):
    gene_list = sorted(set(target_genes) & set(ALL_IMMUNE_GENES))
    if not gene_list:
        return

    cats = list(immune_sets.keys())
    mat  = np.zeros((len(gene_list), len(cats)), dtype=int)

    for j, cat in enumerate(cats):
        for i, gene in enumerate(gene_list):
            if gene in immune_sets[cat]:
                mat[i, j] = 1

    fig, ax = plt.subplots(figsize=(max(5, len(cats) * 0.9), max(4, len(gene_list) * 0.35)))
    ax.imshow(mat, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels(cats, rotation=35, ha="right", fontsize=9)
    ax.set_yticks(range(len(gene_list)))
    ax.set_yticklabels(gene_list, fontsize=8)
    ax.set_title("p53 Immune Targets × Category Membership")
    fig.tight_layout()
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)


def plot_survival_immune(immune_targets, survival_genes, fig_path):
    if not survival_genes:
        return

    sets = {
        "Immune\nTargets":   set(immune_targets),
        "Survival\nSig.":    set(survival_genes),
    }
    both = sets["Immune\nTargets"] & sets["Survival\nSig."]

    sizes  = [len(s) for s in sets.values()]
    labels = list(sets.keys())

    fig, axes_arr = plt.subplots(1, 2, figsize=(10, 4), squeeze=False)
    ax0 = axes_arr[0, 0]
    ax1 = axes_arr[0, 1]

    ax0.bar(labels, sizes, color=["#2ECC71", "#E74C3C"], width=0.5, edgecolor="white")
    ax0.set_ylabel("Gene count")
    ax0.set_title("Immune vs. Survival-significant targets")

    ax1.axis("off")
    if both:
        ax1.text(0.5, 0.5,
                 "Immune + Survival significant:\n\n" + "\n".join(sorted(both)),
                 ha="center", va="center", fontsize=10,
                 bbox={"boxstyle": "round", "facecolor": "#FAD7A0", "alpha": 0.7})
    else:
        ax1.text(0.5, 0.5, "No overlap", ha="center", va="center", fontsize=11)

    fig.suptitle("Immune Targets with Survival Significance", fontsize=12)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)


def main():
    args = parse_args()

    out_dir = Path("results") / "tables" / "immune_targets"
    fig_dir = Path("results") / "figures" / "immune_targets"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    target_genes   = load_conserved(args.conserved)
    survival_genes = load_survival(args.survival, args.fdr_survival)

    print(f"[targets]  Conserved p53 targets: {len(target_genes)}")
    print(f"[immune]   Immune gene database size: {len(ALL_IMMUNE_GENES)}")

    classified_df = classify_targets(target_genes, IMMUNE_GENE_SETS)
    immune_targets = classified_df[classified_df["is_immune"]]["gene"].tolist()

    print(f"[overlap]  Immune targets: {len(immune_targets)} / {len(target_genes)}")

    classified_df.to_csv(out_dir / "target_immune_classification.tsv", sep="\t", index=False)

    enrichment = enrichment_test(target_genes, ALL_IMMUNE_GENES, BACKGROUND_SIZE)
    with open(out_dir / "global_enrichment.json", "w") as fh:
        json.dump(enrichment, fh, indent=2)

    print(f"[enrichment] Overlap={enrichment['overlap']}, "
          f"FE={enrichment['fold_enrichment']}, "
          f"p={enrichment['pvalue_fisher']:.4g}")

    cat_df = per_category_enrichment(target_genes, IMMUNE_GENE_SETS, BACKGROUND_SIZE)
    cat_df.to_csv(out_dir / "category_enrichment.tsv", sep="\t", index=False)

    if survival_genes is not None:
        both = set(immune_targets) & set(survival_genes)
        priority_df = classified_df[classified_df["gene"].isin(both)].copy()
        priority_df["survival_significant"] = True
        priority_df.to_csv(out_dir / "priority_immune_survival_targets.tsv", sep="\t", index=False)
        print(f"[priority] Immune + survival-significant: {len(both)}")

    plot_category_enrichment(cat_df, fig_dir / "category_enrichment.pdf")
    plot_overlap_heatmap(target_genes, IMMUNE_GENE_SETS, fig_dir / "immune_target_heatmap.pdf")
    plot_survival_immune(immune_targets, survival_genes, fig_dir / "immune_survival_overlap.pdf")

    print(f"[done] Tables  → {out_dir}")
    print(f"[done] Figures → {fig_dir}")


if __name__ == "__main__":
    main()