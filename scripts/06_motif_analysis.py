#!/usr/bin/env python3

import argparse
import math
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

BASES = ["A", "C", "G", "T"]
BASE_IDX = {b: i for i, b in enumerate(BASES)}

P53_CONSENSUS = "RRRCWWGYYY" * 2


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("chip_dataset",
                   help="Dataset name matching results/annotation/<name>/")
    p.add_argument("--genome-fasta", required=True)
    p.add_argument("--peak-window", type=int, default=100)
    p.add_argument("--motif-width", type=int, default=20)
    p.add_argument("--n-motifs", type=int, default=5)
    p.add_argument("--gibbs-iterations", type=int, default=500)
    p.add_argument("--gibbs-restarts", type=int, default=10)
    p.add_argument("--top-peaks", type=int, default=500)
    p.add_argument("--meme-bin", default=None,
                   help="Path to MEME suite bin/ directory (optional)")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def find_annotation_files(chip_dataset):
    anno_dir = Path("pipeline_outputs") / "chipseq" / chip_dataset / "annotation"
    if not anno_dir.exists():
        sys.exit(f"Annotation directory not found: {anno_dir}\nRun 07_peak_annotation.R first.")
    files = list(anno_dir.rglob("*_annotation.tsv"))
    if not files:
        sys.exit(f"No annotation TSV files found under {anno_dir}")
    return files


def load_peaks(anno_files, top_n):
    frames = []
    for f in anno_files:
        df = pd.read_csv(f, sep="\t")
        df["source"] = f.parent.name
        frames.append(df)
    peaks = pd.concat(frames, ignore_index=True)

    score_col = next((c for c in ["V5", "score", "signalValue"] if c in peaks.columns), None)
    if score_col:
        peaks = peaks.sort_values(score_col, ascending=False)

    return peaks.head(top_n).reset_index(drop=True)


def extract_sequences(peaks, genome_fasta, window):
    try:
        import pysam
        fasta = pysam.FastaFile(genome_fasta)
        use_pysam = True
    except Exception:
        use_pysam = False
        fasta = None

    chr_col   = next((c for c in ["seqnames", "chr", "chrom", "Chr"] if c in peaks.columns), None)
    start_col = next((c for c in ["start", "Start", "chromStart"] if c in peaks.columns), None)
    end_col   = next((c for c in ["end", "End", "chromEnd"] if c in peaks.columns), None)

    if not all([chr_col, start_col, end_col]):
        sys.exit("Could not identify chr/start/end columns in peak annotation file.")

    sequences = []
    for _, row in peaks.iterrows():
        chrom = str(row[chr_col])
        mid   = int((int(row[start_col]) + int(row[end_col])) / 2)
        s     = max(0, mid - window)
        e     = mid + window

        if use_pysam:
            assert fasta is not None
            if chrom in fasta.references:
                seq = fasta.fetch(chrom, s, e).upper()
            else:
                seq = _fetch_samtools(genome_fasta, chrom, s, e)
        else:
            seq = _fetch_samtools(genome_fasta, chrom, s, e)

        if seq and len(seq) >= 10:
            sequences.append(seq)

    if use_pysam and fasta is not None:
        fasta.close()

    return sequences


def _fetch_samtools(fasta, chrom, start, end):
    region = f"{chrom}:{start + 1}-{end}"
    try:
        result = subprocess.run(
            ["samtools", "faidx", fasta, region],
            capture_output=True, text=True, check=True
        )
        lines = result.stdout.strip().split("\n")
        return "".join(lines[1:]).upper() if len(lines) > 1 else None
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def sequences_to_matrix(seqs):
    filtered = [s for s in seqs if all(b in BASE_IDX for b in s)]
    return filtered


def build_pwm(seqs, width, pseudocount=0.1):
    counts = np.full((4, width), pseudocount)
    valid  = 0
    for seq in seqs:
        if len(seq) < width:
            continue
        for j, b in enumerate(seq[:width]):
            if b in BASE_IDX:
                counts[BASE_IDX[b], j] += 1
        valid += 1
    freqs = counts / counts.sum(axis=0, keepdims=True)
    return freqs, valid


def score_pwm(seq, pwm):
    width = pwm.shape[1]
    if len(seq) < width:
        return -np.inf, 0
    bg = 0.25
    best_score = -np.inf
    best_pos   = 0
    for i in range(len(seq) - width + 1):
        s = 0.0
        for j, b in enumerate(seq[i:i + width]):
            if b not in BASE_IDX:
                s -= 2.0
            else:
                p = pwm[BASE_IDX[b], j]
                s += math.log2(max(p, 1e-9) / bg)
        if s > best_score:
            best_score = s
            best_pos   = i
    return best_score, best_pos


def gibbs_sampler(seqs, width, n_iter, rng):
    n = len(seqs)
    if n < 2:
        return None, -np.inf

    positions = [rng.randint(0, max(1, len(s) - width)) for s in seqs]
    best_pwm   = None
    best_score = -np.inf

    for _ in range(n_iter):
        idx = rng.randint(0, n - 1)

        other_seqs = [
            seqs[i][positions[i]:positions[i] + width]
            for i in range(n) if i != idx and len(seqs[i]) >= positions[i] + width
        ]

        if not other_seqs:
            continue

        pwm, _ = build_pwm(other_seqs, width)

        seq = seqs[idx]
        probs = []
        for p in range(max(1, len(seq) - width + 1)):
            sub = seq[p:p + width]
            if len(sub) < width:
                probs.append(1e-9)
                continue
            log_p = sum(
                math.log(max(pwm[BASE_IDX[b], j], 1e-9))
                for j, b in enumerate(sub) if b in BASE_IDX
            )
            probs.append(math.exp(log_p))

        total = sum(probs)
        if total == 0:
            continue
        probs = [p / total for p in probs]
        positions[idx] = rng.choices(range(len(probs)), weights=probs)[0]

        motif_seqs = [
            seqs[i][positions[i]:positions[i] + width]
            for i in range(n) if len(seqs[i]) >= positions[i] + width
        ]
        current_pwm, _ = build_pwm(motif_seqs, width)
        current_score  = sum(score_pwm(s, current_pwm)[0] for s in seqs)

        if current_score > best_score:
            best_score = current_score
            best_pwm   = current_pwm.copy()

    return best_pwm, best_score


def run_gibbs(seqs, width, n_iter, n_restarts, seed):
    rng        = random.Random(seed)
    best_pwm   = None
    best_score = -np.inf

    for _ in range(n_restarts):
        pwm, score = gibbs_sampler(seqs, width, n_iter, rng)
        if pwm is not None and score > best_score:
            best_score = score
            best_pwm   = pwm

    return best_pwm, best_score


def pwm_to_df(pwm):
    return pd.DataFrame(pwm, index=BASES)


def information_content(pwm):
    bg   = 0.25
    ic   = np.sum(pwm * np.log2(np.maximum(pwm, 1e-9) / bg), axis=0)
    return ic


def plot_logo(pwm, title, path):
    ic    = information_content(pwm)
    width = pwm.shape[1]
    colors = {"A": "#2ECC71", "C": "#3498DB", "G": "#F39C12", "T": "#E74C3C"}

    fig, ax = plt.subplots(figsize=(width * 0.5 + 1, 3))
    for pos in range(width):
        col_ic = ic[pos]
        if col_ic <= 0:
            continue
        sorted_bases = sorted(BASES, key=lambda b, p=pos: pwm[BASE_IDX[b], p])
        bottom = 0.0
        for base in sorted_bases:
            height = pwm[BASE_IDX[base], pos] * col_ic
            ax.bar(pos, height, bottom=bottom, color=colors[base],
                   width=0.8, linewidth=0)
            bottom += height

    patches = [mpatches.Patch(color=colors[b], label=b) for b in BASES]
    ax.legend(handles=patches, loc="upper right", fontsize=8, ncol=4)
    ax.set_xlim(-0.5, width - 0.5)
    ax.set_ylim(0, 2.0)
    ax.set_xticks(range(width))
    ax.set_xticklabels([str(i) for i in range(1, width + 1)], fontsize=7)
    ax.set_ylabel("bits")
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def score_all_peaks(seqs, pwm):
    scores = [score_pwm(s, pwm)[0] for s in seqs]
    return np.array(scores)


def plot_score_distribution(scores, title, path):
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.hist(scores[np.isfinite(scores)], bins=40, color="#8E44AD", edgecolor="white")
    ax.set_xlabel("PWM Score")
    ax.set_ylabel("Peaks")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run_meme(sequences, out_dir, meme_bin, width, n_motifs):
    meme_exe = os.path.join(meme_bin, "meme") if meme_bin else "meme"
    fasta_path = out_dir / "peaks_for_meme.fa"

    with open(fasta_path, "w") as fh:
        for i, seq in enumerate(sequences):
            fh.write(f">peak_{i}\n{seq}\n")

    meme_out = out_dir / "meme_out"
    cmd = [
        meme_exe, str(fasta_path),
        "-oc", str(meme_out),
        "-dna", "-mod", "zoops",
        "-w", str(width),
        "-nmotifs", str(n_motifs),
        "-minw", str(max(6, width - 4)),
        "-maxw", str(width + 4),
        "-revcomp",
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            print(f"[meme] Non-zero exit: {result.stderr[:300]}", file=sys.stderr)
            return False
        print(f"[meme] Output written to {meme_out}")
        return True
    except FileNotFoundError:
        print("[meme] MEME not found — skipping external motif search", file=sys.stderr)
        return False
    except subprocess.TimeoutExpired:
        print("[meme] MEME timed out", file=sys.stderr)
        return False


def main():
    args = parse_args()
    random.seed(args.seed)
    np.random.seed(args.seed)

    out_dir = Path("results") / "tables" / "chipseq" / "motifs" / args.chip_dataset
    fig_dir = Path("results") / "figures" / "chipseq" / "motifs" / args.chip_dataset
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"[setup] Dataset     : {args.chip_dataset}")
    print(f"[setup] Genome FASTA: {args.genome_fasta}")

    anno_files = find_annotation_files(args.chip_dataset)
    print(f"[peaks] Annotation files: {len(anno_files)}")

    peaks = load_peaks(anno_files, args.top_peaks)
    print(f"[peaks] Top peaks selected: {len(peaks)}")

    peaks.to_csv(out_dir / "top_peaks_used.tsv", sep="\t", index=False)

    print(f"[seqs] Extracting ±{args.peak_window}bp around summit")
    sequences = extract_sequences(peaks, args.genome_fasta, args.peak_window)
    sequences = sequences_to_matrix(sequences)
    print(f"[seqs] Valid sequences: {len(sequences)}")

    if len(sequences) < 10:
        sys.exit("Too few valid sequences for motif analysis. Check genome FASTA and peak coordinates.")

    with open(out_dir / "peak_sequences.fa", "w") as fh:
        for i, seq in enumerate(sequences):
            fh.write(f">peak_{i}\n{seq}\n")

    print(f"[gibbs] Running Gibbs sampler (width={args.motif_width}, "
          f"{args.gibbs_iterations} iter × {args.gibbs_restarts} restarts)")

    pwm, score = run_gibbs(
        sequences, args.motif_width,
        args.gibbs_iterations, args.gibbs_restarts, args.seed
    )

    if pwm is None:
        sys.exit("Gibbs sampler failed to converge.")

    print(f"[gibbs] Best score: {score:.2f}")

    pwm_df = pwm_to_df(pwm)
    pwm_df.to_csv(out_dir / "gibbs_pwm.tsv", sep="\t")

    ic = information_content(pwm)
    ic_df = pd.DataFrame({"position": range(1, len(ic) + 1), "ic_bits": ic})
    ic_df["total_ic"] = ic.sum()
    ic_df.to_csv(out_dir / "gibbs_pwm_ic.tsv", sep="\t", index=False)

    plot_logo(pwm, f"Gibbs motif — {args.chip_dataset}", fig_dir / "gibbs_motif_logo.pdf")

    peak_scores = score_all_peaks(sequences, pwm)
    score_df = pd.DataFrame({"peak_index": range(len(peak_scores)), "pwm_score": peak_scores})
    score_df.to_csv(out_dir / "peak_pwm_scores.tsv", sep="\t", index=False)

    plot_score_distribution(
        peak_scores, f"PWM score distribution — {args.chip_dataset}",
        fig_dir / "pwm_score_distribution.pdf"
    )

    print(f"[gibbs] Total IC = {ic.sum():.2f} bits")
    print(f"[gibbs] Median peak score = {np.nanmedian(peak_scores):.2f}")

    if args.meme_bin is not None or shutil.which("meme") is not None:
        print("[meme] Running external MEME")
        run_meme(sequences, out_dir, args.meme_bin, args.motif_width, args.n_motifs)

    print(f"[done] Results → {out_dir}")
    print(f"[done] Figures → {fig_dir}")


if __name__ == "__main__":
    main()