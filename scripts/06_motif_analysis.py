#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("metadata_tsv", help="TSV updated by 05_peak_calling.sh")
    p.add_argument("--genome-fasta", required=True)
    p.add_argument("--peak-window",  type=int, default=100)
    p.add_argument("--motif-width",  type=int, default=20)
    p.add_argument("--n-motifs",     type=int, default=5)
    p.add_argument("--top-peaks",    type=int, default=500)
    p.add_argument("--meme-bin",     default=None)
    return p.parse_args()


def find_peak_files(metadata_tsv: Path) -> list[str]:
    if not metadata_tsv.exists():
        sys.exit(f"Metadata TSV not found: {metadata_tsv}")

    meta = pd.read_csv(metadata_tsv, sep="\t", dtype=str).fillna("")

    missing = {"Sample_Type", "Peak_File"} - set(meta.columns)
    if missing:
        sys.exit(f"Missing TSV columns: {', '.join(sorted(missing))}")

    ip_rows = meta[
        (meta["Sample_Type"].str.lower() == "ip") &
        (meta["Peak_File"].str.len() > 0)
    ]

    if ip_rows.empty:
        sys.exit("No IP rows with Peak_File found. Run 05_peak_calling.sh first.")

    peak_files = list(dict.fromkeys(ip_rows["Peak_File"].tolist()))

    missing_files = [f for f in peak_files if not Path(f).exists()]
    if missing_files:
        sys.exit("Peak files not found:\n" + "\n".join(f"  {f}" for f in missing_files))

    return peak_files


def load_peaks(peak_files: list[str], top_n: int) -> pd.DataFrame:
    frames = []

    for f in peak_files:
        df = pd.read_csv(f, sep="\t", header=None, comment="#")

        if df.shape[1] < 3:
            print(f"[warn] Skipping invalid peak file: {f}", file=sys.stderr)
            continue

        cols = ["chr", "start", "end", "name", "score", "strand",
                "signalValue", "pValue", "qValue", "peak"]
        df = df.iloc[:, :min(df.shape[1], len(cols))]
        df.columns = cols[:df.shape[1]]
        df["peak_file"] = f
        frames.append(df)

    if not frames:
        sys.exit("No valid peak files loaded.")

    peaks = pd.concat(frames, ignore_index=True)

    score_col = next((c for c in ["signalValue", "score"] if c in peaks.columns), None)
    if score_col:
        peaks[score_col] = pd.to_numeric(peaks[score_col], errors="coerce")
        peaks = peaks.sort_values(score_col, ascending=False, na_position="last")

    return peaks.head(top_n).reset_index(drop=True)


def get_summit(row) -> int:
    start = int(row["start"])
    end   = int(row["end"])

    if "peak" in row.index:
        try:
            offset = int(float(row["peak"]))
            if offset >= 0:
                return start + offset
        except (TypeError, ValueError):
            pass

    return (start + end) // 2


def extract_sequences(peaks: pd.DataFrame, genome_fasta: str, window: int) -> list[str]:
    try:
        import pysam
        fasta = pysam.FastaFile(genome_fasta)
    except Exception:
        fasta = None

    sequences = []

    for _, row in peaks.iterrows():
        chrom = str(row["chr"])
        mid   = get_summit(row)
        start = max(0, mid - window)
        end   = mid + window

        if fasta is not None and chrom in fasta.references:
            seq = fasta.fetch(chrom, start, end).upper()
        else:
            region = f"{chrom}:{start + 1}-{end}"
            try:
                result = subprocess.run(
                    ["samtools", "faidx", genome_fasta, region],
                    capture_output=True, text=True, check=True
                )
                lines = result.stdout.strip().split("\n")
                seq = "".join(lines[1:]).upper() if len(lines) > 1 else ""
            except (subprocess.CalledProcessError, FileNotFoundError):
                seq = ""

        if len(seq) >= 10:
            sequences.append(seq)

    if fasta is not None:
        fasta.close()

    return sequences


def write_fasta(sequences: list[str], path: Path) -> None:
    with open(path, "w") as fh:
        for i, seq in enumerate(sequences):
            fh.write(f">peak_{i}\n{seq}\n")


def run_meme(fasta_path: Path, out_dir: Path, meme_bin: str | None,
             width: int, n_motifs: int) -> None:
    meme_exe = str(Path(meme_bin) / "meme") if meme_bin else "meme"
    meme_out = out_dir / "meme_out"

    cmd = [
        meme_exe, str(fasta_path),
        "-oc",      str(meme_out),
        "-dna",
        "-mod",     "zoops",
        "-w",       str(width),
        "-nmotifs", str(n_motifs),
        "-minw",    str(max(6, width - 4)),
        "-maxw",    str(width + 4),
        "-revcomp",
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        if result.returncode != 0:
            print(f"[meme] Non-zero exit:\n{result.stderr[:500]}", file=sys.stderr)
            return

        print(f"[meme] Output → {meme_out}")

    except FileNotFoundError:
        sys.exit("[meme] MEME not found in PATH. Check environment or pass --meme-bin.")
    except subprocess.TimeoutExpired:
        sys.exit("[meme] MEME timed out after 600s.")


def main():
    args = parse_args()

    metadata_tsv  = Path(args.metadata_tsv)
    chip_dataset  = metadata_tsv.stem

    out_dir = Path("results") / "tables"  / "chipseq" / "motifs" / chip_dataset
    fig_dir = Path("results") / "figures" / "chipseq" / "motifs" / chip_dataset
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"[setup] Dataset      : {chip_dataset}")
    print(f"[setup] Genome FASTA : {args.genome_fasta}")

    peak_files = find_peak_files(metadata_tsv)
    print(f"[peaks] Peak files   : {len(peak_files)}")

    peaks = load_peaks(peak_files, args.top_peaks)
    print(f"[peaks] Top peaks    : {len(peaks)}")
    peaks.to_csv(out_dir / "top_peaks_used.tsv", sep="\t", index=False)

    print(f"[seqs] Extracting ±{args.peak_window}bp around summit")
    sequences = extract_sequences(peaks, args.genome_fasta, args.peak_window)
    print(f"[seqs] Valid sequences: {len(sequences)}")

    if len(sequences) < 10:
        sys.exit("Too few valid sequences. Check genome FASTA and peak coordinates.")

    fasta_path = out_dir / "peak_sequences.fa"
    write_fasta(sequences, fasta_path)

    print(f"[meme] Running MEME (width={args.motif_width}, nmotifs={args.n_motifs})")
    run_meme(fasta_path, out_dir, args.meme_bin, args.motif_width, args.n_motifs)

    print(f"[done] Tables  → {out_dir}")
    print(f"[done] Figures → {fig_dir}")


if __name__ == "__main__":
    main()