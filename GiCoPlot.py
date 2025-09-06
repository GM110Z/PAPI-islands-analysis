#!/usr/bin/env python3
"""
Island combo analysis from an existing presence/absence matrix.

Input
-----
A TSV/CSV where the first column is a genome/assembly ID (e.g., assembly/gcf/gca/id/accession)
and remaining columns are island indicators (0/1). Non-zero values are treated as 1.

Outputs
-------
* presence_absence_matrix.normalized.csv   – sanitized copy of the input matrix
* island_prevalence.csv                    – per-island presence counts
* island_combination_counts.csv            – all unique island combos + counts
* island_upsetplot_top{N}.pdf/.svg         – UpSet of the top-N multi-island combos
* island_richness_barplot.pdf/.svg         – bar chart of islands-per-genome

Usage
-----
python island_combo_from_matrix.py matrix.tsv --top-n 50 --drop-version \
  --base-font 11 --figsize 8 11 --orientation horizontal
"""

from __future__ import annotations
import argparse, re, sys
from pathlib import Path
from collections import Counter
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# ---------- Keep text editable in Illustrator ----------
mpl.rcParams['pdf.fonttype'] = 42           # TrueType (Type 42), NOT Type 3 outlines
mpl.rcParams['ps.fonttype']  = 42
mpl.rcParams['svg.fonttype'] = 'none'       # keep SVG text as text (not paths)
mpl.rcParams['text.usetex']  = False
mpl.rcParams['font.family']  = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

try:
    from upsetplot import UpSet, from_indicators
except Exception:
    sys.exit("ERROR: This script requires 'upsetplot'. Install with: pip install upsetplot")

ID_CANDIDATES = {"assembly","gcf","gca","genome","id","accession"}

def norm_id(s: str, drop_version: bool=False) -> str:
    if not isinstance(s, str):
        return s
    s = s.strip().strip("\"'").replace(" ", "_").upper()
    if drop_version:
        s = re.sub(r"^(G[CA]F_\d+)\.\d+$", r"\1", s)
    return s

def read_matrix(path: Path, id_col: str|None, drop_version: bool) -> tuple[pd.DataFrame, list[str]]:
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    df = pd.read_csv(path, sep=sep, dtype=str)
    # choose ID column
    if id_col is None:
        matches = [c for c in df.columns if c.strip().lower() in ID_CANDIDATES]
        id_col = matches[0] if matches else df.columns[0]
    df = df.rename(columns={id_col: "assembly"})
    # normalize IDs
    df["assembly"] = df["assembly"].map(lambda x: norm_id(x, drop_version=drop_version))
    # coerce island columns to 0/1
    island_cols = [c for c in df.columns if c != "assembly"]
    if not island_cols:
        sys.exit("ERROR: No island columns found (matrix must have >=2 columns).")
    for c in island_cols:
        v = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
        df[c] = (v > 0).astype(int)
    # aggregate duplicate assemblies (max)
    df = df.groupby("assembly", as_index=False)[island_cols].max()
    return df, island_cols

def write_normalized_matrix(df: pd.DataFrame, island_cols: list[str], outprefix: str):
    out = Path(f"{outprefix}presence_absence_matrix.normalized.csv")
    df[["assembly"] + island_cols].to_csv(out, index=False)
    print(f"✔ wrote {out}")

def write_prevalence(df: pd.DataFrame, island_cols: list[str], outprefix: str):
    prev = df[island_cols].sum(axis=0).sort_values(ascending=False)
    out = Path(f"{outprefix}island_prevalence.csv")
    prev.rename_axis("island").to_frame("n_present").to_csv(out)
    print(f"✔ wrote {out}")

def compute_combo_counts(df: pd.DataFrame, island_cols: list[str]) -> pd.DataFrame:
    bool_df = df[island_cols].astype(bool)
    combo_counter = Counter()
    for _, row in bool_df.iterrows():
        combo = tuple(sorted([isle for isle, present in row.items() if present]))
        if combo:  # skip empty
            combo_counter[combo] += 1
    rows = []
    for combo, count in combo_counter.items():
        rows.append({
            "combo": " + ".join(combo) if combo else "",
            "count": int(count),
            "combo_size": len(combo)
        })
    combo_df = pd.DataFrame(rows).sort_values("count", ascending=False)
    return combo_df

def plot_upset_from_combos(combo_df: pd.DataFrame, top_n: int, exclude_singletons: bool,
                           base_font: int, figsize: tuple[float,float], orientation: str,
                           outprefix: str):
    df = combo_df.copy()
    if exclude_singletons:
        df = df[df["combo_size"] > 1]
    if df.empty:
        print("! Skipping UpSet: no multi-island combos found.")
        return
    df = df.nlargest(top_n, "count")

    # Build indicator matrix and repeat rows by count
    rows = []
    for combo, count in zip(df["combo"], df["count"]):
        inds = {isle: True for isle in combo.split(" + ")}
        inds["count"] = int(count)
        rows.append(inds)
    ind_df = pd.DataFrame(rows).fillna(False)
    rep_df = ind_df.loc[ind_df.index.repeat(ind_df["count"])].drop(columns=["count"]).astype(bool)

    # Bigger, readable defaults: avoid having to scale down in Illustrator
    plt.rcParams.update({
        "font.size": base_font,
        "axes.labelsize": base_font,
        "axes.titlesize": base_font + 1,
        "xtick.labelsize": base_font,
        "ytick.labelsize": base_font,
        "legend.fontsize": base_font,
        "savefig.bbox": "tight"
    })

    fig = plt.figure(figsize=figsize)
    upset = UpSet(from_indicators(rep_df),
                  orientation=orientation,
                  subset_size="count",
                  show_counts=True,
                  sort_by="cardinality")
    upset.plot()

    pdf = Path(f"{outprefix}island_upsetplot_top{top_n}.pdf")
    svg = Path(f"{outprefix}island_upsetplot_top{top_n}.svg")
    plt.savefig(pdf)
    plt.savefig(svg)
    plt.close(fig)
    print(f"✔ wrote {pdf}")
    print(f"✔ wrote {svg}")

def plot_richness(df: pd.DataFrame, island_cols: list[str],
                  base_font: int, figsize: tuple[float,float], outprefix: str):
    counts = df[island_cols].astype(int).sum(axis=1)
    plt.rcParams.update({
        "font.size": base_font,
        "axes.labelsize": base_font,
        "axes.titlesize": base_font + 1,
        "xtick.labelsize": base_font,
        "ytick.labelsize": base_font,
        "legend.fontsize": base_font,
        "savefig.bbox": "tight"
    })
    fig = plt.figure(figsize=figsize)
    ax = counts.value_counts().sort_index().plot.bar()
    ax.set_xlabel("Islands per genome")
    ax.set_ylabel("Number of genomes")
    ax.set_title("Island richness across assemblies")
    pdf = Path(f"{outprefix}island_richness_barplot.pdf")
    svg = Path(f"{outprefix}island_richness_barplot.svg")
    plt.tight_layout()
    plt.savefig(pdf)
    plt.savefig(svg)
    plt.close(fig)
    print(f"✔ wrote {pdf}")
    print(f"✔ wrote {svg}")

def main():
    ap = argparse.ArgumentParser(description="Analyse island combinations from a presence/absence matrix.")
    ap.add_argument("matrix", type=Path, help="Presence/absence matrix (TSV or CSV). First col = genome/assembly ID.")
    ap.add_argument("--id-col", default=None, help="Name of the ID column (defaults to first or one of assembly/gcf/gca/id/accession).")
    ap.add_argument("--drop-version", action="store_true", help="Drop .version from GCF_/GCA_ (treat .1/.2 as same assembly).")
    ap.add_argument("--top-n", type=int, default=50, help="Top-N multi-island combos to display in the UpSet plot.")
    ap.add_argument("--keep-singletons", action="store_true", help="Include single-island combos in the UpSet (default: exclude).")
    ap.add_argument("--base-font", type=int, default=11, help="Base font size used in figures (larger helps Illustrator).")
    ap.add_argument("--figsize", nargs=2, type=float, default=(8.0, 11.0),
                    help="Figure size (inches): WIDTH HEIGHT. Use tall sizes for readable labels.")
    ap.add_argument("--orientation", choices=["horizontal","vertical"], default="horizontal",
                    help="UpSet orientation (horizontal puts combo labels on Y axis).")
    ap.add_argument("--out-prefix", default="", help="Optional output filename prefix.")
    args = ap.parse_args()

    df, island_cols = read_matrix(args.matrix, args.id_col, args.drop_version)

    # Save normalized inputs & prevalence
    write_normalized_matrix(df, island_cols, args.out_prefix)
    write_prevalence(df, island_cols, args.out_prefix)

    # Combo counts
    combo_df = compute_combo_counts(df, island_cols)
    combo_out = Path(f"{args.out_prefix}island_combination_counts.csv")
    combo_df.to_csv(combo_out, index=False)
    print(f"✔ wrote {combo_out} ({len(combo_df)} combos)")

    # UpSet (multi-island by default)
    plot_upset_from_combos(combo_df,
                           top_n=args.top_n,
                           exclude_singletons=not args.keep_singletons,
                           base_font=args.base_font,
                           figsize=tuple(args.figsize),
                           orientation=args.orientation,
                           outprefix=args.out_prefix)

    # Richness
    # (use a landscape size by default for the bar plot)
    plot_richness(df, island_cols,
                  base_font=max(10, args.base_font),
                  figsize=(max(args.figsize[0], 8.0), 4.5),
                  outprefix=args.out_prefix)

if __name__ == "__main__":
    main()
