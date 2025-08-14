#!/usr/bin/env python3
"""
Plot anti-phage systems in multiple hotspot files vs elsewhere.

Outputs:
--------
1) Stacked bar plot:
   - Each hotspot file is its own colored segment + Elsewhere (light grey)
   - Publication-quality styling (high DPI, top/right frame removed, axes kept)
   - Annotated with total counts per system (even in fraction mode)

2) Summary table (TSV):
   - total_count
   - counts for each hotspot file
   - Elsewhere
   - hotspot_fraction
   - elsewhere_fraction

3) Optional plot matrix (TSV) for Prism via --plot-data:
   - Exactly the values used for the bars (counts or fractions)

Example:
%run plot_defence_hotspots_split.py \
  --all all_systems.tsv \
  --hot hot1.tsv hot2.tsv hot3.tsv \
  --out defence_split.png \
  --summary defence_split.tsv \
  --plot-data defence_plotmatrix.tsv \
  --normalise fraction \
  --min-occ 5 \
  --top 50 \
  --count-mode per_seqid
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


def read_two_col(path: str, col1="seqid", col2="system") -> pd.DataFrame:
    """Read a TSV; return just [seqid, system], cleaned. Tolerates extra columns or headerless two-col files."""
    df = pd.read_csv(path, sep="\t", header=0, dtype=str, comment="#")
    if col1 in df.columns and col2 in df.columns:
        df = df[[col1, col2]].copy()
    else:
        if df.shape[1] < 2:
            raise ValueError(f"{path}: expected at least 2 columns; found {df.shape[1]}")
        df = df.iloc[:, :2].copy()
        df.columns = [col1, col2]
    df[col1] = df[col1].astype(str).str.strip()
    df[col2] = df[col2].astype(str).str.strip()
    df = df[(df[col1] != "") & (df[col2] != "")]
    return df


def maybe_dedupe(df: pd.DataFrame, count_mode: str) -> pd.DataFrame:
    """Avoid PADLOC component over-counting by default."""
    if count_mode == "per_seqid":
        return df.drop_duplicates(subset=["seqid", "system"], keep="first")
    return df  # rows


def build_counts_split(
    df_all: pd.DataFrame,
    hotspot_dfs: list[pd.DataFrame],
    hotspot_labels: list[str],
    count_mode: str
) -> tuple[pd.DataFrame, list[str]]:
    """Build counts per system for each hotspot file and Elsewhere."""
    df_all = maybe_dedupe(df_all, count_mode)
    total = df_all.groupby("system").size().rename("total_count").to_frame()

    for label, df_hot in zip(hotspot_labels, hotspot_dfs):
        df_hot = maybe_dedupe(df_hot, count_mode)
        counts = df_hot.groupby("system").size().rename(label).to_frame()
        total = total.join(counts, how="left")

    for col in total.columns:
        total[col] = total[col].fillna(0).astype(int)

    hotspot_cols = hotspot_labels
    total["Elsewhere"] = total["total_count"] - total[hotspot_cols].sum(axis=1)
    total["Elsewhere"] = total["Elsewhere"].clip(lower=0)

    total["hotspot_total"] = total[hotspot_cols].sum(axis=1)
    total["hotspot_fraction"] = np.where(total["total_count"] > 0,
                                         total["hotspot_total"] / total["total_count"], 0.0)
    total["elsewhere_fraction"] = 1.0 - total["hotspot_fraction"]

    return total, hotspot_cols


def export_plot_matrix(df: pd.DataFrame, hotspot_cols: list[str], out_path: str, mode: str):
    """Write the exact values used for the bars (counts or fractions) for import into Prism."""
    cols = hotspot_cols + ["Elsewhere"]
    plot_df = df[cols].copy()
    if mode == "fraction":
        denom = df["total_count"].replace(0, np.nan)
        plot_df = plot_df.div(denom, axis=0).fillna(0)
    plot_df = plot_df.reset_index(names="system")
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    plot_df.to_csv(out_path, sep="\t", index=False)


def plot_stack_split(df: pd.DataFrame, hotspot_cols: list[str], outfile: str, mode: str = "none"):
    """Stacked bar plot with per-hotspot segments + Elsewhere (publication style)."""
    mpl.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "legend.fontsize": 9
    })

    # Color-blind friendly palette (Wong) + fixed light grey for Elsewhere
    cb_palette = [
        "#E69F00", "#56B4E9", "#009E73", "#F0E442",
        "#0072B2", "#D55E00", "#CC79A7", "#999999",
        "#A6761D", "#666666"
    ]
    palette = {col: cb_palette[i % len(cb_palette)] for i, col in enumerate(hotspot_cols)}
    palette["Elsewhere"] = "#bdbdbd"

    all_cols = hotspot_cols + ["Elsewhere"]

    # Values to plot
    if mode == "fraction":
        plot_df = df[all_cols].div(df["total_count"].replace(0, np.nan), axis=0).fillna(0)
        xlabel = "Fraction of occurrences"
    else:
        plot_df = df[all_cols]
        xlabel = "Number of occurrences"

    systems = df.index.tolist()
    H = len(systems)
    y = np.arange(H)

    fig_h = max(4.0, 0.35 * H + 1.5)
    fig, ax = plt.subplots(figsize=(10, fig_h))

    left = np.zeros(H, dtype=float)
    for col in all_cols:
        ax.barh(y, plot_df[col].values, left=left, label=col, color=palette[col])
        left += plot_df[col].values.astype(float)

    ax.set_yticks(y)
    ax.set_yticklabels(systems)
    ax.set_xlabel(xlabel)
    ax.set_title("Anti-phage systems per hotspot vs elsewhere", pad=10)
    ax.legend(loc="upper right", frameon=False, ncol=2)

    # Keep axes (left/bottom), remove only top/right frame lines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Annotate total counts (always raw totals)
    vmax = float(max(left)) if H else 1.0
    for i, tot in enumerate(df["total_count"].values):
        xpos = float(left[i])
        ax.text(xpos + 0.01 * (vmax + 1), i, str(int(tot)), va="center", fontsize=8, color="black")

    fig.tight_layout()
    Path(outfile).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outfile, dpi=600)  # high-res for publication
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Plot anti-phage systems: split per-hotspot vs elsewhere.")
    ap.add_argument("--all", required=True, help="TSV: all systems genome-wide (seqid, system)")
    ap.add_argument("--hot", required=True, nargs="+", help="One or more TSVs: hotspot systems (seqid, system)")
    ap.add_argument("--out", required=True, help="Output figure path (png/pdf/svg)")
    ap.add_argument("--summary", required=True, help="Output TSV with per-system stats")
    ap.add_argument("--plot-data", help="Optional TSV to save the exact matrix used for plotting (for Prism)")
    ap.add_argument("--top", type=int, default=None, help="Top N systems by total occurrences")
    ap.add_argument("--min-occ", type=int, default=1, help="Minimum total occurrences to include")
    ap.add_argument("--normalise", choices=["none", "fraction"], default="none",
                    help="How to scale values in the plot")
    ap.add_argument("--count-mode", choices=["per_seqid", "rows"], default="per_seqid",
                    help="per_seqid = unique (seqid, system); rows = raw rows")
    args = ap.parse_args()

    # Load
    df_all = read_two_col(args.all, "seqid", "system")
    hotspot_dfs = [read_two_col(f, "seqid", "system") for f in args.hot]
    hotspot_labels = [Path(f).stem for f in args.hot]  # legend & summary labels

    # Build counts
    df, hotspot_cols = build_counts_split(df_all, hotspot_dfs, hotspot_labels, args.count_mode)

    # Filter & sort
    df = df[df["total_count"] >= args.min_occ]
    df = df.sort_values(["hotspot_fraction", "total_count"], ascending=[False, False])
    if args.top is not None and args.top > 0:
        df = df.head(args.top)

    # Summary table
    cols_for_summary = ["total_count"] + hotspot_cols + ["Elsewhere", "hotspot_fraction", "elsewhere_fraction"]
    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    df[cols_for_summary].to_csv(args.summary, sep="\t", index=True, header=True)

    # Plot-data matrix (exact bar heights)
    if args.plot_data:
        export_plot_matrix(df, hotspot_cols, args.plot_data, mode=args.normalise)

    # Plot
    plot_stack_split(df, hotspot_cols, args.out, mode=args.normalise)


if __name__ == "__main__":
    main()

