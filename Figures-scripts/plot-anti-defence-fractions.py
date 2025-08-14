#!/usr/bin/env python3
"""
Antidefence: hotspots vs rest (per hotspot file), with explicit de-duplication.

What it does
------------
- Reads two-column TSVs with: subtype, protein_in_syst
- EXPLICITLY de-duplicates identical (subtype, protein_in_syst) rows in both --all and --hot
  before counting (i.e., unique pairs are counted once)
- Plots a horizontal stacked bar: each hotspot file is a separate segment + Elsewhere (light grey)
- Exports:
    1) summary TSV with totals & fractions
    2) optional plot-matrix TSV with the exact bar heights (counts or fractions) for Prism

Usage example
-------------
%run antidefence_hotspots_split.py \
  --all antidef_all.tsv \
  --hot hotA.tsv hotB.tsv hotC.tsv \
  --out antidef_split.png \
  --summary antidef_summary.tsv \
  --plot-data antidef_plotmatrix.tsv \
  --normalise fraction \
  --min-occ 5 \
  --top 50
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

REQ_COLS = ["subtype", "protein_in_syst"]

# ---------- IO & cleaning ----------

def read_antidef(path: str) -> pd.DataFrame:
    """Read TSV and return only the required columns (subtype, protein_in_syst), cleaned."""
    df = pd.read_csv(path, sep="\t", header=0, dtype=str, comment="#")
    if all(c in df.columns for c in REQ_COLS):
        df = df[REQ_COLS].copy()
    else:
        if df.shape[1] < 2:
            raise ValueError(f"{path}: expected at least 2 columns; found {df.shape[1]}")
        df = df.iloc[:, :2].copy()
        df.columns = REQ_COLS
    for c in REQ_COLS:
        df[c] = df[c].astype(str).str.strip()
    df = df[(df["subtype"] != "") & (df["protein_in_syst"] != "")]
    # EXPLICIT de-duplication of identical pairs
    df = df.drop_duplicates(subset=REQ_COLS, keep="first")
    return df


# ---------- counting ----------

def group_key(df: pd.DataFrame, group_by: str) -> pd.Series:
    if group_by == "subtype":
        return df["subtype"]
    elif group_by == "pair":
        return df["subtype"] + " | " + df["protein_in_syst"]
    else:
        raise ValueError(f"Unknown group_by: {group_by}")

def build_counts_split(
    df_all: pd.DataFrame,
    hotspot_dfs: list[pd.DataFrame],
    hotspot_labels: list[str],
    group_by: str
) -> tuple[pd.DataFrame, list[str]]:
    """
    Build counts per group (subtype or pair) for each hotspot file,
    plus Elsewhere = total - sum(hotspots).
    Note: since no seqid is available, identical pairs across genomes are indistinguishable
    and are counted once after de-duplication.
    """
    g_all = group_key(df_all, group_by)
    total = g_all.value_counts().rename_axis("group").rename("total_count").to_frame()

    for label, df_hot in zip(hotspot_labels, hotspot_dfs):
        g_hot = group_key(df_hot, group_by)
        counts = g_hot.value_counts().rename_axis("group").rename(label).to_frame()
        total = total.join(counts, how="left")

    for col in total.columns:
        total[col] = total[col].fillna(0).astype(int)

    hotspot_cols = hotspot_labels
    total["Elsewhere"] = total["total_count"] - total[hotspot_cols].sum(axis=1)
    total["Elsewhere"] = total["Elsewhere"].clip(lower=0)

    total["hotspot_total"] = total[hotspot_cols].sum(axis=1)
    total["hotspot_fraction"] = np.where(
        total["total_count"] > 0,
        total["hotspot_total"] / total["total_count"],
        0.0
    )
    total["elsewhere_fraction"] = 1.0 - total["hotspot_fraction"]

    total.index.name = "group"
    return total, hotspot_cols


# ---------- exports ----------

def export_plot_matrix(df: pd.DataFrame, hotspot_cols: list[str], out_path: str, mode: str):
    """Write the exact values used for the bars (counts or fractions) for import into Prism."""
    cols = hotspot_cols + ["Elsewhere"]
    plot_df = df[cols].copy()
    if mode == "fraction":
        denom = df["total_count"].replace(0, np.nan)
        plot_df = plot_df.div(denom, axis=0).fillna(0)
    plot_df = plot_df.reset_index()  # 'group' column
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    plot_df.to_csv(out_path, sep="\t", index=False)


# ---------- plotting ----------

def plot_stack_split(
    df: pd.DataFrame,
    hotspot_cols: list[str],
    outfile: str,
    mode: str = "none",
    title: str = "Antidefence per hotspot vs elsewhere"
):
    """Stacked bar plot, per-hotspot segments + Elsewhere (color-blind palette, axes kept, no top/right frame)."""
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

    groups = df.index.tolist()
    H = len(groups)
    y = np.arange(H)

    fig_h = max(4.0, 0.35 * H + 1.5)
    fig, ax = plt.subplots(figsize=(10, fig_h))

    left = np.zeros(H, dtype=float)
    for col in all_cols:
        ax.barh(y, plot_df[col].values, left=left, label=col, color=palette[col])
        left += plot_df[col].values.astype(float)

    ax.set_yticks(y)
    ax.set_yticklabels(groups)
    ax.set_xlabel(xlabel)
    ax.set_title(title, pad=10)
    ax.legend(loc="upper right", frameon=False, ncol=2)

    # Keep axes, remove top/right box
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Annotate raw totals at bar ends
    vmax = float(max(left)) if H else 1.0
    for i, tot in enumerate(df["total_count"].values):
        xpos = float(left[i])
        ax.text(xpos + 0.01 * (vmax + 1), i, str(int(tot)), va="center", fontsize=8, color="black")

    fig.tight_layout()
    Path(outfile).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outfile, dpi=600)
    plt.close(fig)


# ---------- CLI ----------

def main():
    ap = argparse.ArgumentParser(description="Antidefence: hotspots vs rest (per hotspot file), with de-duplication.")
    ap.add_argument("--all", required=True, help="TSV: all antidefence hits (subtype, protein_in_syst)")
    ap.add_argument("--hot", required=True, nargs="+", help="One or more TSVs: hotspot hits (same columns)")
    ap.add_argument("--out", required=True, help="Output figure (png/pdf/svg)")
    ap.add_argument("--summary", required=True, help="Output TSV with per-group stats")
    ap.add_argument("--plot-data", help="Optional TSV to save the exact matrix used for plotting (for Prism)")
    ap.add_argument("--top", type=int, default=None, help="Top N groups by total occurrences")
    ap.add_argument("--min-occ", type=int, default=1, help="Minimum total occurrences to include")
    ap.add_argument("--normalise", choices=["none", "fraction"], default="none", help="Plot counts or fractions")
    ap.add_argument("--group-by", choices=["subtype", "pair"], default="subtype",
                    help="Bar grouping: 'subtype' or exact 'pair' (subtype|protein)")
    args = ap.parse_args()

    # Load & de-duplicate
    df_all = read_antidef(args.all)
    hotspot_dfs = [read_antidef(f) for f in args.hot]
    hotspot_labels = [Path(f).stem for f in args.hot]

    # Build counts
    df, hotspot_cols = build_counts_split(
        df_all=df_all,
        hotspot_dfs=hotspot_dfs,
        hotspot_labels=hotspot_labels,
        group_by=args.group_by
    )

    # Filter & sort
    df = df[df["total_count"] >= args.min_occ]
    df = df.sort_values(["hotspot_fraction", "total_count"], ascending=[False, False])
    if args.top is not None and args.top > 0:
        df = df.head(args.top)

    # Summary
    cols_for_summary = ["total_count"] + hotspot_cols + ["Elsewhere", "hotspot_fraction", "elsewhere_fraction"]
    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    df[cols_for_summary].to_csv(args.summary, sep="\t", index=True, header=True)

    # Plot-matrix (exact bar values)
    if args.plot_data:
        export_plot_matrix(df, hotspot_cols, args.plot_data, mode=args.normalise)

    # Plot
    title = "Antidefence per hotspot vs elsewhere" + (" (fraction)" if args.normalise == "fraction" else "")
    plot_stack_split(df, hotspot_cols, args.out, mode=args.normalise, title=title)


if __name__ == "__main__":
    main()
