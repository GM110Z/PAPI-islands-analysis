#!/usr/bin/env python3
"""
Unified island‑analysis pipeline derived from three earlier scripts.

Outputs
=======
* **presence_absence_matrix.csv** – binary table (assemblies × islands)
* **island_combination_counts.csv** – all unique island combos + counts (from *Secondmat.py*)
* **island_upsetplot_top100.pdf** – UpSet of the 100 most frequent combos **excluding single‑island cases** (from *third.py*, modified)
* **island_richness_barplot.pdf** – bar chart showing how many genomes carry 1, 2, 3… islands (from *matrixupset.py*)

Optional: input TSV copies with an extra “island” column (use `--skip-add-column`).
"""
from __future__ import annotations
import argparse, csv, re
from pathlib import Path
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import pandas as pd
from upsetplot import UpSet, from_indicators

ISLAND_PATTERN = re.compile(r"((?:PA[GP]I|LESGI)\d+)", re.IGNORECASE)
TOP_N = 100  # number of combos for UpSet

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

def extract_label(fname:str)->str:
    m = ISLAND_PATTERN.search(fname)
    if not m:
        raise ValueError(f"No island label in {fname}")
    return m.group(1).upper()

def annotate_file(path:Path, label:str):
    out = path.with_name(f"{path.stem}_with_island.tsv")
    with path.open() as fin, out.open("w", newline="") as fout:
        rdr, wtr = csv.reader(fin, delimiter="\t"), csv.writer(fout, delimiter="\t")
        for row in rdr:
            if row:
                wtr.writerow(row + [label])

def build_matrix(mapping:dict[str,set[str]]):
    islands = sorted(mapping)
    assemblies = sorted({a for s in mapping.values() for a in s})
    table = [["assembly"] + islands]
    for asm in assemblies:
        table.append([asm] + ["1" if asm in mapping[i] else "0" for i in islands])
    return table

# -----------------------------------------------------------------------------
# Analysis & plotting
# -----------------------------------------------------------------------------

def analyse(matrix_csv:Path):
    df = pd.read_csv(matrix_csv)
    island_cols = df.columns[1:]
    bool_df = df[island_cols].astype(bool)

    # ------------ Full combination table  ------------
    combo_counter = Counter()
    for _, row in bool_df.iterrows():
        combo = tuple(sorted(isle for isle, present in row.items() if present))
        if combo:
            combo_counter[combo] += 1

    readable_combo_counts = Counter()
    for combo, count in combo_counter.items():
        combo_name = " + ".join(combo)
        readable_combo_counts[combo_name] = count

    combo_df = pd.DataFrame(
        [(combo, count, combo.count(" + ") + 1) for combo, count in readable_combo_counts.items()],
        columns=["combo", "count", "combo_size"]
    )
    combo_df.sort_values("count", ascending=False, inplace=True)
    combo_df.to_csv("island_combination_counts.csv", index=False)
    print(f"✔ island_combination_counts.csv saved ({len(combo_df)} combos)")

    # ------------ UpSet plot (exclude single‑island)  ------------
    multi_combos_df = combo_df[combo_df["combo_size"] > 1].nlargest(TOP_N, "count")

    rows = []
    for combo, count in zip(multi_combos_df["combo"], multi_combos_df["count"]):
        row_dict = {isle: True for isle in combo.split(" + ")}
        row_dict["count"] = count
        rows.append(row_dict)

    indicator_df = pd.DataFrame(rows).fillna(False)
    repeated_df = indicator_df.loc[indicator_df.index.repeat(indicator_df["count"])]
    upset_data = from_indicators(repeated_df.drop(columns=["count"]).astype(bool))

    plt.figure(figsize=(12, 6))
    UpSet(upset_data, subset_size="count", show_counts=True).plot()
    plt.savefig(f"island_upsetplot_top{TOP_N}.pdf")
    plt.close()

    # ------------ Richness bar plot (islands per genome) ------------
    df["Island_Count"] = bool_df.sum(axis=1)
    plt.figure(figsize=(8, 4))
    df["Island_Count"].value_counts().sort_index().plot.bar()
    plt.xlabel("Islands per genome")
    plt.ylabel("Number of genomes")
    plt.title("Island richness across assemblies")
    plt.tight_layout()
    plt.savefig("island_richness_barplot.pdf")
    plt.close()

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main(folder:Path, output_csv:Path, skip_ann:bool):
    mapping:dict[str,set[str]] = defaultdict(set)
    for p in folder.iterdir():
        if p.suffix.lower() not in {".tsv", ".txt"}:
            continue
        try:
            lbl = extract_label(p.name)
        except ValueError as e:
            print("[warning]", e)
            continue
        if not skip_ann:
            annotate_file(p, lbl)
        with p.open() as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if len(row) >= 3 and row[2]:
                    mapping[lbl].add(row[2].strip())

    if not mapping:
        raise SystemExit("No island files processed – aborting.")

    matrix = build_matrix(mapping)
    with output_csv.open("w", newline="") as fh:
        csv.writer(fh).writerows(matrix)
    print("✔ presence_absence_matrix.csv generated")

    analyse(output_csv)

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Island presence/absence + combination analysis")
    parser.add_argument("folder", type=Path, help="Folder with island TSV files")
    parser.add_argument("-o", "--output", type=Path, default=Path("presence_absence_matrix.csv"))
    parser.add_argument("--skip-add-column", action="store_true")
    args = parser.parse_args()
    main(args.folder, args.output, args.skip_add_column)
