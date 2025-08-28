#!/usr/bin/env python3
import argparse, re, os, csv, sys, glob
from collections import defaultdict, OrderedDict

ACC_RE = re.compile(r'(G[CA]F_\d+)(?:\.(\d+))?')  # accession with optional .version

def parse_args():
    ap = argparse.ArgumentParser(
        description="Build presence/absence matrix from GI hit files.")
    ap.add_argument("--inputs", nargs="+", required=True,
                    help="Files and/or globs (e.g. islands/*.tsv islands/*.txt).")
    ap.add_argument("--output", required=True, help="Output TSV (presence/absence).")
    ap.add_argument("--drop-version", action="store_true",
                    help="Normalize IDs to GCF_XXXXXX (drop .version).")
    ap.add_argument("--id-list", default=None,
                    help="Optional file with genome IDs (one per line) to include/ordering (e.g., tree tips).")
    ap.add_argument("--island-col", type=int, default=None,
                    help="0-based column index to take as island name (default: auto = last non-numeric token).")
    ap.add_argument("--minlen", type=int, default=0,
                    help="Optional: ignore hits shorter than this length (0 = keep all).")
    return ap.parse_args()

def norm_id(s, drop_version=False):
    m = ACC_RE.search(s)
    if not m: return None
    acc = m.group(1) if drop_version else (m.group(1) + ("" if m.group(2) is None else f".{m.group(2)}"))
    return acc

def smart_split(line):
    # robust to tabs / multiple spaces
    return re.split(r'\s+', line.strip())

def pick_island(cols, prefer_idx=None, fallback_name=None):
    cand = None
    if prefer_idx is not None and 0 <= prefer_idx < len(cols):
        cand = cols[prefer_idx]
    # If nothing, try last token that looks non-numeric
    if not cand:
        for tok in reversed(cols):
            t = tok.strip().strip("'\"")
            if t and not re.fullmatch(r'[0-9.]+', t):
                cand = t
                break
    if not cand:
        cand = fallback_name or "UNKNOWN"
    # compress weird chars to safe label
    cand = re.sub(r'[^A-Za-z0-9._+-]+', '_', cand)
    return cand

def main():
    a = parse_args()

    files = []
    for pat in a.inputs:
        hits = glob.glob(pat)
        if not hits and os.path.isfile(pat):
            hits = [pat]
        files.extend(hits)
    files = sorted(set(files))
    if not files:
        sys.exit("No input files matched.")

    presence = defaultdict(set)  # genome -> set(island)
    islands = OrderedDict()      # preserve discovery order

    for path in files:
        base = os.path.basename(path)
        stem, _ = os.path.splitext(base)
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                cols = smart_split(line)
                acc = None
                # Find accession anywhere in the line
                for c in cols:
                    acc = norm_id(c, drop_version=a.drop_version)
                    if acc: break
                if not acc:
                    continue  # no accession on this line

                # Optionally enforce min length if start/end present
                if a.minlen and len(cols) >= 5:
                    try:
                        start = int(cols[2]); end = int(cols[3])
                        if abs(end - start) < a.minlen:
                            continue
                    except Exception:
                        pass

                isl = pick_island(cols, prefer_idx=a.island_col, fallback_name=stem)
                islands.setdefault(isl, None)
                presence[acc].add(isl)

    # Make row order
    if a.id_list:
        with open(a.id_list, "r", encoding="utf-8", errors="ignore") as f:
            order = [l.strip() for l in f if l.strip()]
    else:
        order = sorted(presence.keys())

    # Column order
    island_cols = list(islands.keys())
    # Write matrix
    with open(a.output, "w", newline="", encoding="utf-8") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["gcf"] + island_cols)
        for gid in order:
            rowset = presence.get(gid, set())
            w.writerow([gid] + [1 if c in rowset else 0 for c in island_cols])

    # Stats to stderr
    sys.stderr.write(f"[INFO] files: {len(files)}\n")
    sys.stderr.write(f"[INFO] genomes with >=1 island: {len(presence)}\n")
    sys.stderr.write(f"[INFO] islands (columns): {len(island_cols)} -> {', '.join(island_cols[:10])}{' ...' if len(island_cols)>10 else ''}\n")
    sys.stderr.write(f"[INFO] wrote: {a.output}\n")

if __name__ == "__main__":
    main()

