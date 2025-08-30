# Robustness analysis for co-occurrence before vs after de-replication (cluster collapse)
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, spearmanr
from itertools import combinations

# ====== CONFIG / INPUT PATHS ======
MATRIX = "presence_absence.tsv"              # genome-level 0/1 matrix (index col = genome_id)
CLUST  = "cluster_map.tsv"                   # columns include genome_id, cluster_id
OUT_SY = "Table_Sy_cooccurrence_robustness.tsv"
OUT_LOST = "pairs_lost_after_cluster.tsv"
OUT_GAIN = "pairs_gained_after_cluster.tsv"

FDR_ALPHA = 0.05
# Prevalence filter (skip ultra-rare/ubiquitous features to avoid degenerate 2x2 tables)
MIN_PREV, MAX_PREV = 0.01, 0.99

# ====== LOAD ======
df_full = pd.read_csv(MATRIX, sep="\t", index_col=0)
# Coerce to ints (0/1)
df_full = df_full.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
df_full.index = df_full.index.astype(str).str.strip()

# Load cluster_map, try to find the right two columns
cm_raw = pd.read_csv(CLUST, sep="\t", dtype=str)
cm_raw = cm_raw.applymap(lambda x: x.strip() if isinstance(x, str) else x)

if {"genome_id","cluster_id"}.issubset(set(cm_raw.columns)):
    cluster_map = cm_raw[["genome_id","cluster_id"]].copy()
else:
    # fallback: take first two columns
    cluster_map = cm_raw.iloc[:, :2].copy()
    cluster_map.columns = ["genome_id","cluster_id"]
cluster_map = cluster_map.drop_duplicates().reset_index(drop=True)
cluster_map["genome_id"] = cluster_map["genome_id"].astype(str).str.strip()
cluster_map["cluster_id"] = cluster_map["cluster_id"].astype(str).str.strip()

# ====== HANDLE MISSING GENOMES: assign singleton clusters ======
missing = sorted(set(df_full.index) - set(cluster_map["genome_id"]))
if missing:
    print(f"Warning: {len(missing)} genomes absent in cluster_map; assigning singleton clusters.")
    existing = set(cluster_map["cluster_id"])
    base = "singleton_"
    next_i = 1
    add_rows = []
    for g in missing:
        cid = f"{base}{next_i}"
        while cid in existing:
            next_i += 1
            cid = f"{base}{next_i}"
        add_rows.append({"genome_id": g, "cluster_id": cid})
        existing.add(cid)
        next_i += 1
    cluster_map = pd.concat([cluster_map, pd.DataFrame(add_rows)], ignore_index=True)

# ====== COLLAPSE TO CLUSTER LEVEL (presence = any member) ======
dfc = (
    df_full.reset_index(names="genome_id")
           .merge(cluster_map, on="genome_id", how="left")
           .groupby("cluster_id")
           .max(numeric_only=True)
)

# ====== FEATURE LIST (consistent in both, optional prevalence filter) ======
def valid_features(df):
    prev = df.mean(axis=0)
    return prev[(prev>=MIN_PREV) & (prev<=MAX_PREV)].index.tolist()

features = sorted(set(valid_features(df_full)).intersection(valid_features(dfc)))
if len(features) < 2:
    raise ValueError("Too few features after filtering; relax MIN_PREV/MAX_PREV or check matrix.")

# ====== HELPERS ======
def pair_counts(df, f1, f2):
    x = df[f1].astype(bool); y = df[f2].astype(bool)
    a = int((x & y).sum())        # 1,1
    b = int((x & ~y).sum())       # 1,0
    c = int((~x & y).sum())       # 0,1
    d = int((~x & ~y).sum())      # 0,0
    return a,b,c,d

def haldane_or(a,b,c,d):
    # Haldane–Anscombe correction
    a2,b2,c2,d2 = a+0.5, b+0.5, c+0.5, d+0.5
    return (a2*d2)/(b2*c2)

def fisher_p(a,b,c,d):
    _, p = fisher_exact([[a,b],[c,d]], alternative="two-sided")
    return p

def bh_fdr(pvals):
    p = np.asarray(pvals, float)
    n = p.size
    order = np.argsort(p)
    ranked = np.empty(n, float)
    ranked[order] = p[order] * n / (np.arange(n)+1)
    # enforce monotonicity
    for i in range(n-2, -1, -1):
        ranked[order[i]] = min(ranked[order[i]], ranked[order[i+1]])
    return np.clip(ranked, 0, 1)

# ====== COMPUTE PAIRWISE CO-OCCURRENCES (FULL + CLUSTER) ======
from itertools import combinations
pairs = list(combinations(features, 2))
rows = []
for f1, f2 in pairs:
    a,b,c,d = pair_counts(df_full, f1, f2)
    ORf = haldane_or(a,b,c,d); pf = fisher_p(a,b,c,d)

    a2,b2,c2,d2 = pair_counts(dfc, f1, f2)
    ORc = haldane_or(a2,b2,c2,d2); pc = fisher_p(a2,b2,c2,d2)

    rows.append({
        "f1":f1,"f2":f2,
        "a_full":a,"b_full":b,"c_full":c,"d_full":d,
        "OR_full":ORf,"p_full":pf,
        "a_clust":a2,"b_clust":b2,"c_clust":c2,"d_clust":d2,
        "OR_cluster":ORc,"p_cluster":pc
    })
res = pd.DataFrame(rows)
res["q_full"]    = bh_fdr(res["p_full"].values)
res["q_cluster"] = bh_fdr(res["p_cluster"].values)
res["delta_log2OR"] = np.log2(res["OR_cluster"]) - np.log2(res["OR_full"])
res["sig_full"]    = res["q_full"]    < FDR_ALPHA
res["sig_cluster"] = res["q_cluster"] < FDR_ALPHA
res["concordant_sig"] = res["sig_full"] & res["sig_cluster"]

# ====== CONSISTENCY SUMMARY ======
full_sig = res[res.sig_full]
union_sig = res[res.sig_full | res.sig_cluster]

retained = (res.sig_full & res.sig_cluster).sum()
full_sig_n = int(full_sig.shape[0])
union_sig_n = int(union_sig.shape[0])

retention_pct = (100 * retained / full_sig_n) if full_sig_n > 0 else np.nan
jaccard = (retained / union_sig_n) if union_sig_n > 0 else np.nan

def sgn(x):
    x = np.asarray(x, float)
    return np.sign(np.log2(x))

same_dir_all = float((sgn(res['OR_full']) == sgn(res['OR_cluster'])).mean() * 100)
same_dir_sig = float((sgn(full_sig['OR_full']) == sgn(full_sig['OR_cluster'])).mean() * 100) if full_sig_n > 0 else np.nan

abs_delta = res['delta_log2OR'].abs().replace([np.inf, -np.inf], np.nan).dropna()
median_abs = float(abs_delta.median()) if len(abs_delta) else np.nan
p95_abs = float(abs_delta.quantile(0.95)) if len(abs_delta) else np.nan

mask = np.isfinite(res['OR_full']) & np.isfinite(res['OR_cluster']) & (res['OR_full']>0) & (res['OR_cluster']>0)
r_logOR = float(spearmanr(np.log2(res.loc[mask,'OR_full']), np.log2(res.loc[mask,'OR_cluster']))[0]) if mask.any() else np.nan
r_mlogp = float(spearmanr(-np.log10(res['p_full']), -np.log10(res['p_cluster']))[0])

summary = {
  "retention_pct": "NA" if pd.isna(retention_pct) else round(retention_pct,1),
  "jaccard": "NA" if pd.isna(jaccard) else round(jaccard,3),
  "same_dir_all_%": round(same_dir_all,1),
  "same_dir_sig_%": "NA" if pd.isna(same_dir_sig) else round(same_dir_sig,1),
  "median_|Δlog2OR|": "NA" if pd.isna(median_abs) else round(median_abs,3),
  "p95_|Δlog2OR|": "NA" if pd.isna(p95_abs) else round(p95_abs,3),
  "spearman_log2OR": "NA" if pd.isna(r_logOR) else round(r_logOR,3),
  "spearman_-log10p": "NA" if pd.isna(r_mlogp) else round(r_mlogp,3),
  "n_pairs_tested": len(res)
}

print("Consistency summary:", summary)

# ====== SAVE TABLES ======
res.sort_values(["q_full","q_cluster"]).to_csv(OUT_SY, sep="\t", index=False)

lost = res[ res.sig_full & (~res.sig_cluster) ][["f1","f2","OR_full","q_full","OR_cluster","q_cluster","delta_log2OR"]]
gained = res[ (~res.sig_full) & res.sig_cluster ][["f1","f2","OR_full","q_full","OR_cluster","q_cluster","delta_log2OR"]]
lost.to_csv(OUT_LOST, sep="\t", index=False)
gained.to_csv(OUT_GAIN, sep="\t", index=False)

print(f"Significant pairs (full): {res.sig_full.sum()} | retained after clustering: {retained}")
print(f"Lost after clustering: {lost.shape[0]} | Gained after clustering: {gained.shape[0]}")
