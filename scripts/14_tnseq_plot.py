"""
Visualise Tn-seq DESeq2 output for Paper I (conditional essentiality: HI Serum vs BHI).

Produces three plots:
  tnseq_plot_count_histogram.pdf  — distribution of raw counts per gene
  tnseq_plot_PCA.pdf              — PCA of samples on log2-normalised counts
  tnseq_plot_volcano.pdf          — volcano plot (conditionally essential genes)
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

REPO      = Path(__file__).resolve().parents[1]
DESEQ_DIR = REPO / "data" / "tmp" / "13_tnseq_deseq"
HTSEQ_DIR = REPO / "data" / "tmp" / "12_tnseq_htseq"

RESULTS_CSV = DESEQ_DIR / "Tnseq_DESeq2_results_HI_Serum_vs_BHI.csv"
NORM_CSV    = DESEQ_DIR / "Tnseq_DESeq2_normalized_counts.csv"

SAMPLES = {
    "HI_Serum_ERR1801009": "HI_Serum",
    "HI_Serum_ERR1801010": "HI_Serum",
    "HI_Serum_ERR1801011": "HI_Serum",
    "BHI_ERR1801012":      "BHI",
    "BHI_ERR1801013":      "BHI",
    "BHI_ERR1801014":      "BHI",
}

res  = pd.read_csv(RESULTS_CSV, index_col=0)
norm = pd.read_csv(NORM_CSV,    index_col=0)

def load_htseq(sample: str) -> pd.Series:
    path = HTSEQ_DIR / f"{sample}.counts.tsv"
    df = pd.read_csv(path, sep="\t", header=None, index_col=0)
    df = df[~df.index.str.startswith("__")]
    return df.iloc[:, 0].rename(sample)

raw = pd.concat([load_htseq(s) for s in SAMPLES], axis=1).fillna(0).astype(int)

# Q22 — histogram of total raw counts per gene
total_raw = raw.sum(axis=1)
expressed = total_raw[total_raw > 0]

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.hist(np.log10(expressed + 1), bins=60, color="#3a86ff", edgecolor="white", linewidth=0.3)
ax.axvline(np.log10(10 + 1), color="#e63946", linestyle="--", linewidth=1.2,
           label="count = 10 (rough expression threshold)")
ax.set_xlabel("log10(total raw count + 1)", fontsize=11)
ax.set_ylabel("Number of genes", fontsize=11)
ax.set_title("Tn-seq: Distribution of raw counts per gene (all samples summed)", fontsize=12)
ax.legend(fontsize=9)
sns.despine(ax=ax)
fig.tight_layout()
fig.savefig(DESEQ_DIR / "tnseq_plot_count_histogram.pdf")
plt.close(fig)
print("Wrote tnseq_plot_count_histogram.pdf")

# Q24 — PCA on log2-normalised counts
log_norm = np.log2(norm + 1)
X = StandardScaler().fit_transform(log_norm.T)

pca = PCA(n_components=2)
coords = pca.fit_transform(X)
pct = pca.explained_variance_ratio_ * 100

conditions = [SAMPLES[s] for s in norm.columns]
palette = {"BHI": "#2196f3", "HI_Serum": "#e91e63"}
offsets = {"BHI": (4, 4), "HI_Serum": (-4, 4)}

fig, ax = plt.subplots(figsize=(6, 5))
for cond, color in palette.items():
    idx = [i for i, c in enumerate(conditions) if c == cond]
    ax.scatter(coords[idx, 0], coords[idx, 1], color=color, s=90,
               label=cond, zorder=3)

for i, sample in enumerate(norm.columns):
    cond = conditions[i]
    dx, dy = offsets.get(cond, (4, 4))
    ax.annotate(sample, (coords[i, 0], coords[i, 1]),
                xytext=(dx, dy), textcoords="offset points",
                fontsize=7, ha="center")

ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
ax.axvline(0, color="grey", linewidth=0.5, linestyle="--")
ax.set_xlabel(f"PC1 ({pct[0]:.1f}% variance)", fontsize=11)
ax.set_ylabel(f"PC2 ({pct[1]:.1f}% variance)", fontsize=11)
ax.set_title("Tn-seq: PCA of samples (log2-normalised counts)", fontsize=12)
ax.legend(title="Condition", fontsize=9)
sns.despine(ax=ax)
fig.tight_layout()
fig.savefig(DESEQ_DIR / "tnseq_plot_PCA.pdf")
plt.close(fig)
print("Wrote tnseq_plot_PCA.pdf")

# Q25 — volcano plot: conditionally essential genes
df = res.copy()
df["neg_log10_pval"] = -np.log10(df["pvalue"].clip(lower=1e-300))

def classify(row):
    if pd.isna(row["padj"]) or pd.isna(row["log2FoldChange"]):
        return "Filtered / NA"
    if row["padj"] < 0.05 and row["log2FoldChange"] <= -1:
        return "Depleted in Serum (essential)"
    if row["padj"] < 0.05 and row["log2FoldChange"] >= 1:
        return "Enriched in Serum"
    return "Not significant"

df["category"] = df.apply(classify, axis=1)

cat_colors = {
    "Depleted in Serum (essential)": "#e63946",
    "Enriched in Serum":             "#1d6996",
    "Not significant":               "#adb5bd",
    "Filtered / NA":                 "#dee2e6",
}
cat_sizes = {
    "Depleted in Serum (essential)": 8,
    "Enriched in Serum":             8,
    "Not significant":               4,
    "Filtered / NA":                 3,
}

fig, ax = plt.subplots(figsize=(8, 6))
for cat in ["Filtered / NA", "Not significant", "Enriched in Serum", "Depleted in Serum (essential)"]:
    sub = df[df["category"] == cat]
    ax.scatter(sub["log2FoldChange"], sub["neg_log10_pval"],
               color=cat_colors[cat], s=cat_sizes[cat], alpha=0.55,
               linewidths=0, label=f"{cat} (n={len(sub)})", rasterized=True)

ax.axvline(-1, color="black", linestyle="--", linewidth=0.8)
ax.axvline( 1, color="black", linestyle="--", linewidth=0.8)
ax.axhline(-np.log10(0.05), color="black", linestyle=":", linewidth=0.8)
ax.set_xlabel("log2 fold change  (HI_Serum / BHI)", fontsize=11)
ax.set_ylabel("-log10(p-value)", fontsize=11)
ax.set_title("Tn-seq volcano: HI Serum vs BHI  (sorted by padj)", fontsize=12)
ax.legend(fontsize=8, markerscale=1.5, framealpha=0.7)
sns.despine(ax=ax)
fig.tight_layout()
fig.savefig(DESEQ_DIR / "tnseq_plot_volcano.pdf")
plt.close(fig)
print("Wrote tnseq_plot_volcano.pdf")

print(f"\nAll Tn-seq plots written to {DESEQ_DIR}")
