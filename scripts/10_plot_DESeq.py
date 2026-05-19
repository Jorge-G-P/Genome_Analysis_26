"""
Visualise DESeq2 output for Paper I (Zhang et al., E. faecium Serum vs BH).

Produces three plots (course questions 22, 24, 25):
  plot_count_histogram.pdf  -- distribution of raw counts per gene   (Q22)
  plot_PCA.pdf              -- PCA of samples on log2-normalised counts (Q24)
  plot_volcano.pdf          -- volcano plot sorted by significance     (Q25)
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

# ---------------------------------------------------------------------------
# paths
# ---------------------------------------------------------------------------
REPO     = Path(__file__).resolve().parents[1]
DESEQ_DIR = REPO / "data" / "tmp" / "9_deseq"
HTSEQ_DIR = REPO / "data" / "tmp" / "8_htseq_counts"

RESULTS_CSV = DESEQ_DIR / "DESeq2_results_Serum_vs_BH.csv"
NORM_CSV    = DESEQ_DIR / "DESeq2_normalized_counts.csv"

SAMPLES = {
    "BH_ERR1797972":    "BH",
    "BH_ERR1797973":    "BH",
    "BH_ERR1797974":    "BH",
    "Serum_ERR1797969": "Serum",
    "Serum_ERR1797970": "Serum",
    "Serum_ERR1797971": "Serum",
}

# ---------------------------------------------------------------------------
# load data
# ---------------------------------------------------------------------------
res  = pd.read_csv(RESULTS_CSV, index_col=0)
norm = pd.read_csv(NORM_CSV,    index_col=0)

def load_htseq(sample: str) -> pd.Series:
    path = HTSEQ_DIR / f"{sample}.counts.tsv"
    df = pd.read_csv(path, sep="\t", header=None, index_col=0)
    df = df[~df.index.str.startswith("__")]
    return df.iloc[:, 0].rename(sample)

raw = pd.concat([load_htseq(s) for s in SAMPLES], axis=1).fillna(0).astype(int)

# ---------------------------------------------------------------------------
# Q22 — histogram of total raw counts per gene
# ---------------------------------------------------------------------------
total_raw = raw.sum(axis=1)
expressed = total_raw[total_raw > 0]

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.hist(np.log10(expressed + 1), bins=60, color="#3a86ff", edgecolor="white", linewidth=0.3)
ax.axvline(np.log10(10 + 1), color="#e63946", linestyle="--", linewidth=1.2,
           label="count = 10 (rough expression threshold)")
ax.set_xlabel("log10(total raw count + 1)", fontsize=11)
ax.set_ylabel("Number of genes", fontsize=11)
ax.set_title("Distribution of raw counts per gene (all samples summed)", fontsize=12)
ax.legend(fontsize=9)
sns.despine(ax=ax)
fig.tight_layout()
fig.savefig(DESEQ_DIR / "plot_count_histogram.pdf")
plt.close(fig)
print("Wrote plot_count_histogram.pdf")

# ---------------------------------------------------------------------------
# Q24 — PCA on log2-normalised counts
# ---------------------------------------------------------------------------
log_norm = np.log2(norm + 1)
X = StandardScaler().fit_transform(log_norm.T)   # samples x genes

pca = PCA(n_components=2)
coords = pca.fit_transform(X)
pct = pca.explained_variance_ratio_ * 100

conditions = [SAMPLES[s] for s in norm.columns]
palette = {"BH": "#2196f3", "Serum": "#e91e63"}
offsets = {"BH": (4, 4), "Serum": (-4, 4)}   # x/y text nudge in display units

fig, ax = plt.subplots(figsize=(6, 5))
for cond, color in palette.items():
    idx = [i for i, c in enumerate(conditions) if c == cond]
    ax.scatter(coords[idx, 0], coords[idx, 1], color=color, s=90,
               label=cond, zorder=3)

for i, sample in enumerate(norm.columns):
    cond = conditions[i]
    dx, dy = offsets[cond]
    ax.annotate(sample, (coords[i, 0], coords[i, 1]),
                xytext=(dx, dy), textcoords="offset points",
                fontsize=7, ha="center")

ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
ax.axvline(0, color="grey", linewidth=0.5, linestyle="--")
ax.set_xlabel(f"PC1 ({pct[0]:.1f}% variance)", fontsize=11)
ax.set_ylabel(f"PC2 ({pct[1]:.1f}% variance)", fontsize=11)
ax.set_title("PCA of samples (log2-normalised counts)", fontsize=12)
ax.legend(title="Condition", fontsize=9)
sns.despine(ax=ax)
fig.tight_layout()
fig.savefig(DESEQ_DIR / "plot_PCA.pdf")
plt.close(fig)
print("Wrote plot_PCA.pdf")

# ---------------------------------------------------------------------------
# Q25 — volcano plot (sorted by significance, coloured by direction)
# ---------------------------------------------------------------------------
df = res.copy()
df["neg_log10_pval"] = -np.log10(df["pvalue"].clip(lower=1e-300))

def classify(row):
    if pd.isna(row["padj"]) or pd.isna(row["log2FoldChange"]):
        return "Filtered / NA"
    if row["padj"] < 0.05 and row["log2FoldChange"] >= 1:
        return "Up in Serum"
    if row["padj"] < 0.05 and row["log2FoldChange"] <= -1:
        return "Down in Serum"
    return "Not significant"

df["category"] = df.apply(classify, axis=1)

cat_colors = {
    "Up in Serum":     "#e63946",
    "Down in Serum":   "#1d6996",
    "Not significant": "#adb5bd",
    "Filtered / NA":   "#dee2e6",
}
cat_sizes = {
    "Up in Serum": 8, "Down in Serum": 8,
    "Not significant": 4, "Filtered / NA": 3,
}

fig, ax = plt.subplots(figsize=(8, 6))
for cat in ["Filtered / NA", "Not significant", "Down in Serum", "Up in Serum"]:
    sub = df[df["category"] == cat]
    ax.scatter(sub["log2FoldChange"], sub["neg_log10_pval"],
               color=cat_colors[cat], s=cat_sizes[cat], alpha=0.55,
               linewidths=0, label=f"{cat} (n={len(sub)})", rasterized=True)

ax.axvline(-1, color="black", linestyle="--", linewidth=0.8)
ax.axvline( 1, color="black", linestyle="--", linewidth=0.8)
ax.axhline(-np.log10(0.05), color="black", linestyle=":", linewidth=0.8)
ax.set_xlabel("log2 fold change  (Serum / BH)", fontsize=11)
ax.set_ylabel("-log10(p-value)", fontsize=11)
ax.set_title("Volcano plot: Serum vs BH  (sorted by padj)", fontsize=12)
ax.legend(fontsize=8, markerscale=1.5, framealpha=0.7)
sns.despine(ax=ax)
fig.tight_layout()
fig.savefig(DESEQ_DIR / "plot_volcano.pdf")
plt.close(fig)
print("Wrote plot_volcano.pdf")

print(f"\nAll plots written to {DESEQ_DIR}")
