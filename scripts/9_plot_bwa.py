"""
Parse samtools flagstat and coverage output files and produce two plots:

  1. bwa_mapping_rates.png  — % primary mapped reads per sample (flagstat)
  2. coverage_per_contig.png — mean depth per contig for one representative
                               sample (BH_ERR1797972, samtools coverage)

Input:
  results/6_bwa/<sample>.flagstat.txt   — samtools flagstat output
  results/6_bwa/<sample>.coverage.txt   — samtools coverage output (TSV)

Output:
  results/plots/bwa_mapping_rates.png
  results/plots/coverage_per_contig.png
"""

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# ---------------------------------------------------------------------------
# paths
# ---------------------------------------------------------------------------
REPO     = Path(__file__).resolve().parents[1]
BWA_DIR  = REPO / "results" / "6_bwa"
OUT_DIR  = REPO / "results" / "plots"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Sample order for mapping-rate bar chart
# ---------------------------------------------------------------------------
SAMPLES = [
    ("BH_ERR1797972",    "BHI rep 1",   "BHI"),
    ("BH_ERR1797973",    "BHI rep 2",   "BHI"),
    ("BH_ERR1797974",    "BHI rep 3",   "BHI"),
    ("Serum_ERR1797969", "Serum rep 1", "Serum"),
    ("Serum_ERR1797970", "Serum rep 2", "Serum"),
    ("Serum_ERR1797971", "Serum rep 3", "Serum"),
]

CONDITION_COLOR = {"BHI": "#2196f3", "Serum": "#e91e63"}

# ---------------------------------------------------------------------------
# 1. Parse flagstat files  →  primary mapped %
# ---------------------------------------------------------------------------
# samtools flagstat line of interest:
#   27124574 + 0 primary mapped (98.59% : N/A)
FLAGSTAT_RE = re.compile(r"primary mapped \(([\d.]+)%")

labels   = []
pct_map  = []
colors   = []

for sample_id, label, cond in SAMPLES:
    fpath = BWA_DIR / f"{sample_id}.flagstat.txt"
    if not fpath.exists():
        print(f"WARNING: {fpath} not found — skipping")
        continue
    text = fpath.read_text()
    m = FLAGSTAT_RE.search(text)
    if not m:
        print(f"WARNING: could not parse mapping % from {fpath}")
        continue
    labels.append(label)
    pct_map.append(float(m.group(1)))
    colors.append(CONDITION_COLOR[cond])

# Plot 1 — mapping rates
fig, ax = plt.subplots(figsize=(9, 5))
bars = ax.bar(range(len(labels)), pct_map, color=colors, edgecolor="white", linewidth=0.5)

# Add value labels on bars
for bar, pct in zip(bars, pct_map):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.05,
            f"{pct:.2f}%", ha="center", va="bottom", fontsize=8.5)

ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=9)
ax.set_ylabel("Primary mapped reads (%)", fontsize=11)
ax.set_title("BWA-MEM: primary mapping rate per RNA-seq sample", fontsize=12)
ax.set_ylim(95, 100)

# Manual legend for conditions
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=c, label=cond)
                   for cond, c in CONDITION_COLOR.items()]
ax.legend(handles=legend_elements, fontsize=9)
sns.despine(ax=ax)
fig.tight_layout()
out1 = OUT_DIR / "bwa_mapping_rates.png"
fig.savefig(out1, dpi=150)
plt.close(fig)
print(f"Wrote {out1}")

# ---------------------------------------------------------------------------
# 2. Parse coverage file for representative sample  →  mean depth per contig
# ---------------------------------------------------------------------------
# samtools coverage columns:
# #rname  startpos  endpos  numreads  covbases  coverage  meandepth  meanbaseq  meanmapq
REP_SAMPLE = "BH_ERR1797972"
cov_path   = BWA_DIR / f"{REP_SAMPLE}.coverage.txt"

cov = pd.read_csv(cov_path, sep="\t", comment=None)
# Drop header-comment char if present
cov.columns = [c.lstrip("#") for c in cov.columns]

fig, ax = plt.subplots(figsize=(9, 5))
bar_colors = plt.cm.tab10.colors[:len(cov)]
ax.bar(cov["rname"], cov["meandepth"], color=bar_colors, edgecolor="white", linewidth=0.5)
ax.set_xlabel("Contig", fontsize=11)
ax.set_ylabel("Mean read depth (×)", fontsize=11)
ax.set_title(f"Mean sequencing depth per contig — {REP_SAMPLE} (BHI rep 1)", fontsize=12)
ax.tick_params(axis="x", rotation=30)

# Annotate bars with contig length
for idx, row in cov.iterrows():
    kb = row["endpos"] / 1000
    ax.text(idx, row["meandepth"] + 5,
            f"{kb:.0f} kb", ha="center", va="bottom", fontsize=7.5)

sns.despine(ax=ax)
fig.tight_layout()
out2 = OUT_DIR / "coverage_per_contig.png"
fig.savefig(out2, dpi=150)
plt.close(fig)
print(f"Wrote {out2}")
