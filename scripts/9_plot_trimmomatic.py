"""
Parse Trimmomatic SLURM stderr logs and produce a stacked bar chart showing
read-pair survival per sample.

Input:  results/2_trimmomatic/*.err  (SLURM stderr from Trimmomatic jobs)
Output: results/plots/trimmomatic_survival.png

Each Trimmomatic stderr contains lines like:
  Input Read Pairs: 27078884 Both Surviving: 14397580 (53.17%)
  Forward Only Surviving: 11413645 (42.15%) Reverse Only Surviving: 0 (0.00%)
  Dropped: 1267659 (4.68%)

The script collects all such lines across all err files, maps them to named
samples by accession ID, and draws a 100%-stacked bar chart:
  - Both surviving (paired)
  - Forward only surviving (orphan)
  - Dropped
"""

import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------------------------------------------------------------------
# paths
# ---------------------------------------------------------------------------
REPO      = Path(__file__).resolve().parents[1]
TRIM_DIR  = REPO / "results" / "2_trimmomatic"
OUT_DIR   = REPO / "results" / "plots"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_FILE  = OUT_DIR / "trimmomatic_survival.png"

# ---------------------------------------------------------------------------
# Sample order and labels — fixed so the chart is reproducible
# ---------------------------------------------------------------------------
SAMPLE_ORDER = [
    ("ERR1797972", "BHI rep 1"),
    ("ERR1797973", "BHI rep 2"),
    ("ERR1797974", "BHI rep 3"),
    ("ERR1797969", "Serum rep 1"),
    ("ERR1797970", "Serum rep 2"),
    ("ERR1797971", "Serum rep 3"),
]

# ---------------------------------------------------------------------------
# Parse all .err files for Trimmomatic summary lines
# ---------------------------------------------------------------------------
# Pattern matches the single summary line Trimmomatic prints per run:
#   Input Read Pairs: NNN Both Surviving: NNN (xx.xx%) Forward Only Surviving: NNN
#   (xx.xx%) Reverse Only Surviving: NNN (xx.xx%) Dropped: NNN (xx.xx%)
PATTERN = re.compile(
    r"Input Read Pairs:\s+(\d+)\s+"
    r"Both Surviving:\s+(\d+)\s+\([\d.]+%\)\s+"
    r"Forward Only Surviving:\s+(\d+)\s+\([\d.]+%\)\s+"
    r"Reverse Only Surviving:\s+(\d+)\s+\([\d.]+%\)\s+"
    r"Dropped:\s+(\d+)"
)

# We search the entire stderr text, which also contains the command line that
# includes the accession ID.  We match acc IDs from the SAMPLE_ORDER list.
ACC_PATTERN = re.compile(r"(ERR\d{7})")

raw_stats = {}  # acc -> (total, both, fwd_only, rev_only, dropped)

for err_file in sorted(TRIM_DIR.glob("*.err")):
    text = err_file.read_text(errors="replace")
    # Find each Trimmomatic summary line and the accession that precedes it
    # (the command line always appears before the summary line in the log)
    lines = text.splitlines()
    current_acc = None
    for line in lines:
        acc_match = ACC_PATTERN.search(line)
        if acc_match:
            current_acc = acc_match.group(1)
        m = PATTERN.search(line)
        if m and current_acc:
            total    = int(m.group(1))
            both     = int(m.group(2))
            fwd_only = int(m.group(3))
            rev_only = int(m.group(4))
            dropped  = int(m.group(5))
            # Keep only if this accession is one we care about and not yet seen
            if current_acc not in raw_stats:
                raw_stats[current_acc] = (total, both, fwd_only + rev_only, dropped)

if len(raw_stats) < len(SAMPLE_ORDER):
    found = sorted(raw_stats.keys())
    needed = [acc for acc, _ in SAMPLE_ORDER]
    missing = [a for a in needed if a not in raw_stats]
    print(f"WARNING: parsed {len(raw_stats)} samples; missing: {missing}", file=sys.stderr)

# ---------------------------------------------------------------------------
# Build percentage arrays in SAMPLE_ORDER
# ---------------------------------------------------------------------------
labels    = []
pct_both  = []
pct_orphan = []
pct_drop  = []

for acc, label in SAMPLE_ORDER:
    if acc not in raw_stats:
        print(f"  Skipping {acc} — no data found", file=sys.stderr)
        continue
    total, both, orphan, dropped = raw_stats[acc]
    labels.append(label)
    pct_both.append(100 * both   / total)
    pct_orphan.append(100 * orphan / total)
    pct_drop.append(100 * dropped / total)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
x = range(len(labels))
fig, ax = plt.subplots(figsize=(9, 5))

bar_both   = ax.bar(x, pct_both,   color="#2196f3", label="Both pairs surviving")
bar_orphan = ax.bar(x, pct_orphan, bottom=pct_both, color="#ff9800", label="Forward only (orphan)")
bar_drop   = ax.bar(x, pct_drop,
                    bottom=[b + o for b, o in zip(pct_both, pct_orphan)],
                    color="#e53935", label="Dropped")

ax.set_xticks(list(x))
ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=9)
ax.set_ylabel("Percentage of input read pairs (%)", fontsize=11)
ax.set_title("Trimmomatic read-pair survival per sample", fontsize=12)
ax.set_ylim(0, 105)
ax.legend(loc="lower right", fontsize=9)
sns.despine(ax=ax)
fig.tight_layout()
fig.savefig(OUT_FILE, dpi=150)
plt.close(fig)
print(f"Wrote {OUT_FILE}")
