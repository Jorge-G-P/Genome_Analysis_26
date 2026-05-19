#!/bin/bash
#SBATCH --job-name=quast
#SBATCH --account=uppmax2026-1-61
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/quast-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/quast-%j.err

set -euo pipefail

REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
cd "$REPO_ROOT" || exit 1

module load QUAST/5.3.0-gfbf-2024a

ASSEMBLY=data/tmp/2_canu_assembly/paper1_e745.contigs.fasta
OUT_DIR=data/tmp/3_quast
mkdir -p "$OUT_DIR"

if [ ! -f "$ASSEMBLY" ]; then
  echo "Assembly not found: ${ASSEMBLY}"
  echo "Run scripts/2_assembly_nogrid.sh first."
  exit 1
fi

# Run without a reference — evaluates contiguity metrics only (N50, L50, etc.).
# If you download the E745 reference genome from NCBI, add: -r <reference.fasta>
quast.py \
  --output-dir "$OUT_DIR" \
  --threads 1 \
  "$ASSEMBLY"

echo "QUAST done. Report in ${OUT_DIR}/report.html"
