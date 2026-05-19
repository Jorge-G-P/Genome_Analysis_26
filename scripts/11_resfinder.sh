#!/bin/bash
#SBATCH --job-name=resfinder
#SBATCH --account=uppmax2026-1-61
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/resfinder-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/resfinder-%j.err

set -euo pipefail

REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
cd "$REPO_ROOT" || exit 1

REF=data/tmp/2_canu_assembly/paper1_e745.contigs.fasta
OUT_DIR=data/tmp/11_resfinder
mkdir -p "$OUT_DIR"

module load ResFinder/4.4.2-gompi-2024b 2>/dev/null || module load ResFinder 2>/dev/null || {
  echo "ResFinder module not found. Run: module spider resfinder"
  echo "Then update the module load line in this script."
  exit 1
}

python3 "$RESFINDER_PATH/resfinder.py" -i "$REF" -o "$OUT_DIR" -b "$OUT_DIR/blast"

cp "$REF" "${OUT_DIR}/reference.fasta"
echo "ResFinder analysis complete. Results in ${OUT_DIR}"
