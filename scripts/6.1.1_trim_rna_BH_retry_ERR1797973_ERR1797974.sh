#!/bin/bash
#SBATCH --job-name=trim_bh_rerun
#SBATCH --account=uppmax2026-1-61
#SBATCH --time=06:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/trim_rna-bh-rerun-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/trim_rna-bh-rerun-%j.err

set -euo pipefail

REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
MAIN_SH="${REPO_ROOT}/scripts/6_trim_rna_BH.sh"
[[ -f "$MAIN_SH" ]] || { echo "Missing: $MAIN_SH" >&2; exit 1; }

export TRIM_ACC_LIST="ERR1797973 ERR1797974"
exec bash "$MAIN_SH"
