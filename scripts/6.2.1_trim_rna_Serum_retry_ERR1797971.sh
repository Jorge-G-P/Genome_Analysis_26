#!/bin/bash
#SBATCH --job-name=trim_serum_rerun
#SBATCH --account=uppmax2026-1-61
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/trim_rna-serum-rerun-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/trim_rna-serum-rerun-%j.err

set -euo pipefail

# Slurm may run a *copy* of this script from /var/spool/slurmd/...; do not use
# dirname "${BASH_SOURCE[0]}" to find the workload script (companion file is not copied).
REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
MAIN_SH="${REPO_ROOT}/scripts/6.2_trim_rna_Serum.sh"
[[ -f "$MAIN_SH" ]] || { echo "Missing: $MAIN_SH" >&2; exit 1; }

export TRIM_ACC_LIST=ERR1797971
exec bash "$MAIN_SH"
