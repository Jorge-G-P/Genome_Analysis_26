#!/bin/bash
#SBATCH --job-name=plot_deseq
#SBATCH --account=uppmax2026-1-61
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/plot_deseq-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/plot_deseq-%j.err

set -euo pipefail

REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
cd "$REPO_ROOT" || exit 1

module load matplotlib/3.9.2-gfbf-2024a
module load Seaborn/0.13.2-gfbf-2024a
module load scikit-learn/1.6.1-gfbf-2024a

python3 scripts/10_plot_DESeq.py
