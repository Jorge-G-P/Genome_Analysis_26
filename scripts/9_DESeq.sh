#!/bin/bash
#SBATCH --job-name=deseq2
#SBATCH --account=uppmax2026-1-61
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/deseq2-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/deseq2-%j.err

set -x
set -euo pipefail

REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
cd "$REPO_ROOT" || exit 1

# DESeq2 via Bioconductor bundle (Student Manual — differential expression after counting).
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2

Rscript scripts/9_DESeq.R
