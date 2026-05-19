#!/bin/bash
#SBATCH --job-name=prokka_short
#SBATCH --account=uppmax2026-1-61
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/prokka-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/prokka-%j.err

set -x
set -euo pipefail

REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
cd "$REPO_ROOT" || exit 1

module load prokka/1.14.5-gompi-2024a

# Choose the names of the output files
prokka --outdir data/tmp/4_prokka_ann/ --prefix paper1 data/tmp/2_canu_assembly/paper1_e745.contigs.fasta

# Visualize it in Artemis
#art mydir/mygenome.gff