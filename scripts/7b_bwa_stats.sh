#!/bin/bash
#SBATCH --job-name=bwa_stats
#SBATCH --account=uppmax2026-1-61
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/bwa_stats-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/bwa_stats-%j.err

set -euo pipefail

REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
cd "$REPO_ROOT" || exit 1

module load SAMtools/1.22.1-GCC-13.3.0

BAM_DIR=data/tmp/7_bwa_rna
STATS_DIR=${BAM_DIR}/stats
mkdir -p "$STATS_DIR"

for bam in "${BAM_DIR}"/*.sorted.bam; do
  base=$(basename "$bam" .sorted.bam)
  echo "Processing ${base}..."

  # Q19 — % reads mapped
  samtools flagstat "$bam" > "${STATS_DIR}/${base}.flagstat.txt"

  # Q20 — per-contig coverage depth and breadth
  samtools coverage "$bam" > "${STATS_DIR}/${base}.coverage.txt"

  echo "Done ${base}"
done

echo "Stats written to ${STATS_DIR}"
