#!/bin/bash
#SBATCH --job-name=fastqc_illumina
#SBATCH --account=uppmax2026-1-61
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/fastqc_illumina-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/fastqc_illumina-%j.err

set -euo pipefail

REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
cd "$REPO_ROOT" || exit 1

module load FastQC/0.12.1-Java-11

OUT_DIR=data/tmp/1_fastqc
mkdir -p "$OUT_DIR"

BH_DIR=data/1_Zhang_2017/transcriptomics_data/RNA-Seq_BH/raw
SERUM_DIR=data/1_Zhang_2017/transcriptomics_data/RNA-Seq_Serum/raw

# BH replicates (paired-end):    ERR1797972, ERR1797973, ERR1797974
# Serum replicates (paired-end): ERR1797969, ERR1797970, ERR1797971
fastqc \
  "${BH_DIR}/ERR1797972_1.fastq.gz"    "${BH_DIR}/ERR1797972_2.fastq.gz" \
  "${BH_DIR}/ERR1797973_1.fastq.gz"    "${BH_DIR}/ERR1797973_2.fastq.gz" \
  "${BH_DIR}/ERR1797974_1.fastq.gz"    "${BH_DIR}/ERR1797974_2.fastq.gz" \
  "${SERUM_DIR}/ERR1797969_1.fastq.gz" "${SERUM_DIR}/ERR1797969_2.fastq.gz" \
  "${SERUM_DIR}/ERR1797970_1.fastq.gz" "${SERUM_DIR}/ERR1797970_2.fastq.gz" \
  "${SERUM_DIR}/ERR1797971_1.fastq.gz" "${SERUM_DIR}/ERR1797971_2.fastq.gz" \
  --outdir "$OUT_DIR" \
  --threads "$SLURM_CPUS_PER_TASK"

echo "FastQC complete. Reports in ${OUT_DIR}"
