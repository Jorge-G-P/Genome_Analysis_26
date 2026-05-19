#!/bin/bash
#SBATCH --job-name=trim_rna_serum
#SBATCH --account=uppmax2026-1-61
#SBATCH --time=06:30:30
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/trim_rna-serum-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/trim_rna-serum-%j.err

set -x
set -euo pipefail

REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
cd "$REPO_ROOT" || exit 1

module load Trimmomatic/0.39-Java-17

RAW_DIR=data/1_Zhang_2017/transcriptomics_data/RNA-Seq_Serum/raw
OUT_DIR=data/tmp/5_trim_rna/serum
mkdir -p "$OUT_DIR"

ADAPT="${EBROOTTRIMMOMATIC}/adapters/TruSeq3-PE.fa"

# Space-separated SRA ERR ids (default: all three). Override, e.g.:
#   sbatch --export=ALL,TRIM_ACC_LIST=ERR1797971 scripts/6.2_trim_rna_Serum.sh
TRIM_ACC_LIST=${TRIM_ACC_LIST:-ERR1797969 ERR1797970 ERR1797971}

# Paired-end Illumina RNA libraries in RNA-Seq_Serum/raw.
for acc in $TRIM_ACC_LIST; do
  trimmomatic PE -threads 4 -phred33 \
    "${RAW_DIR}/${acc}_1.fastq.gz" "${RAW_DIR}/${acc}_2.fastq.gz" \
    "${OUT_DIR}/${acc}_1P.fq.gz" "${OUT_DIR}/${acc}_1U.fq.gz" \
    "${OUT_DIR}/${acc}_2P.fq.gz" "${OUT_DIR}/${acc}_2U.fq.gz" \
    ILLUMINACLIP:"${ADAPT}":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 \
    MINLEN:36
done
