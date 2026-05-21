#!/bin/bash
#SBATCH --job-name=trim_serum_ERR1797971
#SBATCH --account=uppmax2026-1-61
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/trim_serum_ERR1797971-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/trim_serum_ERR1797971-%j.err

set -euo pipefail

REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
cd "$REPO_ROOT" || exit 1

module load Trimmomatic/0.39-Java-17

RAW_DIR=data/1_Zhang_2017/transcriptomics_data/RNA-Seq_Serum/raw
OUT_DIR=data/tmp/5_trim_rna/serum
mkdir -p "$OUT_DIR"

ADAPT="${EBROOTTRIMMOMATIC}/adapters/TruSeq3-PE.fa"

echo "Trimming ERR1797971..."
trimmomatic PE -threads 1 -phred33 \
  "${RAW_DIR}/ERR1797971_1.fastq.gz" "${RAW_DIR}/ERR1797971_2.fastq.gz" \
  "${OUT_DIR}/ERR1797971_1P.fq.gz" "${OUT_DIR}/ERR1797971_1U.fq.gz" \
  "${OUT_DIR}/ERR1797971_2P.fq.gz" "${OUT_DIR}/ERR1797971_2U.fq.gz" \
  ILLUMINACLIP:"${ADAPT}":2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 \
  MINLEN:36
echo "Done ERR1797971"
