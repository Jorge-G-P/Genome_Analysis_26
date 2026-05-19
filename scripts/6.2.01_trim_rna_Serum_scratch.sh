#!/bin/bash
#SBATCH --job-name=trim_rna_serum_scr
#SBATCH --account=uppmax2026-1-61
#SBATCH --time=06:30:30
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/trim_rna-serum-scr-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/trim_rna-serum-scr-%j.err

set -x
set -euo pipefail

REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
cd "$REPO_ROOT" || exit 1

module load Trimmomatic/0.39-Java-17

RAW_DIR=data/1_Zhang_2017/transcriptomics_data/RNA-Seq_Serum/raw
OUT_DIR=data/tmp/5_trim_rna/serum
mkdir -p "$OUT_DIR"

ADAPT="${EBROOTTRIMMOMATIC}/adapters/TruSeq3-PE.fa"

# Space-separated SRA ERR ids (default: all three). Example rerun one library:
#   sbatch --export=ALL,TRIM_ACC_LIST=ERR1797971 scripts/6.2_trim_rna_Serum_scratch.sh
TRIM_ACC_LIST=${TRIM_ACC_LIST:-ERR1797969 ERR1797970 ERR1797971}

if [[ -n "${SNIC_TMP:-}" && -d "$SNIC_TMP" && -w "$SNIC_TMP" ]]; then
  WORKDIR="${SNIC_TMP}/trim_rna_serum_${SLURM_JOB_ID}"
  mkdir -p "$WORKDIR"
else
  WORKDIR=""
  echo "WARN: SNIC_TMP not set or not writable; writing outputs directly to ${OUT_DIR}" >&2
fi

for acc in $TRIM_ACC_LIST; do
  if [[ -n "$WORKDIR" ]]; then
    trimmomatic PE -threads 4 -phred33 \
      "${RAW_DIR}/${acc}_1.fastq.gz" "${RAW_DIR}/${acc}_2.fastq.gz" \
      "${WORKDIR}/${acc}_1P.fq.gz" "${WORKDIR}/${acc}_1U.fq.gz" \
      "${WORKDIR}/${acc}_2P.fq.gz" "${WORKDIR}/${acc}_2U.fq.gz" \
      ILLUMINACLIP:"${ADAPT}":2:30:10 \
      MINLEN:36
    cp -f "${WORKDIR}/${acc}_1P.fq.gz" "${WORKDIR}/${acc}_1U.fq.gz" \
      "${WORKDIR}/${acc}_2P.fq.gz" "${WORKDIR}/${acc}_2U.fq.gz" \
      "$OUT_DIR/"
    rm -f "${WORKDIR}/${acc}_1P.fq.gz" "${WORKDIR}/${acc}_1U.fq.gz" \
      "${WORKDIR}/${acc}_2P.fq.gz" "${WORKDIR}/${acc}_2U.fq.gz"
  else
    trimmomatic PE -threads 4 -phred33 \
      "${RAW_DIR}/${acc}_1.fastq.gz" "${RAW_DIR}/${acc}_2.fastq.gz" \
      "${OUT_DIR}/${acc}_1P.fq.gz" "${OUT_DIR}/${acc}_1U.fq.gz" \
      "${OUT_DIR}/${acc}_2P.fq.gz" "${OUT_DIR}/${acc}_2U.fq.gz" \
      ILLUMINACLIP:"${ADAPT}":2:30:10 \
      MINLEN:36
  fi
done

if [[ -n "$WORKDIR" ]]; then
  rmdir "$WORKDIR" 2>/dev/null || true
fi
