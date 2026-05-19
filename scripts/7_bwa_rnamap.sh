#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH --account=uppmax2026-1-61
#SBATCH --time=04:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/bwa-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/bwa-%j.err

set -x
set -euo pipefail

REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
cd "$REPO_ROOT" || exit 1

module load BWA/0.7.19-GCCcore-13.3.0
module load SAMtools/1.22.1-GCC-13.3.0

REF=data/tmp/2_canu_assembly/paper1_e745.contigs.fasta
OUT_DIR=data/tmp/7_bwa_rna
mkdir -p "$OUT_DIR"

[ -f "${REF}.bwt" ] || bwa index "$REF"

# Paired libraries only (trim_paired_*); trim_single_* are Trimmomatic orphan reads — do not mix into PE bwa mem.

# BH condition — three independent replicates.
BH_TRIM=data/1_Zhang_2017/transcriptomics_data/RNA-Seq_BH/trimmed
for acc in ERR1797972 ERR1797973 ERR1797974; do
  bwa mem -t 4 "$REF" \
    "${BH_TRIM}/trim_paired_${acc}_pass_1.fastq.gz" \
    "${BH_TRIM}/trim_paired_${acc}_pass_2.fastq.gz" \
  | samtools sort -@ 4 -o "${OUT_DIR}/BH_${acc}.sorted.bam" -
  samtools index "${OUT_DIR}/BH_${acc}.sorted.bam"
done

# Serum condition — three independent replicates.
SERUM_TRIM=data/1_Zhang_2017/transcriptomics_data/RNA-Seq_Serum/trimmed
for acc in ERR1797969 ERR1797970 ERR1797971; do
  bwa mem -t 4 "$REF" \
    "${SERUM_TRIM}/trim_paired_${acc}_pass_1.fastq.gz" \
    "${SERUM_TRIM}/trim_paired_${acc}_pass_2.fastq.gz" \
  | samtools sort -@ 4 -o "${OUT_DIR}/Serum_${acc}.sorted.bam" -
  samtools index "${OUT_DIR}/Serum_${acc}.sorted.bam"
done
