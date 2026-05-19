#!/bin/bash
#SBATCH --job-name=htseq
#SBATCH --account=uppmax2026-1-61
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/htseq-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/htseq-%j.err

set -x
set -euo pipefail

REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
cd "$REPO_ROOT" || exit 1

module load HTSeq/2.1.2-gfbf-2024a

BAM_DIR=data/tmp/7_bwa_rna
OUT_DIR=data/tmp/8_htseq_counts
mkdir -p "$OUT_DIR"

GFF_SRC=data/tmp/4_prokka_ann/paper1.gff
GFF="${OUT_DIR}/paper1_no_fast.gff"

# Student Manual: .gff must not contain nucleotide sequences (Prokka adds ##FASTA).
awk '/^##FASTA/{exit} {print}' "$GFF_SRC" > "$GFF"

# Manual: BAM sorted by position. -s no: HTSeq default is stranded; use yes/reverse if your library is stranded.
# Prokka: count CDS by ID. -m union is HTSeq default; kept explicit.
shopt -s nullglob
BAMS=( "${BAM_DIR}"/*.sorted.bam )
if [ "${#BAMS[@]}" -eq 0 ]; then
  echo "No *.sorted.bam under ${BAM_DIR}/"
  exit 1
fi

for bam in "${BAMS[@]}"; do
  base=$(basename "$bam" .sorted.bam)
  htseq-count -r pos -s no -m union -t CDS -i ID "$bam" "$GFF" > "${OUT_DIR}/${base}.counts.tsv"
done
