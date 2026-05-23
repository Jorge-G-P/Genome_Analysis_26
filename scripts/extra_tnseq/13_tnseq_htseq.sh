#!/bin/bash
#SBATCH --job-name=tnseq_htseq
#SBATCH --account=uppmax2026-1-61
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/tnseq_htseq-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/tnseq_htseq-%j.err

set -euo pipefail

REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
cd "$REPO_ROOT" || exit 1

module load HTSeq/2.1.2-gfbf-2024a

BAM_DIR=data/tmp/12_tnseq_bowtie
OUT_DIR=data/tmp/13_tnseq_htseq
mkdir -p "$OUT_DIR"

GFF_SRC=data/tmp/4_prokka_ann/paper1.gff
GFF="${OUT_DIR}/paper1_no_fasta.gff"

# Strip the ##FASTA section that HTSeq cannot parse.
[ -f "$GFF" ] || awk '/^##FASTA/{exit} {print}' "$GFF_SRC" > "$GFF"

# Confirm the GFF actually contains CDS features; fall back to 'gene' if not.
FEATURE_TYPE=CDS
if ! grep -q $'\t'"CDS"$'\t' "$GFF"; then
  echo "WARNING: No CDS features found in GFF — falling back to feature type 'gene'"
  FEATURE_TYPE=gene
fi

shopt -s nullglob
BAMS=( "${BAM_DIR}"/*.sorted.bam )
if [ "${#BAMS[@]}" -eq 0 ]; then
  echo "No *.sorted.bam under ${BAM_DIR}/"
  echo "Run scripts/12_tnseq_bwa.sh first."
  exit 1
fi

# -r pos     : BAMs are position-sorted
# -s no      : unstranded — correct for Tn-seq; insertion orientation is not meaningful here
# -m union   : count a read if it overlaps any part of a feature
# -t CDS/-i ID : Prokka annotates coding sequences as CDS with ID= locus-tag attributes
# --additional-attr gene : attach gene names to output for easier DESeq2 interpretation
for bam in "${BAMS[@]}"; do
  base=$(basename "$bam" .sorted.bam)
  htseq-count \
    -r pos \
    -s no \
    -m union \
    -t "$FEATURE_TYPE" \
    -i ID \
    "$bam" "$GFF" \
    > "${OUT_DIR}/${base}.counts.tsv"
  echo "Counted ${base}"
done

echo "Tn-seq counts written to ${OUT_DIR}"
