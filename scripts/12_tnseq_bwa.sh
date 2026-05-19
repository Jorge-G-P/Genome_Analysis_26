#!/bin/bash
#SBATCH --job-name=tnseq_bwa
#SBATCH --account=uppmax2026-1-61
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/tnseq_bwa-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/tnseq_bwa-%j.err

set -euo pipefail

REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
cd "$REPO_ROOT" || exit 1

module load BWA/0.7.19-GCCcore-13.3.0
module load SAMtools/1.22.1-GCC-13.3.0

THREADS=${SLURM_CPUS_PER_TASK:-4}

REF=data/tmp/2_canu_assembly/paper1_e745.contigs.fasta
OUT_DIR=data/tmp/12_tnseq_bwa
mkdir -p "$OUT_DIR"

[ -f "${REF}.bwt" ] || bwa index "$REF"

# Tn-seq: single-end reads mapped to assembled genome.
# Conditionally essential screen: HI Serum vs BHI (3 replicates each).
# Adjust TNSEQ_DIR and sample IDs if Tn-seq data is at a different path.

TNSEQ_DIR=data/1_Zhang_2017/transcriptomics_data/Tn-seq

# Prefer trimmed reads; fall back to raw.
if [ -d "${TNSEQ_DIR}/trimmed" ]; then
  FASTQ_DIR="${TNSEQ_DIR}/trimmed"
elif [ -d "${TNSEQ_DIR}/raw" ]; then
  FASTQ_DIR="${TNSEQ_DIR}/raw"
else
  FASTQ_DIR="${TNSEQ_DIR}"
fi

echo "Tn-seq FASTQ directory: ${FASTQ_DIR}"

# HI Serum samples (heat-inactivated serum, 3 replicates).
HI_SERUM_LIST=${HI_SERUM_LIST:-"ERR1801009 ERR1801010 ERR1801011"}
# BHI samples (3 replicates used in conditional screen).
BHI_LIST=${BHI_LIST:-"ERR1801012 ERR1801013 ERR1801014"}

map_tnseq() {
  local acc=$1 label=$2
  # Use -print -quit so find exits after the first match; avoids SIGPIPE from head.
  local fq
  fq=$(find "$FASTQ_DIR" -name "${acc}*.fastq.gz" -print -quit 2>/dev/null)
  # Also accept .fq.gz if .fastq.gz not found.
  if [ -z "$fq" ]; then
    fq=$(find "$FASTQ_DIR" -name "${acc}*.fq.gz" -print -quit 2>/dev/null)
  fi
  if [ -z "$fq" ]; then
    echo "WARNING: No FASTQ for ${acc} in ${FASTQ_DIR}/ — skipping"
    return
  fi
  # Include @RG header so BAMs carry sample identity for downstream tools.
  bwa mem -t "$THREADS" \
    -R "@RG\tID:${acc}\tSM:${label}\tLB:${label}\tPL:ILLUMINA" \
    "$REF" "$fq" \
    | samtools sort -@ "$THREADS" -o "${OUT_DIR}/${label}_${acc}.sorted.bam" -
  samtools index "${OUT_DIR}/${label}_${acc}.sorted.bam"
  echo "Mapped ${acc} → ${label}"
}

for acc in $HI_SERUM_LIST; do map_tnseq "$acc" "HI_Serum"; done
for acc in $BHI_LIST;       do map_tnseq "$acc" "BHI";        done

echo "Tn-seq BWA mapping complete. BAM files in ${OUT_DIR}"
