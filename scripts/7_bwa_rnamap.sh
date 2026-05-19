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
# Added for Q19 and Q20. Coverage plot.

STATS_DIR=${OUT_DIR}/stats
mkdir -p "$STATS_DIR"

# refactor into a function, it was redundant.
map_and_qc() {
  local label=$1 acc=$2 fq1=$3 fq2=$4
  local bam="${OUT_DIR}/${label}_${acc}.sorted.bam"

  bwa mem -t 4 "$REF" "$fq1" "$fq2" \
    | samtools sort -@ 4 -o "$bam" -
  samtools index "$bam"

  # Q19 — % reads mapped
  samtools flagstat "$bam" > "${STATS_DIR}/${label}_${acc}.flagstat.txt"

  # Q20 — per-contig coverage depth
  samtools coverage "$bam" > "${STATS_DIR}/${label}_${acc}.coverage.txt"
}

# BH condition — three independent replicates.
BH_TRIM=data/1_Zhang_2017/transcriptomics_data/RNA-Seq_BH/trimmed
for acc in ERR1797972 ERR1797973 ERR1797974; do
  map_and_qc BH "$acc" \
    "${BH_TRIM}/trim_paired_${acc}_pass_1.fastq.gz" \
    "${BH_TRIM}/trim_paired_${acc}_pass_2.fastq.gz"
done

# Serum condition — three independent replicates.
SERUM_TRIM=data/1_Zhang_2017/transcriptomics_data/RNA-Seq_Serum/trimmed
for acc in ERR1797969 ERR1797970 ERR1797971; do
  map_and_qc Serum "$acc" \
    "${SERUM_TRIM}/trim_paired_${acc}_pass_1.fastq.gz" \
    "${SERUM_TRIM}/trim_paired_${acc}_pass_2.fastq.gz"
done

echo "Flagstat and coverage reports in ${STATS_DIR}"
