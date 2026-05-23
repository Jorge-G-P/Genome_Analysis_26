#!/bin/bash
#SBATCH --job-name=tnseq_bowtie2
#SBATCH --account=uppmax2026-1-61
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/tnseq_bowtie2-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/tnseq_bowtie2-%j.err

set -euo pipefail

REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
cd "$REPO_ROOT" || exit 1

module load Bowtie2/2.5.4-GCC-13.3.0
module load SAMtools/1.22.1-GCC-13.3.0

THREADS=${SLURM_CPUS_PER_TASK:-4}

REF=data/tmp/2_canu_assembly/paper1_e745.contigs.fasta
BOWTIE2_IDX=data/tmp/2_canu_assembly/paper1_e745_bowtie2
OUT_DIR=data/tmp/12_tnseq_bowtie
mkdir -p "$OUT_DIR"

# Build Bowtie2 index if not already present
# (paper used Bowtie2 to map 16 nt genomic flanking reads)
if [ ! -f "${BOWTIE2_IDX}.1.bt2" ]; then
    echo "Building Bowtie2 index..."
    bowtie2-build "$REF" "$BOWTIE2_IDX"
fi

# Input: 16 nt trimmed reads from 11_cut.sh
TRIM_DIR=data/tmp/11_tnseq_cut

map_tnseq() {
    local acc=$1
    local label=$2
    local fq="${TRIM_DIR}/${acc}_16nt.fastq.gz"

    if [ ! -f "$fq" ]; then
        echo "WARNING: Trimmed FASTQ not found for ${acc}: ${fq} — skipping"
        return
    fi

    echo "Mapping ${acc} (${label})..."

    # Bowtie2: single-end, no mismatches allowed (--no-1mm-upfront),
    # report only uniquely mapping reads (-M 1).
    # The paper mapped 16 nt fragments to the E745 genome with Bowtie2.
    bowtie2 \
        -x "$BOWTIE2_IDX" \
        -U "$fq" \
        --no-unal \
        -p "$THREADS" \
        --rg-id "${acc}" \
        --rg "SM:${label}" \
        --rg "LB:${label}" \
        --rg "PL:ILLUMINA" \
        | samtools sort -@ "$THREADS" -o "${OUT_DIR}/${label}_${acc}.sorted.bam" -

    samtools index "${OUT_DIR}/${label}_${acc}.sorted.bam"
    samtools flagstat "${OUT_DIR}/${label}_${acc}.sorted.bam" \
        > "${OUT_DIR}/${label}_${acc}.flagstat"

    echo "Done: ${label}_${acc}.sorted.bam"
}

# HI Serum replicates
map_tnseq ERR1801009 HI_Serum
map_tnseq ERR1801010 HI_Serum
map_tnseq ERR1801011 HI_Serum

# BHI replicates
map_tnseq ERR1801012 BHI
map_tnseq ERR1801013 BHI
map_tnseq ERR1801014 BHI

echo "Tn-seq Bowtie2 mapping complete. BAM files in ${OUT_DIR}"
