#!/bin/bash
#SBATCH --job-name=tnseq_cut
#SBATCH --account=uppmax2026-1-61
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/tnseq_cut-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/tnseq_cut-%j.err

set -euo pipefail

module load cutadapt/5.0

# Transposon sequence to trim from 3' end (Magellan6 mariner transposon)
TRANSPOSON="AACAGGTTGGATGATAAGTCCCCGGTCTTC"

# Input directories (pre-trimmed Tn-seq FASTQ files)
HSERUM_DIR=/proj/uppmax2026-1-61/Genome_Analysis/1_Zhang_2017/transcriptomics_data/Tn-Seq_HSerum
BHI_DIR=/proj/uppmax2026-1-61/Genome_Analysis/1_Zhang_2017/transcriptomics_data/Tn-Seq_BHI

# Output directory
OUT_DIR=data/tmp/11_tnseq_cut
mkdir -p "$OUT_DIR"

# Process a single sample:
#   1. Remove 3' transposon sequence
#   2. Remove 6 nt barcode from 5' end and keep exactly 16 nt genomic insert
cut_sample() {
    local acc=$1
    local indir=$2
    local input="${indir}/trim_${acc}_pass.fastq.gz"
    local tmp="${OUT_DIR}/${acc}_notransposon.fastq.gz"
    local final="${OUT_DIR}/${acc}_16nt.fastq.gz"

    echo "=== Processing ${acc} ==="
    echo "  Input: ${input}"

    echo "  Step 1: Trim 3' transposon sequence"
    cutadapt -a "$TRANSPOSON" -O 10 \
        --cores 2 \
        -o "$tmp" "$input"

    echo "  Step 2: Remove 6 nt barcode, keep 16 nt genomic insert"
    cutadapt -u 6 --length 16 --minimum-length 16 \
        --cores 2 \
        -o "$final" "$tmp"

    rm "$tmp"
    echo "  Done: ${final}"
}

# HI Serum replicates
cut_sample ERR1801009 "$HSERUM_DIR"
cut_sample ERR1801010 "$HSERUM_DIR"
cut_sample ERR1801011 "$HSERUM_DIR"

# BHI replicates
cut_sample ERR1801012 "$BHI_DIR"
cut_sample ERR1801013 "$BHI_DIR"
cut_sample ERR1801014 "$BHI_DIR"

echo "=== All 6 Tn-seq samples trimmed. Output in ${OUT_DIR} ==="
