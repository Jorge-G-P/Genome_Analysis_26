#!/bin/bash
#SBATCH --job-name=tnseq_cut
#SBATCH --account=uppmax2026-1-94
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/tnseq_cut-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/tnseq_cut-%j.err

module load cutadapt/5.0

# INPUT
INPUT=/proj/uppmax2026-1-61/Genome_Analysis/1_Zhang_2017/transcriptomics_data/Tn-Seq_BHI/trim_ERR1801014_pass.fastq.gz

# CONSTANT SEQUENCES
TRANSPOSON="AACAGGTTGGATGATAAGTCCCCGGTCTTC"

#OUTPUT
OUT1=./step1_notransposon.fastq
OUTFINAL=./control_14_no_transposon.fastq

echo "=== Step 1: Trim 3' magellan6 transposon end ==="
cutadapt -a $TRANSPOSON -O 10 -o $OUT1 $INPUT

echo "=== Step 2: Remove first 6 bases (barcode) and keep 16 nt ==="
cutadapt -u 6 --length 16 --minimum-length 16 -o $OUTFINAL $OUT1