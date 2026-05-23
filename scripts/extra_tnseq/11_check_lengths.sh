#!/bin/bash
#SBATCH --job-name=tnseq_cut
#SBATCH --account=uppmax2026-1-61
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/tnseq_cut-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/tnseq_cut-%j.err

TRIM_DIR=data/tmp/11_tnseq_cut

for fq in "$TRIM_DIR"/*_16nt.fastq.gz; do
    sample=$(basename "$fq" _16nt.fastq.gz)
    # Read lengths: sample every 4th line (sequence lines), get unique lengths
    lengths=$(zcat "$fq" | awk 'NR%4==2 {print length($0)}' | sort -u | tr '\n' ' ')
    echo "${sample}: read lengths = ${lengths}"
done
