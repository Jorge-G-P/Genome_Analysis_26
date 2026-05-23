#!/bin/bash
# Quick check: verify all 6 trimmed Tn-seq files have 16 nt reads
# Run on UPPMAX after 11_cut.sh finishes:
#   bash scripts/extra_tnseq/11_check_lengths.sh

TRIM_DIR=data/tmp/11_tnseq_cut

for fq in "$TRIM_DIR"/*_16nt.fastq.gz; do
    sample=$(basename "$fq" _16nt.fastq.gz)
    # Read lengths: sample every 4th line (sequence lines), get unique lengths
    lengths=$(zcat "$fq" | awk 'NR%4==2 {print length($0)}' | sort -u | tr '\n' ' ')
    echo "${sample}: read lengths = ${lengths}"
done
