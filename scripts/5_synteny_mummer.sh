#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 1
#SBATCH -t 00:15:00
#SBATCH -J mummerplot
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out
module load MUMmer
REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
cd "$REPO_ROOT" || exit 1

# Reference obtained from ncbi database. First blast was run, then the closest sequence was downloaded as 
# a reference. (It was communicated that the reference in uppmax was not valid). 
REF=data/E0019_full.fasta
QUERY=data/tmp/2_canu_assembly/paper1_e745.contigs.fasta
OUT_DIR=data/tmp/5_synteny_mummer
mkdir -p "$OUT_DIR"
nucmer -p "${OUT_DIR}/mummerplot" "$REF" "$QUERY"
cd "$OUT_DIR"
mummerplot --png --layout -p mummerplot mummerplot.delta
echo "MUMmerplot done. Output in ${OUT_DIR}/"