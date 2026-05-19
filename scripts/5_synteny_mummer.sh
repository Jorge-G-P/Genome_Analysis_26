#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 1
#SBATCH -t 00:15:00
#SBATCH -J mummerplot
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

module load MUMmer

# fasta obtained from ncbi website (ran homology there)
nucmer E.faecium-strainVRE095.fasta E.faecium.contigs.fasta

mummerplot --png --layout -p mummerplot out.delta