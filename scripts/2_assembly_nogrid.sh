#!/bin/bash
#SBATCH --job-name=canu_e745_1n
#SBATCH --account=uppmax2026-1-61
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/canu_nogrid-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/canu_nogrid-%j.err

set -euo pipefail

REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
cd "$REPO_ROOT" || exit 1

module load canu
module load SAMtools/1.22.1-GCC-13.3.0

shopt -s nullglob
READS=(data/1_Zhang_2017/genomics_data/PacBio/*.fastq.gz)
if [ "${#READS[@]}" -eq 0 ]; then
  echo "No PacBio FASTQ.gz files under data/1_Zhang_2017/genomics_data/PacBio/"
  exit 1
fi

# useGrid=false: all stages run in this allocation only (no nested sbatch). Tune #SBATCH --time above.
# avoid child job creation that was giving trouble with course quota.
THREADS="${SLURM_CPUS_PER_TASK:-8}"
canu -p paper1_e745 -d data/tmp/2_canu_assembly_nogrid \
  genomeSize=3.1m \
  useGrid=false \
  "maxThreads=${THREADS}" \
  -pacbio-raw "${READS[@]}"