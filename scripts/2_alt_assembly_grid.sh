#!/bin/bash
#SBATCH --job-name=canu_e745
#SBATCH --account=uppmax2026-1-94
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
# Log paths are relative to the submit cwd unless absolute; pin both cwd and logs to REPO_ROOT.
#SBATCH --chdir=/home/gorgis/projects/genome/Genome_Analysis_26
#SBATCH --output=/home/gorgis/projects/genome/Genome_Analysis_26/outputs/canu_driver-%j.out
#SBATCH --error=/home/gorgis/projects/genome/Genome_Analysis_26/errors/canu_driver-%j.err

# set -e exit on first failing command; set -u error on unset variables; pipefail treat a pipeline as failed if any stage fails.
set -euo pipefail

REPO_ROOT=/home/gorgis/projects/genome/Genome_Analysis_26
cd "$REPO_ROOT" || exit 1

# module load bioinfo-tools 2>/dev/null || true
module load canu
module load SAMtools/1.22.1-GCC-13.3.0

# turns on Bash’s nullglob option, so if there's no match we dont pass the literal string name "*.fastaqc"
shopt -s nullglob
READS=(data/1_Zhang_2017/genomics_data/PacBio/*.fastq.gz)
if [ "${#READS[@]}" -eq 0 ]; then
  echo "No PacBio FASTQ.gz files under data/1_Zhang_2017/genomics_data/PacBio/"
  exit 1
fi

# Zhang 2017 E745: bacterial genome ~3 Mbp (adjust if your manual specifies otherwise).
# Input names contain .subreads -> CLR-style reads, not HiFi.
canu -p paper1_e745 -d data/tmp/2_canu_assembly \
  genomeSize=3.1m \
  useGrid=true \
  gridEngine=slurm \
  gridOptions="--account=uppmax2026-1-94 --time=03:00:00" \
  -pacbio-raw "${READS[@]}"



# genomeSize --- note:
# If you set it too small (e.g. true ~3 Mb, you say 1m)
# Canu overestimates read depth (same reads divided by a smaller genome).
# Effects can include: stricter or wrong decisions in correction/overlap steps, more aggressive filtering, higher resource requests for some stages, or weirder stopping conditions. In bad cases you can get poor correction, incomplete assembly, or failure if the pipeline thinks coverage is absurdly high or inconsistent.
# Usually you notice bad stats (claimed coverage off the charts) or quality/contiguity problems, not a literal “3 Mb genome assembled as 1 Mb.”
# If you set it too large (e.g. true ~3 Mb, you say 10m)
# Canu underestimates coverage.
# Risk: steps that need minimum coverage may under-trigger, merge less, or stop early with “not enough data,” or you get a fragmented / low-quality assembly because correction/overlap logic assumed you had much more data than you do.
# If it is moderately wrong (e.g. 2.8m vs 3.2m for a ~3 Mb genome)
# Often little changes in the final contigs, because true depth is still in a sane band. This is the usual case for bacteria when you’re within ~20–30%.