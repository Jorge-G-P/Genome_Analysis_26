# Genome Analysis 26 — Paper I - Jorge

Reproduction of the Zhang et al. 2017 analysis of *Enterococcus faecium* E745:
genome assembly, annotation, RNA-seq differential expression, AMR profiling, and Tn-seq.


## Repository layout

```
scripts/          # SLURM batch scripts (.sh) and analysis code (.py, .R)
scripts/extra_tnseq/  # Tn-seq specific scripts
results/          # All interesting output files (gitignored raw data in data/)
docs/             # wiki.md, analysis plan, questions for grade 4 and 5
errors/           # slurm error logs, kept for log history
outputs/          # slurm output logs, kept for log history
```


## Steps Overview

| Step | Tool |
|---|---|---|
| PacBio assembly | Canu |
| Assembly QC | QUAST |
| Annotation | Prokka |
| Synteny | MUMmer |
| RNA-seq trim | Trimmomatic |
| RNA-seq mapping | BWA mem |
| Read counts | HTSeq-count |
| Diff. expression | DESeq2 |
| AMR (extra) | ResFinder (web) |
| Tn-seq (extra) | BWA + HTSeq + DESeq2 |


## Data

Raw data is under `data/1_Zhang_2017/` (gitignored). 

