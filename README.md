# Genome Analysis 26 — Paper I - Jorge

Reproduction of the Zhang et al. 2017 analysis of *Enterococcus faecium* E745:
genome assembly, annotation, RNA-seq differential expression, AMR profiling, and Tn-seq.

## Overview

| Step | Tool | Output |
|---|---|---|
| PacBio assembly | Canu | 3.15 Mb, 9 contigs, N50 2.77 Mb |
| Assembly QC | QUAST | `results/3.5_quast/` |
| Annotation | Prokka | 3,126 CDS, 85 tRNA, 44% hypothetical |
| Synteny | MUMmer | `results/5_synteny/` |
| RNA-seq trim | Trimmomatic | `results/2_trimmomatic/` |
| RNA-seq mapping | BWA mem | ~98.3–98.6% mapped |
| Read counts | HTSeq-count | `results/7_HTSeq/` |
| Diff. expression | DESeq2 | 1,252 DEGs (padj<0.05, \|LFC\|>1) |
| AMR | ResFinder (web) | `results/E2_resfinder/` |
| Tn-seq | BWA + HTSeq + DESeq2 | `results/E1_TnSeq_DESeq/` |

## Repository layout

```
scripts/          # SLURM batch scripts (.sh) and analysis code (.py, .R)
scripts/extra_tnseq/  # Tn-seq specific scripts
results/          # All interesting output files (gitignored raw data in data/)
docs/             # wiki.md, analysis plan, questions for grade 4 and 5
```

## Data

Raw data is under `data/1_Zhang_2017/` (gitignored). 

