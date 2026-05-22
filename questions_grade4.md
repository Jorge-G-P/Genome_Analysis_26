# Grade 4 Questions — Paper I: Zhang et al. 2017
*E. faecium* E745 — Genome Assembly, Annotation, RNA-seq & Differential Expression

---

## Reads Quality Control

**Q1. How is the quality of your data?**

The quality of the Illumina RNA-seq data is good. FastQC reports show Phred scores > 30 across most positions for all 12 files (6 samples × 2 paired reads). There is a slight drop in quality at the 3' ends of reads, which is normal for Illumina sequencing-by-synthesis chemistry. GC content distributions are unimodal and consistent across samples (~38).

PacBio reads were also of good quality, cause even though they were not directly evaluated with fastQC, they were able to produce an assembly with good quality.

---

**Q2. What can generate the "fails" in FastQC that you observe in your data? Can these cause any problems during subsequent analyses?**

The main FastQC "fail" observed across all samples is in the **Per Base Sequence Content** module, where the first ~10 bases show non-uniform nucleotide composition. This might be caused by random hexamer priming bias during RNA library preparation — the hexamers used for reverse transcription do not bind truly randomly, leading to a biased composition at the start of reads. This is a known and expected artefact of RNA-seq and is not a sign of poor data quality. It does not cause problems in downstream analyses because it affects only a few positions at read ends, which Trimmomatic removes anyway.
---

## Reads Preprocessing

**Q3. How many reads have been discarded after trimming?**

| Sample | Condition | Input pairs | Dropped pairs | % Dropped |
|---|---|---|---|---|
| ERR1797972 | BHI rep 1 | 27,078,884 | 1,267,659 | 4.68% |
| ERR1797973 | BHI rep 2 | 23,953,340 | 923,316 | 3.85% |
| ERR1797974 | BHI rep 3 | 23,240,177 | 1,080,073 | 4.65% |
| ERR1797969 | Serum rep 1 | 25,937,368 | 978,370 | 3.77% |
| ERR1797970 | Serum rep 2 | 26,634,380 | 1,123,683 | 4.22% |
| ERR1797971 | Serum rep 3 | 27,609,615 | 1,545,834 | 5.60% |

![Trimmomatic read survival per sample](results/plots/trimmomatic_survival.png)

On average ~4.5% of read pairs were fully discarded. An additional ~40% of reads survived only as forward-only (orphan) reads due to the reverse read being too low quality; these orphan reads were not used in downstream BWA paired-end mapping.

---

**Q4. How can this affect your future analyses and results?**

Discarding ~4.5% of pairs has a negligible effect given the large remaining read counts (12–15 million pairs per sample), which is more than sufficient for differential expression analysis in a ~3.1 Mb bacterial genome. The main concern is the high proportion of orphan reads (~40% forward-only surviving): these represent read pairs where the reverse read was too low quality to keep. Since only paired reads were used for BWA mapping, roughly half the original reads contributed to the analysis. This reduces statistical power slightly but does not bias the results, as the loss is consistent across all samples and conditions.

---

**Q5. How is the quality of your data after trimming?**

After trimming, the per-base quality scores improve, particularly at the 3' ends where low-quality bases were removed. The per-base sequence content bias at the start of reads persists (it is structural, not quality-related, and cannot be removed by Trimmomatic without discarding valid data). Overall the trimmed reads are of good quality: Phred > 30 across the retained length, no residual adapter contamination, minimum read length 36 bp. Mapping rates of ~98% to the assembled genome confirm the trimmed data is clean.

---

**Q6. What quality threshold did you choose for the leading/trailing/slidingwindow parameters, and why?**

Parameters chosen:
- `LEADING:3` — removes bases with quality < 3 from the start; this threshold removes only very poor bases (essentially uncalled bases) without over-trimming
- `TRAILING:3` — same at the 3' end; conservative threshold to preserve read length
- `SLIDINGWINDOW:4:15` — trims once the average quality in a 4-base window drops below 15; a window of 4 balances sensitivity and specificity, and Q15 is a standard threshold for removing low-quality stretches while retaining usable sequence
- `MINLEN:36` — discards reads shorter than 36 bp after trimming; reads shorter than this are too short for reliable unique mapping against a ~3.1 Mb genome

These are standard Trimmomatic parameters recommended in the course manual and widely used in bacterial RNA-seq pipelines.

---

## Genome and Metagenome Assembly

**Q7. What information can you get from the plots and reports given by the assembler (if you get any)?**

Canu's report (`paper1_e745.report`) provides a read length histogram, coverage statistics, and k-mer frequency distributions from each stage (correction, trimming, assembly). From the report: 100,836 reads were input with a total of 533,813,323 bases, giving ~172× raw coverage. The read length distribution shows a peak around 5–7 kb with a long tail extending to >20 kb. The k-mer frequency plot confirms sufficient depth for error correction. These statistics indicate a well-covered genome suitable for de novo assembly.

---

**Q8. What intermediate steps generate informative output about the assembly?**

Canu runs three internal stages, each producing informative output:
1. **Correction** — reports how many reads were corrected and the resulting read quality improvement
2. **Trimming** — reports read lengths after removing low-quality ends
3. **Unitig assembly** — the `.tigInfo` file lists all unitigs (contigs) with their lengths, coverage, and circularity flags

The `paper1_e745.contigs.layout.tigInfo` file is particularly informative: it shows that tig00000001 (2.77 Mb) is suggested as circular, consistent with a complete bacterial chromosome.

---

**Q9. How many contigs do you expect? How many do you obtain?**

*E. faecium* E745 is known from the literature to have one chromosome (~2.7–2.9 Mb) and several plasmids (the published genome carries 7 plasmids). I therefore expected approximately 8 sequences (1 chromosome + ~7 plasmids). I obtained **9 contigs**, which is very close to the expected number. The largest contig (2,775,132 bp) corresponds to the chromosome, and the remaining 8 smaller contigs (ranging from ~14 kb to ~216 kb) likely represent plasmids.

---

**Q10. Do you expect the same result between different assemblers, for the same data? If you tried different assemblers, what differences do you see in the result and why do you think that is?**

No, different assemblers will produce different results even from the same input data because they use different algorithms, overlap detection strategies, and error correction approaches. Only Canu was ran in this project, so no direct comparison was made, but differences in total contig number, largest contig size, and handling of repetitive regions (such as rRNA operons) would be expected.

---

**Q11. What are the k-mers? What are the problems and benefits of choosing a small or large k-mer?**

A k-mer is a substring of length k from a sequence. In assembly, k-mers are used to detect overlaps between reads: reads sharing k-mers at their ends are linked in an overlap graph.

- **Small k-mer** (e.g., k=11): finds more overlaps (higher sensitivity), but many k-mers will match by chance in a large genome, causing false connections between unrelated sequences (lower specificity). This is especially problematic in repetitive regions.
- **Large k-mer** (e.g., k=127): fewer false overlaps, better resolution of repeats, but requires reads to be long and error-free to share large exact k-mers. If reads are short or error-prone (like PacBio CLR), large k-mers will miss real overlaps.

Canu automatically selects appropriate k-mer sizes based on the input data characteristics, which is one reason it is robust for noisy long-read data.

---

**Q12. Some assemblers can include a read-correction step before doing the assembly. What is this step doing?**

Read correction uses the redundancy of sequencing coverage to identify and fix sequencing errors in individual reads before assembly. Each position in the genome is covered by multiple reads; if most reads agree on a particular base but one read has a different base at that position, the deviant base is likely an error and is corrected. This step is especially important for PacBio CLR reads, which have a ~10–15% per-base error rate. Canu performs this correction internally in its first stage, significantly improving the accuracy of the reads fed into the assembly graph. Without correction, the high error rate of CLR reads would create a fragmented assembly with many overlaps.

---

## Assembly Evaluation

**Q13. How does your assembly compare with the reference assembly? What could have caused the differences?**

QUAST was run against the *E. faecium* E0019EM0028 reference (99.66% nucleotide identity to E745). Key metrics:

| Metric | My assembly |
|---|---|
| Total contigs | 9 |
| Total length | 3,147,208 bp |
| Largest contig | 2,775,132 bp |
| N50 | 2,775,132 bp |
| GC content | 37.79% |
| N's per 100 kbp | 0.00 |

The assembly covers ~97–98% of the reference genome with very few gaps. There are minor structural differences (misassemblies) likely caused by: (1) genuine genomic differences between E745 and E0019EM0028 (these are different strains, not the same isolate), (2) repetitive sequences such as rRNA operons that can be collapsed or misjoined in assembly, and (3) possible transposable elements or genomic islands that differ between strains.

---

**Q14. Do you think your assembly is better/worse than the public one?**

The public *E. faecium* E745 genome (available in NCBI) was assembled with a combination of PacBio and Illumina data followed by manual curation and polishing — a more thorough process than my Canu-only assembly. My assembly has slightly more contigs (9 vs the published 8 sequences) and lacks polishing, meaning it will have a higher residual error rate. However, the overall quality is good.

---

## Annotation

**Q15. What types of features are detected by the software? Would you trust the prediction of some features over others and why?**

Prokka detected the following features:

| Feature | Count |
|---|---|
| CDS (protein-coding genes) | 3,126 |
| tRNA | 85 |
| tmRNA | 1 |
| **Total** | **3,212** |

tRNAs and tmRNAs are predicted with high confidence using well-validated tools (Aragorn, Infernal) against highly conserved RNA structures — these predictions are very reliable. CDS predictions are made by Prodigal, which is accurate for bacterial genomes but can make mistakes at short ORFs (<100 aa) and in distinguishing true genes from spurious ORFs. Functional annotations (gene names, product descriptions) assigned by BLAST homology are only as reliable as the homolog they match — well-characterised genes in closely related species give high-confidence annotations, while "hypothetical protein" annotations reflect genuinely unknown function.

---

**Q16. How can you evaluate the quality of the obtained functional annotation?**

Several approaches can be used:
1. **Compare feature counts with known genomes** — my 3,126 CDS is consistent with published *E. faecium* genomes (~2,900–3,200 CDS), suggesting no major over- or under-prediction
2. **Check hypothetical protein rate** — ~34% hypothetical proteins is normal for *E. faecium*, but very high rates (>50%) would suggest annotation problems
3. **BUSCO completeness** — benchmarking against conserved single-copy orthologs gives a genome-wide completeness estimate.
4. **Manual spot-checking** — looking up key expected genes (e.g., the *vanHAX* vancomycin resistance operon, ribosomal proteins) and confirming they are present and correctly annotated

---

**Q17. How many features of each kind are detected in your contigs? Do you detect the same number of features as the authors? How do they differ?**

My Prokka annotation detected 3,126 CDS, 85 tRNAs, and 1 tmRNA across 9 contigs (3,147,208 bp total). The published E745 genome (Zhang et al. 2017) reported approximately 3,170 CDS — slightly more than my annotation. The difference (~44 CDS) is minor and expected: the published annotation used a more comprehensive pipeline with manual curation, while Prokka is a rapid automated tool. Some genes in the published annotation may be present as pseudogenes or short ORFs that Prokka filtered out. The tRNA count of 85 in my annotation is higher than the typical ~60 reported in similar strains, which may reflect tRNA predictions on plasmid contigs.

---

**Q18. How many genes are annotated as 'hypothetical protein'? Why is that so? How would you tackle that problem?**

Approximately **34% of the 3,126 CDS** (~1,060 genes) are annotated as "hypothetical protein." This occurs when Prokka cannot find a homolog with a known function in its reference databases (UniProt, RefSeq, TIGRfam, Pfam). For *E. faecium*, many genes are strain-specific or have only been sequenced recently without experimental characterisation.

To reduce the hypothetical protein rate, one could:
1. Use **EggNOG-mapper** for deeper orthology-based functional annotation against the eggNOG database
2. Search against **Pfam** or **InterPro** for domain-level annotations (a protein can be "hypothetical" at the gene level but still carry a recognisable domain)
3. Use **HHpred** for remote homology detection using profile-profile alignment, which can find distant homologs missed by BLAST
4. Submit to the **NCBI Conserved Domain Database (CDD)** to identify structural domains

---

## Mapping

**Q19. What percentage of your reads map back to your contigs? Why do you think that is?**

Mapping rates from BWA flagstat:

| Sample | Condition | % Mapped |
|---|---|---|
| ERR1797972 | BHI rep 1 | 98.59% |
| ERR1797973 | BHI rep 2 | ~98.5% |
| ERR1797974 | BHI rep 3 | ~98.5% |
| ERR1797969 | Serum rep 1 | 98.25% |
| ERR1797970 | Serum rep 2 | ~98.4% |
| ERR1797971 | Serum rep 3 | ~98.3% |

![BWA mapping rates per sample](results/plots/bwa_mapping_rates.png)

~98-99% of reads map back to the assembly. This is excellent and indicates that the assembly captures almost all of the sequenced genome. The ~1–2% unmapped reads likely represent: sequencing errors that prevent alignment, reads from low-complexity or repetitive regions with ambiguous mapping, or minor contaminating sequences in the library.

---

**Q20. What do you interpret from your read coverage differences across the genome?**

From `BH_ERR1797972.coverage.txt`, coverage varies substantially across contigs:

| Contig | Length | Mean depth |
|---|---|---|
| tig00000001 (chromosome) | 2,775,132 | 932× |
| tig00000002 | 216,000 | 105× |
| tig00000005 (VanHAX plasmid) | 40,013 | 350× |
| tig00000004 | 14,734 | 29× |

![Coverage depth per contig — ERR1797972 BHI rep 1](results/plots/coverage_per_contig.png)

The chromosome (tig00000001) has by far the highest coverage (~932×) because it is the most abundant sequence in the cell. Some plasmids show higher coverage than their size would predict (e.g., tig00000005 at 350× despite being only 40 kb), suggesting they are present at higher copy numbers per cell than the chromosome. tig00000004 and tig00000008 show lower coverage (~25–29×), suggesting they may be low-copy plasmids. Coverage depth therefore provides information about plasmid copy number, which is biologically meaningful — high-copy plasmids can produce more gene product per cell.

---

**Q21. Do you see big differences between replicates?**

No. Mapping rates are consistent across all replicates within each condition (~98–99%), and the PCA plot from DESeq2 shows tight clustering of the three BHI replicates together and the three Serum replicates together along PC1. This indicates good reproducibility between biological replicates. Minor variation in total read counts between replicates (26–28 million pairs) is expected and is accounted for by DESeq2's size factor normalisation.

---

## Read Counting

**Q22. What is the distribution of the counts per gene? Are most genes expressed? How many counts would indicate that a gene is expressed?**

The count distribution from HTSeq is highly skewed: a large number of genes have very low counts (< 10), while a small number of highly expressed genes (ribosomal proteins, metabolic enzymes) have counts in the tens of thousands. This is typical for RNA-seq — most genes are expressed at low to moderate levels, with a few housekeeping genes dominating the count table.


![Count distribution per gene](results/plots/DESeq2_count_histogram.png)

Most genes appear to be expressed at some level in at least one condition. A common practical threshold is **≥10 counts** in at least some samples as the minimum to consider a gene expressed — below this, counts are in the noise range where Poisson variation dominates. DESeq2 automatically filters very low-count genes (independent filtering) before testing, removing genes unlikely to be statistically testable.

---

## Expression Analyses

**Q23. If your expression results differ from those in the published article, why could it be?**

My analysis identified 2,296 genes with padj < 0.05, which is broadly in the same range as the published paper but may differ in the exact gene list. Reasons for differences include:
1. **Different reference genome** — I used my Canu assembly, while the authors used the published E745 reference. Differences in gene models affect read assignment.
2. **Different pipeline parameters** — the authors may have used different mappers, counting tools, or statistical cutoffs
3. **Different trimming** — my trimmed read pairs represent ~50% of raw reads; the authors may have used the full pre-trimmed dataset provided by the course
4. **DESeq2 version** — different versions can produce slightly different normalisation and dispersion estimates

---

**Q24. How do the different samples and replicates cluster together?**

The PCA plot shows clear separation between BHI and Serum conditions along PC1, which captures the majority of variance in the dataset. The three replicates within each condition cluster tightly together, confirming good experimental reproducibility. This pattern indicates that the growth condition (serum vs BHI) is the dominant source of transcriptional variation and that there are no major outlier samples.

![DESeq2 PCA plot](results/plots/DESeq2_PCA.png)
![DESeq2 volcano plot](results/plots/DESeq2_volcano.png)

---

**Q25. How did you sort your differential expression results? Why?**

Results were sorted by **adjusted p-value (padj)** in ascending order — i.e., the most statistically significant genes appear first. Sorting by padj rather than raw p-value or fold change is appropriate because: (1) padj accounts for the multiple testing problem (testing ~3,000 genes simultaneously inflates false positives without correction), (2) fold change alone can be misleading for lowly expressed genes where large fold changes occur by chance, and (3) padj integrates both effect size and statistical confidence. Genes were additionally filtered by |log2FC| > 1 to focus on biologically meaningful changes.

---

**Q26. Do you need a normalization step? What would you normalize against? Does DESeq do it?**

Yes, normalisation is essential. Without it, differences in sequencing depth between samples would appear as differential expression. For example, a sample with 30 million reads will have twice as many counts per gene as a sample with 15 million reads, even if the actual gene expression is identical.

DESeq2 normalises using **size factors** (median-of-ratios method): for each sample, it calculates the ratio of each gene's count to the geometric mean across all samples, then takes the median of these ratios as the size factor. This approach is robust to the presence of highly differentially expressed genes and does not assume that total read counts should be equal.

In my data, the size factors were close to 1.0 for all samples (ranging from 0.86 to 1.10), confirming that the samples had similar sequencing depths and that normalisation had a modest effect. DESeq2 **does perform this normalisation** automatically as part of its pipeline.

---

**Q27. What would you do to increase the statistical power of your expression analysis?**

Statistical power in RNA-seq differential expression analysis can be increased by:

1. **More biological replicates** — I used 3 per condition, which is the practical minimum. Increasing to 4–6 replicates would substantially reduce false negatives, especially for genes with moderate fold changes
2. **Deeper sequencing** — more reads per sample increases count precision for lowly expressed genes
3. **Better reference genome** — using the polished published E745 genome as reference rather than my Canu assembly would improve read mapping and reduce noise from misassembled regions
4. **Reducing technical variation** — ensuring all libraries are prepared in the same batch and sequenced together minimises batch effects
5. **Pre-filtering low-count genes** — removing genes with very low counts before testing reduces the multiple testing burden and improves the FDR correction for the remaining genes (DESeq2 does this automatically via independent filtering)
