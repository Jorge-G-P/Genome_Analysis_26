# Grade 5 Questions — Paper I: Zhang et al. 2017
Extra analyses: ResFinder (AMR) + Tn-seq (conditional essentiality)

---

## Extra Analysis 6 — Evaluating Antibiotic Resistance Potential (ResFinder)

**Q6.1. How many antibiotics is the strain predicted to be resistant to? And sensitive?**

E. faecium E745 is predicted to be **resistant to 6 antibiotics**: vancomycin, teicoplanin, ampicillin, ciprofloxacin, gentamicin, and erythromycin.

It is predicted to have **no resistance to 6 antibiotics**: fosfomycin, quinupristin+dalfopristin, tigecycline, tetracycline, linezolid, and chloramphenicol.

---

**Q6.2. Which mutation in the gene gyrA confers resistance to nalidixic acid and ciprofloxacin? And in the gene parC?**

**gyrA:** mutation **p.E87G** (GAG → GGG, Glu87Gly at codon 87) confers resistance to nalidixic acid and ciprofloxacin.

**parC:** mutation **p.S80I** (AGC → ATC, Ser80Ile at codon 80) confers resistance to nalidixic acid and ciprofloxacin.

Both GyrA and ParC are targets of fluoroquinolone antibiotics. These substitutions alter the drug-binding pocket of each enzyme, reducing fluoroquinolone binding affinity. The presence of mutations in both genes confers high-level quinolone resistance.

---

**Q6.3. Which mapping method uses ResFinder?**

ResFinder uses **BLAST** when an assembled genome in FASTA format is submitted. When raw reads (FASTQ) are submitted instead, it uses **KMA** (K-Mer Alignment), a faster k-mer based aligner suited for unassembled data. In this analysis, BLAST was used since the Canu assembly was submitted.

---

**Q6.4. Why is it important to evaluate the antibiotic resistance potential of a bacterial strain?**

It is key to know how to defend patients from bacterial infections. In this specific case:

1. E745 is **resistant to vancomycin** — the standard last-resort antibiotic for Gram-positive infections. Knowing this is essential to select an effective alternative.

2. **Infection control:** VRE (*vancomycin-resistant Enterococcus*) spreads easily in hospital settings. Characterising resistance profiles informs isolation and decontamination protocols.

3. **Epidemiology:** Identifying specific resistance mechanisms enables tracking of resistance spread between strains and hospitals.

4. Understanding the full resistance landscape of clinical strains in order to guide rational antibiotic use policies, avoiding risks.

---

## Extra Analysis 7 — Identify Essential Genes for Growth in Human Serum (Tn-seq)

**Q7.1. How did the authors get the data from the Tn-seq analysis? What type of data is it?**

The authors created a saturating transposon insertion library of *E. faecium* E745 by randomly inserting a mariner transposon throughout the genome. Each insertion disrupts the gene it lands in. This library was then grown competitively in two conditions: heat-inactivated human serum and BHI medium. After growth, the bacteria were harvested and the genomic DNA was sequenced (Illumina, single-end 50 nt reads) to quantify how many reads map to each transposon insertion site per gene. Samples used: ERR1801009/010/011 (HI Serum) and ERR1801012/013/014 (BHI), 3 biological replicates each.

It is **short-read Illumina sequencing data** (single-end, 50 nt) targeting transposon-flanking sequences to measure insertion abundance per gene across the genome. Each read has the structure: 6 nt barcode + genomic flanking sequence + transposon sequence. Before mapping, the transposon sequence must be removed to recover the 16 nt genomic insert that identifies the insertion site.

---

**Q7.2. What is the goal of the Tn-seq analysis?**

The goal is to identify genes that are **conditionally essential** for growth specifically in human serum. If a gene is required for survival in serum, mutants with a transposon disrupting that gene will be selectively lost from the serum-grown population but will remain present in the BHI-grown population (where the gene is not needed). By comparing the abundance of transposon insertions per gene between serum and BHI using DESeq2, genes that are depleted in serum are identified as conditionally essential for serum growth. These represent potential virulence and survival factors specific to the bloodstream environment.

---

**Q7.3. Which genes seem to be important for E. faecium to grow in human serum? Attach a plot that supports your conclusion, analyze it and explain briefly your workflow.**

DESeq2 (contrast: HI_Serum vs BHI, reference = BHI) identified **15 genes depleted in serum** (padj < 0.05, log2FC < -1), meaning these genes are conditionally essential for growth in human serum:

| Gene ID | Gene | log2FC | padj | Function |
|---|---|---|---|---|
| LPCHMCBP_01329 | *sorA_1* | -15.04 | 2.3e-31 | PTS sorbose-specific EIIC (carbohydrate uptake) |
| LPCHMCBP_00663 | *purA* | -14.98 | 4.2e-31 | Adenylosuccinate synthetase (purine biosynthesis) |
| LPCHMCBP_02968 | — | -11.14 | 1.4e-07 | Putative PTS IIB component |
| LPCHMCBP_01716 | *rpoN1* | -10.89 | 5.6e-15 | RNA polymerase sigma-54 factor |
| LPCHMCBP_00498 | *msr(C)* | -10.73 | 2.2e-04 | ABC-F ribosomal protection protein |
| LPCHMCBP_02408 | *cmpD* | -10.05 | 9.4e-04 | Bicarbonate transport ATP-binding protein |
| LPCHMCBP_02294 | *yidA_3* | -9.99 | 1.3e-03 | Sugar phosphatase |
| LPCHMCBP_02587 | — | -10.43 | 4.0e-02 | Hypothetical protein |
| LPCHMCBP_01116 | *yhdG* | -9.95 | 2.2e-02 | Putative amino acid permease |
| LPCHMCBP_02043 | — | -9.81 | 5.9e-03 | Hypothetical protein |
| LPCHMCBP_02424 | *pyrC* | -9.52 | 1.4e-02 | Dihydroorotase (pyrimidine biosynthesis) |
| LPCHMCBP_01326 | *dgaR_2* | -5.51 | 3.9e-219 | Transcriptional regulator of sugar uptake |
| LPCHMCBP_01330 | *manZ_3* | -5.21 | 3.8e-29 | PTS mannose-specific EIID (carbohydrate uptake) |
| LPCHMCBP_00636 | *guaB* | -5.73 | 1.6e-03 | IMP dehydrogenase (purine biosynthesis) |
| LPCHMCBP_02996 | — | -1.55 | 1.8e-02 | Hypothetical protein |

An additional **9 genes were enriched in serum** (log2FC > 1, padj < 0.05), meaning their disruption confers a fitness advantage in serum:

| Gene ID | Gene | log2FC | padj | Function |
|---|---|---|---|---|
| LPCHMCBP_02560 | *ytfQ* | +10.41 | 2.7e-05 | ABC transporter periplasmic-binding protein |
| LPCHMCBP_02100 | *arcA* | +9.51 | 2.3e-05 | Arginine deiminase |
| LPCHMCBP_00030 | *pyk* | +9.24 | 2.8e-02 | Pyruvate kinase |
| LPCHMCBP_00070 | — | +3.79 | 4.3e-04 | Alpha-monoglucosyldiacylglycerol synthase |
| LPCHMCBP_01846 | *dacA* | +3.22 | 6.7e-73 | D-alanyl-D-alanine carboxypeptidase |
| LPCHMCBP_02588 | *iolU_2* | +3.18 | 1.6e-17 | scyllo-inositol 2-dehydrogenase |
| LPCHMCBP_01763 | *clsA_1* | +2.96 | 2.5e-03 | Major cardiolipin synthase |
| LPCHMCBP_00437 | — | +2.53 | 3.6e-20 | Hypothetical protein |
| LPCHMCBP_02061 | — | +1.40 | 1.7e-02 | Hypothetical protein |

**Workflow summary:**
1. Remove transposon sequence from raw 50 nt reads with cutadapt, retaining 16 nt genomic insert (`scripts/extra_tnseq/11_cut.sh`)
2. Map 16 nt trimmed reads to Canu assembly with Bowtie2 (`scripts/extra_tnseq/13_tnseq_bowtie.sh`)
3. Count transposon insertions per gene with HTSeq-count using Prokka GFF (`scripts/extra_tnseq/13_tnseq_htseq.sh`)
4. Differential abundance analysis with DESeq2, contrast HI_Serum vs BHI (`scripts/extra_tnseq/14_tnseq_deseq.R`)
5. Visualisation with volcano plot, PCA, and count histogram (`scripts/extra_tnseq/15_tnseq_plot.py`)

The supporting plot is the **volcano plot** below. Each point is a gene. Points to the left (negative log2FC) with high -log10(padj) represent genes depleted in serum — conditionally essential genes. Points to the right represent genes enriched in serum (disruption is advantageous). The two most significant hits (*purA* and *sorA_1*) are the leftmost points.

![Tn-seq volcano plot — HI_Serum vs BHI](../results/plots/TnSeq_volcano.png)
![Tn-seq PCA plot](../results/plots/TnSeq_PCA.png)

The PCA plot shows separation between HI Serum and BHI replicates, confirming that the growth condition drives the main source of variation in insertion abundance.

The two dominant biological themes among conditionally essential genes are:

1. **Nucleotide biosynthesis** — *purA* (adenylosuccinate synthetase), *guaB* (IMP dehydrogenase), and *pyrC* (dihydroorotase) are required for de novo synthesis of purines and pyrimidines. Human serum contains very low concentrations of free nucleotides, forcing *E. faecium* to synthesise them from scratch. This is fully consistent with the paper's key finding: *purA*, *guaB*, *purD*, *purH*, *pyrF*, and *pyrK_2* were identified as essential in serum by Zhang et al. (2017).

2. **Carbohydrate uptake (PTS transporters)** — *manZ_3*, *sorA_1*, *dgaR_2*, and a putative PTS IIB component are all depleted in serum. Glucose is the only carbohydrate freely available in blood; efficient glucose import via phosphotransferase systems is therefore critical. The paper specifically identified *manZ_3* and *manY_2* as among the most essential PTS genes for serum growth.

The enriched genes (whose disruption is advantageous in serum) include cell wall and membrane remodelling genes (*dacA*, *clsA_1*), consistent with the paper's observation that non-essential cell wall biosynthesis imposes a metabolic cost in the nutrient-limited serum environment.
