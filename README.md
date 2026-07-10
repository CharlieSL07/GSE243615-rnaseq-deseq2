# RNA-seq Differential Expression Analysis of Single and Combination Drug Treatment (GSE243615)
Differential gene expression analysis of MZ1, PDD and GSK single drug treatments and MZ1+PDD and MZ1+GSK combination treatments on HCT116 cells using DESeq2 and R.

## Introduction
This project is an independent reanalysis of RNA-seq data originally 
collected and published by Mori et al. (2024):

> Mori Y, Akizuki Y, Honda R, Takao M, et al. & Ohtake F (2024).
> Intrinsic signaling pathways modulate targeted protein degradation.
> *Nature Communications*. https://doi.org/10.1038/s41467-024-49519-z

Raw count data were obtained from the NCBI Gene Expression Omnibus 
(GEO accession: [GSE243615](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243615)).

The original authors examined the transcriptional effects of the PROTAC BET
degrader MZ1 in combination with the PARG (poly(ADP-ribose) glycohydrolase)
inhibitor PDD00017273 and PERK (Protein kinase R-like endoplasmic reticulum
kinase) inhibitor GSK2606414 in HCT116 (RRID:CVCL_0291), an adherent human
colorectal carcinoma cell line widely used in cancer biology, particularly
for studying cell cycle regulation, DNA damage response, and drug
sensitivity. Here I independently reproduce their differential expression
analysis using DESeq2 and extend it with GO pathway enrichment analysis using
clusterProfiler.

The primary study found that PDD alone did not induce significant alterations
in gene expression but enhanced the transcriptional alterations induced by
MZ1. Unlike PDD, GSK altered a distinct subset of genes irrespective of MZ1
co-treatment, while also enhancing MZ1's effects.

MZ1 is a PROTAC-based BET protein degrader, this molecule binds to an E3
ligase and a target protein in the BET protein family such as BRD4, BRD3 and
BRD2. Interaction between the ligase and target protein result in
the ubiquitination of the target protein and subsequent degradation via the
proteasome. Specifically, BET proteins regulate gene transcription and
chromatin architecture. MZ1 most potently degrades BRD4, while BRD3 and BRD2
degradation is less efficient and may be enhanced by co-treatment with
inhibitors such as PDD and GSK (Mori et al., 2024). PDD and GSK exploit
distinct cellular stress pathways that may converge with BRD4-dependent
transcription. I hypothesised that co-treatment with either compound would
amplify MZ1-induced transcriptional changes, particularly in pathways 
associated with chromatin remodelling, MAPK/ERK signalling, and cell
migration, producing a broader disruption of gene expression than MZ1 alone.

## Methods
Raw counts were downloaded from GEO (GSE243615) and imported into R. One GSK
replicate was excluded from the counts matrix by the original authors,
presumably due to QC failure. GSK group therefore has 2 replicates whilst
other conditions have 3. Differential expression analysis was performed using
DESeq2 with a six-level factor design formula comprising vehicle (control), 
PDD, MZ1, GSK, MZ1+PDD, and MZ1+GSK treatment groups. Log2 fold changes were
shrunk using the ashr method. Genes were called significant at padj < 0.05
and |log2FC| > 1. GO enrichment was performed using clusterProfiler with
Entrez gene IDs against the org.Hs.eg.db annotation database, with
Benjamini-Hochberg (BH) multiple testing correction at a significance
threshold of padj < 0.05, in order to control the false discovery rate.

## Results

### Quality Control plots:
#### Principal Component Analysis plot:
![PCA Plot](outputs/pca_plot.png)
#### Sample-sample distances heat map:
![sample distances heat map](outputs/sample_distances.png)

### Differential Expression plots:
#### Summary:
Genes called significant at padj < 0.05 and |log2FC| > 1

|                   |comparison         | total_genes|   up| down| total_sig|
|:------------------|:------------------|-----------:|----:|----:|---------:|
|PDD vs vehicle     |PDD vs vehicle     |       15045|    0|    1|         1|
|MZ1 vs vehicle     |MZ1 vs vehicle     |       21317|  803| 1594|      2397|
|GSK vs vehicle     |GSK vs vehicle     |       17733|  212|  190|       402|
|MZ1+PDD vs vehicle |MZ1+PDD vs vehicle |       22661| 2004| 2712|      4716|
|MZ1+GSK vs vehicle |MZ1+GSK vs vehicle |       22661| 2085| 2708|      4793|
|MZ1+PDD vs MZ1     |MZ1+PDD vs MZ1     |       20421|  673|  494|      1167|
|MZ1+PDD vs PDD     |MZ1+PDD vs PDD     |       22213| 1919| 2491|      4410|
|MZ1+GSK vs MZ1     |MZ1+GSK vs MZ1     |       20869|  853|  642|      1495|
|MZ1+GSK vs GSK     |MZ1+GSK vs GSK     |       22213| 1846| 2463|      4309|

#### Volcano plots:
![outputs/volcano_plots.pdf](outputs/volcano_plots.pdf)
#### Heatmaps:
![outputs/heatmaps.pdf](outputs/heatmaps.pdf)
#### GO Enrichment dot plots:
![outputs/go_enrichment_dotplots.pdf](outputs/go_enrichment_dotplots.pdf)

### Key findings:

## Discussion


## References
