# RNA-seq Differrential Expression Analysis of Single and Combination Drug Treatment (GSE243615)
Differential gene expression analysis of MZ1, PDD and GSK single drug treatments and MZ1+PDD and MZ1+GSK combination treatments on HCT116 cells using DESeq2 and R.

## Introduction

## Methods

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


## Discussion
