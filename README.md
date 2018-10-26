# DENDRO & DENDROplan

Single cell RNA-sequencing (scRNA-seq) technologies have made it possible to study tumor transcriptomes at single-cell resolution, allowing better cell type classification and unprecedented study of intratumor transcriptomic heterogeneity.  Yet, the underlining mechanisms for tumor initiation, progression, metastasis and relapse are often driven by genetic heterogeneity, that is, subclonal evolution involving genetically distinct subpopulations of cells.  Detection of these genetically distinct subclones and understanding the evolutionary dynamics of the tumor is essential for accurate prognosis and effective treatment.  A tumor’s genetic heterogeneity underlies its transcriptome heterogeneity, and thus, intra-tumor transcriptomic variation needs to be quantified within the context of intra-tumor clonal variation. There is yet no method that can derive the genetic mutation profile and gene expression level for the same cell with acceptable accuracy, high throughput, and low cost. Here we describe a statistical method, **DENDRO** (**D**na based **E**volutio**N**ary tree pre**D**iction by sc**R**na-seq techn**O**logy), which takes scRNA-seq data for a tumor and accurately reconstructs its phylogeny, assigning each single cell from the scRNA-seq data to a subclone (Figure 1).

Before conducting a single cell RNA-seq experiment on a tumor sample, it is important to project how subclone detection power depends on the number of cells sequenced and the coverage per cell. To facilitate experiment design, we developed a tool, **DENDROplan** (Figure 2), that  predicts the expected clustering accuracy by DENDRO given sequencing parameters.  Given a tree structure and a target accuracy, **DENDROplan** computes the necessary read depth and number of cells needed based on the spike-in procedure described above. 


## Manuscript

XXX ([link](https://doi.org/10.1093/bioinformatics/bty057))


## Questions & Problems

If you have any questions or problems when using DENDRO or DENDROplan, please feel free to open a new issue [here](https://github.com/zhouzilu/DENDRO/issues). You can also email the maintainers of the corresponding packages -- the contact information is shown under Developers & Maintainers.


## Installation

Install to R/RStudio
Install all packages in the latest version of [R](https://www.r-project.org/).
```r
devtools::install_github("zhouzilu/DENDRO")
```

## Pipeline overview

This DENDRO package includes two type of analysis: (1) **DENDRO** analysis on real dataset such as tumor and hematopoesis, and (2) **DENDROplan** which analyze the accuracy of DENDRO cluster given inferred clonal tree structure, cell number and sequencing depth. Overall pipelines are shown below.

### DENDRO pipeline

<p align="center">
  <img src='https://raw.githubusercontent.com/zhouzilu/DENDRO/tree/master/figure/Pkg_FIG-01.jpg' width='500' height='300'>
  </p>

  **Figure 1.** A flowchart outlining the procedures of DENDRO. DENDRO starts from scRNA-seq raw data. We recommend STAR 2-pass method for mapping because it is more robust with splicing junction. SNA detection was applied to mapped BAM files. Both counts of total allele reads and counts of alternative allele reads for each cell c at mutation position g are collected. In the next step, a cell-to-cell genetic divergence matrix is calculated using a genetic divergence evaluation function. DENDRO further clusters the cells and polls cells from same cluster together and re-estimate SNA profiles. Based on the re-estimated SNA profiles, DENDRO generates a parsimony tree which shows the evolution relationship between subclones.

### DENDROplan pipeline

<p align="center">
  <img src='https://raw.githubusercontent.com/zhouzilu/DENDRO/tree/master/figure/Pkg_FIG-02.jpg' width='500' height='300'>
  </p>

  **Figure 2.** The overall DENDROplan pipeline. The analysis starts with a designed tree with an interested clade (purple clade in the example). Based on the tree model, number of cells, sequencing depth and sequencing error rate, we simulate single cell mutation profile. scRNA data was sampled from a reference scRNA-seq dataset given expression level in bulk. A phylogeny computed by DENDRO is further compared with underlining truth, which measured by three statistics - adjust Rand index (global accuracy statistics), capture rate (subclone specific statistic) and purity (subclone specific statistic). 


## Running DENDRO

  **R notebook** with step-by-step demonstration and rich display is available [***here***](https://raw.githubusercontent.com/zhouzilu/DENDRO/tree/master/vignette/DENDRO_vignette.html). Corresponding **Rmd script** is available [***here***](https://github.com/zhouzilu/DENDRO/blob/master/vignette/DENDRO_vignette.Rmd).

## Running DENDROplan

  **R notebook** with step-by-step demonstration and rich display is available [***here***](https://raw.githubusercontent.com/zhouzilu/DENDRO/tree/master/vignette/DENDROplan_vignette.html). Corresponding **Rmd script** is available [***here***](https://github.com/zhouzilu/DENDRO/blob/master/vignette/DENDROplan_vignette.Rmd).


## Citation

Please cite DENDRO.

* **DENDRO**: [no link yet](https://doi.org/10.1093/bioinformatics/bty057)
<br>
  Genetic Heterogeneity Profiling by Single Cell RNA Sequencing ([GitHub](https://github.com/zhouzilu/DENDRO))

## Developers & Maintainers

* [Zilu Zhou](https://statistics.wharton.upenn.edu/profile/zhouzilu/) (zhouzilu at pennmedicine dot upenn dot edu)
  <br>
  Genomics and Computational Biology Graduate Group, University of Pennsylvania

* [Nancy R. Zhang](https://statistics.wharton.upenn.edu/profile/nzh/) (nzh at wharton dot upenn dot edu)
  <br>
  Department of Statistics, University of Pennsylvania
