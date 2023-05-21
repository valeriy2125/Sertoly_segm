# Search for markers of the Sertoli-like cell population we discovered
# Introduction
The main cells ensuring the course of spermatogenesis in vivo in mammals are Sertoli cells (SC). They are located in the spermatogenic epithelium of the seminal tubules of the testis and provide development and protection of germ cells. Such an important role makes SCs a promising target for the development of in vitro spermatogenesis systems or male infertility treatment techniques, but the proliferation of this cell type in culture is very difficult.

Recently, our laboratory discovered a new cell population in the testis network area, which we named Sertoli-like cells (SLC). They possess similar markers to SC genes, but active proliferation in culture. In addition, under 3D culture conditions with neonatal SCs, SLCs are able to support the development of germ cells up to the initial stages of meiosis. Thus, SLCs could be an alternative to neonatal SCs in in vitro spermatogenesis technology, but more accurate assessment of their transcriptome and identification of marker genes is required.
# Aim and tasks

Aim:
Define markers of the Sertoli-like cells

Objectives:
1) Search for differentially expressed SC and SLC genes
2) Perform a GO enrichment analysis
# Methods

We used assembled transcriptomes from culture of SLC and SC (three replicates for each point; Mus musculus)
SLC and SC cultures were obtained from different parts of the testis (Fig. 1)
RNA was isolated from each culture and sequenced. As a result, we got illumina single-end reads for each culture

1) FastQC (v0.11.9)
2) Samtools (v1.17)
3) HISAT2 (v2.2.1) alignment
4) featureCounts (v2.0.1) (getting counts table)
5) DESeq2 (v1.26.0) (gene expression analysis in Rstudio)
6) also we use Seurat3 for scRNAseq dataa






![image](https://github.com/valeriy2125/Sertoly_segm/assets/101557211/e64f54dd-8d81-4f20-81c4-fc1c0c95bc8f)

Fig. 1. Scheme of isolation of cultures of SLCs and SCs
