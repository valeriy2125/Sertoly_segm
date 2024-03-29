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
6) also we use Seurat(v.4.1.1) for scRNAseq data
R script is in the ...... files.

![image](https://github.com/valeriy2125/Sertoly_segm/assets/101557211/e64f54dd-8d81-4f20-81c4-fc1c0c95bc8f)

Fig. 1. Scheme of isolation of cultures of SLCs and SCs

# Results
Data analysis showed, that the transcriptomes of SLCs and SCs cultures are different.
862 of 20711 genes are differentially expressed. 549 are overexpressed, 313 are underexpressed (Fig. 2)

![image](https://github.com/valeriy2125/Sertoly_segm/assets/101557211/bf27c54b-1444-4381-bad1-5dbc0d53013c)

Fig. 2. Volcano plot differentially expressed genes. Gray - all genes, blue - Padj less than 0.05, LFC greater than 1 modulo, red - average TPM greater than 1.

The top 10 overexpressed and underexpressed genes in SLC and SC cultures are presented (Fig. 3). Pax8 and Sox17 genes are selected as potential markers of SLCs for further IGH studies. GO analysis showed differences between SLCs and SCs cultures in extensive term such as cytosol, membrane, and extracellular matrix. It is worth noting that an enrichment of the term "male gonad development" was shown for SLCs culture, which may indicate incomplete differentiation of SLCs (Fig. 3).

![image](https://github.com/valeriy2125/Sertoly_segm/assets/101557211/b9de36a5-983b-4ffe-b012-0ccd1a6fb2d2)

Fig. 3. The top 10 overexpressed and underexpressed genes for SLCs and SCs (ges potentially suitable for further IGH are marked in red). The top 5 enriched terms for SLCs and SCs

The SLCs and SCs clusters were identified from open scRNAseq data. The population of SLCs is rather small, which may indirectly confirm the hypothesis that there is a niche of Sertoli stem cells in the local region of the transient zone. 

![image](https://github.com/valeriy2125/Sertoly_segm/assets/101557211/6bb9cc38-f794-499e-bab9-e7fc21f00a08)

Fig. 4. UMAP based on the testis scRNAseq data. The arrows indicate the SLC and SC cluster.

Analysis of differentially expressed genes confirms data from bulk RNA seq data.

![image](https://github.com/valeriy2125/Sertoly_segm/assets/101557211/dd834563-e7e5-4aeb-9537-721b892dbdd2)

Fig. 5. The top 10 overexpressed and underexpressed genes for SLCs and SCs from scRNAseq data (genes correlated with bulk RNA seq data are marked in red)

# Conclusions

Identified SLC marker genes – Pax8, Sox17

