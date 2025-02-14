# Comparison of data quality from 10X Ge-nomics-based single-nucleus ATAC-seq technologies

This report aims to compare two approaches of single nucleus assay for transposase accessible chro-matin using sequencing (snATAC-seq) by 10x Genomics: Monomodal snATAC-seq and multimodal single nucleus multiome ATAC + Gene expression sequencing. These methods provide insight into the structure of chromatin in single cells, but the multimodal method gives an additional layer of in-formation by simultaneously performing RNA sequencing on the same cell. However, multimodal snATAC-seq is both double the cost of monomodal snATAC-seq, and more time consuming.  For this reason, when deciding which of the two methods to use, one needs to consider if this extra layer of information is worth the added cost, and if the quality of the snATAC-seq data is compromised for more information. In this report, the quality of monomodal and multimodal snATAC-seq data, is com-pared, using different quality control metrics. These include: Fragment count, Transcriptional star site (TSS) enrichment scores, blacklist ratio, frequency of reads in peaks (FRiP), and proportion of poor quality nuclei. Based on this comparison, it appears that monomodal snATAC-seq produces data of slightly better quality, compared to multimodal snATAC-seq, however, overall, both methods produce high quality snATAC-seq data. Possible errors and reasons for inconsistencies in the samples are dis-cussed in the report.

`ISA_final_script.R`: this file contains all code used for the analysis 

![illustration of work](https://github.com/frejaaugustesen/ISA_2024/blob/main/illustration/ISA_workflow.png)


