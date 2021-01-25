# Longitudinal omics in Syrian hamsters integrated with human data unravel cellular effector responses to moderate COVID-19
Single cell sequencing is a powerful tool to investigate cellular mechanisms of disease pathogenesis. The golden Hamster (Mesocricetus auratus) is a valuable model for SARS-CoV-2 infection as molecular mechanisms seem comparable to humans.

Abstract:

In COVID-19, immune responses are key in determining disease severity. However, cellular mechanisms at the onset of inflammatory lung injury in SARS-CoV-2 infection, particularly involving endothelial cells, remain ill-defined. Using Syrian hamsters as model for moderate COVID-19, we conducted a detailed longitudinal analysis of systemic and pulmonary cellular responses, and corroborated it with datasets from COVID-19 patients. Monocyte-derived macrophages in lungs exerted the earliest and strongest transcriptional response to infection, including induction of pro-inflammatory genes, while epithelial cells showed weak alterations. Without evidence for productive infection, endothelial cells reacted, depending on cell subtypes, by strong and early expression of anti-viral, pro-inflammatory, and T cell recruiting genes. Recruitment of cytotoxic T cells as well as emergence of IgM antibodies preceded viral clearance at day 5 post infection. Investigating SARS-CoV-2 infected Syrian hamsters can thus identify cell type-specific effector functions, provide detailed insights into pathomechanisms of COVID-19, and inform therapeutic strategies.

![Infection scheme]
(https://github.com/Berlin-Hamster-Single-Cell-Consortium/Single-cell-sequencing-of-COVID-19-pathogenesis-in-golden-Hamsters-/blob/main/scheme.jpg)

Infection Scheme Figure was created with Biorender.com

# table of content / data processing

This repository contains the code that was used create to analyze the data an create the figures in the manuscript available on bioRxiv (https://www.biorxiv.org/content/10.1101/2020.12.18.423524v1). 

# Getting started / installation

For running cell ranger, please refer to https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references for requirements.

We have used R version 3.6 to run the code. All required libraries are listed at the beginning of the individual R files and are freely available through Bionconductor (https://www.bioconductor.org) or for direct installation using install.packages(). Installing library and dependencies typically requires about one hour. For a quick start, we recommend starting from the provided Seurat objects provided on http://www.mdc-berlin.de/singlecell-SARSCoV2, which contain all cell type annotations and embeddings that were used in the manuscript mentioned above.

The Seurat object is several gigabytes in size, we therefore recommend to use a computer with a least 16 GB memory.

# Cell ranger
Raw fastq files can be downloaded from GEO, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162208
To create a cell ranger reference, see the description of the mkref command here: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references

We have used a combined hamster/virus fasta/gtf file. For creating and downloading the gtf file, see the https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162208 (data processing entries describe how the downloadable file ma1_genes_longer.gtf was generated).

# Data analysis in R
Processing cell ranger output: hamster_merging.R contains the code to create the various Seurat objects. For further processing of the blood sample object (seu_blood_combined_mtfilt_blood500.rds) refer to blood_preprocessing_clustering.R. For the lung object (seu_lung_combined_mtfilt.rds), proceed with the integration code in lung_hamster_scRNAseq_integrate.R.

Annotation of lung cell types from tabula muris / Travaglini et al.: lung_hamster_annotation.Rmd

Processing and figures for Fig. 1 and Fig. S3: hamster_scRNAseq.R contains the code to annote the cell types and create the panels in Fig. 1 and Fig. S3. The Seurat object ma_int.rds can be downloaded via http://www.mdc-berlin.de/singlecell-SARSCoV2. There is also a downsampled data set available (ma_int_red.rds, 3000 cells) for an initial look at the data, however this is not suitable for any statistical analysis due to the small number of cells. 

Pseudobulk analysis: hamster_scRNAseq_pseudobulkDE.R contains the code for the differential expression dotblots in Fig. 3, 4, 5, S5, S7. The Seurat object ma_int.rds as well as the data table pseudobulk.txt and the enrichment tables (KEGG and process) can be downloaded via http://www.mdc-berlin.de/singlecell-SARSCoV2

The combined blood/lung analysis for 2 dpi (Fig. 4C, Fig. S5C) is detailed in lung_hamster_scRNAseq_bloodcomparison_2dpi.R

