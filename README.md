# Longitudinal omics in Syrian hamsters integrated with human data unravel complexity of moderate immune responses to SARS-CoV-
Single cell sequencing is a powerful tool to investigate cellular mechanisms of disease pathogenesis. The golden Hamster (Mesocricetus auratus) is a valuable model for SARS-CoV-2 infection as molecular mechanisms seem comparable to humans.

Abstract:

In COVID-19, the immune response largely determines disease severity and is key to therapeutic strategies. Cellular mechanisms contributing to inflammatory lung injury and tissue repair in SARS-CoV-2 infection, particularly endothelial cell involvement, remain ill-defined. We performed detailed spatiotemporal analyses of cellular and molecular processes in SARS-CoV-2 infected Syrian hamsters. Comparison of 3 hamster single-cell sequencing and proteomics with data sets from COVID-19 patients demonstrated inter-species concordance of cellular and molecular host-pathogen interactions. In depth vascular and pulmonary compartment analyses (i) supported the hypothesis that monocyte-derived macrophages dominate inflammation, (ii) revealed endothelial inflammation status and T-cell attraction, and (iii) showed that CD4+ and CD8+ cytotoxic T-cell responses precede viral elimination. Using the Syrian hamster model of self-limited moderate COVID-19, we defined the specific roles of endothelial and epithelial cells, among other myeloid and non-myeloid lung cell subtypes, for determining the disease course.

Keywords: coronavirus disease 2019, COVID-19; Syrian hamster; SARS-CoV2, inflammatory monocytes, macrophages, cytotoxic T cells, endothelial cells 

Infection Scheme Figure was created with Biorender.com

# table of content / data processing

Running cell ranger: Raw fastq files can be downloaded from GEO, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162208
To create a cell ranger reference, see the description of the mkref command here: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references

We have used a combined hamster/virus fasta/gtf file. For creating and downloading the gtf file, see the https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162208 (data processing entries describe how the downloadable file ma1_genes_longer.gtf was generated).

Processing cell ranger output: hamster_merging.R contains the code to create the various Seurat objects.

Annotation of lung cell types from tabula muris / Travaglini et al.: lung_hamster_annotation.Rmd

Processing and figures for Fig. 1 and Fig. S3: hamster_scRNAseq.R contains the  code to annote the cell types and create the panels in Fig. 1 and Fig. S3. For convenience, the working object ma_int.rds can be downloaded via http://www.mdc-berlin.de/singlecell-SARSCoV2

Pseudobulk analysis: 

