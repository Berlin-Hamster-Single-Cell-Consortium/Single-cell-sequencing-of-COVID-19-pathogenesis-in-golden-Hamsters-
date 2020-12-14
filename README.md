# -Single-cell-sequencing-of-COVID-19-pathogenesis-in-golden-Hamsters-
Single cell sequencing is a powerful tool to investigate cellular mechanisms of disease pathogenesis. The golden Hamster (Mesocricetus auratus) is a valuable model for SARS-CoV-2 infection as molecular mechanisms seem comparable to humans.

Abstract:

Coronavirus disease 2019 - COVID-19 - a multi-facetted disease initiated by respiratory tract infection with severe acute respiratory syndrome corona virus 2 (SARS-CoV-2), is partially driven by immune responses which are an important determinant of disease severity and may be key to therapeutic strategies. However, exact pulmonary cellular and molecular mechanisms contributing to inflammatory response, lung injury and repair are difficult to investigate in patients. Syrian hamsters develop spontaneously resolving moderate disease with notable pneumonia and systemic clinical symptoms upon SARS-CoV-2 infection. However, a lack of immunological reagents impedes mechanistic analyses. We thus investigated the spatiotemporal pathogen-host interaction in SARS-CoV-2 infected Syrian hamsters by means of single-cell sequencing technologies and ultra-high through-put proteomics. Thereby, details of immune cell recruitment, immune effector functions at the site of infection, cellular subtype-specific responses and their reaction to viral intra- and extracellular exposure as well as initiation of tissue repair were dissected. Direct comparison with according data on human COVID-19 pathogenesis revealed strong concordance regarding pathophysiologic sequences. We here show that the course of disease in SARS-CoV-2 infected Syrian hamsters closely resembles moderate human COVID-19 and yields deep insight into pathogen-host interaction. Our data enhance understanding of COVID-19 pathophysiology and support SARS-CoV-2 infected Syrian hamsters as pertinent model for research on COVID-19.

Keywords: coronavirus disease 2019, COVID-19; Syrian hamster; SARS-CoV2, inflammatory monocytes, macrophages, cytotoxic T cells, endothelial cells 

Infection Scheme Figure was created with Biorender.com

# table of content / data processing

Running cell ranger
Raw fastq files can be downloaded from GEO, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162208
To create a cell ranger reference, see the description of the mkref command here: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references
We have used a combined hamster/virus fasta/gtf file. For creating and downloading the gtf file, see the https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162208 (data processing entries describe how the downloadable file ma1_genes_longer.gtf was generated).

Processing cell ranger output
File hamster_merging.R contains the code to create the various Seurat objects.

Annotation of lung cell types from tabula muris / Travaglini et al.
lung_hamster_annotation.Rmd

http://www.mdc-berlin.de/singlecell-SARSCoV2
