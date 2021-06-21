library(plyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)
library(Seurat)
library(HDF5Array)
library(rhdf5)
library(hdf5r)
library(DelayedArray)
library(DelayedMatrixStats)
require(gridExtra)
library(tidyr)
library("dendextend")
library(DESeq2)
source("~/Documents/Largescale-data/notes and scripts/smooth_DimPlot.R")
library(pheatmap)
#see http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)
source("~/Documents/Largescale-data/notes and scripts/summarySE.R")

##############
#start out with the annotated object ma_int.rds described in lung_hamster_scRNAseq.R
hamster <- readRDS("./ma_int.rds")
Idents(hamster) <- hamster@meta.data$celltype

#################
#Prepare pseudobulk table

expr <- list()
for (cluster in unique(Idents(hamster))) {
  for (sample in unique(hamster@meta.data$orig.ident)) {
    cells <- Cells(hamster)[(hamster@meta.data$orig.ident==sample) & (Idents(hamster)==cluster)]
    if (length(cells) > 20) { 
      expr[[paste0(cluster,'_',sample)]] <- rowSums(hamster@assays$RNA@counts[,cells])
    }
  }
}
for (sample in unique(hamster@meta.data$orig.ident)) {
  cells <- Cells(hamster)[(hamster@meta.data$orig.ident==sample)]
  expr[[paste0('all_',sample)]] <- rowSums(hamster@assays$RNA@counts[,cells])
}

colData <- data.frame(cell.type=factor(gsub('_.*','',names(expr))),
                      donor=gsub('[^_]*_ma_([de][0-9]*)_lung_([0-9]*)','hamster_\\2',names(expr)),
                      day=gsub('[^_]*_ma_([de][0-9]*)_lung_([0-9]*)','\\1',names(expr)),
                      group=names(expr),
                      row.names=names(expr))

counts <- do.call(cbind,expr)

clusters_to_check=c('Tcells', 'Bcells', 'AlveolarMacrophages', 'Endothelial', 'Fibroblasts', 'MonocyticMacrophages', 'AT2', 'InterstitialMacrophages', 'MacrophagesTreml4+', 'AT1', 'Neutrophils')

res <- list()
for (cluster in c('all',clusters_to_check)) {
  take_row <- rowSums(counts) > 0
  take_col <- (colSums(counts) > 0) & (colData[,'cell.type']==cluster)
  try({
    dds <- DESeqDataSetFromMatrix(countData=counts[take_row,take_col],
                                  colData=colData[take_col,,drop=FALSE],
                                  design=~day)
    if (cluster!='all')
      dds <- estimateSizeFactors(dds, type='poscounts')
    dds <- DESeq(dds)
    res[[paste0(cluster,'_d2vsd0')]] <- lfcShrink(dds,
                                                  contrast=c('day','d2','d0'),
                                                  method='normal',
                                                  format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene_name') %>%
      dplyr::mutate(cluster=cluster,contrast='d2vsd0')
    res[[paste0(cluster,'_d3vsd0')]] <- lfcShrink(dds,
                                                  contrast=c('day','d3','d0'),
                                                  method='normal',
                                                  format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene_name') %>%
      dplyr::mutate(cluster=cluster,contrast='d3vsd0')
    res[[paste0(cluster,'_d3vsd2')]] <- lfcShrink(dds,
                                                  contrast=c('day','d3','d2'),
                                                  method='normal',
                                                  format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene_name') %>%
      dplyr::mutate(cluster=cluster,contrast='d3vsd2')
    res[[paste0(cluster,'_d5vsd0')]] <- lfcShrink(dds,
                                                  contrast=c('day','d5','d0'),
                                                  method='normal',
                                                  format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene_name') %>%
      dplyr::mutate(cluster=cluster,contrast='d5vsd0')
    res[[paste0(cluster,'_e14vsd0')]] <- lfcShrink(dds,
                                                   contrast=c('day','e14','d0'),
                                                   method='normal',
                                                   format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene_name') %>%
      dplyr::mutate(cluster=cluster,contrast='e14vsd0')
  })
}

pseudobulk <- do.call(rbind,res)

write.table(pseudobulk, "./pseudobulk.txt", row.names = TRUE, sep="\t", quote=FALSE)
pseudobulk <- read.table(".pseudobulk.txt", sep="\t", header=TRUE)
######################################
#Dotplots for pseudobulk, e.g. Fig. S5

clusters=c('AlveolarMacrophages', 'InterstitialMacrophages', 'MonocyticMacrophages', 'MacrophagesTreml4+', 'Neutrophils', 'Tcells', 'Bcells', 'AT2', 'AT1', 'Fibroblasts', 'Endothelial')

for (contrasts_to_use in c("d2vsd0", "d3vsd0", "d5vsd0", "e14vsd0")) {
#for (contrasts_to_use in c("d2vsd0")) {

  #Select top genes
  #For the small version in Fig. 2 use n_4 in slice_min and save as pseudobulk_small
  genes <- pseudobulk %>%
    dplyr::filter(cluster %in% clusters) %>%
    dplyr::filter(contrast %in% contrasts_to_use) %>%
    dplyr::filter(!grepl('SCoV2',gene_name)) %>%
    dplyr::filter(!(is.na(padj)) & (padj < .01)) %>%
    dplyr::group_by(contrast, cluster) %>%
    dplyr::slice_min(order_by=padj,n=15,with_ties=FALSE) %>%
    dplyr::pull(gene_name)
  
  #cluster values to get proper order for genes (rows) and cell types (columns)
  df <- pseudobulk %>%
    dplyr::filter(cluster %in% clusters) %>%
    dplyr::filter(contrast %in% contrasts_to_use) %>%
    dplyr::filter(gene_name %in% genes) %>%
    dplyr::mutate(group=paste0(cluster,'_',contrast))  %>%
    dplyr::select(group,gene_name,log2FoldChange) %>%
    spread(group,log2FoldChange) %>%
    tibble::column_to_rownames('gene_name')
  
  df[is.na(df)] <- 0
  hc <- hclust(dist(df))
  gene.order <- row.names(df)[order.hclust(hc)]
  hc <- hclust(dist(t(df)))
  group.order <- colnames(df)[order.hclust(hc)]
  group.order = sub("(.+)_.*", "\\1", group.order)

  ggplot(pseudobulk %>%
         dplyr::filter(cluster %in% clusters) %>%
         dplyr::filter(contrast %in% contrasts_to_use) %>%
         dplyr::filter(gene_name %in% genes) %>%
         dplyr::mutate(padj=ifelse(padj < 1E-20, 1E-20, padj)) %>%
         dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
         dplyr::mutate(gene=factor(gene_name,levels=gene.order)) %>%
         dplyr::mutate(cluster=factor(cluster,levels=group.order)),
         aes(y=gene,x=cluster,size=-log10(padj),color=log2FoldChange)) +
    geom_point(shape=16) +
    scale_size_continuous(name='adj. p-value',
                          breaks=c(5,10,15),
                          labels=c("1e-5","1e-10", "<1e-15"),
                          limits=c(0,20)) +
    scale_color_gradient2(low='blue3',mid='gray',high='red3',
                          limits=c(-6,6),oob=scales::squish) +
    facet_wrap(~contrast) +
    coord_cartesian(clip = 'off') +
    theme_minimal(base_size=10) +
    guides(color=guide_colorbar(barheight=4)) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
    labs(x='',y='')
  
  #save 8x12 inches
  h=length(gene.order) / 5
  ggsave(paste("pseudobulk", contrasts_to_use, "pdf", sep="."), width=6, height = h+2, units="in", useDingbats=FALSE)

}

###############
#Plot with inflammatory mediators


cytokines <- read.table("./cytokines_ma_v2.txt", header=TRUE, sep="\t")


contrasts_to_use <- "d2vsd0"
#for (contrasts_to_use in c("d2vsd0", "d3vsd0", "d5vsd0", "e14vsd0")) {
contrasts_to_use <- c("d2vsd0", "d3vsd0", "d5vsd0", "e14vsd0")
  #Select cytokines
  genes <- unique(cytokines$gene)
  
  df <- pseudobulk %>%
    dplyr::filter(cluster %in% clusters) %>%
    dplyr::filter(contrast %in% contrasts_to_use) %>%
    dplyr::filter(!(is.na(padj)) & (padj < .0001) & (abs(log2FoldChange) >= 1)) %>%
    dplyr::filter(gene_name %in% genes) %>%
    dplyr::mutate(group=paste0(cluster,'_',contrast))  %>%
    dplyr::select(group,gene_name,log2FoldChange) %>%
    spread(group,log2FoldChange) %>%
    tibble::column_to_rownames('gene_name')
  
  df[is.na(df)] <- 0
  hc <- hclust(dist(df))
  gene.order <- row.names(df)[order.hclust(hc)]
  hc <- hclust(dist(t(df)))
  group.order <- colnames(df)[order.hclust(hc)]
  group.order = sub("(.+)_.*", "\\1", group.order)
  
  genes <- row.names(df)
  
  ggplot(pseudobulk %>%
           dplyr::filter(cluster %in% clusters) %>%
           dplyr::filter(contrast %in% contrasts_to_use) %>%
           dplyr::filter(gene_name %in% genes) %>%
           dplyr::mutate(padj=ifelse(padj < 1E-20, 1E-20, padj)) %>%
           dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
           dplyr::mutate(gene=factor(gene_name,levels=gene.order)) %>%
           dplyr::mutate(cluster=factor(cluster,levels=clusters)),
         aes(y=gene,x=cluster,size=-log10(padj),color=log2FoldChange)) +
    geom_point(shape=16) +
    scale_size_continuous(name='adj. p-value',
                          breaks=c(2,4,6),
                          labels=c("1e-2","1e-4", "<1e-6"),
                          limits=c(0,NA)) +
    scale_color_gradient2(low='blue3',mid='gray',high='red3',
                          limits=c(-3,3),oob=scales::squish) +
    facet_wrap(~contrast, nrow=1) +
    coord_cartesian(clip = 'off') +
    theme_minimal(base_size=10) +
    guides(color=guide_colorbar(barheight=4)) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
    labs(x='',y='')
  
  #save 8x12 inches
  ggsave(paste("pseudobulk_cytokines_all", "pdf", sep="."), width=10, height = length(genes)/6+2, units="in", useDingbats=FALSE)
  
#}

#Make with reduced set of cytokines
genes <- c("Tgfb3", "Tnfsf14", "Tnfsf10", "S100a9", "Il10", "Il6", "Il1b", "Il1a", "Csf1", "Cxcl13", "Cxcl11", "Cxcl10", "Cxcl5", "ENSMAUG00000016910", "Ccl12", "Ccl8", "Ccl7", "Ccl5", "Ccl4", "Ccl3", "Ccl2")

the_celltypes = c("Tcells",
                  "MonocyticMacrophages",
                  "Endothelial",
                  "AlveolarMacrophages",
                  "AT2",
                  "Bcells",
                  "MacrophagesTreml4+",
                  "Fibroblasts",
                  "MyeloidDendritic",
                  "AT1",
                  "Neutrophils",
                  "SmoothMuscle",
                  "InterstitialMacrophages",
                  "PlasmacytoidDendritic",
                  "Ciliated",
                  "Unclear",
                  "Myofibroblast")


for (contrasts_to_use in c("d2vsd0", "d3vsd0", "d5vsd0", "e14vsd0")) {
  
 
  df <- pseudobulk %>%
    dplyr::filter(cluster %in% clusters) %>%
    dplyr::filter(contrast %in% contrasts_to_use) %>%
    dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
    dplyr::filter(gene_name %in% genes) %>%
    dplyr::mutate(group=paste0(cluster,'_',contrast))  %>%
    dplyr::select(group,gene_name,log2FoldChange) %>%
    spread(group,log2FoldChange) %>%
    tibble::column_to_rownames('gene_name')
  
  # df[is.na(df)] <- 0
  # hc <- hclust(dist(df))
  # gene.order <- row.names(df)[order.hclust(hc)]
  # hc <- hclust(dist(t(df)))
  # group.order <- colnames(df)[order.hclust(hc)]
  # group.order = sub("(.+)_.*", "\\1", group.order)
  # 
  # genes <- row.names(df)
  
  ggplot(pseudobulk %>%
           dplyr::filter(cluster %in% clusters) %>%
           dplyr::filter(contrast %in% contrasts_to_use) %>%
           dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
           dplyr::filter(gene_name %in% genes) %>%
           dplyr::mutate(padj=ifelse(padj < 1E-20, 1E-20, padj)) %>%
           dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
           dplyr::mutate(gene=factor(gene_name,levels=genes))%>%
           dplyr::mutate(cluster=factor(cluster,levels=the_celltypes)),
         aes(y=gene,x=cluster,size=-log10(padj),color=log2FoldChange)) +
    geom_point() +
    scale_size_continuous(name='adj. p-value',
                          breaks=c(2,4,6),
                          labels=c("1e-2","1e-4", "<1e-6"),
                          limits=c(0,NA)) +
    scale_color_gradient2(low='blue3',mid='gray',high='red3',
                          limits=c(-3,3),oob=scales::squish) +
    facet_wrap(~contrast) +
    coord_cartesian(clip = 'off') +
    theme_minimal(base_size=10) +
    guides(color=guide_colorbar(barheight=4)) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
    labs(x='',y='')
  
  #save 8x12 inches
  ggsave(paste("pseudobulk_cytokines_reduced", contrasts_to_use, "pdf", sep="."), width=5, height = length(genes)/6+2, units="in", useDingbats=FALSE)
  
}

####################
#Do a plot with only Cxcl10 but all time points in one

ggplot(pseudobulk %>%
         dplyr::filter(cluster %in% clusters) %>%
         dplyr::filter(!(contrast == "d3vsd2")) %>%
         dplyr::filter(gene_name %in% c("Cxcl10")) %>%
         dplyr::mutate(padj=ifelse(padj < 1E-20, 1E-20, padj)) %>%
         dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
         mutate(contrast = forcats::fct_relevel(contrast, c("e14vsd0", "d5vsd0", "d3vsd0", "d2vsd0"))) %>%
         dplyr::mutate(gene=factor(gene_name,levels=gene.order)),
       aes(y=contrast,x=cluster,size=-log10(padj),color=log2FoldChange)) +
  geom_point() +
  scale_size_continuous(name='adj. p-value',
                        breaks=c(2,4,6),
                        labels=c("1e-2","1e-4", "<1e-6"),
                        limits=c(0,NA)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-3,3),oob=scales::squish) +
  coord_cartesian(clip = 'off') +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
  labs(x='',y='')
ggsave(paste("pseudobulk_cytokines", "Cxcl10", "pdf", sep="."), width=5, height = 3, units="in", useDingbats=FALSE)



################################
#Pseudobulk by comparing with/without virus at d2, d3, d5 (Fig. 6)

Idents(hamster) <- hamster@meta.data$celltype


expr <- list()
for (cluster in unique(Idents(hamster))) {
    for (sample in c("ma_d2_lung_1", "ma_d2_lung_2", "ma_d2_lung_3", "ma_d3_lung_1", "ma_d3_lung_2", "ma_d3_lung_3", "ma_d5_lung_1", "ma_d5_lung_2", "ma_d5_lung_3")) {
      cells <- Cells(hamster)[(hamster@meta.data$orig.ident==sample) & (Idents(hamster)==cluster) & (hamster@meta.data$SCoV2_load > 0)]
      if (length(cells) > 10) { 
        expr[[paste0(cluster,'_',sample,'_withvirus')]] <- rowSums(hamster@assays$RNA@counts[,cells])
      }
      cells <- Cells(hamster)[(hamster@meta.data$orig.ident==sample) & (Idents(hamster)==cluster) & (hamster@meta.data$SCoV2_load == 0)]
      if (length(cells) > 10) { 
        expr[[paste0(cluster,'_',sample,'_novirus')]] <- rowSums(hamster@assays$RNA@counts[,cells])
      }
    }
}


colData <- data.frame(cell.type=factor(gsub('([^_]*)_([^_]*)_([^_]*)_([^_]*)_([^_]*)_([^_]*)','\\1',names(expr))),
                      virus=gsub('([^_]*)_([^_]*)_([^_]*)_([^_]*)_([^_]*)_([^_]*)','\\6',names(expr)),
                      day=gsub('([^_]*)_([^_]*)_([^_]*)_([^_]*)_([^_]*)_([^_]*)','\\3',names(expr)),
                      group=names(expr),
                      row.names=names(expr))

counts <- do.call(cbind,expr)

clusters_to_check=c('MonocyticMacrophages')

res <- list()
for (cluster in c(clusters_to_check)) {
  for (timepoint in c("d2", "d3", "d5")) {
    take_row <- rowSums(counts) > 0
    take_col <- (colSums(counts) > 0) & (colData[,'cell.type']==cluster) & (colData[,'day']==timepoint)
    try({
      dds <- DESeqDataSetFromMatrix(countData=counts[take_row,take_col],
                                    colData=colData[take_col,,drop=FALSE],
                                    design=~virus)
      dds <- DESeq(dds)
      res[[paste(cluster,timepoint,'withnovirus', sep="_")]] <- lfcShrink(dds,
                                                    contrast=c('virus','withvirus','novirus'),
                                                    method='normal',
                                                    format="DataFrame") %>%
        as.data.frame() %>%
        tibble::rownames_to_column('gene_name') %>%
        dplyr::mutate(cluster=cluster,contrast=paste(timepoint,'withnovirus', sep="_"))
    })
  }
}

pseudobulk <- do.call(rbind,res)

genes <- pseudobulk %>%
  dplyr::filter(!grepl('SCoV2',gene_name)) %>%
  dplyr::filter(!(is.na(padj)) & (padj < .01)) %>%
  dplyr::group_by(contrast, cluster) %>%
  dplyr::slice_min(order_by=padj,n=100,with_ties=FALSE) %>%
  dplyr::pull(gene_name)

#For reduced figure use:
genes <- c("Ccl2", "Ccl8", "Ccl12", "Cxcl10", "Cxcl11", "Ccl3", "Ccl4", "Mx2", "Ifit2", "Isg15", "Ccr5", "Siglec1")

df <- pseudobulk %>%
  dplyr::filter(gene_name %in% genes) %>%
  dplyr::mutate(group=contrast)  %>%
  dplyr::select(group,gene_name,log2FoldChange) %>%
  spread(group,log2FoldChange) %>%
  tibble::column_to_rownames('gene_name')

df[is.na(df)] <- 0
hc <- hclust(dist(df))
gene.order <- row.names(df)[order.hclust(hc)]
hc <- hclust(dist(t(df)))
group.order <- colnames(df)[order.hclust(hc)]
group.order = sub("(.+)_.*", "\\1", group.order)

ggplot(pseudobulk %>%
         dplyr::filter(gene_name %in% genes) %>%
         dplyr::mutate(padj=ifelse(padj < 1E-20, 1E-20, padj)) %>%
         dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
         dplyr::mutate(gene=factor(gene_name,levels=gene.order)),
       aes(y=gene,x=contrast,size=-log10(padj),color=log2FoldChange)) +
  geom_point() +
  scale_size_continuous(name='adj. p-value',
                        breaks=c(5,10,15),
                        labels=c("1e-5","1e-10", "<1e-15"),
                        limits=c(0,NA)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-6,6),oob=scales::squish) +
  coord_cartesian(clip = 'off') +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
  labs(x='',y='')

#save 8x12 inches
h=length(gene.order) / 6
ggsave(paste("pseudobulk_MoMa_withnovirus", "pdf", sep="."), width=6, height = h+2, units="in", useDingbats=FALSE)

#Read in exports from string-db.org
#For Fig. S7
#for downregulated no significant KEGG or GO BP was found
#remove leading # and ' before reading in
KEGG <- read.table("./enrichment.KEGG.strongup.tsv", sep="\t", header=TRUE) 
GO <- read.table("./enrichment.Process.strongup.tsv", sep="\t", header=TRUE) 
strongup <- rbind.data.frame(GO, KEGG)
strongup$cluster <- "strongup"
KEGG <- read.table("./enrichment.KEGG.weakup.tsv", sep="\t", header=TRUE) 
GO <- read.table("./enrichment.Process.weakup.tsv", sep="\t", header=TRUE) 
weakup <- rbind.data.frame(GO, KEGG)
weakup$cluster <- "weakup"

KEGG_GO <- rbind.data.frame(strongup, weakup)

terms <- KEGG_GO %>%
  dplyr::filter(!(is.na(false.discovery.rate)) & (false.discovery.rate < .00001)) %>%
  dplyr::pull(term.ID)


df <- KEGG_GO %>%
  dplyr::filter(term.ID %in% terms) %>%
  dplyr::mutate(group=cluster)  %>%
  dplyr::select(group,term.ID,strength) %>%
  spread(group,strength) %>%
  tibble::column_to_rownames('term.ID')

df[is.na(df)] <- 0
hc <- hclust(dist(df))
term.order <- row.names(df)[order.hclust(hc)]

ggplot(KEGG_GO %>%
         dplyr::filter(term.ID %in% terms) %>%
         dplyr::mutate(false.discovery.rate=ifelse(false.discovery.rate < 1E-20, 1E-20, false.discovery.rate)) %>%
         dplyr::mutate(false.discovery.rate=ifelse(is.na(false.discovery.rate), 1, false.discovery.rate)) %>%
         dplyr::mutate(term.ID=factor(term.ID,levels=term.order)),
       aes(y=term.description,x=cluster,size=-log10(false.discovery.rate),color=strength)) +
  geom_point() +
  scale_size_continuous(name='false discovery rate',
                        breaks=c(2,4,6),
                        labels=c("1e-2","1e-4", "<1e-6"),
                        limits=c(0,NA)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-2,2),oob=scales::squish) +
  coord_cartesian(clip = 'off') +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
  labs(x='',y='')
h=length(term.order) / 6
ggsave(paste("strongup_weakup_KEGG_GO", "pdf", sep="."), width=6, height = h+2, units="in", useDingbats=FALSE)
