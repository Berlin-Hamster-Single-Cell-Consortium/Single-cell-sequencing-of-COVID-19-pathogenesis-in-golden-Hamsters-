
library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)
library(ggplot2)
library(DESeq2)
#-------------------------------------------------------------------------------------------

#### PERFORM DIFFERENTIAL GENE EXPRESSION ANALYSIS WITH DESEQ2 PSEUDOBULK METHOD ####

Idents(hamster) <- hamster@meta.data$celltype
unique(hamster$celltype)

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
                      donor=gsub('[^_]*_ma_([de][0-9]*)_blood_([0-9]*)','hamster_\\2',names(expr)),
                      day=gsub('[^_]*_ma_([de][0-9]*)_blood_([0-9]*)','\\1',names(expr)),
                      group=names(expr),
                      row.names=names(expr))

counts <- do.call(cbind,expr)


clusters_to_check=c('T', "Activated T",   'B', "Classical monocyte", "Non-classical monocyte", "Neutrophil", "Immature neutrophil", "NK", "mDC", "pDC", "Platelet", "Mixed" )#, 'AlveolarMacrophages', 'Endothelial', 'Fibroblasts', 'Monocytes', 'SmoothMuscle', 'MyeloidDendritic', 'AT2', 'Macrophages', 'Ciliated', 'MonocytesAce+', 'AT1', 'Neutrophils', 'ProliferatingMacrophages', 'PlasmacytoidDendritic', 'Myofibroblast')

res <- list()
for (cluster in c('all',clusters_to_check)) {
  #for (cluster in c('all')) {
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
                                                  type='normal',
                                                  format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene') %>%
      dplyr::mutate(cluster=cluster,contrast='d2vsd0')
    res[[paste0(cluster,'_d3vsd0')]] <- lfcShrink(dds,
                                                  contrast=c('day','d3','d0'),
                                                  type='normal',
                                                  format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene') %>%
      dplyr::mutate(cluster=cluster,contrast='d3vsd0')
    res[[paste0(cluster,'_d3vsd2')]] <- lfcShrink(dds,
                                                  contrast=c('day','d3','d2'),
                                                  type='normal',
                                                  format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene') %>%
      dplyr::mutate(cluster=cluster,contrast='d3vsd2')
    res[[paste0(cluster,'_d5vsd0')]] <- lfcShrink(dds,
                                                  contrast=c('day','d5','d0'),
                                                  type='normal',
                                                  format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene') %>%
      dplyr::mutate(cluster=cluster,contrast='d5vsd0')
    res[[paste0(cluster,'_e14vsd0')]] <- lfcShrink(dds,
                                                   contrast=c('day','e14','d0'),
                                                   type='normal',
                                                   format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene') %>%
      dplyr::mutate(cluster=cluster,contrast='e14vsd0')
  })
}

pseudobulk <- do.call(rbind,res)



#### CREATE DOTPLOTS TO VISUALISE RESULTS ####

## Rename celltypes 
pseudobulk$cluster <- ifelse(pseudobulk$cluster == "T", "Tcells",
                             ifelse(pseudobulk$cluster == "B", "Bcells", 
                                    ifelse(pseudobulk$cluster == "Classical monocyte", "Classical monocytes",
                                           ifelse(pseudobulk$cluster == "Non-classical monocyte","Non-classical monocytes",
                                                  ifelse(pseudobulk$cluster == "Neutrophil", "Neutrophils",
                                                         ifelse(pseudobulk$cluster == "Immature neutrophil", "Immature neutrophils",
                                                                ifelse(pseudobulk$cluster == "pDC", "PlasmacytoidDendritic",
                                                                       ifelse(pseudobulk$cluster == "mDC", "MyeloidDendritic",
                                                                              ifelse(pseudobulk$cluster == "Activated T", "Activated Tcells",
                                                                                     ifelse(pseudobulk$cluster == "Mixed", "Unclear",
                                                                                            ifelse(pseudobulk$cluster == "Platelet", "Platelet", 
                                                                                                   ifelse(pseudobulk$cluster == "NK", "NK",
                                                                                                          ifelse(pseudobulk$cluster =="all", "all", "none")))))))))))))



#### MAKE DOTPLOTS ####

##change "gene column" to "gene_name"
head(pseudobulk)
pseudobulk <- pseudobulk %>% dplyr::rename(gene_name = gene)
names(pseudobulk)

## CYTOKINES

#get list of cytokines
cytokines <- read.table("C:/path/to/cytokines_ma_v2.txt", header=TRUE, sep="\t")

####


for (contrasts_to_use in c("d2vsd0", "d3vsd0", "d5vsd0", "e14vsd0")) {
  
  #Select cytokines
  genes <- unique(cytokines$gene)
  clusters <- the_celltypes
  
  
  if (contrasts_to_use != "e14vsd0"){
    
    
    df <- pseudobulk %>%
      dplyr::filter(cluster %in% clusters) %>%
      dplyr::filter(contrast %in% contrasts_to_use) %>%
      dplyr::filter(!(is.na(padj)) & (padj < .01)) %>%
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
             dplyr::filter(cluster %in% the_celltypes) %>%
             dplyr::filter(contrast %in% contrasts_to_use) %>%
             dplyr::filter(!(is.na(padj)) & (padj < .01)) %>%
             dplyr::filter(gene_name %in% genes) %>%
             dplyr::mutate(padj=ifelse(padj < 1E-20, 1E-20, padj)) %>%
             dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
             dplyr::mutate(gene=factor(gene_name,levels=gene.order))%>%
             dplyr::mutate(cluster=factor(cluster,levels=group.order)),
           aes(y=gene,x=cluster,size=-log10(padj),color=log2FoldChange)) +
      geom_point() +
      scale_size_continuous(name='adj. p-value',
                            breaks=c(2,4,6),
                            labels=c("1e-2","1e-4", "<1e-6"),
                            limits=c(0,NA)) +
      scale_color_gradient2(low='blue3',mid='gray',high='red3',
                            limits=c(-3,3),oob=scales::squish) +
      #facet_wrap(~contrast) +
      #coord_cartesian(clip = 'off') +
      theme_minimal(base_size=10) +
      guides(color=guide_colorbar(barheight=4)) +
      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
      labs(x='',y='')
    
    #save 8x12 inches
    ggsave(paste("pseudobulk_cytokines", contrasts_to_use, "pdf", sep="."), width=6, height = length(genes)/6+2, units="in", useDingbats=FALSE)
    ggsave(paste("pseudobulk_cytokines", contrasts_to_use, "png", sep="."), width=6, height = length(genes)/6+2, units="in") 
    
  } else {
    
    df <- pseudobulk %>%
      dplyr::filter(cluster %in% the_celltypes) %>%
      dplyr::filter(contrast %in% contrasts_to_use) %>%
      #dplyr::filter(!(is.na(padj)) & (padj < .01)) %>%
      dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
      dplyr::filter(gene_name %in% genes) %>%
      dplyr::mutate(group=paste0(cluster,'_',contrast))  %>%
      dplyr::select(group,gene_name,log2FoldChange) %>%
      spread(group,log2FoldChange) %>%
      tibble::column_to_rownames('gene_name')
    
    
    df
    
    df[is.na(df)] <- 0
    hc <- hclust(dist(df))
    gene.order <- row.names(df)[order.hclust(hc)]
    
    
    genes <- row.names(df)
    
    ggplot(pseudobulk %>%
             dplyr::filter(cluster %in% clusters) %>%
             dplyr::filter(contrast %in% contrasts_to_use) %>%
             #dplyr::filter(!(is.na(padj) & (padj < .01))) %>%
             dplyr::filter(gene_name %in% genes) %>%
             dplyr::mutate(padj=ifelse(padj < 1E-20, 1E-20, padj)) %>%
             dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
             dplyr::mutate(gene=factor(gene_name,levels=gene.order)),
           aes(y=gene_name,x=cluster,size=-log10(padj),color=log2FoldChange)) +
      geom_point() +
      scale_size_continuous(name='adj. p-value',
                            breaks=c(2,4,6),
                            labels=c("1e-2","1e-4", "<1e-6"),
                            limits=c(0,NA)) +
      scale_color_gradient2(low='blue3',mid='gray',high='red3',
                            limits=c(-3,3),oob=scales::squish) +
      #facet_wrap(~contrast) +
      #coord_cartesian(clip = 'off') +
      theme_minimal(base_size=10) +
      guides(color=guide_colorbar(barheight=4)) +
      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
      labs(x='',y='')
    
    
    
    #save 8x12 inches
    ggsave(paste("pseudobulk_cytokines", contrasts_to_use, "pdf", sep="."), width=6, height = length(genes)/6+2, units="in", useDingbats=FALSE)
    ggsave(paste("pseudobulk_cytokines", contrasts_to_use, "png", sep="."), width=6, height = length(genes)/6+2, units="in") 
    
    
  }
  
}


#Make with reduced set of cytokines


unique(pseudobulk$cluster)

the_celltypes = c("Classical monocytes",
                  "Non-classical monocytes",
                  "Neutrophils",
                  "Immature neutrohils",
                  "MyeloidDendritic",
                  "NK",
                  "Tcells",
                  "Activated Tcells",
                  "Bcells")

genes <- c("Tgfb3", "Tnfsf14", "Tnfsf10", "Il10", "Il6", "Il1b", "Il1a", "Csf1", "Cxcl13", "Cxcl11", "Cxcl10", "Ccl12", "Ccl8", "Ccl7", "Ccl5", "Ccl4", "Ccl3", "Ccl2")

for (contrasts_to_use in c("d2vsd0", "d3vsd0", "d5vsd0", "e14vsd0")) {
  
  
  df <- pseudobulk %>%
    dplyr::filter(cluster %in% the_celltypes) %>%
    dplyr::filter(contrast %in% contrasts_to_use) %>%
    dplyr::filter(!(is.na(padj)) & (padj < .01)) %>%
    dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
    dplyr::filter(gene_name %in% genes) %>%
    dplyr::mutate(group=paste0(cluster,'_',contrast))  %>%
    dplyr::select(group,gene_name,log2FoldChange) %>%
    spread(group,log2FoldChange) %>%
    tibble::column_to_rownames('gene_name')
  
  

  ggplot(pseudobulk %>%
           dplyr::filter(cluster %in% the_celltypes) %>%
           dplyr::filter(contrast %in% contrasts_to_use) %>%
           dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
           dplyr::filter(gene_name %in% genes) %>%
           dplyr::mutate(padj=ifelse(padj < 1E-20, 1E-20, padj)) %>%
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
  ggsave(paste("pseudobulk_cytokines_reduced", contrasts_to_use, "pdf", sep="."), width=6, height = length(genes)/6+2, units="in", useDingbats=FALSE)
  ggsave(paste("pseudobulk_cytokines_reduced", contrasts_to_use, "png", sep="."), width=6, height = length(genes)/6+2, units="in")
  
}


## GLOBAL DIFFERENTIAL GENE EXPRESSION

###top 15 per celltype


for (contrasts_to_use in c("d2vsd0", "d3vsd0", "d5vsd0", "e14vsd0")) {
  
  #Select for each contrast and each cell type the 15 most changing genes (derived vom p-padj)
  
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
  
  
  #Adjust data as following: all padj less than 1E-20 are corrected to 1E-20, all NA values to 1. Since these values are -log10
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
                          limits=c(-3,3),oob=scales::squish) +
    facet_wrap(~contrast) +
    coord_cartesian(clip = 'off') +
    theme_minimal(base_size=10) +
    guides(color=guide_colorbar(barheight=4)) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
    labs(x='',y='')
  
  #save 8x12 inches
  h=length(gene.order) / 5
  ggsave(paste("pseudobulk_top15", contrasts_to_use, "pdf", sep="."), width=6, height = h+2, units="in", useDingbats=FALSE)
  ggsave(paste("pseudobulk_top15", contrasts_to_use, "png", sep="."), width=6, height = h+2, units="in")
  
}



###top 4 per celltype

for (contrasts_to_use in c("d2vsd0", "d3vsd0", "d5vsd0", "e14vsd0")) {
  
  #Select for each contrast and each cell type the 15 most changing genes (derived vom p-padj)
  
  #For the small version in Fig. 2 use n_4 in slice_min and save as pseudobulk_small
  genes <- pseudobulk %>%
    dplyr::filter(cluster %in% clusters) %>%
    dplyr::filter(contrast %in% contrasts_to_use) %>%
    dplyr::filter(!grepl('SCoV2',gene_name)) %>%
    dplyr::filter(!(is.na(padj)) & (padj < .01)) %>%
    dplyr::group_by(contrast, cluster) %>%
    dplyr::slice_min(order_by=padj,n=4,with_ties=FALSE) %>%
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
  
  
  #Adjust data as following: all padj less than 1E-20 are corrected to 1E-20, all NA values to 1. Since these values are -log10
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
                          limits=c(-3,3),oob=scales::squish) +
    facet_wrap(~contrast) +
    coord_cartesian(clip = 'off') +
    theme_minimal(base_size=10) +
    guides(color=guide_colorbar(barheight=4)) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
    labs(x='',y='')
  
  #save 8x12 inches
  h=length(gene.order) / 5
  ggsave(paste("pseudobulk_small", contrasts_to_use, "pdf", sep="."), width=6, height = h+2, units="in", useDingbats=FALSE)
  ggsave(paste("pseudobulk_small", contrasts_to_use, "png", sep="."), width=6, height = h+2, units="in")
  
}








