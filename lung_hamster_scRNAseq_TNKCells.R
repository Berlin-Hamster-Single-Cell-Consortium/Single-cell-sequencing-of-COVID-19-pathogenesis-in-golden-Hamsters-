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

library(pheatmap)
source("./summarySE.R")
library(DESeq2)



###################
#Subclustering and analysis of T/NK cells (Fig. 7)

seu_T <- subset(hamster.integrated, subset = (celltype == "Tcells"))
DefaultAssay(seu_T) <- "integrated"                  
seu_T <- RunPCA(seu_T, verbose = FALSE)
seu_T <- RunUMAP(seu_T, dims = 1:30)
seu_T <- FindNeighbors(seu_T, dims = 1:30)
seu_T <- FindClusters(seu_T, resolution = 0.9)
UMAPPlot(seu_T, label=TRUE)
DefaultAssay(seu_T) <- "SCT"                  

a <- seu_T@meta.data %>% group_by(seurat_clusters, infection) %>% tally()

FeaturePlot(seu_T, "Cd4")
FeaturePlot(seu_T, "ENSMAUG00000000153")
FeaturePlot(seu_T, "Ifng")
FeaturePlot(seu_T, "Gzma")
FeaturePlot(seu_T, "Cd3e")


seu_T@meta.data$UMAP_1 <- NULL
seu_T@meta.data$UMAP_2 <- NULL
seu_T@meta.data = cbind(seu_T@meta.data, seu_T@reductions$umap@cell.embeddings)


ggplot()+geom_point(data=seu_T@meta.data, aes(x=UMAP_1, y=UMAP_2, color=seurat_clusters), size=0.5)+theme_bw()
ggplot()+geom_point(data=seu_T@meta.data, aes(x=UMAP_1, y=UMAP_2, color=infection), size=0.5, alpha=0.4)+theme_bw()
#
SCoV2_rawcounts <- FetchData(seu_T, grep("SCoV2", seu_T@assays$RNA@counts@Dimnames[[1]], value="TRUE"), slot="counts")
seu_T@meta.data$SCoV2_sum <- rowSums(SCoV2_rawcounts)
SCoV2_rawcounts$SCoV2sum <- rowSums(SCoV2_rawcounts)
seu_T@meta.data$SCoV2_load <- SCoV2_rawcounts$SCoV2sum/seu_T@meta.data$nCount_RNA*100
seu_T@meta.data$withvirus = ifelse(seu_T@meta.data$SCoV2_load>0,"withvirus","novirus")

saveRDS(seu_T, "/Applications/stuff/hamster/seu_red_Tcells.rds")

seu_T <- readRDS("/Applications/stuff/hamster/seu_red_Tcells.rds")


#Percent positive for a gene in clusters
gene <- "Cd3e"
df1 <- cbind.data.frame(FetchData(seu_T, gene),
                        seu_T@meta.data$seurat_clusters,
                        str_replace_all(gene, "-", "_"),
                        seu_T@meta.data$celltype,
                        seu_T@meta.data$infection)
colnames(df1) <- c("expression", "cluster", "gene", "celltype", "infection")
a = df1 %>% group_by(cluster, infection) %>% tally(name="tot") 
b = df1 %>% group_by(cluster, infection) %>% filter(expression>0) %>% tally(name="pos") 
barplotvalues =
  left_join(a , b , by = c('cluster', 'infection')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
ggplot()+
  geom_bar(data=barplotvalues, aes(x=cluster, y=fraction, fill=infection), stat="identity", position=position_dodge())+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("percent positiv", gene, "in", "clusters", sep=" "))



#Clusters:
Idents(seu_T) <- seu_T@meta.data$seurat_clusters
markers <- c("Cd3e", "Cd4", "ENSMAUG00000000153", "Gzma", "Nkg7")
markers1 <- c("Cd3e" ="Cd3e", "Cd4"="Cd4","ENSMAUG00000000153"= "Cd8a", "Gzma" = "Gzma", "Nkg7" = "Nkg7")
avg <- AverageExpression(seu_T, features = markers, return.seurat = T, assays = "SCT") 

DoHeatmap(avg, assay="SCT", size=5, features=markers, label=F, group.bar=F, group.bar.height = 0,
          draw.lines = F, hjust = 0)+ 
  scale_fill_gradientn(colors = c("blue", "white","red"), na.value = "white") +
  NoLegend()+
  theme_bw()+
  scale_y_discrete(labels=markers1, position = "right" )+
  scale_x_discrete( limits = rev(levels(seu_T$seurat_clusters)))+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, vjust = -2, hjust = -0.15))

ggsave("Tcells_ctype_markers_hmap.pdf")

seu_T <- RenameIdents(seu_T, 
                        '0'='CD4_Tcells',
                        '1'='CD4_Tcells',
                        '2'='NKcells',
                        '3'='NKcells',
                        '4'='CD4_Tcells',
                        '5'='CD4_Tcells',
                        '6'='ILC',
                        '7'='CD8_Tcells',
                        '8'='CD8_Tcells',
                        '9'='NKcells',
                        '10'='CD4_Tcells',
                        '11'='CD8_Tcells',
                        '12'='CD8_Tcells',
                        '13'='NKcells',
                        '14'='ILC')
seu_T@meta.data$celltype <- Idents(seu_T)

cluster_marks <- FindAllMarkers(seu_T)
write.table(cluster_markers, "Tcells_celltypes_markers.txt", sep="\t", quote=FALSE, row.names = FALSE)
top10 <- cluster_marks %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DotPlot(seu_T, group.by='celltype',
        features = unique(top10$gene, assay='RNA')) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))
ggsave("Tcells_celltypemarkers.pdf", width=12, height = 8, units="in", useDingbats=FALSE)


####################
#Plot with cell types
#Tcells in Fig. 1: #368F8B
celltypecolors = c("CD4_Tcells" = "#287F8B",
                   "CD8_Tcells" = "#40C9A2",
                   "NKcells" = "#A3F7B5",
                   "ILC" = "#81F0E5")
means <- seu_T@meta.data %>%
  group_by(celltype) %>%
  summarise(mean_U1 = mean(UMAP_1),
            mean_U2 = mean(UMAP_2))
ggplot()+
  geom_point(data=seu_T@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=celltype), size=0.5, stroke=0)+
  geom_text(data=subset(means, celltype!="Unclear"), aes(x=mean_U1, y=mean_U2, label=celltype))+
  coord_fixed(ratio=1)+
  ggtitle("celltype")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  scale_colour_manual(values = celltypecolors)
ggsave("Tcells_celltypes.pdf", useDingbats=FALSE)


#Make percentages over time
celltypecolors <- as.data.frame(celltypecolors)
celltypecolors <- celltypecolors[order(row.names(celltypecolors)),]
bright_colors = as.vector(celltypecolors)

#colors brighter with time
expr <- list()
for (cbright in bright_colors) {
  r=(col2rgb(cbright)+40)[[1]]
  r=ifelse(r>255,255,r)
  g=(col2rgb(cbright)+40)[[2]]
  g=ifelse(g>255,255,g)
  b=(col2rgb(cbright)+40)[[3]]
  b=ifelse(b>255,255,b)
  cdark = rgb(r, g, b, maxColorValue = 255)
  expr[[cbright]] <- scales::seq_gradient_pal(cbright, cdark, "Lab")(seq(0,1,length.out=5))
}
the_colors = as.vector(unlist(expr))


a = seu_T@meta.data %>% group_by(infection, hamster) %>% tally(name="tot") 
b = seu_T@meta.data %>% group_by(infection, celltype, hamster) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('infection', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
write.table(c, paste("Tcells_subtypes", "csv", sep="."), sep=",", quote=FALSE, row.names = FALSE)
tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("infection", "celltype"))


ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, infection, sep="_")))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.9))+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("Tcell subtype proportion", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  theme(legend.position = "none")+
  scale_fill_manual(values=the_colors)
ggsave("Tcells_percentages_celltypes.pdf", useDingbats=FALSE)

############################
#only d0 and d5

the_colors_d0d5 = c("#81F0E5", "#63D2C7", "#A3F7B5", "#85D997", "#40C9A2", "#24AB84", "#72E0AC", "#54C28E", "#287F8B", "#0C616D")

ggplot(subset(tgc, infection %in% c("d0", "d5")), aes(x=celltype, y=fraction, fill=paste(celltype, infection, sep="_")))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.9))+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("Tcell subtype proportion", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
#  theme(legend.position = "none")+
  scale_fill_manual(values=the_colors_d0d5)
ggsave("Tcells_percentages_celltypes_d0d5.pdf", useDingbats=FALSE)

#######################


#positive for gene in celltype
for (gene in c("Ifng", "Top2a", "Zap70", "Gzma")) {
  df1 <- cbind.data.frame(FetchData(seu_T, gene), seu_T@meta.data$seurat_clusters, seu_T@meta.data$celltype, str_replace_all(gene, "-", "_"), seu_T@meta.data$SCoV2_load, seu_T@meta.data$infection, seu_T@meta.data$hamster)
  colnames(df1) <- c("expression", "cluster", "celltype", "gene", "load", "day", "hamster")
  a = df1 %>% group_by(day, celltype, hamster) %>% tally(name="tot") 
  b = df1 %>% group_by(day, celltype, hamster) %>% filter(expression>0) %>% tally(name="pos") 
  c =
    left_join(a , b , by = c('day', 'celltype', 'hamster')) %>% 
    replace(., is.na(.), 0) %>%
    mutate(fraction = pos / tot)
  #%>%
  #  mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))
  tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("day", "celltype"))
  write.table(c, paste("Tcells_subtypes", gene, "csv", sep="."), sep=",", quote=FALSE, row.names = FALSE)

  
  #ggplot(subset(tgc, day %in% c("d0", "d5")), aes(x=celltype, y=fraction, fill=paste(celltype, day, sep="_")))+
  ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, day, sep="_")))+
    geom_bar(position=position_dodge(), stat="identity")+
    geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.9))+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    ggtitle(paste("fraction positive", gene, "in virus levels mean and se of three animals", sep=" "))+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
#    scale_fill_manual(values=the_colors_d0d5)
    scale_fill_manual(values=the_colors)

#  ggsave(paste("Tcells_percposbycelltype_d0d5", gene,  "pdf", sep="."), useDingbats=FALSE)
  ggsave(paste("Tcells_percposbycelltype", gene,  "pdf", sep="."), useDingbats=FALSE)
}


###############
#pseudobulk

hamster <- seu_T
Idents(hamster) <- seu_T@meta.data$celltype

clusters_to_check=c('CD4_Tcells', 'CD8_Tcells', 'NKcells', 'ILC')
expr <- list()
for (cluster in clusters_to_check) {
  for (sample in unique(hamster@meta.data$orig.ident)) {
    cells <- Cells(hamster)[(hamster@meta.data$orig.ident==sample) & (Idents(hamster)==cluster)]
    print(sample)
    print(cluster)
    print(length(cells))
    if (length(cells) > 0) { 
      expr[[paste0(cluster,'_',sample)]] <- rowSums(hamster@assays$RNA@counts[,cells])
    }
  }
}

colData <- data.frame(cell.type=factor(gsub('(.*)_ma_([de][0-9]*)_lung_([0-9]*)','\\1',names(expr))),
                      donor=gsub('.*_ma_([de][0-9]*)_lung_([0-9]*)','hamster_\\2',names(expr)),
                      day=gsub('.*_ma_([de][0-9]*)_lung_([0-9]*)','\\1',names(expr)),
                      group=names(expr),
                      row.names=names(expr))

counts <- do.call(cbind,expr)

res <- list()
for (cluster in c(clusters_to_check)) {
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

pseudobulk_T <- do.call(rbind,res)

write.table(pseudobulk_T, "/Applications/stuff/hamster/pseudobulk_T.txt", row.names = TRUE, sep="\t", quote=FALSE)
pseudobulk_T <- read.table("/Applications/stuff/hamster/pseudobulk_T.txt", sep="\t", header=TRUE)

genes <- c("Gata3", "Cxcr3", "Ms4a4a", "Tbx21", "Ifng", "Tnf", "Il2", "Gzma", "Gzmb", "Gzmk", "Prf1", "Ccl5", "Ftla4", "Il10", "Il2ra", "Foxp3", "Klrg1")
genes <- grep("Tgf", row.names(hamster.integrated), value = TRUE)

clusters = clusters_to_check

contrasts_to_use <- c("d2vsd0", "d3vsd0", "d5vsd0", "e14vsd0")
df <- pseudobulk_T %>%
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

ggplot(pseudobulk_T %>%
         dplyr::filter(cluster %in% clusters) %>%
         dplyr::filter(contrast %in% contrasts_to_use) %>%
         dplyr::filter(gene_name %in% genes) %>%
         dplyr::mutate(padj=ifelse(padj < 1E-20, 1E-20, padj)) %>%
         dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
         dplyr::mutate(gene=factor(gene_name,levels=gene.order))%>%
         dplyr::mutate(cluster=factor(cluster,levels=clusters)),
       aes(y=gene,x=cluster,size=-log10(padj),color=log2FoldChange)) +
  geom_point(shape=16) +
  scale_size_continuous(name='adj. p-value',
                        breaks=c(3,6,9),
                        labels=c("1e-3","1e-6", "<1e-9"),
                        limits=c(0,NA)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-6,6),oob=scales::squish) +
  facet_wrap(~contrast, nrow=1) +
  coord_cartesian(clip = 'off') +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
  labs(x='',y='')

ggsave(paste("pseudobulk_tgfs_Tcells", "pdf", sep="."), width=6, height = 4, units="in", useDingbats=FALSE)

ggsave(paste("pseudobulk_Teffectors", "pdf", sep="."), width=6, height = 4, units="in", useDingbats=FALSE)



for (contrasts_to_use in c("d2vsd0", "d3vsd0", "d5vsd0", "e14vsd0")) {
  #for (contrasts_to_use in c("d2vsd0")) {
  
  #Select top genes
  #For the small version in Fig. 2 use n_4 in slice_min and save as pseudobulk_small
  genes <- pseudobulk_T %>%
    dplyr::filter(cluster %in% clusters) %>%
    dplyr::filter(contrast %in% contrasts_to_use) %>%
    dplyr::filter(!grepl('SCoV2',gene_name)) %>%
    dplyr::filter(!(is.na(padj)) & (padj < .01)) %>%
    dplyr::group_by(contrast, cluster) %>%
    dplyr::slice_min(order_by=padj,n=40,with_ties=FALSE) %>%
    dplyr::pull(gene_name)
  
  #cluster values to get proper order for genes (rows) and cell types (columns)
  df <- pseudobulk_T %>%
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
  
  ggplot(pseudobulk_T %>%
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
  ggsave(paste("pseudobulk_T", contrasts_to_use, "pdf", sep="."), width=4, height = h+2, units="in", useDingbats=FALSE)
  
}
