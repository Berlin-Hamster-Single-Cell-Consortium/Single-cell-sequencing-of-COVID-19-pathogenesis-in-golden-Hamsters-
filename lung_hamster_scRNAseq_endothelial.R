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
#Subclustering and analysis of endothelial cells (Fig. 6)




#Subclustering
hamster.integrated <- readRDS("./ma_int.rds")
seu_endo_ma <- subset(hamster.integrated, subset = (celltype == "Endothelial"))
DefaultAssay(seu_endo_ma) <- "integrated"                  
seu_endo_ma <- RunPCA(seu_endo_ma, verbose = FALSE)
seu_endo_ma <- RunUMAP(seu_endo_ma, dims = 1:30)
seu_endo_ma <- FindNeighbors(seu_endo_ma, dims = 1:30)
seu_endo_ma <- FindClusters(seu_endo_ma, resolution = 0.6)
UMAPPlot(seu_endo_ma, label=TRUE)
DefaultAssay(seu_endo_ma) <- "SCT"                  



seu_endo_ma@meta.data$UMAP_1 <- NULL
seu_endo_ma@meta.data$UMAP_2 <- NULL
seu_endo_ma@meta.data = cbind(seu_endo_ma@meta.data, seu_endo_ma@reductions$umap@cell.embeddings)


ggplot()+geom_point(data=seu_endo_ma@meta.data, aes(x=UMAP_1, y=UMAP_2, color=seurat_clusters), size=0.5)+theme_bw()
ggplot()+geom_point(data=seu_endo_ma@meta.data, aes(x=UMAP_1, y=UMAP_2, color=infection), size=0.5, alpha=0.4)+theme_bw()
#
SCoV2_rawcounts <- FetchData(seu_endo_ma, grep("SCoV2", seu_endo_ma@assays$RNA@counts@Dimnames[[1]], value="TRUE"), slot="counts")
seu_endo_ma@meta.data$SCoV2_sum <- rowSums(SCoV2_rawcounts)
SCoV2_rawcounts$SCoV2sum <- rowSums(SCoV2_rawcounts)
seu_endo_ma@meta.data$SCoV2_load <- SCoV2_rawcounts$SCoV2sum/seu_endo_ma@meta.data$nCount_RNA*100
seu_endo_ma@meta.data$withvirus = ifelse(seu_endo_ma@meta.data$SCoV2_load>0,"withvirus","novirus")

saveRDS(seu_endo_ma, "./seu_red_endothelial.rds")
seu_endo_ma <- readRDS("./seu_red_endothelial.rds")


################
#For cell type annotation, check the genes from the lung atlas, Travaglini et al. and other endothelial marker genes
#Percent positive for a gene in clusters
for (gene in c("Ccl7", "Nkg7", "Ccl5", "Klre1", "Cldn5", "Plvap", "Nr2f2", "Bmx", "Pecam1", "Ednrb", "Myc", "Hspg2", "Mpzl2", "Vwa1", "Spry1", "Dkk2", "Ackr1", "Lyve1", "Ccl21")) {
  df1 <- cbind.data.frame(FetchData(seu_endo_ma, gene),
                          seu_endo_ma@meta.data$seurat_clusters,
                          str_replace_all(gene, "-", "_"),
                          seu_endo_ma@meta.data$celltype,
                          seu_endo_ma@meta.data$infection)
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
  ggsave(paste("percentpositive", gene, "pdf", sep="."), useDingbats=FALSE)
}


############
#Look for cluster markers
cluster_marks <- FindAllMarkers(seu_endo_ma)
top30 <- cluster_marks %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seu_endo_ma, features = top30$gene) + NoLegend()


DotPlot(seu_endo_ma, group.by='seurat_clusters',
        features = unique(top10$gene, assay='RNA')) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))
ggsave("Endothelial_celltypemarkers.pdf", width=12, height = 8, units="in", useDingbats=FALSE)


#Clusters:
Idents(seu_endo_ma) <- seu_endo_ma@meta.data$seurat_clusters

markers <- c("Ccl7", "Nkg7", "Ccl5", "Klre1", "Cldn5", "Plvap", "Nr2f2", "Bmx", "Pecam1", "Ednrb", "Myc", "Hspg2", "Mpzl2", "Vwa1", "Spry1", "Dkk2", "Ackr1", "Lyve1", "Ccl21")
avg <- AverageExpression(seu_endo_ma, features = markers, return.seurat = T, assays = "SCT") 

DoHeatmap(avg, assay="SCT", size=5, features=markers, draw.lines = F, hjust = 0)+ 
  scale_fill_gradientn(colors = c("blue", "white","red"), na.value = "white")
+
  NoLegend()+
  scale_y_discrete(labels=markers, position = "right" )

DoHeatmap(avg, assay="SCT", size=5, features=markers, label=F, group.bar=F, group.bar.height = 0,
          draw.lines = F, hjust = 0)+ 
  scale_fill_gradientn(colors = c("blue", "white","red"), na.value = "white") +
  NoLegend()+
  theme_bw()+
  scale_y_discrete(labels=markers1, position = "right" )+
  scale_x_discrete( limits = rev(levels(seu_endo_ma$seurat_clusters)))+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, vjust = -2, hjust = -0.15))

ggsave("endo_ctype_markers_hmap.pdf")


#################
#annotate clusters with cell types as following:

seu_endo_ma <- RenameIdents(seu_endo_ma, 
                        '0'='Bronchial',
                        '1'='Capillary',
                        '2'='Bronchial',
                        '3'='Bronchial',
                        '4'='Artery',
                        '5'='Bronchial',
                        '6'='Bronchial',
                        '7'='Vein',
                        '8'='Bronchial',
                        '9'='Capillary',
                        '10'='Unclear1',
                        '11'='Lympathic',
                        '12'='Unclear2',
                        '13'='Unclear3')
seu_endo_ma@meta.data$celltype <- Idents(seu_endo_ma)

cluster_marks <- FindAllMarkers(seu_endo_ma)
write.table(cluster_markers, "Endothelialcells_celltypes_markers.txt", sep="\t", quote=FALSE, row.names = FALSE)
top10 <- cluster_marks %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


DotPlot(seu_endo_ma, group.by='celltype',
        features = unique(top10$gene, assay='RNA')) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))
ggsave("Endothelialcells_celltypemarkers.pdf", width=16, height = 8, units="in", useDingbats=FALSE)


####################
#Plot with cell types
#Endothelial: 0D3B66
celltypecolors = c("Artery" = "#5D3F99",
                   "Bronchial" = "#0090D2",
                   "Capillary" = "#528DCA",
                   "Lympathic" = "#88C9EF",
                   "Unclear1" = "#B2B2B2",
                   "Unclear2" = "#B2B2B2",
                   "Unclear3" = "#B2B2B2",
                   "Vein" = "#2C479E")
means <- seu_endo_ma@meta.data %>%
  group_by(celltype) %>%
  summarise(mean_U1 = mean(UMAP_1),
            mean_U2 = mean(UMAP_2))
ggplot()+
  geom_point(data=seu_endo_ma@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=celltype), size=0.5, shape=16)+
  geom_text(data=subset(means, celltype!="Unclear"), aes(x=mean_U1, y=mean_U2, label=celltype))+
  coord_fixed(ratio=1)+
  ggtitle("celltype")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  scale_colour_manual(values = celltypecolors)
ggsave("Endothelialcells_celltypes.pdf", useDingbats=FALSE)


#Make percentages over time, without unclear cell types
celltypecolors = c("Artery" = "#5D3F99",
                   "Bronchial" = "#0090D2",
                   "Capillary" = "#528DCA",
                   "Lympathic" = "#88C9EF",
                   "Vein" = "#2C479E")
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


a = seu_endo_ma@meta.data %>% group_by(infection, hamster) %>% tally(name="tot") 
b = seu_endo_ma@meta.data %>% group_by(infection, celltype, hamster) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('infection', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  filter(!(celltype %in% c("Unclear1", "Unclear2", "Unclear3"))) %>%
  mutate(fraction = pos / tot)
tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("infection", "celltype")) %>% replace(., is.na(.), 0) %>% replace(., "NaN", 0)

ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, infection, sep="_")))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.9))+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("Endothelial cells subtype proportion", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  theme(legend.position = "none")+
  scale_fill_manual(values=the_colors)
ggsave("Endothelialcells_percentages_celltypes.pdf", useDingbats=FALSE)

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
ggsave("Endothelialcells_percentages_celltypes_d0d5.pdf", useDingbats=FALSE)

#######################


#positive for gene in celltype
for (gene in c("Igfbp5")) {
  df1 <- cbind.data.frame(FetchData(seu_endo_ma, gene), seu_endo_ma@meta.data$seurat_clusters, seu_endo_ma@meta.data$celltype, str_replace_all(gene, "-", "_"), seu_endo_ma@meta.data$SCoV2_load, seu_endo_ma@meta.data$infection, seu_endo_ma@meta.data$hamster)
  colnames(df1) <- c("expression", "cluster", "celltype", "gene", "load", "day", "hamster")
  a = df1 %>% group_by(day, celltype, hamster) %>% tally(name="tot") 
  b = df1 %>% group_by(day, celltype, hamster) %>% filter(expression>0) %>% tally(name="pos") 
  c =
    left_join(a , b , by = c('day', 'celltype', 'hamster')) %>% 
    replace(., is.na(.), 0) %>%
    mutate(fraction = pos / tot)
  tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("day", "celltype"))
  
  ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, day, sep="_")))+
    geom_bar(position=position_dodge(), stat="identity")+
    geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.9))+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    ggtitle(paste("fraction positive", gene, "in virus levels mean and se of three animals", sep=" "))+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))  
  ggsave(paste("Endothelialcells_percposbycelltype", gene,  "pdf", sep="."), useDingbats=FALSE)
}

for (gene in c("Angpt2", "Tek")) {
  df1 <- cbind.data.frame(FetchData(seu_endo_ma, gene), seu_endo_ma@meta.data$seurat_clusters, seu_endo_ma@meta.data$celltype, str_replace_all(gene, "-", "_"), seu_endo_ma@meta.data$SCoV2_load, seu_endo_ma@meta.data$infection, seu_endo_ma@meta.data$hamster)
  colnames(df1) <- c("expression", "cluster", "celltype", "gene", "load", "day", "hamster")

  ggplot()+
    geom_boxplot(data=subset(df1, expression>0), aes(x=celltype, y=expression, fill=day))+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    ggtitle(paste("expression of", gene, sep=" "))+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))#+
  #    scale_fill_manual(values=the_colors_d0d5)
  #scale_fill_manual(values=the_colors)
  
  #  ggsave(paste("Endothelialcells_percposbycelltype_d0d5", gene,  "pdf", sep="."), useDingbats=FALSE)
  ggsave(paste("Endothelialcells_expression", gene,  "pdf", sep="."), useDingbats=FALSE)
}

#########################
#UMAPs

for (gene in c("Ccl7", "Nkg7", "Ccl5", "Klre1", "Cldn5", "Plvap", "Nr2f2", "Bmx", "Pecam1", "Ednrb", "Myc", "Hspg2", "Mpzl2", "Vwa1", "Spry1", "Dkk2", "Ackr1", "Lyve1", "Ccl21")) {
  res <- try(df <- FetchData(seu_endo_ma, gene))
  if(inherits(res, "try-error")) {
    print(paste(gene, " does not exist"))
  }
  else {
    gene <- str_replace_all(gene, "-", "_")
    colnames(df) <- gene
    df = cbind(df, seu_endo_ma@reductions$umap@cell.embeddings)
    ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", shape=16, size=0.5)+
      geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2, colour=eval(parse(text = gene))), stroke=0, size=0.5)+
      scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
      theme_bw()+
      theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
      coord_fixed(ratio=1)+
      ggtitle(gene)
    ggsave(paste(gene, "endo.pdf", sep="."), useDingbats=FALSE)
  }
}




#########################
#pseudobulk
hamster <- seu_endo_ma
Idents(hamster) <- seu_endo_ma@meta.data$celltype

clusters_to_check=c('Bronchial', 'Capillary', 'Artery', 'Vein')
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

#write/read table
write.table(pseudobulk, "./endothelial_pseudobulk.txt", row.names = TRUE, sep="\t", quote=FALSE)
pseudbulk_endo <- read.table("./endothelial_pseudobulk.txt", sep="\t", header=TRUE)


################
#Differential expression dotblot with cytokines or other gene sets

genes <- c("Ccl2", "Ccl3", "Ccl4", "Ccl7", "Ccl8", "Cxcl10", "Cxcl11", "Icam1", "Icam2", "Tnfsf10", "Vcam1")
#ENSMAUG00000020409 is Akap12
genes <- c("Thbs1", "Mmp14", "Hgf", "Vegfa", "Vegfb", "Vegfc", "Vegfd", "Ca4", "Ube2c", "Top2a", "Mki67", "Gadd45a", "Ppp1r15a", "Cdkn1a", "Bbc3", "Nrarp", "ENSMAUG00000020409", "Nr4a1", "Nr2f2")
genes <- grep("Tgf", row.names(hamster.integrated), value = TRUE)
genes <- c("Igfbp5")

clusters = c('Bronchial', 'Capillary', 'Artery', 'Vein')
contrasts_to_use <- c("d2vsd0", "d3vsd0", "d5vsd0", "e14vsd0")

df <- pseudbulk_endo %>%
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

ggplot(pseudbulk_endo %>%
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

ggsave(paste("pseudobulk_tgfs_endo", "pdf", sep="."), width=6.1, height = 4, units="in", useDingbats=FALSE)

ggsave(paste("pseudobulk_someendogenes_endo", "pdf", sep="."), width=6.1, height = 4, units="in", useDingbats=FALSE)



#save 8x12 inches
#ggsave(paste("pseudobulk_cytokines", contrasts_to_use, "pdf", sep="."), width=4, height = 4, units="in", useDingbats=FALSE)
ggsave(paste("pseudobulk_cytokines_endothelial", "pdf", sep="."), width=4, height = 3, units="in", useDingbats=FALSE)

#This is for Figure S5  
for (contrasts_to_use in c("d2vsd0", "d3vsd0", "d5vsd0", "e14vsd0")) {
  #for (contrasts_to_use in c("d2vsd0")) {
  
  #Select top genes
  genes <- pseudbulk_endo %>%
    dplyr::filter(cluster %in% clusters) %>%
    dplyr::filter(contrast %in% contrasts_to_use) %>%
    dplyr::filter(!grepl('SCoV2',gene_name)) %>%
    dplyr::filter(!(is.na(padj)) & (padj < .01)) %>%
    dplyr::group_by(contrast, cluster) %>%
    dplyr::slice_min(order_by=padj,n=50,with_ties=FALSE) %>%
    dplyr::pull(gene_name)
  
  #cluster values to get proper order for genes (rows) and cell types (columns)
  df <- pseudbulk_endo %>%
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
  
  ggplot(pseudbulk_endo %>%
           dplyr::filter(cluster %in% clusters) %>%
           dplyr::filter(contrast %in% contrasts_to_use) %>%
           dplyr::filter(gene_name %in% genes) %>%
           dplyr::mutate(padj=ifelse(padj < 1E-20, 1E-20, padj)) %>%
           dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
           dplyr::mutate(gene=factor(gene_name,levels=gene.order)),
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
  ggsave(paste("pseudobulk_endo", contrasts_to_use, "pdf", sep="."), width=4, height = h+2, units="in", useDingbats=FALSE)
  
}

#For figure 6b

contrasts_to_use <- c("d2vsd0", "d3vsd0")

#Select top genes
genes <- pseudbulk_endo %>%
  dplyr::filter(cluster %in% clusters) %>%
  dplyr::filter(contrast %in% contrasts_to_use) %>%
  dplyr::filter(!grepl('SCoV2',gene_name)) %>%
  dplyr::filter(!(is.na(padj)) & (padj < .01)) %>%
  dplyr::group_by(contrast, cluster) %>%
  dplyr::slice_min(order_by=padj,n=10,with_ties=FALSE) %>%
  dplyr::pull(gene_name)

#cluster values to get proper order for genes (rows) and cell types (columns)
df <- pseudbulk_endo %>%
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

ggplot(pseudbulk_endo %>%
         dplyr::filter(cluster %in% clusters) %>%
         dplyr::filter(contrast %in% contrasts_to_use) %>%
         dplyr::filter(gene_name %in% genes) %>%
         dplyr::mutate(padj=ifelse(padj < 1E-20, 1E-20, padj)) %>%
         dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
         dplyr::mutate(gene=factor(gene_name,levels=gene.order)) %>%
         dplyr::mutate(cluster = paste(cluster, contrast, sep="_")) %>%
         dplyr::mutate(cluster=factor(cluster,levels=c("Bronchial_d2vsd0", "Capillary_d2vsd0", "Artery_d2vsd0", "Vein_d2vsd0", "Bronchial_d3vsd0", "Capillary_d3vsd0", "Artery_d3vsd0", "Vein_d3vsd0"))),
       aes(y=gene,x=cluster,size=-log10(padj),color=log2FoldChange)) +
  geom_point(shape=16) +
  scale_size_continuous(name='adj. p-value',
                        breaks=c(5,10,15),
                        labels=c("1e-5","1e-10", "<1e-15"),
                        limits=c(0,20)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-6,6),oob=scales::squish) +
  coord_cartesian(clip = 'off') +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
  labs(x='',y='')

#save 8x12 inches
h=length(gene.order) / 5
ggsave("pseudobulk_endo_small.pdf", width=6, height = h+1, units="in", useDingbats=FALSE)



  
