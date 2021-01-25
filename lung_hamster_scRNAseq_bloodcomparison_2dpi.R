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
library(ggpatern)


#For preprocessing:
#First script hamster_merging.R is run to read in and merge the files after mt filtering
#Second Dylan runs the integration

hamster.lungblood.d2 <- readRDS("./seu_d2_combined_mtfilt.integrated.rds")
DefaultAssay(hamster.lungblood.d2) <- "integrated"
hamster.lungblood.d2 <- FindNeighbors(hamster.lungblood.d2,verbose=FALSE, assay='integrated') %>%
  FindClusters(resolution=1.0, verbose=FALSE, assay='integrated')
hamster.lungblood.d2$hamster <- gsub("ma_([de][0-5]{1,2})_.*(_[0-5])","\\1\\2",hamster.lungblood.d2@meta.data$orig.ident)
hamster.lungblood.d2$infection <- gsub("_[0-5]","",hamster.lungblood.d2@meta.data$hamster)
hamster.lungblood.d2$tissue <- gsub("ma_[de][0-5]{1,2}_(.*)_[0-5]","\\1",hamster.lungblood.d2@meta.data$orig.ident)
DefaultAssay(hamster.lungblood.d2) <- "RNA"
SCoV2_rawcounts <- FetchData(hamster.lungblood.d2, grep("SCoV2", hamster.lungblood.d2@assays$RNA@counts@Dimnames[[1]], value="TRUE"), slot="counts")
DefaultAssay(hamster.lungblood.d2) <- "SCT"
hamster.lungblood.d2@meta.data$SCoV2_sum <- rowSums(SCoV2_rawcounts)
SCoV2_rawcounts$SCoV2sum <- rowSums(SCoV2_rawcounts)
hamster.lungblood.d2@meta.data$SCoV2_load <- SCoV2_rawcounts$SCoV2sum/hamster.lungblood.d2@meta.data$nCount_RNA*100
hamster.lungblood.d2@meta.data = cbind(hamster.lungblood.d2@meta.data, hamster.lungblood.d2@reductions$umap@cell.embeddings)
DimPlot(hamster.lungblood.d2, label=T, reduction = "umap", group.by = "seurat_clusters")


#Verify cluster / cell type identities
means <- hamster.lungblood.d2@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(mean_U1 = mean(UMAP_1),
            mean_U2 = mean(UMAP_2))
ggplot()+
  geom_point(data=hamster.lungblood.d2@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=seurat_clusters), size=0.5, shape=16)+
  geom_text(data=means, aes(x=mean_U1, y=mean_U2, label=seurat_clusters))+
  coord_fixed(ratio=1)+
  ggtitle("clusters")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("lungblood.d2.seurat_clusters.pdf", useDingbats=FALSE)
ggplot()+
  geom_point(data=hamster.lungblood.d2@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=tissue), size=0.5, shape=16)+
  coord_fixed(ratio=1)+
  ggtitle("clusters")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("lungblood.d2.tissue.pdf", useDingbats=FALSE)
for (gene in c("C1qb", "Marco", "Treml4", "Cd14", "Cd68")) {
  res <- try(df <- FetchData(hamster.lungblood.d2, gene))
  if(inherits(res, "try-error")) {
    print(paste(gene, " does not exist"))
  }
  else {
    gene <- str_replace_all(gene, "-", "_")
    colnames(df) <- gene
    df = cbind(df, hamster.lungblood.d2@reductions$umap@cell.embeddings)
    ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", stroke=0, size=0.5)+
      geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2, colour=eval(parse(text = gene))), shape=16, size=0.5)+
      scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
      theme_bw()+
      theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
      coord_fixed(ratio=1)+
      ggtitle(gene)
    ggsave(paste(gene, "d2", "pdf", sep="."), useDingbats=FALSE)
  }
}  

#############################
#Same with d0
hamster.lungblood.d0 <- readRDS("./seu_d0_combined_mtfilt.integrated.rds")
DefaultAssay(hamster.lungblood.d0) <- "integrated"
hamster.lungblood.d0 <- FindNeighbors(hamster.lungblood.d0,verbose=FALSE, assay='integrated') %>%
  FindClusters(resolution=1.0, verbose=FALSE, assay='integrated')
hamster.lungblood.d0$hamster <- gsub("ma_([de][0-5]{1,2})_.*(_[0-5])","\\1\\2",hamster.lungblood.d0@meta.data$orig.ident)
hamster.lungblood.d0$infection <- gsub("_[0-5]","",hamster.lungblood.d0@meta.data$hamster)
hamster.lungblood.d0$tissue <- gsub("ma_[de][0-5]{1,2}_(.*)_[0-5]","\\1",hamster.lungblood.d0@meta.data$orig.ident)
DefaultAssay(hamster.lungblood.d0) <- "RNA"
SCoV2_rawcounts <- FetchData(hamster.lungblood.d0, grep("SCoV2", hamster.lungblood.d0@assays$RNA@counts@Dimnames[[1]], value="TRUE"), slot="counts")
DefaultAssay(hamster.lungblood.d0) <- "SCT"
hamster.lungblood.d0@meta.data$SCoV2_sum <- rowSums(SCoV2_rawcounts)
SCoV2_rawcounts$SCoV2sum <- rowSums(SCoV2_rawcounts)
hamster.lungblood.d0@meta.data$SCoV2_load <- SCoV2_rawcounts$SCoV2sum/hamster.lungblood.d0@meta.data$nCount_RNA*100
hamster.lungblood.d0@meta.data = cbind(hamster.lungblood.d0@meta.data, hamster.lungblood.d0@reductions$umap@cell.embeddings)
DimPlot(hamster.lungblood.d0, label=T, reduction = "umap", group.by = "seurat_clusters")

#Verify cluster / cell type identities
means <- hamster.lungblood.d0@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(mean_U1 = mean(UMAP_1),
            mean_U2 = mean(UMAP_2))
ggplot()+
  geom_point(data=hamster.lungblood.d0@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=seurat_clusters), size=0.5, shape=16)+
  geom_text(data=means, aes(x=mean_U1, y=mean_U2, label=seurat_clusters))+
  coord_fixed(ratio=1)+
  ggtitle("clusters")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("lungblood.d0.seurat_clusters.pdf", useDingbats=FALSE)
ggplot()+
  geom_point(data=hamster.lungblood.d0@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=tissue), size=0.5, shape=16)+
  coord_fixed(ratio=1)+
  ggtitle("clusters")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("lungblood.d0.tissue.pdf", useDingbats=FALSE)

for (gene in c("C1qb", "Marco", "Treml4", "Cd14", "Cd68")) {
  res <- try(df <- FetchData(hamster.lungblood.d0, gene))
  if(inherits(res, "try-error")) {
    print(paste(gene, " does not exist"))
  }
  else {
    gene <- str_replace_all(gene, "-", "_")
    colnames(df) <- gene
    df = cbind(df, hamster.lungblood.d0@reductions$umap@cell.embeddings)
    ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", stroke=0, size=0.5)+
      geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2, colour=eval(parse(text = gene))), shape=16, size=0.5)+
      scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
      theme_bw()+
      theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
      coord_fixed(ratio=1)+
      ggtitle(gene)
    ggsave(paste(gene, "d0", "pdf", sep="."), useDingbats=FALSE)
  }
}  



###############
#Prepare data for scatter plot in Figure S5C

#values for macrophages without interstitial macrophages but with Treml4+ Macrophages
#celltype_clusters <- c(3,5,6,7,14,20,24)
celltype_clusters <- c(0,4,6,7,17)
celltype <- "MonocytesMacrophages"
t <- 0.4
cond1 <- paste(celltype, "blood", sep="_")
cond2 <- paste(celltype, "lung", sep="_")
plottitle <- paste("average expression d2", celltype, "cells", sep=" ")
filename <- paste("diff_d2", celltype, cond1, cond2, "txt", sep=".")

hamster.lungblood.d2@meta.data$cell_id <- row.names(hamster.lungblood.d2@meta.data)
hamster.lungblood.d2 <- SetIdent(object = hamster.lungblood.d2, value=hamster.lungblood.d2@meta.data$orig.ident)
hamster.lungblood.d2 <- SetIdent(object = hamster.lungblood.d2, cells = subset(hamster.lungblood.d2@meta.data, (seurat_clusters %in% celltype_clusters) & tissue == "blood")$cell_id, value = cond1)
hamster.lungblood.d2 <- SetIdent(object = hamster.lungblood.d2, cells = subset(hamster.lungblood.d2@meta.data, (seurat_clusters %in% celltype_clusters) & tissue == "lung")$cell_id, value = cond2)
avg.cluster_x <- log1p(AverageExpression(subset(hamster.lungblood.d2, idents = cond1), verbose = FALSE)$SCT)
avg.cluster_y <- log1p(AverageExpression(subset(hamster.lungblood.d2, idents = cond2), verbose = FALSE)$SCT)
avg.cluster_x$gene <- rownames(avg.cluster_x)
df <- cbind(avg.cluster_x, avg.cluster_y)
df$diff <- abs(df[[cond1]] - df[[cond2]])
ggplot()+geom_point(data=df, aes(x=eval(parse(text = cond1)), y=eval(parse(text = cond2))), size=1, shape=16)+
  ggtitle(plottitle)+
  ggrepel::geom_text_repel(data=subset(df, diff>t), aes(x=eval(parse(text = cond1)), y=eval(parse(text = cond2)), label=gene), size=2, colour="red")+
  ggrepel::geom_text_repel(data=subset(df, gene %in% c("Isg15", "Ifit2")), aes(x=eval(parse(text = cond1)), y=eval(parse(text = cond2)), label=gene), size=2, colour="red")+
  theme_bw()+
  xlab(cond1)+
  ylab(cond2)

diffgenes.d2 <- FindMarkers(hamster.lungblood.d2, ident.1= cond2, ident.2= cond1, min.pct=0.001, min.diff.pct=0, logfc.threshold=0)
diffgenes.d2$gene <- row.names(diffgenes.d2)
write.table(diffgenes, filename, sep="\t", quote=FALSE, row.names = FALSE)

#values for macrophages without interstitial macrophages but with Treml4+ Macrophages
#celltype_clusters <- c(4, 5, 6, 8, 26)
celltype_clusters <- c(4, 5, 6, 8, 10, 19)
celltype <- "MonocytesMacrophages"
t <- 0.4
cond1 <- paste(celltype, "blood", sep="_")
cond2 <- paste(celltype, "lung", sep="_")
plottitle <- paste("average expression d0", celltype, "cells", sep=" ")
filename <- paste("diff_d0", celltype, cond1, cond2, "txt", sep=".")

hamster.lungblood.d0@meta.data$cell_id <- row.names(hamster.lungblood.d0@meta.data)
hamster.lungblood.d0 <- SetIdent(object = hamster.lungblood.d0, value=hamster.lungblood.d0@meta.data$orig.ident)
hamster.lungblood.d0 <- SetIdent(object = hamster.lungblood.d0, cells = subset(hamster.lungblood.d0@meta.data, (seurat_clusters %in% celltype_clusters) & tissue == "blood")$cell_id, value = cond1)
hamster.lungblood.d0 <- SetIdent(object = hamster.lungblood.d0, cells = subset(hamster.lungblood.d0@meta.data, (seurat_clusters %in% celltype_clusters) & tissue == "lung")$cell_id, value = cond2)
avg.cluster_x <- log1p(AverageExpression(subset(hamster.lungblood.d0, idents = cond1), verbose = FALSE)$SCT)
avg.cluster_y <- log1p(AverageExpression(subset(hamster.lungblood.d0, idents = cond2), verbose = FALSE)$SCT)
avg.cluster_x$gene <- rownames(avg.cluster_x)
df <- cbind(avg.cluster_x, avg.cluster_y)
df$diff <- abs(df[[cond1]] - df[[cond2]])
ggplot()+geom_point(data=df, aes(x=eval(parse(text = cond1)), y=eval(parse(text = cond2))), size=1, shape=16)+
  ggtitle(plottitle)+
  ggrepel::geom_text_repel(data=subset(df, diff>t), aes(x=eval(parse(text = cond1)), y=eval(parse(text = cond2)), label=gene), size=2, colour="red")+
  ggrepel::geom_text_repel(data=subset(df, gene %in% c("Isg15", "Ifit2")), aes(x=eval(parse(text = cond1)), y=eval(parse(text = cond2)), label=gene), size=2, colour="red")+
  theme_bw()+
  xlab(cond1)+
  ylab(cond2)

diffgenes.d0 <- FindMarkers(hamster.lungblood.d0, ident.1= cond2, ident.2= cond1, min.pct=0.001, min.diff.pct=0, logfc.threshold=0)
diffgenes.d0$gene <- row.names(diffgenes.d0)
write.table(diffgenes.d0, filename, sep="\t", quote=FALSE, row.names = FALSE)

diffgenes <- merge(diffgenes.d0, diffgenes.d2, by="gene")

ggplot()+
  geom_point(data=diffgenes, aes(x=avg_logFC.x, y=avg_logFC.y), size=1.5, shape=16)+
  ggrepel::geom_text_repel(data=subset(diffgenes, (p_val_adj.x>0.01 & p_val_adj.y<0.01) & (abs(avg_logFC.x) > 0.4 | abs(avg_logFC.y) > 0.4)), aes(x=avg_logFC.x, y=avg_logFC.y, label=gene), colour="red", min.segment.length = unit(0, 'lines'))+
  ggrepel::geom_text_repel(data=subset(diffgenes, (p_val_adj.x<0.01 & p_val_adj.y>0.01) & (abs(avg_logFC.x) > 0.4 | abs(avg_logFC.y) > 0.4)), aes(x=avg_logFC.x, y=avg_logFC.y, label=gene), colour="green", min.segment.length = unit(0, 'lines'))+
  xlim(-1.5, 1.5)+
  ylim(-1.5, 1.5)+
  coord_fixed(ratio=1)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("MonoMacro blood/lung")+
  xlab("lung/blood d0")+
  ylab("lung/blood d2")
  
ggsave("blood_lung_d0_d2.pdf", useDingbats=FALSE)

###################
#Boxplot with Cxcl10 for Fig. 4c

#d2
celltype_clusters_d2 <- c(0,4,6,7,17)
celltype_clusters_d0 <- c(4, 5, 6, 8, 10, 19)
celltype <- "MonocytesMacrophages"
the_colors = c("#cc658c", "#d37b9d", "#b7245c", "#be396c")
the_colors = c("#b7245c", "#d37b9d")

for (gene in c("Cxcl10")) {
  df0 <- cbind.data.frame(FetchData(hamster.lungblood.d0, gene), hamster.lungblood.d0@meta.data$seurat_clusters, str_replace_all(gene, "-", "_"), hamster.lungblood.d0@meta.data$infection, hamster.lungblood.d0@meta.data$hamster, hamster.lungblood.d0@meta.data$tissue)
  colnames(df0) <- c("expression", "cluster", "gene", "day", "hamster", "tissue")
  df0 <- subset(df0, cluster %in% celltype_clusters_d0)
  df2 <- cbind.data.frame(FetchData(hamster.lungblood.d2, gene), hamster.lungblood.d2@meta.data$seurat_clusters, str_replace_all(gene, "-", "_"), hamster.lungblood.d2@meta.data$infection, hamster.lungblood.d2@meta.data$hamster, hamster.lungblood.d2@meta.data$tissue)
  colnames(df2) <- c("expression", "cluster", "gene", "day", "hamster", "tissue")
  df2 <- subset(df2, cluster %in% celltype_clusters_d2)
  df <- rbind.data.frame(df0, df2)

  a = df %>% group_by(day, tissue, hamster) %>% tally(name="tot") 
  b = df %>% group_by(day, tissue, hamster) %>% filter(expression>0) %>% tally(name="pos") 
  c =
    left_join(a , b , by = c('day', 'tissue', 'hamster')) %>% 
    replace(., is.na(.), 0) %>%
    mutate(fraction = pos / tot)
  tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("day", "tissue"))

  expr <- list()
  #p-value calculation
  for (d in c("d0", "d2")) {
    df2 <- subset(df, day == d)
    df2$expr <- ifelse(df2$expression>0, TRUE, FALSE)
    df2$load <- NULL
    df2$expression <- NULL
    df2$cluster <- NULL
    df2$day <- NULL
    df2$gene <- NULL
    res <- try(p <- summary(glmer(df2$expr ~ df2$tissue + (1 | df2$hamster), family=binomial))$coefficients[2,4])
    if(inherits(res, "try-error")) {
      p <- "NA"
    }
    expr[[d]] <- p
  }
  pvalues <- as.data.frame(do.call(rbind, expr))
  write.table(pvalues, paste("./pvalues_bloodlungd0d2_barplots", gene, "txt", sep="."), quote=FALSE, sep="\t", col.names=FALSE)
  
  expr <- list()
  #p-value calculation
  for (d in c("d0", "d2")) {
    df2 <- subset(df, day == d)
    res <- try(p <- wilcox.test(x=subset(df2, tissue == "blood")$expression, y=subset(df2, tissue == "lung")$expression)$p.value)
    if(inherits(res, "try-error")) {
      p <- "NA"
    }
    expr[[paste(d)]] <- p
  }
  pvalues <- as.data.frame(do.call(rbind, expr))
  write.table(pvalues, paste("./pvalues_bloodlungd0d2_boxplots", gene, "txt", sep="."), quote=FALSE, sep="\t", col.names=FALSE)
  
   
  p0 <- ggplot(tgc, aes(x=tissue, y=fraction, fill=day))+
    geom_bar(position=position_dodge(.8), stat="identity", width = 0.7)+
    geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.8))+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    ggtitle(paste("fraction positive", gene, "in virus levels mean and sd of three animals", sep=" "))+
    scale_fill_manual(values=the_colors)
  
  p1 <- ggplot(tgc, aes(x=tissue, y=fraction, fill=day, pattern_angle=tissue))+
    ggpattern::geom_bar_pattern(position=position_dodge(.8), width = 0.7, stat="identity", pattern='placeholder', pattern_type='kitten')+
    geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.8))+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    ggtitle(paste("fraction positive", gene, "in virus levels mean and sd of three animals", sep=" "))+
    scale_fill_manual(values=the_colors)
  
  p2 <- ggplot()+
    geom_boxplot(data=subset(df, expression>0), aes(x=tissue, y=expression, fill=day), outlier.size=0.2, position=position_dodge(.8), width = 0.7, lwd=0.25, fatten=4)+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    ggtitle(paste("boxplot_quintiles", gene,  sep=" "))+
    ylim(0, NA)+
    scale_fill_manual(values=the_colors)
  pdf(file=paste(gene, "percentpositive", "boxplotgt0", "bloodlung_v2", "pdf", sep="."), width=6, height=8, useDingbats=FALSE)
  grid.arrange(p0, p2, ncol=1)
  dev.off()
}
