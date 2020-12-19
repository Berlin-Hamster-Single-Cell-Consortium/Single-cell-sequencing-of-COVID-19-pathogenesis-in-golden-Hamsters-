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



#For preprocessing:
#First script hamster_merging.R is run to read in and merge the files after mt filtering
#Second scintegrate workflow
#Third Tabula Muris/Travglini cell type annotation (lung_hamster_annotation.Rmd)

hamster.integrated <- readRDS("./lung_hamster_annotated.rds")
hamster.integrated$hamster <- gsub("ma_([de][0-5]{1,2})_lung(_[0-5])","\\1\\2",hamster.integrated@meta.data$orig.ident)
hamster.integrated$infection <- gsub("_[0-5]","",hamster.integrated@meta.data$hamster )
SCoV2_rawcounts <- FetchData(hamster.integrated, grep("SCoV2", hamster.integrated@assays$RNA@counts@Dimnames[[1]], value="TRUE"), slot="counts")
hamster.integrated@meta.data$SCoV2_sum <- rowSums(SCoV2_rawcounts)
SCoV2_rawcounts$SCoV2sum <- rowSums(SCoV2_rawcounts)
hamster.integrated@meta.data$SCoV2_load <- SCoV2_rawcounts$SCoV2sum/hamster.integrated@meta.data$nCount_RNA*100
hamster.integrated@meta.data = cbind(hamster.integrated@meta.data, hamster.integrated@reductions$umap@cell.embeddings)
DimPlot(hamster.integrated, label=F, reduction = "umap", group.by = "hamster")
DimPlot(hamster.integrated, label=F, reduction = "umap", group.by = "infection")
DimPlot(hamster.integrated, label=F, reduction = "umap", group.by = "orig.ident")
DimPlot(hamster.integrated, label=T, reduction = "umap", group.by = "seurat_clusters")
DefaultAssay(hamster.integrated) <- "SCT"
UMAPPlot(object=seu_red, do.label = TRUE)

#Blood object
hamster.integrated <- readRDS("./hamster.integrated_noery_pca_annot_ctype.RDS")


#save/read the working object
hamster.integrated <- readRDS("./ma_int.rds")
saveRDS(hamster.integrated, "./ma_int.rds")

#Random subsampling to 10000 cells for shiny iSEE
seu_red <- subset(hamster.integrated, cells = sample(Cells(hamster.integrated), 10000))

#Basic plots
means <- hamster.integrated@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(mean_U1 = mean(UMAP_1),
            mean_U2 = mean(UMAP_2))
ggplot()+
  geom_point(data=hamster.integrated@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=seurat_clusters), size=0.5, stroke=0)+
  geom_text(data=means, aes(x=mean_U1, y=mean_U2, label=seurat_clusters))+
  coord_fixed(ratio=1)+
  ggtitle("clusters")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("seurat_clusters.pdf", useDingbats=FALSE)

ggplot()+
  geom_point(data=hamster.integrated@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=infection), size=0.5, stroke=0)+
  coord_fixed(ratio=1)+
  ggtitle("clusters")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("timepoints.pdf", useDingbats=FALSE)


means <- hamster.integrated@meta.data %>%
  group_by(predicted.id.tabula_muris) %>%
  summarise(mean_U1 = mean(UMAP_1),
            mean_U2 = mean(UMAP_2))
ggplot()+
  geom_point(data=hamster.integrated@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=predicted.id.tabula_muris), size=0.5, stroke=0)+
  geom_text(data=means, aes(x=mean_U1, y=mean_U2, label=predicted.id.tabula_muris))+
  coord_fixed(ratio=1)+
  ggtitle("predicted.id.tabula_muris")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("predicted_tabulamuris.pdf", useDingbats=FALSE)
means <- hamster.integrated@meta.data %>%
  group_by(predicted.id.travaglini) %>%
  summarise(mean_U1 = mean(UMAP_1),
            mean_U2 = mean(UMAP_2))
ggplot()+
  geom_point(data=hamster.integrated@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=predicted.id.travaglini), size=0.5, stroke=0)+
  geom_text(data=means, aes(x=mean_U1, y=mean_U2, label=predicted.id.travaglini))+
  coord_fixed(ratio=1)+
  ggtitle("predicted.id.tabula_muris")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("predicted_travaglini.pdf", useDingbats=FALSE)



#Think about which cluster is which celltype, this is the conclusion:
Idents(hamster.integrated) <- hamster.integrated@meta.data$seurat_clusters
hamster.integrated <- RenameIdents(hamster.integrated, 
                        '0'='Tcells',
                        '1'='AlveolarMacrophages',
                        '2'='Bcells',
                        '3'='AT2',
                        '4'='AlveolarMacrophages',
                        '5'='MonocyticMacrophages',
                        '6'='Endothelial',
                        '7'='Tcells',
                        '8'='MonocyticMacrophages',
                        '9'='MacrophagesTreml4+',
                        '10'='MonocyticMacrophages',
                        '11'='Tcells',
                        '12'='MonocyticMacrophages',
                        '13'='InterstitialMacrophages',
                        '14'='Fibroblasts',
                        '15'='Tcells',
                        '16'='MyeloidDendritic',
                        '17'='Bcells',
                        '18'='AT1',
                        '19'='Endothelial',
                        '20'='Endothelial',
                        '21'='Tcells',
                        '22'='SmoothMuscle',
                        '23'='MonocyticMacrophages',
                        '24'='Neutrophils',
                        '25'='PlasmacytoidDendritic',
                        '26'='AlveolarMacrophages',
                        '27'='Neutrophils',
                        '28'='Tcells',
                        '29'='Ciliated',
                        '30'='Unclear',
                        '31'='Endothelial',
                        '32'='Endothelial',
                        '33'='InterstitialMacrophages',
                        '34'='Myofibroblast',
                        '35'='Endothelial',
                        '36'='Endothelial',
                        '37'='Unclear',
                        '38'='AlveolarMacrophages')
hamster.integrated@meta.data$celltype <- Idents(hamster.integrated)


celltypecolors = c(
  "Tcells" = "#368F8B",
  "MonocyticMacrophages" = "#B7245C",
  "Endothelial" = "#0D3B66",
  "AlveolarMacrophages" = "#DFACC4",
  "AT2" = "#F97E44",
  "Bcells" = "#62C370",
  "MacrophagesTreml4+" = "#3E2F5B",
  "Fibroblasts" = "#B2675E",
  "MyeloidDendritic" = "#4F6D7A",
  "AT1" = "#F7C548",
  "Neutrophils" = "#0081AF",
  "SmoothMuscle" = "#644536",
  "InterstitialMacrophages" = "#B97C9D",
  "PlasmacytoidDendritic" = "#7C6A0A",
  "Ciliated" = "#FB3640",
  "Unclear" = "#CAD2C5",
  "Myofibroblast" = "#C4A381"
)


means <- hamster.integrated@meta.data %>%
  group_by(celltype) %>%
  summarise(mean_U1 = mean(UMAP_1),
            mean_U2 = mean(UMAP_2))
ggplot()+
  geom_point(data=hamster.integrated@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=celltype), size=0.5, shape=16)+
  geom_text(data=subset(means, celltype!="Unclear"), aes(x=mean_U1, y=mean_U2, label=celltype))+
  coord_fixed(ratio=1)+
  ggtitle("celltype")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  scale_colour_manual(values = celltypecolors)
ggsave("celltypes_2.pdf", useDingbats=FALSE)


#Order cell types by abundance as in Fig. 1
the_celltypes = c("AlveolarMacrophages",
                  "InterstitialMacrophages",
                  "MonocyticMacrophages",
                  "MacrophagesTreml4+",
                  "Neutrophils",
                  "MyeloidDendritic",
                  "PlasmacytoidDendritic",
                  "Tcells",
                  "Bcells",
                  "AT1",
                  "AT2",
                  "Fibroblasts",
                  "Ciliated",
                  "Endothelial",
                  "Myofibroblast",
                  "SmoothMuscle",
                  "Unclear")


############
#Coloring schemes
celltypecolors <- as.data.frame(celltypecolors)
celltypecolors <- celltypecolors[order(row.names(celltypecolors)),]
bright_colors = as.vector(celltypecolors)

#Make colors darker with time
expr <- list()
for (cbright in bright_colors) {
  r=(col2rgb(cbright)-40)[[1]]
  r=ifelse(r<0,0,r)
  g=(col2rgb(cbright)-40)[[2]]
  g=ifelse(g<0,0,g)
  b=(col2rgb(cbright)-40)[[3]]
  b=ifelse(b<0,0,b)
  cdark = rgb(r, g, b, maxColorValue = 255)
  expr[[cbright]] <- scales::seq_gradient_pal(cbright, cdark, "Lab")(seq(0,1,length.out=5))
}

#Make colors brighter with time
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


#############
#Barplot with percentages


a = hamster.integrated@meta.data %>% group_by(infection, hamster) %>% tally(name="tot") 
b = hamster.integrated@meta.data %>% group_by(infection, celltype, hamster) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('infection', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("infection", "celltype"))
tgc = tgc %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))


ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, infection, sep="_")))+
  geom_bar(position=position_dodge(), stat="identity")+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.9))+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("celltype proportions", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  theme(legend.position = "none")+
  scale_fill_manual(values=the_colors)



#Make Plot with cluster markers
cluster_markers <- FindAllMarkers(hamster.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(cluster_markers, "lung_celltypes_markers.txt", sep="\t", quote=FALSE, row.names = FALSE)

#To use only d0 for this, do
cluster_markers_d0 <- FindAllMarkers(subset(hamster.integrated, infection == "d0"), only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top3 <- cluster_markers_d0 %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)


top3 <- cluster_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
DotPlot(hamster.integrated, group.by='celltype',
        features = unique(top3$gene, assay='RNA')) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))
ggsave("lung_celltypemarkers_onlyd0.pdf", width=12, height = 8, units="in", useDingbats=FALSE)


#Percent positive for a gene in clusters
gene <- "Ifng"
df1 <- cbind.data.frame(FetchData(hamster.integrated, gene),
                        hamster.integrated@meta.data$seurat_clusters,
                        str_replace_all(gene, "-", "_"),
                        hamster.integrated@meta.data$celltype)
colnames(df1) <- c("expression", "cluster", "gene", "celltype")
a = df1 %>% group_by(cluster) %>% tally(name="tot") 
b = df1 %>% filter(expression>0) %>% group_by(cluster) %>% tally(name="pos") 
barplotvalues =
  left_join(a , b , by = 'cluster') %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
ggplot()+
  geom_bar(data=barplotvalues, aes(x=cluster, y=fraction), stat="identity")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("percent positiv", gene, "in", "clusters", sep=" "))

#Percent positive for a gene in celltypes
#positive for gene in celltype
for (gene in c("Cxcl10")) {
  df1 <- cbind.data.frame(FetchData(hamster.integrated, gene), hamster.integrated@meta.data$seurat_clusters, hamster.integrated@meta.data$celltype, str_replace_all(gene, "-", "_"), hamster.integrated@meta.data$SCoV2_load, hamster.integrated@meta.data$infection, hamster.integrated@meta.data$hamster)
  colnames(df1) <- c("expression", "cluster", "celltype", "gene", "load", "day", "hamster")
  a = df1 %>% group_by(day, celltype, hamster) %>% tally(name="tot")
  #b = df1 %>% group_by(day, celltype, hamster) %>% filter(load>0) %>% tally(name="pos") 
  b = df1 %>% group_by(day, celltype, hamster) %>% filter(expression>0) %>% tally(name="pos") 
  c =
    left_join(a , b , by = c('day', 'celltype', 'hamster')) %>% 
    replace(., is.na(.), 0) %>%
    mutate(fraction = pos / tot)
  write.table(c, paste("cells", gene, "csv", sep="."), sep=",", quote=FALSE, row.names = FALSE)
  #%>%
  #  mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))
  tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("day", "celltype"))
  
  ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, day, sep="_")))+
    geom_bar(position=position_dodge(), stat="identity")+
    geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.9))+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    ggtitle(paste("fraction positive", gene, "in cell types mean and se of three animals", sep=" "))+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
    #scale_fill_manual(values=the_colors)+
    theme(legend.position = "none")

  
  #  ggsave(paste("Endothelialcells_percposbycelltype_d0d5", gene,  "pdf", sep="."), useDingbats=FALSE)
  ggsave(paste("celltypes_percposbycelltype", gene,  "pdf", sep="."), useDingbats=FALSE)
  ggplot()+
    geom_boxplot(data=subset(df1, expression>0), aes(x=celltype, y=expression, fill=day))+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    ggtitle(paste("expression gt0 of ", gene, "in cell types", sep=" "))+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
    #scale_fill_manual(values=the_colors)+
    theme(legend.position = "none")
  ggsave(paste("celltypes_expression_gt0", gene,  "pdf", sep="."), useDingbats=FALSE)
}

#############################
#Statistics


#Make tables with numbers
expr <- list()
for (animal in unique(hamster.integrated@meta.data$orig.ident)) {
  expr[[animal]] = hamster.integrated@meta.data %>% filter(orig.ident == animal) %>% group_by(celltype) %>% tally()
  expr[[animal]] = expr[[animal]] %>% add_row(celltype="all", n=sum(expr[[animal]]$n)) %>% add_column(hamster=animal)
}
values <- reshape2::melt(do.call(rbind, expr), id.vars=c("hamster", "celltype"))
values$variable <- NULL
values <- reshape2::dcast(values, celltype ~ hamster)
write.table(values, "hamster_celltypes.txt", sep="\t", quote=FALSE, row.names = FALSE)

#Same with adding plus/minus virus information DOES NOT WORK YET
a <- hamster.integrated@meta.data %>% group_by(orig.ident, celltype) %>% tally() %>% mutate(key=paste(orig.ident, celltype, sep="___"))
b <- hamster.integrated@meta.data %>% filter(SCoV2_load > 0) %>% group_by(orig.ident, celltype) %>% tally(name="pos")  %>% mutate(key=paste(orig.ident, celltype, sep="___"))
barplotvalues = left_join(a , b , by = 'key') %>% replace(., is.na(.), 0) %>% mutate(fraction = pos / n)

a = hamster.integrated@meta.data %>% group_by(infection, celltype, hamster) %>% tally(name="tot") 
b = hamster.integrated@meta.data %>% group_by(infection, celltype, hamster) %>% filter(SCoV2_load>0) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('infection', 'celltype', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
write.table(c, "celltypes_with_virus.csv", sep=",", quote=FALSE, row.names = FALSE)

a <- hamster.integrated@meta.data %>% group_by(infection, celltype) %>% tally() %>% mutate(key=paste(infection, celltype, sep="___"))
a$infection <- NULL
a$celltype <- NULL
b <- hamster.integrated@meta.data %>% filter(SCoV2_load > 0) %>% group_by(infection, celltype) %>% tally(name="pos")  %>% mutate(key=paste(infection, celltype, sep="___"))
b$infection <- NULL
b$celltype <- NULL
barplotvalues = left_join(a , b , by = 'key') %>% replace(., is.na(.), 0) %>% mutate(fraction = pos / n) %>% separate(key, c("day", "celltype"), sep="___")

ggplot()+
  geom_bar(data=barplotvalues, aes(x=celltype, y=fraction, fill=day), position=position_dodge2(), stat="identity")+
  theme_bw()+
  ggtitle(paste("percent with virus"))+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))

ggsave("celltype_fraction_with_virus.pdf", useDingbats=FALSE)

############################
#smooth dimplots
smooth_DimPlot(subset(hamster.integrated, infection %in% c("d0", "d2")),
               group.by='infection',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave("density_d2_vs_d0.pdf", useDingbats=FALSE)
smooth_DimPlot(subset(hamster.integrated, infection %in% c("d0", "d3")),
               group.by='infection',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave("density_d3_vs_d0.pdf", useDingbats=FALSE)
smooth_DimPlot(subset(hamster.integrated, infection %in% c("d0", "d5")),
               group.by='infection',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave("density_d5_vs_d0.pdf", useDingbats=FALSE)
smooth_DimPlot(subset(hamster.integrated, infection %in% c("d0", "e14")),
               group.by='infection',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave("density_d14_vs_d0.pdf", useDingbats=FALSE)
smooth_DimPlot(subset(hamster.integrated, infection %in% c("d2", "d3")),
               group.by='infection',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave("density_d3_vs_d2.pdf", useDingbats=FALSE)



######################
#UMAP SCoV2 load
ggplot()+geom_point(data=subset(hamster.integrated@meta.data, SCoV2_load<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", stroke=0, size=0.5)+
  geom_point(data=subset(hamster.integrated@meta.data, SCoV2_load>0), aes(x=UMAP_1, y=UMAP_2, colour=log10(SCoV2_load)), stroke=0, size=0.5)+
  scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  coord_fixed(ratio=1)+
  ggtitle("log10 SCoV2 load")
ggsave(paste("log10SCoV2_load", "pdf", sep="."), useDingbats=FALSE)


#Density plot
the_celltypes = c("AT2", "Monocytes")
ggplot()+
  geom_density(data=subset(hamster.integrated@meta.data, (SCoV2_load>0) & (celltype %in% the_celltypes) &(infection =="d5")), aes(x=log10(SCoV2_load), colour=celltype))+
  theme_bw()

#Boxplot viral load
ggplot()+
  geom_boxplot(data=subset(hamster.integrated@meta.data, SCoV2_load>=0 & infection %in% c("d2", "d3", "d5")), aes(x=celltype, y=log10(SCoV2_load), fill=infection))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))

#############################
#UMAPs for genes

for (gene in c("Cxcl10")) {
res <- try(df <- FetchData(hamster.integrated, gene))
  if(inherits(res, "try-error")) {
    print(paste(gene, " does not exist"))
  }
  else {
    gene <- str_replace_all(gene, "-", "_")
    colnames(df) <- gene
    df = cbind(df, hamster.integrated@reductions$umap@cell.embeddings)
    ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", shape=16, size=0.5)+
      geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2, colour=eval(parse(text = gene))), stroke=0, size=0.5)+
      scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
      theme_bw()+
      theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
      coord_fixed(ratio=1)+
      ggtitle(gene)
    ggsave(paste(gene, "pdf", sep="."), useDingbats=FALSE)
  }
}


###############################
#UMAPs with two genes

gene <- "Log10_SCoV2_percent"
secondgene <- "Ace2"
df <- cbind.data.frame(hamster.integrated@meta.data$SCoV2_load, FetchData(hamster.integrated, secondgene))
colnames(df) <- c(gene, secondgene)
df = cbind(df, hamster.integrated@reductions$umap@cell.embeddings)
ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", size=0.5, shape = 16)+
  geom_point(data=subset(df, eval(parse(text = secondgene)) > 0), aes(x=UMAP_1, y=UMAP_2), colour="deepskyblue3", size=0.5, shape = 16)+
  geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2, colour=log10(eval(parse(text = gene)))), alpha=0.8, size=0.5, shape = 16)+
  scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  coord_fixed(ratio=1)+
  ggtitle(paste(gene, secondgene, sep=" and "))
ggsave(paste(gene, secondgene, "pdf", sep="."), useDingbats=FALSE)


gene <- "Tmprss2"
secondgene <- "Ace2"
df <- cbind.data.frame(FetchData(hamster.integrated, gene), FetchData(hamster.integrated, secondgene))
colnames(df) <- c(gene, secondgene)
df = cbind(df, hamster.integrated@reductions$umap@cell.embeddings)
ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0 & eval(parse(text = secondgene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", stroke=0, size=0.5, shape = 16)+
  geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2), alpha=0.5, colour="tomato1", stroke=0, size=0.8, shape = 16)+
  geom_point(data=subset(df, eval(parse(text = secondgene)) > 0), aes(x=UMAP_1, y=UMAP_2), alpha=0.7, colour="deepskyblue3", stroke=0, size=0.8, shape = 16)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  coord_fixed(ratio=1)+
  ggtitle(paste(gene, secondgene, sep=" and "))
ggsave(paste(gene, secondgene, "pdf", sep="."), useDingbats=FALSE)



