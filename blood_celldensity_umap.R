
library(patchwork)
library(ggplot2)
library(writexl)
library(dplyr)
library(Seurat)
#----------------------------------------------------------------------------------------------------------------------------------------------------

#### CREATE CELL DENSITY UMAP PROJECTION ####

#Load smooth dimplot function
source("./smooth_DimPlot.R")

hamster.integrated <- readRDS("./hamster.integrated_noery_pca_annot_ctype.RDS")


## rename celltypes
hamster.integrated$celltype <- ifelse(hamster.integrated$celltype == "T", "Tcells",
                                      ifelse(hamster.integrated$celltype == "B", "Bcells", 
                                             ifelse(hamster.integrated$celltype == "Classical monocyte", "Classical monocytes",
                                                    ifelse(hamster.integrated$celltype == "Non-classical monocyte","Non-classical monocytes",
                                                           ifelse(hamster.integrated$celltype == "Neutrophil", "Neutrophils",
                                                                  ifelse(hamster.integrated$celltype == "Immature neutrophil", "Immature neutrophils",
                                                                         ifelse(hamster.integrated$celltype == "pDC", "PlasmacytoidDendritic",
                                                                                ifelse(hamster.integrated$celltype == "mDC", "MyeloidDendritic",
                                                                                       ifelse(hamster.integrated$celltype == "Activated T", "Activated Tcells",
                                                                                              ifelse(hamster.integrated$celltype == "Mixed", "Unclear",
                                                                                                     ifelse(hamster.integrated$celltype == "Platelet", "Platelet", 
                                                                                                            ifelse(hamster.integrated$celltype == "NK", "NK",
                                                                                                                   "none"))))))))))))


hamster.integrated <- subset(hamster.integrated, celltype != "Unclear")
Idents(hamster.integrated) <- "celltype"

## CREATE PLOTS

##D2 vs D0

smooth_DimPlot(subset(hamster.integrated, infection %in% c("d0", "d2")),
               group.by='infection',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)+
  scale_color_gradient2(low='blue',mid='gray', high="red" , limits=c(-3,5),oob=scales::squish)

ggsave(filename = "./smoothdim/density_d2vs0_umap.pdf", useDingbats=FALSE )
ggsave(filename = "./smoothdim/density_d2vs0_umap.png" )



##D3 vs D0
smooth_DimPlot(subset(hamster.integrated, infection %in% c("d0", "d3")),
               group.by='infection',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)+
  scale_color_gradient2(low='blue',mid='gray', high="red" , limits=c(-3,5),oob=scales::squish)

ggsave(filename = "./smoothdim/density_d3vs0_umap.pdf", useDingbats=FALSE )
ggsave(filename = "./smoothdim/density_d3vs0_umap.png")



##D5 vs D0
smooth_DimPlot(subset(hamster.integrated, infection %in% c("d0", "d5")),
               group.by='infection',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)+
  scale_color_gradient2(low='blue',mid='gray', high="red" , limits=c(-3,5),oob=scales::squish)

ggsave(filename = "./smoothdim/density_d5vs0_umap.pdf", useDingbats=FALSE)
ggsave(filename = "./smoothdim/density_d5vs0_umap.png" )


##D14 vs D0
smooth_DimPlot(subset(hamster.integrated, infection %in% c("d0", "d14")),
               group.by='infection',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)+
  scale_color_gradient2(low='blue',mid='gray', high="red" , limits=c(-3,5),oob=scales::squish)

ggsave(filename = "./smoothdim/density_d14vs0_umap.pdf", useDingbats=FALSE)
ggsave(filename = "./smoothdim/density_d14vs0_umap.png" )



