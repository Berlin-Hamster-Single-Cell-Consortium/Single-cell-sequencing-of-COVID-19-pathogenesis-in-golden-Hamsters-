library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)
library(ggplot2)
library(ggpubr)
#------------------------------------------------------------------------------------------

#### HEATMAP MARKER GENES ####

Idents(hamster.integrated) <- "seurat_clusters"
avg <- AverageExpression(hamster.integrated, features = markers, return.seurat = T, assays = "SCT") 
avg

## determine average expression of marker genes per cluster
avg <- AverageExpression(hamster.integrated, features = markers, return.seurat = T, assays = "SCT") 
avg

## plot heatmap
DoHeatmap(avg, assay="SCT", size=5, features=markers, label=F, group.bar=F, group.bar.height = 0,
               draw.lines = F, hjust = 0)+ 
  scale_fill_gradientn(colors = c("blue", "white","red"), na.value = "white") +
  NoLegend()+
  theme_bw()+
  scale_y_discrete(labels=markers1, position = "right" )+
  scale_x_discrete( limits = rev(levels(hamster.integrated$seurat_clusters)))+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, vjust = -2, hjust = -0.15))


#### ASSIGN CELLTYPE ANNOTATION TO CLUSTERS ####

hamster.integrated <- RenameIdents(hamster.integrated, `0` = "Neutrophil", `1` = "B", `2` = "T", 
                                   `3` = "Classical monocyte", `4` = "Classical monocyte", `5` = "Neutrophil", `6` = "Immature neutrophil", 
                                   `7` = "Immature neutrophil", `8` = "B", '9'="NK", '10'="Neutrophil", '11'="Neutrophil", 
                                   '12'="T", '13'="Activated T", '14'="Non-classical monocyte", '15' = "B", '16' = "B", "17" = "Classical monocyte",
                                   "18"= "T", "19" = "mDC", "20"="Platelet", "21"="B", '22' = "pDC", '23' ="Mixed", '24'="Classical monocyte")
Idents(hamster.integrated) <- factor(Idents(hamster.integrated), levels = c( "Neutrophil", "Immature neutrophil",   "Classical monocyte",   "Non-classical monocyte", "NK", "T", "Activated T", "B", "Plasma cell", "mDC", "pDC", "Platelet" , "Mixed"))
hamster.integrated$celltype <- Idents(hamster.integrated)


## get counts for celltypes per hamster
tbl <- (table(hamster.integrated$infection, hamster.integrated$celltype))

# convert to percentage

tbl.perc <- round(tbl/rowSums(tbl) *100, 2)
tbl.perc <- data.frame(tbl.perc)
colnames(tbl.perc) <- c("infection", "celltype", "percentage")


## determine mean percentage of celltype composition + std dev 

tbl2 <- (table(hamster.integrated$hamster, hamster.integrated$celltype))
tbl2.perc <- round(tbl2/rowSums(tbl2) *100, 2)
tbl2.perc <- data.frame(tbl2.perc)
colnames(tbl2.perc) <- c("hamster", "celltype", "percentage")
tbl2.perc$infection <- gsub("(d[0-5]{1,2})_[0-5]", "\\1", tbl2.perc$hamster)

##function to calc mean and std dev
data_summary <- function(data, varname, groupname){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=T), sd = sd(x[[col]], na.rm=T))
  }
  data_sum <- ddply(data, groupname, .fun=summary_func, varname )
  data_sum <- rename(data_sum, c("mean"=varname))
  return(data_sum)
}

df2 <- data_summary(tbl2.perc, varname="percentage", groupname=c("celltype","infection"))
df2$infection <- factor(df2$infection, levels = c("d0", "d1", "d2", "d3", "d5", "d14"))
