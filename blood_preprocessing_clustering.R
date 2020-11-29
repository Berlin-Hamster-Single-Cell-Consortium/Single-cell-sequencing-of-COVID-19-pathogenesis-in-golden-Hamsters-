library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)
library(ggplot2)
library(ggpubr)
#------------------------------------------------------------------------------------------

#### PREPROCESS DATA ####


##Read in data from CellRanger output

hamster_blood <- Read10X(data.dir = "../path/to/CellRanger/output/")
# Initialize the Seurat object with the raw (non-normalized data).
hamster_all <- CreateSeuratObject(counts = hamster_blood, project = "pbmc3k", min.cells = 3, min.features = 200)
# Filter object to cells with at least 500 genes per cell
hamster_all <- subset(pbmc, subset = nFeature_RNA >500)
# determine percentage mitochondrial reads and subset
hamster_all[["percent.mt"]] <- PercentageFeatureSet(hamster_all, features =  c("COX1", "CYTB", "ND1", "ND2",  "ND4", "ND5", "ND6"))
hamster_all <- subset(hamster_all, percent.mt < 18)


DefaultAssay(hamster_all) <- "RNA"


## integrate by hamster to remove batch effects

hamster.list <- SplitObject(hamster_all, split.by = "orig.ident")
for (i in 1:length(hamster.list)){
  hamster.list[[i]] <- SCTransform(hamster.list[[i]], verbose =T)
}

hamster.features <- SelectIntegrationFeatures(object.list = hamster.list, nfeatures = 3000)
hamster.list <- PrepSCTIntegration(object.list = hamster.list, anchor.features = hamster.features, 
                                   verbose = T)

hamster.anchors <- FindIntegrationAnchors(object.list = hamster.list, normalization.method = "SCT", 
                                          anchor.features = hamster.features, verbose = T)
hamster.integrated <- IntegrateData(anchorset = hamster.anchors, normalization.method = "SCT", 
                                    verbose = T)

## run dimensional reductions
#   PCA
hamster.integrated<- RunPCA(hamster.integrated, verbose = FALSE)
#   UMAP
hamster.integrated<- RunUMAP(hamster.integrated, dims = 1:30, verbose = FALSE)


#### CLUSTERING AND ANNOTATION ####

##cluster
hamster.integrated <- FindNeighbors(hamster.integrated, dims = 1:30)
hamster.integrated <- FindClusters(hamster.integrated, resolution = 0.5)

##annotate by treatment
hamster.integrated$hamster <- gsub("ma_([de][0-5]{1,2})_blood(_[0-5]{1,2})","\\1\\2",hamster.integrated@meta.data$orig.ident, )
hamster.integrated$hamster <- gsub("e", "d", hamster.integrated$hamster)
hamster.integrated$infection <- gsub("_[0-5]{1,2}","",hamster.integrated@meta.data$hamster )
unique(hamster.integrated$hamster)
unique(hamster.integrated$infection)


##annotate libraries and barcodes
unique(hamster.integrated@meta.data$orig.ident)
hamster.integrated$library <- ifelse(hamster.integrated$orig.ident == "ma_d0_blood_1", "-1",
                                     ifelse(hamster.integrated$orig.ident == "ma_d0_blood_2", "-2",
                                            ifelse(hamster.integrated$orig.ident == "ma_d0_blood_3", "-3",
                                                   ifelse(hamster.integrated$orig.ident == "ma_d2_blood_1", "-4",
                                                          ifelse(hamster.integrated$orig.ident == "ma_d2_blood_2", "-5",
                                                                 ifelse(hamster.integrated$orig.ident == "ma_d2_blood_3", "-6",
                                                                        ifelse(hamster.integrated$orig.ident == "ma_d3_blood_1", "-7",
                                                                               ifelse(hamster.integrated$orig.ident == "ma_d3_blood_2", "-8",
                                                                                      ifelse(hamster.integrated$orig.ident == "ma_d3_blood_3", "-9", 
                                                                                             ifelse(hamster.integrated$orig.ident == "ma_d5_blood_1", "-10",
                                                                                                    ifelse(hamster.integrated$orig.ident == "ma_d5_blood_2", "-11",
                                                                                                           ifelse(hamster.integrated$orig.ident == "ma_d5_blood_3", "-12", 
                                                                                                                  ifelse(hamster.integrated$orig.ident == "ma_e14_blood_1", "-13",
                                                                                                                         ifelse(hamster.integrated$orig.ident =="ma_e14_blood_2", "-14",
                                                                                                                                ifelse(hamster.integrated$orig.ident =="ma_e14_blood_3", "-15", "none")))))))))))))))





hamster.integrated$barcode <- gsub("^ma_.*_[0-9]{1,2}_", "", row.names(hamster.integrated@meta.data))
hamster.integrated$barcode <- paste(hamster.integrated$barcode, hamster.integrated$library, sep="")


##find top 10 marker genes of each cluster
cluster_marks <- FindAllMarkers(hamster.integrated)
cluster_marks.top <- cluster_marks %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC) %>% write.csv(file = "top10_clustmarks.csv")



#### FILTER ERYTHROCYTES AND DEAD CELLS ####

#Cells filtered out if clusters have erythrocyte marker genes:
  # Fam46C, Snca, Alas2
#as well as significantly high expression of mitochondrial genes (see above)


hamster.integrated <- subset(hamster.integrated, subset = seurat_clusters != 8 & seurat_clusters != 17 & seurat_clusters != 21 & seurat_clusters != 23 & seurat_clusters != 20)

##repeat above integration and clustering workflow
hamster.list <- SplitObject(hamster.integrated, split.by = "orig.ident")
for (i in 1:length(hamster.list)){
  hamster.list[[i]] <- SCTransform(hamster.list[[i]], verbose =T)
}
hamster.features <- SelectIntegrationFeatures(object.list = hamster.list, nfeatures = 3000)
hamster.list <- PrepSCTIntegration(object.list = hamster.list, anchor.features = hamster.features, 
                                   verbose = T)
hamster.anchors <- FindIntegrationAnchors(object.list = hamster.list, normalization.method = "SCT", 
                                          anchor.features = hamster.features, verbose = T)
hamster.integrated <- IntegrateData(anchorset = hamster.anchors, normalization.method = "SCT", 
                                    verbose = T)
hamster.integrated<- RunPCA(hamster.integrated, verbose = FALSE)
hamster.integrated<- RunUMAP(hamster.integrated, dims = 1:30, verbose = FALSE)
hamster.integrated <- FindNeighbors(hamster.integrated, dims = 1:30)
hamster.integrated <- FindClusters(hamster.integrated, resolution = 0.5)

##check cluster markers again
cluster_marks.top <- cluster_marks %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC) %>% write.csv(file = "top10_clustmarks_first_filter.csv")

##remove two clusters still expressing high levels of erythrocyte marker genes
hamster.integrated <- subset(hamster.integrated, subset = seurat_clusters != 19 & seurat_clusters != 24)

##repeat above integration and clustering workflow
hamster.list <- SplitObject(hamster.integrated, split.by = "orig.ident")
for (i in 1:length(hamster.list)){
  hamster.list[[i]] <- SCTransform(hamster.list[[i]], verbose =T)
}
hamster.features <- SelectIntegrationFeatures(object.list = hamster.list, nfeatures = 3000)
hamster.list <- PrepSCTIntegration(object.list = hamster.list, anchor.features = hamster.features, 
                                   verbose = T)
hamster.anchors <- FindIntegrationAnchors(object.list = hamster.list, normalization.method = "SCT", 
                                          anchor.features = hamster.features, verbose = T)
hamster.integrated <- IntegrateData(anchorset = hamster.anchors, normalization.method = "SCT", 
                                    verbose = T)
hamster.integrated<- RunPCA(hamster.integrated, verbose = FALSE)
hamster.integrated<- RunUMAP(hamster.integrated, dims = 1:30, verbose = FALSE)
hamster.integrated <- FindNeighbors(hamster.integrated, dims = 1:30)
hamster.integrated <- FindClusters(hamster.integrated, resolution = 0.5)






