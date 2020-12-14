library(ggplot2)
library(stringr)
library(Seurat)
library(HDF5Array)
library(DelayedArray)
library(DelayedMatrixStats)
require(gridExtra)

data_dir <- './ma-d0-lung-1/ma-d0-lung-1/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d0_lung_1 = CreateSeuratObject(counts = data, project="ma_d0_lung_1", min.cells=5, min.features=1000)

data_dir <- './ma-d0-lung-2/ma-d0-lung-2/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d0_lung_2 = CreateSeuratObject(counts = data, project="ma_d0_lung_2", min.cells=5, min.features=1000)

data_dir <- './ma-d0-lung-3/ma-d0-lung-3/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d0_lung_3 = CreateSeuratObject(counts = data, project="ma_d0_lung_3", min.cells=5, min.features=1000)

data_dir <- './ma-d2-lung-1/ma-d2-lung-1/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d2_lung_1 = CreateSeuratObject(counts = data, project="ma_d2_lung_1", min.cells=5, min.features=1000)

data_dir <- './ma-d2-lung-2/ma-d2-lung-2/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d2_lung_2 = CreateSeuratObject(counts = data, project="ma_d2_lung_2", min.cells=5, min.features=1000)

data_dir <- './ma-d2-lung-3/ma-d2-lung-3/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d2_lung_3 = CreateSeuratObject(counts = data, project="ma_d2_lung_3", min.cells=5, min.features=1000)

data_dir <- './ma-d3-lung-1/ma-d3-lung-1/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d3_lung_1 = CreateSeuratObject(counts = data, project="ma_d3_lung_1", min.cells=5, min.features=1000)

data_dir <- './ma-d3-lung-2/ma-d3-lung-2/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d3_lung_2 = CreateSeuratObject(counts = data, project="ma_d3_lung_2", min.cells=5, min.features=1000)

data_dir <- './ma-d3-lung-3/ma-d3-lung-3/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d3_lung_3 = CreateSeuratObject(counts = data, project="ma_d3_lung_3", min.cells=5, min.features=1000)

data_dir <- './ma-d5-lung-1/ma-d5-lung-1/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d5_lung_1 = CreateSeuratObject(counts = data, project="ma_d5_lung_1", min.cells=5, min.features=1000)

data_dir <- './ma-d5-lung-2/ma-d5-lung-2/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d5_lung_2 = CreateSeuratObject(counts = data, project="ma_d5_lung_2", min.cells=5, min.features=1000)

data_dir <- './ma-d5-lung-3/ma-d5-lung-3/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d5_lung_3 = CreateSeuratObject(counts = data, project="ma_d5_lung_3", min.cells=5, min.features=1000)

data_dir <- './ma-e14-lung-1/ma-e14-lung-1/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_e14_lung_1 = CreateSeuratObject(counts = data, project="ma_e14_lung_1", min.cells=5, min.features=1000)

data_dir <- './ma-e14-lung-2/ma-e14-lung-2/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_e14_lung_2 = CreateSeuratObject(counts = data, project="ma_e14_lung_2", min.cells=5, min.features=1000)

data_dir <- './ma-e14-lung-3/ma-e14-lung-3/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_e14_lung_3 = CreateSeuratObject(counts = data, project="ma_e14_lung_3", min.cells=5, min.features=1000)

data_dir <- './ma-d0-blood-1/ma-d0-blood-1/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d0_blood_1 = CreateSeuratObject(counts = data, project="ma_d0_blood_1", min.cells=5, min.features=500)

data_dir <- './ma-d0-blood-2/ma-d0-blood-2/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d0_blood_2 = CreateSeuratObject(counts = data, project="ma_d0_blood_2", min.cells=5, min.features=500)

data_dir <- './ma-d0-blood-3/ma-d0-blood-3/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d0_blood_3 = CreateSeuratObject(counts = data, project="ma_d0_blood_3", min.cells=5, min.features=500)

data_dir <- './ma-d2-blood-1/ma-d2-blood-1/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d2_blood_1 = CreateSeuratObject(counts = data, project="ma_d2_blood_1", min.cells=5, min.features=500)

data_dir <- './ma-d2-blood-2/ma-d2-blood-2/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d2_blood_2 = CreateSeuratObject(counts = data, project="ma_d2_blood_2", min.cells=5, min.features=500)

data_dir <- './ma-d2-blood-3/ma-d2-blood-3/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d2_blood_3 = CreateSeuratObject(counts = data, project="ma_d2_blood_3", min.cells=5, min.features=500)

data_dir <- './ma-d3-blood-1/ma-d3-blood-1/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d3_blood_1 = CreateSeuratObject(counts = data, project="ma_d3_blood_1", min.cells=5, min.features=500)

data_dir <- './ma-d3-blood-2/ma-d3-blood-2/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d3_blood_2 = CreateSeuratObject(counts = data, project="ma_d3_blood_2", min.cells=5, min.features=500)

data_dir <- './ma-d3-blood-3/ma-d3-blood-3/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d3_blood_3 = CreateSeuratObject(counts = data, project="ma_d3_blood_3", min.cells=5, min.features=500)

data_dir <- './ma-d5-blood-1/ma-d5-blood-1/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d5_blood_1 = CreateSeuratObject(counts = data, project="ma_d5_blood_1", min.cells=5, min.features=500)

data_dir <- './ma-d5-blood-2/ma-d5-blood-2/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d5_blood_2 = CreateSeuratObject(counts = data, project="ma_d5_blood_2", min.cells=5, min.features=500)

data_dir <- './ma-d5-blood-3/ma-d5-blood-3/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d5_blood_3 = CreateSeuratObject(counts = data, project="ma_d5_blood_3", min.cells=5, min.features=500)

data_dir <- './ma-e14-blood-1/ma-e14-blood-1/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_e14_blood_1 = CreateSeuratObject(counts = data, project="ma_e14_blood_1", min.cells=5, min.features=500)

data_dir <- './ma-e14-blood-2/ma-e14-blood-2/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_e14_blood_2 = CreateSeuratObject(counts = data, project="ma_e14_blood_2", min.cells=5, min.features=500)

data_dir <- './ma-e14-blood-3/ma-e14-blood-3/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_e14_blood_3 = CreateSeuratObject(counts = data, project="ma_e14_blood_3", min.cells=5, min.features=500)

seu_red <- merge(ma_d0_blood_1, y = c(ma_d0_blood_2, ma_d0_blood_3, ma_d2_blood_1, ma_d2_blood_2, ma_d2_blood_3, ma_d3_blood_1, ma_d3_blood_2, ma_d3_blood_3, ma_d5_blood_1, ma_d5_blood_2, ma_d5_blood_3, ma_e14_blood_1, ma_e14_blood_2, ma_e14_blood_3), add.cell.ids = c("ma_d0_blood_1", "ma_d0_blood_2", "ma_d0_blood_3", "ma_d2_blood_1", "ma_d2_blood_2", "ma_d2_blood_3", "ma_d3_blood_1", "ma_d3_blood_2", "ma_d3_blood_3", "ma_d5_blood_1", "ma_d5_blood_2", "ma_d5_blood_3", "ma_e14_blood_1", "ma_e14_blood_2", "ma_e14_blood_3"), project = "ma_blood")
seu_red[["percent.mt"]] <- PercentageFeatureSet(seu_red, features =  c("COX1", "CYTB", "ND1", "ND2",  "ND4", "ND5", "ND6"))
seu_red <- subset(seu_red, subset = percent.mt < 7)
saveRDS(seu_red, "./seu_blood_combined_mtfilt_blood500.rds")

seu_red <- merge(ma_d0_lung_1, y = c(ma_d0_lung_2, ma_d0_lung_3, ma_d2_lung_1, ma_d2_lung_2, ma_d2_lung_3, ma_d3_lung_1, ma_d3_lung_2, ma_d3_lung_3, ma_d5_lung_1, ma_d5_lung_2, ma_d5_lung_3, ma_e14_lung_1, ma_e14_lung_2, ma_e14_lung_3), add.cell.ids = c("ma_d0_lung_1", "ma_d0_lung_2", "ma_d0_lung_3", "ma_d2_lung_1", "ma_d2_lung_2", "ma_d2_lung_3", "ma_d3_lung_1", "ma_d3_lung_2", "ma_d3_lung_3", "ma_d5_lung_1", "ma_d5_lung_2", "ma_d5_lung_3", "ma_e14_lung_1", "ma_e14_lung_2", "ma_e14_lung_3"), project = "ma_lung")
seu_red[["percent.mt"]] <- PercentageFeatureSet(seu_red, features =  c("COX1", "CYTB", "ND1", "ND2",  "ND4", "ND5", "ND6"))
seu_red <- subset(seu_red, subset = percent.mt < 7)
saveRDS(seu_red, "./seu_lung_combined_mtfilt.rds")

seu_red <- merge(ma_d2_blood_1, y = c(ma_d2_blood_2, ma_d2_blood_3, ma_d2_lung_1, ma_d2_lung_2, ma_d2_lung_3), add.cell.ids = c("ma_d2_blood_1", "ma_d2_blood_2", "ma_d2_blood_3", "ma_d2_lung_1", "ma_d2_lung_2", "ma_d2_lung_3"), project = "ma_combined_d2")
seu_red[["percent.mt"]] <- PercentageFeatureSet(seu_red, features =  c("COX1", "CYTB", "ND1", "ND2",  "ND4", "ND5", "ND6"))
seu_red <- subset(seu_red, subset = percent.mt < 7)
saveRDS(seu_red, "./seu_d2_combined_mtfilt_blood500.rds")

seu_red <- merge(ma_d3_blood_1, y = c(ma_d3_blood_2, ma_d3_blood_3, ma_d3_lung_1, ma_d3_lung_2, ma_d3_lung_3), add.cell.ids = c("ma_d3_blood_1", "ma_d3_blood_2", "ma_d3_blood_3", "ma_d3_lung_1", "ma_d3_lung_2", "ma_d3_lung_3"), project = "ma_combined_d3")
seu_red[["percent.mt"]] <- PercentageFeatureSet(seu_red, features =  c("COX1", "CYTB", "ND1", "ND2",  "ND4", "ND5", "ND6"))
seu_red <- subset(seu_red, subset = percent.mt < 7)
saveRDS(seu_red, "./seu_d3_combined_mtfilt_blood500.rds")

seu_red <- merge(ma_d5_blood_1, y = c(ma_d5_blood_2, ma_d5_blood_3, ma_d5_lung_1, ma_d5_lung_2, ma_d5_lung_3), add.cell.ids = c("ma_d5_blood_1", "ma_d5_blood_2", "ma_d5_blood_3", "ma_d5_lung_1", "ma_d5_lung_2", "ma_d5_lung_3"), project = "ma_combined_d5")
seu_red[["percent.mt"]] <- PercentageFeatureSet(seu_red, features =  c("COX1", "CYTB", "ND1", "ND2",  "ND4", "ND5", "ND6"))
seu_red <- subset(seu_red, subset = percent.mt < 7)
saveRDS(seu_red, "./seu_d5_combined_mtfilt_blood500.rds")

seu_red <- merge(ma_d0_blood_1, y = c(ma_d0_blood_2, ma_d0_blood_3, ma_d0_lung_1, ma_d0_lung_2, ma_d0_lung_3), add.cell.ids = c("ma_d0_blood_1", "ma_d0_blood_2", "ma_d0_blood_3", "ma_d0_lung_1", "ma_d0_lung_2", "ma_d0_lung_3"), project = "ma_combined_d0")
seu_red[["percent.mt"]] <- PercentageFeatureSet(seu_red, features =  c("COX1", "CYTB", "ND1", "ND2",  "ND4", "ND5", "ND6"))
seu_red <- subset(seu_red, subset = percent.mt < 7)
saveRDS(seu_red, "./seu_d0_combined_mtfilt_blood500.rds")

seu_red <- merge(ma_e14_blood_1, y = c(ma_e14_blood_2, ma_e14_blood_3, ma_e14_lung_1, ma_e14_lung_2, ma_e14_lung_3), add.cell.ids = c("ma_e14_blood_1", "ma_e14_blood_2", "ma_e14_blood_3", "ma_e14_lung_1", "ma_e14_lung_2", "ma_e14_lung_3"), project = "ma_blood_d14")
seu_red[["percent.mt"]] <- PercentageFeatureSet(seu_red, features =  c("COX1", "CYTB", "ND1", "ND2",  "ND4", "ND5", "ND6"))
seu_red <- subset(seu_red, subset = percent.mt < 7)
saveRDS(seu_red, "./seu_d14_combined_mtfilt_blood500.rds")

data_dir <- './ma-d0-blood-1/ma-d0-blood-1/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d0_blood_1 = CreateSeuratObject(counts = data, project="ma_d0_blood_1", min.cells=5, min.features=1000)

data_dir <- './ma-d0-blood-2/ma-d0-blood-2/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d0_blood_2 = CreateSeuratObject(counts = data, project="ma_d0_blood_2", min.cells=5, min.features=1000)

data_dir <- './ma-d0-blood-3/ma-d0-blood-3/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d0_blood_3 = CreateSeuratObject(counts = data, project="ma_d0_blood_3", min.cells=5, min.features=1000)

data_dir <- './ma-d2-blood-1/ma-d2-blood-1/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d2_blood_1 = CreateSeuratObject(counts = data, project="ma_d2_blood_1", min.cells=5, min.features=1000)

data_dir <- './ma-d2-blood-2/ma-d2-blood-2/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d2_blood_2 = CreateSeuratObject(counts = data, project="ma_d2_blood_2", min.cells=5, min.features=1000)

data_dir <- './ma-d2-blood-3/ma-d2-blood-3/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d2_blood_3 = CreateSeuratObject(counts = data, project="ma_d2_blood_3", min.cells=5, min.features=1000)

data_dir <- './ma-d3-blood-1/ma-d3-blood-1/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d3_blood_1 = CreateSeuratObject(counts = data, project="ma_d3_blood_1", min.cells=5, min.features=1000)

data_dir <- './ma-d3-blood-2/ma-d3-blood-2/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d3_blood_2 = CreateSeuratObject(counts = data, project="ma_d3_blood_2", min.cells=5, min.features=1000)

data_dir <- './ma-d3-blood-3/ma-d3-blood-3/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d3_blood_3 = CreateSeuratObject(counts = data, project="ma_d3_blood_3", min.cells=5, min.features=1000)

data_dir <- './ma-d5-blood-1/ma-d5-blood-1/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d5_blood_1 = CreateSeuratObject(counts = data, project="ma_d5_blood_1", min.cells=5, min.features=1000)

data_dir <- './ma-d5-blood-2/ma-d5-blood-2/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d5_blood_2 = CreateSeuratObject(counts = data, project="ma_d5_blood_2", min.cells=5, min.features=1000)

data_dir <- './ma-d5-blood-3/ma-d5-blood-3/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_d5_blood_3 = CreateSeuratObject(counts = data, project="ma_d5_blood_3", min.cells=5, min.features=1000)

data_dir <- './ma-e14-blood-1/ma-e14-blood-1/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_e14_blood_1 = CreateSeuratObject(counts = data, project="ma_e14_blood_1", min.cells=5, min.features=1000)

data_dir <- './ma-e14-blood-2/ma-e14-blood-2/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_e14_blood_2 = CreateSeuratObject(counts = data, project="ma_e14_blood_2", min.cells=5, min.features=1000)

data_dir <- './ma-e14-blood-3/ma-e14-blood-3/outs/raw_feature_bc_matrix/'
data <- Read10X(data.dir = data_dir)
ma_e14_blood_3 = CreateSeuratObject(counts = data, project="ma_e14_blood_3", min.cells=5, min.features=1000)

seu_red <- merge(ma_d2_blood_1, y = c(ma_d2_blood_2, ma_d2_blood_3, ma_d2_lung_1, ma_d2_lung_2, ma_d2_lung_3), add.cell.ids = c("ma_d2_blood_1", "ma_d2_blood_2", "ma_d2_blood_3", "ma_d2_lung_1", "ma_d2_lung_2", "ma_d2_lung_3"), project = "ma_combined_d2")
seu_red[["percent.mt"]] <- PercentageFeatureSet(seu_red, features =  c("COX1", "CYTB", "ND1", "ND2",  "ND4", "ND5", "ND6"))
seu_red <- subset(seu_red, subset = percent.mt < 7)
saveRDS(seu_red, "./seu_d2_combined_mtfilt.rds")

seu_red <- merge(ma_d3_blood_1, y = c(ma_d3_blood_2, ma_d3_blood_3, ma_d3_lung_1, ma_d3_lung_2, ma_d3_lung_3), add.cell.ids = c("ma_d3_blood_1", "ma_d3_blood_2", "ma_d3_blood_3", "ma_d3_lung_1", "ma_d3_lung_2", "ma_d3_lung_3"), project = "ma_combined_d3")
seu_red[["percent.mt"]] <- PercentageFeatureSet(seu_red, features =  c("COX1", "CYTB", "ND1", "ND2",  "ND4", "ND5", "ND6"))
seu_red <- subset(seu_red, subset = percent.mt < 7)
saveRDS(seu_red, "./seu_d3_combined_mtfilt.rds")

seu_red <- merge(ma_d5_blood_1, y = c(ma_d5_blood_2, ma_d5_blood_3, ma_d5_lung_1, ma_d5_lung_2, ma_d5_lung_3), add.cell.ids = c("ma_d5_blood_1", "ma_d5_blood_2", "ma_d5_blood_3", "ma_d5_lung_1", "ma_d5_lung_2", "ma_d5_lung_3"), project = "ma_combined_d5")
seu_red[["percent.mt"]] <- PercentageFeatureSet(seu_red, features =  c("COX1", "CYTB", "ND1", "ND2",  "ND4", "ND5", "ND6"))
seu_red <- subset(seu_red, subset = percent.mt < 7)
saveRDS(seu_red, "./seu_d5_combined_mtfilt.rds")

seu_red <- merge(ma_d0_blood_1, y = c(ma_d0_blood_2, ma_d0_blood_3, ma_d0_lung_1, ma_d0_lung_2, ma_d0_lung_3), add.cell.ids = c("ma_d0_blood_1", "ma_d0_blood_2", "ma_d0_blood_3", "ma_d0_lung_1", "ma_d0_lung_2", "ma_d0_lung_3"), project = "ma_combined_d0")
seu_red[["percent.mt"]] <- PercentageFeatureSet(seu_red, features =  c("COX1", "CYTB", "ND1", "ND2",  "ND4", "ND5", "ND6"))
seu_red <- subset(seu_red, subset = percent.mt < 7)
saveRDS(seu_red, "./seu_d0_combined_mtfilt.rds")

seu_red <- merge(ma_e14_blood_1, y = c(ma_e14_blood_2, ma_e14_blood_3, ma_e14_lung_1, ma_e14_lung_2, ma_e14_lung_3), add.cell.ids = c("ma_e14_blood_1", "ma_e14_blood_2", "ma_e14_blood_3", "ma_e14_lung_1", "ma_e14_lung_2", "ma_e14_lung_3"), project = "ma_combined_d14")
seu_red[["percent.mt"]] <- PercentageFeatureSet(seu_red, features =  c("COX1", "CYTB", "ND1", "ND2",  "ND4", "ND5", "ND6"))
seu_red <- subset(seu_red, subset = percent.mt < 7)
saveRDS(seu_red, "./seu_d14_combined_mtfilt.rds")


