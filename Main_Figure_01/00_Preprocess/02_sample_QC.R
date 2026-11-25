setwd('Data11/')
#GZ8------------------------
GZ8 <- Read10X(data.dir = "../../runs/GZ0008/outs/multi/count/raw_feature_bc_matrix/")
GZ8 <- CreateSeuratObject(counts = GZ8, project = "GZ8", min.cells = 3, min.features = 200)
GZ8 <- RenameCells(GZ8, add.cell.id = 'GZ8')
GZ8[["percent.mt"]] <- PercentageFeatureSet(GZ8, pattern = "^MT-")
VlnPlot(GZ8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(GZ8,'GZ8QC前.rds')#12006
GZ8 <- subset(GZ8, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(GZ8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#11667
saveRDS(GZ8,'GZ8QC后.rds')
head(GZ8@meta.data)
#YZ8------------------------
YZ8 <- Read10X(data.dir = "../../runs/YZ0008/outs/multi/count/raw_feature_bc_matrix/")
YZ8 <- CreateSeuratObject(counts = YZ8, project = "YZ8", min.cells = 3, min.features = 200)
YZ8 <- RenameCells(YZ8, add.cell.id = 'YZ8')
YZ8[["percent.mt"]] <- PercentageFeatureSet(YZ8, pattern = "^MT-")
VlnPlot(YZ8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(YZ8,'YZ8QC前.rds')#6922
YZ8 <- subset(YZ8, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(YZ8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#6745
saveRDS(YZ8,'YZ8QC后.rds')
#YZ8A------------------------
YZ8A <- Read10X(data.dir = "../runs/CellRanger/YZ0008a/outs/multi/count/raw_feature_bc_matrix/")
YZ8A <- CreateSeuratObject(counts = YZ8A, project = "YZ8A", min.cells = 3, min.features = 200)
YZ8A <- RenameCells(YZ8A, add.cell.id = 'YZ8A')
YZ8A[["percent.mt"]] <- PercentageFeatureSet(YZ8A, pattern = "^MT-")
VlnPlot(YZ8A, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(YZ8A,'YZ8AQC前.rds')#11197
YZ8A <- subset(YZ8A, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(YZ8A, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#10833
saveRDS(YZ8A,'YZ8AQC后.rds')

#YZ8C------------------------
YZ8C <- Read10X(data.dir = "../../runs/CellRanger/YZ0008c/outs/multi/count/raw_feature_bc_matrix/")
YZ8C <- CreateSeuratObject(counts = YZ8C, project = "YZ8C", min.cells = 3, min.features = 200)
YZ8C <- RenameCells(YZ8C, add.cell.id = 'YZ8C')
YZ8C[["percent.mt"]] <- PercentageFeatureSet(YZ8C, pattern = "^MT-")
p<- VlnPlot(YZ8C, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("YZ8C.png", plot = p, width = 6, height = 4, units = "in", dpi = 300)
saveRDS(YZ8C,'YZ8CQC前.rds')#6833
YZ8C <- subset(YZ8C, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(YZ8C, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#6745
saveRDS(YZ8C,'YZ8CQC后.rds')
head(YZ8C@meta.data)

#YZ8D------------------------
YZ8D <- Read10X(data.dir = "../../runs/CellRanger/Y29D/outs/multi/count/raw_feature_bc_matrix/")
YZ8D <- CreateSeuratObject(counts = YZ8D, project = "YZ8D", min.cells = 3, min.features = 200)
YZ8D <- RenameCells(YZ8D, add.cell.id = 'YZ8D')
YZ8D[["percent.mt"]] <- PercentageFeatureSet(YZ8D, pattern = "^MT-")
VlnPlot(YZ8D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(YZ8D,'YZ8DQC前.rds')#4867
YZ8D <- subset(YZ8D, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(YZ8D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#4840
saveRDS(YZ8D,'YZ8DQC后.rds')
head(YZ8D@meta.data)

#GZ5------------------------
GZ5 <- Read10X(data.dir = "../../runs/CellRanger/GZ0005/outs/multi/count/raw_feature_bc_matrix/")
GZ5 <- CreateSeuratObject(counts = GZ5, project = "GZ5", min.cells = 3, min.features = 200)
GZ5 <- RenameCells(GZ5, add.cell.id = 'GZ5')
GZ5[["percent.mt"]] <- PercentageFeatureSet(GZ5, pattern = "^MT-")
VlnPlot(GZ5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(GZ5,'GZ5QC前.rds')#3267
GZ5 <- subset(GZ5, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(GZ5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#3140
saveRDS(GZ5,'GZ5QC后.rds')
head(GZ5@meta.data)
#YZ5------------------------
YZ5 <- Read10X(data.dir = "../../runs/YZ0005/outs/multi/count/raw_feature_bc_matrix/")
YZ5 <- CreateSeuratObject(counts = YZ5, project = "YZ5", min.cells = 3, min.features = 200)
YZ5 <- RenameCells(YZ5, add.cell.id = 'YZ5')
YZ5[["percent.mt"]] <- PercentageFeatureSet(YZ5, pattern = "^MT-")
VlnPlot(YZ5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(YZ5,'YZ5QC前.rds')#8215
YZ5 <- subset(YZ5, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(YZ5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#7528
saveRDS(YZ5,'YZ5QC后.rds')
head(YZ5@meta.data)
#YZ5A------------------------
YZ5A <- Read10X(data.dir = "../../runs/CellRanger/YZ0005a/outs/multi/count/raw_feature_bc_matrix/")
YZ5A <- CreateSeuratObject(counts = YZ5A, project = "YZ5A", min.cells = 3, min.features = 200)
YZ5A <- RenameCells(YZ5A, add.cell.id = 'YZ5A')
YZ5A[["percent.mt"]] <- PercentageFeatureSet(YZ5A, pattern = "^MT-")
VlnPlot(YZ5A, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(YZ5A,'YZ5AQC前.rds')#9022
YZ5A <- subset(YZ5A, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(YZ5A, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#8827
saveRDS(YZ5A,'YZ5AQC后.rds')
head(YZ5A@meta.data)
#YZ5C------------------------
YZ5C <- Read10X(data.dir = "../../runs/CellRanger/YZ5C/outs/multi/count/raw_feature_bc_matrix/")
YZ5C <- CreateSeuratObject(counts = YZ5C, project = "YZ5C", min.cells = 3, min.features = 200)
YZ5C <- RenameCells(YZ5C, add.cell.id = 'YZ5C')
YZ5C[["percent.mt"]] <- PercentageFeatureSet(YZ5C, pattern = "^MT-")
VlnPlot(YZ5C, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(YZ5C,'YZ5CQC前.rds')#12325
YZ5C <- subset(YZ5C, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(YZ5C, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#11899
saveRDS(YZ5C,'YZ5CQC后.rds')
head(YZ5C@meta.data)
#YZ5D------------------------
YZ5D <- Read10X(data.dir = "../../runs/CellRanger/Y25D/outs/multi/count/raw_feature_bc_matrix/")
YZ5D <- CreateSeuratObject(counts = YZ5D, project = "YZ5D", min.cells = 3, min.features = 200)
YZ5D <- RenameCells(YZ5D, add.cell.id = 'YZ5D')
YZ5D[["percent.mt"]] <- PercentageFeatureSet(YZ5D, pattern = "^MT-")
VlnPlot(YZ5D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(YZ5D,'YZ5DQC前.rds')#7687
YZ5D <- subset(YZ5D, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(YZ5D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#7451
saveRDS(YZ5D,'YZ5DQC后.rds')
head(YZ5D@meta.data)

#GZ6------------------------
GZ6 <- Read10X(data.dir = "../../runs/GZ0006/outs/multi/count/raw_feature_bc_matrix/")
GZ6 <- CreateSeuratObject(counts = GZ6, project = "GZ6", min.cells = 3, min.features = 200)
GZ6 <- RenameCells(GZ6, add.cell.id = 'GZ6')
GZ6[["percent.mt"]] <- PercentageFeatureSet(GZ6, pattern = "^MT-")
VlnPlot(GZ6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(GZ6,'GZ6QC前.rds')#10565
GZ6 <- subset(GZ6, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(GZ6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#10176
saveRDS(GZ6,'GZ6QC后.rds')
head(GZ6@meta.data)
#YZ6------------------------
YZ6 <- Read10X(data.dir = "../../runs/YZ0006/outs/multi/count/raw_feature_bc_matrix/")
YZ6 <- CreateSeuratObject(counts = YZ6, project = "YZ6", min.cells = 3, min.features = 200)
YZ6 <- RenameCells(YZ6, add.cell.id = 'YZ6')
YZ6[["percent.mt"]] <- PercentageFeatureSet(YZ6, pattern = "^MT-")
VlnPlot(YZ6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(YZ6,'YZ6QC前.rds')#10121
YZ6 <- subset(YZ6, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(YZ6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#9978
saveRDS(YZ6,'YZ6QC后.rds')
head(YZ6@meta.data)
#YZ6A------------------------
YZ6A <- Read10X(data.dir = "../../runs/YZ0006a/outs/multi/count/raw_feature_bc_matrix/")
YZ6A <- CreateSeuratObject(counts = YZ6A, project = "YZ6A", min.cells = 3, min.features = 200)
YZ6A <- RenameCells(YZ6A, add.cell.id = 'YZ6A')
YZ6A[["percent.mt"]] <- PercentageFeatureSet(YZ6A, pattern = "^MT-")
VlnPlot(YZ6A, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(YZ6A,'YZ6AQC前.rds')#10065
YZ6A <- subset(YZ6A, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(YZ6A, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#9426
saveRDS(YZ6A,'YZ6AQC后.rds')
head(YZ6A@meta.data)
#YZ6C------------------------
YZ6C <- Read10X(data.dir = "../../runs/YZ0006c/outs/multi/count/raw_feature_bc_matrix/")
YZ6C <- CreateSeuratObject(counts = YZ6C, project = "YZ6C", min.cells = 3, min.features = 200)
YZ6C <- RenameCells(YZ6C, add.cell.id = 'YZ6C')
YZ6C[["percent.mt"]] <- PercentageFeatureSet(YZ6C, pattern = "^MT-")
VlnPlot(YZ6C, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(YZ6C,'YZ6CQC前.rds')#11075
YZ6C <- subset(YZ6C, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(YZ6C, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#1083
saveRDS(YZ6C,'YZ6CQC后.rds')
head(YZ6C@meta.data)
#YZ6D------------------------
YZ6D <- Read10X(data.dir = "../../runs/CellRanger/YZ0006d/outs/multi/count/raw_feature_bc_matrix/")
YZ6D <- CreateSeuratObject(counts = YZ6D, project = "YZ6D", min.cells = 3, min.features = 200)
YZ6D <- RenameCells(YZ6D, add.cell.id = 'YZ6D')
YZ6D[["percent.mt"]] <- PercentageFeatureSet(YZ6D, pattern = "^MT-")
VlnPlot(YZ6D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(YZ6D,'YZ6DQC前.rds')#11346
YZ6D <- subset(YZ6D, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(YZ6D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#10818
saveRDS(YZ6D,'YZ6DQC后.rds')
head(YZ6D@meta.data)

#GZ7------------------------
GZ7 <- Read10X(data.dir = "../../runs/GZ0007/outs/multi/count/raw_feature_bc_matrix/")
GZ7 <- CreateSeuratObject(counts = GZ7, project = "GZ7", min.cells = 3, min.features = 200)
GZ7 <- RenameCells(GZ7, add.cell.id = 'GZ7')
GZ7[["percent.mt"]] <- PercentageFeatureSet(GZ7, pattern = "^MT-")
VlnPlot(GZ7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(GZ7,'GZ7QC前.rds')#8535
GZ7 <- subset(GZ7, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(GZ7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#8359
saveRDS(GZ7,'GZ7QC后.rds')
head(GZ7@meta.data)
#YZ7------------------------
YZ7 <- Read10X(data.dir = "../../runs/YZ0007/outs/multi/count/raw_feature_bc_matrix/")
YZ7 <- CreateSeuratObject(counts = YZ7, project = "YZ7", min.cells = 3, min.features = 200)
YZ7 <- RenameCells(YZ7, add.cell.id = 'YZ7')
YZ7[["percent.mt"]] <- PercentageFeatureSet(YZ7, pattern = "^MT-")
VlnPlot(YZ7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(YZ7,'YZ7QC前.rds')#10121
YZ7 <- subset(YZ7, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(YZ7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#9181
saveRDS(YZ7,'YZ7QC后.rds')
head(YZ7@meta.data)
#YZ7A------------------------
YZ7A <- Read10X(data.dir = "../../runs/YZ0007a/outs/multi/count/raw_feature_bc_matrix/")
YZ7A <- CreateSeuratObject(counts = YZ7A, project = "YZ7A", min.cells = 3, min.features = 200)
YZ7A <- RenameCells(YZ7A, add.cell.id = 'YZ7A')
YZ7A[["percent.mt"]] <- PercentageFeatureSet(YZ7A, pattern = "^MT-")
VlnPlot(YZ7A, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(YZ7A,'YZ7AQC前.rds')#6880
YZ7A <- subset(YZ7A, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(YZ7A, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#6880
saveRDS(YZ7A,'YZ7AQC后.rds')
head(YZ7A@meta.data)
#YZ7C------------------------
YZ7C <- Read10X(data.dir = "../../runs/CellRanger/YZ0007C/outs/multi/count/raw_feature_bc_matrix/")
YZ7C <- CreateSeuratObject(counts = YZ7C, project = "YZ7C", min.cells = 3, min.features = 200)
YZ7C <- RenameCells(YZ7C, add.cell.id = 'YZ7C')
YZ7C[["percent.mt"]] <- PercentageFeatureSet(YZ7C, pattern = "^MT-")
VlnPlot(YZ7C, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(YZ7C,'YZ7CQC前.rds')#8302
YZ7C <- subset(YZ7C, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(YZ7C, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#7972
saveRDS(YZ7C,'YZ7CQC后.rds')
head(YZ7C@meta.data)
#YZ7D------------------------
YZ7D <- Read10X(data.dir = "../../runs/CellRanger/YZ0007D/outs/multi/count/raw_feature_bc_matrix/")
YZ7D <- CreateSeuratObject(counts = YZ7D, project = "YZ7D", min.cells = 3, min.features = 200)
YZ7D <- RenameCells(YZ7D, add.cell.id = 'YZ7D')
YZ7D[["percent.mt"]] <- PercentageFeatureSet(YZ7D, pattern = "^MT-")
VlnPlot(YZ7D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(YZ7D,'YZ7DQC前.rds')#7893
YZ7D <- subset(YZ7D, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(YZ7D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#7753
saveRDS(YZ7D,'YZ7DQC后.rds')
head(YZ7D@meta.data)
#GZ9------------------------
GZ9 <- Read10X(data.dir = "../../runs/CellRanger/GZ0009/outs/multi/count/raw_feature_bc_matrix/")
GZ9 <- CreateSeuratObject(counts = GZ9, project = "GZ9", min.cells = 3, min.features = 200)
GZ9 <- RenameCells(GZ9, add.cell.id = 'GZ9')
GZ9[["percent.mt"]] <- PercentageFeatureSet(GZ9, pattern = "^MT-")
VlnPlot(GZ9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(GZ9,'GZ9QC前.rds')#6526
GZ9 <- subset(GZ9, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(GZ9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#6353
saveRDS(GZ9,'GZ9QC后.rds')
head(GZ9@meta.data)
GZ9
#YZ9------------------------
YZ9 <- Read10X(data.dir = "../../runs/CellRanger/YZ0009/outs/multi/count/raw_feature_bc_matrix/")
YZ9 <- CreateSeuratObject(counts = YZ9, project = "YZ9", min.cells = 3, min.features = 200)
YZ9 <- RenameCells(YZ9, add.cell.id = 'YZ9')
YZ9[["percent.mt"]] <- PercentageFeatureSet(YZ9, pattern = "^MT-")
VlnPlot(YZ9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(YZ9,'YZ9QC前.rds')#6660
YZ9 <- subset(YZ9, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(YZ9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#6339
saveRDS(YZ9,'YZ9QC后.rds')
head(YZ9@meta.data)
#YZ9A------------------------
YZ9A <- Read10X(data.dir = "../../runs/CellRanger/YZ0009a/outs/multi/count/raw_feature_bc_matrix/")
YZ9A <- CreateSeuratObject(counts = YZ9A, project = "YZ9A", min.cells = 3, min.features = 200)
YZ9A <- RenameCells(YZ9A, add.cell.id = 'YZ9A')
YZ9A[["percent.mt"]] <- PercentageFeatureSet(YZ9A, pattern = "^MT-")
VlnPlot(YZ9A, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(YZ9A,'YZ9AQC前.rds')#4388
YZ9A <- subset(YZ9A, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 18000 & percent.mt < 5)
VlnPlot(YZ9A, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)#4364
saveRDS(YZ9A,'YZ9AQC后.rds')
head(YZ9A@meta.data)
