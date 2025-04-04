# 导入Seurat包
library(Seurat)
library(harmony)
library(dplyr)


Seurat_object <- readRDS("all_mono.RDS")
Seurat_object@meta.data$sample <- Idents(Seurat_object)

# 计算线粒体基因比例
# 该操作会在数据里面增加一列叫做percent.mt
Seurat_object[["percent.mt"]] <- PercentageFeatureSet(Seurat_object, pattern = "^MT-")
VlnPlot(Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 正式筛选，筛选的是细胞，最终细胞减少
Seurat_object <- subset(Seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 15)

Seurat_object <- Seurat_object %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>%
  ScaleData() %>%
  RunPCA()

ElbowPlot(Seurat_object)

DimPlot(Seurat_object, reduction = "pca")

# 去除批次效应
Seurat_object <- RunHarmony(Seurat_object, reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony")
harmony_embeddings <- Embeddings(Seurat_object, 'harmony')


Seurat_object <- RunUMAP(Seurat_object, reduction = "harmony", dims = 1:30, reduction.name = "umap")


# 转置表达量数据，使得样本为行，基因为列
scale_data <- t(as.matrix(Seurat_object@assays$RNA@scale.data))

# 创建一个数据框（DataFrame）用于保存到 CSV
scale_data_df <- as.data.frame(scale_data)

# 保存到 CSV，同时保留行名和列名
write.csv(scale_data_df, file = "scale_data.csv")

# 提示信息
cat("Expression data saved to scale_data.csv\n")


# 创建一个数据框（DataFrame）用于保存到 CSV
harmony_embeddings_df <- as.data.frame(harmony_embeddings)

# 拼接列名
new_rownames <- paste(Seurat_object@meta.data$sample, rownames(harmony_embeddings_df), sep = "_")

# 设置新列名
rownames(harmony_embeddings_df) <- new_rownames

# 保存到 CSV，同时保留行名和列名
write.csv(harmony_embeddings_df, file = "harmony_embeddings.csv")

# 提示信息
cat("Expression data saved to harmony_embeddings.csv\n")