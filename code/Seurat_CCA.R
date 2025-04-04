library(Seurat)
library(dplyr)
library(ggplot2)

# 1. Perform batch effect correction using CCA
# Split the data into multiple lists, each corresponding to a sample
mono_list <- SplitObject(mono, split.by = "group2")

# Normalize and find variable features for each sample
mono_list <- lapply(mono_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

# Perform batch effect correction using CCA
mono_anchors <- FindIntegrationAnchors(object.list = mono_list, dims = 1:50)
mono_combined <- IntegrateData(anchorset = mono_anchors, dims = 1:50)

# Set the default assay to the integrated data
DefaultAssay(mono_combined) <- "integrated"

# 2. Output the average expression of each gene per sample
Idents(object= mono_combined)<-"sample"
av <-AverageExpression(mono_combined, assays = "RNA")
av=as.data.frame(av$RNA)
write.csv(av, file="gene_mono_sample.csv")
