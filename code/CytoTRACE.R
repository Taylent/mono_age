
rm(list = ls())
# devtools::install_local("CytoTRACE.tar.gz")
# library(reticulate)
# conda_create("cytoTRACE",python_version = '3.7')
# use_condaenv("cytoTRACE")
# conda_install("cytoTRACE", "numpy") 
library(Seurat)
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(reshape2)
library(ggsci)
library(ggpubr)
library(data.table)
library(clusterProfiler)
library(CytoTRACE)
set.seed(12345)
dir.create("cytoTRACE")
setwd("./cytoTRACE")

#load data

Idents(mono)=mono$ cell_type # Assign cell type to the Idents attribute; the cell types have already been annotated
DimPlot(mono,split.by = "group")
mono=SplitObject(mono, split.by = "group")
seurat<-mono
table(seurat$group)  
table(seurat$cell_type)
seurat
rm(mono)
table(Idents(seurat))


# Select the cell types for which to construct the cell differentiation trajectory (subset extracts the names of interest, which have already been prepared above)

exp <- seurat@assays$RNA@counts
exp <- as.matrix(exp)
phe <- seurat$cell_type
phe = as.character(phe)
names(phe) <- rownames(seurat@meta.data)
emb <- seurat@reductions$umap@cell.embeddings

results <- CytoTRACE(mat = exp)

plotCytoGenes(results, numOfGenes = 30)
plotCytoTRACE(results, phenotype = phe, emb=emb)

save(results, file = 'CytoTRACE_results.Rdata')

