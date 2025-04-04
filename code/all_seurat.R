library(dplyr)
library(Seurat)
library(patchwork)
library(gridExtra) #grid.arrange(plot1,plot2,ncol = 2, nrow = 1)
library(DOSE)
library(topGO)
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(readr)
library(readtext)
library(cowplot)
library(harmony)

ALL <-readRDS("C:\\Users\\Administrator\\Data\\All.rds")
ALL <- NormalizeData(ALL, normalization.method = "LogNormalize", scale.factor = 10000)
dim(Liver)
ALL <- FindVariableFeatures(object = ALL, selection.method = "vst", nfeatures = 3000)
#PCA
all.genes <- rownames(ALL)
ALL <- ScaleData(ALL, features = all.genes)
ALL=RunPCA(object= ALL,npcs = 30,pc.genes=VariableFeatures(object = ALL))     
pdf(file="05.pcaHeatmap.pdf",width=10,height=20)
DimHeatmap(object = ALL, dims = 1:50, cells = 500, balanced = TRUE,nfeatures = 50,ncol=2)
dev.off()
#Harmony
ALL <- ALL %>% 
  RunHarmony("sample", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(ALL, 'harmony')
#Cluster
pcSelect=50
ALL <- FindNeighbors(object = ALL, dims = 1:pcSelect,reduction = "harmony")               
ALL <- FindClusters(object = ALL, resolution = 0.6) 
ALL <- RunUMAP(object = ALL, dims = 1:pcSelect,reduction = "harmony") 
pdf(file="06.UMAP1.pdf",width=7,height=6)
UMAPPlot(object =ALL, pt.size = 0.6, label = TRUE)
dev.off()

