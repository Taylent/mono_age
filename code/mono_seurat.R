library(dplyr)
library(Seurat)
library(patchwork)
library(gridExtra) 
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
MONO <- FindVariableFeatures(object = MONO, selection.method = "vst", nfeatures = 3000)
# Output the feature variance plot
# PCA
all.genes <- rownames(MONO)
MONO <- ScaleData(MONO, features = all.genes)
MONO=RunPCA(object= MONO,npcs = 50,pc.genes=VariableFeatures(object = MONO))     # PCA analysis
pdf(file="all.cellcycle.age.pdf",width=30,height=10)
DimHeatmap(object = MONO, dims = 1:50, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()
MONO <- MONO %>%
  RunHarmony("sample", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(MONO, 'harmony')
pcSelect=30
MONO <- FindNeighbors(object = MONO, dims = 1:pcSelect,reduction = "harmony")                # Calculate adjacency distance
MONO <- FindClusters(object = MONO, resolution = 0.6)
MONO <- RunUMAP(object = MONO, dims = 1:pcSelect,reduction = "harmony")
UMAPPlot(object =MONO, pt.size= 0.6, label = TRUE)
# Check for contamination of other cell types
FeaturePlot(object = MONO, features = c('CD14','FCGR3A','CD3E','NKG7','CD19','ITGAX','CD3D'),cols = c("ivory2","ivory","tomato","tomato1" ,"tomato2","tomato3","tomato4"))
# Remove contaminated cells
MONO<- subset(MONO, ident=c("0",'1','2','3','4','5','6','7','9','11','13','14','16','17' ))
UMAPPlot(object =MONO, pt.size= 0.6, label = TRUE)
# Re-cluster
MONO <- FindClusters(object = MONO, resolution = 0.8)
MONO <- RunUMAP(object = MONO, dims = 1:pcSelect,reduction = "harmony")
UMAPPlot(object =MONO, pt.size= 0.6, label = TRUE)
new.cluster.ids <- c( 'c1','c2','c3','c4','c5','c6','c7','c6','c8','c6','c6','c4','c9')
names(new.cluster.ids) <- levels(MONO)
MONO <- RenameIdents(MONO, new.cluster.ids)
UMAPPlot(object =MONO, pt.size= 0.6, label = TRUE)
UMAPPlot(object =MONO, pt.size= 0.6, label = TRUE,split.by='group')
ClusterPercentagessample<-prop.table(table(Idents(MONO), MONO$sample), margin = 2)
write.csv(ClusterPercentagessample, file="all.mono.csv")
saveRDS(MONO, file = "all_mono.rds")
