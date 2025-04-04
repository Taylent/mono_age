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
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(AUCell)
library(GSEABase)
cluAll <- readRDS("C:/Users/Administrator/Documents/all_mono2023.12.4.rds")
setwd("C:/Users/Administrator/Documents/allnewmono/")
c5 <- read.gmt("ALLNEWMONO.gmt")
head(c5$term)
head(c5$gene)
geneSets <- lapply(unique(c5$term), function(x){print(x);c5$gene[c5$term == x]})
names(geneSets) <- unique(c5$term)
cluAll<- subset(cluAll, ident=c("c1",'c2','c3','c4','c5','c7','c8'))
cells_rankings1 <- AUCell_buildRankings(cluAll@assays$RNA@data)
## Select the top 5% highly expressed genes to calculate AUC
cells_AUC1 <- AUCell_calcAUC(geneSets, cells_rankings1, aucMaxRank=nrow(cells_rankings1)*0.05)
geneSet <- "AUC01"
aucs <- as.numeric(getAUC(cells_AUC1)[geneSet, ])
cluAll$AUC01 <- aucs
geneSet <- "AUC02"
aucs <- as.numeric(getAUC(cells_AUC1)[geneSet, ])
cluAll$AUC02 <- aucs
geneSet <- "AUC03"
aucs <- as.numeric(getAUC(cells_AUC1)[geneSet, ])
# Middle part omitted
geneSet <- "AUC16"
aucs <- as.numeric(getAUC(cells_AUC1)[geneSet, ])
cluAll$AUC70 <- aucs
df<- data.frame(cluAll@meta.data, cluAll@reductions$umap@cell.embeddings)
class_avg <- df %>%
  group_by(seurat_clusters) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
write.csv(df,'mergeIM_AUC_10.csv')
p <-DotPlot(object = cluAll, features = c('AUC01','AUC02','AUC03','AUC04','AUC05','AUC06','AUC07','AUC08','AUC09','AUC10','AUC11','AUC12','AUC13','AUC14','AUC15','AUC16')) + coord_flip()
write.csv(p[["data"]],'allmonoaucnew2.csv')

