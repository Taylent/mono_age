library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# 1. Clear the environment, set the output path, and working path ------------------------------------------------------

rm(list=ls())
getwd() 
dir.create("./1-QC")
setwd("./1-QC") 
getwd()
cell.all=readRDS("../cell.all_raw.rds")



# Quality Control

# Calculate gene proportions???
# 01 Calculate mitochondrial gene proportions??? 
# str_to_upper() converts to uppercase; str_to_title() converts to title case???

mito_genes=rownames(cell.all)[grep("^MT-", rownames(cell.all))] 
mito_genes 
cell.all=PercentageFeatureSet(cell.all, "^MT-", col.name = "pMT")

# 02 Calculate ribosomal gene proportions
ribo_genes=rownames(cell.all)[grep("^RP[SL]", rownames(cell.all))] 
ribo_genes 
cell.all=PercentageFeatureSet(cell.all, "^RP[SL]", col.name = "pRP") 


# Visualize Quality Control

feats <- c("nFeature_RNA", "nCount_RNA", "pMT", "pRP")
feats <- c("nFeature_RNA", "nCount_RNA")
p1=VlnPlot(cell.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
  NoLegend() 
p1
library(ggplot2)
ggsave(filename="Vlnplot1.pdf",plot=p1, height = 5, width = 8) ## Adjust size as needed???


feats <- c("pMT", "pRP")
p2=VlnPlot(cell.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, same.y.lims=T) + 
  scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend() 
p2	
ggsave(filename="Vlnplot2.pdf",plot=p2, height = 5, width = 8)


p3=FeatureScatter(cell.all, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
p3
ggsave(filename="Scatterplot.pdf",plot=p3,height = 4, width = 6)



# Filter low-quality cells and genes

## Filter criterion 1: Cells with abnormal gene counts & genes with abnormal cell counts
nFeature_lower <- 200
nFeature_upper <- 4500
nCount_lower <- 200
nCount_upper <- 45000
selected_c <- WhichCells(cell.all, expression = nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper
                         &nCount_RNA > nCount_lower & nCount_RNA < nCount_upper) 
selected_f <- rownames(cell.all)[Matrix::rowSums(cell.all@assays$RNA@counts > 0 ) > 3] # Each gene must be expressed in at least 3 cells
cell.all.filt <- subset(cell.all, features = selected_f, cells = selected_c) # subset to keep only qualified cells (filtering)
dim(cell.all) 
dim(cell.all.filt) 
table(cell.all@meta.data$orig.ident) 
table(cell.all.filt@meta.data$orig.ident) 


# Filter criterion 2: Mitochondrial genes
pMT_lower <- 0
pMT_upper <- 15
selected_mito <- WhichCells(cell.all.filt, expression = pMT < pMT_upper & pMT >= pMT_lower)

length(selected_mito)

cell.all.filt <- subset(cell.all.filt, cells = selected_mito)
dim(cell.all.filt)
