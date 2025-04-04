library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# Cell cycle

cell.all.filt = NormalizeData(cell.all.filt)
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
cell.all.filt=CellCycleScoring(object = cell.all.filt, 
                               s.features = s.genes, 
                               g2m.features = g2m.genes, 
                               set.ident = TRUE)
p4=VlnPlot(cell.all.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", 
           ncol = 2, pt.size = 0)
p4 
ggsave(filename="Vlnplot4_cycle.pdf",plot=p4,height = 5,width = 8)
cell.all.filt@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal() 
ggsave(filename="cycle_details.pdf",height = 6,width = 8)


cell.all.filt = FindVariableFeatures(cell.all.filt, selection.method = "vst")
cell.all.filt <- ScaleData(cell.all.filt, features = rownames(cell.all.filt))
cell.all.filt <- RunPCA(cell.all.filt, features = VariableFeatures(cell.all.filt), ndims.print = 6:10, nfeatures.print = 10)

cell.all.filt <- RunPCA(cell.all.filt, features = c(s.genes, g2m.genes))
DimPlot(cell.all.filt)


if(F){
  
  cell.all.filt <- ScaleData(cell.all.filt, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(cell.all.filt))
  cell.all.filt <- RunPCA(cell.all.filt, features = VariableFeatures(cell.all.filt), nfeatures.print = 10)
  cell.all.filt <- RunPCA(cell.all.filt, features = c(s.genes, g2m.genes))
  DimPlot(cell.all.filt)
}
