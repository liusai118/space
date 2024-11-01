library(ha5d)
library(Seurat)
library(SeuratDisk)
library(rhdf5)
library(dior)
library(dplyr)
library(reticulate)
library(dior)
library(anndata)
library(reticulate)
library(dior)
Convert("AO.h5ad", dest="h5seurat", assay = "RNA",overwrite=T)
seurat_object <- LoadH5Seurat("AO.h5seurat",meta.data = FALSE, misc = FALSE)
scRNA_harmony <- NormalizeData(seurat_object) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:30)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:30) %>% FindClusters()
p1 <- DimPlot(scRNA_harmony,group.by = "CELL",raster=FALSE)
system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "Batch")})
meta <- read.csv("all_pbmcs_metadata.csv")
meta2 <- meta[meta$X %in% colnames(scRNA_harmony),]
rownames(meta2) <- meta2[,1]
meta2 <- meta2[,-1]
scRNA_harmony <- AddMetaData(scRNA_harmony,metadata = meta2)
meta2 <- meta[which(rownames(meta) %in% name),]  
metaa <- scRNA_harmony@meta.data
scRNA_harmony <- readRDS("b.rds")
meta <- read.csv("all_pbmcs_metadata.csv")
meta2 <- meta[meta$X %in% colnames(scRNA_harmony),]
rownames(meta2) <- meta2[,1]
meta2 <- meta2[,-1]
scRNA_harmony <- AddMetaData(scRNA_harmony,metadata = meta2)
metadata <- scRNA_harmony@meta.data
m2 <- metadata[!duplicated(metadata$Donor_id),]
write.csv(m2,"m2.csv")
system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "Batch")})
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:30)
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:30)
options(future.globals.maxSize = 2 * 1024^3)  # 2 GB
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:30) %>% FindClusters()
all.markers  <- FindAllMarkers(scRNA_harmony, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.75)
write.csv(all.markers,"mark.csv")
meta_need <- scRNA_harmony@meta.data
DimPlot(scRNA_harmony,label = T,group.by = "seurat_clusters")
###############################
new.cluster.ids <- c(
  "0"="CD4+ T cells",
  "1"="CD4+ T cells",
  "2"="Natural Killer Cell",
  "3"="CD14+ Monocyte",
  "4"="CD8+ T cells",
  "5"="CD8+ T cells",
  "6"="CD14+ Monocyte",
  "7"="CD4+ T cells",
  "8"="CD8+ T cells",
  "9"="CD16+ Monocyte",
  "10"="CD4+ T cells",
  "11"="Other T Cell",
  "12"="B Cell",
  "13"="B Cell",
  "14"="CD4+ T cells",
  "15"="Other T Cell",
  "16"="CD4+ T cells",
  "17"="Natural Killer Cell",
  "18"="other",
  "19"="CD8+ T cells",
  "20"="CD4+ T cells",
  "21"="Dendritic Cell",
  "22"="other",
  "23"="CD4+ T cells",
  "24"="CD4+ T cells",
  "25"="CD4+ T cells"
)
Idents(scRNA_harmony) <- "seurat_clusters"
scRNA_harmony <- RenameIdents(scRNA_harmony, new.cluster.ids)   
scRNA_harmony$celltype1 <- scRNA_harmony@active.ident
DimPlot(scRNA_harmony, group.by = "celltype1",label= T,reduction = "umap")
######################
new.cluster.ids <- c(
  "0"="CD4 T",
  "1"="CD4 T",
  "2"="NK",
  "3"="Mono",
  "4"="CD8 T",
  "5"="CD8 T",
  "6"="Mono",
  "7"="CD4 T",
  "8"="CD8 T",
  "9"="Mono",
  "10"="CD4 T",
  "11"="other T",
  "12"="B Cell",
  "13"="B Cell",
  "14"="CD4 T",
  "15"="other T",
  "16"="CD4 T",
  "17"="NK",
  "18"="other",
  "19"="CD8 T",
  "20"="CD4 T",
  "21"="DC",
  "22"="other",
  "23"="CD4 T",
  "24"="CD4 T",
  "25"="CD4 T"
)
Idents(scRNA_harmony) <- "seurat_clusters"
scRNA_harmony <- RenameIdents(scRNA_harmony, new.cluster.ids)   
scRNA_harmony$celltype2 <- scRNA_harmony@active.ident
DimPlot(scRNA_harmony, group.by = "celltype2",label= T,reduction = "umap")
###############################
Idents(scRNA_harmony) <- "Age_group"
degs = lapply(unique(scRNA_harmony$celltype1), function(x){
  FindMarkers(scRNA_harmony[,scRNA_harmony$celltype1==x],ident.1 = 'E',
              ident.2 = 'A',verbose = T, test.use = 'wilcox',min.pct = 0.1)
})
deg_total <- data.frame()
a <- unique(scRNA_harmony$celltype1)
for (i in seq_along(a)){
  degs2 <- as.data.frame(degs[[i]])
  degs2$cluster <- c(rep(a[i],nrow(degs2)))
  degs2$gene <- rownames(degs2)
  deg_total<- rbind(deg_total,degs2)
}


write.csv(deg_total,"deg_aging.csv")
