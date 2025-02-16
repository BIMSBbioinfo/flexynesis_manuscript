# download and prepare data for flexynesis 

# 1. scrnaseq data with cell type labels 
# See https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis

library(Seurat)
# devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(dplyr)
library(ggplot2)
# InstallData("bmcite")
bm <- LoadData(ds = "bmcite")
# randomly select 5K cells per train/test 
N <- 1000
in_train <- sample(1:ncol(bm), N)
in_test <- sample(setdiff(1:ncol(bm), in_train), N)
bm_train <- bm[,in_train]
bm_test <- bm[,in_test]   

# separately plot for RNA
DefaultAssay(bm_train) <- 'RNA'
bm_train <- NormalizeData(bm_train) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
bm_train <- RunUMAP(bm_train, reduction = 'pca', dims = 1:30)
p1 <- DimPlot(bm_train, group.by = 'celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) 

DefaultAssay(bm_test) <- 'RNA'
bm_test <- NormalizeData(bm_test) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
bm_test <- RunUMAP(bm_test, reduction = 'pca', dims = 1:30)
p2 <- DimPlot(bm_test, group.by = 'celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) 

DefaultAssay(bm_train) <- 'ADT'
bm_train <- NormalizeData(bm_train) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = 25)
bm_train <- RunUMAP(bm_train, reduction = 'pca', dims = 1:10)
p3 <- DimPlot(bm_train, group.by = 'celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) 

DefaultAssay(bm_test) <- 'ADT'
bm_test <- NormalizeData(bm_test) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = 25)
bm_test <- RunUMAP(bm_test, reduction = 'pca', dims = 1:10)
p4 <- DimPlot(bm_test, group.by = 'celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) 

plots <- list('rna_train' = p1, 'rna_test' = p2, 'adt_train' = p3, 'adt_test' = p4)
p <- cowplot::plot_grid(plotlist = plots, labels = names(plots))
ggsave("bm_umap_plots.pdf", plot = p, height = 12, width = 20)
# Export train/test data 
# obj: seurat object
# split: train/test
# datatype: e.g. RNA/ADT
export_data <- function(obj, split, datatypes) {
  # print raw counts
  if(!dir.exists(split)) { dir.create(split)}
  for (datatype in datatypes) {
    M <- as.matrix(obj[[datatype]]@counts)
    write.table(M, file.path(split, paste0(datatype, ".csv")), 
                sep = ',', row.names = T, col.names = T, quote = F)
    
  }
  # print cell labels
  df <- obj@meta.data
  colnames(df) <- gsub("\\.", "_", colnames(df))
  write.table(df, file.path(split, paste0("clin.csv")),
              sep = ',', row.names = T, col.names = T, quote = F)
}

# export raw count tables and cell metadata to folder 
export_data(bm_train, 'train', c('RNA', 'ADT'))
export_data(bm_test, 'test', c('RNA', 'ADT'))


