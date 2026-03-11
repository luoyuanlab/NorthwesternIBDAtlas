library(rjson)
library(data.table)
library(Matrix)
library(Seurat)
library(viridis)
library(ggplot2)
options(Seurat.object.assay.version = "v5") # Ensure using Seurat V5 assays
library(patchwork)
library(progressr)

grepr <- function(pattern, v){
  return(v[grep(pattern, v)])
}

sample_list <- fread("~/CosMx-vs-Xenium-16.txt")

for (i in 1:nrow(sample_list)) {
  # i <- 1
  sample_id <- sample_list$sample_id[i]
  
  smi <- readRDS(paste0("~/comparison/", as.character(sample_id), "_smi_so.RDS"))
  xen <- readRDS(paste0("~/comparison/", as.character(sample_id), "_xen_so.RDS"))
  xen_unFilter <- readRDS(paste0("~/comparison_Cenfu/", as.character(sample_id), "_xen_so_unQVfiltered.RDS"))

  # Normalize data
  smi <- NormalizeData(smi, normalization.method = "LogNormalize", scale.factor = 10000)
  xen <- NormalizeData(xen, normalization.method = "LogNormalize", scale.factor = 10000)
  xen_unFilter <- NormalizeData(xen_unFilter, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Find variable features and dimensionality reduction
  smi <- FindVariableFeatures(smi, selection.method = "vst", nfeatures = 300)
  smi <- ScaleData(smi)
  smi <- RunPCA(smi)
  smi <- RunUMAP(smi, dims = 1:40)
  xen <- FindVariableFeatures(xen, selection.method = "vst", nfeatures = nrow(xen))
  xen <- ScaleData(xen)
  xen <- RunPCA(xen)
  xen <- RunUMAP(xen, dims = 1:40)
  xen_unFilter <- FindVariableFeatures(xen_unFilter, selection.method = "vst", nfeatures = nrow(xen_unFilter))
  xen_unFilter <- ScaleData(xen_unFilter)
  xen_unFilter <- RunPCA(xen_unFilter)
  xen_unFilter <- RunUMAP(xen_unFilter, dims = 1:40)
  
  saveRDS(smi, paste0("~/", sample_id, "_smi_so_norm.RDS"))
  saveRDS(xen, paste0("~/", sample_id, "_xen_so_norm.RDS"))
  saveRDS(xen_unFilter, paste0("~/", sample_id, "_xen_so_unQVfiltered_norm.RDS"))
  
}