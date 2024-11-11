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

sample_list <- fread("CosMx-vs-Xenium-16.txt")

for (i in 1:nrow(sample_list)) {
	# i <- 1
	sample_id <- sample_list$sample_id[i]
	
	smi <- readRDS(paste0("comparison/", as.character(sample_id), "_smi_so.RDS"))
	xen <- readRDS(paste0("comparison/", as.character(sample_id), "_xen_so.RDS"))
	
	# Filter CosMx
	use_cells <- colnames(smi)[smi$nCount_Nanostring >= 20]
	smi <- smi[, use_cells]
	
	# Filter Xenium
	use_cells <- colnames(xen)[xen$nCount_Xenium >= 5]
	xen <- xen[, use_cells]
	
	# Normalize data
	smi <- NormalizeData(smi, normalization.method = "LogNormalize", scale.factor = 10000)
	xen <- NormalizeData(xen, normalization.method = "LogNormalize", scale.factor = 10000)
	
	# Find variable features and dimensionality reduction
	smi <- FindVariableFeatures(smi, selection.method = "vst", nfeatures = 300)
	smi <- ScaleData(smi)
	smi <- RunPCA(smi)
	smi <- RunUMAP(smi, dims = 1:40)
	xen <- FindVariableFeatures(xen, selection.method = "vst", nfeatures = nrow(xen))
	xen <- ScaleData(xen)
	xen <- RunPCA(xen)
	xen <- RunUMAP(xen, dims = 1:40)
	
	saveRDS(smi, paste0("comparison/", sample_id, "_smi_so_qcnorm.RDS"))
	saveRDS(xen, paste0("comparison/", sample_id, "_xen_so_qcnorm.RDS"))
}