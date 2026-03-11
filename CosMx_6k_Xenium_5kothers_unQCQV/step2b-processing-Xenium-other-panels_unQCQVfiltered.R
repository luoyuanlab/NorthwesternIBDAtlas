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


sample_list <- fread("~/CosMx-vs-Xenium-4.txt")


for (i in 1:nrow(sample_list)) {
	# i <- 1
	sample_id <- sample_list$sample_id[i]
	
	
	xen <- readRDS(paste0("~/", as.character(sample_id), "_xen_377_with_seg_so.RDS"))
	xen <- NormalizeData(xen, normalization.method = "LogNormalize", scale.factor = 10000)
	xen <- FindVariableFeatures(xen, selection.method = "vst", nfeatures = 3000)
	xen <- ScaleData(xen)
	xen <- RunPCA(xen)
	xen <- RunUMAP(xen, dims = 1:40)
	saveRDS(xen, paste0("~/", sample_id, "_xen_377_with_seg_so_norm.RDS"))
	rm(xen)
	gc()
	#QV unfiltered
	xen <- readRDS(paste0("~/", as.character(sample_id), "_xen_377_with_seg_so_unQVfiltered.RDS"))
	xen <- NormalizeData(xen, normalization.method = "LogNormalize", scale.factor = 10000)
	xen <- FindVariableFeatures(xen, selection.method = "vst", nfeatures = 3000)
	xen <- ScaleData(xen)
	xen <- RunPCA(xen)
	xen <- RunUMAP(xen, dims = 1:40)
	saveRDS(xen, paste0("~/", sample_id, "_xen_377_with_seg_so_unQVfiltered_norm.RDS"))
	rm(xen)
	gc()
	
	print("checkpoint with seg")
	#----
	
	xen <- readRDS(paste0("~/comparison/", as.character(sample_id), "_xen_377_no_seg_so.RDS"))
	xen <- NormalizeData(xen, normalization.method = "LogNormalize", scale.factor = 10000)
	xen <- FindVariableFeatures(xen, selection.method = "vst", nfeatures = 3000)
	xen <- ScaleData(xen)
	xen <- RunPCA(xen)
	xen <- RunUMAP(xen, dims = 1:40)
	saveRDS(xen, paste0("~/", sample_id, "_xen_377_no_seg_so_norm.RDS"))
	rm(xen)
	gc()
	#qv unfiltered
	xen <- readRDS(paste0("~/", as.character(sample_id), "_xen_377_no_seg_so_unQVfiltered.RDS"))
	xen <- NormalizeData(xen, normalization.method = "LogNormalize", scale.factor = 10000)
	xen <- FindVariableFeatures(xen, selection.method = "vst", nfeatures = 3000)
	xen <- ScaleData(xen)
	xen <- RunPCA(xen)
	xen <- RunUMAP(xen, dims = 1:40)
	saveRDS(xen, paste0("~/", sample_id, "_xen_377_no_seg_so_unQVfiltered_norm.RDS"))
	rm(xen)
	gc()
	
	print("checkpoint no seg")
	
	#---
	
	xen <- readRDS(paste0("~/comparison/", as.character(sample_id), "_xen_colon_with_seg_so.RDS"))
	xen <- NormalizeData(xen, normalization.method = "LogNormalize", scale.factor = 10000)
	xen <- FindVariableFeatures(xen, selection.method = "vst", nfeatures = 3000)
	xen <- ScaleData(xen)
	xen <- RunPCA(xen)
	xen <- RunUMAP(xen, dims = 1:40)
	saveRDS(xen, paste0("~/", sample_id, "_xen_colon_with_seg_so_norm.RDS"))
	rm(xen)
	gc()
	#qv unfiltered
	xen <- readRDS(paste0("~/", as.character(sample_id), "_xen_colon_with_seg_so_unQVfiltered.RDS"))
	xen <- NormalizeData(xen, normalization.method = "LogNormalize", scale.factor = 10000)
	xen <- FindVariableFeatures(xen, selection.method = "vst", nfeatures = 3000)
	xen <- ScaleData(xen)
	xen <- RunPCA(xen)
	xen <- RunUMAP(xen, dims = 1:40)
	saveRDS(xen, paste0("~/", sample_id, "_xen_colon_with_seg_so_unQVfiltered_norm.RDS"))
	rm(xen)
	gc()
}
