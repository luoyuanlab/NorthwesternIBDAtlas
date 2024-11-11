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

source("./ReadXenium.FIX.R")
sample_list <- fread("CosMx-vs-Xenium-4.txt")

### Calculate the number of before/after QC cells using different QC thresholds
cell_counts <- list()
j <- 1
for (i in 1:nrow(sample_list)) {
	# i <- 1
	sample_id <- sample_list$sample_id[i]
	
	smi <- readRDS(paste0("comparison/", as.character(sample_id), "_smi_so.RDS"))
	xen <- readRDS(paste0("comparison/", as.character(sample_id), "_xen_so.RDS"))
	
	cell_counts[[j]] <- data.table(sample_id = sample_id, tech = "CosMx", type = "Before QC", n_cell = ncol(smi))
	cell_counts[[j+1]] <- data.table(sample_id = sample_id, tech = "Xenium", type = "Before QC", n_cell = ncol(xen))
	j <- j + 2
	
	for (threshold in seq(10, 200, 10)) {
		# Filter SMI
		use_cells <- colnames(smi)[smi$nCount_Nanostring >= threshold]
		smi <- smi[, use_cells]
	
		# Filter Xenium
		use_cells <- colnames(xen)[xen$nCount_Xenium >= threshold]
		xen <- xen[, use_cells]
		
		cell_counts[[j]] <- data.table(sample_id = sample_id, tech = "CosMx", type = paste0("QC threshold = ", as.character(threshold)), n_cell = ncol(smi))
		cell_counts[[j+1]] <- data.table(sample_id = sample_id, tech = "Xenium", type = paste0("QC threshold = ", as.character(threshold)), n_cell = ncol(xen))
		j <- j + 2
	
		smi <- NormalizeData(smi, normalization.method = "LogNormalize", scale.factor = 10000)
		smi <- FindVariableFeatures(smi, selection.method = "vst", nfeatures = 3000)
		smi <- ScaleData(smi)
		smi <- RunPCA(smi)
		smi <- RunUMAP(smi, dims = 1:40)
		
		xen <- NormalizeData(xen, normalization.method = "LogNormalize", scale.factor = 10000)
		xen <- FindVariableFeatures(xen, selection.method = "vst", nfeatures = 3000)
		xen <- ScaleData(xen)
		xen <- RunPCA(xen)
		xen <- RunUMAP(xen, dims = 1:40)
		
		saveRDS(smi, paste0("comparison_multi_thr/", sample_id, "_smi_so_qcnorm_thr", as.character(threshold), ".RDS"))
		saveRDS(xen, paste0("comparison_multi_thr/", sample_id, "_xen_so_qcnorm_thr", as.character(threshold), ".RDS"))
	}
}
cell_counts <- rbindlist(cell_counts)
fwrite(cell_counts, "comparison_multi_thr/cell_counts.txt", sep = "\t")



### Visualize the number of cells before/after QC under different thresholds

library(data.table)
library(stringr)
library(ggpubr)

both <- fread("comparison_multi_thr/cell_counts.txt")
both$type <- ifelse(both$type == "Before QC", 0, str_replace(both$type, "QC threshold = ", ""))
both$type <- as.integer(both$type)
both$sample_id <- as.character(both$sample_id)

cosmx <- both[both$tech == "CosMx",]
xen <- both[both$tech == "Xenium",]
min(both$n_cell)
max(both$n_cell)

pdf("comparison_multi_thr/cell_counts.pdf", width = 6)
ggline(cosmx, x = "type", y = "n_cell", color = "sample_id", add = "mean_se", xlab = "nCount QC threshold", ylab = "Number of cells after QC", palette = "npg", title = "CosMx")
ggline(xen, x = "type", y = "n_cell", color = "sample_id", add = "mean_se", xlab = "nCount QC threshold", ylab = "Number of cells after QC", palette = "npg", title = "Xenium")
ggline(cosmx, x = "type", y = "n_cell", color = "sample_id", add = "mean_se", xlab = "nCount QC threshold", ylab = "Number of cells after QC", palette = "npg", ylim = c(0, 40000), title = "CosMx")
ggline(xen, x = "type", y = "n_cell", color = "sample_id", add = "mean_se", xlab = "nCount QC threshold", ylab = "Number of cells after QC", palette = "npg", ylim = c(0, 40000), title = "Xenium")
dev.off()



### For downstream analysis, use nCount = 50 as the threshold for both CosMx 6K and Xenium 5K

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

source("./ReadXenium.FIX.R")

sample_list <- fread("CosMx-vs-Xenium-4.txt")
threshold <- 50

for (i in 1:nrow(sample_list)) {
	# i <- 1
	sample_id <- sample_list$sample_id[i]
	
	smi <- readRDS(paste0("comparison/", as.character(sample_id), "_smi_so.RDS"))
	xen <- readRDS(paste0("comparison/", as.character(sample_id), "_xen_so.RDS"))
	
	# Filter CosMx
	use_cells <- colnames(smi)[smi$nCount_Nanostring >= threshold]
	smi <- smi[, use_cells]

	# Filter Xenium
	use_cells <- colnames(xen)[xen$nCount_Xenium >= threshold]
	xen <- xen[, use_cells]
	
	smi <- NormalizeData(smi, normalization.method = "LogNormalize", scale.factor = 10000)
	smi <- FindVariableFeatures(smi, selection.method = "vst", nfeatures = 3000)
	smi <- ScaleData(smi)
	smi <- RunPCA(smi)
	smi <- RunUMAP(smi, dims = 1:40)
	
	xen <- NormalizeData(xen, normalization.method = "LogNormalize", scale.factor = 10000)
	xen <- FindVariableFeatures(xen, selection.method = "vst", nfeatures = 3000)
	xen <- ScaleData(xen)
	xen <- RunPCA(xen)
	xen <- RunUMAP(xen, dims = 1:40)
	
	saveRDS(smi, paste0("comparison/", sample_id, "_smi_so_qcnorm.RDS"))
	saveRDS(xen, paste0("comparison/", sample_id, "_xen_so_qcnorm.RDS"))
}
