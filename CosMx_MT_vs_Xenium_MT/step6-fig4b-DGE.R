library(data.table)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)

subsetList <- function(myList, elementNames) {
	sapply(elementNames, FUN=function(x) myList[[x]])
}

sample_list <- fread("sample_list_16.txt")
sample_list$category <- paste0(sample_list$tissue_type, ", ", sample_list$disease_state, " ", sample_list$disease)

### Merging of CosMx data
so_list <- list()
for (i in 1:nrow(sample_list)) {
	# i <- 1
	sample_id <- sample_list$sample_id[i]
	so <- readRDS(paste0("comparison/", as.character(sample_id), "_smi_so_qcnorm.RDS"))
	so@meta.data$sample_id <- sample_id
	so@meta.data$block <- sample_list$block[i]
	so@meta.data$tissue_type <- sample_list$tissue_type[i]
	so@meta.data$disease <- sample_list$disease[i]
	so@meta.data$disease_state <- sample_list$disease_state[i]
	so@meta.data$category <- sample_list$category[i]
	
	so <- SCTransform(so, assay = "Nanostring", return.only.var.genes = FALSE, clip.range = c(-10, 10), verbose = FALSE)
	so_list[[sample_id]] <- so
}
so <- merge(x = so_list[[1]], y = c(subsetList(so_list, seq(2, length(so_list)))), add.cell.ids = as.character(1:16), project = "CosMx")
so <- PrepSCTFindMarkers(object = so)
VariableFeatures(so) <- rownames(so) # set all the genes to be variable because we want to use all the genes
so <- RunPCA(so, assay = "SCT", verbose = FALSE)
so <- RunUMAP(so, reduction = "pca", dims = 1:40, verbose = FALSE)
so <- FindNeighbors(so, reduction = "pca", dims = 1:40, verbose = FALSE)

### CosMx sample category DGE
annotations <- so[["category"]]$category
names(annotations) <- rownames(so[["category"]])
Idents(so) <- annotations
category_DGE <- FindAllMarkers(so, assay = "SCT")
fwrite(category_DGE, "comparison/cosmx_merged_category_DGE.txt", sep = "\t")

#####

### Merging of Xenium data
so_list <- list()
for (i in 1:nrow(sample_list)) {
	# i <- 1
	sample_id <- sample_list$sample_id[i]
	so <- readRDS(paste0("comparison/", as.character(sample_id), "_xen_so_qcnorm.RDS"))
	so@meta.data$sample_id <- sample_id
	so@meta.data$block <- sample_list$block[i]
	so@meta.data$tissue_type <- sample_list$tissue_type[i]
	so@meta.data$disease <- sample_list$disease[i]
	so@meta.data$disease_state <- sample_list$disease_state[i]
	so@meta.data$category <- sample_list$category[i]
	
	so <- SCTransform(so, assay = "Xenium", return.only.var.genes = FALSE, clip.range = c(-10, 10), verbose = FALSE)
	so_list[[sample_id]] <- so
}

so <- merge(x = so_list[[1]], y = c(subsetList(so_list, seq(2, length(so_list)))), add.cell.ids = as.character(1:16), project = "Xenium")
so <- PrepSCTFindMarkers(object = so)
VariableFeatures(so) <- rownames(so) # set all the genes to be variable because we want to use all the genes
so <- RunPCA(so, assay = "SCT", verbose = FALSE)
so <- RunUMAP(so, reduction = "pca", dims = 1:40, verbose = FALSE)
so <- FindNeighbors(so, reduction = "pca", dims = 1:40, verbose = FALSE)

### Xenium sample category DGE
annotations <- so[["category"]]$category
names(annotations) <- rownames(so[["category"]])
Idents(so) <- annotations
category_DGE <- FindAllMarkers(so, assay = "SCT")
fwrite(category_DGE, "comparison/xenium_merged_category_DGE.txt", sep = "\t")