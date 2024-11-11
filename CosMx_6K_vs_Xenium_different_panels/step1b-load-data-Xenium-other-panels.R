library(rjson)
library(data.table)
library(Matrix)
library(Seurat)
library(viridis)
library(ggplot2)
options(Seurat.object.assay.version = "v5") # Ensure using Seurat V5 assays

grepr <- function(pattern, v){
	return(v[grep(pattern, v)])
}

library(progressr)
source("ReadXenium.FIX.R")

sample_list <- fread("CosMx-vs-Xenium-4.txt")

### Xenium multi-tissue, with segmentation

for (i in 1:nrow(sample_list)) {
	# i <- 4
	sample_id <- sample_list$sample_id[i]
	
	data_path <- sample_list$xenium_377_with_seg[i]
	xen <- LoadXenium_FIX(data_path)
	
	# Read metadata file
	md <- fread(file.path(data_path, 'cells.csv.gz'))
	md <- as.data.frame(md)
	rownames(md) <- md$cell_id
	
	# Add cell area information to the Seurat object
	xen$cell_area <- md[colnames(xen), 'cell_area']
	xen$nucleus_area <- md[colnames(xen), 'nucleus_area']

	# Add cell positions to Seurat object
	xen$x_centroid <- md[colnames(xen), 'x_centroid']
	xen$y_centroid <- md[colnames(xen), 'y_centroid']

	# Read cluster assignments and add to the Seurat object
	clust_methods <- c('gene_expression_graphclust', paste0('gene_expression_kmeans_', 4:10, '_clusters'))
	for(clust_method_i in clust_methods){
		clusts <- fread(file.path(data_path, 'analysis', 'clustering', clust_method_i, 'clusters.csv'))
		clusts <- as.data.frame(clusts)
		rownames(clusts) <- clusts$Barcode
		xen[[clust_method_i]] <- as.factor(clusts[colnames(xen), 'Cluster'])
	}
	
	saveRDS(xen, paste0("comparison/", sample_id, "_xen_377_with_seg_so.RDS"))
	gc()
}

### Xenium multi-tissue, no segmentation

for (i in 1:nrow(sample_list)) {
	# i <- 1
	sample_id <- sample_list$sample_id[i]
	
	data_path <- sample_list$xenium_377_no_seg[i]
	xen <- LoadXenium_FIX(data_path)
	
	# Read metadata file
	md <- fread(file.path(data_path, 'cells.csv.gz'))
	md <- as.data.frame(md)
	rownames(md) <- md$cell_id

	# Add cell area information to the Seurat object
	xen$cell_area <- md[colnames(xen), 'cell_area']
	xen$nucleus_area <- md[colnames(xen), 'nucleus_area']

	# Add cell positions to Seurat object
	xen$x_centroid <- md[colnames(xen), 'x_centroid']
	xen$y_centroid <- md[colnames(xen), 'y_centroid']

	# Read cluster assignments and add to the Seurat object
	clust_methods <- c('gene_expression_graphclust', paste0('gene_expression_kmeans_', 4:10, '_clusters'))
	for(clust_method_i in clust_methods){
		clusts <- fread(file.path(data_path, 'analysis', 'clustering', clust_method_i, 'clusters.csv'))
		clusts <- as.data.frame(clusts)
		rownames(clusts) <- clusts$Barcode
		xen[[clust_method_i]] <- as.factor(clusts[colnames(xen), 'Cluster'])
	}
	
	saveRDS(xen, paste0("comparison/", sample_id, "_xen_377_no_seg_so.RDS"))
	gc()
}

### Xenium colon-specific, 322 genes

for (i in 1:nrow(sample_list)) {
	# i <- 1
	sample_id <- sample_list$sample_id[i]
	
	data_path <- sample_list$xenium_folder_colon[i]
	xen <- LoadXenium_FIX(data_path)
	
	# Read metadata file
	md <- fread(file.path(data_path, 'cells.csv.gz'))
	md <- as.data.frame(md)
	rownames(md) <- md$cell_id

	# Add cell area information to the Seurat object
	xen$cell_area <- md[colnames(xen), 'cell_area']
	xen$nucleus_area <- md[colnames(xen), 'nucleus_area']

	# Add cell positions to Seurat object
	xen$x_centroid <- md[colnames(xen), 'x_centroid']
	xen$y_centroid <- md[colnames(xen), 'y_centroid']

	# Read cluster assignments and add to the Seurat object
	clust_methods <- c('gene_expression_graphclust', paste0('gene_expression_kmeans_', 4:10, '_clusters'))
	for(clust_method_i in clust_methods){
		clusts <- fread(file.path(data_path, 'analysis', 'clustering', clust_method_i, 'clusters.csv'))
		clusts <- as.data.frame(clusts)
		rownames(clusts) <- clusts$Barcode
		xen[[clust_method_i]] <- as.factor(clusts[colnames(xen), 'Cluster'])
	}
	
	saveRDS(xen, paste0("comparison/", sample_id, "_xen_colon_with_seg_so.RDS"))
	gc()
}
