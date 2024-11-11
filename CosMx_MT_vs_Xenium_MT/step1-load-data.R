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
source("LoadNanostring.FIX.R")

sample_list <- fread("CosMx-vs-Xenium-16.txt")

for (i in 1:nrow(sample_list)) {
	# i <- 1
	sample_id <- sample_list$sample_id[i]
	panel <- list(smi = list(), xen = list())
	
	### Xenium
	data_path <- sample_list$xenium_folder[i]
	xen_panel <- file.path(data_path, 'gene_panel.json')
	xen_panel <- fromJSON(file = xen_panel)
	xen_panel <- sapply(1:length(xen_panel$payload$targets), function(i) xen_panel$payload$targets[[i]]$type$data$name)
	
	# Determine the negprobes and RNA targets
	panel$xen$negprobes <- grepr('NegControlProbe', xen_panel)
	panel$xen$rna <- xen_panel[!(xen_panel %in% panel$xen$negprobes)]
	
	# Read Xenium data and metadata file
	xen <- LoadXenium(data_path)
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
	
	###########################
	
	### CosMx
	flow_cell_name <- sample_list$cosmx_flow_cell[i]
	fov_start <- sample_list$fov_start[i]
	fov_end <- sample_list$fov_end[i]
	data_path2 <- paste0(flow_cell_name, "/flat")
	
	# Create temporary folder including sample-specific FoVs only
	dir.create("tmp2")
	tmp_exprMat <- fread(paste0(data_path2, "/", flow_cell_name, "_exprMat_file.csv"))
	tmp_fovPos <- fread(paste0(data_path2, "/", flow_cell_name, "_fov_positions_file.csv"))
	tmp_metadata <- fread(paste0(data_path2, "/", flow_cell_name, "_metadata_file.csv"))
	tmp_polygon <- fread(paste0(data_path2, "/", flow_cell_name, "-polygons.csv"))
	tmp_txfile <- fread(paste0(data_path2, "/", flow_cell_name, "_tx_file.csv"))
	tmp_exprMat <- tmp_exprMat[tmp_exprMat$fov <= fov_end & tmp_exprMat$fov >= fov_start,]
	fwrite(tmp_exprMat, paste0("tmp2/", flow_cell_name, "_exprMat_file.csv"), sep = ",")
	tmp_fovPos <- tmp_fovPos[tmp_fovPos$FOV <= fov_end & tmp_fovPos$FOV >= fov_start,]
	fwrite(tmp_fovPos, paste0("tmp2/", flow_cell_name, "_fov_positions_file.csv"), sep = ",")
	tmp_metadata <- tmp_metadata[tmp_metadata$fov <= fov_end & tmp_metadata$fov >= fov_start,]
	fwrite(tmp_metadata, paste0("tmp2/", flow_cell_name, "_metadata_file.csv"), sep = ",")
	tmp_polygon <- tmp_polygon[tmp_polygon$fov <= fov_end & tmp_polygon$fov >= fov_start,]
	fwrite(tmp_polygon, paste0("tmp2/", flow_cell_name, "-polygons.csv"), sep = ",")
	tmp_txfile <- tmp_txfile[tmp_txfile$fov <= fov_end & tmp_txfile$fov >= fov_start,]
	fwrite(tmp_txfile, paste0("tmp2/", flow_cell_name, "_tx_file.csv"), sep = ",")
	smi <- LoadNanostring.FIX("tmp2", fov = 'smi')
	
	panel$smi$negprobes <- grepr('Negative', rownames(smi))
	panel$smi$falsecode <- grepr('SystemControl', rownames(smi))
	exclude_vars <- c(grepr('Negative', rownames(smi)), grepr('SystemControl', rownames(smi)))
	panel$smi$rna <- rownames(smi)[!(rownames(smi) %in% exclude_vars)]
	
	# Extract negprobes and select assay data
	smi_negs <- smi[panel$smi$negprobes,]
	smi_negs <- smi_negs@assays$Nanostring
	
	# Extract falsecodes and select assay data
	smi_falsecode <- smi[panel$smi$falsecode,]
	smi_falsecode <- smi_falsecode@assays$Nanostring
	
	# Keep only RNA targets in the assay
	smi <- smi[panel$smi$rna,]
	Key(object = smi[["Nanostring"]]) <- "rna_"
	
	# Add negprobes back as the second assay
	smi_negs@key <- 'negprobes_'
	smi[['negprobes']] <- smi_negs
	
	# Add falsecodes back as the third assay
	smi_falsecode@key <- 'falsecode_'
	smi[['falsecode']] <- smi_falsecode
	
	# Add metadata
	metaSMI <- read.csv(paste0(data_path2, "/", flow_cell_name, "_metadata_file.csv"))
	# Select cells in the seurat object
	rownames(metaSMI) <- paste0(metaSMI$cell_ID, "_", metaSMI$fov)
	smi@meta.data <- cbind(smi@meta.data, metaSMI[rownames(smi@meta.data),])
	
	saveRDS(panel, paste0("comparison/", sample_id, "_panel.RDS"))
	saveRDS(smi, paste0("comparison/", sample_id, "_smi_so.RDS"))
	saveRDS(xen, paste0("comparison/", sample_id, "_xen_so.RDS"))
	
	unlink("tmp2", recursive = TRUE)
	gc()
}
