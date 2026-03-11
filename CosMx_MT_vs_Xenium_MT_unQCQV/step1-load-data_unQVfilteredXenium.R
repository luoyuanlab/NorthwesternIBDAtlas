# cd /share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex
# conda activate R4
# R


library(rjson)
library(data.table)
library(Matrix)
library(Seurat)
#library(viridis)
library(dplyr)
library(ggplot2)
library(arrow)
options(Seurat.object.assay.version = "v5") # Ensure using Seurat V5 assays



sample_list <- fread("~/CosMx-vs-Xenium-16.txt") # Under /share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex/

for (i in 1:nrow(sample_list)) {
	#i <- 1
	sample_id <- sample_list$sample_id[i]
	print(paste0("working on ",sample_list$xenium_folder[i]))
	#panel <- list(smi = list(), xen = list())
	
	### Xenium
	data_path <- paste0("~/", sample_list$xenium_folder[i])
	parquet_path <- file.path(data_path, "transcripts.parquet")
	
	#read in the parquet unfiltered transcripts file
	parquet_file <- arrow::read_parquet(parquet_path, as_data_frame = TRUE)
	
	#remove the unassigned transcripts and non-gene probes, note there is a QV column in the df, 
	#we don't use the default creteria QV >=20 to filter it 
	parquet_file_assign <- parquet_file %>% filter(cell_id != "UNASSIGNED") #%>% filter(is_gene == TRUE) only for 5k
	
	#need to remove all the negprobe and unassiged gene in the features
	bad_pat <- "(?i)^(NegControl(Codeword|Probe)|UnassignedCodeword)"  # case-insensitive
	keep <- unique(parquet_file_assign$feature_name)[!grepl(bad_pat, unique(parquet_file_assign$feature_name))]
	parquet_file_assign <- parquet_file_assign %>% filter(feature_name %in% keep)
	
	# Create a factor for cell_id to ensure all cells are included, even those with zero counts later
	all_cell_ids <- unique(parquet_file_assign$cell_id) # Use all detected cell barcodes for a truly unfiltered matrix
	
	# Build a sparse matrix: rows are genes, columns are cell barcodes
	count_data <- parquet_file_assign %>%
	  count(feature_name, cell_id)
	
	# Now, convert feature_name and cell_id to factors with levels
	count_data$i <- as.integer(factor(count_data$feature_name, levels = unique(parquet_file_assign$feature_name)))
	count_data$j <- as.integer(factor(count_data$cell_id, levels = all_cell_ids))
	
	# Now, create the sparse matrix
	raw_matrix <- sparseMatrix(
	  i = count_data$i,
	  j = count_data$j,
	  x = count_data$n,  # This is the count
	  dims = c(length(unique(parquet_file_assign$feature_name)), length(all_cell_ids)),
	  dimnames = list(genes = unique(parquet_file_assign$feature_name), cells = all_cell_ids)
	)
	
	seurat_obj_unfiltered <- CreateSeuratObject(
	  counts = raw_matrix,
	  assay = "Xenium",
	  project = "Xenium_Unfiltered"
	)
	
	#add coordinates
	# Read the cell metadata which contains centroid coordinates
	cor_path <- file.path(data_path, "cells.parquet")
	cells_metadata <- arrow::read_parquet(cor_path)
	
	# Ensure the order of cells in the metadata matches the order in the Seurat object
	# This is critical for adding coordinates correctly
	metadata_ordered <- cells_metadata[match(colnames(seurat_obj_unfiltered), cells_metadata$cell_id), ]
	
	# Add spatial coordinates as a dimensional reduction object
	seurat_obj_unfiltered$x_centroid <- cells_metadata$x_centroid
	seurat_obj_unfiltered$y_centroid <- cells_metadata$y_centroid
	
	saveRDS(seurat_obj_unfiltered, paste0("~/", sample_id, "_xen_so_unQVfiltered.RDS"))
}
