library(Matrix)
library(mclust)
library(pheatmap)
library(RColorBrewer)
library(Seurat)
library(ggplot2)
library(ggsci)
library(data.table)

sample_list <- fread("/share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex/CosMx-vs-Xenium-16.txt")
sample_list <- sample_list[sample_list$disease == "CD",]

### Crohn's disease reference -- https://pubmed.ncbi.nlm.nih.gov/34497389/
sc <- readRDS("/share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex/scRNA-seq/CD/sc.RDS")
cell_map <- readRDS("/share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex/scRNA-seq/CD/cell_map.RDS")
ref_profile <- readRDS("/share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex/scRNA-seq/CD/minor_profile.RDS")

for (i in 1:nrow(sample_list)) {
	# i <- 1
	sample_id <- sample_list$sample_id[i]
	
	smi <- readRDS(paste0("/share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex/comparison_Cenfu/", as.character(sample_id), "_smi_so_norm.RDS"))
	xen <- readRDS(paste0("/share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex/comparison_Cenfu/", as.character(sample_id), "_xen_so_norm.RDS"))
	xen_unfiltered <- readRDS(paste0("/share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex/comparison_Cenfu/", as.character(sample_id), "_xen_so_unQVfiltered_norm.RDS"))
	
	
	### Using Seurat to perform cell typing
	# Find anchors exceeds future.globals.maxSize, so raise it before running
	options(future.globals.maxSize = 3000*1024^2)
	
	### CosMx (RNA only)
	# Find anchors 
	anchors <- FindTransferAnchors(reference = sc, query = smi)
	# Pull predictions
	predictions <- TransferData(anchorset = anchors, refdata = sc$celltype_minor, weight.reduction = "pcaproject")
	# Add predictions to seurat
	if (all(colnames(smi) == rownames(predictions))) {
		smi@meta.data$ct_minor <- predictions$predicted.id
	}
	rm(anchors, predictions)
	
	### Xenium
	# Find anchors 
	anchors <- FindTransferAnchors(reference = sc, query =xen)
	# Pull predictions
	predictions <-TransferData(anchorset = anchors, refdata = sc$celltype_minor, weight.reduction = "pcaproject")
	# Add predictions to seurat
	if (all(colnames(xen) == rownames(predictions))) {
		xen@meta.data$ct_minor <- predictions$predicted.id 
	}
	rm(anchors, predictions)
	
	### Xenium un QV filtered
	# Find anchors 
	anchors <- FindTransferAnchors(reference = sc, query =xen_unfiltered)
	# Pull predictions
	predictions <-TransferData(anchorset = anchors, refdata = sc$celltype_minor, weight.reduction = "pcaproject")
	# Add predictions to seurat
	if (all(colnames(xen_unfiltered) == rownames(predictions))) {
	  xen_unfiltered@meta.data$ct_minor <- predictions$predicted.id 
	}
	rm(anchors, predictions)
	
	smi$ct_major <- unname(cell_map[smi$ct_minor])
	xen$ct_major <- unname(cell_map[xen$ct_minor])
	xen_unfiltered$ct_major <- unname(cell_map[xen_unfiltered$ct_minor])
	saveRDS(smi, paste0("/share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex/comparison_Cenfu/", sample_id, "_smi_so_norm.RDS"))
	saveRDS(xen, paste0("/share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex/comparison_Cenfu/", sample_id, "_xen_so_norm.RDS"))
	saveRDS(xen_unfiltered, paste0("/share/fsmresfiles/UC/AtoMx/UC_serial_cuts/1000plex/comparison_Cenfu/", sample_id, "_xen_so_unQVfiltered_norm.RDS"))
	
	gc()
}




