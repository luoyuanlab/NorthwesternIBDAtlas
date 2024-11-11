library(data.table)
library(ggrepel)
library(ggpubr)
library(Seurat)
library(Matrix)
library(ggsci)
library(stringr)

sample_list <- fread("CosMx-vs-Xenium-4.txt")

##### Get transcript diversity, transcript total counts, and total negprobes in the full/overlapping panel

results_uncorrected_full <- list()
j <- 1
results_uncorrected_overlap <- list()
k <- 1

for (i in 1:nrow(sample_list)) {
	# i <- 1
	sample_id <- sample_list$sample_id[i]
	
	smi <- readRDS(paste0("comparison/", as.character(sample_id), "_smi_so_qcnorm.RDS"))
	xen <- readRDS(paste0("comparison/", as.character(sample_id), "_xen_so_qcnorm.RDS"))
	panel <- readRDS(paste0("comparison/", as.character(sample_id), "_panel_full.RDS"))
	
	### Retrieve negprobes
	smi_neg <- GetAssayData(smi, "negprobes")
	xen_neg <- GetAssayData(xen, "ControlProbe")
	
	### Retrieve counts (full panel)
	smi_counts <- smi@assays$Nanostring@layers$counts
	rownames(smi_counts) <- rownames(smi@assays$Nanostring@features)
	colnames(smi_counts) <- rownames(smi@assays$Nanostring@cells)
	xen_counts <- xen@assays$Xenium$counts
	
	### CosMx results (full panel)
	total_count <- colSums(smi_counts)
	smi_counts_pos <- smi_counts > 0
	diversity <- colSums(smi_counts_pos)
	total_neg_count <- colSums(smi_neg)
	smi_neg_pos <- smi_neg > 0
	n_neg_detected <- colSums(smi_neg_pos)
	results_uncorrected_full[[j]] <- data.table(sample_id = rep(sample_id, ncol(smi_counts)), tech = rep("CosMx", ncol(smi_counts)), cell_id = colnames(smi_counts), n_cell = rep(ncol(smi_counts), ncol(smi_counts)), n_gene = rep(nrow(smi_counts), ncol(smi_counts)), total_count = as.integer(total_count), diversity = as.integer(diversity), total_neg_count = as.integer(total_neg_count), n_neg_detected = as.integer(n_neg_detected))
	j <- j + 1
	rm(total_count, diversity, total_neg_count, n_neg_detected)
	
	### Xenium results (full panel)
	total_count <- colSums(xen_counts)
	xen_counts_pos <- xen_counts > 0
	diversity <- colSums(xen_counts_pos)
	total_neg_count <- colSums(xen_neg)
	xen_neg_pos <- xen_neg > 0
	n_neg_detected <- colSums(xen_neg_pos)
	results_uncorrected_full[[j]] <- data.table(sample_id = rep(sample_id, ncol(xen_counts)), tech = rep("Xenium", ncol(xen_counts)), cell_id = colnames(xen_counts), n_cell = rep(ncol(xen_counts), ncol(xen_counts)), n_gene = rep(nrow(xen_counts), ncol(xen_counts)), total_count = as.integer(total_count), diversity = as.integer(diversity), total_neg_count = as.integer(total_neg_count), n_neg_detected = as.integer(n_neg_detected))
	j <- j + 1
	rm(total_count, diversity, total_neg_count, n_neg_detected)
	
	### Retrieve counts (overlapping panel)
	common_genes <- intersect(rownames(smi_counts), rownames(xen_counts))
	smi_counts_overlap <- smi_counts[common_genes,]
	xen_counts_overlap <- xen_counts[common_genes,]

	### CosMx results (overlapping panel)
	total_count <- colSums(smi_counts_overlap)
	smi_counts_overlap_pos <- smi_counts_overlap > 0
	diversity <- colSums(smi_counts_overlap_pos)
	results_uncorrected_overlap[[k]] <- data.table(sample_id = rep(sample_id, ncol(smi_counts_overlap)), tech = rep("CosMx", ncol(smi_counts_overlap)), cell_id = colnames(smi_counts_overlap), n_cell = rep(ncol(smi_counts_overlap), ncol(smi_counts_overlap)), n_gene = rep(nrow(smi_counts_overlap), ncol(smi_counts_overlap)), total_count = as.integer(total_count), diversity = as.integer(diversity))
	k <- k + 1

	### Xenium results (overlapping panel)
	total_count <- colSums(xen_counts_overlap)
	xen_counts_overlap_pos <- xen_counts_overlap > 0
	diversity <- colSums(xen_counts_overlap_pos)
	results_uncorrected_overlap[[k]] <- data.table(sample_id = rep(sample_id, ncol(xen_counts_overlap)), tech = rep("Xenium", ncol(xen_counts_overlap)), cell_id = colnames(xen_counts_overlap), n_cell = rep(ncol(xen_counts_overlap), ncol(xen_counts_overlap)), n_gene = rep(nrow(xen_counts_overlap), ncol(xen_counts_overlap)), total_count = as.integer(total_count), diversity = as.integer(diversity))
	k <- k + 1
	
	gc()
}

results_uncorrected_full <- rbindlist(results_uncorrected_full)
results_uncorrected_overlap <- rbindlist(results_uncorrected_overlap)

fwrite(results_uncorrected_full, "comparison/results_uncorrected_full.txt", sep = "\t")
fwrite(results_uncorrected_overlap, "comparison/results_uncorrected_overlap.txt", sep = "\t")

##### Pull cell areas

results_areas <- list()

for (i in 1:nrow(sample_list)) {
	# i <- 1
	sample_id <- sample_list$sample_id[i]
	
	smi <- readRDS(paste0("comparison/", as.character(sample_id), "_smi_so_qcnorm.RDS"))
	xen <- readRDS(paste0("comparison/", as.character(sample_id), "_xen_so_qcnorm.RDS"))
	
	smi$cell_area_um2 <- smi$Area * 0.12 * 0.12
	results_areas[[i]] <- data.table(sample_id = rep(sample_id, ncol(smi) + ncol(xen)), tech = c(rep("CosMx", ncol(smi)), rep("Xenium", ncol(xen))), cell_id = c(colnames(smi), colnames(xen)), cell_area = c(smi$cell_area_um2, xen$cell_area), cell_type_minor = c(smi$ct_minor, as.character(xen$ct_minor)), cell_type_major = c(smi$ct_major, xen$ct_major))
}

results_areas <- rbindlist(results_areas)
fwrite(results_areas, "comparison/results_areas.txt", sep = "\t")

##### Standardize quality metrics by cell area

library(data.table)

results_uncorrected_full <- fread("comparison/results_uncorrected_full.txt")
results_uncorrected_overlap <- fread("comparison/results_uncorrected_overlap.txt")
results_areas <- fread("comparison/results_areas.txt")

results_uncorrected_full <- merge(results_uncorrected_full, results_areas, by = c("sample_id", "tech", "cell_id"))
results_uncorrected_overlap <- merge(results_uncorrected_overlap, results_areas, by = c("sample_id", "tech", "cell_id"))

results_uncorrected_full$diversity_by_area <- results_uncorrected_full$diversity / results_uncorrected_full$cell_area
results_uncorrected_full$total_count_by_area <- results_uncorrected_full$total_count / results_uncorrected_full$cell_area
results_uncorrected_full$total_neg_count_by_area <- results_uncorrected_full$total_neg_count / results_uncorrected_full$cell_area
results_uncorrected_full$n_neg_detected_by_area <- results_uncorrected_full$n_neg_detected / results_uncorrected_full$cell_area

results_uncorrected_overlap$diversity_by_area <- results_uncorrected_overlap$diversity / results_uncorrected_overlap$cell_area
results_uncorrected_overlap$total_count_by_area <- results_uncorrected_overlap$total_count / results_uncorrected_overlap$cell_area

fwrite(results_uncorrected_full, "comparison/results_combined_full_panel.txt", sep = "\t")
fwrite(results_uncorrected_overlap, "comparison/results_combined_overlap_panel.txt", sep = "\t")
