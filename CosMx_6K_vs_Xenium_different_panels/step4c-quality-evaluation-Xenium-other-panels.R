library(data.table)
library(ggrepel)
library(ggpubr)
library(Seurat)
library(Matrix)
library(ggsci)
library(stringr)

sample_list <- fread("CosMx-vs-Xenium-4.txt")

results_uncorrected_full <- list()
j <- 1

results_uncorrected_overlap <- list()
k <- 1

file_paths <- c("_xen_377_with_seg_so.RDS", "_xen_377_no_seg_so.RDS", "_xen_colon_with_seg_so.RDS", "_xen_377_with_seg_so_qcnorm.RDS", "_xen_377_no_seg_so_qcnorm.RDS", "_xen_colon_with_seg_so_qcnorm.RDS")
file_descs <- c("Xenium 377, with seg, no QC", "Xenium 377, no seg, no QC", "Xenium colon-specific, with seg, no QC", "Xenium 377, with seg, with QC", "Xenium 377, no seg, with QC", "Xenium colon-specific, with seg, with QC")

i <- 1
sample_id <- sample_list$sample_id[i]
xen_5K <- readRDS(paste0("comparison/", sample_id, "_xen_so.RDS"))
xen_377 <- readRDS(paste0("comparison/", sample_id, "_xen_377_with_seg_so.RDS"))
xen_cs <- readRDS(paste0("comparison/", sample_id, "_xen_colon_with_seg_so.RDS"))
common_genes <- intersect(rownames(xen_5K), rownames(xen_377))
common_genes <- intersect(common_genes, rownames(xen_cs)) # 48

for (i in 1:nrow(sample_list)) {
	# i <- 1
	sample_id <- sample_list$sample_id[i]
	
	xen_5K <- readRDS(paste0("comparison/", sample_id, "_xen_so.RDS"))
	xen_5K_counts <- xen_5K@assays$Xenium$counts
	xen_5K_neg <- GetAssayData(xen_5K, "ControlProbe")
	xen_5K_counts <- xen_5K@assays$Xenium$counts
	total_count <- colSums(xen_5K_counts)
	xen_5K_counts_pos <- xen_5K_counts > 0
	diversity <- colSums(xen_5K_counts_pos)
	total_neg_count <- colSums(xen_5K_neg)
	xen_5K_neg_pos <- xen_5K_neg > 0
	n_neg_detected <- colSums(xen_5K_neg_pos)
	
	results_uncorrected_full[[j]] <- data.table(sample_id = rep(sample_id, ncol(xen_5K_counts)), tech = rep("Xenium 5K, with seg, no QC", ncol(xen_5K_counts)), cell_id = colnames(xen_5K_counts), n_cell = rep(ncol(xen_5K_counts), ncol(xen_5K_counts)), n_gene = rep(nrow(xen_5K_counts), ncol(xen_5K_counts)), total_count = as.integer(total_count), diversity = as.integer(diversity), total_neg_count = as.integer(total_neg_count), n_neg_detected = as.integer(n_neg_detected), cell_area = xen_5K$cell_area)
	j <- j + 1
	rm(total_count, diversity, total_neg_count, n_neg_detected)
	
	xen_5K <- readRDS(paste0("comparison/", sample_id, "_xen_so_qcnorm.RDS"))
	xen_5K_counts_qc <- xen_5K@assays$Xenium$counts
	xen_5K_neg <- GetAssayData(xen_5K, "ControlProbe")
	xen_5K_counts_qc <- xen_5K@assays$Xenium$counts
	total_count <- colSums(xen_5K_counts_qc)
	xen_5K_counts_qc_pos <- xen_5K_counts_qc > 0
	diversity <- colSums(xen_5K_counts_qc_pos)
	total_neg_count <- colSums(xen_5K_neg)
	xen_5K_neg_pos <- xen_5K_neg > 0
	n_neg_detected <- colSums(xen_5K_neg_pos)
	
	results_uncorrected_full[[j]] <- data.table(sample_id = rep(sample_id, ncol(xen_5K_counts_qc)), tech = rep("Xenium 5K, with seg, with QC", ncol(xen_5K_counts_qc)), cell_id = colnames(xen_5K_counts_qc), n_cell = rep(ncol(xen_5K_counts_qc), ncol(xen_5K_counts_qc)), n_gene = rep(nrow(xen_5K_counts_qc), ncol(xen_5K_counts_qc)), total_count = as.integer(total_count), diversity = as.integer(diversity), total_neg_count = as.integer(total_neg_count), n_neg_detected = as.integer(n_neg_detected), cell_area = xen_5K$cell_area)
	j <- j + 1
	rm(total_count, diversity, total_neg_count, n_neg_detected)
	
	for (m in 1:6) {
		file_path <- file_paths[m]
		file_desc <- file_descs[m]
		
		xen <- readRDS(paste0("comparison/", sample_id, file_path))
		### Retrieve negprobes and counts (full)
		xen_neg <- GetAssayData(xen, "ControlProbe")
		xen_counts <- xen@assays$Xenium$counts
		total_count <- colSums(xen_counts)
		xen_counts_pos <- xen_counts > 0
		diversity <- colSums(xen_counts_pos)
		total_neg_count <- colSums(xen_neg)
		xen_neg_pos <- xen_neg > 0
		n_neg_detected <- colSums(xen_neg_pos)
		
		results_uncorrected_full[[j]] <- data.table(sample_id = rep(sample_id, ncol(xen_counts)), tech = rep(file_desc, ncol(xen_counts)), cell_id = colnames(xen_counts), n_cell = rep(ncol(xen_counts), ncol(xen_counts)), n_gene = rep(nrow(xen_counts), ncol(xen_counts)), total_count = as.integer(total_count), diversity = as.integer(diversity), total_neg_count = as.integer(total_neg_count), n_neg_detected = as.integer(n_neg_detected), cell_area = xen$cell_area)
		j <- j + 1
		rm(total_count, diversity, total_neg_count, n_neg_detected)
	
		### Retrieve counts (overlapped)
		xen_counts_overlap <- xen_counts[common_genes,]
		total_count <- colSums(xen_counts_overlap)
		xen_counts_overlap_pos <- xen_counts_overlap > 0
		diversity <- colSums(xen_counts_overlap_pos)
		results_uncorrected_overlap[[k]] <- data.table(sample_id = rep(sample_id, ncol(xen_counts_overlap)), tech = rep(file_desc, ncol(xen_counts_overlap)), cell_id = colnames(xen_counts_overlap), n_cell = rep(ncol(xen_counts_overlap), ncol(xen_counts_overlap)), n_gene = rep(nrow(xen_counts_overlap), ncol(xen_counts_overlap)), total_count = as.integer(total_count), diversity = as.integer(diversity), cell_area = xen$cell_area)
		k <- k + 1
		rm(total_count, diversity)
		
		if (m == 1) {
			xen_5K_counts_overlap <- xen_5K_counts[common_genes,]
			total_count <- colSums(xen_5K_counts_overlap)
			xen_5K_counts_overlap_pos <- xen_5K_counts_overlap > 0
			diversity <- colSums(xen_5K_counts_overlap_pos)
			results_uncorrected_overlap[[k]] <- data.table(sample_id = rep(sample_id, ncol(xen_5K_counts_overlap)), tech = rep("Xenium 5K, with seg, no QC", ncol(xen_5K_counts_overlap)), cell_id = colnames(xen_5K_counts_overlap), n_cell = rep(ncol(xen_5K_counts_overlap), ncol(xen_5K_counts_overlap)), n_gene = rep(nrow(xen_5K_counts_overlap), ncol(xen_5K_counts_overlap)), total_count = as.integer(total_count), diversity = as.integer(diversity), cell_area = xen_5K$cell_area)
			k <- k + 1
			rm(total_count, diversity)
			
			xen_5K_counts_qc_overlap <- xen_5K_counts_qc[common_genes,]
			total_count <- colSums(xen_5K_counts_qc_overlap)
			xen_5K_counts_qc_overlap_pos <- xen_5K_counts_qc_overlap > 0
			diversity <- colSums(xen_5K_counts_qc_overlap_pos)
			results_uncorrected_overlap[[k]] <- data.table(sample_id = rep(sample_id, ncol(xen_5K_counts_qc_overlap)), tech = rep("Xenium 5K, with seg, with QC", ncol(xen_5K_counts_qc_overlap)), cell_id = colnames(xen_5K_counts_qc_overlap), n_cell = rep(ncol(xen_5K_counts_qc_overlap), ncol(xen_5K_counts_qc_overlap)), n_gene = rep(nrow(xen_5K_counts_qc_overlap), ncol(xen_5K_counts_qc_overlap)), total_count = as.integer(total_count), diversity = as.integer(diversity), cell_area = xen_5K$cell_area)
			k <- k + 1
			rm(total_count, diversity)
		}
		gc()
	}
}

results_uncorrected_full <- rbindlist(results_uncorrected_full)
results_uncorrected_overlap <- rbindlist(results_uncorrected_overlap)

results_uncorrected_full$total_count_by_area <- results_uncorrected_full$total_count / results_uncorrected_full$cell_area
results_uncorrected_full$diversity_by_area <- results_uncorrected_full$diversity / results_uncorrected_full$cell_area
results_uncorrected_full$total_neg_count_by_area <- results_uncorrected_full$total_neg_count / results_uncorrected_full$cell_area
results_uncorrected_full$n_neg_detected_by_area <- results_uncorrected_full$n_neg_detected / results_uncorrected_full$cell_area

results_uncorrected_overlap$total_count_by_area <- results_uncorrected_overlap$total_count / results_uncorrected_overlap$cell_area
results_uncorrected_overlap$diversity_by_area <- results_uncorrected_overlap$diversity / results_uncorrected_overlap$cell_area

fwrite(results_uncorrected_full, "comparison/xen_results_combined_full_panel.txt", sep = "\t")
fwrite(results_uncorrected_overlap, "comparison/xen_results_combined_overlap_panel.txt", sep = "\t")