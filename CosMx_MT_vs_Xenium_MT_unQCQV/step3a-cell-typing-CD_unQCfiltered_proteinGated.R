### CosMx data cell type annotation - protein-gated

library(Matrix)
library(mclust)
library(pheatmap)
library(RColorBrewer)
library(Seurat)
library(ggplot2)
library(ggsci)
library(data.table)

mm_scale <- function(x) {
	return((x - min(x)) / (max(x) - min(x)))
}

sample_list <- fread("~/CosMx-vs-Xenium-16.txt")
sample_list <- sample_list[sample_list$disease == "CD",]

### Crohn's disease reference -- https://pubmed.ncbi.nlm.nih.gov/34497389/
### CD

sc_imm <- readRDS("~/scRNA-seq/CD/sc_imm.RDS")
cell_map_imm <- readRDS("~/scRNA-seq/CD/cell_map_imm.RDS")

sc_epi <- readRDS("~/scRNA-seq/CD/sc_epi.RDS")
cell_map_epi <- readRDS("~/scRNA-seq/CD/cell_map_epi.RDS")

sc_others <- readRDS("~/scRNA-seq/CD/sc_others.RDS")
cell_map_others <- readRDS("~/scRNA-seq/CD/cell_map_others.RDS")
gc()

summary <- list()

for (i in 1:nrow(sample_list)) {
	# i <- 1
	sample_id <- sample_list$sample_id[i]
	smi <- readRDS(paste0("~/", as.character(sample_id), "_smi_so_norm.RDS"))

	### Determine cell category (epithelial, immune, others)
	smi$Max.PanCK.scaled <- mm_scale(smi@meta.data$Max.PanCK)
	smi$Max.CD45.scaled <- mm_scale(smi@meta.data$Max.CD45)
	anchor_value <- max(smi@meta.data$Max.PanCK.scaled[which(abs(smi@meta.data$Max.PanCK - 5000) == min(abs(smi@meta.data$Max.PanCK - 5000)))])
	smi$PanCK_pos <- smi@meta.data$Max.PanCK > 5000
	smi$CD45_pos <- smi@meta.data$Max.CD45.scaled > anchor_value
	smi$category <- ifelse(smi@meta.data$PanCK_pos, "Epithelial", ifelse(smi@meta.data$CD45_pos, "Immune", "Others"))
	summary[[i]] <- data.table(sample_id = sample_id, disease = "CD", n_cell = nrow(smi@meta.data), n_cell_epithelial = sum(smi@meta.data$category == "Epithelial"), n_cell_immune = sum(smi@meta.data$category == "Immune"), n_cell_others = sum(smi@meta.data$category == "Others"), scaled_thr = anchor_value)
	
	# Find anchors exceeds future.globals.maxSize, so raise it before running
	options(future.globals.maxSize = 16000*1024^2)
	
	### Annotate epithelial cells
	smi_epi <- smi[,smi@meta.data$category == "Epithelial"]	
	anchors <- FindTransferAnchors(reference = sc_epi, query = smi_epi)
	predictions <- TransferData(anchorset = anchors, refdata = sc_epi$celltype_minor, weight.reduction = "pcaproject")
	smi_epi$ct_minor <- predictions$predicted.id
	smi_epi$ct_major <- unname(cell_map_epi[smi_epi@meta.data$ct_minor])
	meta_epi <- smi_epi@meta.data[,c("cell", "ct_minor", "ct_major")]
	rm(anchors, predictions, smi_epi)
	gc()
	
	### Annotate immune cells
	smi_imm <- smi[,smi@meta.data$category == "Immune"]	
	anchors <- FindTransferAnchors(reference = sc_imm, query = smi_imm)
	predictions <- TransferData(anchorset = anchors, refdata = sc_imm$celltype_minor, weight.reduction = "pcaproject")
	smi_imm$ct_minor <- predictions$predicted.id
	smi_imm$ct_major <- unname(cell_map_imm[smi_imm@meta.data$ct_minor])
	meta_imm <- smi_imm@meta.data[,c("cell", "ct_minor", "ct_major")]
	rm(anchors, predictions, smi_imm)
	gc()
	
	### Annotate other cells
	smi_others <- smi[,smi@meta.data$category == "Others"]	
	anchors <- FindTransferAnchors(reference = sc_others, query = smi_others)
	predictions <- TransferData(anchorset = anchors, refdata = sc_others$celltype_minor, weight.reduction = "pcaproject")
	smi_others$ct_minor <- predictions$predicted.id
	smi_others$ct_major <- unname(cell_map_others[smi_others@meta.data$ct_minor])
	meta_others <- smi_others@meta.data[,c("cell", "ct_minor", "ct_major")]
	rm(anchors, predictions, smi_others)
	gc()
	
	### Merge back the stratified cell type annotations
	meta_new <- rbind(meta_epi, meta_imm, meta_others)
	meta <- smi@meta.data
	meta$index <- 1:nrow(meta)
	meta <- meta[,c("cell", "index")]
	meta <- merge(meta, meta_new, by = "cell")
	meta <- meta[order(meta$index),]
	smi@meta.data$ct_minor_new <- meta$ct_minor
	smi@meta.data$ct_major_new <- meta$ct_major
	
	saveRDS(smi, paste0("~/", sample_id, "_smi_so_norm.RDS"))
	gc()
}

summary <- rbindlist(summary)
summary$p_epithelial <- summary$n_cell_epithelial / summary$n_cell
summary$p_immune <- summary$n_cell_immune / summary$n_cell
summary$p_others <- summary$n_cell_others / summary$n_cell
fwrite(summary, "~/CosMx16_CD_ct_prop_summary.txt", sep = "\t")
