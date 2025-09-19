# Fibroblast Subpopulation Analysis


# Load required libraries
library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)
library(DoubletFinder)
library(clustree)
library(ggplot2)
library(ggsci)

# Clean environment and set working directory
rm(list = ls())
gc()
setwd("D:/data/scRNA/")

# Load annotated Seurat object
scRNA_data <- readRDS("manually_annotated_seurat.rds")

# Extract fibroblast population
fibroblast_data <- subset(scRNA_data, subset = cell_type == "Fibroblast")

# Save fibroblast subset
saveRDS(fibroblast_data, "fibroblast_raw.rds")

# Doublet detection and removal
# ------------------------------
fibroblast_clean <- fibroblast_data

# DoubletFinder function
find_doublets <- function(seurat_obj) {
  # Find optimal pK value
  sweep.res.list <- paramSweep_v3(seurat_obj, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  optimal_pk <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC == max(bcmvn$MeanBC), ]$pK))
  
  # Calculate expected doublet rate
  doublet_rate <- ncol(seurat_obj) * 8 * 1e-6
  nExp_poi <- round(doublet_rate * ncol(seurat_obj))
  
  # Adjust for homotypic doublets
  homotypic_prop <- modelHomotypic(seurat_obj@meta.data$seurat_clusters)
  nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))
  
  # Run DoubletFinder
  seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:20, pN = 0.25, pK = optimal_pk, 
                                 nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(seurat_obj@meta.data)[ncol(seurat_obj@meta.data)] <- "doublet_low"
  
  seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:20, pN = 0.25, pK = optimal_pk, 
                                 nExp = nExp_poi_adj, reuse.pANN = FALSE, sct = FALSE)
  colnames(seurat_obj@meta.data)[ncol(seurat_obj@meta.data)] <- "doublet_high"
  
  return(seurat_obj)
}

# Apply doublet detection
fibroblast_clean <- find_doublets(fibroblast_clean)

# Classify doublet confidence
fibroblast_clean@meta.data$doublet_class <- fibroblast_clean@meta.data$doublet_low
fibroblast_clean@meta.data$doublet_class[
  which(fibroblast_clean@meta.data$doublet_class == "Doublet" & 
          fibroblast_clean@meta.data$doublet_high == "Singlet")
] <- "Doublet-Low"
fibroblast_clean@meta.data$doublet_class[
  which(fibroblast_clean@meta.data$doublet_class == "Doublet")
] <- "Doublet-High"

# Visualize doublets
pdf("fibroblast_doublet_detection.pdf", width = 10, height = 8)
DimPlot(fibroblast_clean, reduction = "umap", group.by = "doublet_class",
        cols = c("red", "gold", "#1bb3b6"), pt.size = 1.0) +
  ggtitle("DoubletFinder Results - Fibroblasts")
dev.off()

# Remove doublets and save clean data
fibroblast_clean <- subset(fibroblast_clean, subset = doublet_class == "Singlet")
saveRDS(fibroblast_clean, "fibroblast_clean.rds")

# Re-clustering of fibroblast subpopulations
# ------------------------------------------
fibroblast_clean <- NormalizeData(fibroblast_clean, verbose = FALSE)
fibroblast_clean <- FindVariableFeatures(fibroblast_clean, selection.method = "vst", 
                                         nfeatures = 2000, verbose = FALSE)
fibroblast_clean <- ScaleData(fibroblast_clean, verbose = FALSE)
fibroblast_clean <- RunPCA(fibroblast_clean, verbose = FALSE)

# Determine optimal number of PCs
pdf("fibroblast_elbow_plot.pdf", width = 8, height = 6)
ElbowPlot(fibroblast_clean, ndims = 30)
dev.off()

# Use 8 PCs based on elbow plot
fibroblast_clean <- FindNeighbors(fibroblast_clean, reduction = "pca", dims = 1:8)

# Test multiple resolutions
resolutions <- c(0.05, 0.08, 0.1, 0.12, 0.14, 0.2, 0.3)
for (res in resolutions) {
  fibroblast_clean <- FindClusters(fibroblast_clean, resolution = res, algorithm = 1)
}

# Visualize resolution selection
pdf("fibroblast_clustree.pdf", width = 12, height = 10)
clustree(fibroblast_clean@meta.data, prefix = "RNA_snn_res.")
dev.off()

# Select optimal resolution (0.08)
fibroblast_clean <- FindClusters(fibroblast_clean, resolution = 0.08)
Idents(fibroblast_clean) <- "RNA_snn_res.0.08"

# Dimensionality reduction
fibroblast_clean <- RunUMAP(fibroblast_clean, dims = 1:8)
fibroblast_clean <- RunTSNE(fibroblast_clean, dims = 1:8)

# Visualize clustering
pdf("fibroblast_subclustering.pdf", width = 16, height = 8)
p1 <- DimPlot(fibroblast_clean, reduction = "umap", label = TRUE, pt.size = 1.5) +
  ggtitle("UMAP - Fibroblast Subclusters")
p2 <- DimPlot(fibroblast_clean, reduction = "tsne", label = TRUE, pt.size = 1.5) +
  ggtitle("t-SNE - Fibroblast Subclusters")
p1 + p2
dev.off()

# Find cluster markers
fibroblast_markers <- FindAllMarkers(fibroblast_clean, 
                                     only.pos = TRUE, 
                                     min.pct = 0.25, 
                                     logfc.threshold = 0.25)

write.csv(fibroblast_markers, "fibroblast_cluster_markers.csv", row.names = FALSE)

# Top markers for each cluster
top_markers <- fibroblast_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

write.csv(top_markers, "fibroblast_top_markers.csv", row.names = FALSE)

# Sample grouping for comparative analysis
fibroblast_clean@meta.data$sample_group <- case_when(
  fibroblast_clean$orig.ident %in% c('GSM5573484_sample19', 'GSM5573485_sample20', 'GSM5573503_sample38') ~ 'Peritoneum_Diffuse',
  fibroblast_clean$orig.ident == 'SeuratProject' ~ 'Peritoneum_Intestinal',
  TRUE ~ 'Primary_gastric'
)

# Find sample-specific markers
Idents(fibroblast_clean) <- "sample_group"
sample_markers <- FindAllMarkers(fibroblast_clean,
                                 only.pos = TRUE,
                                 min.pct = 0.25,
                                 logfc.threshold = 0.25)

write.csv(sample_markers, "fibroblast_sample_markers.csv", row.names = FALSE)

# Composition analysis
# --------------------
composition_data <- data.frame(
  sample = fibroblast_clean$orig.ident,
  cluster = Idents(fibroblast_clean),
  sample_group = fibroblast_clean$sample_group
)

# Calculate cell counts
cell_counts <- composition_data %>%
  group_by(sample_group, cluster) %>%
  summarise(count = n(), .groups = 'drop')

# Stacked bar plots
pdf("fibroblast_composition_plots.pdf", width = 14, height = 10)

# Absolute counts
p1 <- ggplot(cell_counts, aes(x = sample_group, y = count, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_nejm() +
  theme_classic() +
  labs(x = "Sample Group", y = "Cell Count", title = "Fibroblast Subpopulation Composition - Absolute Counts")

# Percentage
p2 <- ggplot(cell_counts, aes(x = sample_group, y = count, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_nejm() +
  theme_classic() +
  labs(x = "Sample Group", y = "Percentage", title = "Fibroblast Subpopulation Composition - Percentage") +
  scale_y_continuous(labels = scales::percent)

# Horizontal view
p3 <- ggplot(cell_counts, aes(x = count, y = sample_group, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_nejm() +
  theme_classic() +
  labs(x = "Percentage", y = "Sample Group", title = "Fibroblast Composition - Horizontal View") +
  scale_x_continuous(labels = scales::percent)

(p1 / p2) | p3

dev.off()

# Save final fibroblast object
saveRDS(fibroblast_clean, "fibroblast_analyzed.rds")

# Generate analysis summary
cat("Fibroblast Subpopulation Analysis Summary\n")
cat("========================================\n\n")
cat("Initial fibroblast cells:", ncol(fibroblast_data), "\n")
cat("After doublet removal:", ncol(fibroblast_clean), "\n")
cat("Number of subclusters:", length(unique(fibroblast_clean$RNA_snn_res.0.08)), "\n")
cat("Sample groups:", paste(unique(fibroblast_clean$sample_group), collapse = ", "), "\n\n")

cat("Subcluster Distribution:\n")
print(table(fibroblast_clean$RNA_snn_res.0.08))

cat("\nSample Group Distribution:\n")
print(table(fibroblast_clean$sample_group))

# Save session information
writeLines(capture.output(sessionInfo()), "fibroblast_analysis_session_info.txt")

cat("\nFibroblast subpopulation analysis completed successfully!\n")