# Dimensionality Reduction, Clustering, and Cell Type Annotation


# Load required libraries
library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)
library(celldex)
library(SingleR)
library(clustree)
library(ggplot2)
library(ggsci)

# Set random seed for reproducibility
set.seed(123)

# Load processed Seurat object
seurat_obj <- readRDS("processed_scRNA_data.rds")

# Principal Component Analysis
# ----------------------------

# Visualize top 20 PCs
pdf("PC_heatmap_top20.pdf", width = 7.5, height = 9)
DimHeatmap(seurat_obj, dims = 1:20, cells = 1000, balanced = TRUE)
dev.off()

# JackStraw and Elbow plots for PC selection
seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20)

pdf("PC_selection_plots.pdf", width = 12, height = 6)
p1 <- JackStrawPlot(seurat_obj, dims = 1:20)
p2 <- ElbowPlot(seurat_obj, ndims = 30)
p1 + p2
dev.off()

# Determine optimal number of PCs
# Based on visual inspection of JackStraw and Elbow plots
optimal_pcs <- 18

# Clustering Analysis
# ------------------

# Find neighbors and perform clustering at multiple resolutions
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:optimal_pcs, reduction = "harmony")

# Test multiple resolutions
resolution_range <- c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1, 1.2)
for (res in resolution_range) {
  seurat_obj <- FindClusters(seurat_obj, resolution = res, algorithm = 1)
}

# Visualize clustering tree
pdf("clustree_resolution_selection.pdf", width = 12, height = 10)
clustree_plot <- clustree(seurat_obj@meta.data, prefix = "RNA_snn_res.")
print(clustree_plot)
dev.off()

# Select optimal resolution based on clustree analysis
optimal_resolution <- 0.4
seurat_obj <- FindClusters(seurat_obj, resolution = optimal_resolution)
Idents(seurat_obj) <- "seurat_clusters"

# Find cluster markers
cluster_markers <- FindAllMarkers(seurat_obj,
                                  only.pos = TRUE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25)

write.csv(cluster_markers, "cluster_markers.csv", row.names = FALSE)

# Visualize top markers
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# Create color palette
color_palette <- c(pal_npg()(9), pal_jco()(9), pal_jama()(7), pal_nejm()(8))

pdf("cluster_markers_heatmap.pdf", width = 22, height = 16)
DoHeatmap(seurat_obj,
          features = top_markers$gene,
          group.colors = color_palette) +
  scale_fill_gradient2(low = '#0099CC', mid = 'white', high = '#CC0033',
                       name = 'Z-score')
dev.off()

# Dimensionality Reduction Visualization
# --------------------------------------

seurat_obj <- RunUMAP(seurat_obj, dims = 1:optimal_pcs, reduction = "harmony")
seurat_obj <- RunTSNE(seurat_obj, dims = 1:optimal_pcs, reduction = "harmony")

# UMAP plots
pdf("clustering_umap_plots.pdf", width = 15, height = 12)
p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE,
              label.size = 3.5, pt.size = 0.5) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5),
        legend.position = "right")

p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident") +
  theme_classic()

p3 <- DimPlot(seurat_obj, reduction = "umap", split.by = "orig.ident") +
  theme_classic()

p1 / (p2 | p3)
dev.off()

# t-SNE plots
pdf("clustering_tsne_plots.pdf", width = 15, height = 12)
p1 <- DimPlot(seurat_obj, reduction = "tsne", label = TRUE,
              label.size = 3.5, pt.size = 0.5) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5),
        legend.position = "right")

p2 <- DimPlot(seurat_obj, reduction = "tsne", group.by = "orig.ident") +
  theme_classic()

p1 | p2
dev.off()

# Save clustered data
saveRDS(seurat_obj, "clustered_data.rds")
