# Manual Cell Type Annotation Pipeline


# Load required libraries
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggsci)
library(clustree)
library(SCpubr)

# Set random seed for reproducibility
set.seed(123)

# Load clustered Seurat object
seurat_obj <- readRDS("clustered_data.rds")

# First-round Manual Annotation
# -----------------------------

# Find neighbors and test multiple resolutions for initial broad annotation
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:18, reduction = "harmony")

# Test multiple resolutions for optimal broad clustering
resolutions <- c(0.01, 0.05, 0.1, 1, 1.5, 2, 2.5, 3, 3.5, 4)
for (res in resolutions) {
  seurat_obj <- FindClusters(seurat_obj, 
                             graph.name = "RNA_snn", 
                             resolution = res, 
                             algorithm = 1)
}

# Visualize resolution selection
pdf("resolution_selection_clustree.pdf", width = 12, height = 10)
clustree_plot <- clustree(seurat_obj@meta.data, prefix = "RNA_snn_res.")
print(clustree_plot)
dev.off()

# Select optimal resolution for initial annotation
seurat_obj <- FindClusters(seurat_obj, resolution = 0.4)
Idents(seurat_obj) <- "seurat_clusters"

# Define broad cell type marker genes
cell_type_markers <- list(
  immune = c("PTPRC"),           # Pan-immune marker
  epithelial = c("EPCAM"),       # Epithelial cells
  stromal = c("MME", "PECAM1")   # Stromal cells (mesenchymal and endothelial)
)

# Visualize broad cell type markers
pdf("broad_cell_type_dotplot.pdf", width = 20, height = 15)
dot_plot <- SCpubr::do_DotPlot(
  sample = seurat_obj,
  features = cell_type_markers,
  dot.scale = 12,
  colors.use = c("white", "yellow", "red"),
  legend.framewidth = 2,
  font.size = 10
)
print(dot_plot)
dev.off()

# Detailed Marker Gene Visualization
# ----------------------------------

# Fibroblast markers
fibroblast_markers <- c("DCN", "COL1A1", "COL1A2")
pdf("fibroblast_markers.pdf", width = 12, height = 8)
FeaturePlot(seurat_obj,
            reduction = "tsne",
            features = fibroblast_markers,
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            raster = FALSE)
dev.off()

# Endothelial markers
endothelial_markers <- c("ENG", "VWF")
pdf("endothelial_markers.pdf", width = 12, height = 8)
FeaturePlot(seurat_obj,
            reduction = "tsne",
            features = endothelial_markers,
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            raster = FALSE)
dev.off()

# Epithelial markers
epithelial_markers <- c("EPCAM", "KRT8", "TFF1")
pdf("epithelial_markers.pdf", width = 12, height = 8)
FeaturePlot(seurat_obj,
            reduction = "tsne",
            features = epithelial_markers,
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            raster = FALSE)
dev.off()

# B cell markers
b_cell_markers <- c("CD79A", "MS4A1")
pdf("b_cell_markers.pdf", width = 12, height = 8)
FeaturePlot(seurat_obj,
            reduction = "tsne",
            features = b_cell_markers,
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            raster = FALSE)
dev.off()

# T cell markers
t_cell_markers <- c("CD3D", "CD3E", "CD2")
pdf("t_cell_markers.pdf", width = 12, height = 8)
FeaturePlot(seurat_obj,
            reduction = "tsne",
            features = t_cell_markers,
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            raster = FALSE)
dev.off()

# Myeloid cell markers
myeloid_markers <- c("CD14", "CD163", "CD68", "AIF1")
pdf("myeloid_markers.pdf", width = 12, height = 8)
FeaturePlot(seurat_obj,
            reduction = "tsne",
            features = myeloid_markers,
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            raster = FALSE)
dev.off()

# Endocrine cell markers
endocrine_markers <- c("GHRL", "GAST")
pdf("endocrine_markers.pdf", width = 12, height = 8)
FeaturePlot(seurat_obj,
            reduction = "tsne",
            features = endocrine_markers,
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            raster = FALSE)
dev.off()

# Mast cell markers
mast_cell_markers <- c("TPSAB1", "CPA3")
pdf("mast_cell_markers.pdf", width = 12, height = 8)
FeaturePlot(seurat_obj,
            reduction = "tsne",
            features = mast_cell_markers,
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            raster = FALSE)
dev.off()

# Manual Cell Type Assignment
# ---------------------------

Idents(seurat_obj) <- "seurat_clusters"

# Assign cell types based on marker expression patterns
seurat_obj <- RenameIdents(seurat_obj,
                           "0" = "T_NK_cell",
                           "1" = "T_NK_cell", 
                           "2" = "Epithelial_cell", 
                           "3" = "B_cell", 
                           "4" = "Endothelial_cell", 
                           "5" = "Myeloid_cell",
                           "6" = "Fibroblast", 
                           "7" = "Myeloid_cell", 
                           "8" = "Epithelial_cell",
                           "9" = "Fibroblast", 
                           "10" = "Myeloid_cell", 
                           "11" = "Mast_cell", 
                           "12" = "Endocrine_cell",
                           "13" = "B_cell", 
                           "14" = "Endothelial_cell", 
                           "15" = "B_cell",
                           "16" = "Fibroblast"
)

# Add cell type information to metadata
seurat_obj@meta.data$cell_type <- Idents(seurat_obj)

# Visualize annotated clusters
pdf("manual_annotation_results.pdf", width = 15, height = 12)
p1 <- DimPlot(seurat_obj, 
              reduction = "tsne", 
              group.by = "cell_type",
              label = TRUE, 
              pt.size = 0.5) +
  ggtitle("t-SNE - Manual Cell Type Annotation") +
  theme_classic()

p2 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "cell_type",
              label = TRUE, 
              pt.size = 0.5) +
  ggtitle("UMAP - Manual Cell Type Annotation") +
  theme_classic()

p1 | p2
dev.off()

# Save cell type distribution statistics
cell_type_distribution <- table(seurat_obj@meta.data$cell_type)
write.csv(as.data.frame(cell_type_distribution), 
          "cell_type_distribution.csv")

# Generate marker expression violin plots for validation
validation_markers <- c("PTPRC", "EPCAM", "MME", "CD3D", "CD79A", "VWF", "COL1A1")
pdf("validation_markers_violin.pdf", width = 14, height = 10)
VlnPlot(seurat_obj, features = validation_markers, group.by = "cell_type", 
        pt.size = 0, ncol = 4)
dev.off()

# Save annotated object
saveRDS(seurat_obj, "manually_annotated_seurat.rds")

# Generate annotation summary
cat("Manual Annotation Summary:\n")
cat("=========================\n")
cat("Total clusters annotated:", length(unique(seurat_obj$seurat_clusters)), "\n")
cat("Cell types identified:", length(unique(seurat_obj$cell_type)), "\n")
cat("Total cells:", ncol(seurat_obj), "\n\n")

cat("Cell Type Distribution:\n")
print(cell_type_distribution)

# Save session information
writeLines(capture.output(sessionInfo()), "annotation_session_info.txt")

cat("Manual annotation completed successfully!\n")