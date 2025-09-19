# Single-cell RNA Sequencing Analysis Pipeline

# Load required libraries
library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)
library(data.table)
library(Matrix)

# Clean environment and set working directory
rm(list = ls())
gc()
memory.limit(size = 99990000000)
setwd("D:/data/scRNA/work/")

# Function to create Seurat objects from different file formats
create_seurat_object <- function(dir_path, file_pattern = NULL, 
                                 sample_names = NULL, data_type = "csv") {
  
  if (data_type == "10X") {
    # For 10X Genomics data
    samples <- list.files(dir_path)
    seurat_list <- lapply(samples, function(sample) {
      cat("Processing 10X sample:", sample, "\n")
      counts <- Read10X(file.path(dir_path, sample))
      sce <- CreateSeuratObject(counts = counts,
                                project = sample,
                                min.cells = 5,
                                min.features = 300)
      return(sce)
    })
    names(seurat_list) <- samples
    
  } else if (data_type == "csv") {
    # For CSV/GZ files
    samples <- list.files(dir_path, pattern = file_pattern)
    seurat_list <- lapply(samples, function(pro) {
      cat("Processing CSV sample:", pro, "\n")
      ct <- fread(file.path(dir_path, pro), data.table = FALSE)
      rownames(ct) <- ct[, 1]
      ct <- ct[, -1]
      
      # Create project name
      project_name <- gsub('\\.csv\\.gz$', '', gsub('^GSM[0-9]*_', '', pro))
      
      sce <- CreateSeuratObject(counts = ct,
                                project = project_name,
                                min.cells = 5,
                                min.features = 300)
      return(sce)
    })
    names(seurat_list) <- gsub('\\.csv\\.gz$', '', gsub('^GSM[0-9]*_', '', samples))
    
  } else if (data_type == "txt") {
    # For TXT/GZ files
    samples <- list.files(dir_path, pattern = file_pattern)
    seurat_list <- lapply(samples, function(pro) {
      cat("Processing TXT sample:", pro, "\n")
      ct <- fread(file.path(dir_path, pro), data.table = FALSE)
      rownames(ct) <- ct[, 1]
      ct <- ct[, -1]
      
      project_name <- gsub('_scgex\\.txt\\.gz$', '', gsub('^GSM[0-9]*_', '', pro))
      
      sce <- CreateSeuratObject(counts = ct,
                                project = project_name,
                                min.cells = 5,
                                min.features = 300)
      return(sce)
    })
    names(seurat_list) <- gsub('_scgex\\.txt\\.gz$', '', gsub('^GSM[0-9]*_', '', samples))
  }
  
  return(seurat_list)
}

# Function to merge Seurat objects with proper naming
merge_seurat_objects <- function(seurat_list, output_name) {
  if (length(seurat_list) == 0) {
    stop("Empty Seurat object list")
  }
  
  # Ensure all objects have proper names
  if (is.null(names(seurat_list))) {
    names(seurat_list) <- paste0("Sample_", seq_along(seurat_list))
  }
  
  # Merge objects
  merged_seurat <- seurat_list[[1]]
  if (length(seurat_list) > 1) {
    for (i in 2:length(seurat_list)) {
      merged_seurat <- merge(merged_seurat, y = seurat_list[[i]])
    }
  }
  
  # Save merged object
  saveRDS(merged_seurat, file = paste0(output_name, ".rds"))
  
  # Print summary
  cat("Merged object dimensions:", dim(merged_seurat), "\n")
  cat("Cell counts per sample:\n")
  print(table(merged_seurat@meta.data$orig.ident))
  
  return(merged_seurat)
}

# Process Primary Tumor Datasets
# ------------------------------

# GSE183904 - CSV format
gse183904_dir <- "D:/data/scRNA/datas/D/GSE183904"
gse183904_data <- create_seurat_object(gse183904_dir, 
                                       file_pattern = "\\.csv\\.gz$", 
                                       data_type = "csv")

# GSE163558 - 10X format
gse163558_dir <- "D:/data/scRNA/datas/D/GSE163558/GT"
gse163558_data <- create_seurat_object(gse163558_dir, 
                                       data_type = "10X")

# GSE206785 - TXT format
gse206785_dir <- "D:/data/scRNA/datas/D/GSE206785"
gse206785_data <- create_seurat_object(gse206785_dir, 
                                       file_pattern = "\\.txt\\.gz$", 
                                       data_type = "txt")

# Process Metastatic Datasets
# ---------------------------

# P directory - CSV format
p_dir <- "D:/data/scRNA/datas/P"
p_data <- create_seurat_object(p_dir, 
                               file_pattern = "\\.csv\\.gz$", 
                               data_type = "csv")

# Manual naming for P samples
names(p_data) <- c("GSM5573484_s19_P1", "GSM5573485_s20_P2", "GSM5573503_s38_P3")

# GSM5004187 - 10X format
gsm5004187_dir <- "D:/data/scRNA/datas/P/GSM5004187_P1"
gsm5004187_counts <- Read10X(data.dir = gsm5004187_dir)
gsm5004187_data <- CreateSeuratObject(counts = gsm5004187_counts,
                                      project = "GSM5004187_P1",
                                      min.cells = 5,
                                      min.features = 300)
names(gsm5004187_data) <- "GSM5004187_P1"

# Li directory - 10X format
li_dir <- "D:/data/scRNA/datas/P/Li"
li_data <- create_seurat_object(li_dir, data_type = "10X")

# LN directory - 10X format
ln_dir <- "D:/data/scRNA/datas/P/LN"
ln_data <- create_seurat_object(ln_dir, data_type = "10X")

# Data Integration
# ----------------

# Group datasets
primary_data <- c(gse183904_data, gse163558_data, gse206785_data)
metastatic_data <- c(p_data, list(gsm5004187_data), li_data, ln_data)

# Save individual groups
saveRDS(primary_data, "primary_data_list.rds")
saveRDS(metastatic_data, "metastatic_data_list.rds")

# Merge primary data
primary_merged <- merge_seurat_objects(primary_data, "primary_merged")

# Merge metastatic data
metastatic_merged <- merge_seurat_objects(metastatic_data, "metastatic_merged")

# Merge all data
all_data <- c(list(primary_merged), list(metastatic_merged))
all_merged <- merge_seurat_objects(all_data, "all_data_merged")

# Add dataset origin metadata
all_merged@meta.data$dataset_origin <- ifelse(
  all_merged@meta.data$orig.ident %in% names(primary_merged@meta.data$orig.ident),
  "Primary",
  "Metastatic"
)

# Add sample group metadata
all_merged@meta.data$sample_group <- case_when(
  all_merged@meta.data$orig.ident %in% c('GSM5573484_s19_P1', 'GSM5573485_s20_P2', 'GSM5573503_s38_P3') ~ 'Peritoneum_Diffuse',
  all_merged@meta.data$orig.ident == 'SeuratProject' ~ 'Peritoneum_Intestinal',
  TRUE ~ 'Primary_Gastric'
)

# Quality control metrics
all_merged[["percent.mt"]] <- PercentageFeatureSet(all_merged, pattern = "^MT-")
all_merged[["percent.rb"]] <- PercentageFeatureSet(all_merged, pattern = "^RP[SL]")
all_merged[["percent.hb"]] <- PercentageFeatureSet(all_merged, pattern = "^HB[^(P)]")

# Save final integrated dataset
saveRDS(all_merged, "integrated_scRNA_data.rds")

# Generate integration report
cat("Data Integration Summary Report\n")
cat("==============================\n\n")
cat("Primary datasets:", length(primary_data), "\n")
cat("Metastatic datasets:", length(metastatic_data), "\n")
cat("Total samples:", length(unique(all_merged@meta.data$orig.ident)), "\n")
cat("Total cells:", ncol(all_merged), "\n")
cat("Total genes:", nrow(all_merged), "\n\n")

merged_seurat <- all_merged

# QC Visualization
qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.hb")

pdf("QC_metrics.pdf", width = 12, height = 8)
VlnPlot(merged_seurat, group.by = "orig.ident", features = qc_features, 
        pt.size = 0.1, ncol = 3) + NoLegend()
VlnPlot(merged_seurat, group.by = "dataset_origin", features = qc_features, 
        pt.size = 0.1, ncol = 3) + NoLegend()
dev.off()

# Filter cells based on QC metrics
filtered_seurat <- subset(merged_seurat,
                          subset = nFeature_RNA > 200 & 
                            nFeature_RNA < 6000 & 
                            percent.mt < 20)

# Data Processing
filtered_seurat <- NormalizeData(filtered_seurat)
filtered_seurat <- FindVariableFeatures(filtered_seurat, 
                                        selection.method = "vst", 
                                        nfeatures = 2000)

# Integration using Harmony
filtered_seurat <- ScaleData(filtered_seurat, features = VariableFeatures(filtered_seurat))
filtered_seurat <- RunPCA(filtered_seurat, features = VariableFeatures(filtered_seurat))

# Batch correction
filtered_seurat <- RunHarmony(filtered_seurat, 
                              reduction = "pca",
                              group.by.vars = "orig.ident",
                              reduction.save = "harmony")

# Dimensionality reduction
filtered_seurat <- RunUMAP(filtered_seurat, reduction = "harmony", dims = 1:30)
filtered_seurat <- RunTSNE(filtered_seurat, reduction = "harmony", dims = 1:30)

# Save processed data
saveRDS(filtered_seurat, "processed_scRNA_data.rds")

# Clustering
filtered_seurat <- FindNeighbors(filtered_seurat, reduction = "harmony", dims = 1:30)
filtered_seurat <- FindClusters(filtered_seurat, resolution = 0.4)

# Find cluster markers
cluster_markers <- FindAllMarkers(filtered_seurat, 
                                  only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  logfc.threshold = 0.25)

write.csv(cluster_markers, "cluster_markers.csv")

# Visualization
pdf("clustering_results.pdf", width = 12, height = 10)
DimPlot(filtered_seurat, reduction = "umap", label = TRUE, group.by = "seurat_clusters")
DimPlot(filtered_seurat, reduction = "umap", label = TRUE, group.by = "dataset_origin")
DimPlot(filtered_seurat, reduction = "tsne", label = TRUE, group.by = "seurat_clusters")
DimPlot(filtered_seurat, reduction = "tsne", label = TRUE, group.by = "dataset_origin")
dev.off()

# Save final results
saveRDS(filtered_seurat, "final_clustered_data.rds")

cat("Analysis completed successfully!\n")
cat("Dataset composition:\n")
print(table(filtered_seurat@meta.data$dataset_origin))