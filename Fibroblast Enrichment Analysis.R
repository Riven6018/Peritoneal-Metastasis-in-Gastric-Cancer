# Fibroblast Subpopulation Enrichment Analysis

# Load required libraries
library(Seurat)
library(dplyr)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(GSVA)
library(GSEABase)
library(DOSE)
library(ggplot2)
library(pheatmap)
library(ggsci)

# Load fibroblast data
fibroblast_data <- readRDS("fibroblast_analyzed.rds")

# Set default assay to RNA
DefaultAssay(fibroblast_data) <- "RNA"

# Find marker genes for all clusters
all_markers <- FindAllMarkers(fibroblast_data,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)

# Get top 10 markers per cluster
top10_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Calculate average expression
average_expression <- AverageExpression(fibroblast_data,
                                        features = top10_markers$gene,
                                        group.by = 'RNA_snn_res.0.08',
                                        slot = 'data')$RNA

# Correlation heatmap of top variable genes
variable_genes <- names(tail(sort(apply(average_expression, 1, sd)), 1000))
pdf("fibroblast_correlation_heatmap.pdf", width = 8, height = 6)
pheatmap(cor(average_expression[variable_genes, ]),
         main = "Correlation Heatmap of Fibroblast Subclusters")
dev.off()

# GSEA Analysis for Specific Cluster
# ----------------------------------
cluster_of_interest <- "4"  # Change this to your cluster of interest

# Prepare gene list for GSEA
gsea_markers <- all_markers %>%
  dplyr::filter(cluster == cluster_of_interest) %>%
  arrange(desc(avg_log2FC))

gene_list <- gsea_markers$avg_log2FC
names(gene_list) <- rownames(gsea_markers)
gene_list <- sort(gene_list, decreasing = TRUE)

# Load gene sets
msigdb_df <- msigdbr(species = 'Homo sapiens', category = "C2")

# Perform GSEA
gsea_results <- GSEA(gene_list,
                     TERM2GENE = msigdb_df[, c('gs_name', 'gene_symbol')],
                     pvalueCutoff = 0.05,
                     verbose = FALSE)

# Save GSEA results
write.csv(gsea_results@result, "fibroblast_gsea_results.csv")

# GSEA Visualization
pdf("fibroblast_gsea_plots.pdf", width = 12, height = 10)

# Dot plot
print(dotplot(gsea_results, showCategory = 15,
              title = paste("GSEA - Cluster", cluster_of_interest)))

# Split by sign
print(dotplot(gsea_results, split = '.sign') +
        facet_wrap(~.sign, scales = 'free') +
        ggtitle("GSEA - Activated vs Suppressed Pathways"))

# Top pathways by NES
gsea_tidy <- gsea_results@result %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(p.adjust < 0.05)

if (nrow(gsea_tidy) > 0) {
  top_pathways <- head(gsea_tidy, 20)
  print(
    ggplot(top_pathways, aes(reorder(Description, NES), NES)) +
      geom_col(aes(fill = NES)) +
      coord_flip() +
      scale_fill_gsea() +
      labs(x = "Pathway", y = "Normalized Enrichment Score",
           title = paste("Top GSEA Pathways - Cluster", cluster_of_interest)) +
      theme_minimal()
  )
}

dev.off()

# GSVA Analysis for All Clusters
# ------------------------------
# Get Reactome gene sets
reactome_genesets <- msigdbr(species = "Homo sapiens",
                             category = "C2",
                             subcategory = "CP:REACTOME")

# Prepare expression matrix
expression_matrix <- GetAssayData(fibroblast_data, slot = "data")
expression_matrix <- as.matrix(expression_matrix)

# Filter genes with sufficient expression
expressed_genes <- rownames(expression_matrix)[
  rowMeans(expression_matrix) > 0.1
]
expression_matrix <- expression_matrix[expressed_genes, ]

# Perform GSVA
gsva_results <- gsva(expression_matrix,
                     reactome_genesets[, c('gs_name', 'gene_symbol')],
                     min.sz = 20,
                     verbose = FALSE)

# Heatmap of GSVA results
pdf("fibroblast_gsva_heatmap.pdf", width = 10, height = 8)
pheatmap(gsva_results,
         show_rownames = FALSE,
         main = "GSVA Enrichment Scores - Reactome Pathways")
dev.off()

# Find cluster-specific pathways
pathway_scores <- gsva_results
specific_pathways <- data.frame()

for (cluster in colnames(pathway_scores)) {
  cluster_score <- pathway_scores[, cluster]
  other_scores <- rowMeans(pathway_scores[, setdiff(colnames(pathway_scores), cluster)])
  specificity_score <- cluster_score - other_scores
  
  specific_pathways <- rbind(specific_pathways,
                             data.frame(pathway = rownames(pathway_scores),
                                        cluster = cluster,
                                        specificity = specificity_score))
}

# Top specific pathways per cluster
top_specific_pathways <- specific_pathways %>%
  group_by(cluster) %>%
  top_n(5, specificity) %>%
  ungroup()

# Heatmap of top specific pathways
top_pathways_matrix <- pathway_scores[unique(top_specific_pathways$pathway), ]
rownames(top_pathways_matrix) <- gsub("REACTOME_", "", rownames(top_pathways_matrix))
rownames(top_pathways_matrix) <- substr(rownames(top_pathways_matrix), 1, 30)

pdf("fibroblast_specific_pathways_heatmap.pdf", width = 12, height = 10)
pheatmap(top_pathways_matrix,
         main = "Top Specific Pathways per Fibroblast Subcluster")
dev.off()

# GO and KEGG Enrichment Analysis
# -------------------------------
# Get top 200 markers per cluster
top200_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 200, wt = avg_log2FC)

top200_genes <- split(top200_markers$gene, top200_markers$cluster)

# Convert to ENTREZ IDs
top200_entrez <- lapply(top200_genes, function(genes) {
  entrez_ids <- bitr(genes,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = "org.Hs.eg.db")
  return(entrez_ids$ENTREZID)
})

# GO Enrichment Analysis
go_results <- compareCluster(top200_entrez,
                             fun = "enrichGO",
                             OrgDb = "org.Hs.eg.db",
                             ont = "ALL",
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.2)

# Save GO results
write.csv(go_results@result, "fibroblast_go_enrichment.csv")

# GO Dot Plots
pdf("fibroblast_go_enrichment.pdf", width = 15, height = 12)

# Biological Process
bp_plot <- dotplot(go_results, showCategory = 5, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("GO Enrichment - Biological Process")

print(bp_plot)

# Separate plots for each ontology
cc_plot <- dotplot(go_results %>% filter(ONTOLOGY == "CC"),
                   showCategory = 5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("GO Enrichment - Cellular Component")

mf_plot <- dotplot(go_results %>% filter(ONTOLOGY == "MF"),
                   showCategory = 5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("GO Enrichment - Molecular Function")

print(cc_plot)
print(mf_plot)

dev.off()

# KEGG Enrichment Analysis
kegg_results <- compareCluster(top200_entrez,
                               fun = "enrichKEGG",
                               organism = "hsa",
                               pvalueCutoff = 0.05)

# Save KEGG results
write.csv(kegg_results@result, "fibroblast_kegg_enrichment.csv")

# KEGG Dot Plot
pdf("fibroblast_kegg_enrichment.pdf", width = 12, height = 8)
kegg_plot <- dotplot(kegg_results, showCategory = 5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("KEGG Pathway Enrichment")
print(kegg_plot)
dev.off()

# Cluster-specific KEGG analysis
cluster_specific_kegg <- enrichKEGG(gene = top200_entrez[[cluster_of_interest]],
                                    organism = "hsa",
                                    pvalueCutoff = 0.05)

if (!is.null(cluster_specific_kegg)) {
  cluster_specific_kegg <- setReadable(cluster_specific_kegg,
                                       OrgDb = org.Hs.eg.db,
                                       keyType = "ENTREZID")
  
  pdf(paste0("fibroblast_cluster_", cluster_of_interest, "_kegg.pdf"),
      width = 10, height = 6)
  print(dotplot(cluster_specific_kegg, showCategory = 15) +
          ggtitle(paste("KEGG Enrichment - Cluster", cluster_of_interest)))
  dev.off()
  
  write.csv(cluster_specific_kegg@result,
            paste0("fibroblast_cluster_", cluster_of_interest, "_kegg.csv"))
}

# Save analysis results
saveRDS(list(
  gsea_results = gsea_results,
  go_results = go_results,
  kegg_results = kegg_results,
  gsva_results = gsva_results
), "fibroblast_enrichment_results.rds")

# Generate analysis summary
cat("Fibroblast Enrichment Analysis Summary\n")
cat("=====================================\n\n")
cat("Total clusters analyzed:", length(unique(all_markers$cluster)), "\n")
cat("Total marker genes:", nrow(all_markers), "\n")
cat("GSEA significant pathways:", sum(gsea_results@result$p.adjust < 0.05), "\n")
cat("GO enriched terms:", sum(go_results@result$p.adjust < 0.05), "\n")
cat("KEGG enriched pathways:", sum(kegg_results@result$p.adjust < 0.05), "\n\n")

# Save session information
writeLines(capture.output(sessionInfo()), "enrichment_analysis_session_info.txt")

cat("Enrichment analysis completed successfully!\n")