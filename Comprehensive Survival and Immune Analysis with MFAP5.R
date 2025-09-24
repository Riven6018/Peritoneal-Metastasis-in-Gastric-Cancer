
# Load required libraries
library(IOBR)
library(survminer)
library(survival)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(pheatmap)
library(tibble)
library(patchwork)

# Set working directory
setwd("D:/data/diffuse/")

# Step 1: Data Preparation and MFAP5 Survival Analysis
# ----------------------------------------------------

prepare_survival_data <- function(expression_matrix, clinical_data, gene_name = "MFAP5") {
  # Clean sample names in expression matrix
  colnames(expression_matrix) <- gsub("(^TCGA-\\w{2}-\\w{4})-.*", "\\1", colnames(expression_matrix))
  
  # Extract MFAP5 expression
  if (!gene_name %in% rownames(expression_matrix)) {
    stop(paste("Gene", gene_name, "not found in expression matrix"))
  }
  
  mfap5_expression <- as.numeric(expression_matrix[gene_name, ])
  names(mfap5_expression) <- colnames(expression_matrix)
  
  # Find common samples
  common_samples <- intersect(names(mfap5_expression), clinical_data$Sample)
  
  if (length(common_samples) == 0) {
    stop("No common samples found between expression and clinical data")
  }
  
  # Filter data
  mfap5_filtered <- mfap5_expression[common_samples]
  clinical_filtered <- clinical_data %>% filter(Sample %in% common_samples)
  
  # Create survival dataframe
  surv_data <- data.frame(
    Sample = clinical_filtered$Sample,
    time = clinical_filtered$OS,
    status = clinical_filtered$Status,
    expression = mfap5_filtered
  )
  
  return(surv_data)
}

# Prepare survival data
surv_data <- prepare_survival_data(drawdata, clinical_indexed_clean, "MFAP5")

# Determine optimal cutpoint for MFAP5
calculate_optimal_cutpoint <- function(surv_data) {
  cutpoint <- surv_cutpoint(
    surv_data,
    time = "time",
    event = "status",
    variables = "expression"
  )
  
  if (is.null(cutpoint$cutpoint)) {
    # Fallback to median if optimal cutpoint cannot be determined
    median_value <- median(surv_data$expression, na.rm = TRUE)
    warning("Optimal cutpoint not found, using median: ", median_value)
    return(median_value)
  }
  
  cutoff_value <- cutpoint$cutpoint$cutpoint
  cat("Optimal cutpoint value:", cutoff_value, "\n")
  return(cutoff_value)
}

cutoff_value <- calculate_optimal_cutpoint(surv_data)

# Create groups
surv_data$group <- ifelse(surv_data$expression >= cutoff_value, "High", "Low")
surv_data$group <- factor(surv_data$group, levels = c("Low", "High"))

# Merge with clinical data
clin_final <- merge(clinical_indexed_clean, surv_data[, c("Sample", "group")], by = "Sample", all.x = TRUE)

# Step 2: Survival Analysis Functions
# -----------------------------------
perform_survival_analysis <- function(clin_data, time_var, event_var, group_var, 
                                      title_suffix = "", output_prefix = "") {
  
  # Create survival object
  surv_object <- Surv(time = clin_data[[time_var]], event = clin_data[[event_var]])
  
  # Fit survival curve
  fit <- survfit(surv_object ~ clin_data[[group_var]])
  
  # Create survival plot
  surv_plot <- ggsurvplot(
    fit,
    data = clin_data,
    pval = TRUE,
    risk.table = TRUE,
    conf.int = TRUE,
    ggtheme = theme_minimal(),
    palette = c("#E7B800", "#2E9FDF"),
    title = paste("Survival Analysis", title_suffix)
  )
  
  # Save plots
  ggsave(paste0(output_prefix, "_survival_plot.png"), surv_plot$plot, width = 8, height = 6)
  ggsave(paste0(output_prefix, "_survival_plot.pdf"), surv_plot$plot, width = 8, height = 6)
  
  return(surv_plot)
}

# Perform OS and DFS survival analysis
os_plot <- perform_survival_analysis(clin_final, "OS", "Status", "group", 
                                     "for Overall Survival", "OS_MFAP5")

dfs_plot <- perform_survival_analysis(clin_final, "DFS", "Status", "group",
                                      "for Disease-Free Survival", "DFS_MFAP5")

# Step 3: Immune Infiltration Analysis
# ------------------------------------
perform_comprehensive_immune_analysis <- function(expression_matrix) {
  methods_list <- list()
  
  analysis_methods <- c("cibersort", "mcpcounter", "xcell", "epic",
                        "estimate", "timer", "quantiseq", "ips")
  
  for (method in analysis_methods) {
    cat("Performing", method, "analysis...\n")
    
    tryCatch({
      if (method == "timer") {
        result <- deconvo_tme(eset = expression_matrix, method = method,
                              group_list = rep("sample_type", ncol(expression_matrix)))
      } else if (method == "quantiseq") {
        result <- deconvo_tme(eset = expression_matrix, method = method,
                              arrays = TRUE, tumor = TRUE, scale_mrna = TRUE)
      } else {
        result <- deconvo_tme(eset = expression_matrix, method = method, arrays = TRUE)
      }
      
      # Ensure ID column exists
      if (!"ID" %in% colnames(result)) {
        result$ID <- rownames(result)
      }
      
      methods_list[[toupper(method)]] <- result
      
    }, error = function(e) {
      cat("Error in", method, "analysis:", e$message, "\n")
    })
  }
  
  return(methods_list)
}

# Perform immune analysis
immune_results <- perform_comprehensive_immune_analysis(GSE62254)

# Prepare grouping data
group_df <- clin_final[, c("Sample", "group")]
colnames(group_df) <- c("ID", "MFAP5_group")

# Step 4: Immune Infiltration Visualization
# -----------------------------------------
prepare_immune_long_data <- function(immune_results, group_df) {
  prepare_method_data <- function(df, method_name) {
    if (!"ID" %in% colnames(df)) {
      df$ID <- rownames(df)
    }
    
    df_long <- df %>%
      pivot_longer(cols = -ID, names_to = "CellType", values_to = "Fraction") %>%
      mutate(Method = method_name)
    
    return(df_long)
  }
  
  # Combine all methods
  all_data_long <- map2_dfr(immune_results, names(immune_results), 
                            ~ prepare_method_data(.x, .y))
  
  # Merge with grouping information
  all_data_clean <- all_data_long %>%
    inner_join(group_df, by = "ID") %>%
    filter(!is.na(Fraction), !is.na(MFAP5_group), !is.na(CellType))
  
  return(all_data_clean)
}

# Prepare data for visualization
immune_data_clean <- prepare_immune_long_data(immune_results, group_df)

# Create combined immune infiltration plot
create_combined_immune_plot <- function(immune_data) {
  ggplot(immune_data, aes(x = CellType, y = Fraction, fill = MFAP5_group)) +
    geom_boxplot(position = position_dodge(0.8), alpha = 0.7, outlier.shape = NA) +
    facet_wrap(~ Method, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("Low" = "#4DBBD5", "High" = "#E64B35")) +
    theme_minimal() +
    labs(title = "Immune Cell Composition by Analysis Method",
         x = "Immune Cell Type", 
         y = "Cell Fraction",
         fill = "MFAP5 Group") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "top") +
    stat_compare_means(aes(group = MFAP5_group), method = "t.test", 
                       label = "p.signif", size = 3)
}

# Generate and save combined plot
combined_immune_plot <- create_combined_immune_plot(immune_data_clean)
ggsave("immune_infiltration_combined.pdf", combined_immune_plot, width = 14, height = 10)

# Step 5: Individual Method Plots
# -------------------------------
generate_individual_immune_plots <- function(immune_data, output_dir = ".") {
  methods <- unique(immune_data$Method)
  
  for (method in methods) {
    method_data <- immune_data %>% filter(Method == method)
    
    # Adjust width based on number of cell types
    n_celltypes <- length(unique(method_data$CellType))
    plot_width <- ifelse(method == "XCELL", 20, min(12, ceiling(n_celltypes / 3) * 4))
    
    p <- ggboxplot(
      method_data,
      x = "CellType",
      y = "Fraction",
      color = "MFAP5_group",
      palette = c("Low" = "#4DBBD5", "High" = "#E64B35"),
      ylab = "Immune Infiltration Score",
      xlab = "",
      title = paste("Immune Cell Infiltration -", method)
    ) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      stat_compare_means(aes(group = MFAP5_group), method = "wilcox.test", 
                         label = "p.signif")
    
    # Save plots
    ggsave(paste0(output_dir, "/immune_infiltration_", method, ".png"), 
           p, width = plot_width, height = 8, dpi = 300)
    ggsave(paste0(output_dir, "/immune_infiltration_", method, ".pdf"), 
           p, width = plot_width, height = 8)
  }
}

# Generate individual plots
generate_individual_immune_plots(immune_data_clean)

# Step 6: Significant Features Heatmap
# ------------------------------------
create_immune_heatmap <- function(immune_results, group_df, p_cutoff = 0.05) {
  # Combine all immune results
  all_immune <- reduce(immune_results, full_join, by = "ID")
  
  # Merge with grouping information
  merged_data <- inner_join(all_immune, group_df, by = "ID")
  
  # Calculate p-values for each immune feature
  immune_features <- setdiff(colnames(merged_data), c("ID", "MFAP5_group"))
  
  p_values <- map_dbl(immune_features, function(feature) {
    high_group <- merged_data %>% filter(MFAP5_group == "High") %>% pull(feature)
    low_group <- merged_data %>% filter(MFAP5_group == "Low") %>% pull(feature)
    
    if (length(high_group) > 1 && length(low_group) > 1) {
      return(t.test(high_group, low_group)$p.value)
    }
    return(NA)
  })
  
  # Select significant features
  sig_features <- immune_features[p_values < p_cutoff & !is.na(p_values)]
  
  if (length(sig_features) == 0) {
    cat("No significant features found at p <", p_cutoff, "\n")
    return(NULL)
  }
  
  # Prepare heatmap data
  heatmap_data <- merged_data %>%
    select(ID, all_of(sig_features), MFAP5_group) %>%
    arrange(MFAP5_group)
  
  rownames(heatmap_data) <- heatmap_data$ID
  heatmap_matrix <- as.matrix(heatmap_data[, sig_features])
  heatmap_matrix_scaled <- scale(heatmap_matrix)
  
  # Create annotation
  annotation_df <- data.frame(MFAP5_Group = heatmap_data$MFAP5_group)
  rownames(annotation_df) <- heatmap_data$ID
  
  # Create heatmap
  pheatmap(heatmap_matrix_scaled,
           annotation_row = annotation_df,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = FALSE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = paste("Significant Immune Features (p <", p_cutoff, ")"),
           fontsize_col = 8)
}

# Generate heatmap
pdf("significant_immune_features_heatmap.pdf", width = 10, height = 8)
create_immune_heatmap(immune_results, group_df)
dev.off()

# Step 7: TIDE Analysis
# ---------------------
analyze_tide_data <- function(tide_file, group_df) {
  # Load TIDE data
  tide_data <- read.csv(tide_file, header = TRUE)
  
  # Merge with grouping information
  tide_merged <- merge(tide_data, group_df, by.x = "Patient", by.y = "ID")
  
  # Create comparison list
  my_comparisons <- list(c("Low", "High"))
  
  # Create TIDE violin plots
  tide_metrics <- c("TIDE", "Dysfunction", "Exclusion", "MSI")
  plot_list <- list()
  
  for (metric in tide_metrics) {
    if (metric %in% colnames(tide_merged)) {
      p <- ggviolin(tide_merged, x = "MFAP5_group", y = metric, fill = "MFAP5_group",
                    palette = c("#F7A7A8", "#7E9BC7"),
                    add = "boxplot", add.params = list(fill = "white")) +
        stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
                           method = "t.test") +
        labs(title = paste(metric, "by MFAP5 Expression"))
      
      plot_list[[metric]] <- p
    }
  }
  
  # Combine plots
  combined_tide_plot <- wrap_plots(plot_list, ncol = 2)
  ggsave("TIDE_analysis_MFAP5.pdf", combined_tide_plot, width = 16, height = 12)
  
  return(combined_tide_plot)
}

# Perform TIDE analysis (if TIDE data available)
# tide_plot <- analyze_tide_data("GSE62254_TIDE.csv", group_df)

# Step 8: Save Results and Summary
# --------------------------------
# Save grouping information
write.csv(group_df, "MFAP5_optimal_cutpoint_groups.csv", row.names = FALSE)

# Generate summary report
generate_analysis_summary <- function(surv_data, immune_results, group_df) {
  cat("MFAP5 Survival and Immune Analysis Summary\n")
  cat("==========================================\n\n")
  
  cat("Survival Analysis:\n")
  cat("Total samples:", nrow(surv_data), "\n")
  cat("MFAP5 High group:", sum(surv_data$group == "High"), "\n")
  cat("MFAP5 Low group:", sum(surv_data$group == "Low"), "\n")
  cat("Optimal cutpoint:", cutoff_value, "\n\n")
  
  cat("Immune Analysis:\n")
  cat("Methods completed:", paste(names(immune_results), collapse = ", "), "\n")
  cat("Total immune features analyzed:", ncol(reduce(immune_results, full_join, by = "ID")) - 1, "\n")
}

# Generate summary
generate_analysis_summary(surv_data, immune_results, group_df)

# Save workspace
save.image("MFAP5_survival_immune_analysis.RData")

cat("Analysis completed successfully!\n")