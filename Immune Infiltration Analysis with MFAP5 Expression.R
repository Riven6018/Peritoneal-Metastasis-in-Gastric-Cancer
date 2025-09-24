# Immune Infiltration Analysis with MFAP5 Expression


# Load required libraries
library(IOBR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(pheatmap)
library(tibble)

# Set working directory
setwd("D:/data/diffuse/")

# Step 1: Data Preparation and MFAP5 Grouping
# -------------------------------------------
prepare_mfap5_groups <- function(expression_matrix, gene_name = "MFAP5") {
  # Extract MFAP5 expression
  if (gene_name %in% rownames(expression_matrix)) {
    mfap5_expression <- as.numeric(expression_matrix[gene_name, ])
  } else {
    stop(paste("Gene", gene_name, "not found in expression matrix"))
  }
  
  # Create groups based on median expression
  median_expression <- median(mfap5_expression, na.rm = TRUE)
  groups <- ifelse(mfap5_expression > median_expression, "High", "Low")
  
  # Create grouping dataframe
  group_df <- data.frame(
    ID = colnames(expression_matrix),
    MFAP5_expression = mfap5_expression,
    MFAP5_group = factor(groups, levels = c("Low", "High"))
  )
  
  return(group_df)
}

# Apply grouping
group_df <- prepare_mfap5_groups(GSE62254, "MFAP5")

# Step 2: Immune Infiltration Analysis
# ------------------------------------
perform_immune_analysis <- function(expression_matrix) {
  methods_list <- list()
  
  # Define analysis methods
  analysis_methods <- c("cibersort", "mcpcounter", "xcell", "epic", 
                        "estimate", "timer", "quantiseq", "ips")
  
  for (method in analysis_methods) {
    cat("Performing", method, "analysis...\n")
    
    tryCatch({
      if (method == "timer") {
        # TIMER requires group_list parameter
        result <- deconvo_tme(eset = expression_matrix, method = method,
                              group_list = rep("sample_type", ncol(expression_matrix)))
      } else if (method == "quantiseq") {
        # quanTIseq specific parameters
        result <- deconvo_tme(eset = expression_matrix, method = method,
                              arrays = TRUE, tumor = TRUE, scale_mrna = TRUE)
      } else {
        # Standard analysis
        result <- deconvo_tme(eset = expression_matrix, method = method,
                              arrays = TRUE)
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

# Perform immune infiltration analysis
immune_results <- perform_immune_analysis(TPMexp)

# Step 3: Visualization Functions
# -------------------------------
create_infiltration_plot <- function(data_long, method_name) {
  ggplot(data_long, aes(x = MFAP5_group, y = Fraction, fill = MFAP5_group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 0.5) +
    facet_wrap(~CellType, scales = "free_y", ncol = 4) +
    theme_minimal() +
    scale_fill_manual(values = c("Low" = "#4DBBD5", "High" = "#E64B35")) +
    labs(title = paste("Immune Cell Fractions -", method_name),
         x = "MFAP5 Expression Group",
         y = "Cell Fraction") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    stat_compare_means(method = "t.test", label = "p.signif", size = 3)
}

prepare_long_data <- function(df, method_name) {
  if (!"ID" %in% colnames(df)) {
    df$ID <- rownames(df)
  }
  
  df_long <- df %>%
    pivot_longer(cols = -ID, names_to = "CellType", values_to = "Fraction") %>%
    mutate(Method = method_name)
  
  return(df_long)
}

# Step 4: Generate Individual Method Plots
# ----------------------------------------
generate_individual_plots <- function(immune_results, group_df, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (method_name in names(immune_results)) {
    cat("Processing", method_name, "...\n")
    
    # Prepare data
    method_data <- immune_results[[method_name]]
    data_long <- prepare_long_data(method_data, method_name)
    data_long <- inner_join(data_long, group_df, by = "ID")
    
    # Create plot
    p <- create_infiltration_plot(data_long, method_name)
    
    # Adjust dimensions based on number of cell types
    n_celltypes <- length(unique(data_long$CellType))
    plot_width <- max(10, ceiling(n_celltypes / 4) * 3)
    plot_height <- 8
    
    # Save plots
    ggsave(file.path(output_dir, paste0(method_name, "_immune_infiltration.pdf")),
           plot = p, width = plot_width, height = plot_height)
    
    ggsave(file.path(output_dir, paste0(method_name, "_immune_infiltration.png")),
           plot = p, width = plot_width, height = plot_height, dpi = 300)
  }
}

# Generate individual plots
output_dir <- "D:/data/diffuse/immune_plots"
generate_individual_plots(immune_results, group_df, output_dir)

# Step 5: Combined Analysis
# -------------------------
# Combine all immune data
all_data_long <- map2_dfr(immune_results, names(immune_results), 
                          ~ prepare_long_data(.x, .y))

# Merge with grouping information
all_data_clean <- all_data_long %>%
  inner_join(group_df, by = "ID") %>%
  filter(!is.na(Fraction), !is.na(MFAP5_group), !is.na(CellType))

# Create combined plot
combined_plot <- ggplot(all_data_clean, 
                        aes(x = CellType, y = Fraction, fill = MFAP5_group)) +
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

# Save combined plot
ggsave("combined_immune_infiltration.pdf", combined_plot, width = 16, height = 12)
ggsave("combined_immune_infiltration.png", combined_plot, width = 16, height = 12, dpi = 300)

# Step 6: Significant Features Heatmap
# ------------------------------------
create_significance_heatmap <- function(immune_results, group_df, p_cutoff = 0.05) {
  # Merge all results
  all_results <- reduce(immune_results, function(x, y) {
    full_join(x, y, by = "ID", suffix = c("", ""))
  })
  
  # Merge with grouping information
  merged_df <- inner_join(all_results, group_df, by = "ID")
  
  # Remove ID and group columns for analysis
  immune_features <- setdiff(colnames(merged_df), c("ID", "MFAP5_group", "MFAP5_expression"))
  
  # Calculate p-values for each feature
  p_values <- map_dbl(immune_features, function(feature) {
    high_group <- merged_df %>% filter(MFAP5_group == "High") %>% pull(feature)
    low_group <- merged_df %>% filter(MFAP5_group == "Low") %>% pull(feature)
    
    if (length(high_group) > 1 && length(low_group) > 1) {
      test_result <- t.test(high_group, low_group)
      return(test_result$p.value)
    } else {
      return(NA)
    }
  })
  
  # Create results dataframe
  significance_df <- data.frame(
    Feature = immune_features,
    p_value = p_values,
    significant = p_values < p_cutoff
  ) %>% filter(!is.na(p_value))
  
  # Select significant features
  sig_features <- significance_df %>% 
    filter(significant) %>% 
    pull(Feature)
  
  if (length(sig_features) == 0) {
    cat("No significant features found at p <", p_cutoff, "\n")
    return(NULL)
  }
  
  # Prepare heatmap data
  heatmap_data <- merged_df %>%
    select(ID, all_of(sig_features), MFAP5_group) %>%
    arrange(MFAP5_group)
  
  rownames(heatmap_data) <- heatmap_data$ID
  heatmap_matrix <- as.matrix(heatmap_data[, sig_features, drop = FALSE])
  
  # Scale the matrix
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
           show_colnames = TRUE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = paste("Significant Immune Features (p <", p_cutoff, ")"),
           fontsize_col = 8,
           fontsize_row = 6)
}

# Generate significance heatmap
pdf("significant_immune_features_heatmap.pdf", width = 10, height = 8)
create_significance_heatmap(immune_results, group_df, p_cutoff = 0.05)
dev.off()

# Step 7: Summary Statistics
# --------------------------
generate_summary_report <- function(immune_results, group_df) {
  cat("Immune Infiltration Analysis Summary\n")
  cat("====================================\n\n")
  
  cat("Sample Information:\n")
  cat("Total samples:", nrow(group_df), "\n")
  cat("MFAP5 High group:", sum(group_df$MFAP5_group == "High"), "\n")
  cat("MFAP5 Low group:", sum(group_df$MFAP5_group == "Low"), "\n\n")
  
  cat("Analysis Methods Completed:\n")
  cat(paste(names(immune_results), collapse = ", "), "\n\n")
  
  # Calculate significant features
  all_results <- reduce(immune_results, function(x, y) {
    full_join(x, y, by = "ID", suffix = c("", ""))
  })
  
  merged_df <- inner_join(all_results, group_df, by = "ID")
  immune_features <- setdiff(colnames(merged_df), c("ID", "MFAP5_group", "MFAP5_expression"))
  
  sig_count <- sum(map_dbl(immune_features, function(feature) {
    high <- merged_df %>% filter(MFAP5_group == "High") %>% pull(feature)
    low <- merged_df %>% filter(MFAP5_group == "Low") %>% pull(feature)
    if (length(high) > 1 && length(low) > 1) {
      p_val <- t.test(high, low)$p.value
      return(p_val < 0.05)
    }
    return(FALSE)
  }))
  
  cat("Significant immune features (p < 0.05):", sig_count, "/", length(immune_features), "\n")
}

# Generate summary report
generate_summary_report(immune_results, group_df)

# Save workspace
save.image("immune_infiltration_analysis.RData")

cat("Analysis completed successfully!\n")
cat("Results saved to:", output_dir, "\n")