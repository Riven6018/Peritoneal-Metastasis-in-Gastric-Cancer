

# Clean environment and set working directory
gc()
rm(list = ls())
options(stringsAsFactors = FALSE)
setwd("D:/data/diffuse/")

# Load required libraries
library(GEOquery)
library(limma)
library(AnnoProbe)
library(dplyr)
library(tidyr)
library(tibble)
library(survival)
library(survminer)
library(ggplot2)
library(gridExtra)
library(data.table)

# Step 1: GEO Data Download and Processing
# ----------------------------------------

download_geo_data <- function(geo_id, dest_dir = ".") {
  cat("Downloading GEO dataset:", geo_id, "\n")
  
  tryCatch({
    # Download GEO data
    geo_data <- getGEO(geo_id, 
                       destdir = dest_dir, 
                       getGPL = FALSE)
    
    # Extract expression matrix
    if (length(geo_data) > 0) {
      expr_matrix <- exprs(geo_data[[1]])
      return(expr_matrix)
    } else {
      stop("No data found for GEO ID: ", geo_id)
    }
  }, error = function(e) {
    cat("Error downloading", geo_id, ":", e$message, "\n")
    return(NULL)
  })
}

# Process expression matrix with annotation
process_expression_matrix <- function(expr_matrix, platform = "GPL570") {
  cat("Processing expression matrix with platform:", platform, "\n")
  
  # Get probe annotation
  probe_annotation <- AnnoProbe::idmap(platform, destdir = tempdir())
  probe_annotation <- dplyr::rename(probe_annotation, ID_REF = probe_id, symbol = symbol)
  
  # Convert expression matrix to data frame
  expr_df <- as.data.frame(expr_matrix)
  expr_df$ID_REF <- rownames(expr_df)
  
  # Merge with annotation
  annotated_expr <- merge(probe_annotation, expr_df, by = "ID_REF")
  
  # Remove probe ID column and aggregate by gene symbol
  annotated_expr <- annotated_expr[, -1]  # Remove ID_REF
  
  # Aggregate by gene symbol (keep maximum expression for each gene)
  aggregated_expr <- annotated_expr %>%
    group_by(symbol) %>%
    summarise_all(max, na.rm = TRUE) %>%
    filter(symbol != "") %>%  # Remove empty symbols
    column_to_rownames("symbol")
  
  return(aggregated_expr)
}

# Download and process GSE62254
gse62254_expr <- download_geo_data("GSE62254")
if (!is.null(gse62254_expr)) {
  gse62254_processed <- process_expression_matrix(gse62254_expr, "GPL570")
  write.csv(gse62254_processed, "GSE62254_matrix_processed.csv")
  cat("GSE62254 processed and saved. Dimensions:", dim(gse62254_processed), "\n")
}

# Download and process GSE66229
gse66229_expr <- download_geo_data("GSE66229")
if (!is.null(gse66229_expr)) {
  gse66229_processed <- process_expression_matrix(gse66229_expr, "GPL570")
  write.csv(gse66229_processed, "GSE66229_matrix_processed.csv")
  cat("GSE66229 processed and saved. Dimensions:", dim(gse66229_processed), "\n")
}

# Step 2: Clinical Data Processing
# --------------------------------

process_clinical_data <- function(clin_file) {
  cat("Processing clinical data from:", clin_file, "\n")
  
  clin_data <- read.csv(clin_file, stringsAsFactors = FALSE)
  
  # Process clinical data
  clin_processed <- clin_data %>%
    select(GEO_ID, sex, age, Death, OS.m, DFS.m, Stage, `stage(TNM)`, 
           Lauren, Pathology, 
           `EBV ISH\n0: negative\n1: positive\nNA: Not available (no available tumor blocks)`, 
           pnode, `# of positive node (+)`) %>%
    mutate(
      OS = OS.m * 30,   # Convert months to days
      DFS = DFS.m * 30, # Convert months to days
      EBV_Status = case_when(
        `EBV ISH\n0: negative\n1: positive\nNA: Not available (no available tumor blocks)` == "1" ~ "Positive",
        `EBV ISH\n0: negative\n1: positive\nNA: Not available (no available tumor blocks)` == "0" ~ "Negative",
        TRUE ~ "Not Available"
      ),
      Lauren_Type = recode(Lauren, 
                           "1" = "Intestinal", 
                           "2" = "Diffuse", 
                           "3" = "Mixed"),
      Positive_Lymph_Node_Ratio = ifelse(pnode != 0, `# of positive node (+)` / pnode, NA)
    ) %>%
    select(
      Sample = GEO_ID, Gender = sex, Age = age, OS, DFS, Status = Death,
      Stage, pTNM = `stage(TNM)`, Lauren = Lauren_Type, EBV_Status, 
      Positive_Lymph_Node_Ratio
    )
  
  # Remove rows with missing critical data
  clin_clean <- clin_processed %>%
    filter(!is.na(OS), !is.na(DFS), !is.na(Status))
  
  cat("Clinical data processed. Samples:", nrow(clin_clean), "\n")
  return(clin_clean)
}

# Process clinical data
clin_final <- process_clinical_data("GSE62254_clin.csv")
write.csv(clin_final, "GSE62254_clin_processed.csv")

# Step 3: Survival Analysis Functions
# -----------------------------------

perform_gene_survival_analysis <- function(expr_matrix, clin_data, genes, 
                                           low_quantile = 0.4, high_quantile = 0.6) {
  
  results <- list()
  survival_plots <- list()
  
  for (gene in genes) {
    cat("Analyzing survival for gene:", gene, "\n")
    
    if (gene %in% rownames(expr_matrix)) {
      # Extract gene expression
      gene_expr <- as.numeric(expr_matrix[gene, ])
      names(gene_expr) <- colnames(expr_matrix)
      
      # Add gene expression to clinical data
      clin_with_gene <- clin_data %>%
        mutate(Gene_Expression = gene_expr[match(Sample, colnames(expr_matrix))])
      
      # Remove samples with missing expression
      clin_filtered <- clin_with_gene %>% filter(!is.na(Gene_Expression))
      
      if (nrow(clin_filtered) > 10) {  # Ensure sufficient samples
        # Calculate expression thresholds
        low_threshold <- quantile(clin_filtered$Gene_Expression, low_quantile, na.rm = TRUE)
        high_threshold <- quantile(clin_filtered$Gene_Expression, high_quantile, na.rm = TRUE)
        
        # Create expression groups
        clin_filtered <- clin_filtered %>%
          mutate(
            Expression_Group = case_when(
              Gene_Expression > high_threshold ~ "High",
              Gene_Expression < low_threshold ~ "Low",
              TRUE ~ "Intermediate"
            )
          ) %>%
          filter(Expression_Group %in% c("High", "Low"))
        
        if (length(unique(clin_filtered$Expression_Group)) == 2) {
          # Perform survival analysis
          survival_results <- perform_survival_analysis(clin_filtered, gene)
          results[[gene]] <- survival_results
          survival_plots <- c(survival_plots, survival_results$plots)
        } else {
          warning("Insufficient groups for gene: ", gene)
        }
      } else {
        warning("Insufficient samples for gene: ", gene)
      }
    } else {
      warning("Gene not found in expression matrix: ", gene)
    }
  }
  
  return(list(results = results, plots = survival_plots))
}

perform_survival_analysis <- function(clin_data, gene_name) {
  # Ensure proper factor levels
  clin_data$Expression_Group <- factor(clin_data$Expression_Group, 
                                       levels = c("Low", "High"))
  
  # Group counts
  group_counts <- table(clin_data$Expression_Group)
  
  # OS Analysis
  os_fit <- survfit(Surv(OS, Status) ~ Expression_Group, data = clin_data)
  os_test <- survdiff(Surv(OS, Status) ~ Expression_Group, data = clin_data)
  os_pvalue <- 1 - pchisq(os_test$chisq, df = 1)
  
  # DFS Analysis
  dfs_fit <- survfit(Surv(DFS, Status) ~ Expression_Group, data = clin_data)
  dfs_test <- survdiff(Surv(DFS, Status) ~ Expression_Group, data = clin_data)
  dfs_pvalue <- 1 - pchisq(dfs_test$chisq, df = 1)
  
  # Create plots
  os_plot <- create_survival_plot(os_fit, clin_data, 
                                  paste0(gene_name, " - Overall Survival"),
                                  "Overall Survival Probability", group_counts)
  
  dfs_plot <- create_survival_plot(dfs_fit, clin_data,
                                   paste0(gene_name, " - Disease-Free Survival"),
                                   "Disease-Free Survival Probability", group_counts)
  
  # Save individual plots
  ggsave(paste0("results/", gene_name, "_OS.png"), os_plot, width = 8, height = 6, dpi = 300)
  ggsave(paste0("results/", gene_name, "_DFS.png"), dfs_plot, width = 8, height = 6, dpi = 300)
  
  return(list(
    os = list(fit = os_fit, pvalue = os_pvalue),
    dfs = list(fit = dfs_fit, pvalue = dfs_pvalue),
    plots = list(OS = os_plot, DFS = dfs_plot),
    group_counts = group_counts
  ))
}

create_survival_plot <- function(surv_fit, clin_data, title, ylab, group_counts) {
  ggsurvplot(
    surv_fit,
    data = clin_data,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.3,
    tables.y.text = FALSE,
    title = title,
    xlab = "Time (Days)",
    ylab = ylab,
    legend.title = "Expression",
    legend.labs = c(paste0("Low (n=", group_counts["Low"], ")"),
                    paste0("High (n=", group_counts["High"], ")")),
    palette = c("#2E9FDF", "#E7B800"),
    pval.coord = c(0, 0.1),
    ggtheme = theme_minimal()
  )$plot
}

# Step 4: Execute Survival Analysis
# ---------------------------------

# Create results directory
if (!dir.exists("results")) {
  dir.create("results")
}

# Define genes of interest
target_genes <- c("MFAP5", "CXCL14", "TFF2", "SYT4", "IGFN1")

# Perform survival analysis
if (exists("gse62254_processed") && exists("clin_final")) {
  survival_results <- perform_gene_survival_analysis(
    gse62254_processed, clin_final, target_genes
  )
  
  # Create combined plot
  if (length(survival_results$plots) > 0) {
    combined_plot <- grid.arrange(grobs = survival_results$plots, ncol = 2)
    ggsave("results/combined_survival_plots.png", combined_plot, 
           width = 16, height = 12, dpi = 300)
    ggsave("results/combined_survival_plots.pdf", combined_plot, 
           width = 16, height = 12)
  }
}

# Step 5: Generate Summary Report
# -------------------------------

generate_summary_report <- function(survival_results, clin_data) {
  cat("\n" + strrep("=", 70) + "\n")
  cat("GEO SURVIVAL ANALYSIS SUMMARY REPORT\n")
  cat(strrep("=", 70) + "\n\n")
  
  cat("STUDY POPULATION:\n")
  cat("Total clinical samples:", nrow(clin_data), "\n")
  cat("Samples with complete survival data:", 
      sum(!is.na(clin_data$OS) & !is.na(clin_data$DFS) & !is.na(clin_data$Status)), "\n\n")
  
  cat("GENES ANALYZED:\n")
  cat(paste(target_genes, collapse = ", "), "\n\n")
  
  if (!is.null(survival_results$results)) {
    cat("SURVIVAL ANALYSIS RESULTS:\n")
    for (gene in names(survival_results$results)) {
      result <- survival_results$results[[gene]]
      cat("•", gene, ":\n")
      cat("  OS P-value:", format.pval(result$os$pvalue, digits = 3), "\n")
      cat("  DFS P-value:", format.pval(result$dfs$pvalue, digits = 3), "\n")
      cat("  Samples (Low/High):", paste(result$group_counts, collapse = "/"), "\n\n")
    }
  }
  
  cat("OUTPUT FILES:\n")
  cat("✓ Individual gene survival plots in 'results/' directory\n")
  cat("✓ Combined survival plots: combined_survival_plots.png/pdf\n")
  cat("✓ Processed expression matrices: GSE*_matrix_processed.csv\n")
  cat("✓ Processed clinical data: GSE62254_clin_processed.csv\n\n")
}

# Generate report
generate_summary_report(survival_results, clin_final)

# Save workspace
save.image("geo_survival_analysis_workspace.RData")

cat("Analysis completed successfully!\n")
cat("Workspace saved to: geo_survival_analysis_workspace.RData\n")