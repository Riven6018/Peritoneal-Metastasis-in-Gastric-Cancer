# MFAP5-CAFs Survival Analysis

# Load required libraries
library(data.table)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)

# Read and prepare data
# ---------------------
GSE62254 <- fread("GSE62254_matrix_去重后.csv", encoding = "UTF-8")
CAFs <- fread("CIBERSORTx_Job40_Results.txt", encoding = "UTF-8")
clin <- fread("GSE62254_clin_整理后.csv", stringsAsFactors = FALSE)

# Check data structure
cat("Data Structure Overview:\n")
cat("GSE62254 dimensions:", dim(GSE62254), "\n")
cat("CAFs dimensions:", dim(CAFs), "\n")
cat("Clinical data dimensions:", dim(clin), "\n\n")

# Prepare clinical data column names
setnames(clin, c("V1", "Sample", "Gender", "Age", "OS_days", "DFS_days", 
                 "Death_Status", "Stage", "pTNM", "Lauren", "EBV_Status"))

# Merge CAFs data with clinical data
merged_data <- merge(
  CAFs[, .(Sample = Mixture, MFAP5_CAFs = `MFAP5-CAFs`)], 
  clin, 
  by = "Sample"
)

cat("Merged data dimensions:", dim(merged_data), "\n")
cat("Sample overlap:", nrow(merged_data), "samples\n\n")

# Data Quality Check
# ------------------
cat("Data Quality Check:\n")
cat("MFAP5-CAFs summary statistics:\n")
print(summary(merged_data$MFAP5_CAFs))

cat("\nMissing values:\n")
cat("MFAP5-CAFs NA:", sum(is.na(merged_data$MFAP5_CAFs)), "\n")
cat("OS_days NA:", sum(is.na(merged_data$OS_days)), "\n")
cat("DFS_days NA:", sum(is.na(merged_data$DFS_days)), "\n")
cat("Death_Status NA:", sum(is.na(merged_data$Death_Status)), "\n")

# Remove rows with missing survival data
initial_rows <- nrow(merged_data)
merged_data <- merged_data[!is.na(OS_days) & !is.na(DFS_days) & !is.na(Death_Status)]
cat("\nRows removed due to missing survival data:", initial_rows - nrow(merged_data), "\n")

# Survival Analysis Functions
# ---------------------------
create_survival_plot <- function(surv_object, data, group_var, title, ylab, 
                                 group_counts, output_file) {
  
  # Create survival plot
  surv_plot <- ggsurvplot(
    surv_object,
    data = data,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = FALSE,
    tables.theme = theme_cleantable(),
    title = title,
    xlab = "Time (Days)",
    ylab = ylab,
    legend.title = "MFAP5-CAFs",
    legend.labs = c(
      paste0("Positive (n=", group_counts["Expression"], ")"), 
      paste0("Negative (n=", group_counts["No Expression"], ")")
    ),
    palette = c("#E7B800", "#2E9FDF"),
    pval.coord = c(0, 0.1),
    pval.size = 4,
    font.main = 14,
    font.x = 12,
    font.y = 12,
    font.legend = 10,
    font.tickslab = 10
  )
  
  # Customize risk table
  surv_plot$table <- surv_plot$table + 
    theme(
      plot.title = element_text(size = 12),
      axis.title.x = element_text(size = 10),
      axis.text.x = element_text(size = 8)
    )
  
  # Save plot
  ggsave(output_file, plot = surv_plot$plot, width = 8, height = 7, dpi = 300)
  
  return(surv_plot)
}

perform_survival_analysis <- function(merged_data) {
  # Define groups based on MFAP5-CAFs expression
  merged_data[, Group := factor(
    ifelse(MFAP5_CAFs == 0, "No Expression", "Expression"),
    levels = c("No Expression", "Expression")
  )]
  
  # Check group distribution
  group_counts <- table(merged_data$Group)
  cat("\nGroup Distribution:\n")
  print(group_counts)
  
  # Ensure sufficient samples in each group
  if (any(group_counts < 5)) {
    warning("One or more groups have very few samples (<5). Results may be unreliable.")
  }
  
  # Overall Survival (OS) Analysis
  # ------------------------------
  cat("\n=== Overall Survival Analysis ===\n")
  
  # Create survival object for OS
  surv_os <- Surv(time = merged_data$OS_days, event = merged_data$Death_Status)
  fit_os <- survfit(surv_os ~ Group, data = merged_data)
  
  # Log-rank test for OS
  os_test <- survdiff(surv_os ~ Group, data = merged_data)
  os_pvalue <- 1 - pchisq(os_test$chisq, df = length(os_test$n) - 1)
  
  cat("OS Log-rank test:\n")
  cat("Chi-squared:", round(os_test$chisq, 3), "\n")
  cat("Degrees of freedom:", length(os_test$n) - 1, "\n")
  cat("P-value:", format.pval(os_pvalue, digits = 3), "\n")
  
  # Create OS plot
  os_plot <- create_survival_plot(
    surv_object = fit_os,
    data = merged_data,
    group_var = "Group",
    title = "Overall Survival by MFAP5-CAFs Expression",
    ylab = "Overall Survival Probability",
    group_counts = group_counts,
    output_file = "results/MFAP5_OS_Survival.pdf"
  )
  
  # Disease-Free Survival (DFS) Analysis
  # ------------------------------------
  cat("\n=== Disease-Free Survival Analysis ===\n")
  
  # Create survival object for DFS
  surv_dfs <- Surv(time = merged_data$DFS_days, event = merged_data$Death_Status)
  fit_dfs <- survfit(surv_dfs ~ Group, data = merged_data)
  
  # Log-rank test for DFS
  dfs_test <- survdiff(surv_dfs ~ Group, data = merged_data)
  dfs_pvalue <- 1 - pchisq(dfs_test$chisq, df = length(dfs_test$n) - 1)
  
  cat("DFS Log-rank test:\n")
  cat("Chi-squared:", round(dfs_test$chisq, 3), "\n")
  cat("Degrees of freedom:", length(dfs_test$n) - 1, "\n")
  cat("P-value:", format.pval(dfs_pvalue, digits = 3), "\n")
  
  # Create DFS plot
  dfs_plot <- create_survival_plot(
    surv_object = fit_dfs,
    data = merged_data,
    group_var = "Group",
    title = "Disease-Free Survival by MFAP5-CAFs Expression",
    ylab = "Disease-Free Survival Probability",
    group_counts = group_counts,
    output_file = "results/MFAP5_DFS_Survival.pdf"
  )
  
  # Additional Summary Statistics
  # -----------------------------
  cat("\n=== Additional Statistics ===\n")
  
  # Median survival times
  os_median <- surv_median(fit_os)
  dfs_median <- surv_median(fit_dfs)
  
  cat("Median OS times:\n")
  print(os_median)
  
  cat("\nMedian DFS times:\n")
  print(dfs_median)
  
  # Return comprehensive results
  results <- list(
    data = merged_data,
    group_counts = group_counts,
    os_analysis = list(
      fit = fit_os,
      test = os_test,
      p_value = os_pvalue,
      median_times = os_median,
      plot = os_plot
    ),
    dfs_analysis = list(
      fit = fit_dfs,
      test = dfs_test,
      p_value = dfs_pvalue,
      median_times = dfs_median,
      plot = dfs_plot
    ),
    summary_stats = list(
      total_samples = nrow(merged_data),
      mfap5_summary = summary(merged_data$MFAP5_CAFs),
      os_summary = summary(merged_data$OS_days),
      dfs_summary = summary(merged_data$DFS_days)
    )
  )
  
  return(results)
}

# Create results directory
if (!dir.exists("results")) {
  dir.create("results")
}

# Perform survival analysis
analysis_results <- perform_survival_analysis(merged_data)

# Save numerical results
# ----------------------
# Save group information
fwrite(merged_data[, .(Sample, MFAP5_CAFs, Group)], 
       "results/MFAP5_group_assignments.csv")

# Save statistical results
stats_summary <- data.frame(
  Analysis = c("Overall Survival", "Disease-Free Survival"),
  ChiSquared = c(analysis_results$os_analysis$test$chisq, 
                 analysis_results$dfs_analysis$test$chisq),
  PValue = c(analysis_results$os_analysis$p_value, 
             analysis_results$dfs_analysis$p_value),
  Samples_Positive = analysis_results$group_counts["Expression"],
  Samples_Negative = analysis_results$group_counts["No Expression"]
)

fwrite(stats_summary, "results/survival_analysis_statistics.csv")

# Generate Summary Report
# -----------------------
cat("\n" + strrep("=", 60) + "\n")
cat("MFAP5-CAFs SURVIVAL ANALYSIS SUMMARY REPORT\n")
cat(strrep("=", 60) + "\n\n")

cat("STUDY POPULATION:\n")
cat("Total analyzed samples:", analysis_results$summary_stats$total_samples, "\n")
cat("MFAP5-CAFs positive:", analysis_results$group_counts["Expression"], "\n")
cat("MFAP5-CAFs negative:", analysis_results$group_counts["No Expression"], "\n\n")

cat("SURVIVAL OUTCOMES:\n")
cat("Overall Survival - P-value:", 
    format.pval(analysis_results$os_analysis$p_value, digits = 3), "\n")
cat("Disease-Free Survival - P-value:", 
    format.pval(analysis_results$dfs_analysis$p_value, digits = 3), "\n\n")

cat("DATA QUALITY:\n")
cat("MFAP5-CAFs range:", 
    round(min(merged_data$MFAP5_CAFs, na.rm = TRUE), 3), "to",
    round(max(merged_data$MFAP5_CAFs, na.rm = TRUE), 3), "\n")
cat("OS follow-up range:", 
    round(min(merged_data$OS_days, na.rm = TRUE), 1), "to",
    round(max(merged_data$OS_days, na.rm = TRUE), 1), "days\n")
cat("DFS follow-up range:", 
    round(min(merged_data$DFS_days, na.rm = TRUE), 1), "to",
    round(max(merged_data$DFS_days, na.rm = TRUE), 1), "days\n\n")

cat("OUTPUT FILES:\n")
cat("✓ MFAP5_OS_Survival.pdf - Overall survival plot\n")
cat("✓ MFAP5_DFS_Survival.pdf - Disease-free survival plot\n")
cat("✓ MFAP5_group_assignments.csv - Sample group assignments\n")
cat("✓ survival_analysis_statistics.csv - Statistical results\n\n")

cat("Analysis completed successfully!\n")

# Save workspace for future reference
saveRDS(analysis_results, "results/survival_analysis_results.rds")
save.image("MFAP5_survival_analysis_workspace.RData")

cat("Workspace saved to: MFAP5_survival_analysis_workspace.RData\n")