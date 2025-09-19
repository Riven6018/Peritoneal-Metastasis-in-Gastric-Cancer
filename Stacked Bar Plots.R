# Cell Type Composition Stacked Bar Plots
# Author: [Your Name]
# Date: [Date]

# Load required libraries
library(Seurat)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(scales)

# Load annotated Seurat object
seurat_obj <- readRDS("manually_annotated_seurat.rds")

# Prepare data for stacked bar plots
stack_data <- data.frame(
  sample = seurat_obj$orig.ident,
  cell_type = Idents(seurat_obj)
)

head(stack_data)

# Calculate cell counts per sample and cell type
cell_type_counts <- stack_data %>%
  group_by(sample, cell_type) %>%
  summarise(count = n(), .groups = 'drop')

head(cell_type_counts)

# Create sample grouping based on sample names
cell_type_counts <- cell_type_counts %>%
  mutate(sample_group = case_when(
    sample %in% c('GSM5573484_sample19', 'GSM5573485_sample20', 'GSM5573503_sample38') ~ 'Peritoneum_Diffuse',
    sample == 'SeuratProject' ~ 'Peritoneum_Intestinal',
    TRUE ~ 'Primary_gastric'
  ))

# Plot 1: Stacked bar plot (absolute counts)
plot_absolute <- ggplot(cell_type_counts, aes(x = sample, y = count, fill = cell_type)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_npg(alpha = 0.6) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 0.9),
    legend.position = "right"
  ) +
  labs(
    x = "",
    y = "Cell Count",
    title = "Cell Type Composition - Absolute Counts",
    fill = "Cell Type"
  ) +
  scale_y_continuous(labels = comma)

# Plot 2: Stacked bar plot (percentage)
plot_percentage <- ggplot(cell_type_counts, aes(x = sample, y = count, fill = cell_type)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_npg(alpha = 0.6) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 0.9),
    legend.position = "right"
  ) +
  labs(
    x = "",
    y = "Percentage",
    title = "Cell Type Composition - Percentage",
    fill = "Cell Type"
  ) +
  scale_y_continuous(labels = percent_format())

# Plot 3: Grouped by sample type (percentage)
plot_grouped <- ggplot(cell_type_counts, aes(x = sample_group, y = count, fill = cell_type)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_npg(alpha = 0.6) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 0.9),
    legend.position = "right"
  ) +
  labs(
    x = "Sample Type",
    y = "Percentage",
    title = "Cell Type Composition by Sample Group",
    fill = "Cell Type"
  ) +
  scale_y_continuous(labels = percent_format())

# Plot 4: Horizontal stacked bar plot (percentage)
plot_horizontal <- ggplot(cell_type_counts, aes(x = count, y = sample_group, fill = cell_type)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_npg(alpha = 0.6) +
  theme_classic() +
  theme(
    legend.position = "right"
  ) +
  labs(
    x = "Percentage",
    y = "Sample Group",
    title = "Cell Type Composition - Horizontal View",
    fill = "Cell Type"
  ) +
  scale_x_continuous(labels = percent_format())

# Save all plots
pdf("cell_type_composition_plots.pdf", width = 12, height = 8)
print(plot_absolute)
print(plot_percentage)
print(plot_grouped)
print(plot_horizontal)
dev.off()

# Additional detailed analysis
# ----------------------------

# Calculate summary statistics
summary_stats <- cell_type_counts %>%
  group_by(sample_group, cell_type) %>%
  summarise(
    total_cells = sum(count),
    mean_percentage = mean(count / sum(count) * 100),
    .groups = 'drop'
  )

write.csv(summary_stats, "cell_type_composition_summary.csv", row.names = FALSE)

# Create a more detailed grouped plot with facets
plot_detailed <- ggplot(cell_type_counts, aes(x = sample_group, y = count, fill = cell_type)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_npg(alpha = 0.6) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 4) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 0.9),
    legend.position = "none"
  ) +
  labs(
    x = "Sample Group",
    y = "Percentage",
    title = "Detailed Cell Type Distribution by Sample Group"
  ) +
  scale_y_continuous(labels = percent_format())

pdf("detailed_cell_type_distribution.pdf", width = 14, height = 10)
print(plot_detailed)
dev.off()

# Create a table of cell type percentages
percentage_table <- cell_type_counts %>%
  group_by(sample_group, cell_type) %>%
  summarise(total = sum(count), .groups = 'drop') %>%
  group_by(sample_group) %>%
  mutate(percentage = total / sum(total) * 100) %>%
  select(-total) %>%
  pivot_wider(names_from = sample_group, values_from = percentage, values_fill = 0)

write.csv(percentage_table, "cell_type_percentage_table.csv", row.names = FALSE)

# Generate a comprehensive report
cat("Cell Type Composition Analysis Report\n")
cat("====================================\n\n")
cat("Total samples:", length(unique(cell_type_counts$sample)), "\n")
cat("Sample groups:", paste(unique(cell_type_counts$sample_group), collapse = ", "), "\n")
cat("Cell types identified:", length(unique(cell_type_counts$cell_type)), "\n\n")

cat("Cell Type Distribution Summary:\n")
for (group in unique(cell_type_counts$sample_group)) {
  group_data <- cell_type_counts %>% filter(sample_group == group)
  total_cells <- sum(group_data$count)
  cat("\n", group, "(", total_cells, "cells):\n")
  
  for (ctype in unique(group_data$cell_type)) {
    ctype_count <- sum(group_data$count[group_data$cell_type == ctype])
    percentage <- ctype_count / total_cells * 100
    cat("  ", ctype, ": ", ctype_count, " cells (", round(percentage, 1), "%)\n", sep = "")
  }
}

# Save session information
writeLines(capture.output(sessionInfo()), "stacked_plot_session_info.txt")

cat("\nAnalysis completed successfully! Plots saved to cell_type_composition_plots.pdf\n")