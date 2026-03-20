#!/usr/bin/env Rscript

set.seed(123)

# =============================================================================
# SIMPLE beta-diversity analysis using microeco::trans_beta
# Main analysis:
#   1. Keep only the groups you want
#   2. Build Bray-Curtis distance
#   3. PCoA plot
#   4. PERMANOVA
#   5. PERMDISP
#   6. Within-group Bray-Curtis distances
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Load packages
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(microeco)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(readr)
  library(tibble)
  library(vegan)
})

# -----------------------------------------------------------------------------
# 2. Helper function to save plots
# -----------------------------------------------------------------------------
save_plot <- function(p, filename_base, width = 7, height = 6, dpi = 600) {
  ggsave(
    filename = paste0(filename_base, ".png"),
    plot = p,
    width = width,
    height = height,
    dpi = dpi
  )
  ggsave(
    filename = paste0(filename_base, ".pdf"),
    plot = p,
    width = width,
    height = height
  )
}

# -----------------------------------------------------------------------------
# 3. Input and output
# Change here if your file path is different
# -----------------------------------------------------------------------------
#input_rds  <- "Results/Results_02C_contaminant_removal_from_blanks/microeco_cleaned_after_blank_contaminants.rds"
input_rds  <- "Results/Results_Step_03_A_rarefaction_Skin_for_alpha_diet/microeco_skin_rarefied_20000.rds"
output_dir <- "results/Results_Step_06_beta_diversity_Skin_samples"
fig_dir    <- file.path(output_dir, "figures")
tab_dir    <- file.path(output_dir, "tables")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_rds)) {
  stop("Input file not found: ", input_rds)
}

# -----------------------------------------------------------------------------
# 4. Read microeco object
# -----------------------------------------------------------------------------
mt <- readRDS(input_rds)

if (!inherits(mt, "microtable")) {
  stop("Input file is not a microeco::microtable object.")
}

# -----------------------------------------------------------------------------
# 5. Check metadata
# We use only one metadata column: SampleTypes
# -----------------------------------------------------------------------------
mt$sample_table <- as.data.frame(mt$sample_table)
colnames(mt$sample_table) <- make.unique(colnames(mt$sample_table))

if (!"SampleTypes" %in% colnames(mt$sample_table)) {
  stop("Column 'SampleTypes' not found in sample_table.")
}

mt$sample_table$SampleTypes <- trimws(as.character(mt$sample_table$SampleTypes))

cat("Metadata columns:\n")
print(colnames(mt$sample_table))

cat("\nSampleTypes counts before filtering:\n")
print(table(mt$sample_table$SampleTypes, useNA = "ifany"))

# -----------------------------------------------------------------------------
# 6. Keep only the three target groups
# Change only these group names if needed
# -----------------------------------------------------------------------------
target_groups <- c("Control", "Dose_01%", "Dose_2_2%")

keep_samples <- rownames(mt$sample_table)[mt$sample_table$SampleTypes %in% target_groups]

if (length(keep_samples) == 0) {
  stop("No samples found for target groups: ", paste(target_groups, collapse = ", "))
}

# Keep only samples that are also in OTU table
keep_samples <- intersect(keep_samples, colnames(mt$otu_table))

if (length(keep_samples) == 0) {
  stop("No matching sample IDs found between sample_table and otu_table.")
}

mt_sub <- mt$clone(deep = TRUE)
mt_sub$sample_table <- mt_sub$sample_table[keep_samples, , drop = FALSE]
mt_sub$otu_table    <- mt_sub$otu_table[, keep_samples, drop = FALSE]
mt_sub$tidy_dataset()

cat("\nSampleTypes counts after filtering:\n")
print(table(mt_sub$sample_table$SampleTypes, useNA = "ifany"))

# -----------------------------------------------------------------------------
# 7. Set group order
# This controls plot order and legend order
# -----------------------------------------------------------------------------
mt_sub$sample_table$SampleTypes <- factor(
  mt_sub$sample_table$SampleTypes,
  levels = target_groups
)
mt_sub$sample_table$SampleTypes <- droplevels(mt_sub$sample_table$SampleTypes)

# -----------------------------------------------------------------------------
# 8. Define labels and colors
# Change here if you want different display names or colors
# -----------------------------------------------------------------------------
display_labels <- c(
  "Control"   = "F1",
  "Dose_01%"  = "F2",
  "Dose_2_2%" = "F3"
)

group_colors <- c(
  "Control"   = "#1b9e77",
  "Dose_01%"  = "#d95f02",
  "Dose_2_2%" = "#7570b3"
)

group_order <- c("Control", "Dose_01%", "Dose_2_2%")

# -----------------------------------------------------------------------------
# 9. Build Bray-Curtis distance matrix
# trans_beta uses distance matrices stored in mt_sub$beta_diversity
# -----------------------------------------------------------------------------
otu_mat <- as.matrix(mt_sub$otu_table)

# Make otu table columns match metadata row order
otu_mat <- otu_mat[, rownames(mt_sub$sample_table), drop = FALSE]

# Convert to samples x taxa
sample_by_taxa <- t(otu_mat)

if (is.null(mt_sub$beta_diversity)) {
  mt_sub$beta_diversity <- list()
}

mt_sub$beta_diversity$bray <- as.matrix(
  vegdist(sample_by_taxa, method = "bray")
)

# Save Bray matrix
write_tsv(
  as.data.frame(mt_sub$beta_diversity$bray) %>% rownames_to_column("SampleID"),
  file.path(tab_dir, "beta_bray_matrix.tsv")
)

# -----------------------------------------------------------------------------
# 10. Create trans_beta object
# This is the main beta-diversity object
# -----------------------------------------------------------------------------
tb <- trans_beta$new(
  dataset = mt_sub,
  group   = "SampleTypes",
  measure = "bray"
)

# -----------------------------------------------------------------------------
# 11. PCoA ordination
# -----------------------------------------------------------------------------
tb$cal_ordination(method = "PCoA")

# Save ordination scores
write_tsv(
  as.data.frame(tb$res_ordination$scores),
  file.path(tab_dir, "beta_bray_pcoa_scores.tsv")
)

# -----------------------------------------------------------------------------
# 12. Plot PCoA
# -----------------------------------------------------------------------------
p_pcoa <- tb$plot_ordination(
  plot_color = "SampleTypes",
  plot_type = c("point", "ellipse"),
  plot_group_order = group_order,
  point_size = 3.5,
  point_alpha = 0.85,
  ellipse_level = 0.95,
  ellipse_chull_alpha = 0.15
) +
  scale_color_manual(values = group_colors, labels = display_labels, drop = FALSE) +
  scale_fill_manual(values = group_colors, labels = display_labels, drop = FALSE) +
  labs(
    x = "PCoA1",
    y = "PCoA2",
    color = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

save_plot(
  p_pcoa,
  file.path(fig_dir, "beta_bray_pcoa"),
  width = 7,
  height = 6
)

# -----------------------------------------------------------------------------
# 13. PERMANOVA
# This tests whether overall community composition differs among groups
# -----------------------------------------------------------------------------
tb$cal_manova(
  manova_all = TRUE,
  permutations = 999
)

cat("\nPERMANOVA result:\n")
print(tb$res_manova)

write_tsv(
  as.data.frame(tb$res_manova),
  file.path(tab_dir, "beta_bray_permanova.tsv")
)

# -----------------------------------------------------------------------------
# 14. PERMDISP
# This checks whether group spread is different
# -----------------------------------------------------------------------------
tb$cal_betadisper()

cat("\nPERMDISP result:\n")
print(tb$res_betadisper)

capture.output(
  print(tb$res_betadisper),
  file = file.path(tab_dir, "beta_bray_betadisper.txt")
)

# -----------------------------------------------------------------------------
# 15. Within-group Bray-Curtis distances
# This shows how variable samples are inside each group
# -----------------------------------------------------------------------------
tb$cal_group_distance(within_group = TRUE)

write_tsv(
  as.data.frame(tb$res_group_distance),
  file.path(tab_dir, "beta_bray_within_group_distances.tsv")
)

# -----------------------------------------------------------------------------
# 16. Statistics for within-group distances
# -----------------------------------------------------------------------------
tb$cal_group_distance_diff(
  method = "KW_dunn",
  p_adjust_method = "holm",
  KW_dunn_letter = TRUE
)

write_tsv(
  as.data.frame(tb$res_group_distance_diff),
  file.path(tab_dir, "beta_bray_within_group_distance_stats.tsv")
)

# -----------------------------------------------------------------------------
# 17. Plot within-group distances
# -----------------------------------------------------------------------------
p_dist <- tb$plot_group_distance(
  plot_group_order = group_order,
  plot_type = "ggboxplot",
  add = "jitter",
  xtext_angle = 20,
  xtext_size = 12,
  ytitle_size = 14
) +
  scale_fill_manual(values = group_colors, labels = display_labels, drop = FALSE) +
  scale_color_manual(values = group_colors, labels = display_labels, drop = FALSE) +
  scale_x_discrete(labels = display_labels, drop = FALSE) +
  labs(
    x = NULL,
    y = "Within-group Bray-Curtis distance"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none"
  )

save_plot(
  p_dist,
  file.path(fig_dir, "beta_bray_within_group_distance"),
  width = 7,
  height = 5
)

# -----------------------------------------------------------------------------
# 18. Save the sample metadata used
# -----------------------------------------------------------------------------
write_tsv(
  mt_sub$sample_table %>%
    as.data.frame() %>%
    rownames_to_column("SampleID"),
  file.path(tab_dir, "beta_metadata_used.tsv")
)

# -----------------------------------------------------------------------------
# 19. Done
# -----------------------------------------------------------------------------
message("Simple beta-diversity analysis finished")
message("Results saved in: ", normalizePath(output_dir))