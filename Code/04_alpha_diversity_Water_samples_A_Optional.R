#!/usr/bin/env Rscript

set.seed(123)

#"Note Only exicute this code if you have more than three water samples from a single tank or a cage

# =============================================================================
# Alpha-diversity analysis of WATER samples using microeco::trans_alpha
#
# Corrected metadata structure:
# - subset SampleType == "Water"
# - group using SampleTypes = F1_Water, F2_Water, F3_Water
# - save water-specific outputs
# - use matching colours and labels
# =============================================================================

# -------------------- Load packages --------------------
suppressPackageStartupMessages({
  library(microeco)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(readr)
  library(tibble)
  library(grid)
})

# -------------------- Input and output --------------------
input_rds  <- "Results/Step_05_A_rarefaction_Skin_for_alpha_diet/microeco_skin_rarefied_20000.rds"
output_dir <- "results/Step_05_C_alpha_diversity_Water"
fig_dir    <- file.path(output_dir, "figures")
tab_dir    <- file.path(output_dir, "tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_rds)) {
  stop("Input file not found: ", input_rds)
}

# -------------------- User-defined metadata settings --------------------
# Change these only if your metadata column names are different
sample_class_col <- "SampleType"   # column that says Fish / Water / Blank
group_col        <- "SampleTypes"  # column that says F1_Water / F2_Water / F3_Water

target_sample_class <- "Water"
target_groups <- c("F1_Water", "F2_Water", "F3_Water")

# -------------------- Read microeco object --------------------
mt <- readRDS(input_rds)

if (!inherits(mt, "microtable")) {
  stop("Input file is not a microeco::microtable object")
}

# -------------------- Check metadata --------------------
cat("Metadata columns:\n")
print(colnames(mt$sample_table))

sample_df <- as.data.frame(mt$sample_table)
colnames(sample_df) <- make.unique(colnames(sample_df))

if (!(sample_class_col %in% colnames(sample_df))) {
  stop("Column not found in sample_table: ", sample_class_col)
}

if (!(group_col %in% colnames(sample_df))) {
  stop("Column not found in sample_table: ", group_col)
}

sample_df[[sample_class_col]] <- trimws(as.character(sample_df[[sample_class_col]]))
sample_df[[group_col]]        <- trimws(as.character(sample_df[[group_col]]))

cat("\nCounts in SampleType column:\n")
print(table(sample_df[[sample_class_col]], useNA = "ifany"))

cat("\nCounts in SampleTypes column:\n")
print(table(sample_df[[group_col]], useNA = "ifany"))

# -------------------- Keep only target water samples --------------------
keep_samples <- rownames(sample_df)[
  sample_df[[sample_class_col]] == target_sample_class &
    sample_df[[group_col]] %in% target_groups
]

if (length(keep_samples) == 0) {
  stop(
    "No samples found after filtering.\n",
    "Check whether:\n",
    "  ", sample_class_col, " == '", target_sample_class, "'\n",
    "and ", group_col, " contains: ",
    paste(target_groups, collapse = ", ")
  )
}

# keep only samples that are also present in otu_table
keep_samples <- intersect(keep_samples, colnames(mt$otu_table))

if (length(keep_samples) == 0) {
  stop("No matching sample IDs found between sample_table and otu_table after filtering.")
}

mt_water <- mt$clone(deep = TRUE)
mt_water$sample_table <- sample_df[keep_samples, , drop = FALSE]
mt_water$otu_table    <- mt_water$otu_table[, keep_samples, drop = FALSE]
mt_water$tidy_dataset()

cat("\nRetained samples after filtering:\n")
print(table(mt_water$sample_table[[group_col]], useNA = "ifany"))

# -------------------- Set group order --------------------
mt_water$sample_table[[group_col]] <- factor(
  mt_water$sample_table[[group_col]],
  levels = target_groups
)

# -------------------- Define x-axis labels --------------------
display_labels <- c(
  "F1_Water" = "F1_Water",
  "F2_Water" = "F2_Water",
  "F3_Water" = "F3_Water"
)

# -------------------- Check phylogenetic tree --------------------
if (is.null(mt_water$phylo_tree)) {
  stop("phylo_tree is missing, so Faith's PD cannot be calculated")
}

# -------------------- Define selected alpha-diversity indices --------------------
alpha_measures <- c("Observed", "Coverage", "Shannon", "PD", "Pielou", "Chao1")

# -------------------- Calculate alpha diversity --------------------
mt_water$cal_alphadiv(
  measures = c("Observed", "Coverage", "Shannon", "Pielou", "Chao1"),
  PD = TRUE
)

cat("\nAlpha diversity columns after calculation:\n")
print(colnames(mt_water$alpha_diversity))

missing_measures <- setdiff(alpha_measures, colnames(mt_water$alpha_diversity))
if (length(missing_measures) > 0) {
  stop("Missing alpha-diversity measures: ", paste(missing_measures, collapse = ", "))
}

# -------------------- Save alpha values --------------------
alpha_table <- mt_water$alpha_diversity %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(
    mt_water$sample_table %>%
      as.data.frame() %>%
      rownames_to_column("SampleID"),
    by = "SampleID"
  )

write_tsv(alpha_table, file.path(tab_dir, "alpha_values_water.tsv"))

# -------------------- Create trans_alpha object --------------------
t1 <- trans_alpha$new(
  dataset = mt_water,
  group   = group_col
)

write_tsv(as.data.frame(t1$data_alpha), file.path(tab_dir, "trans_alpha_data_alpha_water.tsv"))
write_tsv(as.data.frame(t1$data_stat),  file.path(tab_dir, "trans_alpha_data_stat_water.tsv"))

# -------------------- Statistical test: Kruskal-Wallis + Dunn --------------------
t1$cal_diff(
  measure = alpha_measures,
  method = "KW_dunn",
  p_adjust_method = "holm",
  KW_dunn_letter = TRUE
)

write_tsv(as.data.frame(t1$res_diff), file.path(tab_dir, "alpha_stats_KW_Dunn_water.tsv"))

# -------------------- Colors --------------------
group_colors <- c(
  "F1_Water" = "#1b9e77",
  "F2_Water" = "#d95f02",
  "F3_Water" = "#7570b3"
)

# -------------------- Plot each alpha metric --------------------
plot_list <- list()

for (m in alpha_measures) {
  
  p <- t1$plot_alpha(
    measure           = m,
    plot_type         = "ggboxplot",
    add               = "jitter",
    add_sig           = TRUE,
    color_values      = unname(group_colors),
    xtext_angle       = 20,
    xtext_size        = 16,
    ytitle_size       = 18,
    point_size        = 2.8,
    point_alpha       = 0.9,
    add_sig_text_size = 5
  ) +
    scale_x_discrete(labels = display_labels) +
    labs(x = NULL, y = m) +
    theme_classic(base_size = 16) +
    theme(
      text = element_text(face = "bold", colour = "black"),
      axis.title = element_text(face = "bold", size = 18, colour = "black"),
      axis.text = element_text(face = "bold", size = 16, colour = "black"),
      axis.line = element_line(linewidth = 1.2, colour = "black"),
      axis.ticks = element_line(linewidth = 1.2, colour = "black"),
      axis.ticks.length = grid::unit(0.22, "cm"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
      legend.position = "none",
      plot.margin = margin(12, 12, 12, 12)
    )
  
  # make layers thicker
  for (i in seq_along(p$layers)) {
    if ("GeomBoxplot" %in% class(p$layers[[i]]$geom)) {
      p$layers[[i]]$aes_params$linewidth <- 1.1
    }
    if ("GeomPoint" %in% class(p$layers[[i]]$geom)) {
      p$layers[[i]]$aes_params$size <- 2.8
      p$layers[[i]]$aes_params$stroke <- 0.8
    }
    if ("GeomSegment" %in% class(p$layers[[i]]$geom)) {
      p$layers[[i]]$aes_params$linewidth <- 1.1
    }
    if ("GeomErrorbar" %in% class(p$layers[[i]]$geom)) {
      p$layers[[i]]$aes_params$linewidth <- 1.1
    }
  }
  
  plot_list[[m]] <- p
  
  ggsave(
    filename = file.path(fig_dir, paste0("alpha_", m, "_water.png")),
    plot = p,
    width = 9,
    height = 9,
    dpi = 600
  )
  
  ggsave(
    filename = file.path(fig_dir, paste0("alpha_", m, "_water.pdf")),
    plot = p,
    width = 9,
    height = 9
  )
}

# -------------------- Make composite figure --------------------
p_combined <- ggarrange(
  plotlist = plot_list,
  ncol = 3,
  nrow = 2,
  labels = c("A", "B", "C", "D", "E", "F"),
  font.label = list(size = 18, face = "bold", color = "black"),
  common.legend = FALSE
)

ggsave(
  filename = file.path(fig_dir, "alpha_diversity_combined_water.png"),
  plot = p_combined,
  width = 24,
  height = 16,
  dpi = 900
)

ggsave(
  filename = file.path(fig_dir, "alpha_diversity_combined_water.pdf"),
  plot = p_combined,
  width = 24,
  height = 16
)

# -------------------- Done --------------------
message("Alpha diversity analysis for water samples finished")
message("Results saved in: ", output_dir)