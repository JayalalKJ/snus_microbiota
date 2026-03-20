#!/usr/bin/env Rscript

set.seed(123)

# =============================================================================
# Alpha-diversity analysis of skin fish samples using microeco::trans_alpha
# Final version:
# - uses Observed, Coverage, Shannon, PD, Pielou, Chao1
# - avoids duplicate fill/colour scales
# - uses bold text and thicker borders/lines
# - creates a 3 x 2 composite figure
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
input_rds  <- "Results/Results_Step_03_A_rarefaction_Skin_for_alpha_diet/microeco_skin_rarefied_20000.rds"
output_dir <- "results/Results_Step_03_B_alpha_diversity_Skin_simple_Skin_samples"
fig_dir    <- file.path(output_dir, "figures")
tab_dir    <- file.path(output_dir, "tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_rds))

# -------------------- Read microeco object --------------------
mt <- readRDS(input_rds)

if (!inherits(mt, "microtable")) {
  stop("Input file is not a microeco::microtable object")
}

# -------------------- Check metadata --------------------
cat("Metadata columns:\n")
print(colnames(mt$sample_table))

cat("\nSampleType counts:\n")
print(table(mt$sample_table$SampleType, useNA = "ifany"))

cat("\nSampleTypes counts:\n")
print(table(mt$sample_table$SampleTypes, useNA = "ifany"))

# -------------------- Keep only fish samples --------------------
fish_samples <- rownames(mt$sample_table)[mt$sample_table$SampleType == "Fish"]

if (length(fish_samples) == 0) {
  stop("No Fish samples found in SampleType column")
}

mt_fish <- mt$clone(deep = TRUE)
mt_fish$sample_table <- mt_fish$sample_table[fish_samples, , drop = FALSE]
mt_fish$otu_table    <- mt_fish$otu_table[, fish_samples, drop = FALSE]
mt_fish$tidy_dataset()

# -------------------- Set group order --------------------
mt_fish$sample_table$SampleTypes <- factor(
  mt_fish$sample_table$SampleTypes,
  levels = c("Control", "Dose_01%", "Dose_2_2%")
)

# -------------------- Define x-axis labels --------------------
display_labels <- c(
  "Control"   = "F1",
  "Dose_01%"  = "F2",
  "Dose_2_2%" = "F3"
)

# -------------------- Check phylogenetic tree --------------------
if (is.null(mt_fish$phylo_tree)) {
  stop("phylo_tree is missing, so Faith's PD cannot be calculated")
}

# -------------------- Define selected alpha-diversity indices --------------------
alpha_measures <- c("Observed", "Coverage", "Shannon", "PD", "Pielou", "Chao1")

# -------------------- Calculate alpha diversity --------------------
# PD is calculated separately with PD = TRUE
mt_fish$cal_alphadiv(
  measures = c("Observed", "Coverage", "Shannon", "Pielou", "Chao1"),
  PD = TRUE
)

cat("\nAlpha diversity columns after calculation:\n")
print(colnames(mt_fish$alpha_diversity))

# Safety check
missing_measures <- setdiff(alpha_measures, colnames(mt_fish$alpha_diversity))
if (length(missing_measures) > 0) {
  stop("These selected alpha measures are missing after cal_alphadiv(): ",
       paste(missing_measures, collapse = ", "))
}

# -------------------- Save alpha values --------------------
alpha_table <- mt_fish$alpha_diversity %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(
    mt_fish$sample_table %>%
      as.data.frame() %>%
      rownames_to_column("SampleID"),
    by = "SampleID"
  )

write_tsv(alpha_table, file.path(tab_dir, "alpha_values_fish.tsv"))

# -------------------- Create trans_alpha object --------------------
t1 <- trans_alpha$new(
  dataset = mt_fish,
  group   = "SampleTypes"
)

write_tsv(as.data.frame(t1$data_alpha), file.path(tab_dir, "trans_alpha_data_alpha.tsv"))
write_tsv(as.data.frame(t1$data_stat),  file.path(tab_dir, "trans_alpha_data_stat.tsv"))

# -------------------- Statistical test: Kruskal-Wallis + Dunn --------------------
t1$cal_diff(
  measure = alpha_measures,
  method  = "KW_dunn",
  p_adjust_method = "holm",
  KW_dunn_letter  = TRUE
)

write_tsv(as.data.frame(t1$res_diff), file.path(tab_dir, "alpha_stats_KW_Dunn.tsv"))

# -------------------- Colors --------------------
group_colors <- c(
  "Control"   = "#1b9e77",
  "Dose_01%"  = "#d95f02",
  "Dose_2_2%" = "#7570b3"
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
  
  # Thicken existing layers created by plot_alpha()
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
    filename = file.path(fig_dir, paste0("alpha_", m, ".png")),
    plot = p,
    width = 9,
    height = 9,
    dpi = 600
  )
  
  ggsave(
    filename = file.path(fig_dir, paste0("alpha_", m, ".pdf")),
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
  filename = file.path(fig_dir, "alpha_diversity_combined.png"),
  plot = p_combined,
  width = 24,
  height = 16,
  dpi = 900
)

ggsave(
  filename = file.path(fig_dir, "alpha_diversity_combined.pdf"),
  plot = p_combined,
  width = 24,
  height = 16
)

# -------------------- Done --------------------
message("Alpha diversity analysis finished")
message("Results saved in: ", output_dir)