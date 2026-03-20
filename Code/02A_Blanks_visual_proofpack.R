#!/usr/bin/env Rscript

## 0. Packages ----
suppressPackageStartupMessages({
  library(microeco)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)      # for unit()
  library(cowplot)
  library(stringr)
  library(tools)
  library(viridis)   
})

## 0a. Heatmap helpers (palette + tiles) ----
make_heat_scale <- function(key = "current") {
  if (key == "viridis") {
    viridis::scale_fill_viridis(trans = "log", name = "% Relative Abundance", na.value = "grey80")
  } else if (key == "magma") {
    viridis::scale_fill_viridis(option = "A", trans = "log", name = "% Relative Abundance", na.value = "grey80")
  } else {
    scale_fill_gradientn(
      colors = c("#A0A0A0", "#4B4B4B", "#000000", "#FF0000"),
      trans  = "log",
      name   = "% Relative Abundance",
      na.value = "grey80"
    )
  }
}

widen_tiles <- function(p, width = 1.2) {
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomTile")) {
      p$layers[[i]]$aes_params$width <- width
    }
  }
  p
}

## 1. Load microeco object ----
dataset_path <- "Results/Results_01_import_mothur_to_microeco/microeco_raw.rds"
if (!file.exists(dataset_path)) stop("❌ Input RDS not found: ", dataset_path)

mt <- readRDS(dataset_path)
if (!inherits(mt, "microtable")) stop("❌ Loaded object is not a microeco::microtable")
message("✅ Loaded microeco object: ", dataset_path)

## 2. Rank resolver (plural/singular smart pick) ------------------------------
tt <- mt$tax_table
avail_tax_cols <- colnames(tt)

rank_groups <- list(
  Phylum  = c("Phylum", "Phyla"),
  Class   = c("Class", "Classes"),
  Order   = c("Order", "Orders"),
  Family  = c("Family", "Families"),
  Genus   = c("Genus_level_OTUs", "Genera_level_OTUs")
)

pick_rank_label <- function(group_key, available_cols) {
  opts <- rank_groups[[group_key]]
  hit  <- opts[opts %in% available_cols]
  if (length(hit)) return(hit[1])
  NA_character_
}

ranks_to_use <- vapply(
  names(rank_groups),
  pick_rank_label,
  FUN.VALUE = character(1),
  available_cols = avail_tax_cols
)
ranks_to_use <- ranks_to_use[!is.na(ranks_to_use)]
message("• Will plot ranks: ", paste(ranks_to_use, collapse = ", "))

## Facets (only those present will be used) ----
#facet_cols <- c("Tank", "SampleTypes", "Labs")
facet_cols <- c("Tanks", "SampleTypes")
facet_cols <- facet_cols[facet_cols %in% colnames(mt$sample_table)]
if (!length(facet_cols)) facet_cols <- NULL

## 3. Config + single output folder + saver ==================================
bar_ntaxa   <- 25
heat_ntaxa  <- 100
fig_w       <- 28
fig_h       <- 16
tile_width  <- 1.2
palette_key <- "current" # "current" | "viridis" | "magma"

out_dir <- file.path("results", "Resuls_02A_Blanks_Visual_proofpack", "Full_dataset")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_png_pdf <- function(file_base, plot, width, height, dpi = 600) {
  png_file <- paste0(file_base, ".png")
  pdf_file <- paste0(file_base, ".pdf")
  
  ggsave(filename = png_file, plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(filename = pdf_file, plot = plot, width = width, height = height, device = cairo_pdf)
  
  message("✓ Saved: ", png_file)
  message("✓ Saved: ", pdf_file)
}

## 4. Export sample table (info) ---------------------------------------------
try({
  write.csv(
    mt$sample_table,
    file = file.path(out_dir, "FULL_sample_table.csv"),
    row.names = TRUE
  )
}, silent = TRUE)

## 5. FULL barplots (Top N) ---------------------------------------------------
# groupmean expects a column name in sample_table (your script uses SampleTypes)
groupmean_col <- if ("SampleTypes" %in% colnames(mt$sample_table)) "SampleTypes" else NULL
if (is.null(groupmean_col)) {
  warning("⚠️ 'SampleTypes' column not found in sample_table. Barplots will still run, but groupmean may fail.")
}

for (tax_rank in ranks_to_use) {
  
  ntop <- min(bar_ntaxa, nrow(mt$tax_table))  # guard
  
  ta <- trans_abund$new(
    dataset   = mt,
    taxrank   = tax_rank,
    ntaxa     = ntop,
    groupmean = groupmean_col
  )
  
  gp <- ta$plot_bar(others_color = "grey70", xtext_angle = 30) +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text.x  = element_text(size = 10, angle = 65, hjust = 1),
      axis.text.y  = element_text(size = 10),
      strip.text   = element_text(size = 11, face = "bold"),
      legend.position  = "right",
      legend.key.size  = unit(0.5, "cm"),
      legend.title     = element_text(size = 9, face = "bold"),
      legend.text      = element_text(size = 8),
      plot.title       = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.margin      = margin(10, 40, 10, 10)
    ) +
    guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
    labs(
      fill  = tax_rank,
      title = paste0("Top ", ntop, " taxa at ", tax_rank, " (group = ", ifelse(is.null(groupmean_col), "N/A", groupmean_col), ") — FULL"),
      x     = "Sample Types",
      y     = "Relative Abundance (%)"
    )
  
  save_png_pdf(
    file_base = file.path(out_dir, paste0("FULL_group_mean_", tax_rank, "_barplot")),
    plot      = gp,
    width     = 12,
    height    = 7
  )
  
  if (!is.null(ta$data_abund)) {
    out_csv <- file.path(out_dir, paste0("FULL_bar_", tax_rank, "_data.csv"))
    try(write.csv(ta$data_abund, file = out_csv, row.names = FALSE), silent = TRUE)
  }
}

## 6. FULL heatmap (Top 100 Genus_level_OTUs) --------------------------------
genus_rank_for_heat <- pick_rank_label("Genus", colnames(mt$tax_table))

if (is.na(genus_rank_for_heat)) {
  warning("⚠️ No Genus_level_OTUs/Genera_level_OTUs column found; skipping heatmap.")
} else {
  
  ntop_hm <- min(heat_ntaxa, sum(!is.na(mt$tax_table[[genus_rank_for_heat]])))
  heat_data <- trans_abund$new(dataset = mt, taxrank = genus_rank_for_heat, ntaxa = ntop_hm)
  
  breaks_vec <- c(0.001, 0.01, 0.1, 1, 10)
  
  hm <- try(
    heat_data$plot_heatmap(
      facet = facet_cols,
      xtext_keep = TRUE,
      withmargin = TRUE,
      plot_breaks = breaks_vec
    ),
    silent = TRUE
  )
  
  if (inherits(hm, "try-error")) {
    hm <- heat_data$plot_heatmap(
      facet = facet_cols,
      xtext_keep = TRUE,
      withmargin = TRUE
    )
  }
  
  hm <- hm +
    make_heat_scale(palette_key) +
    scale_x_discrete(expand = c(0, 0)) +
    theme(panel.spacing.x = unit(0.2, "lines"))
  
  hm <- widen_tiles(hm, tile_width)
  
  hm <- hm +
    theme_classic() +
    theme(
      axis.title   = element_text(size = 14, face = "bold"),
      axis.text.x  = element_text(size = 10, angle = 35, hjust = 1),
      axis.text.y  = element_text(size = 10, face = "italic"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text  = element_text(size = 10),
      plot.title   = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.margin  = margin(10, 40, 10, 10)
    ) +
    labs(
      title = paste0("Heatmap of Top ", ntop_hm, " ", genus_rank_for_heat, " — FULL"),
      x = "Samples",
      y = genus_rank_for_heat
    )
  
  save_png_pdf(
    file_base = file.path(out_dir, paste0("FULL_heatmap_", genus_rank_for_heat, "_top", ntop_hm)),
    plot      = hm,
    width     = fig_w,
    height    = fig_h
  )
  
  if (!is.null(heat_data$data_abund)) {
    try(
      write.csv(
        heat_data$data_abund,
        file = file.path(out_dir, "FULL_heatmap_data.csv"),
        row.names = FALSE
      ),
      silent = TRUE
    )
  }
}

message("FULL dataset plots + info saved under: ", normalizePath(out_dir))
