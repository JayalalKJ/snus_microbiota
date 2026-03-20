#!/usr/bin/env Rscript
# =============================================================================


set.seed(123)

suppressPackageStartupMessages({
  library(microeco)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(grid)
  library(plotly)
  library(htmlwidgets)
  library(readr)
})

# ---- 1) I/O -----------------------------------------------------------------
rds_in   <- "results/Results_01_import_mothur_to_microeco/microeco_raw.rds"
out_base <- "results/Results_02B_Negative_control_sample_composition_Plots_MultiRank"

dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(rds_in)) stop("Input RDS not found: ", rds_in)
mt <- readRDS(rds_in)
if (!inherits(mt, "microtable")) stop("Input object is not a microeco::microtable")

# ---- 2) Column resolvers ----------------------------------------------------
pick_col <- function(cands, available) {
  hit <- cands[cands %in% available]
  if (length(hit)) hit[1] else NA_character_
}

st_cols <- colnames(mt$sample_table)
tt_cols <- colnames(mt$tax_table)

sampletype_candidates <- "SampleTypes"
sampletype_col <- pick_col(sampletype_candidates, st_cols)
if (is.na(sampletype_col)) {
  stop(
    "sample_table must contain SampleTypes.\n",
    "Tried: {", paste(sampletype_candidates, collapse = ", "), "}\n",
    "Found: ", paste(st_cols, collapse = ", ")
  )
}
mt$sample_table[[sampletype_col]] <- trimws(as.character(mt$sample_table[[sampletype_col]]))

rank_map <- list(
  Phyla             = c("Phyla", "Phylum"),
  Classes           = c("Classes", "Class"),
  Orders            = c("Orders", "Order"),
  Families          = c("Families", "Family"),
  Genera            = c("Genera", "Genus"),
  Genera_level_OTUs = c("Genera_level_OTUs", "Genus_level_OTUs", "Genera_level_OTU", "Genus_level_OTU")
)
rank_cols <- lapply(rank_map, pick_col, available = tt_cols)

# ---- 3) Plot sets -----------------------------------------------------------
plot_sets <- list(
  All_samples = NULL,
  Blank_only  = "Blank",
  Water_only  = "Water"
)

# ---- 4) Subsetting helper ---------------------------------------------------
subset_microtable_by_sampletypes <- function(m, sampletype_col, keep_types = NULL) {
  m2 <- m$clone(deep = TRUE)
  st <- as.data.frame(m2$sample_table)
  
  if (is.null(keep_types)) {
    m2$tidy_dataset()
    return(m2)
  }
  
  keep_mask <- as.character(st[[sampletype_col]]) %in% keep_types
  st2 <- st[keep_mask, , drop = FALSE]
  if (nrow(st2) == 0) return(NULL)
  
  keep_ids <- rownames(st2)
  m2$sample_table <- st2
  m2$otu_table <- m2$otu_table[, keep_ids, drop = FALSE]
  m2$tidy_dataset()
  m2
}

# ---- 5) Plot builder --------------------------------------------------------
make_barplot <- function(m, taxrank_col, ntaxa, facet_var = NULL) {
  st <- as.data.frame(m$sample_table)
  
  if ("Sample_ID" %in% colnames(st)) {
    rownames(st) <- as.character(st$Sample_ID)
    m$sample_table <- st
    
    if (all(rownames(st) %in% colnames(m$otu_table))) {
      m$otu_table <- m$otu_table[, rownames(st), drop = FALSE]
    }
    m$tidy_dataset()
  }
  
  t1 <- trans_abund$new(dataset = m, taxrank = taxrank_col, ntaxa = ntaxa)
  
  p <- t1$plot_bar(
    others_color = "grey70",
    facet = facet_var,
    xtext_keep = TRUE,
    legend_text_italic = FALSE
  )
  
  p +
    theme_classic(base_size = 13) +
    theme(
      axis.title  = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 6, angle = 65,  vjust = 1, hjust = 1),
      axis.ticks.x = element_line(linewidth = 0.2),
      strip.text  = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text  = element_text(size = 10),
      legend.key.size = unit(0.35, "cm"),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.margin = margin(6, 6, 6, 6)
    ) +
    labs(
      x = "Samples",
      y = "Relative abundance (%)"
    )
}

# ---- 6) Export helper (PNG + PDF + HTML only) -------------------------------
save_all_formats <- function(p, out_dir, prefix) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  w <- 16; h <- 9; dpi <- 450
  
  ggsave(
    filename = file.path(out_dir, paste0(prefix, ".png")),
    plot = p, width = w, height = h, dpi = dpi, bg = "white"
  )
  
  ggsave(
    filename = file.path(out_dir, paste0(prefix, ".pdf")),
    plot = p, width = w, height = h, device = cairo_pdf, bg = "white"
  )
  
  html_file <- file.path(out_dir, paste0(prefix, ".html"))
  widget <- plotly::ggplotly(p)
  
  has_pandoc <- requireNamespace("rmarkdown", quietly = TRUE) &&
    isTRUE(rmarkdown::pandoc_available())
  
  if (has_pandoc) {
    htmlwidgets::saveWidget(widget, file = html_file, selfcontained = TRUE)
  } else {
    libdir <- file.path(out_dir, paste0(prefix, "_files"))
    dir.create(libdir, recursive = TRUE, showWarnings = FALSE)
    htmlwidgets::saveWidget(widget, file = html_file, selfcontained = FALSE, libdir = libdir)
  }
}

# ---- 7) ntaxa per rank ------------------------------------------------------
ntaxa_by_rank <- c(
  Phyla = 8,
  Classes = 12,
  Orders = 15,
  Families = 20,
  Genera = 30,
  Genera_level_OTUs = 50
)

# ---- 8) Main loop -----------------------------------------------------------
available_types <- sort(unique(mt$sample_table[[sampletype_col]]))
message("SampleTypes: ", paste(available_types, collapse = ", "))

readr::write_tsv(
  tibble(
    SampleTypes = names(table(mt$sample_table[[sampletype_col]])),
    n = as.integer(table(mt$sample_table[[sampletype_col]]))
  ) %>% arrange(desc(n)),
  file.path(out_base, "SampleTypes_counts.tsv")
)

for (rank_name in names(rank_cols)) {
  taxrank_col <- rank_cols[[rank_name]]
  
  if (is.na(taxrank_col)) {
    warning("Skipping rank '", rank_name, "': no matching column found in tax_table.")
    next
  }
  
  ntaxa <- unname(ntaxa_by_rank[[rank_name]])
  if (is.na(ntaxa) || !is.finite(ntaxa)) ntaxa <- 20
  
  for (set_name in names(plot_sets)) {
    keep_types <- plot_sets[[set_name]]
    
    mt_sub <- subset_microtable_by_sampletypes(mt, sampletype_col, keep_types = keep_types)
    if (is.null(mt_sub)) {
      warning("No samples for ", set_name, " (rank ", rank_name, "). Skipping.")
      next
    }
    
    facet_var <- if (set_name == "All_samples") sampletype_col else NULL
    
    p <- make_barplot(mt_sub, taxrank_col = taxrank_col, ntaxa = ntaxa, facet_var = facet_var) +
      ggtitle(
        if (set_name == "All_samples") {
          paste0(rank_name, " composition (Top ", ntaxa, ") — All samples (", sampletype_col, ")")
        } else {
          paste0(rank_name, " composition (Top ", ntaxa, ") — ", set_name)
        }
      )
    
    out_dir <- file.path(out_base, rank_name, set_name)
    prefix  <- paste0("barplot_", rank_name, "_top", ntaxa, "_", set_name)
    
    save_all_formats(p, out_dir, prefix)
    message("Saved: ", file.path(out_dir, prefix))
  }
}

message("Done. Output: ", normalizePath(out_base, winslash = "/", mustWork = FALSE))