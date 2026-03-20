#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

## If needed, install first:
## install.packages(c("paletteer", "ggsci", "ggthemes", "RColorBrewer"))

## 0 │ Libraries --------------------------------------------------------------
suppressPackageStartupMessages({
  library(microeco)
  library(dplyr)
  library(ggplot2)
  library(ggh4x)
  library(readr)
  library(tibble)
  library(stringr)
  library(grid)
  library(tools)
  library(paletteer)
})

## 0.1 │ Logging --------------------------------------------------------------
stamp <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = "")))
}

## 1 │ I/O -------------------------------------------------------------------
dataset_path <- "Results/Results_02C_contaminant_removal_from_blanks/microeco_cleaned_after_blank_contaminants.rds"
output_root  <- "results/Step_11_Taxonomic_composition_Skin_bySampleTypes"
output_dir   <- output_root

## 2 │ Config ----------------------------------------------------------------
# Remove these tanks from the whole dataset
drop_tanks <- c("SC")

# Remove these SampleTypes only from the group-mean dataset
drop_sampletypes_groupmean <- c("Water")

# Ranks to plot
rank_cfg <- list(
  Phyla             = list(taxrank = "Phyla",             ntaxa = 8),
  Classes           = list(taxrank = "Classes",           ntaxa = 10),
  Orders            = list(taxrank = "Orders",            ntaxa = 12),
  Families          = list(taxrank = "Families",          ntaxa = 15),
  Genera            = list(taxrank = "Genera",            ntaxa = 25),
  Genera_level_OTUs = list(taxrank = "Genera_level_OTUs", ntaxa = 35)
)

# Different palette for each rank
rank_palettes <- list(
  Phyla             = "RColorBrewer::Set3",
  Classes           = "RColorBrewer::Paired",
  Orders            = "ggsci::category20_d3",
  Families          = "ggthemes::Classic_20",
  Genera            = "ggsci::category20c_d3",
  Genera_level_OTUs = "ggthemes::Classic_20"
)

# Facets for sample-level plot
facet_vars <- c("Tanks", "SampleTypes")

# Group-mean settings
MAKE_GROUPMEAN_PLOTS <- TRUE
groupmean_var        <- "Tank_SampleType"   # choose: "Tanks", "SampleTypes", or "Tank_SampleType"

# Default plot controls
BAR_WIDTH         <- 17
BAR_HEIGHT        <- 7
GROUPMEAN_WIDTH   <- 12
GROUPMEAN_HEIGHT  <- 8
DPI               <- 600
X_TEXT_ANGLE      <- 65
GROUP_X_ANGLE     <- 35
LEGEND_NCOL       <- 1
OTHERS_COLOR      <- "grey70"
BARWIDTH          <- 0.92
BASE_FAMILY       <- ""

QC_DIR <- file.path(output_root, "QC")

## 2.1 │ Special group-mean settings for ALL ranks ----------------------------
special_groupmean_cfg <- list(
  Phyla = list(
    width             = 12,
    height            = 8,
    barwidth          = 0.90,
    x_angle           = 35,
    x_text_size       = 11,
    legend_ncol       = 1,
    legend_text_size  = 14,
    legend_title_size = 15,
    legend_key_size   = 0.35,
    wrap_x_labels     = TRUE,
    wrap_width        = 56,
    custom_title      = "Phylum-level composition across sample groups (group mean)"
  ),
  Classes = list(
    width             = 12,
    height            = 8,
    barwidth          = 0.90,
    x_angle           = 35,
    x_text_size       = 11,
    legend_ncol       = 1,
    legend_text_size  = 13,
    legend_title_size = 15,
    legend_key_size   = 0.35,
    wrap_x_labels     = TRUE,
    wrap_width        = 56,
    custom_title      = "Class-level composition across sample groups (group mean)"
  ),
  Orders = list(
    width             = 12,
    height            = 8,
    barwidth          = 0.90,
    x_angle           = 40,
    x_text_size       = 10,
    legend_ncol       = 1,
    legend_text_size  = 12,
    legend_title_size = 15,
    legend_key_size   = 0.35,
    wrap_x_labels     = TRUE,
    wrap_width        = 56,
    custom_title      = "Order-level composition across sample groups (group mean)"
  ),
  Families = list(
    width             = 12,
    height            = 8,
    barwidth          = 0.90,
    x_angle           = 45,
    x_text_size       = 10,
    legend_ncol       = 1,
    legend_text_size  = 12,
    legend_title_size = 15,
    legend_key_size   = 0.35,
    wrap_x_labels     = TRUE,
    wrap_width        = 56,
    custom_title      = "Family-level composition across sample groups (group mean)"
  ),
  Genera = list(
    width             = 12,
    height            = 8,
    barwidth          = 0.90,
    x_angle           = 50,
    x_text_size       = 9,
    legend_ncol       = 1,
    legend_text_size  = 12,
    legend_title_size = 15,
    legend_key_size   = 0.35,
    wrap_x_labels     = TRUE,
    wrap_width        = 56,
    custom_title      = "Genus-level composition across sample groups (group mean)"
  ),
  Genera_level_OTUs = list(
    width             = 12,
    height            = 8,
    barwidth          = 0.90,
    x_angle           = 55,
    x_text_size       = 9,
    legend_ncol       = 1,
    legend_text_size  = 14,
    legend_title_size = 15,
    legend_key_size   = 0.35,
    wrap_x_labels     = TRUE,
    wrap_width        = 56,
    custom_title      = "Genera-level OTU composition across sample groups (group mean)"
  )
)

## 3 │ Helpers ----------------------------------------------------------------

load_microeco_object <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  ext <- tolower(file_ext(path))
  
  if (ext == "rds") {
    obj <- readRDS(path)
    if (!inherits(obj, "microtable")) {
      stop("RDS does not contain a microtable. Got: ", paste(class(obj), collapse = ", "))
    }
    return(obj)
  }
  
  if (ext %in% c("rdata", "rda")) {
    env <- new.env(parent = emptyenv())
    objs <- load(path, envir = env)
    for (nm in objs) {
      if (inherits(env[[nm]], "microtable")) return(env[[nm]])
    }
    stop("No microtable found in ", path, ". Objects: ", paste(objs, collapse = ", "))
  }
  
  stop("Unsupported extension: ", ext)
}

save_plot_png <- function(out_dir, name, plot, w, h, dpi = 300, bg = "white") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  png_path <- file.path(out_dir, paste0(name, ".png"))
  ggsave(
    filename  = png_path,
    plot      = plot,
    width     = w,
    height    = h,
    dpi       = dpi,
    limitsize = FALSE,
    bg        = bg
  )
  stamp("Saved: ", png_path)
}

assert_cols <- function(df, cols, where = "sample_table") {
  miss <- setdiff(cols, colnames(df))
  if (length(miss)) {
    stop(
      "Missing columns in ", where, ": ", paste(miss, collapse = ", "),
      "\nAvailable: ", paste(colnames(df), collapse = ", ")
    )
  }
  invisible(TRUE)
}

get_special_groupmean_cfg <- function(rank_name, default_width, default_height,
                                      default_barwidth, default_angle,
                                      default_legend_ncol) {
  cfg <- list(
    width             = default_width,
    height            = default_height,
    barwidth          = default_barwidth,
    x_angle           = default_angle,
    x_text_size       = 10,
    legend_ncol       = default_legend_ncol,
    legend_text_size  = 13,
    legend_title_size = 15,
    legend_key_size   = 0.40,
    wrap_x_labels     = FALSE,
    wrap_width        = 18,
    custom_title      = NULL
  )
  
  if (rank_name %in% names(special_groupmean_cfg)) {
    this_cfg <- special_groupmean_cfg[[rank_name]]
    for (nm in names(this_cfg)) cfg[[nm]] <- this_cfg[[nm]]
  }
  
  cfg
}

get_rank_colors <- function(rank_name, n_taxa) {
  pal_name <- rank_palettes[[rank_name]]
  
  if (is.null(pal_name)) {
    pal_name <- "RColorBrewer::Dark2"
  }
  
  base_cols <- as.character(paletteer::paletteer_d(pal_name))
  
  if (length(base_cols) >= n_taxa) {
    return(base_cols[seq_len(n_taxa)])
  }
  
  grDevices::colorRampPalette(base_cols)(n_taxa)
}

## 4 │ Load data --------------------------------------------------------------
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

stamp("Loading: ", dataset_path)
mt <- load_microeco_object(dataset_path)

if (!inherits(mt, "microtable")) {
  stop("Input is not a microeco::microtable")
}

stamp("INPUT COUNTS:")
stamp("  • Total OTUs   : ", nrow(mt$otu_table))
stamp("  • Total samples: ", ncol(mt$otu_table))

st0 <- as.data.frame(mt$sample_table)
stamp("SAMPLE TABLE: ", nrow(st0), " rows × ", ncol(st0), " cols")

## 5 │ Filter tanks + preserve sample IDs safely ------------------------------
st <- as.data.frame(mt$sample_table) %>%
  rownames_to_column(var = "Sample_ID")

assert_cols(st, facet_vars, where = "sample_table")

st <- st %>%
  mutate(
    across(where(is.character), ~ trimws(.x)),
    Tanks       = trimws(Tanks),
    SampleTypes = trimws(SampleTypes)
  ) %>%
  filter(!Tanks %in% drop_tanks) %>%
  mutate(
    Tank_SampleType = paste(Tanks, SampleTypes, sep = " | ")
  )

assert_cols(
  st,
  c("Sample_ID", "Tanks", "SampleTypes", "Tank_SampleType"),
  where = "filtered sample_table"
)

st <- st %>%
  column_to_rownames(var = "Sample_ID")

st <- as.data.frame(st, stringsAsFactors = FALSE)
st$Tanks <- droplevels(as.factor(st$Tanks))
st$SampleTypes <- droplevels(as.factor(st$SampleTypes))
st$Tank_SampleType <- factor(st$Tank_SampleType, levels = unique(st$Tank_SampleType))

mt$sample_table <- st

keep_samp <- rownames(mt$sample_table)
missing_in_otu <- setdiff(keep_samp, colnames(mt$otu_table))
if (length(missing_in_otu)) {
  stop("These filtered samples are missing from otu_table: ", paste(missing_in_otu, collapse = ", "))
}

mt$otu_table <- mt$otu_table[, keep_samp, drop = FALSE]

if (ncol(mt$otu_table) == 0) {
  stop("No samples left after filtering Tanks. Check drop_tanks names.")
}

stopifnot(identical(rownames(mt$sample_table), colnames(mt$otu_table)))
mt$tidy_dataset()

stamp("AFTER dropping unwanted tanks:")
stamp("  • OTUs   : ", nrow(mt$otu_table))
stamp("  • Samples: ", ncol(mt$otu_table))

write_tsv(
  tibble(Sample_ID = rownames(mt$sample_table)),
  file.path(QC_DIR, "samples_used_after_tank_filter.tsv")
)

write_tsv(
  mt$sample_table %>%
    rownames_to_column(var = "Sample_ID") %>%
    as_tibble(),
  file.path(QC_DIR, "sample_table_after_tank_filter.tsv")
)

## 5c │ Validate groupmean column ---------------------------------------------
if (MAKE_GROUPMEAN_PLOTS) {
  if (!groupmean_var %in% colnames(mt$sample_table)) {
    stop(
      "groupmean_var '", groupmean_var, "' not found in sample_table.\n",
      "Available columns: ", paste(colnames(mt$sample_table), collapse = ", ")
    )
  }
}

## 5d │ TSS normalization + save ----------------------------------------------
tn  <- trans_norm$new(dataset = mt)
TSS <- tn$norm(method = "TSS")

save(
  TSS,
  file = file.path(output_dir, "Composition_TSS.RData"),
  compress = TRUE
)

stamp("Saved: ", file.path(output_dir, "Composition_TSS.RData"))

mt_plot <- TSS
mt_plot$tidy_dataset()
mt_plot$cal_abund()

## 5e │ Prepare separate dataset for group-mean plots -------------------------
mt_groupmean <- NULL

if (MAKE_GROUPMEAN_PLOTS) {
  stamp("Preparing group-mean dataset ...")
  
  mt_groupmean <- mt_plot$clone(deep = TRUE)
  
  st_gm <- as.data.frame(mt_groupmean$sample_table) %>%
    rownames_to_column(var = "Sample_ID") %>%
    mutate(
      Tanks       = trimws(Tanks),
      SampleTypes = trimws(SampleTypes)
    ) %>%
    filter(!SampleTypes %in% drop_sampletypes_groupmean) %>%
    mutate(
      Tank_SampleType = paste(Tanks, SampleTypes, sep = " | ")
    ) %>%
    column_to_rownames(var = "Sample_ID")
  
  st_gm <- as.data.frame(st_gm, stringsAsFactors = FALSE)
  
  if (nrow(st_gm) == 0) {
    stop("No samples left after filtering SampleTypes for group-mean plots.")
  }
  
  st_gm$Tanks <- droplevels(as.factor(st_gm$Tanks))
  st_gm$SampleTypes <- droplevels(as.factor(st_gm$SampleTypes))
  st_gm$Tank_SampleType <- factor(st_gm$Tank_SampleType, levels = unique(st_gm$Tank_SampleType))
  
  keep_gm <- rownames(st_gm)
  missing_gm <- setdiff(keep_gm, colnames(mt_groupmean$otu_table))
  if (length(missing_gm)) {
    stop("These group-mean samples are missing from otu_table: ", paste(missing_gm, collapse = ", "))
  }
  
  mt_groupmean$sample_table <- st_gm
  mt_groupmean$otu_table    <- mt_groupmean$otu_table[, keep_gm, drop = FALSE]
  
  stopifnot(identical(rownames(mt_groupmean$sample_table), colnames(mt_groupmean$otu_table)))
  
  if (!groupmean_var %in% colnames(mt_groupmean$sample_table)) {
    stop(
      "groupmean_var '", groupmean_var, "' not found in group-mean sample_table.\n",
      "Available columns: ", paste(colnames(mt_groupmean$sample_table), collapse = ", ")
    )
  }
  
  mt_groupmean$tidy_dataset()
  mt_groupmean$cal_abund()
  
  write_tsv(
    mt_groupmean$sample_table %>%
      rownames_to_column(var = "Sample_ID") %>%
      as_tibble(),
    file.path(QC_DIR, "sample_table_groupmean_after_sampletype_filter.tsv")
  )
  
  stamp(
    "Group-mean dataset prepared after dropping SampleTypes: ",
    paste(drop_sampletypes_groupmean, collapse = ", ")
  )
}

## 6 │ Plot all ranks ---------------------------------------------------------
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

if (MAKE_GROUPMEAN_PLOTS && is.null(mt_groupmean)) {
  stop("Group-mean dataset was not created properly.")
}

for (nm in names(rank_cfg)) {
  
  taxrank <- rank_cfg[[nm]]$taxrank
  ntaxa   <- rank_cfg[[nm]]$ntaxa
  
  if (!taxrank %in% colnames(mt_plot$tax_table)) {
    stop(
      "taxrank not found in tax_table: ", taxrank,
      "\nAvailable: ", paste(colnames(mt_plot$tax_table), collapse = ", ")
    )
  }
  
  rank_colors <- get_rank_colors(nm, ntaxa)
  
  stamp("Plotting: ", nm, " | taxrank=", taxrank, " | top=", ntaxa)
  
  ## ---- Sample-level barplot ------------------------------------------------
  t1 <- trans_abund$new(
    dataset = mt_plot,
    taxrank = taxrank,
    ntaxa   = ntaxa
  )
  
  p <- t1$plot_bar(
    others_color       = OTHERS_COLOR,
    color_values       = rank_colors,
    facet              = facet_vars,
    xtext_keep         = FALSE,
    legend_text_italic = FALSE,
    barwidth           = BARWIDTH
  ) +
    theme_classic(base_family = BASE_FAMILY) +
    theme(
      strip.text       = element_text(size = 12, face = "bold"),
      axis.title.x     = element_text(size = 15, face = "bold"),
      axis.title.y     = element_text(size = 15, face = "bold"),
      axis.text.y      = element_text(size = 10),
      axis.text.x      = element_text(angle = X_TEXT_ANGLE, hjust = 1, vjust = 1, size = 9),
      legend.position  = "right",
      legend.title     = element_text(size = 15, face = "bold"),
      legend.text      = element_text(size = 13),
      legend.key.size  = unit(0.40, "cm"),
      legend.spacing   = unit(0.05, "cm"),
      plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.margin      = margin(3, 3, 3, 3)
    ) +
    labs(
      title = paste0(nm, " (Top ", ntaxa, ")"),
      x     = "Sample source",
      y     = "Relative abundance (%)",
      fill  = nm
    ) +
    guides(fill = guide_legend(ncol = LEGEND_NCOL, byrow = TRUE), colour = "none")
  
  out_dir_rank <- file.path(output_root, nm)
  out_name     <- paste0("BAR_", nm, "_top", ntaxa, "_TSS")
  
  save_plot_png(
    out_dir = out_dir_rank,
    name    = out_name,
    plot    = p,
    w       = BAR_WIDTH,
    h       = BAR_HEIGHT,
    dpi     = DPI
  )
  
  ## ---- Group-mean barplot --------------------------------------------------
  if (MAKE_GROUPMEAN_PLOTS) {
    
    stamp(
      "Plotting group-mean: ", nm,
      " | taxrank=", taxrank,
      " | top=", ntaxa,
      " | groupmean=", groupmean_var
    )
    
    gm_cfg <- get_special_groupmean_cfg(
      rank_name           = nm,
      default_width       = GROUPMEAN_WIDTH,
      default_height      = GROUPMEAN_HEIGHT,
      default_barwidth    = BARWIDTH,
      default_angle       = GROUP_X_ANGLE,
      default_legend_ncol = LEGEND_NCOL
    )
    
    t1_mean <- trans_abund$new(
      dataset   = mt_groupmean,
      taxrank   = taxrank,
      ntaxa     = ntaxa,
      groupmean = groupmean_var
    )
    
    p_mean <- t1_mean$plot_bar(
      others_color       = OTHERS_COLOR,
      color_values       = rank_colors,
      legend_text_italic = FALSE,
      barwidth           = gm_cfg$barwidth
    )
    
    if (isTRUE(gm_cfg$wrap_x_labels)) {
      p_mean <- p_mean +
        scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = gm_cfg$wrap_width))
    }
    
    this_groupmean_title <- if (!is.null(gm_cfg$custom_title)) {
      gm_cfg$custom_title
    } else {
      paste0(nm, " group-mean composition (Top ", ntaxa, ")")
    }
    
    p_mean <- p_mean +
      theme_classic(base_family = BASE_FAMILY) +
      theme(
        axis.title.x     = element_text(size = 15, face = "bold"),
        axis.title.y     = element_text(size = 15, face = "bold"),
        axis.text.y      = element_text(size = 10),
        axis.text.x      = element_text(angle = gm_cfg$x_angle, hjust = 1, vjust = 1, size = gm_cfg$x_text_size),
        legend.position  = "right",
        legend.title     = element_text(size = gm_cfg$legend_title_size, face = "bold"),
        legend.text      = element_text(size = gm_cfg$legend_text_size),
        legend.key.size  = unit(gm_cfg$legend_key_size, "cm"),
        legend.spacing   = unit(0.05, "cm"),
        plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.margin      = margin(3, 3, 3, 3)
      ) +
      labs(
        title = this_groupmean_title,
        x     = "Sample source groups",
        y     = "Mean relative abundance (%)",
        fill  = nm
      ) +
      guides(fill = guide_legend(ncol = gm_cfg$legend_ncol, byrow = TRUE), colour = "none")
    
    out_dir_groupmean  <- file.path(out_dir_rank, "GroupMean")
    out_name_groupmean <- paste0("BAR_GROUPMEAN_", nm, "_top", ntaxa, "_", groupmean_var, "_TSS")
    
    save_plot_png(
      out_dir = out_dir_groupmean,
      name    = out_name_groupmean,
      plot    = p_mean,
      w       = gm_cfg$width,
      h       = gm_cfg$height,
      dpi     = DPI
    )
  }
}

stamp("✅ Done. Output folder: ", normalizePath(output_root))