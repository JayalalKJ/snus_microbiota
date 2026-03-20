#!/usr/bin/env Rscript

# =============================================================================
# Alpha diversity from skin samples with rarefaction decision step
# Corrected version:
#   - uses cleaned non-rarefied object as input
#   - plots rarefaction curve first
#   - checks read depth
#   - chooses rarefaction depth
#   - rarefies using microtable$rarefy_samples()
#   - calculates alpha diversity from the rarefied object
# =============================================================================

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# 1. Load packages
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(microeco)
  library(magrittr)
  library(mecodev)
  library(dplyr)
  library(tibble)
  library(readr)
  library(cowplot)
})

# -----------------------------------------------------------------------------
# 2. Create output directory
# -----------------------------------------------------------------------------
output_dir <- "Results/Results_Step_03_A_rarefaction_Skin_for_alpha_diet"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("Output directory:\n")
print(normalizePath(output_dir))

# -----------------------------------------------------------------------------
# 3. Load cleaned non-rarefied object
# -----------------------------------------------------------------------------
input_nonrarefied <- "Results/Results_02C_contaminant_removal_from_blanks/microeco_cleaned_after_blank_contaminants.rds"

if (!file.exists(input_nonrarefied)) {
  stop("❌ Input file not found:\n", input_nonrarefied)
}

amplicon_16S_microtable <- readRDS(input_nonrarefied)

if (!inherits(amplicon_16S_microtable, "microtable")) {
  stop("❌ Input object is not a microeco::microtable")
}

cat("\nCleaned non-rarefied object loaded successfully.\n")
cat("OTU table dimensions:\n")
print(dim(amplicon_16S_microtable$otu_table))

# -----------------------------------------------------------------------------
# 4. Keep only skin diet groups
# -----------------------------------------------------------------------------


sampletype_keep <- c("Control", "Dose_01%", "Dose_2_2%")


st0 <- as.data.frame(amplicon_16S_microtable$sample_table)
colnames(st0) <- make.unique(colnames(st0))

if (!"SampleTypes" %in% colnames(st0)) {
  stop("❌ Column 'SampleTypes' not found in sample_table")
}

st0$SampleTypes <- trimws(as.character(st0$SampleTypes))

keep_ids <- rownames(st0)[st0$SampleTypes %in% sampletype_keep]

if (length(keep_ids) == 0) {
  stop("❌ No samples found for SampleTypes: ", paste(sampletype_keep, collapse = ", "))
}

keep_ids <- intersect(keep_ids, colnames(amplicon_16S_microtable$otu_table))

if (length(keep_ids) == 0) {
  stop("❌ No matching sample IDs between sample_table and otu_table")
}

amplicon_16S_microtable$otu_table   <- amplicon_16S_microtable$otu_table[, keep_ids, drop = FALSE]
amplicon_16S_microtable$sample_table <- st0[keep_ids, , drop = FALSE]
amplicon_16S_microtable$tidy_dataset()

cat("\nAfter filtering to skin diet groups:\n")
print(table(amplicon_16S_microtable$sample_table$SampleTypes))

# -----------------------------------------------------------------------------
# 5. Save filtered non-rarefied object
# -----------------------------------------------------------------------------
saveRDS(
  amplicon_16S_microtable,
  file.path(output_dir, "microeco_skin_filtered_nonrarefied.rds")
)

cat("\nSaved filtered non-rarefied object:\n")
print(file.path(output_dir, "microeco_skin_filtered_nonrarefied.rds"))

# -----------------------------------------------------------------------------
# 6. Rarefaction curve based on Shannon diversity
# -----------------------------------------------------------------------------
tmp_rarefy <- trans_rarefy$new(
  amplicon_16S_microtable,
  alphadiv = c("Shannon"),
  depth = c(
    0, 10, 50, 500, 1000, 2000, 4000, 6000,
    10000, 15000, 20000, 25000, 30000, 35000,
    40000, 45000, 50000, 55000, 60000, 65000, 70000
  )
)

g1 <- tmp_rarefy$plot_rarefy(
  show_samplename = TRUE,
  color_values = rep("grey50", 100),
  show_legend = FALSE
)

cowplot::save_plot(
  file.path(output_dir, "AlphaDiv_Rarefactioncurve_Shannon.png"),
  g1,
  base_aspect_ratio = 1.4,
  dpi = 300,
  base_height = 5
)

cat("\nSaved rarefaction curve:\n")
print(file.path(output_dir, "AlphaDiv_Rarefactioncurve_Shannon.png"))

# -----------------------------------------------------------------------------
# 6A. Check sequencing depth after plotting rarefaction curve
# -----------------------------------------------------------------------------
sample_depths <- amplicon_16S_microtable$sample_sums()

depth_range <- range(sample_depths)
min_reads <- depth_range[1]
max_reads <- depth_range[2]

cat("\nRead depth range after filtering:\n")
cat("Minimum reads:", min_reads, "\n")
cat("Maximum reads:", max_reads, "\n")

depth_table <- data.frame(
  SampleID = names(sample_depths),
  Reads = as.numeric(sample_depths),
  stringsAsFactors = FALSE
)

write.csv(
  depth_table,
  file.path(output_dir, "Sample_read_depths.csv"),
  row.names = FALSE
)

cat("\nSaved sample read-depth table:\n")
print(file.path(output_dir, "Sample_read_depths.csv"))

# -----------------------------------------------------------------------------
# 6B. Decide rarefaction depth
# -----------------------------------------------------------------------------
rarefy_depth <- 20000

cat("\nChosen rarefaction depth:\n")
print(rarefy_depth)

if (min_reads < rarefy_depth) {
  stop(
    "❌ Chosen rarefaction depth is too high.\n",
    "Minimum reads in filtered samples = ", min_reads, "\n",
    "Chosen rarefaction depth = ", rarefy_depth, "\n",
    "Please lower rarefy_depth after checking the rarefaction curve."
  )
}

cat("\n✅ Rarefaction depth is valid for all filtered samples.\n")

# -----------------------------------------------------------------------------
# 6C. Rarefy the filtered object
# Why:
# Use rarefy_samples() because this is a microtable object.
# -----------------------------------------------------------------------------
tmp_microtable <- clone(amplicon_16S_microtable)

amplicon_16S_microtable_rarefy <- tmp_microtable$rarefy_samples(
  sample.size = rarefy_depth
)

cat("\nRarefaction completed successfully.\n")
cat("Rarefied OTU table dimensions:\n")
print(dim(amplicon_16S_microtable_rarefy$otu_table))

saveRDS(
  amplicon_16S_microtable_rarefy,
  file.path(output_dir, paste0("microeco_skin_rarefied_", rarefy_depth, ".rds"))
)

cat("\nSaved rarefied object:\n")
print(file.path(output_dir, paste0("microeco_skin_rarefied_", rarefy_depth, ".rds")))

# -----------------------------------------------------------------------------
# 7. Calculate alpha diversity from rarefied abundance table
# -----------------------------------------------------------------------------
tmp_microtable <- clone(amplicon_16S_microtable_rarefy)

if (is.null(tmp_microtable$phylo_tree)) {
  stop("❌ Faith's PD cannot be calculated because phylo_tree is missing in the microtable.")
}

tmp_microtable$cal_alphadiv(
  measures = c(
    "Observed",
    "Coverage",
    "Chao1",
    "ACE",
    "Shannon",
    "Simpson",
    "InvSimpson",
    "Fisher",
    "Pielou"
  ),
  PD = TRUE
)

cat("\nAlpha diversity calculated.\n")
cat("Alpha diversity columns:\n")
print(colnames(tmp_microtable$alpha_diversity))

cat("\nFirst rows of alpha_diversity table:\n")
print(head(tmp_microtable$alpha_diversity))

# -----------------------------------------------------------------------------
# 8. Export alpha diversity table
# -----------------------------------------------------------------------------
write.csv(
  tmp_microtable$alpha_diversity,
  file.path(output_dir, "AlphaDiv_metrics.csv"),
  row.names = TRUE
)

cat("\nSaved alpha diversity table:\n")
print(file.path(output_dir, "AlphaDiv_metrics.csv"))

# -----------------------------------------------------------------------------
# 9. Save microtable object with alpha diversity
# -----------------------------------------------------------------------------
saveRDS(
  tmp_microtable,
  file.path(output_dir, "microeco_Skin_rarefied_with_alphadiv.rds")
)

cat("\nSaved alpha-diversity microtable object:\n")
print(file.path(output_dir, "microeco_Skin_rarefied_with_alphadiv.rds"))

cat("\n✅ Alpha diversity workflow completed successfully.\n")