#!/usr/bin/env Rscript

# =============================================================================
# Step_01_import_mothur_to_microeco.R
# Purpose : Import mothur shared + taxonomy + tree + metadata into phyloseq
#           and convert to microeco object with QC reports
# Author  : Jayalal K Jayanthan
# Updated : 2026-03-11
# =============================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(phyloseq)
  library(microeco)
  library(file2meco)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ape)
})

# -----------------------------------------------------------------------------
# 0) Logging
# -----------------------------------------------------------------------------
stamp <- function(...) {
  message(sprintf("[%s] %s",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  paste0(..., collapse = "")))
}

# -----------------------------------------------------------------------------
# 1) I/O
# -----------------------------------------------------------------------------
INPUT_DIR  <- "Data/97_metadata"
OUTPUT_DIR <- "results/Results_01_import_mothur_to_microeco"
REPORT_DIR <- file.path(OUTPUT_DIR, "reports")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

FILES <- list(
  metadata = "metadata.csv",
  shared   = "AlgOpti1.skin.trim.contigs.renamed.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.shared",
  taxonomy = "AlgOpti1.skin.trim.contigs.renamed.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.0.03.cons.taxonomy",
  tree     = "AlgOpti1.skin.trim.contigs.renamed.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.0.03.rep.phylip.tre"
)

meta_path   <- file.path(INPUT_DIR, FILES$metadata)
shared_path <- file.path(INPUT_DIR, FILES$shared)
tax_path    <- file.path(INPUT_DIR, FILES$taxonomy)
tree_path   <- file.path(INPUT_DIR, FILES$tree)

all_paths <- c(meta_path, shared_path, tax_path, tree_path)
missing_files <- all_paths[!file.exists(all_paths)]

if (length(missing_files) > 0) {
  stop("❌ Missing input file(s):\n", paste(" -", missing_files, collapse = "\n"))
}

# -----------------------------------------------------------------------------
# 2) Helpers
# -----------------------------------------------------------------------------
trim_chr <- function(x) {
  trimws(as.character(x))
}

safe_write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(as.data.frame(df), path)
}

clean_rank_prefix <- function(x) {
  x <- trim_chr(x)
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x <- stringr::str_replace_all(x, "^g__", "")
  x <- stringr::str_replace_all(x, "^g_", "")
  trim_chr(x)
}

is_unclassified_token <- function(x_clean) {
  x0 <- tolower(trim_chr(x_clean))
  x0 %in% c("", "na", "n/a", "unknown", "unclassified", "uncultured", "unidentified", "metagenome") |
    grepl("uncultured|unidentified|unclassified|unknown|metagenome", x0)
}

# -----------------------------------------------------------------------------
# 3) Metadata
# -----------------------------------------------------------------------------
stamp("Reading metadata ...")

metadata <- readr::read_csv(meta_path, show_col_types = FALSE)

if (!"Sample_ID" %in% names(metadata)) {
  stop("❌ Metadata file must contain a column named 'Sample_ID'")
}

metadata <- metadata %>%
  dplyr::mutate(Sample_ID = trim_chr(Sample_ID)) %>%
  dplyr::filter(!is.na(Sample_ID), Sample_ID != "") %>%
  dplyr::distinct(Sample_ID, .keep_all = TRUE)

stamp("✅ metadata: ", nrow(metadata), " rows × ", ncol(metadata), " cols")

# -----------------------------------------------------------------------------
# 4) Shared -> OTU table
# -----------------------------------------------------------------------------
stamp("Reading mothur shared file ...")

otu_df <- readr::read_tsv(shared_path, show_col_types = FALSE)

req_cols <- c("label", "Group", "numOtus")
miss <- setdiff(req_cols, names(otu_df))
if (length(miss) > 0) {
  stop("❌ mothur shared file missing required columns: ", paste(miss, collapse = ", "))
}

otu_counts <- otu_df %>%
  dplyr::select(-dplyr::all_of(c("label", "numOtus"))) %>%
  dplyr::rename(Sample_ID = Group) %>%
  dplyr::mutate(Sample_ID = trim_chr(Sample_ID))

otu_cols <- setdiff(names(otu_counts), "Sample_ID")

otu_counts_num <- otu_counts %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(otu_cols),
      ~ suppressWarnings(as.integer(.x))
    )
  )

# detect coercion issues only where original was not NA/blank
bad_count_cols <- otu_cols[
  vapply(otu_cols, function(col) {
    original <- otu_counts[[col]]
    converted <- otu_counts_num[[col]]
    any(!is.na(original) & trim_chr(original) != "" & is.na(converted))
  }, logical(1))
]

if (length(bad_count_cols) > 0) {
  safe_write_tsv(
    tibble::tibble(OTU_column = bad_count_cols),
    file.path(REPORT_DIR, "shared_columns_with_integer_coercion_issues.tsv")
  )
  stamp("⚠️ Some OTU count columns produced NA during integer coercion; report written.")
}

otu_mat <- otu_counts_num %>%
  tibble::column_to_rownames("Sample_ID") %>%
  as.matrix()

storage.mode(otu_mat) <- "integer"

if (any(is.na(otu_mat))) {
  stop("❌ OTU matrix contains NA values after coercion. Check shared file and coercion report.")
}

stamp("✅ shared OTU table: ", nrow(otu_mat), " samples × ", ncol(otu_mat), " OTUs")

# -----------------------------------------------------------------------------
# 5) Taxonomy
# -----------------------------------------------------------------------------
stamp("Reading taxonomy file ...")

tax_raw <- readr::read_tsv(
  tax_path,
  col_types = readr::cols(
    OTU      = readr::col_character(),
    Taxonomy = readr::col_character(),
    Size     = readr::col_double()
  )
)

req_tax_cols <- c("OTU", "Taxonomy")
miss_tax <- setdiff(req_tax_cols, names(tax_raw))
if (length(miss_tax) > 0) {
  stop("❌ taxonomy file missing required columns: ", paste(miss_tax, collapse = ", "))
}

tax_raw <- tax_raw %>%
  dplyr::rename(
    otu      = OTU,
    taxonomy = Taxonomy,
    size     = Size
  ) %>%
  dplyr::mutate(
    taxonomy = taxonomy %>%
      stringr::str_remove_all("\\(\\d+\\)") %>%
      stringr::str_remove_all(";$")
  ) %>%
  tidyr::separate(
    taxonomy,
    into  = c("Kingdoms", "Phyla", "Classes", "Orders", "Families", "Genera"),
    sep   = ";",
    fill  = "right",
    extra = "drop"
  ) %>%
  dplyr::mutate(
    dplyr::across(
      c("Kingdoms", "Phyla", "Classes", "Orders", "Families", "Genera"),
      ~ dplyr::na_if(stringr::str_squish(.x), "")
    )
  ) %>%
  dplyr::distinct(otu, .keep_all = TRUE)

stamp("✅ taxonomy rows: ", nrow(tax_raw), " OTUs")

# genus diagnostics
genus_freq <- tibble::tibble(
  Genera_raw = as.character(tax_raw$Genera)
) %>%
  dplyr::mutate(
    Genera_raw = ifelse(is.na(Genera_raw), "NA", trim_chr(Genera_raw))
  ) %>%
  dplyr::count(Genera_raw, sort = TRUE)

safe_write_tsv(genus_freq, file.path(REPORT_DIR, "genus_raw_frequency.tsv"))

# robust genus cleaning
genus_raw   <- as.character(tax_raw$Genera)
genus_clean <- clean_rank_prefix(genus_raw)
unclassified_mask <- is.na(genus_raw) | is_unclassified_token(genus_clean)

tax_raw <- tax_raw %>%
  dplyr::mutate(
    Genus_clean = ifelse(unclassified_mask, "Unclassified_genus", genus_clean),
    Genera = Genus_clean,
    Genera_level_OTUs = paste0(Genus_clean, "_", otu)
  )

unclassified_log <- tax_raw %>%
  dplyr::transmute(
    OTU         = otu,
    Genera_raw  = genus_raw,
    Genus_clean = Genus_clean,
    Label_used  = Genera_level_OTUs
  ) %>%
  dplyr::filter(Genus_clean == "Unclassified_genus")

safe_write_tsv(unclassified_log, file.path(REPORT_DIR, "unclassified_genus_OTUs.tsv"))
stamp("✅ unclassified genus OTUs: ", nrow(unclassified_log), " (logged)")

tax_raw <- tax_raw %>%
  tibble::column_to_rownames("otu")

# -----------------------------------------------------------------------------
# 6) Harmonize OTUs between shared and taxonomy
# -----------------------------------------------------------------------------
stamp("Harmonizing OTUs between shared and taxonomy ...")

otus_common <- intersect(colnames(otu_mat), rownames(tax_raw))
if (length(otus_common) == 0) {
  stop("❌ No overlapping OTU IDs between shared and taxonomy")
}

dropped_from_otu <- setdiff(colnames(otu_mat), otus_common)
dropped_from_tax <- setdiff(rownames(tax_raw), otus_common)

if (length(dropped_from_otu) > 0) {
  safe_write_tsv(
    tibble::tibble(OTU = dropped_from_otu),
    file.path(REPORT_DIR, "otus_dropped_missing_in_taxonomy.tsv")
  )
}

if (length(dropped_from_tax) > 0) {
  safe_write_tsv(
    tibble::tibble(OTU = dropped_from_tax),
    file.path(REPORT_DIR, "otus_dropped_missing_in_shared.tsv")
  )
}

otu_mat <- otu_mat[, otus_common, drop = FALSE]

tax_df <- tax_raw[
  otus_common,
  c("Kingdoms", "Phyla", "Classes", "Orders", "Families", "Genera", "Genus_clean", "Genera_level_OTUs"),
  drop = FALSE
]

# -----------------------------------------------------------------------------
# 7) Harmonize Samples
# -----------------------------------------------------------------------------
stamp("Harmonizing samples between metadata and shared file ...")

samples_before <- rownames(otu_mat)
meta_ids <- metadata$Sample_ID

common_samples <- samples_before[samples_before %in% meta_ids]
lost_no_meta   <- setdiff(samples_before, meta_ids)
meta_only      <- setdiff(meta_ids, samples_before)

if (length(lost_no_meta) > 0) {
  safe_write_tsv(
    tibble::tibble(Sample_ID = lost_no_meta),
    file.path(REPORT_DIR, "samples_dropped_no_metadata.tsv")
  )
}

if (length(meta_only) > 0) {
  safe_write_tsv(
    tibble::tibble(Sample_ID = meta_only),
    file.path(REPORT_DIR, "metadata_rows_without_counts.tsv")
  )
}

otu_mat <- otu_mat[common_samples, , drop = FALSE]

metadata_ps <- metadata %>%
  dplyr::filter(Sample_ID %in% common_samples) %>%
  dplyr::arrange(match(Sample_ID, rownames(otu_mat))) %>%
  tibble::column_to_rownames("Sample_ID")

if (!identical(rownames(metadata_ps), rownames(otu_mat))) {
  stop("❌ Sample order mismatch between metadata and OTU table after harmonization")
}

stamp("✅ working dataset: ", nrow(otu_mat), " samples × ", ncol(otu_mat), " OTUs")

# -----------------------------------------------------------------------------
# 8) Tree
# -----------------------------------------------------------------------------
stamp("Reading and harmonizing tree ...")

tr <- ape::read.tree(tree_path)
if (is.null(tr$tip.label)) {
  stop("❌ Tree has no tip labels")
}

tips_common         <- intersect(tr$tip.label, colnames(otu_mat))
tips_only_in_tree   <- setdiff(tr$tip.label, colnames(otu_mat))
otus_missing_in_tree <- setdiff(colnames(otu_mat), tr$tip.label)

if (length(tips_only_in_tree) > 0) {
  safe_write_tsv(
    tibble::tibble(OTU = tips_only_in_tree),
    file.path(REPORT_DIR, "tree_tips_dropped_not_in_otu_table.tsv")
  )
}

if (length(otus_missing_in_tree) > 0) {
  safe_write_tsv(
    tibble::tibble(OTU = otus_missing_in_tree),
    file.path(REPORT_DIR, "otus_dropped_missing_in_tree.tsv")
  )
}

if (length(tips_common) == 0) {
  stop("❌ No overlapping OTU IDs between tree tips and OTU table")
}

tr <- ape::drop.tip(tr, setdiff(tr$tip.label, tips_common))
otu_mat <- otu_mat[, tips_common, drop = FALSE]
tax_df  <- tax_df[tips_common, , drop = FALSE]

# force identical order across tree / otu / taxonomy
otu_mat <- otu_mat[, tr$tip.label, drop = FALSE]
tax_df  <- tax_df[tr$tip.label, , drop = FALSE]

if (!identical(colnames(otu_mat), rownames(tax_df))) {
  stop("❌ OTU/taxonomy order mismatch before phyloseq construction")
}

# -----------------------------------------------------------------------------
# 9) Build phyloseq
# -----------------------------------------------------------------------------
stamp("Building phyloseq object ...")

ps <- phyloseq::phyloseq(
  phyloseq::otu_table(t(otu_mat), taxa_are_rows = TRUE),
  phyloseq::tax_table(as.matrix(tax_df)),
  phyloseq::sample_data(metadata_ps),
  phyloseq::phy_tree(tr)
)

saveRDS(ps, file.path(OUTPUT_DIR, "phyloseq_raw.rds"))
stamp("💾 saved: ", file.path(OUTPUT_DIR, "phyloseq_raw.rds"))

# -----------------------------------------------------------------------------
# 10) Convert to microeco
# -----------------------------------------------------------------------------
stamp("Converting phyloseq to microeco ...")

mt <- file2meco::phyloseq2meco(ps)
mt$tidy_dataset()
mt$cal_abund()
print(mt$tidy_dataset())
print(mt$cal_abund())




saveRDS(mt, file.path(OUTPUT_DIR, "microeco_raw.rds"))
stamp("💾 saved: ", file.path(OUTPUT_DIR, "microeco_raw.rds"))

# -----------------------------------------------------------------------------
# 11) QC exports
# -----------------------------------------------------------------------------
stamp("Writing QC reports ...")

qc <- tibble::tibble(
  n_samples    = phyloseq::nsamples(ps),
  n_otus       = phyloseq::ntaxa(ps),
  min_depth    = min(phyloseq::sample_sums(ps)),
  median_depth = as.numeric(stats::median(phyloseq::sample_sums(ps))),
  max_depth    = max(phyloseq::sample_sums(ps))
)
safe_write_tsv(qc, file.path(REPORT_DIR, "import_qc_summary.tsv"))

depths <- tibble::tibble(
  Sample_ID = names(phyloseq::sample_sums(ps)),
  LibSize   = as.integer(phyloseq::sample_sums(ps))
)
safe_write_tsv(depths, file.path(REPORT_DIR, "library_sizes.tsv"))

stamp("Step 01 complete")
stamp("• Output: ", normalizePath(OUTPUT_DIR))