#!/usr/bin/env Rscript
# =============================================================================
# (METADATA-DRIVEN NEGATIVES + BUCKETED OTUs)
# Uses SampleTypes column only:
#   - Negative samples removed: SampleTypes %in% c("Blank","Water")
# Silences contaminant OTUs using BUCKET RULES:
#   - Always silence: kit_water, human, enteric_flags (Zymo-like), hood_dust, rare_unexpected
#   - Conditionally silence: aquatic_ambiguous ONLY if blank-driven (prev/reads thresholds)
# =============================================================================

suppressPackageStartupMessages({
  library(microeco)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

# ---- I/O -------------------------------------------------------------------
input_rds  <- "results/Results_01_import_mothur_to_microeco/microeco_raw.rds"
output_dir <- "results/Results_02C_contaminant_removal_from_blanks"
fig_dir    <- file.path(output_dir, "figures")
rep_dir    <- file.path(output_dir, "reports")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(rep_dir,    recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_rds))

mt <- readRDS(input_rds)
if (!inherits(mt, "microtable")) stop("❌ Input is not a microeco::microtable")

# ---- Config ----------------------------------------------------------------
NEG_TYPE_COL   <- "SampleTypes"
NEG_TYPE_LEVEL <- "Blank"   # IMPORTANT: "Control" is a fish group, not a negative

manual_drop_ids <- character(0)  # optional, e.g. c("SM10")

# Blank-driven rule for "aquatic_ambiguous"
TIER_C_MIN_BLANK_PREV  <- 0.20
TIER_C_MIN_BLANK_READS <- 10

make_qc_plots <- TRUE
top_n <- 30

# ---- Bucketed contaminant TaxLabels (OTU-level strings) ---------------------
contaminant_taxa_kit_water <- c(
  "g__Ralstonia_Otu00173",
  "g__Bradyrhizobium_Otu00077",
  "g__Methylobacterium_Otu00013",
  "g__Methylobacterium_Otu00024",
  "g__Methylobacterium_Otu00052",
  "g__Methylobacterium_Otu00068",
  "g__Methylobacterium_Otu00085",
  "g__Sphingomonas_Otu00003",
  "g__Sphingomonas_Otu00008",
  "g__Sphingomonas_Otu00197",
  "g__Novosphingobium_Otu00006",
  "g__Novosphingobium_Otu00045",
  "g__Sphingobium_Otu00026",
  "g__Brevundimonas_Otu00094",
  "g__Brevundimonas_Otu00086",
  "g__Caulobacter_Otu00036",
  "g__Phyllobacterium_Otu00011"
)

contaminant_taxa_human <- c(
  "g__Cutibacterium_Otu00013",
  "g__Cutibacterium_Otu00025",
  "g__Staphylococcus_Otu00047",
  "g__Micrococcus_Otu00043",
  "g__Enterococcus_Otu00071"
)

# Zymo/mock-like / enteric flags (often bleed/carryover)
contaminant_taxa_enteric_flags <- c(
  "g__Escherichia-Shigella_Otu00019",
  "g__Salmonella_Otu00018",
  "g__Listeria_Otu00076"
)

contaminant_taxa_hood_dust <- c(
  "g__Bacillus_Otu00017"
)

# Aquatic-but-ambiguous: silence ONLY if blank-driven
contaminant_taxa_aquatic_ambiguous <- c(
  "g__Pseudomonas_Otu00039",
  "g__Pseudomonas_Otu00035",
  "g__Pseudomonas_Otu00145",
  "g__Pseudomonas_Otu00103",
  "g__Pseudomonas_Otu00107",
  "g__Acinetobacter_Otu00006",
  "g__Acinetobacter_Otu00009",
  "g__Acinetobacter_Otu00005",
  "g__Acinetobacter_Otu00010",
  "g__Acinetobacter_Otu00084",
  "g__Acidovorax_Otu00029",
  "g__Acidovorax_Otu00048",
  "g__Acidovorax_Otu00289",
  "g__Burkholderia-Caballeronia-Paraburkholderia_Otu00092",
  "g__Hydrogenophaga_Otu00138",
  "g__Curvibacter_Otu00517",
  "g__Microbacterium_Otu00015",
  "g__Cloacibacterium_Otu00250",
  "g__Methyloversatilis_Otu00003",
  "g__Methyloversatilis_Otu00004"
)

contaminant_taxa_rare_unexpected <- c(
  "g__Truepera_Otu00002",
  "g__Meiothermus_Otu00033"
)

contaminant_taxa <- unique(c(
  contaminant_taxa_kit_water,
  contaminant_taxa_human,
  contaminant_taxa_enteric_flags,
  contaminant_taxa_hood_dust,
  contaminant_taxa_aquatic_ambiguous,
  contaminant_taxa_rare_unexpected
))

# ---- Helpers ----------------------------------------------------------------
ensure_dir <- function(p) { dir.create(p, recursive = TRUE, showWarnings = FALSE); invisible(p) }

safe_write_tsv <- function(df, path) {
  ensure_dir(dirname(path))
  readr::write_tsv(as.data.frame(df), path)
}

pick_tax_col <- function(tt_cols) {
  if ("Genera_level_OTUs" %in% tt_cols) return("Genera_level_OTUs")
  if ("Genus_level_OTUs"  %in% tt_cols) return("Genus_level_OTUs")
  if ("Genus"             %in% tt_cols) return("Genus")
  if ("Genera"            %in% tt_cols) return("Genera")
  NA_character_
}

sync_to_sample_id <- function(m) {
  st <- as.data.frame(m$sample_table)
  if (!"Sample_ID" %in% names(st)) return(m)
  
  new_rn <- as.character(st$Sample_ID)
  old_rn <- rownames(st)
  
  rownames(st) <- new_rn
  m$sample_table <- st
  
  otu <- m$otu_table
  if (identical(colnames(otu), old_rn)) colnames(otu) <- new_rn
  
  if (all(new_rn %in% colnames(otu))) {
    otu <- otu[, new_rn, drop = FALSE]
  } else {
    missing <- setdiff(new_rn, colnames(otu))
    stop("❌ Cannot align otu_table columns to Sample_ID. Missing: ", paste(missing, collapse = ", "))
  }
  
  m$otu_table <- otu
  stopifnot(identical(rownames(m$sample_table), colnames(m$otu_table)))
  m
}

get_negative_ids_from_sampletypes <- function(m) {
  st <- as.data.frame(m$sample_table)
  if (!NEG_TYPE_COL %in% names(st)) stop("❌ sample_table missing column: ", NEG_TYPE_COL)
  neg_mask <- as.character(st[[NEG_TYPE_COL]]) %in% NEG_TYPE_LEVEL
  rownames(st)[neg_mask]
}

drop_samples <- function(m, ids, label, out_log) {
  st <- as.data.frame(m$sample_table)
  rm_mask <- rownames(st) %in% ids
  removed <- rownames(st)[rm_mask]
  if (!length(removed)) return(list(mt = m, removed = character(0)))
  
  safe_write_tsv(tibble(Removed_Sample_ID = removed, Reason = label), out_log)
  
  m$sample_table <- m$sample_table[!rm_mask, , drop = FALSE]
  m$otu_table    <- m$otu_table[, rownames(m$sample_table), drop = FALSE]
  stopifnot(identical(rownames(m$sample_table), colnames(m$otu_table)))
  
  list(mt = m, removed = removed)
}

prune_zero_otus <- function(m) {
  keep <- rownames(m$otu_table)[rowSums(m$otu_table) > 0]
  m$otu_table <- m$otu_table[keep, , drop = FALSE]
  m$tax_table <- m$tax_table[keep, , drop = FALSE]
  m
}

plot_top_taxa <- function(m, tax_col, out_stub) {
  tt <- as.data.frame(m$tax_table)
  ok <- !is.na(tt[[tax_col]]) & tt[[tax_col]] != ""
  tt <- tt[ok, , drop = FALSE]
  
  counts <- rowSums(m$otu_table[rownames(tt), , drop = FALSE])
  df <- tibble(Taxon = as.character(tt[[tax_col]]), Reads = as.numeric(counts)) %>%
    group_by(Taxon) %>%
    summarise(Reads = sum(Reads), .groups = "drop") %>%
    arrange(desc(Reads)) %>%
    slice_head(n = top_n)
  
  p <- ggplot(df, aes(x = reorder(Taxon, Reads), y = Reads)) +
    geom_col() +
    coord_flip() +
    theme_classic(base_size = 13) +
    labs(x = NULL, y = "Total reads", title = NULL)
  
  ggsave(paste0(out_stub, ".png"), p, width = 10, height = 7, dpi = 450, bg = "white")
  ggsave(paste0(out_stub, ".pdf"), p, width = 10, height = 7, device = "pdf")
  invisible(p)
}

# ---- FIXED: vectorised '&' instead of '&&' inside mutate/case_when ----------
score_bucketed_contaminants <- function(m, tax_col, rep_dir) {
  bucket_map <- dplyr::bind_rows(
    tibble(TaxLabel = contaminant_taxa_kit_water,         Bucket = "kit_water"),
    tibble(TaxLabel = contaminant_taxa_human,             Bucket = "human"),
    tibble(TaxLabel = contaminant_taxa_enteric_flags,     Bucket = "enteric_flags"),
    tibble(TaxLabel = contaminant_taxa_hood_dust,         Bucket = "hood_dust"),
    tibble(TaxLabel = contaminant_taxa_aquatic_ambiguous, Bucket = "aquatic_ambiguous"),
    tibble(TaxLabel = contaminant_taxa_rare_unexpected,   Bucket = "rare_unexpected")
  ) %>% distinct()
  
  st <- as.data.frame(m$sample_table)
  neg_ids <- rownames(st)[as.character(st[[NEG_TYPE_COL]]) %in% NEG_TYPE_LEVEL]
  pos_ids <- setdiff(rownames(st), neg_ids)
  
  safe_write_tsv(
    tibble(Negative_Sample_ID = neg_ids,
           SampleTypes = as.character(st[neg_ids, NEG_TYPE_COL, drop = TRUE])),
    file.path(rep_dir, "negative_samples_detected.tsv")
  )
  
  tax_df <- as.data.frame(m$tax_table)
  otu_map <- tibble(
    OTU = rownames(tax_df),
    TaxLabel = as.character(tax_df[[tax_col]])
  )
  
  cand <- otu_map %>% inner_join(bucket_map, by = "TaxLabel")
  
  if (!nrow(cand)) {
    safe_write_tsv(bucket_map, file.path(rep_dir, "bucket_map_taxlabels.tsv"))
    stop("❌ No bucketed TaxLabels matched your tax_col: ", tax_col,
         "\n   (Wrote bucket_map_taxlabels.tsv so you can compare labels.)")
  }
  
  otu_mat <- m$otu_table[cand$OTU, , drop = FALSE]
  calc_prev <- function(x) mean(x > 0)
  calc_sum  <- function(x) sum(x)
  
  cand$blank_prev  <- if (length(neg_ids)) apply(otu_mat[, neg_ids, drop = FALSE], 1, calc_prev) else 0
  cand$blank_reads <- if (length(neg_ids)) apply(otu_mat[, neg_ids, drop = FALSE], 1, calc_sum)  else 0
  cand$samp_prev   <- if (length(pos_ids)) apply(otu_mat[, pos_ids, drop = FALSE], 1, calc_prev) else 0
  cand$samp_reads  <- if (length(pos_ids)) apply(otu_mat[, pos_ids, drop = FALSE], 1, calc_sum)  else 0
  
  always_remove_buckets <- c("kit_water","human","enteric_flags","hood_dust","rare_unexpected")
  
  cand <- cand %>%
    mutate(
      Remove = case_when(
        Bucket %in% always_remove_buckets ~ TRUE,
        
        # aquatic_ambiguous: remove ONLY if blank-driven (VECTORISED &)
        Bucket == "aquatic_ambiguous" ~ (
          blank_prev  >= TIER_C_MIN_BLANK_PREV &
            blank_reads >= TIER_C_MIN_BLANK_READS &
            blank_prev  >= samp_prev
        ),
        
        TRUE ~ FALSE
      ),
      Remove_reason = case_when(
        Bucket %in% always_remove_buckets ~ paste0("always_remove_", Bucket),
        Bucket == "aquatic_ambiguous" & Remove ~ "aquatic_blank_driven",
        Bucket == "aquatic_ambiguous" & !Remove ~ "aquatic_keep_not_blank_driven",
        TRUE ~ "keep"
      )
    )
  
  safe_write_tsv(cand, file.path(rep_dir, "candidate_contaminants_bucketed_scored.tsv"))
  
  remove_otus <- cand %>% filter(Remove) %>% pull(OTU)
  safe_write_tsv(
    cand %>% filter(Remove) %>%
      select(OTU, TaxLabel, Bucket, Remove_reason, blank_prev, blank_reads, samp_prev, samp_reads),
    file.path(rep_dir, "otus_to_silence_bucketed.tsv")
  )
  
  list(remove_otus = remove_otus, cand_tbl = cand, neg_ids = neg_ids)
}

# ---- Run --------------------------------------------------------------------
mt <- sync_to_sample_id(mt)

tax_col <- pick_tax_col(colnames(mt$tax_table))
if (is.na(tax_col)) stop("❌ No usable genus label column in mt$tax_table")

st0 <- as.data.frame(mt$sample_table)
safe_write_tsv(
  st0 %>% count(.data[[NEG_TYPE_COL]], name = "n") %>% arrange(desc(.data$n)),
  file.path(rep_dir, "SampleTypes_counts.tsv")
)

qc_before <- tibble(
  n_samples = nrow(mt$sample_table),
  n_otus    = nrow(mt$otu_table),
  min_depth = min(colSums(mt$otu_table)),
  median_depth = as.numeric(stats::median(colSums(mt$otu_table))),
  max_depth = max(colSums(mt$otu_table))
)
safe_write_tsv(qc_before, file.path(rep_dir, "qc_before.tsv"))
if (isTRUE(make_qc_plots)) plot_top_taxa(mt, tax_col, file.path(fig_dir, "top_taxa_before"))

if (length(manual_drop_ids)) {
  res <- drop_samples(
    mt, ids = manual_drop_ids, label = "manual_drop",
    out_log = file.path(rep_dir, "samples_dropped_manual.tsv")
  )
  mt <- res$mt
}

sc <- score_bucketed_contaminants(mt, tax_col, rep_dir)
remove_otus <- sc$remove_otus

safe_write_tsv(
  tibble(
    n_total_unique_taxlabels = length(contaminant_taxa),
    n_remove_otus = length(remove_otus),
    tierC_min_blank_prev = TIER_C_MIN_BLANK_PREV,
    tierC_min_blank_reads = TIER_C_MIN_BLANK_READS
  ),
  file.path(rep_dir, "contaminant_silencing_summary.tsv")
)

if (length(remove_otus)) mt$otu_table[remove_otus, ] <- 0
mt <- prune_zero_otus(mt)

negative_ids <- get_negative_ids_from_sampletypes(mt)
safe_write_tsv(
  tibble(Negative_Sample_ID = negative_ids,
         SampleTypes = as.character(mt$sample_table[negative_ids, NEG_TYPE_COL, drop = TRUE])),
  file.path(rep_dir, "negative_samples_removed_final.tsv")
)

if (length(negative_ids)) {
  res2 <- drop_samples(
    mt, ids = negative_ids, label = "negative_control",
    out_log = file.path(rep_dir, "samples_dropped_negative_controls.tsv")
  )
  mt <- prune_zero_otus(res2$mt)
}

mt$tidy_dataset()
mt$cal_abund()

qc_after <- tibble(
  n_samples = nrow(mt$sample_table),
  n_otus    = nrow(mt$otu_table),
  min_depth = min(colSums(mt$otu_table)),
  median_depth = as.numeric(stats::median(colSums(mt$otu_table))),
  max_depth = max(colSums(mt$otu_table))
)
safe_write_tsv(qc_after, file.path(rep_dir, "qc_after.tsv"))
if (isTRUE(make_qc_plots)) plot_top_taxa(mt, tax_col, file.path(fig_dir, "top_taxa_after"))

out_rds <- file.path(output_dir, "microeco_cleaned_after_blank_contaminants.rds")
saveRDS(mt, out_rds)

message("✅ Step 03 complete (SampleTypes negatives + bucketed contaminant OTUs)")
message("• Input : ", normalizePath(input_rds))
message("• Output: ", normalizePath(out_rds))
message("• Reports: ", normalizePath(rep_dir))
message("• Figures: ", normalizePath(fig_dir))