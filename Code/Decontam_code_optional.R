#!/usr/bin/env Rscript
# =============================================================================
#  Frequency-based decontam did not perform reliably for this dataset; therefore, 
#we used the extraction blank for visual separation and contaminant screening

options(stringsAsFactors = FALSE)
set.seed(100)

suppressPackageStartupMessages({
  library(phyloseq)
  library(decontam)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(file2meco)
  library(microeco)
})

# --------------------------- USER SETTINGS -----------------------------------
rds_in  <- "results/Step_01_import_mothur_to_microeco/phyloseq_raw.rds"
out_dir <- "results/Decontam_abort_version"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Required columns in sample_data(ps)
TYPE_COL   <- "SampleType"     # Fish / Blank / Water
CONC_COL   <- "DNA_ng_ul"      # numeric; may include 0
GROUP_SRC  <- "SampleTypes"    # Control / Dose_01% / Dose_2_2% (used for facet)
TANK_COL   <- "Tanks"          # optional (used only for ordering if you extend later)

# Taxonomy column present in tax_table(ps) from Step_01
# (You created it as: paste0(Genus_clean, "_", OTU))
FEATURE_COL <- "Genera_level_OTUs"

# decontam settings
NEG_LABELS <- c("Blank")       # treat tank water as environment, not negative control
DECONTAM_THRESHOLD <- 0.1
CONC_COL_USE <- "DNA_ng_ul_decontam"

# single-blank heuristic (supporting evidence; transparent)
MIN_BLANK_READS  <- 50
MIN_BLANK_REL    <- 0.005
MAX_FISH_PREV    <- 0.10
MAX_FISH_MAX_REL <- 0.0001

# Composition plots (microeco::trans_abund)
COMP_DIR <- file.path(out_dir, "Composition_BeforeAfter")
dir.create(COMP_DIR, recursive = TRUE, showWarnings = FALSE)

FACET_COL <- "Group"           # we will create mt$sample_table$Group from GROUP_SRC
TOP_PHYLUM <- 8
TOP_GLO    <- 25
OTHERS_COLOR <- "grey70"
XTEXT_KEEP <- FALSE
LEGEND_TEXT_ITALIC <- FALSE

# frequency example plots
N_PLOT_CONTAMS <- 6

stamp <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep="")

# =============================================================================
# 1) Load phyloseq from Step_01
# =============================================================================
stamp("Loading phyloseq: ", rds_in)
ps <- readRDS(rds_in)
if (!inherits(ps, "phyloseq")) stop("Input RDS is not a phyloseq object: ", rds_in)

sd <- as.data.frame(sample_data(ps), stringsAsFactors = FALSE)
for (req in c(TYPE_COL, CONC_COL, GROUP_SRC)) {
  if (!req %in% names(sd)) stop("sample_data(ps) missing required column: ", req)
}

sd[[TYPE_COL]] <- as.character(sd[[TYPE_COL]])
sd[[CONC_COL]] <- suppressWarnings(as.numeric(sd[[CONC_COL]]))
if (anyNA(sd[[CONC_COL]])) stop(CONC_COL, " contains NA/non-numeric values.")

# Create a positive concentration column for frequency method (zeros break it)
pos_vals <- sd[[CONC_COL]][sd[[CONC_COL]] > 0]
if (length(pos_vals) == 0) stop("All DNA concentrations are 0; frequency method cannot run.")
floor_val <- max(min(pos_vals, na.rm = TRUE) / 10, 0.001)

sd[[CONC_COL_USE]] <- sd[[CONC_COL]]
sd[[CONC_COL_USE]][sd[[CONC_COL_USE]] <= 0] <- floor_val

sd$is_neg <- sd[[TYPE_COL]] %in% NEG_LABELS
sample_data(ps) <- sample_data(sd)

fish_ids  <- rownames(sd)[sd[[TYPE_COL]] == "Fish"]
blank_ids <- rownames(sd)[sd[[TYPE_COL]] == "Blank"]
water_ids <- rownames(sd)[sd[[TYPE_COL]] == "Water"]

stamp("Samples: total=", nsamples(ps),
      " | fish=", length(fish_ids),
      " | blank=", length(blank_ids),
      " | water=", length(water_ids))
if (length(fish_ids) < 5) stop("Too few Fish samples for frequency method.")

# =============================================================================
# 2) QC plot: library sizes
# =============================================================================
sd_df <- data.frame(sample_data(ps), check.names = FALSE, stringsAsFactors = FALSE)
sd_df$Sample_ID <- rownames(sd_df)
libvec <- sample_sums(ps)
sd_df$LibrarySize <- as.numeric(libvec[sd_df$Sample_ID])
if (anyNA(sd_df$LibrarySize)) stop("LibrarySize NA (sample name mismatch).")

sd_df <- sd_df[order(sd_df$LibrarySize), , drop = FALSE]
sd_df$Index <- seq_len(nrow(sd_df))

p_ls <- ggplot(sd_df, aes(x = Index, y = LibrarySize, color = .data[[TYPE_COL]])) +
  geom_point() +
  scale_y_continuous(trans = "log10") +
  labs(title = "Library sizes (log10)", x = "Sample (sorted)", y = "Reads (log10)") +
  theme_bw()

ggsave(file.path(out_dir, "QC_library_sizes_log10.png"), p_ls, width = 9, height = 4.5, dpi = 300)

# =============================================================================
# 3) decontam FREQUENCY (Fish only)
# =============================================================================
ps_fish <- prune_samples(sample_data(ps)[[TYPE_COL]] == "Fish", ps)

stamp("Running decontam frequency on Fish only (threshold=", DECONTAM_THRESHOLD, ")")
contam_freq <- decontam::isContaminant(
  ps_fish,
  method    = "frequency",
  conc      = CONC_COL_USE,
  threshold = DECONTAM_THRESHOLD
)

# Map results back to full ps taxa
freq_flag <- rep(FALSE, ntaxa(ps)); names(freq_flag) <- taxa_names(ps)
freq_p    <- rep(NA_real_, ntaxa(ps)); names(freq_p) <- taxa_names(ps)
freq_flag[names(contam_freq$contaminant)] <- contam_freq$contaminant
freq_p[names(contam_freq$p)]              <- contam_freq$p

stamp("decontam frequency flagged: ", sum(freq_flag, na.rm = TRUE), " taxa")

# Example frequency plots (if any flagged)
contam_taxa <- names(which(freq_flag))
if (length(contam_taxa) > 0) {
  pick <- sample(contam_taxa, size = min(N_PLOT_CONTAMS, length(contam_taxa)))
  png(file.path(out_dir, "QC_frequency_plots_examples.png"), width = 2200, height = 1400, res = 200)
  print(decontam::plot_frequency(ps_fish, pick, conc = CONC_COL_USE) +
          xlab("DNA concentration (ng/µL; floored)") +
          ylab("Feature frequency"))
  dev.off()
}

# =============================================================================
# 4) Single-blank heuristic (supporting evidence)
# =============================================================================
otu <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu <- t(otu)

col_sums <- colSums(otu)
rel <- sweep(otu, 2, col_sums, "/")
rel[is.na(rel)] <- 0

blank_reads  <- rep(0, ntaxa(ps)); names(blank_reads) <- taxa_names(ps)
blank_rel   <- rep(0, ntaxa(ps)); names(blank_rel) <- taxa_names(ps)

fish_prev    <- rowMeans(otu[, fish_ids, drop = FALSE] > 0)
fish_max_rel <- apply(rel[, fish_ids, drop = FALSE], 1, max)

if (length(blank_ids) > 0) {
  blank_reads <- rowSums(otu[, blank_ids, drop = FALSE])
  blank_rel   <- rowSums(rel[, blank_ids, drop = FALSE]) / length(blank_ids)
}

blank_flag <- (blank_reads >= MIN_BLANK_READS) &
  (blank_rel   >= MIN_BLANK_REL) &
  (fish_prev   <= MAX_FISH_PREV) &
  (fish_max_rel <= MAX_FISH_MAX_REL)

stamp("blank red-flag flagged: ", sum(blank_flag, na.rm = TRUE), " taxa")

# =============================================================================
# 5) Define remove sets and export reports
# =============================================================================
remove_lenient <- freq_flag
remove_strict  <- freq_flag | blank_flag

tax <- as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)
tax$OTU <- rownames(tax)

report <- tibble(
  OTU = taxa_names(ps),
  decontam_freq_p = as.numeric(freq_p),
  decontam_freq_flag = as.logical(freq_flag),
  blank_reads = as.numeric(blank_reads),
  blank_rel = as.numeric(blank_rel),
  fish_prev = as.numeric(fish_prev),
  fish_max_rel = as.numeric(fish_max_rel),
  blank_flag = as.logical(blank_flag),
  remove_lenient = as.logical(remove_lenient),
  remove_strict  = as.logical(remove_strict)
) %>%
  left_join(tax, by = "OTU") %>%
  arrange(desc(remove_strict), desc(blank_reads), decontam_freq_p)

write_tsv(report, file.path(out_dir, "contaminant_report.tsv"))
write_tsv(filter(report, remove_lenient), file.path(out_dir, "removed_OTUs_lenient.tsv"))
write_tsv(filter(report, remove_strict),  file.path(out_dir, "removed_OTUs_strict.tsv"))

summary_tbl <- tibble(
  n_samples_total = nsamples(ps),
  n_samples_fish  = sum(sample_data(ps)[[TYPE_COL]] == "Fish"),
  n_samples_blank = sum(sample_data(ps)[[TYPE_COL]] == "Blank"),
  n_samples_water = sum(sample_data(ps)[[TYPE_COL]] == "Water"),
  n_otus_total    = ntaxa(ps),
  n_flag_freq     = sum(freq_flag, na.rm = TRUE),
  n_flag_blank    = sum(blank_flag, na.rm = TRUE),
  n_remove_lenient = sum(remove_lenient, na.rm = TRUE),
  n_remove_strict  = sum(remove_strict,  na.rm = TRUE),
  decontam_threshold = DECONTAM_THRESHOLD,
  conc_floor_val = floor_val
)
write_tsv(summary_tbl, file.path(out_dir, "decontam_summary.tsv"))

# =============================================================================
# 6) Prune + save cleaned phyloseq and microeco
# =============================================================================
ps_lenient <- prune_taxa(!remove_lenient, ps)
ps_strict  <- prune_taxa(!remove_strict,  ps)

saveRDS(ps,         file.path(out_dir, "ps_raw.rds"))
saveRDS(ps_lenient, file.path(out_dir, "ps_lenient.rds"))
saveRDS(ps_strict,  file.path(out_dir, "ps_strict.rds"))

# Convert ALL three to microeco for plotting convenience
mt_before  <- file2meco::phyloseq2meco(ps);         mt_before$tidy_dataset();  mt_before$cal_abund()
mt_lenient <- file2meco::phyloseq2meco(ps_lenient); mt_lenient$tidy_dataset(); mt_lenient$cal_abund()
mt_strict  <- file2meco::phyloseq2meco(ps_strict);  mt_strict$tidy_dataset();  mt_strict$cal_abund()

saveRDS(mt_before,  file.path(out_dir, "mt_before.rds"))
saveRDS(mt_lenient, file.path(out_dir, "mt_lenient.rds"))
saveRDS(mt_strict,  file.path(out_dir, "mt_strict.rds"))

# =============================================================================
# 7) Composition barplots using microeco::trans_abund (BEFORE vs AFTER)
# =============================================================================
ensure_group <- function(mt, group_src = GROUP_SRC, group_out = FACET_COL) {
  if (is.null(mt$sample_table) || nrow(mt$sample_table) == 0) stop("microtable sample_table missing/empty.")
  if (!group_src %in% colnames(mt$sample_table)) {
    stop("sample_table missing column '", group_src, "'. Available: ",
         paste(colnames(mt$sample_table), collapse = ", "))
  }
  mt$sample_table[[group_out]] <- as.character(mt$sample_table[[group_src]])
  mt
}

pick_taxrank <- function(mt, preferred, fallback) {
  tt <- mt$tax_table
  if (is.null(tt) || nrow(tt) == 0) stop("microtable tax_table missing/empty.")
  if (preferred %in% colnames(tt)) return(preferred)
  if (fallback %in% colnames(tt)) return(fallback)
  stop("Neither rank found: ", preferred, " / ", fallback, ". Available: ",
       paste(colnames(tt), collapse = ", "))
}

plot_and_save <- function(p, prefix, width = 13, height = 5) {
  ggsave(file.path(COMP_DIR, paste0(prefix, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(COMP_DIR, paste0(prefix, ".pdf")), p, width = width, height = height)
}

make_bar <- function(mt, taxrank, ntaxa, title_txt) {
  mt$tidy_dataset()
  mt$cal_abund()
  t1 <- trans_abund$new(dataset = mt, taxrank = taxrank, ntaxa = ntaxa)
  p <- t1$plot_bar(
    others_color = OTHERS_COLOR,
    facet = FACET_COL,
    xtext_keep = XTEXT_KEEP,
    legend_text_italic = LEGEND_TEXT_ITALIC
  ) +
    ggtitle(title_txt) +
    theme(plot.title = element_text(hjust = 0.5))
  p
}

stamp("Composition plots (microeco::trans_abund) ...")

mt_before  <- ensure_group(mt_before,  group_src = GROUP_SRC, group_out = FACET_COL)
mt_lenient <- ensure_group(mt_lenient, group_src = GROUP_SRC, group_out = FACET_COL)
mt_strict  <- ensure_group(mt_strict,  group_src = GROUP_SRC, group_out = FACET_COL)

# Phylum rank name may be "Phyla" in your parsed mothur taxonomy
RANK_PHYLUM <- pick_taxrank(mt_before, preferred = "Phylum", fallback = "Phyla")
RANK_GLO    <- pick_taxrank(mt_before, preferred = FEATURE_COL, fallback = FEATURE_COL)

# --- Top Phyla (8) ---
p_phy_before  <- make_bar(mt_before,  taxrank = RANK_PHYLUM, ntaxa = TOP_PHYLUM,
                          title_txt = paste0("Top ", TOP_PHYLUM, " Phyla — BEFORE decontam"))
p_phy_lenient <- make_bar(mt_lenient, taxrank = RANK_PHYLUM, ntaxa = TOP_PHYLUM,
                          title_txt = paste0("Top ", TOP_PHYLUM, " Phyla — AFTER decontam (lenient)"))
p_phy_strict  <- make_bar(mt_strict,  taxrank = RANK_PHYLUM, ntaxa = TOP_PHYLUM,
                          title_txt = paste0("Top ", TOP_PHYLUM, " Phyla — AFTER decontam (strict)"))

plot_and_save(p_phy_before,  "Bar_Phylum_top8_before",  width = 13, height = 5)
plot_and_save(p_phy_lenient, "Bar_Phylum_top8_lenient", width = 13, height = 5)
plot_and_save(p_phy_strict,  "Bar_Phylum_top8_strict",  width = 13, height = 5)

# --- Top Genera_level_OTUs (N) ---
p_glo_before  <- make_bar(mt_before,  taxrank = RANK_GLO, ntaxa = TOP_GLO,
                          title_txt = paste0("Top ", TOP_GLO, " Genera_level_OTUs — BEFORE decontam"))
p_glo_lenient <- make_bar(mt_lenient, taxrank = RANK_GLO, ntaxa = TOP_GLO,
                          title_txt = paste0("Top ", TOP_GLO, " Genera_level_OTUs — AFTER decontam (lenient)"))
p_glo_strict  <- make_bar(mt_strict,  taxrank = RANK_GLO, ntaxa = TOP_GLO,
                          title_txt = paste0("Top ", TOP_GLO, " Genera_level_OTUs — AFTER decontam (strict)"))

plot_and_save(p_glo_before,  paste0("Bar_GeneraLevelOTUs_top", TOP_GLO, "_before"),  width = 14, height = 5.5)
plot_and_save(p_glo_lenient, paste0("Bar_GeneraLevelOTUs_top", TOP_GLO, "_lenient"), width = 14, height = 5.5)
plot_and_save(p_glo_strict,  paste0("Bar_GeneraLevelOTUs_top", TOP_GLO, "_strict"),  width = 14, height = 5.5)

stamp("Done.")
stamp("Removed (lenient) = ", sum(remove_lenient, na.rm = TRUE), " OTUs")
stamp("Removed (strict)  = ", sum(remove_strict,  na.rm = TRUE), " OTUs")
stamp("Outputs: ", normalizePath(out_dir))