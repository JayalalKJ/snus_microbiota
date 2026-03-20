#!/usr/bin/env Rscript
# =============================================================================
# Runs the full pipeline (00 → 05) in order, with logs + fail-fast behavior.
#
# Usage:
#   Rscript run_all_pipeline.R
#
# Assumes scripts are in: scripts/
# =============================================================================

options(stringsAsFactors = FALSE)

# ---- Config ----------------------------------------------------------------
scripts_dir <- "Code"
log_dir     <- file.path("results", "_pipeline_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

pipeline <- c(
  #"00_project_setup.R",
  #"01_import_mothur_to_microeco.R",
  #"02A_Blanks_visual_proofpack.R",
  #"02B_Negative_control_samples_Composition.R",
  #"02C_contaminant_removal_from_blanks.R",
  "03_rarefaction_Skin_sample_for_alpha.R",
  "04_alpha_diversity_Skin_samples.R"
  #"05_Beta_diversity_Skin_samples.R",
  #"07_composition_analysis_Skin_water_samples.R"
)

# ---- Helpers ---------------------------------------------------------------
stop_missing <- function(path) {
  if (!file.exists(path)) stop("❌ Missing file: ", path)
}

stamp <- function() format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

run_one <- function(script_name) {
  script_path <- file.path(scripts_dir, script_name)
  stop_missing(script_path)

  tag <- paste0(gsub("\\.R$", "", script_name), "__", stamp())
  out_log <- file.path(log_dir, paste0(tag, ".out.log"))
  err_log <- file.path(log_dir, paste0(tag, ".err.log"))

  message("\n──────────────────────────────────────────────────────────────")
  message("▶ Running: ", script_name)
  message("  out: ", out_log)
  message("  err: ", err_log)

  # Use the same R that is running this driver
  rbin <- file.path(R.home("bin"), "Rscript")
  if (!file.exists(rbin)) rbin <- "Rscript"

  exit_code <- system2(
    command = rbin,
    args    = c("--vanilla", script_path),
    stdout  = out_log,
    stderr  = err_log
  )

  if (!identical(exit_code, 0L)) {
    message("❌ Failed: ", script_name)
    message("   Check logs:\n   - ", out_log, "\n   - ", err_log)
    stop("Pipeline aborted (exit code = ", exit_code, ").")
  }

  message("✅ Done: ", script_name)
  invisible(list(script = script_name, out_log = out_log, err_log = err_log))
}

# ---- Run pipeline ----------------------------------------------------------
results <- lapply(pipeline, run_one)

# ---- Manifest --------------------------------------------------------------
manifest <- data.frame(
  step     = seq_along(pipeline),
  script   = vapply(results, `[[`, character(1), "script"),
  out_log  = vapply(results, `[[`, character(1), "out_log"),
  err_log  = vapply(results, `[[`, character(1), "err_log"),
  stringsAsFactors = FALSE
)

manifest_file <- file.path(log_dir, paste0("pipeline_manifest__", stamp(), ".tsv"))
write.table(manifest, file = manifest_file, sep = "\t", row.names = FALSE, quote = FALSE)

message("\n🎉 Pipeline complete.")
message("• Manifest: ", normalizePath(manifest_file))
message("• Logs dir: ", normalizePath(log_dir))
