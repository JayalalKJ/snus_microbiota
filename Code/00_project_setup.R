# R
# Project setup script: create folders, ensure renv, install required packages,
# snapshot environment, and save session/package version info.

# Configuration ---------------------------------------------------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 1000)               # network timeout (seconds)
install_opts <- list(                 # optional per-package R CMD INSTALL flags
  # xml2 = "--no-lock"
)

dirs <- list(
  scripts = "scripts",
  data    = "data",
  results = "results",
  out     = file.path("results", "Step_00_setup")
)

cran_pkgs <- c(
  "microeco", "file2meco",
  "tidyverse", "readr", "dplyr", "tibble", "tidyr", "stringr", "purrr", "glue",
  "ggplot2", "ggpubr", "cowplot",
  "vegan", "ape", "sessioninfo",
  "svglite"
)

bioc_pkgs <- c("phyloseq", "decontam", "DESeq2", "edgeR")

github_repos <- c(
  "ChiLiubio/microeco",
  "gmteunisse/ggnested",
  "taowenmicro/ggClusterNet"
)

# Helpers ---------------------------------------------------------------------
dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    message("Created: ", path)
  }
}

ensure_package <- function(pkg, lib = .libPaths()[1]) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  invisible(TRUE)
}

install_missing_cran <- function(pkgs) {
  installed <- rownames(installed.packages())
  missing <- setdiff(pkgs, installed)
  if (length(missing)) {
    install.packages(missing, dependencies = TRUE)
  }
  invisible(missing)
}

install_missing_bioc <- function(pkgs) {
  installed <- rownames(installed.packages())
  missing <- setdiff(pkgs, installed)
  if (length(missing)) {
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  }
  invisible(missing)
}

install_missing_github <- function(repos) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  for (repo in repos) {
    pkg <- sub(".*/", "", repo)
    if (!requireNamespace(pkg, quietly = TRUE)) {
      remotes::install_github(repo, dependencies = TRUE, upgrade = "never")
    }
  }
  invisible(TRUE)
}

safe_snapshot <- function() {
  tryCatch({
    renv::snapshot(prompt = FALSE)
    TRUE
  }, error = function(e) {
    warning("renv snapshot failed: ", conditionMessage(e))
    FALSE
  })
}

# Run setup ------------------------------------------------------------------
# Create directories
lapply(dirs, dir_create)

# Ensure minimal base packages for installer tasks
ensure_package("utils")
ensure_package("remotes")    # used for GitHub installs
ensure_package("BiocManager")# used for Bioconductor installs

# renv bootstrap / restore ---------------------------------------------------
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
suppressPackageStartupMessages(library(renv))

# Use an explicit workflow: init if lock missing, else restore
if (!file.exists("renv.lock")) {
  message("Initializing renv (bare)...")
  renv::init(bare = TRUE)
} else {
  message("Activating renv and restoring from renv.lock...")
  renv::activate()
  renv::restore(prompt = FALSE)
}

# Optionally set install-time flags for R CMD INSTALL
if (length(install_opts)) {
  options(install.opts = install_opts)
}

# Install packages -----------------------------------------------------------
install_missing_cran(cran_pkgs)
install_missing_bioc(bioc_pkgs)
install_missing_github(github_repos)

# Snapshot environment and record versions -----------------------------------
safe_snapshot()

# Save session info and installed package versions for reproducibility
ensure_package("sessioninfo")
si <- capture.output(sessioninfo::session_info())
writeLines(si, con = file.path(dirs$out, "session_info.txt"))

pk <- as.data.frame(installed.packages()[, c("Package", "Version")], stringsAsFactors = FALSE)
pk <- pk[order(pk$Package), , drop = FALSE]
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
readr::write_tsv(pk, file.path(dirs$out, "package_versions.tsv"))

message("✅ Step 00 complete")
message("• renv.lock: ", if (file.exists("renv.lock")) normalizePath("renv.lock") else "not found")
message("• session:   ", normalizePath(file.path(dirs$out, "session_info.txt")))
message("• versions:  ", normalizePath(file.path(dirs$out, "package_versions.tsv")))
