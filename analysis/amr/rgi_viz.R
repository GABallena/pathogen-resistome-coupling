#!/usr/bin/env Rscript
# Portfolio-safe script (paths/identifiers generalized; inputs not included).

# ------- Working directory handling -------
# When run via `Rscript`, set working directory to the script's directory.
try({
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if(length(file_arg) == 1){
    script_path <- normalizePath(sub("^--file=", "", file_arg))
    setwd(dirname(script_path))
  }
}, silent = TRUE)
cat("Working directory:", getwd(), "\n")


# rgi_viz.R — Visualization toolkit for CARD RGI outputs (contig-level)
#
# Features
# - Auto-discovers RGI result files in a directory (multiple samples)
# - Robust parser for common RGI tabular outputs; basic JSON support
# - Clean aggregation by Drug Class, Mechanism, and Gene (ARO/Best_Hit_ARO)
# - Publication-ready plots:
#     1) Stacked bar by Drug Class per sample
#     2) Top-N ARG genes across samples
#     3) Heatmap of gene counts (presence/abundance) across samples
# - Exports CSV summaries and high-res PNG/PDF figures
#
# Usage
#   Rscript rgi_viz.R \
#     --input-dir rgi_outputs \
#     --output-dir rgi_viz \
#     --pattern "(.*)\\.(txt|tsv|csv|json)$" \
#     --sample-regex "^(.*?)(?:[._]rgi.*)?$" \
#     --top-n 50 \
#     --ontology-dir card-ontology \
#     --model-filter "protein homolog model" \
#     --identity-min 90 --coverage-min 50 --bitscore-min 400 --cutoff-levels "Perfect,Strict" \
#     --count-mode hit
#
# Notes
# - Expects RGI outputs from `rgi main` (tab-delimited) with columns similar to:
#   Best_Hit_ARO, ARO, Drug Class, Resistance Mechanism, AMR Gene Family, Best_Identities, Contig, ORF_ID, Percentage Length of Reference Sequence
# - If some columns are missing, the script will skip dependent plots and report a warning.

suppressWarnings({})
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(rlang)
  library(forcats)

# ---- Config ----
RGI_TSV <- Sys.getenv("RGI_TSV", "results/tables/rgi_summary.tsv")
CONTIG_META_TSV <- Sys.getenv("CONTIG_META_TSV", "results/tables/contig_metadata.tsv")
SAMPLE_META_TSV <- Sys.getenv("SAMPLE_META_TSV", "metadata/sample_metadata.tsv")
OUT_DIR <- Sys.getenv("OUT_DIR", "results/figures/amr")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

})

# -- Helpers -------------------------------------------------------------------

require_pkg <- function(pkg, install_if_missing = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (isTRUE(install_if_missing)) {
      message(sprintf("Installing missing package: %s", pkg))
      install.packages(pkg, repos = "https://cloud.r-project.org")
    } else {
      warning(sprintf("Package '%s' is not installed. Some features may be limited.", pkg))
      return(FALSE)
    }
  }
  TRUE
}

make_clean <- function(x) {
  # Lowercase, replace non-alnum with underscore, squeeze underscores
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

first_present <- function(nms, candidates) {
  for (c in candidates) {
    if (c %in% nms) return(c)
  }
  NA_character_
}

infer_sample_name <- function(filepath, sample_regex) {
  base <- basename(filepath)
  base <- sub("\\.(txt|tsv|csv|json|gz)$", "", base, ignore.case = TRUE)
  default_pat <- "^(.*?)(?:[._]rgi.*)?$"
  pat <- if (is.character(sample_regex) && length(sample_regex) == 1L && !is.na(sample_regex) && nzchar(sample_regex)) sample_regex else default_pat
  m <- tryCatch(regexec(pat, base, perl = TRUE), error = function(e) regexec(default_pat, base, perl = TRUE))
  matched <- regmatches(base, m)
  if (length(matched) == 1 && length(matched[[1]]) >= 2) {
    return(matched[[1]][2])
  }
  # Fallback: strip common trailing tokens
  sub("[._]rgi.*$", "", base, perl = TRUE)
}

read_rgi_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "json") {
    if (!require_pkg("jsonlite", install_if_missing = FALSE)) {
      stop("jsonlite not available to parse JSON: ", path)
    }
    dat <- jsonlite::fromJSON(path, flatten = TRUE)
    # Try to coerce to tibble sensibly
    if (is.data.frame(dat)) {
      as_tibble(dat)
    } else if (is.list(dat) && length(dat) > 0 && is.data.frame(dat[[1]])) {
      # Some RGI JSONs are lists of records
      bind_rows(lapply(dat, as_tibble))
    } else {
      stop("Unsupported JSON structure in ", path)
    }
  } else if (ext %in% c("txt", "tsv")) {
    suppressWarnings(readr::read_tsv(path, guess_max = 100000, show_col_types = FALSE))
  } else if (ext == "csv") {
    suppressWarnings(readr::read_csv(path, guess_max = 100000, show_col_types = FALSE))
  } else {
    # Default try TSV
    suppressWarnings(readr::read_tsv(path, guess_max = 100000, show_col_types = FALSE))
  }
}

standardize_rgi_cols <- function(df) {
  if (nrow(df) == 0) return(tibble())
  original_names <- names(df)
  clean_names <- make_clean(original_names)
  names(df) <- clean_names

  nm <- names(df)
  col_map <- list(
    contig = first_present(nm, c("contig", "contig_id")),
    orf_id = first_present(nm, c("orf_id", "orf")),
    best_hit_aro = first_present(nm, c("best_hit_aro", "best_hit", "besthitaro")),
    aro = first_present(nm, c("aro")),
  model_type = first_present(nm, c("model_type", "model", "modeltype")),
  cut_off = first_present(nm, c("cut_off", "cutoff")),
  pass_bitscore = first_present(nm, c("pass_bitscore", "passbitscore", "pass_bit_score")),
    drug_class = first_present(nm, c("drug_class", "drugclass")),
    resistance_mechanism = first_present(nm, c("resistance_mechanism", "mechanism", "resistance")),
    amr_gene_family = first_present(nm, c("amr_gene_family", "gene_family", "amr_family")),
    best_identities = first_present(nm, c("best_identities", "identity", "percent_identity", "perc_identity")),
    percent_length_ref = first_present(nm, c("percentage_length_of_reference_sequence", "percent_length_of_reference_sequence", "perc_length_ref", "coverage", "alignment_coverage"))
  )

  # Canonical selection
  to_char <- function(x) { if (is.null(x)) NULL else as.character(x) }
  out <- tibble::tibble(
    contig = if (!is.na(col_map$contig)) df[[col_map$contig]] else NA_character_,
    orf_id = if (!is.na(col_map$orf_id)) df[[col_map$orf_id]] else NA_character_,
    gene = dplyr::coalesce(
      to_char(if (!is.na(col_map$best_hit_aro)) df[[col_map$best_hit_aro]] else NULL),
      to_char(if (!is.na(col_map$aro)) df[[col_map$aro]] else NULL),
      NA_character_
    ),
  model_type = if (!is.na(col_map$model_type)) df[[col_map$model_type]] else NA_character_,
  cut_off = if (!is.na(col_map$cut_off)) df[[col_map$cut_off]] else NA_character_,
  pass_bitscore = if (!is.na(col_map$pass_bitscore)) df[[col_map$pass_bitscore]] else NA,
  aro_accession = to_char(if (!is.na(col_map$aro)) df[[col_map$aro]] else NA_character_),
    drug_class = if (!is.na(col_map$drug_class)) df[[col_map$drug_class]] else NA_character_,
    mechanism = if (!is.na(col_map$resistance_mechanism)) df[[col_map$resistance_mechanism]] else NA_character_,
    gene_family = if (!is.na(col_map$amr_gene_family)) df[[col_map$amr_gene_family]] else NA_character_,
    identity = suppressWarnings(as.numeric(if (!is.na(col_map$best_identities)) df[[col_map$best_identities]] else NA_real_)),
    coverage_ref_pct = suppressWarnings(as.numeric(if (!is.na(col_map$percent_length_ref)) df[[col_map$percent_length_ref]] else NA_real_))
  )

  # Normalize text fields, split multi-valued columns if necessary (take first)
  out <- out %>%
    mutate(across(
  c(.data$gene, .data$model_type, .data$cut_off, .data$drug_class, .data$mechanism, .data$gene_family),
      ~ ifelse(is.na(.x) | .x == "", NA_character_, .x)
    ))

  # RGI often uses '; ' separated for Drug Class etc.; take first for plotting or keep full
  split_first <- function(v) {
    ifelse(
      is.na(v) | v == "",
      NA_character_,
      vapply(strsplit(v, ";\\s*"), function(x) if (length(x)) x[1] else NA_character_, character(1))
    )
  }
  out <- out %>% mutate(
    drug_class_primary = split_first(.data$drug_class),
    mechanism_primary = split_first(.data$mechanism),
    gene_family_primary = split_first(.data$gene_family)
  )

  out
}

save_plot <- function(p, file_base, width = 10, height = 6, dpi = 300, bg = "white") {
  ggsave(paste0(file_base, ".png"), p, width = width, height = height, dpi = dpi, bg = bg, limitsize = FALSE)
  pdf_dev <- if (isTRUE(capabilities("cairo"))) cairo_pdf else "pdf"
  ggsave(paste0(file_base, ".pdf"), p, width = width, height = height, device = pdf_dev, bg = bg, limitsize = FALSE)
}

# -- Palette utilities ---------------------------------------------------------

get_base_palette <- function(name = "economist") {
  name <- tolower(name)
  if (name == "ft" || name == "financial_times") {
    c(
      "#173E6A", "#4899B8", # Blues
      "#B0133D", "#D65A76", "#E393A4", # Reds/Pinks
      "#9DAA62", "#9D8D7B", # Green/Neutral
      "#C9B9A9" # Neutral
    )
  } else {
    # Economist default
    c(
      "#2F6FA1", "#41A9B8", "#3BB3DC", # Blues
      "#E2931B", "#C7B434", # Orange/Yellow
      "#E5251F", "#E9541B", # Reds
      "#907070", "#C4B7A4" # Neutral
    )
  }
}

get_palette <- function(name = "economist", n = 10) {
  # For small n, use curated palettes; for larger n, extend deterministically without viridis
  base <- get_base_palette(name)
  if (n <= length(base)) return(base[seq_len(n)])
  # Extend by repeating base and applying slight lightness adjustments per cycle
  reps <- ceiling(n / length(base))
  cols <- rep(base, reps)[seq_len(n)]
  # Apply subtle lightness shifts to each cycle to improve distinction
  if (requireNamespace("colorspace", quietly = TRUE)) {
    cyc <- rep(seq_len(reps), each = length(base))[seq_len(n)]
    shift <- (cyc - 1) * 0.07
    cols <- colorspace::lighten(cols, amount = pmin(shift, 0.35))
  }
  cols
}

wrap_lab <- function(x, width = 18) {
  stringr::str_replace_all(stringr::str_wrap(x, width = width), "\n", "\n")
}

# Heatmap two-color mapping (Economist dark blue ↔ dark red)
# Low (blue):  #1B5F85  | High (red): #A61C2E | Mid (light neutral): #F1EFEA
get_heatmap_low_high <- function(name = "economist_div") {
  low <- get_opt(args, "heatmap-low"); mid <- get_opt(args, "heatmap-mid"); high <- get_opt(args, "heatmap-high")
  def_low <- "#1B5F85"; def_mid <- "#F1EFEA"; def_high <- "#A61C2E"
  list(low = ifelse(nchar(low) > 0, low, def_low), high = ifelse(nchar(high) > 0, high, def_high), mid = ifelse(nchar(mid) > 0, mid, def_mid))
}

heatmap_scale_layer <- function(norm_mode, pal_name, lims = NULL) {
  cols <- get_heatmap_low_high(pal_name)
  nm <- tolower(norm_mode)
  legend_title <- switch(nm,
                         "row_z" = "Row z-score",
                         "row_prop" = "Row proportion",
                         "col_prop" = "Column proportion",
                         "log1p" = "log1p(count)",
                         "Count")
  if (nm == "row_z") {
    return(scale_fill_gradient2(low = cols$low, mid = cols$mid, high = cols$high, limits = lims, oob = scales::squish, midpoint = 0, name = legend_title))
  }
  scale_fill_gradient(low = cols$low, high = cols$high, limits = lims, oob = scales::squish, name = legend_title)
}

normalize_matrix <- function(mat, mode = "none") {
  m <- as.matrix(mat)
  mode <- tolower(mode)
  if (mode == "row_z") {
    m <- t(scale(t(m)))
    m[is.nan(m)] <- 0
  } else if (mode == "row_prop") {
    rs <- rowSums(m)
    rs[rs == 0] <- 1
    m <- m / rs
  } else if (mode == "col_prop") {
    cs <- colSums(m)
    cs[cs == 0] <- 1
    m <- sweep(m, 2, cs, "/")
  } else if (mode == "log1p") {
    m <- log1p(m)
  }
  m
}

dist_matrix <- function(m, method = "euclidean") {
  method <- tolower(method)
  if (method %in% c("euclidean", "manhattan")) return(dist(m, method = method))
  if (method %in% c("pearson", "spearman")) {
    cor_method <- method
    cmat <- suppressWarnings(stats::cor(t(m), method = cor_method, use = "pairwise.complete.obs"))
    return(as.dist(1 - cmat))
  }
  dist(m)
}

cluster_order <- function(m, distance = "euclidean", linkage = "complete", axis = "rows") {
  if (axis == "cols") m <- t(m)
  d <- dist_matrix(m, distance)
  hc <- hclust(d, method = linkage)
  hc$order
}

# -- CLI -----------------------------------------------------------------------

option_list <- list(
  make_option(c("--config"), type = "character", default = "", help = "Path to TOML config file (CLI overrides config)."),
  make_option(c("-i", "--input-dir"), type = "character", default = "rgi_outputs", help = "Directory containing RGI output files."),
  make_option(c("-o", "--output-dir"), type = "character", default = "rgi_viz", help = "Directory to write figures and summaries."),
  make_option(c("-p", "--pattern"), type = "character", default = "(.*)\\.(txt|tsv|csv|json)$", help = "Regex for filenames to include."),
  make_option(c("-r", "--sample-regex"), type = "character", default = "^(.*?)(?:[._]rgi.*)?$", help = "Regex with one capture group to extract sample name from filename (default: basename minus trailing rgi tokens)."),
  make_option(c("-n", "--top-n"), type = "integer", default = 50L, help = "Top N genes for heatmap and top-genes plot."),
  make_option(c("--install-missing"), type = "logical", default = FALSE, help = "Attempt to install missing packages automatically."),
  make_option(c("--ontology-dir"), type = "character", default = "card-ontology", help = "Directory containing CARD ontology files (expects aro.tsv)."),
  make_option(c("--model-filter"), type = "character", default = "protein homolog model", help = "Only include rows with this RGI model_type. Set empty string to include all."),
  make_option(c("--identity-min"), type = "double", default = NA, help = "Minimum percent identity (e.g., 90)."),
  make_option(c("--coverage-min"), type = "double", default = NA, help = "Minimum percent reference length coverage (e.g., 50)."),
  make_option(c("--bitscore-min"), type = "double", default = NA, help = "Minimum pass bitscore (numeric)."),
  make_option(c("--cutoff-levels"), type = "character", default = "", help = "Comma-separated list of allowed RGI cut_off levels (e.g., 'Perfect,Strict'). Empty for no filter."),
  make_option(c("--count-mode"), type = "character", default = "hit", help = "Counting mode: hit|orf|contig|presence."),
  make_option(c("--drop-ontology-categories"), type = "logical", default = TRUE, help = "Drop non-gene ontology rows (family/group/mechanism/drug/regulator).")
  , make_option(c("--palette"), type = "character", default = "economist", help = "Palette to use: 'economist' or 'ft' [default: %default]")
  , make_option(c("--label-wrap"), type = "integer", default = 18, help = "Wrap width for long labels [default: %default]")
  , make_option(c("--lump-top-n"), type = "integer", default = 12, help = "For categorical stacks, keep top-N and lump the rest into 'Other' [default: %default]")
  , make_option(c("--metadata-file"), type = "character", default = "", help = "Path to metadata file (csv/tsv) with sample annotations.")
  , make_option(c("--metadata-sample-col"), type = "character", default = "", help = "Column name in metadata corresponding to sample IDs.")
  , make_option(c("--metadata-type-col"), type = "character", default = "", help = "Column name in metadata corresponding to sample type/group.")
  , make_option(c("--drop-singletons"), type = "logical", default = FALSE, help = "Drop genes with only one copy per sample from gene-based summaries/plots [default: %default]")
  , make_option(c("--heatmap-normalize"), type = "character", default = "none", help = "Normalization for heatmaps: none|row_z|row_prop|col_prop|log1p [default: %default]")
  , make_option(c("--heatmap-cluster"), type = "logical", default = TRUE, help = "Cluster rows/columns in heatmaps [default: %default]")
  , make_option(c("--heatmap-low"), type = "character", default = "", help = "Override low color for heatmaps (hex, e.g., #1B4F72)")
  , make_option(c("--heatmap-mid"), type = "character", default = "", help = "Override mid color for diverging heatmaps (hex)")
  , make_option(c("--heatmap-high"), type = "character", default = "", help = "Override high color for heatmaps (hex, e.g., #A61C2E)")
  , make_option(c("--cluster-distance"), type = "character", default = "euclidean", help = "Distance: euclidean|manhattan|pearson|spearman [default: %default]")
  , make_option(c("--cluster-linkage"), type = "character", default = "complete", help = "Linkage: complete|average|ward.D2|single [default: %default]")
)

args <- parse_args(OptionParser(option_list = option_list))
# -- Config loading (TOML) ----------------------------------------------------

simple_parse_toml <- function(path) {
  # Minimal TOML subset: [section], key = value (strings, numbers, booleans)
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path, warn = FALSE)
  cfg <- list()
  cur <- "root"
  cfg[[cur]] <- list()
  for (ln in lines) {
    s <- trimws(sub("#.*$", "", ln))
    if (s == "") next
    if (grepl("^\\[.*\\]$", s)) {
      cur <- sub("^\\[(.*)\\]$", "\\1", s)
      cfg[[cur]] <- cfg[[cur]] %||% list()
    } else if (grepl("=", s)) {
      kv <- strsplit(s, "=", fixed = TRUE)[[1]]
      key <- trimws(kv[1])
      val <- trimws(paste(kv[-1], collapse = "="))
      # strip quotes if present
      if (grepl('^".*"$', val) || grepl("^'.*'$", val)) {
        val <- sub('^"', '', sub('"$', '', val))
        val <- sub("^'", '', sub("'$", '', val))
        # unescape common sequences: convert double backslashes to single
        val <- gsub("\\\\\\\\", "\\\\", val)
      } else if (grepl("^(true|false)$", tolower(val))) {
        val <- tolower(val) == "true"
      } else if (suppressWarnings(!is.na(as.numeric(val)))) {
        val <- as.numeric(val)
      }
      cfg[[cur]][[key]] <- val
    }
  }
  cfg
}

`%||%` <- function(a,b) if (!is.null(a)) a else b

get_opt2 <- function(opts, base) {
  cand <- c(base, gsub("-", "_", base, fixed = TRUE), gsub("_", "-", base, fixed = TRUE))
  for (nm in cand) {
    if (!is.null(opts[[nm]])) return(opts[[nm]])
  }
  NULL
}

defaults_map <- list(
  `input-dir` = "rgi_outputs", `output-dir` = "rgi_viz", pattern = "(.*)\\.(txt|tsv|csv|json)$",
  `sample-regex` = "^(.*?)(?:[._]rgi.*)?$", `top-n` = 50L, `install-missing` = FALSE,
  `ontology-dir` = "card-ontology", `model-filter` = "protein homolog model",
  `identity-min` = NA, `coverage-min` = NA, `bitscore-min` = NA,
  `cutoff-levels` = "", `count-mode` = "hit", `drop-ontology-categories` = TRUE,
  palette = "economist", `label-wrap` = 18L, `lump-top-n` = 12L,
  `metadata-file` = "", `metadata-sample-col` = "", `metadata-type-col` = "",
  `heatmap-normalize` = "none", `heatmap-cluster` = TRUE,
  `cluster-distance` = "euclidean", `cluster-linkage` = "complete",
  `drop-singletons` = FALSE
)

set_opt <- function(lst, base, value) {
  lst[[base]] <- value
  lst[[gsub("-", "_", base, fixed = TRUE)]] <- value
  lst[[gsub("_", "-", base, fixed = TRUE)]] <- value
  lst
}

read_cfg <- function(path) {
  cfg <- NULL
  if (nzchar(path) && file.exists(path)) {
    # Try RcppTOML if available, else fallback
    if (requireNamespace("RcppTOML", quietly = TRUE)) {
      cfg <- tryCatch(RcppTOML::parseTOML(path), error = function(e) NULL)
    }
    if (is.null(cfg)) cfg <- simple_parse_toml(path)
  }
  cfg
}

cfg <- read_cfg(get_opt2(args, "config"))
config_note <- if (!is.null(cfg)) normalizePath(get_opt2(args, "config")) else "none"

pick_cfg <- function(section, key, cli_name) {
  cli_val <- get_opt2(args, cli_name)
  def_val <- defaults_map[[cli_name]] %||% NULL
  has_cli <- !is.null(cli_val) && !identical(cli_val, def_val)
  if (has_cli) return(cli_val)
  if (!is.null(cfg) && !is.null(cfg[[section]]) && !is.null(cfg[[section]][[key]])) return(cfg[[section]][[key]])
  cli_val
}

# Apply config to args (CLI overrides config)
args <- set_opt(args, "input-dir", pick_cfg("input", "input_dir", "input-dir"))
args <- set_opt(args, "pattern", pick_cfg("input", "pattern", "pattern"))
args <- set_opt(args, "sample-regex", pick_cfg("input", "sample_regex", "sample-regex"))
args <- set_opt(args, "output-dir", pick_cfg("output", "output_dir", "output-dir"))

args <- set_opt(args, "model-filter", pick_cfg("filters", "model_filter", "model-filter"))
args <- set_opt(args, "identity-min", pick_cfg("filters", "identity_min", "identity-min"))
args <- set_opt(args, "coverage-min", pick_cfg("filters", "coverage_min", "coverage-min"))
args <- set_opt(args, "bitscore-min", pick_cfg("filters", "bitscore_min", "bitscore-min"))
args <- set_opt(args, "cutoff-levels", pick_cfg("filters", "cutoff_levels", "cutoff-levels"))
args <- set_opt(args, "count-mode", pick_cfg("filters", "count_mode", "count-mode"))
args <- set_opt(args, "drop-ontology-categories", pick_cfg("filters", "drop_ontology_categories", "drop-ontology-categories"))

args <- set_opt(args, "top-n", pick_cfg("top", "top_n", "top-n"))

args <- set_opt(args, "palette", pick_cfg("palette", "name", "palette"))
args <- set_opt(args, "label-wrap", pick_cfg("palette", "label_wrap", "label-wrap"))
args <- set_opt(args, "lump-top-n", pick_cfg("palette", "lump_top_n", "lump-top-n"))

args <- set_opt(args, "metadata-file", pick_cfg("metadata", "file", "metadata-file"))
args <- set_opt(args, "metadata-sample-col", pick_cfg("metadata", "sample_col", "metadata-sample-col"))
args <- set_opt(args, "metadata-type-col", pick_cfg("metadata", "type_col", "metadata-type-col"))

args <- set_opt(args, "heatmap-normalize", pick_cfg("heatmap", "normalize", "heatmap-normalize"))
args <- set_opt(args, "heatmap-cluster", pick_cfg("heatmap", "cluster", "heatmap-cluster"))
args <- set_opt(args, "heatmap-low", pick_cfg("heatmap", "low", "heatmap-low"))
args <- set_opt(args, "heatmap-mid", pick_cfg("heatmap", "mid", "heatmap-mid"))
args <- set_opt(args, "heatmap-high", pick_cfg("heatmap", "high", "heatmap-high"))
args <- set_opt(args, "cluster-distance", pick_cfg("heatmap", "cluster_distance", "cluster-distance"))
args <- set_opt(args, "cluster-linkage", pick_cfg("heatmap", "cluster_linkage", "cluster-linkage"))

args <- set_opt(args, "drop-singletons", pick_cfg("genes", "drop_singletons", "drop-singletons"))


normalize_names <- function(x) {
  x <- gsub("[\u00A0\t\n\r]+", " ", x, perl = TRUE)
  x <- gsub("[ .]+", "_", trimws(tolower(x)))
  x
}

sanitize_key <- function(x) {
  # alphanumeric lowercase key without separators
  tolower(gsub("[^A-Za-z0-9]+", "", x))
}

read_metadata <- function(file_path) {
  if (!nzchar(file_path) || is.na(file_path) || !file.exists(file_path)) return(NULL)
  ext <- tolower(tools::file_ext(file_path))
  df <- tryCatch({
    if (ext %in% c("tsv", "txt")) readr::read_tsv(file_path, show_col_types = FALSE)
    else readr::read_csv(file_path, show_col_types = FALSE)
  }, error = function(e) NULL)
  if (is.null(df)) return(NULL)
  names(df) <- normalize_names(names(df))
  df
}

autodetect_metadata <- function(input_dir, explicit_path = "") {
  if (nzchar(explicit_path) && file.exists(explicit_path)) return(explicit_path)
  base_dir <- normalizePath(file.path(input_dir, ".."))
  cand <- list.files(base_dir, pattern = "(?i)(metadata|masterlist|sample).*\\.(csv|tsv|txt)$", full.names = TRUE)
  if (length(cand) > 0) return(cand[[1]])
  ""
}

pick_col <- function(nms, preferred = "", patterns = c()) {
  if (nzchar(preferred) && preferred %in% nms) return(preferred)
  for (pat in patterns) {
    hit <- nms[grepl(pat, nms, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[[1]])
  }
  ""
}
# Helper to fetch options whether defined with hyphen or underscore
get_opt <- function(opts, base) {
  cand <- c(base, gsub("-", "_", base, fixed = TRUE), gsub("_", "-", base, fixed = TRUE))
  for (nm in cand) {
    if (!is.null(opts[[nm]])) return(opts[[nm]])
  }
  NULL
}

# Sanitize path args
input_dir_raw <- get_opt(args, "input-dir")
output_dir_raw <- get_opt(args, "output-dir")
input_dir <- if (!is.null(input_dir_raw) && length(input_dir_raw) >= 1) path.expand(as.character(input_dir_raw)[1]) else NA_character_
output_dir <- if (!is.null(output_dir_raw) && length(output_dir_raw) >= 1) path.expand(as.character(output_dir_raw)[1]) else NA_character_

if (is.na(input_dir) || length(input_dir) != 1 || !is.character(input_dir) || !nzchar(input_dir)) {
  stop("Invalid --input-dir value")
}
if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir)
}
if (is.na(output_dir) || length(output_dir) != 1 || !is.character(output_dir) || !nzchar(output_dir)) {
  stop("Invalid --output-dir value")
}
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Optional: install extra packages for nicer text rendering
if (isTRUE(args$install_missing)) {
  require_pkg("systemfonts", TRUE)
  require_pkg("textshaping", TRUE)
  require_pkg("ragg", TRUE)
}

# -- Discover and ingest -------------------------------------------------------

pattern_str <- args$pattern
if (grepl("\\\\\\\\", pattern_str)) pattern_str <- gsub("\\\\\\\\", "\\\\", pattern_str)
message("Using input dir=", input_dir, ", pattern=", pattern_str)
files <- list.files(input_dir, full.names = TRUE)
bn <- basename(files)
message("Found ", length(files), " files before regex. Examples: ", paste(utils::head(bn, 5), collapse=", "))
hits <- grepl(pattern_str, bn, ignore.case = TRUE, perl = TRUE)
message("Regex hit count: ", sum(hits))
files <- files[hits]

if (length(files) == 0) stop("No files matched in ", args$input_dir, " with pattern ", args$pattern)

message("Found ", length(files), " file(s). Parsing…")

all_df <- purrr::map_dfr(files, function(f) {
  df <- tryCatch(read_rgi_file(f), error = function(e) { warning(e$message); return(NULL) })
  if (is.null(df) || nrow(df) == 0) return(NULL)
  sdf <- standardize_rgi_cols(df)
  if (nrow(sdf) == 0) return(NULL)
  sample <- infer_sample_name(f, args$sample_regex)
  sdf$sample <- sample
  sdf$source_file <- basename(f)
  sdf
})

if (nrow(all_df) == 0) stop("No valid RGI records parsed.")

# Optional model filter (default: protein homolog model)
model_filter <- get_opt(args, "model-filter")
if (is.null(model_filter)) model_filter <- "protein homolog model"
if (!is.null(model_filter) && nzchar(model_filter)) {
  all_df <- all_df %>% filter(is.na(.data$model_type) | .data$model_type == model_filter)
}

# rgi main filters
identity_min <- suppressWarnings(as.numeric(get_opt(args, "identity-min")))
coverage_min <- suppressWarnings(as.numeric(get_opt(args, "coverage-min")))
bitscore_min <- suppressWarnings(as.numeric(get_opt(args, "bitscore-min")))
cutoff_levels <- get_opt(args, "cutoff-levels")
count_mode <- get_opt(args, "count-mode")
if (is.null(count_mode) || !nzchar(count_mode)) count_mode <- "hit"  # hit|orf|contig|presence

# Coerce numerics robustly just before filtering
all_df <- all_df %>% mutate(
  identity = suppressWarnings(as.numeric(.data$identity)),
  coverage_ref_pct = suppressWarnings(as.numeric(.data$coverage_ref_pct)),
  pass_bitscore_num = suppressWarnings(as.numeric(.data$pass_bitscore))
)

if (!is.na(identity_min)) all_df <- all_df %>% filter(is.na(.data$identity) | .data$identity >= identity_min)
if (!is.na(coverage_min)) all_df <- all_df %>% filter(is.na(.data$coverage_ref_pct) | .data$coverage_ref_pct >= coverage_min)
if (!is.na(bitscore_min)) all_df <- all_df %>% filter(is.na(.data$pass_bitscore_num) | .data$pass_bitscore_num >= bitscore_min)
if (!is.null(cutoff_levels) && nzchar(cutoff_levels)) {
  # Keep rows where cut_off is in provided levels (case-insensitive match)
  lev <- trimws(strsplit(cutoff_levels, ",")[[1]])
  lev <- tolower(lev)
  all_df <- all_df %>% mutate(cut_off_low = tolower(.data$cut_off)) %>% filter(is.na(.data$cut_off_low) | .data$cut_off_low %in% lev) %>% select(-cut_off_low)
}

# -- Metadata join (sample type) ---------------------------------------------
meta_path <- autodetect_metadata(input_dir, get_opt(args, "metadata-file"))
sample_type_info <- NULL
meta_note <- "none"
if (nzchar(meta_path)) {
  md <- read_metadata(meta_path)
  if (!is.null(md)) {
    nms <- names(md)
    sample_col <- get_opt(args, "metadata-sample-col")
    type_col <- get_opt(args, "metadata-type-col")
    # If user provided names, try to normalize and match against normalized headers
    if (nzchar(sample_col)) {
      sc_norm <- normalize_names(sample_col)
      if (sc_norm %in% nms) {
        sample_col <- sc_norm
      } else {
        sc_san <- sanitize_key(sample_col)
        cand <- nms[sanitize_key(nms) == sc_san]
        if (length(cand) >= 1) sample_col <- cand[[1]]
      }
    }
    if (nzchar(type_col)) {
      tc_norm <- normalize_names(type_col)
      if (tc_norm %in% nms) {
        type_col <- tc_norm
      } else {
        tc_san <- sanitize_key(type_col)
        cand <- nms[sanitize_key(nms) == tc_san]
        if (length(cand) >= 1) type_col <- cand[[1]]
      }
    }
    # If still not set or not present, fall back to autodetection patterns
    if (!nzchar(sample_col) || !(sample_col %in% nms)) sample_col <- pick_col(nms, patterns = c("^sample$", "sample_id", "^id$", "^sampleid$", "^sample_name$", "^name$", "sample_code", "samplecode"))
    if (!nzchar(type_col) || !(type_col %in% nms)) type_col <- pick_col(nms, patterns = c("sample_type", "^type$", "group", "category", "site", "source", "matrix", "environment"))
    if (nzchar(sample_col) && nzchar(type_col) && sample_col %in% nms && type_col %in% nms) {
      base_key <- function(x) {
        x <- as.character(x)
        x <- sub("_.*$", "", x) # drop suffix after underscore
        x <- gsub("[^A-Za-z0-9]+", "", x)
        tolower(x)
      }
      sample_type_info <- md %>%
        mutate(
          .sample_key = sanitize_key(.data[[sample_col]]),
          .sample_base = base_key(.data[[sample_col]]),
          sample_type = as.character(.data[[type_col]])
        ) %>%
        select(.sample_key, .sample_base, sample_type) %>% distinct()

      all_df <- all_df %>% mutate(.sample_key = sanitize_key(.data$sample), .sample_base = base_key(.data$sample)) %>%
        left_join(sample_type_info %>% select(.sample_key, sample_type), by = ".sample_key")
      # Fallback by base key for those still NA
      left_base <- all_df %>% filter(is.na(sample_type)) %>% select(.sample_base) %>% distinct()
      if (nrow(left_base) > 0) {
        all_df <- all_df %>%
          left_join(sample_type_info %>% select(.sample_base, sample_type_base = sample_type), by = ".sample_base") %>%
          mutate(sample_type = coalesce(sample_type, sample_type_base)) %>%
          select(-sample_type_base)
      }
  all_df <- all_df %>% select(-.sample_key, -.sample_base)
      matched <- sum(!is.na(all_df$sample_type))
      total_samp <- n_distinct(all_df$sample)
      meta_note <- paste0(basename(meta_path), " [sample_col=", sample_col, ", type_col=", type_col, ", matched_types=", matched, "/", total_samp, "]")
    } else {
      all_df <- all_df %>% mutate(sample_type = NA_character_)
      meta_note <- paste0(basename(meta_path), " [columns not found]")
    }
  }
}
if (!"sample_type" %in% names(all_df)) all_df <- all_df %>% mutate(sample_type = NA_character_)
all_df <- all_df %>% mutate(sample_type = ifelse(is.na(sample_type) | !nzchar(sample_type), "Unknown", sample_type))

# -- Ontology enrichment (optional) -------------------------------------------

load_aro_ontology <- function(onto_dir) {
  aro_path <- file.path(onto_dir, "aro.tsv")
  if (!file.exists(aro_path)) return(NULL)
  on <- tryCatch(readr::read_tsv(aro_path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(on) || nrow(on) == 0) return(NULL)
  # Expect columns: Accession, Name, Description, ID, CARD Short Name
  cn <- make_clean(names(on))
  names(on) <- cn
  # Standardize accession
  on <- on %>% mutate(
    accession_std = toupper(trimws(.data$accession)),
    id_numeric = suppressWarnings(as.numeric(.data$id)),
    name = .data$name,
    card_short_name = dplyr::coalesce(.data$card_short_name, .data$`card_short_name`),
    description = .data$description
  )
  on %>% dplyr::select(.data$accession_std, .data$id_numeric, .data$name, .data$card_short_name, .data$description)
}

aro_map <- NULL
if (!is.null(args$ontology_dir) && dir.exists(args$ontology_dir)) {
  aro_map <- load_aro_ontology(args$ontology_dir)
}

# Standardize ARO accession in results and enrich labels
if (!is.null(aro_map)) {
  std_acc <- function(x) {
    y <- toupper(trimws(as.character(x)))
    # Add prefix if missing and looks like ARO numeric
    y <- ifelse(grepl("^ARO:", y), y, ifelse(grepl("^[0-9]+$", y), paste0("ARO:", y), y))
    y
  }
  all_df <- all_df %>% mutate(
    aro_accession_std = std_acc(.data$aro_accession)
  )
  all_df <- all_df %>% left_join(aro_map, by = c("aro_accession_std" = "accession_std")) %>%
    mutate(
      gene_display = dplyr::coalesce(.data$card_short_name, .data$name, .data$gene)
    )
  # Optional ontology category drop (leaf genes only)
  drop_flag <- isTRUE(get_opt(args, "drop-ontology-categories"))
  if (isTRUE(drop_flag)) {
    drop_onto_pat <- "(family|group|antibiotic|efflux|resistance mechanism|target alteration|enzyme class|regulator)"
    all_df <- all_df %>%
      mutate(
        .keep_gene = (
          !is.na(.data$gene_display) & stringr::str_detect(.data$gene_display, "[A-Za-z].*[0-9\\)\\(\\'_-]")
        ) & !(
          (!is.na(.data$gene_display) & stringr::str_detect(.data$gene_display, stringr::regex(drop_onto_pat, ignore_case = TRUE))) |
          (!is.na(.data$gene_family) & stringr::str_detect(.data$gene_family, stringr::regex(drop_onto_pat, ignore_case = TRUE)))
        )
      ) %>%
      filter(.keep_gene) %>% select(-.keep_gene)
  }
} else {
  all_df <- all_df %>% mutate(gene_display = .data$gene)
}

# Persist raw standardized and enriched table
readr::write_csv(all_df, file.path(output_dir, "rgi_standardized_records.csv"))

# -- Summaries -----------------------------------------------------------------

has_gene <- !all(is.na(all_df$gene_display))
drop_singletons <- isTRUE(get_opt(args, "drop-singletons"))
has_drug <- !all(is.na(all_df$drug_class_primary))

class_counts <- if (has_drug) {
  dfc <- all_df %>% filter(!is.na(drug_class_primary))
  if (count_mode == "presence") {
    dfc %>% distinct(sample, drug_class_primary) %>% count(sample, drug_class_primary, name = "n")
  } else if (count_mode == "contig") {
    dfc %>% distinct(sample, contig, drug_class_primary) %>% count(sample, drug_class_primary, name = "n")
  } else if (count_mode == "orf") {
    dfc %>% distinct(sample, orf_id, drug_class_primary) %>% count(sample, drug_class_primary, name = "n")
  } else {
    dfc %>% count(sample, drug_class_primary, name = "n")
  }
} else tibble()

# Also export full Drug Class strings in summaries
if (has_drug) {
  dff <- all_df %>% filter(!is.na(.data$drug_class))
  class_counts_full <- if (count_mode == "presence") {
    dff %>% distinct(sample, drug_class = .data$drug_class) %>% count(sample, drug_class, name = "n")
  } else if (count_mode == "contig") {
    dff %>% distinct(sample, contig, drug_class = .data$drug_class) %>% count(sample, drug_class, name = "n")
  } else if (count_mode == "orf") {
    dff %>% distinct(sample, orf_id, drug_class = .data$drug_class) %>% count(sample, drug_class, name = "n")
  } else {
    dff %>% count(sample, drug_class = .data$drug_class, name = "n")
  }
  readr::write_csv(class_counts_full, file.path(output_dir, "counts_by_sample_and_drug_class_full.csv"))
}

sample_map <- all_df %>% distinct(sample, sample_type)
if (nrow(class_counts) > 0) class_counts <- class_counts %>% left_join(sample_map, by = "sample")
if (exists("class_counts_full")) class_counts_full <- class_counts_full %>% left_join(sample_map, by = "sample")

readr::write_csv(class_counts, file.path(output_dir, "counts_by_sample_and_drug_class.csv"))

# Type-level class counts
class_counts_by_type <- if (has_drug) {
  dfc <- all_df %>% filter(!is.na(drug_class_primary))
  if (count_mode == "presence") {
    dfc %>% distinct(sample_type, sample, drug_class_primary) %>% count(sample_type, drug_class_primary, name = "n")
  } else if (count_mode == "contig") {
    dfc %>% distinct(sample_type, sample, contig, drug_class_primary) %>% count(sample_type, drug_class_primary, name = "n")
  } else if (count_mode == "orf") {
    dfc %>% distinct(sample_type, sample, orf_id, drug_class_primary) %>% count(sample_type, drug_class_primary, name = "n")
  } else {
    dfc %>% count(sample_type, drug_class_primary, name = "n")
  }
} else tibble()
readr::write_csv(class_counts_by_type, file.path(output_dir, "counts_by_type_and_drug_class.csv"))

# Gene counts (use gene column)
if (!has_gene) {
  warning("Gene identifiers (ARO/Best_Hit_ARO) not found. Gene-based plots will be skipped.")
}

gene_counts <- if (has_gene) {
  dfg <- all_df %>% filter(!is.na(gene_display))
  if (drop_singletons) {
    # Pre-aggregate per sample+gene and drop n==1
    pre <- dfg %>% count(sample, gene = .data$gene_display, name = "n")
    pre <- pre %>% filter(n > 1)
    dfg <- dfg %>% semi_join(pre %>% select(sample, gene), by = c("sample", "gene_display" = "gene"))
  }
  if (count_mode == "presence") {
    dfg %>% distinct(sample, gene = .data$gene_display) %>% count(sample, gene, name = "n")
  } else if (count_mode == "contig") {
    dfg %>% distinct(sample, contig, gene = .data$gene_display) %>% count(sample, gene, name = "n")
  } else if (count_mode == "orf") {
    dfg %>% distinct(sample, orf_id, gene = .data$gene_display) %>% count(sample, gene, name = "n")
  } else {
    dfg %>% count(sample, gene = .data$gene_display, name = "n")
  }
} else tibble()

if (nrow(gene_counts) > 0) gene_counts <- gene_counts %>% left_join(sample_map, by = "sample")

readr::write_csv(gene_counts, file.path(output_dir, "counts_by_sample_and_gene.csv"))

gene_counts_by_type <- if (has_gene) {
  dfg <- all_df %>% filter(!is.na(gene_display))
  if (drop_singletons) {
    pre <- dfg %>% count(sample, gene = .data$gene_display, name = "n")
    pre <- pre %>% filter(n > 1)
    dfg <- dfg %>% semi_join(pre %>% select(sample, gene), by = c("sample", "gene_display" = "gene"))
  }
  if (count_mode == "presence") {
    dfg %>% distinct(sample_type, sample, gene = .data$gene_display) %>% count(sample_type, gene, name = "n")
  } else if (count_mode == "contig") {
    dfg %>% distinct(sample_type, sample, contig, gene = .data$gene_display) %>% count(sample_type, gene, name = "n")
  } else if (count_mode == "orf") {
    dfg %>% distinct(sample_type, sample, orf_id, gene = .data$gene_display) %>% count(sample_type, gene, name = "n")
  } else {
    dfg %>% count(sample_type, gene = .data$gene_display, name = "n")
  }
} else tibble()
readr::write_csv(gene_counts_by_type, file.path(output_dir, "counts_by_type_and_gene.csv"))

# Top-N genes overall
overall_gene_counts <- if (has_gene) gene_counts %>% group_by(gene) %>% summarise(n = sum(n), .groups = "drop") %>% arrange(desc(n)) else tibble()

top_n <- min(args$top_n, if (nrow(overall_gene_counts) > 0) nrow(overall_gene_counts) else 0)

top_genes <- if (top_n > 0) overall_gene_counts %>% slice_head(n = top_n) %>% pull(gene) else character(0)

# --- QA log and session info ---
readr::write_lines(c(
  paste("N records after parse:", nrow(all_df)),
  paste(
    "Filters -> model:", model_filter,
    "| id_min:", ifelse(is.na(identity_min), "none", identity_min),
    "cov_min:", ifelse(is.na(coverage_min), "none", coverage_min),
    "bits_min:", ifelse(is.na(bitscore_min), "none", bitscore_min),
    "| cutoff:", ifelse(!is.null(cutoff_levels) && nzchar(cutoff_levels), cutoff_levels, "none"),
    "| leaf-only:", isTRUE(get_opt(args, "drop-ontology-categories"))
  )
  , paste("Metadata:", meta_note)
  , paste("Config:", config_note)
  , paste("drop_singletons:", isTRUE(get_opt(args, "drop-singletons")))
  , paste("Heatmap:", paste0("normalize=", tolower(get_opt(args, "heatmap-normalize")), ", cluster=", isTRUE(get_opt(args, "heatmap-cluster")), ", dist=", tolower(get_opt(args, "cluster-distance")), ", link=", get_opt(args, "cluster-linkage")))
), file.path(output_dir, "rgi_viz_run.log"))

sink(file.path(output_dir, "sessionInfo.txt")); print(sessionInfo()); sink()

# -- Plots ---------------------------------------------------------------------

theme_set(theme_minimal(base_size = 12))
pal_name <- tolower(get_opt(args, "palette"))
wrap_width <- as.integer(get_opt(args, "label-wrap"))
lump_top_n <- as.integer(get_opt(args, "lump-top-n"))
discrete_pal <- function(n) get_palette(pal_name, n)
cont_pal <- function(n) get_palette(pal_name, n)

# Provided palette for stacked bars (prioritize these colors)
provided_stack_cols <- c(
  "#D98586", "#CC5877", "#DC8B9F", "#DCB1AF", # reds/pinks
  "#587F8E", "#4C96B7", "#5BBDCE", "#75CECD", "#75AECB", "#236D89" # blues/teals
)
get_provided_stack_palette <- function(levels_vec) {
  lv <- as.character(levels_vec)
  k <- length(lv)
  base <- provided_stack_cols
  if (k <= length(base)) return(setNames(base[seq_len(k)], lv))
  # If more categories than provided colors, fall back by reducing lump-top-n until it fits
  setNames(rep(base, length.out = k), lv)
}

# 1) Stacked bar by Drug Class per sample
if (nrow(class_counts) > 0) {
  # Order samples by total ARGs
  sample_order <- class_counts %>% group_by(sample) %>% summarise(tot = sum(n), .groups = "drop") %>% arrange(desc(tot)) %>% pull(sample)
  p1_data <- class_counts %>% mutate(drug_class_lumped = forcats::fct_lump_n(
      ifelse(is.na(drug_class_primary) | drug_class_primary == "", "Unknown", drug_class_primary),
      n = lump_top_n, other_level = "Other"
    ))
  n_levels <- nlevels(p1_data$drug_class_lumped)
  # Recompute sample order per type to sort stacks within facets
  ord_by_type <- p1_data %>% group_by(sample_type, sample) %>% summarise(tot = sum(n), .groups = "drop") %>%
    group_by(sample_type) %>% arrange(desc(tot), .by_group = TRUE) %>% mutate(sample_fac = paste(sample_type, sample, sep = "|"))
  p1 <- p1_data %>%
    mutate(sample_fac = paste(sample_type, sample, sep = "|")) %>%
    left_join(ord_by_type %>% select(sample_type, sample_fac), by = c("sample_type", "sample_fac")) %>%
    mutate(sample_fac = factor(sample_fac, levels = unique(ord_by_type$sample_fac)),
           sample_lab = forcats::fct_relabel(factor(sample), ~ wrap_lab(.x, width = wrap_width))) %>%
    ggplot(aes(x = sample_fac, y = n, fill = drug_class_lumped)) +
    geom_col(width = 0.85) +
  scale_fill_manual(values = discrete_pal(n_levels), name = "Drug Class") +
  scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_si(""))) +
    labs(title = "ARG drug classes per sample (RGI on contigs)", x = "Sample", y = "Count") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = if (n_levels > 20) "none" else "right"
    ) +
    facet_wrap(~ sample_type, scales = "free_x") +
    scale_x_discrete(labels = function(v) sub("^[^|]+\\|", "", v))
  save_plot(p1, file.path(output_dir, "bar_drug_class_by_sample"), width = max(8, length(sample_order) * 0.4), height = 6)
  # Aggregated by type
  p1_type_data <- class_counts_by_type %>% mutate(drug_class_lumped = forcats::fct_lump_n(
      ifelse(is.na(drug_class_primary) | drug_class_primary == "", "Unknown", drug_class_primary),
      n = lump_top_n, other_level = "Other"
    ))
  n_levels_type <- nlevels(p1_type_data$drug_class_lumped)
  type_order <- class_counts_by_type %>% group_by(sample_type) %>% summarise(tot = sum(n), .groups = "drop") %>% arrange(desc(tot)) %>% pull(sample_type)
  p1_type <- p1_type_data %>%
    mutate(sample_type = factor(sample_type, levels = type_order)) %>%
    ggplot(aes(x = sample_type, y = n, fill = drug_class_lumped)) +
    geom_col(position = position_stack(), width = 0.8) +
  scale_fill_manual(values = discrete_pal(n_levels_type), name = "Drug Class") +
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_si(""))) +
    labs(title = "ARG drug classes by sample type", x = "Sample Type", y = "Count") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = if (n_levels_type > 20) "none" else "right")
  save_plot(p1_type, file.path(output_dir, "bar_drug_class_by_type"), width = max(8, nlevels(factor(p1_type_data$sample_type)) * 1.8), height = 6)
}

# 2) Top-N genes across samples
if (has_gene && length(top_genes) > 0) {
  gene_counts_top <- gene_counts %>% filter(gene %in% top_genes)
  # Order genes by abundance within each sample_type facet (descending height)
  ord_by_type_gene <- gene_counts_top %>%
    group_by(sample_type, gene) %>% summarise(tot = sum(n), .groups = "drop") %>%
    group_by(sample_type) %>% arrange(desc(tot), .by_group = TRUE) %>%
    mutate(gene_fac = paste(sample_type, gene, sep = "|"))
  n_samples <- nlevels(factor(gene_counts_top$sample))
  p2 <- gene_counts_top %>%
    mutate(gene_fac = paste(sample_type, gene, sep = "|")) %>%
    left_join(ord_by_type_gene %>% select(sample_type, gene_fac), by = c("sample_type", "gene_fac")) %>%
    mutate(gene_fac = factor(gene_fac, levels = unique(ord_by_type_gene$gene_fac))) %>%
    ggplot(aes(x = n, y = gene_fac, fill = sample)) +
    geom_col(position = position_stack(), width = 0.8) +
    scale_fill_manual(values = discrete_pal(n_samples), name = "Sample") +
    scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_si(""))) +
    labs(title = paste0("Top-", top_n, " ARG genes across samples"), x = "Count", y = NULL) +
    theme(legend.position = if (n_samples > 25) "none" else "right") +
    facet_wrap(~ sample_type, scales = "free_y") +
    scale_y_discrete(labels = function(v) stringr::str_replace(v, "^[^|]+\\|", ""))
  save_plot(p2, file.path(output_dir, "bar_top_genes_stacked_by_sample"), width = 10, height = max(6, length(top_genes) * 0.18))

  # Aggregated by type (top-N by overall), sort genes by total height
  gene_counts_by_type_top <- gene_counts_by_type %>% filter(gene %in% top_genes)
  gene_order_type <- gene_counts_by_type_top %>% group_by(gene) %>% summarise(tot = sum(n), .groups = "drop") %>% arrange(desc(tot)) %>% pull(gene)
  p2_type <- gene_counts_by_type_top %>%
    mutate(sample_type = factor(sample_type, levels = type_order), gene = factor(gene, levels = gene_order_type)) %>%
    ggplot(aes(x = n, y = gene, fill = sample_type)) +
    geom_col(position = position_stack(), width = 0.8) +
    scale_fill_manual(values = discrete_pal(nlevels(factor(gene_counts_by_type_top$sample_type))), name = "Sample Type") +
    scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_si(""))) +
    labs(title = paste0("Top-", top_n, " ARG genes by sample type"), x = "Count", y = NULL) +
    theme(legend.position = "right")
  save_plot(p2_type, file.path(output_dir, "bar_top_genes_by_type"), width = 10, height = max(6, length(top_genes) * 0.18))
}

# 3) Heatmap of gene counts across samples (top-N)
if (has_gene && length(top_genes) > 0) {
  mat <- gene_counts %>%
    filter(gene %in% top_genes) %>%
    select(sample, gene, n) %>%
    tidyr::pivot_wider(names_from = sample, values_from = n, values_fill = 0L) %>%
    as.data.frame()
  rownames(mat) <- mat$gene
  mat$gene <- NULL
  # Normalization and clustering
  hm_norm <- tolower(get_opt(args, "heatmap-normalize"))
  do_cluster <- isTRUE(get_opt(args, "heatmap-cluster"))
  dist_m <- tolower(get_opt(args, "cluster-distance"))
  link_m <- get_opt(args, "cluster-linkage")
  m_proc <- normalize_matrix(mat, hm_norm)
  # Orders
  if (do_cluster) {
    row_ord <- cluster_order(m_proc, distance = dist_m, linkage = link_m, axis = "rows")
    col_ord <- cluster_order(m_proc, distance = dist_m, linkage = link_m, axis = "cols")
    mat_ord <- m_proc[row_ord, col_ord, drop = FALSE]
  } else {
    # fallback: order by totals
    row_ord <- order(rowSums(m_proc), decreasing = FALSE)
    col_ord <- order(colSums(m_proc), decreasing = TRUE)
    mat_ord <- m_proc[row_ord, col_ord, drop = FALSE]
  }
  # Convert to long for ggplot heatmap
  hm_long <- tibble::as_tibble(mat_ord, rownames = "gene") %>%
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "count") %>%
    left_join(sample_map, by = "sample")
  # Order axes for readability
  gene_order_hm <- unique(rownames(mat_ord))
  sample_order_hm <- unique(colnames(mat_ord))

  # Choose limits so maxima use full hue
  if (hm_norm == "row_z") {
    lims_hm <- c(min(hm_long$count, na.rm = TRUE), max(hm_long$count, na.rm = TRUE))
  } else {
    lims_hm <- c(0, max(hm_long$count, na.rm = TRUE))
  }

  p3 <- hm_long %>%
    mutate(gene = factor(gene, levels = gene_order_hm), sample = factor(sample, levels = sample_order_hm)) %>%
    ggplot(aes(x = sample, y = gene, fill = count)) +
    geom_tile(color = "grey90", linewidth = 0.1) +
  heatmap_scale_layer(hm_norm, pal_name, lims = lims_hm) +
    labs(title = paste0("ARG gene counts heatmap (Top-", top_n, ")"), x = "Sample", y = NULL) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right"
  ) +
  facet_wrap(~ sample_type, scales = "free_x")
  save_plot(p3, file.path(output_dir, "heatmap_top_genes_by_sample"), width = max(8, length(sample_order_hm) * 0.35), height = max(6, length(gene_order_hm) * 0.18))
  # Also export a single clustered heatmap without facets (for clear clustering view)
  p3c <- tibble::as_tibble(mat_ord, rownames = "gene") %>%
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "count") %>%
    ggplot(aes(x = factor(sample, levels = sample_order_hm), y = factor(gene, levels = gene_order_hm), fill = count)) +
    geom_tile(color = "grey90", linewidth = 0.1) +
  heatmap_scale_layer(hm_norm, pal_name, lims = lims_hm) +
    labs(title = paste0("Clustered ARG gene counts heatmap (Top-", top_n, ")"), x = "Sample", y = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p3c, file.path(output_dir, "heatmap_top_genes_clustered"), width = max(8, length(sample_order_hm) * 0.35), height = max(6, length(gene_order_hm) * 0.18))
  # Heatmap aggregated by type
  hm_type <- gene_counts_by_type %>% filter(gene %in% top_genes) %>%
    tidyr::pivot_wider(names_from = sample_type, values_from = n, values_fill = 0L)
  mat_t <- as.data.frame(hm_type)
  rownames(mat_t) <- mat_t$gene
  mat_t$gene <- NULL
  hm_long_t <- tibble::as_tibble(mat_t, rownames = "gene") %>% tidyr::pivot_longer(-gene, names_to = "sample_type", values_to = "count")
  gene_order_t <- hm_long_t %>% group_by(gene) %>% summarise(tot = sum(count), .groups = "drop") %>% arrange(tot) %>% pull(gene)
  type_order_t <- hm_long_t %>% group_by(sample_type) %>% summarise(tot = sum(count), .groups = "drop") %>% arrange(desc(tot)) %>% pull(sample_type)
  # Limits for type heatmap
  if (hm_norm == "row_z") {
    lims_hm_t <- c(min(hm_long_t$count, na.rm = TRUE), max(hm_long_t$count, na.rm = TRUE))
  } else {
    lims_hm_t <- c(0, max(hm_long_t$count, na.rm = TRUE))
  }
  p3_type <- hm_long_t %>% mutate(gene = factor(gene, levels = gene_order_t), sample_type = factor(sample_type, levels = type_order_t)) %>%
    ggplot(aes(x = sample_type, y = gene, fill = count)) +
    geom_tile(color = "grey90", linewidth = 0.1) +
    heatmap_scale_layer(hm_norm, pal_name, lims = lims_hm_t) +
    labs(title = paste0("ARG gene counts heatmap by type (Top-", top_n, ")"), x = "Sample Type", y = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p3_type, file.path(output_dir, "heatmap_top_genes_by_type"), width = max(8, length(type_order_t) * 0.7), height = max(6, length(gene_order_t) * 0.18))
}

# Bonus: mechanism breakdown if available
if (!all(is.na(all_df$mechanism_primary)) && !all(is.na(all_df$drug_class_primary))) {
  mech_counts <- all_df %>%
    mutate(mechanism_primary = ifelse(is.na(mechanism_primary) | mechanism_primary == "", "Unknown", mechanism_primary)) %>%
    count(drug_class_primary, mechanism_primary, name = "n")
  readr::write_csv(mech_counts, file.path(output_dir, "counts_by_drug_class_and_mechanism.csv"))

  # Rank mechanisms by total across drug classes, pick top-N, and sort y by total height
  top_n_mech <- 20L
  mech_totals <- mech_counts %>% group_by(mechanism_primary) %>% summarise(tot = sum(n), .groups = "drop") %>% arrange(desc(tot))
  top_mechs <- head(mech_totals$mechanism_primary, top_n_mech)
  order_levels <- rev(top_mechs)  # highest at top of the plot

  p4_data <- mech_counts %>% filter(mechanism_primary %in% top_mechs) %>%
    mutate(
      mechanism_primary = factor(mechanism_primary, levels = order_levels),
      drug_class_lumped = forcats::fct_lump_n(
        ifelse(is.na(drug_class_primary) | drug_class_primary == "", "Unknown", drug_class_primary),
        n = lump_top_n, other_level = "Other"
      )
    )
  n_levels_mech <- nlevels(p4_data$drug_class_lumped)
  p4 <- p4_data %>%
    ggplot(aes(x = n, y = mechanism_primary, fill = drug_class_lumped)) +
    geom_col(position = position_stack(), width = 0.8) +
  scale_fill_manual(values = discrete_pal(n_levels_mech), name = "Drug Class") +
    scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_si(""))) +
    labs(title = "Top mechanisms by drug class (overall)", x = "Count", y = NULL) +
    theme(legend.position = if (n_levels_mech > 20) "none" else "right")
  save_plot(p4, file.path(output_dir, "bar_top_mechanisms_overall"), width = 10, height = 8)
}

# 4) Circular packing: mechanism composition per sample -----------------------
if (!all(is.na(all_df$mechanism_primary))) {
  # Ensure required packages
  have_ggraph <- require_pkg("ggraph", install_if_missing = isTRUE(args$install_missing))
  have_ig <- require_pkg("igraph", install_if_missing = isTRUE(args$install_missing))
  have_ggforce <- require_pkg("ggforce", install_if_missing = isTRUE(args$install_missing))
  if (have_ggraph && have_ig) {
    # Build per-sample mechanism counts and lump long tail for readability
    cc_lumped <- all_df %>%
      filter(!is.na(mechanism_primary)) %>%
      mutate(mechanism_lumped = forcats::fct_lump_n(
        ifelse(is.na(mechanism_primary) | mechanism_primary == "", "Unknown", mechanism_primary),
        n = lump_top_n, other_level = "Other"
      )) %>%
      count(sample, sample_type, mechanism_lumped, name = "n")

    # Nodes (root, samples, classes) and edges
    root_id <- "ROOT"
    samples <- cc_lumped %>% distinct(sample, sample_type)
    classes <- cc_lumped %>% mutate(mechanism_label = as.character(mechanism_lumped)) %>%
      mutate(node_id = paste(sample, mechanism_label, sep = "::")) %>% select(sample, sample_type, mechanism_label, node_id, n)
    
    nodes_root <- tibble::tibble(id = root_id, name = "All", type = "root", sample = NA_character_, sample_type = NA_character_, mechanism_label = NA_character_, weight = 0)
    nodes_samples <- samples %>% mutate(id = paste0("S:", sample), name = sample, type = "sample", mechanism_label = NA_character_, weight = 0) %>%
      select(id, name, type, sample, sample_type, mechanism_label, weight)
    nodes_classes <- classes %>% mutate(id = paste0("M:", node_id), name = mechanism_label, type = "class", weight = as.numeric(n)) %>%
      select(id, name, type, sample, sample_type, mechanism_label, weight)
    nodes <- dplyr::bind_rows(nodes_root, nodes_samples, nodes_classes)

  edges_root <- samples %>% transmute(from = root_id, to = paste0("S:", sample))
  edges_sc <- classes %>% transmute(from = paste0("S:", sample), to = paste0("M:", node_id))
    edges <- dplyr::bind_rows(edges_root, edges_sc)

    g <- igraph::graph_from_data_frame(edges, vertices = nodes, directed = TRUE)
    w <- igraph::V(g)$weight
    w[is.na(w)] <- 0

    lay <- ggraph::create_layout(g, layout = "circlepack", weight = w)

    # Palette for classes
  cls_levels <- sort(unique(na.omit(nodes_classes$mechanism_label)))
  pal_vals <- discrete_pal(length(cls_levels))
  names(pal_vals) <- cls_levels
    names(pal_vals) <- cls_levels

    p_cp <- ggplot2::ggplot() +
      ggraph::geom_node_circle(data = subset(lay, type == "sample"), ggplot2::aes(x0 = x, y0 = y, r = r),
                               linewidth = 0.2, color = "grey40", fill = NA) +
  ggraph::geom_node_circle(data = subset(lay, type == "class"), ggplot2::aes(x0 = x, y0 = y, r = r, fill = mechanism_label),
                               linewidth = 0.1, color = "grey85") +
  ggplot2::scale_fill_manual(values = pal_vals, name = "Mechanism") +
      ggplot2::coord_fixed() +
  ggplot2::labs(title = "Circular packing: mechanism composition per sample",
        subtitle = "Each large circle is a sample; inner circles are mechanisms sized by counts",
                    x = NULL, y = NULL) +
      ggplot2::theme_void(base_size = 12) +
      ggplot2::theme(legend.position = "right")

    # Label top-N largest sample circles (use ggrepel if available)
    samp_nodes <- subset(lay, type == "sample")
    if (nrow(samp_nodes) > 0) {
      tot_by_sample <- cc_lumped %>% group_by(sample) %>% summarise(tot = sum(n), .groups = "drop")
      samp_nodes <- dplyr::left_join(samp_nodes, tot_by_sample, by = c("name" = "sample")) %>%
        dplyr::arrange(dplyr::desc(.data$tot)) %>% dplyr::mutate(.keep_lab = seq_len(dplyr::n()) <= 12)
      lab_dat <- subset(samp_nodes, .keep_lab)
      if (require_pkg("ggrepel", install_if_missing = isTRUE(args$install_missing))) {
        p_cp <- p_cp + ggrepel::geom_text_repel(data = lab_dat,
                          ggplot2::aes(x = x, y = y, label = name),
                          size = 3, color = "grey15", seed = 42, max.overlaps = 100,
                          min.segment.length = 0, segment.alpha = 0.5)
      } else {
        p_cp <- p_cp + ggraph::geom_node_text(data = lab_dat,
                          ggplot2::aes(x = x, y = y, label = name), size = 3, color = "grey15")
      }
    }

  save_plot(p_cp, file.path(output_dir, "circlepack_mechanism_per_sample"), width = 12, height = 10)
  } else {
    warning("Skipping circular packing plot: required packages ggraph, igraph not available.")
  }
}

# 5) Radial hub: central mechanisms with outer sample circle-packs connected ----
if (!all(is.na(all_df$mechanism_primary)) && nrow(class_counts) > 0) {
  have_ggraph <- require_pkg("ggraph", install_if_missing = isTRUE(args$install_missing))
  have_ig <- require_pkg("igraph", install_if_missing = isTRUE(args$install_missing))
  if (have_ggraph && have_ig) {
    # Reuse class circlepack at sample->class level
    cc_lumped <- class_counts %>%
      mutate(drug_class_lumped = forcats::fct_lump_n(
        ifelse(is.na(drug_class_primary) | drug_class_primary == "", "Unknown", drug_class_primary),
        n = lump_top_n, other_level = "Other"
      )) %>% group_by(sample, sample_type, drug_class_lumped) %>% summarise(n = sum(n), .groups = "drop")

    # Map each sample+drug_class_lumped to a dominant mechanism for coloring
    class_mech_map <- all_df %>%
      filter(!is.na(drug_class_primary), !is.na(mechanism_primary)) %>%
      mutate(drug_class_lumped = forcats::fct_lump_n(
        ifelse(is.na(drug_class_primary) | drug_class_primary == "", "Unknown", drug_class_primary),
        n = lump_top_n, other_level = "Other"
      )) %>%
      count(sample, sample_type, drug_class_lumped, mechanism_primary, name = "w") %>%
      group_by(sample, drug_class_lumped) %>%
      slice_max(order_by = w, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(sample, sample_type, drug_class_lumped, mech_for_class = mechanism_primary)

    # Build graph for circlepack
    samples <- cc_lumped %>% distinct(sample, sample_type)
    classes <- cc_lumped %>% mutate(drug_class_label = as.character(drug_class_lumped)) %>%
      left_join(class_mech_map, by = c("sample", "sample_type", "drug_class_lumped")) %>%
      mutate(mech_for_class = ifelse(is.na(mech_for_class) | mech_for_class == "", "Unknown", mech_for_class)) %>%
      mutate(node_id = paste(sample, drug_class_label, sep = "::")) %>% select(sample, sample_type, drug_class_label, mech_for_class, node_id, n)
    nodes_root <- tibble::tibble(id = "ROOT", name = "All", type = "root", sample = NA_character_, sample_type = NA_character_, drug_class_label = NA_character_, weight = 0)
    nodes_samples <- samples %>% mutate(id = paste0("S:", sample), name = sample, type = "sample", drug_class_label = NA_character_, weight = 0) %>%
      select(id, name, type, sample, sample_type, drug_class_label, weight)
    nodes_classes <- classes %>% mutate(id = paste0("C:", node_id), name = drug_class_label, type = "class", weight = as.numeric(n)) %>%
      select(id, name, type, sample, sample_type, drug_class_label, mech_for_class, weight)
    nodes <- dplyr::bind_rows(nodes_root, nodes_samples, nodes_classes)
    edges_root <- samples %>% transmute(from = "ROOT", to = paste0("S:", sample))
    edges_sc <- classes %>% transmute(from = paste0("S:", sample), to = paste0("C:", node_id))
    g_pack <- igraph::graph_from_data_frame(dplyr::bind_rows(edges_root, edges_sc), vertices = nodes, directed = TRUE)
    w_pack <- igraph::V(g_pack)$weight; w_pack[is.na(w_pack)] <- 0
    lay <- ggraph::create_layout(g_pack, layout = "circlepack", weight = w_pack)

    # Compute new positions: push each sample pack to outer circle
    samp <- subset(lay, type == "sample")
    cls <- subset(lay, type == "class")
  N <- nrow(samp); if (N < 1) N <- 1
  R <- max(samp$r, na.rm = TRUE) * 8
  # Deterministic alphabetical arrangement of sample packs around the ring
  samp_ord <- samp[order(samp$name), , drop = FALSE]
  theta <- 2 * pi * (seq_len(N) - 1) / N
  samp_ord$x <- R * cos(theta)
  samp_ord$y <- R * sin(theta)
  samp_ord$theta <- atan2(samp_ord$y, samp_ord$x)
  samp_new <- samp_ord[match(samp$name, samp_ord$name), , drop = FALSE]
    # Map old->new centers by sample name
    map_centers <- samp %>% dplyr::select(name, x0 = x, y0 = y) %>% dplyr::left_join(samp_new %>% dplyr::select(name, x1 = x, y1 = y), by = "name")
    cls_new <- cls %>% dplyr::left_join(map_centers, by = c("sample" = "name")) %>%
      dplyr::mutate(x = x + (x1 - x0), y = y + (y1 - y0)) %>% dplyr::select(-x0, -y0, -x1, -y1)

    # Central mechanisms arranged on small inner circle
    mech_counts <- all_df %>% mutate(mechanism_primary = ifelse(is.na(mechanism_primary) | mechanism_primary == "", "Unknown", mechanism_primary))
    mech_tot <- mech_counts %>% count(mechanism_primary, name = "weight") %>% dplyr::arrange(dplyr::desc(weight))
    M <- nrow(mech_tot); rM <- R * 0.25
    # Arrange mechanisms to reduce crossings: align to sample sector order by dominant sample
    # 1) Determine sample order by angle
    samp_angle_map <- samp_ord %>% dplyr::transmute(sample = name, theta = atan2(y, x)) %>% dplyr::arrange(theta) %>% dplyr::mutate(pos = dplyr::row_number())
    # 2) For each mechanism, find dominant sample (max weight)
    mech_dom <- mech_counts %>% dplyr::count(mechanism_primary, sample, name = "w") %>%
      dplyr::group_by(mechanism_primary) %>% dplyr::slice_max(order_by = w, n = 1, with_ties = FALSE) %>% dplyr::ungroup() %>%
      dplyr::left_join(samp_angle_map, by = "sample")
    # 3) Order mechanisms by the position of their dominant samples
    mech_order <- mech_dom %>% dplyr::arrange(pos) %>% dplyr::pull(mechanism_primary)
    mech_tot <- mech_tot %>% dplyr::mutate(mechanism_primary = factor(mechanism_primary, levels = unique(c(mech_order, as.character(mechanism_primary))))) %>% dplyr::arrange(mechanism_primary)
    thM <- 2 * pi * (seq_len(M) - 1) / max(M, 1)
    mech_df <- mech_tot %>% dplyr::mutate(x = rM * cos(thM), y = rM * sin(thM))

    # Edges from mechanisms to sample centers weighted by shared counts
    ms <- mech_counts %>% count(mechanism_primary, sample, name = "w") %>% dplyr::left_join(samp_new %>% dplyr::select(sample = name, xs = x, ys = y), by = "sample") %>%
      dplyr::left_join(mech_df %>% dplyr::select(mechanism_primary, xm = x, ym = y), by = "mechanism_primary")
    # Scale linewidth within each mechanism for local contrast
  ms <- ms %>% dplyr::group_by(mechanism_primary) %>% dplyr::mutate(lw = scales::rescale(w, to = c(0, 1))) %>% dplyr::ungroup()
  # Prepare simple sector-based bundling control points
  bund_bins <- min(24L, max(8L, ceiling(N / 2)))
  theta_s <- atan2(ms$ys, ms$xs)
  sector <- floor((theta_s + pi) / (2 * pi) * bund_bins)
  center_angle <- (sector + 0.5) * (2 * pi / bund_bins) - pi
  r_mid <- R * 0.55
  ms$cx <- r_mid * cos(center_angle)
  ms$cy <- r_mid * sin(center_angle)
  ms$edge_id <- paste(ms$mechanism_primary, ms$sample, sep = "||")
  edge_pts <- dplyr::bind_rows(
    dplyr::transmute(ms, edge_id, mechanism_primary, lw, x = xm, y = ym, ord = 1L),
    dplyr::transmute(ms, edge_id, mechanism_primary, lw, x = cx, y = cy, ord = 2L),
    dplyr::transmute(ms, edge_id, mechanism_primary, lw, x = xs, y = ys, ord = 3L)
  )

    # Palette for mechanisms
    mech_levels <- mech_df$mechanism_primary
    mech_cols <- discrete_pal(length(mech_levels)); names(mech_cols) <- mech_levels

    p_hub <- ggplot2::ggplot() +
      # edges colored by mechanism (bundled when ggforce is available)
      {
        if (have_ggforce) {
          ggforce::geom_bezier(data = edge_pts, ggplot2::aes(x = x, y = y, group = edge_id, color = mechanism_primary, linewidth = lw),
                               alpha = 0.65, show.legend = FALSE, lineend = "round")
        } else {
          ggplot2::geom_curve(data = ms, ggplot2::aes(x = xm, y = ym, xend = xs, yend = ys, color = mechanism_primary, linewidth = lw),
                               curvature = 0.15, alpha = 0.65, lineend = "round", show.legend = FALSE)
        }
      } +
      # central mechanisms colored by mechanism
      ggplot2::geom_point(data = mech_df, ggplot2::aes(x = x, y = y, size = weight, fill = mechanism_primary), shape = 21, color = "#1B3E5A", stroke = 0.3, show.legend = FALSE) +
      # sample pack outlines
      ggraph::geom_node_circle(data = samp_new, ggplot2::aes(x0 = x, y0 = y, r = r), linewidth = 0.25, color = "grey40", fill = NA) +
  # class circles (packed per sample) colored by dominant mechanism
  ggraph::geom_node_circle(data = cls_new, ggplot2::aes(x0 = x, y0 = y, r = r, fill = mech_for_class), linewidth = 0.1, color = "grey85", alpha = 0.65, show.legend = TRUE) +
  # (mechanism labels removed in favor of legend)
      # sample labels for each outer pack
      {
        if (require_pkg("ggrepel", install_if_missing = isTRUE(args$install_missing))) {
          ggrepel::geom_text_repel(data = samp_new, ggplot2::aes(x = x, y = y, label = name), seed = 42, size = 2.8,
                                   max.overlaps = Inf, box.padding = 0.3, point.padding = 0.2,
                                   min.segment.length = 0, segment.alpha = 0.35, segment.size = 0.2, show.legend = FALSE)
        } else {
          ggplot2::geom_text(data = samp_new, ggplot2::aes(x = x, y = y, label = name), size = 2.8, color = "#333333", show.legend = FALSE)
        }
      } +
  ggplot2::scale_color_manual(values = mech_cols, guide = "none") +
  ggplot2::scale_linewidth(range = c(0.02, 2), guide = "none") +
      ggplot2::scale_fill_manual(values = mech_cols, name = "Mechanism") +
      ggplot2::scale_size_continuous(range = c(3, 12), guide = "none") +
      ggplot2::coord_fixed() +
      ggplot2::theme_void(base_size = 12) +
  ggplot2::labs(title = "Mechanism hub with outer sample circle-packs",
        subtitle = "Color encodes mechanism (nodes, edges, and class circles in outer packs via dominant mechanism). Edge width ~ shared mechanism counts (thin range).")

    save_plot(p_hub, file.path(output_dir, "radial_mechanism_hub"), width = 12, height = 12)

    # Faceted hub: one panel per mechanism highlighting that mechanism, greying others
    if (nrow(mech_df) > 0) {
      focus_levels <- as.character(mech_df$mechanism_primary)
      focus_levels <- unique(focus_levels[!is.na(focus_levels)])

      # Build faceted datasets
      # Edges (bundled points), replicate across focus mechanisms
      edges_fac <- edge_pts %>%
        dplyr::mutate(mechanism_primary = factor(mechanism_primary, levels = mech_levels)) %>%
        tidyr::crossing(focus_mech = mech_levels) %>%
        dplyr::mutate(
          is_focus = .data$mechanism_primary == .data$focus_mech,
          color_key = dplyr::if_else(.data$is_focus, as.character(.data$mechanism_primary), "Other"),
          alpha_v = dplyr::if_else(.data$is_focus, 0.65, 0.08),
          lw_v = dplyr::if_else(.data$is_focus, .data$lw, pmax(.data$lw * 0.2, 0.02)),
          group_id = paste(.data$focus_mech, .data$edge_id, sep = "||")
        ) %>%
        dplyr::arrange(.data$group_id, .data$ord)

      # Mechanism center nodes
      mech_fac <- mech_df %>%
        tidyr::crossing(focus_mech = mech_levels) %>%
        dplyr::mutate(
          is_focus = .data$mechanism_primary == .data$focus_mech,
          fill_key = dplyr::if_else(.data$is_focus, as.character(.data$mechanism_primary), "Other"),
          alpha_v = dplyr::if_else(.data$is_focus, 0.95, 0.15)
        )

      # Sample ring and class circles (replicate across facets; render class circles in neutral greys to reduce distraction)
      samp_fac <- samp_new %>% tidyr::crossing(focus_mech = mech_levels)
      cls_fac <- cls_new %>% tidyr::crossing(focus_mech = mech_levels) %>%
        dplyr::mutate(
          fill_cls = dplyr::if_else(.data$mech_for_class == .data$focus_mech, as.character(.data$focus_mech), "Other"),
          alpha_cls = dplyr::if_else(.data$mech_for_class == .data$focus_mech, 0.65, 0.25)
        )

      # Palette with 'Other' grey for non-focused elements
  mech_cols_full <- c(mech_cols, Other = "#878771")

      nF <- length(mech_levels)
      ncol_grid <- ceiling(sqrt(nF))

      p_hub_fac <- ggplot2::ggplot() +
        {
          if (have_ggforce) {
            ggforce::geom_bezier(
              data = edges_fac,
              ggplot2::aes(x = x, y = y, group = group_id, color = color_key, linewidth = lw_v, alpha = alpha_v),
              lineend = "round", show.legend = FALSE
            )
          } else {
            ggplot2::geom_curve(
              data = ms %>% tidyr::crossing(focus_mech = mech_levels) %>% dplyr::mutate(
                is_focus = .data$mechanism_primary == .data$focus_mech,
                color_key = dplyr::if_else(.data$is_focus, as.character(.data$mechanism_primary), "Other"),
                alpha_v = dplyr::if_else(.data$is_focus, 0.65, 0.08),
                lw_v = dplyr::if_else(.data$is_focus, .data$lw, pmax(.data$lw * 0.2, 0.02))
              ),
              ggplot2::aes(x = xm, y = ym, xend = xs, yend = ys, color = color_key, linewidth = lw_v, alpha = alpha_v),
              curvature = 0.15, lineend = "round", show.legend = FALSE
            )
          }
        } +
        ggplot2::geom_point(
          data = mech_fac,
          ggplot2::aes(x = x, y = y, fill = fill_key, alpha = alpha_v, size = weight),
          shape = 21, color = "#1B3E5A", stroke = 0.25, show.legend = FALSE
        ) +
        ggraph::geom_node_circle(data = samp_fac, ggplot2::aes(x0 = x, y0 = y, r = r), linewidth = 0.2, color = "grey65", fill = NA) +
  ggraph::geom_node_circle(data = cls_fac, ggplot2::aes(x0 = x, y0 = y, r = r, fill = fill_cls, alpha = alpha_cls), linewidth = 0.08, color = "grey85", show.legend = FALSE) +
  ggplot2::scale_color_manual(values = mech_cols_full) +
  ggplot2::scale_fill_manual(values = mech_cols_full) +
        ggplot2::scale_alpha_identity() +
        ggplot2::scale_size_continuous(range = c(3, 12)) +
        ggplot2::scale_linewidth(range = c(0.02, 2)) +
        ggplot2::facet_wrap(~focus_mech, ncol = ncol_grid) +
        ggplot2::coord_fixed() +
        ggplot2::theme_void(base_size = 11) +
        ggplot2::theme(strip.text = ggplot2::element_text(color = "#333333", face = "bold")) +
        ggplot2::labs(title = "Mechanism-focused facets of radial hub", subtitle = "In each facet, the selected mechanism is colored; others are greyed out.")

      save_plot(p_hub_fac, file.path(output_dir, "radial_mechanism_hub_facets"), width = 14, height = 12)
    }
  }
}

message("Done. Outputs written to: ", normalizePath(output_dir))
