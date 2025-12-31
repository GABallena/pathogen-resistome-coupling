#!/usr/bin/env Rscript
# -------------------------------------------------------------------------
# Portfolio-ready note:
# - No hard-coded local paths.
# - Run from the project root (or set PROJECT_ROOT env var).
# -------------------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = getwd())
try(setwd(PROJECT_ROOT), silent = TRUE)

suppressPackageStartupMessages({
  library(tidyverse); library(readr); library(stringr); library(janitor)
})

# Inputs
merged_path <- "diversity_results/metaphlan_bracken_merged_both.tsv"
pathogen_list <- "clean_pathogens.txt"

dir.create("results", showWarnings = FALSE)

cat("[1/3] Reading merged taxonomy and pathogen list...\n")
if(!file.exists(merged_path)) stop("Missing ", merged_path)
if(!file.exists(pathogen_list)) stop("Missing ", pathogen_list)

merged <- suppressMessages(read_tsv(merged_path, show_col_types = FALSE)) %>% janitor::clean_names()
# Expecting columns like: sample, taxon, rank, metaphlan, bracken
stopifnot(all(c("sample","taxon","rank") %in% names(merged)))

# Use best available abundance (prefer bracken, fallback metaphlan)
if(!"bracken" %in% names(merged)) merged$bracken <- NA_real_
if(!"metaphlan" %in% names(merged)) merged$metaphlan <- NA_real_
merged <- merged %>% mutate(abund = dplyr::coalesce(bracken, metaphlan, 0))

# Parse pathogen genera from list (first token up to space or underscore)
pat_lines <- readLines(pathogen_list, warn = FALSE)
pat_genera <- tibble(line = pat_lines) %>%
  filter(nchar(trimws(line)) > 0) %>%
  mutate(genus = str_extract(line, "^[A-Za-z]+")) %>%
  filter(!is.na(genus)) %>% distinct(genus) %>% arrange(genus)
cat("Unique pathogen genera:", nrow(pat_genera), "\n")

# Filter merged to genus rank and build wide matrix of genus % per sample (zeros allowed)
genus_tbl <- merged %>% filter(tolower(rank) == "genus") %>%
  transmute(
    sample,
    genus = dplyr::case_when(
      str_detect(taxon, "^g__") ~ str_replace(taxon, "^g__", ""),
      TRUE ~ taxon
    ),
    value = abund
  ) %>%
  group_by(sample, genus) %>% summarise(value = suppressWarnings(max(value, na.rm=TRUE)), .groups='drop') %>%
  mutate(value = ifelse(is.infinite(value) | is.na(value), 0, value))

# Presence/absence just for pathogen genera subset
P <- genus_tbl %>% filter(genus %in% pat_genera$genus) %>%
  mutate(present = value > 0) %>%
  group_by(sample, genus) %>% summarise(present = any(present), value = max(value, na.rm=TRUE), .groups='drop') %>%
  mutate(value = ifelse(is.infinite(value) | is.na(value), 0, value))

# Sample-level metrics: % abundance in pathogen genera and richness (# present genera)
metrics <- P %>% group_by(sample) %>% summarise(
  pathogen_pct = sum(value, na.rm=TRUE),
  pathogen_richness = sum(present, na.rm=TRUE), .groups='drop')

# Wide table for genera abundances (allow zeros)
all_samp <- sort(unique(genus_tbl$sample))
pat_wide <- genus_tbl %>% filter(genus %in% pat_genera$genus) %>%
  tidyr::pivot_wider(names_from = genus, values_from = value, values_fill = 0) %>%
  arrange(sample)

readr::write_tsv(metrics, "results/pathogen_sample_metrics.tsv")
readr::write_tsv(pat_wide, "results/pathogen_genus_wide.tsv")

cat("[3/3] Wrote:\n - results/pathogen_sample_metrics.tsv\n - results/pathogen_genus_wide.tsv\n")
