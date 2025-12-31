#!/usr/bin/env Rscript
# -------------------------------------------------------------------------
# Portfolio-ready note:
# - No hard-coded local paths.
# - Run from the project root (or set PROJECT_ROOT env var).
# -------------------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = getwd())
try(setwd(PROJECT_ROOT), silent = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(glue)
  library(janitor)
  library(patchwork)
})

# Config
path_genus_table <- "diversity_results/metaphlan_bracken_merged_both.tsv"
path_pathogen_list <- "clean_pathogens.txt"
path_class_wide <- "results/class_abundance_cpm_wide.tsv"
out_dir <- "results"
fig_dir <- "figs"
pa_threshold <- 0 # presence if abundance > threshold
min_non_na <- 5   # require at least this many paired samples per correlation

# Helpers
read_pathogen_list <- function(f){
  if(!file.exists(f)) stop("Pathogen list not found: ", f)
  # Expect one taxon per line; keep genus part only when available
  tx <- read_lines(f) %>% str_trim() %>% purrr::discard(~.x=="")
  tibble::tibble(raw = tx) %>%
    dplyr::mutate(genus = str_extract(raw, "^[A-Za-z]+")) %>%
    dplyr::filter(!is.na(genus)) %>% dplyr::distinct(genus)
}

message("[1/5] Loading inputs…")
stopifnot(file.exists(path_genus_table))
stopifnot(file.exists(path_class_wide))
patho <- read_pathogen_list(path_pathogen_list)

# Genus table: keep rank==genus and create presence/absence per genus per sample
message("[2/5] Building presence/absence matrix for pathogen genera…")
gtbl <- suppressMessages(read_tsv(path_genus_table, show_col_types = FALSE))
req_cols <- c("sample","taxon","rank","metaphlan","bracken")
missing <- setdiff(req_cols, names(gtbl))
if(length(missing)>0) stop("Missing columns in genus table: ", paste(missing, collapse=", "))

genus_only <- gtbl %>%
  filter(rank=="genus") %>%
  # standardize sample IDs and genus names
  mutate(sample = as.character(sample), genus = as.character(taxon)) %>%
  # Presence if either profiler reports > 0: use max of bracken/metaphlan abundances
  mutate(abund = pmax(dplyr::coalesce(bracken, 0), dplyr::coalesce(metaphlan, 0))) %>%
  select(sample, genus, abund)

# Filter to pathogen genera list (match by exact genus)
patho_genus <- unique(patho$genus)
message("Unique pathogen genera in list: ", length(patho_genus))

pa_long <- genus_only %>%
  mutate(genus_extracted = str_extract(genus, "^[A-Za-z]+")) %>%
  filter(genus_extracted %in% patho_genus) %>%
  transmute(sample, genus = genus_extracted, present = as.integer(abund > pa_threshold)) %>%
  group_by(sample, genus) %>% summarise(present = as.integer(any(present>0)), .groups="drop")

# Wide presence matrix
pa_wide <- pa_long %>%
  pivot_wider(names_from = genus, values_from = present, values_fill = 0) %>%
  arrange(sample)

message("[3/5] Loading ARG class CPM wide…")
class_wide <- suppressMessages(read_tsv(path_class_wide, show_col_types = FALSE)) %>%
  clean_names()
# Ensure 'sample' column exists
if(!"sample" %in% names(class_wide)){
  stop("Class CPM wide table missing 'sample' column: ", path_class_wide)
}

# Align samples: keep intersection
common_samples <- intersect(pa_wide$sample, class_wide$sample)
pa_wide <- pa_wide %>% filter(sample %in% common_samples) %>% arrange(sample)
class_wide <- class_wide %>% filter(sample %in% common_samples) %>% arrange(sample)
stopifnot(nrow(pa_wide) == nrow(class_wide))

# Identify class columns (everything except 'sample')
class_cols <- setdiff(names(class_wide), "sample")
patho_cols <- setdiff(names(pa_wide), "sample")

message("[4/5] Computing correlations (point-biserial via Pearson; Spearman fallback)…")
results <- map_dfr(patho_cols, function(pg){
  y <- pa_wide[[pg]]
  map_dfr(class_cols, function(ac){
    x <- class_wide[[ac]]
    # require minimal non-NA paired
    ok <- is.finite(x) & is.finite(y)
    if(sum(ok) < min_non_na) return(tibble(genus=pg, arg_class=ac, n=sum(ok), r=NA_real_, p=NA_real_, method=NA_character_))
    # point-biserial = Pearson between numeric x and binary y
    ct <- tryCatch({ suppressWarnings(cor.test(x[ok], y[ok], method="pearson")) }, error=function(e) NULL)
    if(!is.null(ct)){
      return(tibble(genus=pg, arg_class=ac, n=sum(ok), r=unname(ct$estimate), p=ct$p.value, method="pearson"))
    } else {
      ct2 <- tryCatch({ suppressWarnings(cor.test(x[ok], y[ok], method="spearman")) }, error=function(e) NULL)
      if(!is.null(ct2)) return(tibble(genus=pg, arg_class=ac, n=sum(ok), r=unname(ct2$estimate), p=ct2$p.value, method="spearman"))
    }
    tibble(genus=pg, arg_class=ac, n=sum(ok), r=NA_real_, p=NA_real_, method=NA_character_)
  })
})

# Adjust p-values BH
results <- results %>% mutate(p_adj = p.adjust(p, method="BH")) %>% arrange(p_adj, desc(abs(r)))

# Save results
fs::dir_create(out_dir)
write_tsv(results, file.path(out_dir, "pathogen_arg_correlations.tsv"))

# Optional heatmap of significant correlations
sig <- results %>% filter(!is.na(p_adj), p_adj <= 0.1)
if(nrow(sig) > 0){
  # Build matrix for heatmap
  mat <- sig %>% select(genus, arg_class, r) %>% pivot_wider(names_from = arg_class, values_from = r, values_fill = 0)
  mat_df <- as.data.frame(mat)
  rownames(mat_df) <- mat_df$genus; mat_df$genus <- NULL
  # Simple ggplot tile heatmap
  hm_df <- sig %>% mutate(sign = ifelse(r>0, "+", "-"))
  p <- ggplot(hm_df, aes(x = arg_class, y = genus, fill = r)) +
    geom_tile(color = "grey70") +
    scale_fill_gradient2(low = "#3b4cc0", mid = "#f7f7f7", high = "#b40426", midpoint = 0, name = "r") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
    labs(title = "Pathogen genus vs ARG class correlations", subtitle = "Cells shown pass FDR q <= 0.1", x = "ARG class", y = "Pathogen genus")
  fs::dir_create(fig_dir)
  ggsave(file.path(fig_dir, "pathogen_arg_correlation_heatmap.png"), p, width = 12, height = 10, dpi = 200)
}

message("[5/5] Done. Outputs:\n - ", file.path(out_dir, "pathogen_arg_correlations.tsv"), if(nrow(sig)>0) glue("\n - {file.path(fig_dir, 'pathogen_arg_correlation_heatmap.png')}") else "")
