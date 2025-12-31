#!/usr/bin/env Rscript

# ARGs-OAP per-site Sankey (alluvial) plots with Top 10 ARG tables/plots
# - Normalizes counts per GB (and per L when available)
# - Outputs per-site + global CSVs
# - Generates Sankey plots (Untreated -> Treated)
# - Extracts Top 10 ARGs per site and globally
# - Produces Top 10 barplots

suppressPackageStartupMessages({
  req_pkgs <- c("data.table","dplyr","stringr","readr","tools",
                "ggplot2","ggalluvial","RColorBrewer","scales",
                "ggrepel","patchwork")
  to_install <- req_pkgs[!sapply(req_pkgs, function(p) requireNamespace(p, quietly = TRUE))]
  if (length(to_install) > 0) {
    try(suppressWarnings(install.packages(to_install, repos="https://cloud.r-project.org")), silent=TRUE)
  }
  library(data.table)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tools)
  library(ggplot2)
  library(ggalluvial)
  library(RColorBrewer)
  # scales, ggrepel, patchwork loaded
})

root <- "."
samples_dir <- file.path(root,"arg_oap_work","samples")
masterlist_path <- file.path(root,"metadata/masterlist.tsv")
out_dir <- file.path(root,"diversity_results","arg_sankey")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
trimmed_dir <- file.path(root,"trimmed_reads")

# ---- Visualization config ----
TOP_N <- as.integer(Sys.getenv("ARG_SANKEY_TOP_N", unset = "15"))
SHOW_LEGEND_MAX <- as.integer(Sys.getenv("ARG_SANKEY_SHOW_LEGEND_MAX", unset = "25"))
SAVE_PDF <- tolower(Sys.getenv("ARG_SANKEY_SAVE_PDF", unset = "true")) %in% c("1","true","yes")
DELTA_TOP_N <- as.integer(Sys.getenv("ARG_SANKEY_DELTA_TOP_N", unset = "12"))
USE_COLORBLIND <- tolower(Sys.getenv("ARG_SANKEY_USE_COLORBLIND", unset = "false")) %in% c("1","true","yes")

# Deterministic color mapping for ARG types across sites
arg_type_colors <- function(types) {
  types <- unique(types)
  # coerce to character and replace NAs
  types_chr <- as.character(types)
  types_chr[is.na(types_chr)] <- "(NA)"
  if (USE_COLORBLIND && length(types_chr) <= 8) {
    base <- c("#0072B2","#E69F00","#009E73","#D55E00","#CC79A7","#56B4E9","#F0E442","#999999")
    cols <- setNames(rep(base, length.out = length(types_chr)), types_chr)
    cols["Other"] <- "#bdbdbd"
  } else {
    cols <- vapply(types_chr, function(t) {
      if (is.na(t) || t == "(NA)") return("#999999")
      if (t == "Other") return("#bdbdbd")
      h <- (sum(as.integer(charToRaw(as.character(t)))) %% 360)
      grDevices::hcl(h = h, c = 60, l = 65)
    }, character(1))
  }
  structure(cols, names = types_chr)
}

collapse_topN <- function(long_df, top_n = 15) {
  agg <- long_df %>% group_by(ARG_type) %>% summarise(tot = sum(count, na.rm = TRUE), .groups = "drop")
  top_types <- agg %>% arrange(desc(tot)) %>% slice_head(n = top_n) %>% pull(ARG_type)
  long_df %>% mutate(ARG_type = ifelse(ARG_type %in% top_types, ARG_type, "Other")) %>%
    group_by(state, ARG_type) %>% summarise(count = sum(count, na.rm = TRUE), .groups = "drop")
}

# Order ARG_type factor consistently by overall mass (Other last)
order_arg_types <- function(long_df) {
  ord <- long_df %>% group_by(ARG_type) %>% summarise(tot = sum(count, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(tot)) %>% pull(ARG_type)
  ord <- c(setdiff(ord, "Other"), "Other")
  long_df %>% mutate(ARG_type = factor(ARG_type, levels = ord))
}

# Helper: find the sample directory for a sample code
find_sample_dir <- function(code) {
  dirs <- list.dirs(samples_dir, full.names=TRUE, recursive=FALSE)
  m <- dirs[grepl(paste0("^", code,"(?:_|$)"), basename(dirs), perl=TRUE) & 
              !grepl("\\.part_", basename(dirs))]
  if (length(m)==0) return(NA_character_)
  m[which.min(nchar(m))]
}

# Helper: read unnormalized type counts
read_arg_type_counts <- function(sample_dir) {
  f <- file.path(sample_dir,"unnormalized_count.type.txt")
  if (!file.exists(f)) return(NULL)
  dt <- suppressWarnings(fread(f, sep="\t", header=TRUE))
  if (!nrow(dt)) return(NULL)
  num_cols <- setdiff(names(dt), names(dt)[1])
  num_cols <- num_cols[sapply(num_cols,function(x) is.numeric(dt[[x]]))]
  if (length(num_cols)==0) return(NULL)
  sum_count <- rowSums(as.data.frame(dt[,num_cols,with=FALSE]), na.rm=TRUE)
  tibble(ARG_type=dt[[1]], count=sum_count)
}

# Helper: get total trimmed read size (R1+R2) in GB
get_sample_size_gb <- function(trim_prefix) {
  r1 <- file.path(trimmed_dir, paste0(trim_prefix,"_R1_trimmed.fq.gz"))
  r2 <- file.path(trimmed_dir, paste0(trim_prefix,"_R2_trimmed.fq.gz"))
  sizes <- c(0,0)
  if (file.exists(r1)) sizes[1] <- file.info(r1)$size
  if (file.exists(r2)) sizes[2] <- file.info(r2)$size
  total_bytes <- sum(sizes, na.rm=TRUE)
  if (is.na(total_bytes) || total_bytes<=0) return(NA_real_)
  total_bytes/1e9
}

# --- Parse master list ---
ml <- suppressWarnings(read_tsv(masterlist_path, show_col_types=FALSE, col_names=TRUE))
names(ml) <- str_trim(names(ml))
stopifnot(all(c("SAMPLE CODE","SAMPLE DESCRIPTION","SAMPLE TYPE") %in% names(ml)))

map <- ml %>%
  transmute(sample_code=str_trim(`SAMPLE CODE`),
            description=str_trim(`SAMPLE DESCRIPTION`),
            sample_type=str_trim(`SAMPLE TYPE`)) %>%
  filter(!is.na(sample_code), !is.na(description), !is.na(sample_type)) %>%
  mutate(
    site = description %>%
      str_replace("\\s+Untreated.*$","") %>%
      str_replace("\\s+Treated.*$",""),
    state = case_when(
      str_detect(str_to_lower(description),"untreated") ~ "Untreated",
      str_detect(str_to_lower(description),"treated") ~ "Treated",
      TRUE ~ NA_character_
    ),
    replicate_hint = case_when(
      str_detect(description,"\\b2\\b") ~ 2L,
      str_detect(description,"\\b1\\b") ~ 1L,
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(!is.na(state))

sites_with_pairs <- map %>%
  group_by(site) %>%
  summarise(has_treated=any(state=="Treated"),
            has_untreated=any(state=="Untreated"), .groups="drop") %>%
  filter(has_treated & has_untreated) %>%
  pull(site)

map_pairs_resolved <- map %>%
  filter(site %in% sites_with_pairs) %>%
  group_by(site,state) %>%
  arrange(desc(coalesce(replicate_hint,0L))) %>%
  slice_head(n=1) %>%
  ungroup() %>%
  mutate(sample_dir=vapply(sample_code, find_sample_dir, character(1))) %>%
  filter(!is.na(sample_dir) & file.exists(file.path(sample_dir,"unnormalized_count.type.txt")))

if (nrow(map_pairs_resolved)==0) stop("No site pairs with available ARG counts.")

message("Sites with pairs:")
print(map_pairs_resolved %>% select(site,state,sample_code,sample_dir))

all_site_summaries <- list()
top10_site_summaries <- list()

# --- Main per-site loop ---
for (site_name in unique(map_pairs_resolved$site)) {
  sel <- map_pairs_resolved %>% filter(site==site_name)
  if (!all(c("Treated","Untreated") %in% sel$state)) next
  s_unt <- sel %>% filter(state=="Untreated") %>% slice(1)
  s_trt <- sel %>% filter(state=="Treated") %>% slice(1)

  dt_unt <- read_arg_type_counts(s_unt$sample_dir)
  dt_trt <- read_arg_type_counts(s_trt$sample_dir)
  if (is.null(dt_unt) || is.null(dt_trt)) next

  size_gb_unt <- get_sample_size_gb(basename(s_unt$sample_dir))
  size_gb_trt <- get_sample_size_gb(basename(s_trt$sample_dir))

  sum_df <- full_join(dt_unt, dt_trt, by="ARG_type", suffix=c("_unt","_trt")) %>%
    mutate(
      untreated_count_raw = coalesce(count_unt,0),
      treated_count_raw   = coalesce(count_trt,0),
      untreated_count_per_GB = if (!is.na(size_gb_unt) && size_gb_unt>0) untreated_count_raw/size_gb_unt else NA_real_,
      treated_count_per_GB   = if (!is.na(size_gb_trt) && size_gb_trt>0) treated_count_raw/size_gb_trt else NA_real_,
  # 1.0 L for untreated, 3.0 L for treated
  untreated_count_per_GB_per_L = ifelse(!is.na(untreated_count_per_GB), untreated_count_per_GB/1.0, NA_real_),
      treated_count_per_GB_per_L   = ifelse(!is.na(treated_count_per_GB), treated_count_per_GB/3.0, NA_real_)
    ) %>%
    select(ARG_type, untreated_count_raw, treated_count_raw,
           untreated_count_per_GB, treated_count_per_GB,
           untreated_count_per_GB_per_L, treated_count_per_GB_per_L) %>%
    arrange(desc(coalesce(untreated_count_per_GB_per_L,0) + coalesce(treated_count_per_GB_per_L,0)))

  # --- Top 10 ARGs ---
  sum_df <- sum_df %>%
    mutate(total_abundance = rowSums(select(., matches("untreated|treated")), na.rm=TRUE))
  top10 <- sum_df %>% arrange(desc(total_abundance)) %>% slice_head(n=10)
  write_csv(top10, file.path(out_dir, paste0("top10_ARGs_",make.names(site_name),".csv")))

  # Top 10 barplot
  palette <- setNames(colorRampPalette(brewer.pal(12,"Set3"))(nrow(top10)), top10$ARG_type)
  p_top10 <- ggplot(top10, aes(x=reorder(ARG_type,total_abundance),
                               y=total_abundance, fill=ARG_type)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values=palette) +
    labs(title=paste("Top 10 ARG types -",site_name),
         y="Total abundance (norm or raw)", x=NULL) +
    theme_minimal(base_size=12) +
    theme(legend.position="none")
  ggsave(file.path(out_dir,paste0("top10_ARGs_",make.names(site_name),".png")),
         p_top10, width=7, height=5, dpi=300)

  # --- Alluvial plot (normalized if possible) + raw ---
  use_per_gbperl <- all(!is.na(sum_df$untreated_count_per_GB_per_L)) && all(!is.na(sum_df$treated_count_per_GB_per_L))
  use_per_gb     <- all(!is.na(sum_df$untreated_count_per_GB)) && all(!is.na(sum_df$treated_count_per_GB))
  long_df_norm <- if (use_per_gbperl) {
    bind_rows(sum_df %>% transmute(state="Untreated", ARG_type, count=untreated_count_per_GB_per_L),
              sum_df %>% transmute(state="Treated", ARG_type, count=treated_count_per_GB_per_L))
  } else if (use_per_gb) {
    bind_rows(sum_df %>% transmute(state="Untreated", ARG_type, count=untreated_count_per_GB),
              sum_df %>% transmute(state="Treated", ARG_type, count=treated_count_per_GB))
  } else {
    bind_rows(sum_df %>% transmute(state="Untreated", ARG_type, count=untreated_count_raw),
              sum_df %>% transmute(state="Treated", ARG_type, count=treated_count_raw))
  }
  long_df_raw <- bind_rows(
    sum_df %>% transmute(state="Untreated", ARG_type, count=untreated_count_raw),
    sum_df %>% transmute(state="Treated", ARG_type, count=treated_count_raw)
  )
  # Collapse to Top-N + Other, filter zeros, enforce axis order and ARG order
  long_df_norm <- long_df_norm %>% collapse_topN(TOP_N) %>% filter(count>0) %>%
    mutate(state=factor(state, levels=c("Untreated","Treated"))) %>%
    order_arg_types()
  long_df_raw  <- long_df_raw  %>% collapse_topN(TOP_N) %>% filter(count>0) %>%
    mutate(state=factor(state, levels=c("Untreated","Treated"))) %>%
    order_arg_types()

  if (nrow(long_df_norm)>0) {
    arg_types <- unique(long_df_norm$ARG_type)
    palette <- arg_type_colors(arg_types)
    p_norm <- ggplot(long_df_norm,
                     aes(x=state, stratum=ARG_type, alluvium=ARG_type,
                         y=count, fill=ARG_type, label=ARG_type)) +
      geom_flow(stat="alluvium", lode.guidance="frontback", alpha=0.7) +
      geom_stratum(width=0.32, color="grey20") +
      scale_fill_manual(values=palette) +
      scale_y_continuous(labels=scales::comma) +
      scale_x_discrete(limits=c("Untreated","Treated")) +
      labs(title=paste0("ARG type flow (",site_name,")"),
           subtitle=paste0(s_unt$sample_code," (Untreated) → ",s_trt$sample_code," (Treated)"),
           y=if (use_per_gbperl) "Counts per GB per L" else if (use_per_gb) "Counts per GB" else "Counts",
           x=NULL, fill="ARG type") +
      theme_minimal(base_size=13) +
      theme(legend.position=if (length(arg_types) <= SHOW_LEGEND_MAX) "right" else "none",
            legend.key.size=unit(0.4,"cm"))
    out_norm <- file.path(out_dir, paste0("arg_type_alluvial_",make.names(site_name),".png"))
    ggsave(out_norm, p_norm, width=12, height=7, dpi=300, bg="white")
    if (SAVE_PDF) ggsave(sub("\\.png$",".pdf", out_norm), p_norm, width=12, height=7)
  }

  # Raw plot as well
  if (nrow(long_df_raw)>0) {
    arg_types_r <- unique(long_df_raw$ARG_type)
    palette_r <- arg_type_colors(arg_types_r)
    p_raw <- ggplot(long_df_raw,
                    aes(x=state, stratum=ARG_type, alluvium=ARG_type,
                        y=count, fill=ARG_type, label=ARG_type)) +
      geom_flow(stat="alluvium", lode.guidance="frontback", alpha=0.7) +
      geom_stratum(width=0.32, color="grey20") +
      scale_fill_manual(values=palette_r) +
      scale_y_continuous(labels=scales::comma) +
      scale_x_discrete(limits=c("Untreated","Treated")) +
      labs(title=paste0("ARG type raw counts (",site_name,")"),
           subtitle=paste0(s_unt$sample_code," (Untreated) → ",s_trt$sample_code," (Treated)"),
           y="Counts", x=NULL, fill="ARG type") +
      theme_minimal(base_size=13) +
      theme(legend.position=if (length(arg_types_r) <= SHOW_LEGEND_MAX) "right" else "none",
            legend.key.size=unit(0.4,"cm"))
    out_raw <- file.path(out_dir, paste0("arg_type_alluvial_",make.names(site_name),"_raw.png"))
    ggsave(out_raw, p_raw, width=12, height=7, dpi=300, bg="white")
    if (SAVE_PDF) ggsave(sub("\\.png$",".pdf", out_raw), p_raw, width=12, height=7)
  }

  # Collect summaries
  all_site_summaries[[site_name]] <- sum_df %>%
    mutate(site=site_name,
           untreated_sample_code=s_unt$sample_code, treated_sample_code=s_trt$sample_code,
           untreated_size_GB=size_gb_unt, treated_size_GB=size_gb_trt,
            untreated_volume_L=1.0, treated_volume_L=3.0)
  top10_site_summaries[[site_name]] <- top10
}

# --- Global summaries ---
if (length(all_site_summaries)>0) {
  combined <- bind_rows(all_site_summaries) %>% relocate(site, untreated_sample_code, treated_sample_code)
  write_csv(combined, file.path(out_dir,"arg_type_counts_all_sites.csv"))

  top10_global <- combined %>%
    group_by(ARG_type) %>%
    summarise(total_abundance=sum(total_abundance, na.rm=TRUE), .groups="drop") %>%
    arrange(desc(total_abundance)) %>%
    slice_head(n=10)
  write_csv(top10_global, file.path(out_dir,"top10_ARGs_global.csv"))

  # Global Top 10 barplot
  palette <- setNames(colorRampPalette(brewer.pal(12,"Set3"))(nrow(top10_global)), top10_global$ARG_type)
  p_global <- ggplot(top10_global, aes(x=reorder(ARG_type,total_abundance),
                                       y=total_abundance, fill=ARG_type)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values=palette) +
    labs(title="Global Top 10 ARG types", y="Total abundance", x=NULL) +
    theme_minimal(base_size=12) + theme(legend.position="none")
  ggsave(file.path(out_dir,"top10_ARGs_global.png"), p_global, width=7, height=5, dpi=300)

  message("Wrote per-site + global summaries, Sankey plots, and Top10 outputs to ", out_dir)
} else {
  message("No results generated.")
}


# Long format with negative for untreated
top10_global_long <- combined %>%
  group_by(ARG_type) %>%
  summarise(
    untreated = sum(untreated_count_per_GB_per_L, na.rm=TRUE),
    treated   = sum(treated_count_per_GB_per_L, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(total = untreated + treated) %>%
  arrange(desc(total)) %>%
  slice_head(n=10) %>%
  tidyr::pivot_longer(cols=c("untreated","treated"), names_to="state", values_to="count") %>%
  mutate(count = ifelse(state=="untreated", -count, count))
palette2 <- c("untreated"="#d95f02", "treated"="#1b9e77")  # orange vs teal

p_div <- ggplot(top10_global_long,
       aes(x=reorder(ARG_type, abs(count), FUN=max), y=count, fill=state)) +
  geom_col(width=0.7) +
  coord_flip() +
  scale_fill_manual(values=palette2) +
  scale_y_continuous(labels=function(x) scales::comma(abs(x))) +
  labs(
    title="WWTPs Top 10 ARGs (Untreated vs Treated)",
    y="Total normalized abundance (counts per GB per L)",
    x=NULL, fill=NULL
  ) +
  theme_minimal(base_size=13) +
  theme(panel.grid.major.y=element_blank(),
        panel.grid.minor=element_blank())

ggsave(file.path(out_dir,"top10_ARGs_global_diverging.png"),
       p_div, width=8, height=6, dpi=300, bg="white")


# Optional: multi-site grid of normalized Sankeys (if patchwork is available)
if (exists("long_df_norm") && requireNamespace("patchwork", quietly = TRUE)) {
  # Rebuild site-level normalized plots and patch them together
  site_plots <- list()
  for (site_name in unique(map_pairs_resolved$site)) {
    sel <- map_pairs_resolved %>% filter(site==site_name)
    if (!all(c("Treated","Untreated") %in% sel$state)) next
    s_unt <- sel %>% filter(state=="Untreated") %>% slice(1)
    s_trt <- sel %>% filter(state=="Treated") %>% slice(1)
    # read cached per-site CSV to avoid recomputation
    csv <- file.path(out_dir, paste0("arg_type_counts_",make.names(site_name),".csv"))
    if (!file.exists(csv)) next
    df <- suppressWarnings(read_csv(csv, show_col_types = FALSE))
    if (!all(c("untreated_count_per_GB_per_L","treated_count_per_GB_per_L") %in% names(df))) next
    ldf <- bind_rows(
      df %>% transmute(state="Untreated", ARG_type, count=untreated_count_per_GB_per_L),
      df %>% transmute(state="Treated",   ARG_type, count=treated_count_per_GB_per_L)
    ) %>% filter(!is.na(count) & count>0) %>% collapse_topN(TOP_N) %>% order_arg_types() %>%
      mutate(state=factor(state, levels=c("Untreated","Treated")))
    if (nrow(ldf)==0) next
    pal <- arg_type_colors(unique(ldf$ARG_type))
    p <- ggplot(ldf, aes(x=state, stratum=ARG_type, alluvium=ARG_type, y=count, fill=ARG_type))+
      geom_flow(stat="alluvium", alpha=0.65) + geom_stratum(width=0.3,color="grey25")+
      scale_fill_manual(values=pal)+ scale_x_discrete(limits=c("Untreated","Treated"))+ theme_void()+
      labs(title=site_name)
    site_plots[[site_name]] <- p
  }
  if (length(site_plots)>0) {
    grid <- Reduce(`+`, site_plots) + patchwork::plot_layout(ncol = 2, guides = "collect")
    ggsave(file.path(out_dir, "arg_type_alluvial_sites_grid.png"), grid, width=14, height=10, dpi=300, bg="white")
  }
}
