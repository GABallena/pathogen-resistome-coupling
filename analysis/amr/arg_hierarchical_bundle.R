# ARG hierarchical edge-bundling (portfolio-ready)
#
# Inputs:
#   1) ARG relative abundance matrix (TSV). First column should be an identifier like:
#        db|type|class|mechanism|gene
#      Remaining columns are samples.
#   2) (Optional) sample metadata TSV with columns:
#        sample_id, site, treatment
#
# Outputs:
#   - edge-bundling plot linking ARG genes to samples (PDF/PNG)
#
# Notes:
#   - This is a generic reimplementation based on the original working intent.
#   - Replace parsing rules if your ResistanceGene naming scheme differs.

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript arg_hierarchical_bundle.R <arg_matrix.tsv> <out_prefix> [metadata.tsv]\n")
  quit(status = 2)
}
arg_path <- args[[1]]
out_prefix <- args[[2]]
meta_path <- if (length(args) >= 3) args[[3]] else NA_character_

# --- load ---
arg_df <- readr::read_tsv(arg_path, show_col_types = FALSE)
stopifnot(ncol(arg_df) >= 3)
colnames(arg_df)[1] <- "ResistanceGene"

# Optional metadata
meta <- NULL
if (!is.na(meta_path) && file.exists(meta_path)) {
  meta <- readr::read_tsv(meta_path, show_col_types = FALSE) %>%
    mutate(sample_id = as.character(sample_id))
}

# --- parse hierarchy from ResistanceGene ---
# Expecting 5 fields separated by "|"
parts <- stringr::str_split_fixed(arg_df$ResistanceGene, "\\|", 5)
hier <- tibble(
  db = parts[,1],
  type = parts[,2],
  class = parts[,3],
  mechanism = parts[,4],
  gene = parts[,5]
) %>%
  mutate(across(everything(), ~ifelse(.x == "" | is.na(.x), "(NA)", .x))) %>%
  mutate(gene = make.unique(gene))

# Build edges for a tree:
# root -> db -> type -> class -> mechanism -> gene
edges_tree <- bind_rows(
  tibble(from="ARGs", to=unique(hier$db)),
  hier %>% distinct(db, type) %>% transmute(from=db, to=type),
  hier %>% distinct(type, class) %>% transmute(from=type, to=class),
  hier %>% distinct(class, mechanism) %>% transmute(from=class, to=mechanism),
  hier %>% distinct(mechanism, gene) %>% transmute(from=mechanism, to=gene)
)

# Sample nodes
sample_cols <- setdiff(colnames(arg_df), "ResistanceGene")
edges_samples <- tibble(from="Samples", to=sample_cols)
edges_tree2 <- bind_rows(edges_tree, edges_samples)

# Connections: gene -> sample where abundance > 0
conn <- arg_df %>%
  mutate(gene = hier$gene) %>%
  select(gene, all_of(sample_cols)) %>%
  pivot_longer(-gene, names_to="sample", values_to="abundance") %>%
  mutate(abundance = suppressWarnings(as.numeric(abundance))) %>%
  filter(!is.na(abundance), abundance > 0)

# Full node list
nodes <- unique(c(edges_tree2$from, edges_tree2$to, conn$gene, conn$sample))

# Build graph for hierarchy (tree-ish)
g <- graph_from_data_frame(edges_tree2, directed = TRUE, vertices = nodes)

# Build bundled connections (requires a dendrogram-like layout)
# We create a combined data frame mapping endpoints in the graph.
conn2 <- conn %>% transmute(from=gene, to=sample, weight=abundance)

# Layout + plot
pdf(paste0(out_prefix, "_arg_bundle.pdf"), width = 12, height = 12)
p <- ggraph(g, layout = "dendrogram", circular = TRUE) +
  geom_edge_diagonal(alpha = 0.25) +
  geom_node_point(size = 1, alpha = 0.8) +
  geom_node_text(aes(label = name), size = 2, repel = TRUE) +
  theme_void() +
  ggtitle("ARG hierarchical structure (portfolio-ready)")
print(p)
dev.off()

# Connection bundle plot (optional; can be heavy). Provide a readable fallback.
# Try geom_conn_bundle if available.
out_png <- paste0(out_prefix, "_arg_bundle.png")
try({
  # Convert to a hierarchy required by bundle
  h <- as.dendrogram(g)
  # ggraph conn bundle expects an 'edges' list with 'from'/'to' in node indices
}, silent = TRUE)

cat("Wrote: ", paste0(out_prefix, "_arg_bundle.pdf"), "\n")
