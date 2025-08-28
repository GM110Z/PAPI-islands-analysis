tree_file    <- "tree.nwk"   # Newick with GCF_/GCA_ tip labels
matrix_file  <- "presence.tsv"        # TSV/CSV: first col = genome ID; others = 0/1
drop_version <- FALSE                 # TRUE if IDs look like GCF_XXXXXX (no ".1")

# ===== Layout & style =====
open_angle        <- 10    # degrees of opening wedge (try 160–220)
rotate_degrees    <- 90     # rotate entire plot to place the opening
tree_line_size    <- 0.15   # thicker branches → broader looking tree
radius_scale      <- 1.0    # >1 stretches tree radially (e.g., 1.5 or 2.0)
ring_offset       <- 0.0005 # tiny gap from leaf tips to ring (smaller = closer)
ring_width        <- 0.75   # ring thickness
ring_absent_grey  <- "#E8EDF5" # 0's
colname_angle     <- -7     # rotate dataset headers (unused now; we hide colnames)
colname_size      <- 2      # header size (unused now)
colname_push_out  <- 10     # push headers outward from ring (unused now)

# ===== Packages =====
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  for (pkg in c("ape","ggplot2","dplyr","tibble","readr","tidyr"))
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg, repos = "https://cloud.r-project.org")
  if (!requireNamespace("ggtree",   quietly = TRUE)) BiocManager::install("ggtree",   update = FALSE, ask = FALSE)
  if (!requireNamespace("phangorn", quietly = TRUE)) BiocManager::install("phangorn", update = FALSE, ask = FALSE)
})
library(ape); library(ggtree); library(ggplot2)
library(dplyr); library(tibble); library(readr); library(tidyr); library(phangorn)

# ===== Helpers =====
norm_id <- function(x, drop_ver = FALSE){
  x <- gsub("\\s+", "_", x)
  x <- gsub("^['\"]|['\"]$", "", x)
  x <- toupper(x)
  if (drop_ver) x <- sub("^(G[CA]F_\\d+)\\.\\d+$", "\\1", x)
  x
}

# ===== Tree =====
tr <- read.tree(tree_file)
tr$tip.label <- norm_id(tr$tip.label, drop_version)
tr <- phangorn::midpoint(tr) |> ladderize(right = TRUE)
if (radius_scale != 1) tr$edge.length <- tr$edge.length * radius_scale

# ===== Matrix =====
sep <- if (grepl("\\.csv$", matrix_file, ignore.case = TRUE)) "," else "\t"
mat <- read.delim(matrix_file, sep = sep, header = TRUE, check.names = FALSE,
                  quote = "", stringsAsFactors = FALSE)
id_col <- names(mat)[tolower(names(mat)) %in% c("gcf","gca","genome","id","accession")][1]
stopifnot(length(id_col) == 1)
mat <- mat %>% rename(gcf = all_of(id_col)) %>% mutate(gcf = norm_id(gcf, drop_version))
value_cols <- setdiff(names(mat), "gcf"); stopifnot(length(value_cols) > 0)
mat[value_cols] <- lapply(mat[value_cols], function(v){
  v <- suppressWarnings(as.numeric(v)); v[is.na(v)] <- 0; as.integer(v > 0)
})
tip_df <- tibble(gcf = tr$tip.label)
mat2 <- left_join(tip_df, mat, by = "gcf")
mat2[value_cols] <- lapply(mat2[value_cols], function(v){ v[is.na(v)] <- 0L; v })

# Order columns by prevalence
heat <- as.data.frame(mat2[value_cols], stringsAsFactors = FALSE)
rownames(heat) <- mat2$gcf
prev <- colSums(heat, na.rm = TRUE)
heat <- heat[, order(prev, decreasing = TRUE), drop = FALSE]

# === Option A: one color per island (present), single grey for absent ===
# Build a key matrix: "0" for absent; "<island_name>" for present
heat_key <- heat
for (j in seq_along(heat_key)) {
  nm <- colnames(heat_key)[j]
  heat_key[[j]] <- ifelse(heat_key[[j]] > 0, nm, "0")
}
# Ensure consistent levels across columns (so legend shows all islands)
heat_key[] <- lapply(heat_key, factor, levels = c("0", colnames(heat_key)))

# Color palette: one color per island + grey for "absent"
present_pal <- setNames(grDevices::hcl.colors(n = ncol(heat_key), palette = "Spectral"),
                        colnames(heat_key))
fill_vals <- c("0" = ring_absent_grey, present_pal)

# ===== Base circular tree WITH FORCED OPENING =====
p <- ggtree(tr, layout = "circular", size = tree_line_size)
p <- open_tree(p, open_angle)
p <- rotate_tree(p, rotate_degrees)

# Clean canvas
p <- p + theme_void() +
  theme(
    plot.margin      = margin(8, 8, 8, 8),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    legend.position  = "right"
  )

# ===== Tight ring (no column headers; legend carries names) =====
p <- tryCatch(
  gheatmap(
    p, heat_key,
    offset = ring_offset,
    width  = ring_width,
    colnames = FALSE,   # <-- hide ring labels; legend only
    color = NA
  ),
  error = function(e) {
    gheatmap(
      p, heat_key,
      offset = ring_offset, width = ring_width,
      colnames = FALSE,
      color = NA
    )
  }
) + scale_fill_manual(
  values = fill_vals,
  name   = "Island (present)",
  breaks = colnames(heat_key)  # exclude "0" so legend shows only islands
  # To include "Absent", use: breaks = c("0", colnames(heat_key))
)

# ===== Save =====
ggsave("papi_tree_OPEN_forced_TIGHT_multicolorLegendBerlin.pdf", p,
       width = 18, height = 18, units = "in", limitsize = FALSE)
