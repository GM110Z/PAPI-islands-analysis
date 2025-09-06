# ============================ CONFIG ============================
tree_file      <- "New-tree-I-Hope.nwk"   # Newick with GCF_/GCA_ tip labels
matrix_file    <- "matrix.tsv"            # TSV/CSV: first col = genome ID; others = 0/1
drop_version   <- TRUE                    # strip .version (e.g. GCF_XXXXXX.1 -> GCF_XXXXXX)

# ===== Layout & style =====
open_angle        <- 10
rotate_degrees    <- 90
tree_line_size    <- 0.15
radius_scale      <- 1.0
ring_offset       <- 0.0005
ring_width        <- 0.75
ring_absent_grey  <- "#E8EDF5"

# ============================ PACKAGES ============================
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

# ============================ HELPERS ============================
norm_id <- function(x, drop_ver = FALSE){
  x <- gsub("\\s+", "_", x)
  x <- gsub("^['\"]|['\"]$", "", x)
  x <- toupper(x)
  if (drop_ver) x <- sub("^(G[CA]F_\\d+)\\.\\d+$", "\\1", x)
  x
}

read_matrix_safely <- function(path){
  sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  mat <- read.delim(path, sep = sep, header = TRUE, check.names = FALSE,
                    quote = "", stringsAsFactors = FALSE)
  # detect ID column
  id_col <- names(mat)[tolower(names(mat)) %in% c("gcf","gca","genome","id","accession")]
  if (length(id_col) == 0) id_col <- names(mat)[1]
  names(mat)[names(mat) == id_col[1]] <- "gcf"
  # normalise IDs (drop_version global)
  mat$gcf <- norm_id(mat$gcf, drop_version)
  # coerce value cols to 0/1
  value_cols <- setdiff(names(mat), "gcf")
  mat[value_cols] <- lapply(mat[value_cols], function(v){
    v <- suppressWarnings(as.numeric(v)); v[is.na(v)] <- 0; as.integer(v > 0)
  })
  # aggregate duplicated genomes by max (safety)
  mat <- mat %>% group_by(gcf) %>% summarise(across(everything(), max), .groups = "drop")
  list(mat = mat, value_cols = setdiff(names(mat), "gcf"))
}

# ============================ TREE ============================
message(sprintf("Reading tree:   %s", tree_file))
tr_raw <- read.tree(tree_file)             # keep RAW labels (with versions) for plotting
lbl_raw  <- tr_raw$tip.label               # unique rownames for gheatmap (e.g., .1 vs .2)
lbl_norm <- norm_id(lbl_raw, drop_version) # normalized IDs (for joining only)

# Root/ladderize; keep edge lengths
tr <- phangorn::midpoint(tr_raw) |> ladderize(right = TRUE)
if (radius_scale != 1) tr$edge.length <- tr$edge.length * radius_scale

# Warn if multiple raw labels collapse to same normalized ID
if (drop_version) {
  dup_tbl <- tibble(raw = lbl_raw, dropped = lbl_norm) |>
    group_by(dropped) |>
    summarise(n = n(), raws = paste(sort(unique(raw)), collapse="; "), .groups="drop") |>
    filter(n > 1)
  if (nrow(dup_tbl) > 0) {
    message("[warn] Multiple tips collapse to same ID after dropping version:")
    print(dup_tbl, n = nrow(dup_tbl))
  }
}

# ============================ MATRIX ============================
message(sprintf("Reading matrix: %s", matrix_file))
mx <- read_matrix_safely(matrix_file)
mat <- mx$mat
value_cols <- mx$value_cols

# ============================ ALIGN / CHECKS ============================
# Join using normalized IDs but preserve RAW labels for plotting order
tip_df <- tibble(label_raw = lbl_raw, gcf = lbl_norm)
mat2 <- left_join(tip_df, mat, by = "gcf")
mat2[value_cols] <- lapply(mat2[value_cols], function(v){ v[is.na(v)] <- 0L; v })

# Hard-stop diagnostics
matched <- sum(lbl_norm %in% mat$gcf)
message(sprintf("[check] drop_version=%s", as.character(drop_version)))
message(sprintf("[check] tips matched to matrix IDs: %d / %d", matched, length(lbl_norm)))
if (matched < 0.95 * length(lbl_norm)) {
  stop("Fewer than 95% of tips matched matrix IDs — check drop_version and that tree/matrix files are intended.")
}
col_nonzero <- colSums(as.data.frame(mat2[value_cols], stringsAsFactors = FALSE))
nz <- sum(col_nonzero > 0)
message(sprintf("[check] islands with any presence: %d / %d", nz, length(col_nonzero)))
if (nz == 0) stop("All-zero ring — join succeeded but all values are zero. Check island columns in the matrix.")

# Order columns by prevalence
heat <- as.data.frame(mat2[value_cols], stringsAsFactors = FALSE)
rownames(heat) <- tip_df$label_raw   # RAW labels keep rownames unique
prev <- colSums(heat, na.rm = TRUE)
heat <- heat[, order(prev, decreasing = TRUE), drop = FALSE]

# ============================ RING KEY & COLORS ============================
# Build categorical matrix for legend coloring
heat_key <- heat
for (j in seq_along(heat_key)) {
  nm <- colnames(heat_key)[j]
  heat_key[[j]] <- ifelse(heat_key[[j]] > 0, nm, "0")
}
# Consistent levels across columns
heat_key[] <- lapply(heat_key, function(col) factor(col, levels = c("0", colnames(heat_key))))

present_pal <- setNames(grDevices::hcl.colors(n = ncol(heat_key), palette = "Spectral"),
                        colnames(heat_key))
fill_vals <- c("0" = ring_absent_grey, present_pal)

# ============================ PLOT ============================
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

# Tight ring (legend-only column names)
p <- tryCatch(
  gheatmap(
    p, heat_key,
    offset = ring_offset,
    width  = ring_width,
    colnames = FALSE,
    color = NA
  ),
  error = function(e) {
    gheatmap(
      p, heat_key,
      offset = ring_offset, width = ring_width,
      colnames = FALSE, color = NA
    )
  }
) + scale_fill_manual(
  values = fill_vals,
  name   = "Island (present)",
  breaks = colnames(heat_key)  # islands only (omit "Absent" grey)
  # To include "Absent" in legend: breaks = c("0", colnames(heat_key))
)

# ============================ SAVE ============================
out_pdf <- "papi_tree_OPEN_forced_TIGHT_multicolorLegend_fixed.pdf"
ggsave(out_pdf, p, width = 18, height = 18, units = "in", limitsize = FALSE)
message(sprintf("[ok] wrote %s", out_pdf))

