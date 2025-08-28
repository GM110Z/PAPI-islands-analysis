tree_file    <- "New-tree-I-Hope.nwk"   # Newick with GCF_/GCA_ tip labels
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
ring_present_blue <- "#0000FF" # 1's
colname_angle     <- -5     # rotate dataset headers
colname_size      <- 2  # header size
colname_push_out  <- 30   # push headers outward from ring

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

# order dataset columns by prevalence
heat <- as.data.frame(mat2[value_cols], stringsAsFactors = FALSE)
rownames(heat) <- mat2$gcf
prev <- colSums(heat, na.rm = TRUE)
heat <- heat[, order(prev, decreasing = TRUE), drop = FALSE]
heat_disc <- as.data.frame(lapply(heat, function(v) factor(ifelse(v > 0, "1", "0"), levels = c("0","1"))))
rownames(heat_disc) <- rownames(heat)

# ===== Base circular tree WITH FORCED OPENING =====
p <- ggtree(tr, layout = "circular", size = tree_line_size)
p <- open_tree(p, open_angle)       # <-- force the wedge opening, works on all ggtree
p <- rotate_tree(p, rotate_degrees) # position the opening

# clean canvas
p <- p + theme_void() +
  theme(
    plot.margin      = margin(8, 8, 8, 8),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    legend.position  = "right"
  )

# ===== Tight ring (respects the opening) =====
p <- tryCatch(
  gheatmap(
    p, heat_disc,
    offset = ring_offset,
    width  = ring_width,
    colnames = TRUE,
    colnames_position = "top",
      colnames_angle = colname_angle,
    colnames_offset_y = colname_push_out,
    font.size = colname_size,
    hjust = 0.5,
    color = NA
  ),
  error = function(e) {
    gheatmap(
      p, heat_disc,
      offset = ring_offset, width = ring_width,
      colnames = TRUE, colnames_position = "top",
      font.size = colname_size, hjust = 0.5, color = NA
    ) + theme(axis.text.x = element_text(angle = colname_angle, vjust = 0.5, hjust = 0))
  }
) + scale_fill_manual(values = c("0" = ring_absent_grey, "1" = ring_present_blue),
                      na.value = "white", name = "Presence")

# ===== Save =====
ggsave("papi_tree_OPEN_forced_TIGHT_lightGreyRing_SUPERblue.pdf", p,
       width = 18, height = 18, units = "in", limitsize = FALSE)
