#!/usr/bin/env Rscript
# Minimal joint per-cell filtering:
#   - Read ATAC QC TSV (must have `cell`; may have `pass`)
#   - Read GEX QC  TSV (must have `cell`; may have doublet column; may have `pass`)
#   - Harmonize GEX cell IDs (optional: replace second "_" with "#")
#   - Define keep sets and write a single TSV with keep flags
#
# CLI:
#   --atac_qc      path/to/atac_qc.tsv
#   --gex_qc       path/to/gex_qc.tsv
#   --out_tsv      path/to/cell_filters.tsv
#   --gex_fix_second_underscore  TRUE/FALSE (default TRUE)
#   --drop_gex_doublets          TRUE/FALSE (default FALSE)
#   --gex_doublet_col            default "scDblFinder_class"
#   --gex_doublet_label          default "doublet"

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt_list <- list(
  make_option("--atac_qc", type="character"),
  make_option("--gex_qc",  type="character"),
  make_option("--out_tsv", type="character"),
  make_option("--gex_fix_second_underscore", type="logical", default=TRUE),
  make_option("--drop_gex_doublets", type="logical", default=FALSE),
  make_option("--gex_doublet_col", type="character", default="scDblFinder_class"),
  make_option("--gex_doublet_label", type="character", default="doublet")
)
opt <- parse_args(OptionParser(option_list=opt_list))

dir.create(dirname(opt$out_tsv), recursive = TRUE, showWarnings = FALSE)

# ---------- helpers ----------
fix_second_underscore <- function(v) {
  # replace 2nd "_" with "#": "SAMPLE_BARCODE_X" -> "SAMPLE#BARCODE_X"
  sub("^(.*?_.*?)_(.*)$", "\\1#\\2", v, perl = TRUE)
}
as_int01 <- function(x) as.integer(ifelse(x, 1L, 0L))

# ---------- load ----------
atac <- fread(opt$atac_qc, sep = "\t", header = TRUE)
gex  <- fread(opt$gex_qc,  sep = "\t", header = TRUE)

stopifnot("cell" %in% names(atac), "cell" %in% names(gex))

# harmonize GEX cell names if requested
if (isTRUE(opt$gex_fix_second_underscore)) {
  gex$cell <- fix_second_underscore(gex$cell)
}

# define pass flags
atac$pass <- if ("pass" %in% names(atac)) as.integer(atac$pass != 0) else 1L
gex$.__not_doublet__ <- 1L
if (isTRUE(opt$drop_gex_doublets) && opt$gex_doublet_col %in% names(gex)) {
  gex$.__not_doublet__ <- as.integer(gex[[opt$gex_doublet_col]] != opt$gex_doublet_label)
}

gex$pass <- if ("pass" %in% names(gex)) as.integer(gex$pass != 0) else 1L
gex$pass <- as.integer(gex$pass & gex$.__not_doublet__)

# sets
atac_pass_cells <- atac$cell[atac$pass == 1L]
gex_pass_cells  <- gex$cell[gex$pass == 1L]

keep_cells <- intersect(atac_pass_cells, gex_pass_cells)

# annotate and write
dt <- data.table(
  cell = unique(c(atac$cell, gex$cell))
)
dt[, in_atac := as_int01(cell %in% atac$cell)]
dt[, atac_pass := as_int01(cell %in% atac_pass_cells)]
dt[, in_gex := as_int01(cell %in% gex$cell)]
dt[, gex_pass := as_int01(cell %in% gex_pass_cells)]
dt[, keep_joint := as_int01(cell %in% keep_cells)]

# optional: carry Sample if available
if ("Sample" %in% names(atac)) dt[match(atac$cell, dt$cell), Sample_atac := atac$Sample]
if ("Sample" %in% names(gex))  dt[match(gex$cell,  dt$cell), Sample_gex  := gex$Sample]

setcolorder(dt, c("cell","keep_joint","in_atac","atac_pass","in_gex","gex_pass",
                  setdiff(names(dt), c("cell","keep_joint","in_atac","atac_pass","in_gex","gex_pass"))))
fwrite(dt, opt$out_tsv, sep = "\t", quote = FALSE, na = "NA")

message(sprintf("Joint filter: kept %d cells (ATAC pass %d ∩ GEX pass %d).",
                sum(dt$keep_joint), sum(dt$atac_pass), sum(dt$gex_pass)))
