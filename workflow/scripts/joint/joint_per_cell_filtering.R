#!/usr/bin/env Rscript
# Minimal joint per-cell filtering:
#   - Read ATAC QC TSV (must have `cell`)
#   - Read GEX QC  TSV (must have `cell`; may have doublet column)
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
    library(dplyr)
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

# filter gex by seurat clusters
# p.s. originally in GEX_per_cell.Rmd, moved here
non_neuronal_clusters <- gex %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(
        non_neuronal_prop = mean(class == "Non-Neuronal"),
        .groups = "drop"
    ) %>%
    dplyr::filter(non_neuronal_prop > 0.6) %>%
    dplyr::pull(seurat_clusters)

`%ni%` <- Negate(`%in%`)
gex$is_neuron <- (gex$seurat_clusters %ni% non_neuronal_clusters) &
    gex$class.prediction.score.max > 0.4

# --- extra strict neuron/doublet/TSS filter (build per-cell join) ---
merge_dt <- merge(
    gex[, .(cell, scDblFinder_class, is_neuron)],
    atac[, .(cell, TSSEnrichment)],
    by = "cell", all = FALSE
)
kept_extra <- merge_dt[
    scDblFinder_class == "singlet" & as.logical(is_neuron) & TSSEnrichment >= 5.25,
    cell
]

keep_cells <- intersect(atac$cells, gex$cells)

# refine keep set with strict neuron/doublet/TSS condition
keep_cells <- intersect(keep_cells, kept_extra)

# annotate and write
dt <- data.table(
    cell = unique(c(atac$cell, gex$cell))
)
dt[, in_atac := as_int01(cell %in% atac$cell)]
dt[, in_gex := as_int01(cell %in% gex$cell)]
dt[, keep_joint := as_int01(cell %in% keep_cells)]

# optional: carry Sample if available
if ("Sample" %in% names(atac)) dt[match(atac$cell, dt$cell), Sample_atac := atac$Sample]
if ("Sample" %in% names(gex))  dt[match(gex$cell,  dt$cell), Sample_gex  := gex$Sample]

setcolorder(dt, c("cell", "keep_joint", "in_atac", "in_gex",
                  setdiff(names(dt), c("cell", "keep_joint", "in_atac", "in_gex"))))
fwrite(dt, opt$out_tsv, sep = "\t", quote = FALSE, na = "NA")

message(sprintf("Joint filter: kept %d cells (ATAC pass %d ∩ GEX pass %d).",
                sum(dt$keep_joint), sum(dt$in_atac), sum(dt$in_gex)))
