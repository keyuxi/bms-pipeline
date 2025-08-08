#!/usr/bin/env Rscript
# Merge of dataset_annotation + cell_type_labels:
# - Read filtered ArchR project (from atac_per_cell_qc)
# - Read merged Seurat object (from gex_per_cell)
# - Read joint cell filter (keep_joint)
# - Filter both to BMS_clean
# - Copy Seurat labels (class/subclass/neighborhood/cluster/Seurat clusters/Sample) into ArchR colData
# - Save filtered ArchR project + Seurat + annotations TSV
#
# No parameter tuning. Hardcoded behavior. 4-space indents.

suppressPackageStartupMessages({
    library(optparse)
    library(ArchR)
    library(Seurat)
    library(dplyr); library(tidyr)
    library(data.table)
})

opt_list <- list(
    make_option("--archr_in",  type = "character"),   # e.g. results/atac/single_qc/ArchRProject
    make_option("--gex_in",    type = "character"),   # e.g. results/gex/merged/gex_merged.rds
    make_option("--filters",   type = "character"),   # e.g. results/joint/cell_filters.tsv
    make_option("--archr_out", type = "character"),   # e.g. results/joint/BMS_clean/ArchRProject
    make_option("--gex_out",   type = "character"),   # e.g. results/joint/BMS_clean/gex_BMS_clean.rds
    make_option("--anno_out",   type = "character")    # e.g. results/joint/BMS_clean/annotations.tsv
)
opt <- parse_args(OptionParser(option_list = opt_list))

dir.create(opt$archr_out, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$gex_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$ann_out), recursive = TRUE, showWarnings = FALSE)

# ---- load inputs -------------------------------------------------------------
addArchRGenome("mm10")
addArchRThreads(threads = max(1, parallel::detectCores() - 2))

proj <- loadArchRProject(opt$archr_in, showLogo = FALSE)
gex  <- readRDS(opt$gex_in)
filters <- fread(opt$filters, sep = "\t", header = TRUE)

stopifnot("cell" %in% names(filters), "keep_joint" %in% names(filters))
keep_cells <- filters$cell[filters$keep_joint == 1L]

# ---- filter both objects ------------------------------------------------------
proj <- proj[keep_cells, ]
gex  <- subset(gex, cells = keep_cells)

# ---- propagate Seurat annotations into ArchR ----------------------------------
md <- gex@meta.data
ann_cols <- intersect(c("Sample", "class", "subclass", "neighborhood", "cluster", "seurat_clusters"),
                      colnames(md))
ann <- md[keep_archr, ann_cols, drop = FALSE]   # rows by ArchR cells (same IDs)
for (nm in colnames(ann)) {
    proj <- addCellColData(
        ArchRProj = proj,
        data      = ann[[nm]],
        name      = nm,
        cells     = rownames(ann),
        force     = TRUE
    )
}

# ---- write annotations TSV ----------------------------------------------------
ann_out <- data.table(cell = rownames(ann), ann, check.names = FALSE)
fwrite(ann_out, opt$ann_out, sep = "\t", quote = FALSE)

# ---- save filtered objects ----------------------------------------------------
saveArchRProject(proj, outputDirectory = opt$archr_out, load = TRUE)
saveRDS(gex, opt$gex_out)

message(sprintf("BMS_clean: %d ArchR cells, %d GEX cells.",
                nCells(proj), ncol(gex)))
