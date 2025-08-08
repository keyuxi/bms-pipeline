#!/usr/bin/env Rscript
# Minimal per-cell GEX processing + label transfer + merge
# Inputs (CLI):
#   --samples_csv   : CSV with columns sample_id,gex_dir (one row per sample)
#   --allen_ref     : path to Allen Seurat RDS (precomputed reference)
#   --out_rds       : merged Seurat object (RDS)
#   --out_qc        : per-cell QC TSV
#   --out_fig       : tiny PDF (UMAP by class & by sample)
# Optional thresholds:
#   --min_features 500 --max_features 7500 --min_counts 1000 --max_counts 80000 --max_percent_mt 5
#   --anchor_dims 40 --umap_dims 30 --cluster_res 5

suppressPackageStartupMessages({
  library(optparse); library(Seurat); library(ggplot2)
  library(scDblFinder); library(SingleCellExperiment)
})

opt_list <- list(
  make_option("--samples_csv",  type="character"),
  make_option("--aligned_dir",  type="character"),
  make_option("--allen_ref",    type="character"),
  make_option("--out_rds",      type="character"),
  make_option("--out_qc",       type="character"),
  make_option("--out_fig",      type="character"),
  make_option("--min_features", type="integer", default=500),
  make_option("--max_features", type="integer", default=7500),
  make_option("--min_counts",   type="integer", default=1000),
  make_option("--max_counts",   type="integer", default=80000),
  make_option("--max_percent_mt", type="double", default=5),
  make_option("--anchor_dims",  type="integer", default=40),
  make_option("--umap_dims",    type="integer", default=30),
  make_option("--cluster_res",  type="double",  default=5.0)
)
opt <- parse_args(OptionParser(option_list=opt_list))

dir.create(dirname(opt$out_rds), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(opt$out_qc),  recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(opt$out_fig), recursive=TRUE, showWarnings=FALSE)

# ---- util helper -------------------------------------------------
util_path <- "workflow/scripts/utils/load_seurat_object.R"
if (file.exists(util_path)) source(util_path)

# ---- inputs ------------------------------------------------------------------
samples <- read.csv(opt$samples_csv, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot("sample_id" %in% colnames(samples))

allen <- readRDS(opt$allen_ref)
label.cols <- c("cluster","subclass","neighborhood","class")

# ---- per-sample objects ------------------------------------------------------
mk_obj_from_dir <- function(dir_path, sample_id) {
  obj <- load_seurat_object(dir_path, project=sample_id)

  obj$Sample <- sample_id
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern="^mt-")
  obj <- subset(obj, subset =
                  nFeature_RNA > opt$min_features & nFeature_RNA < opt$max_features &
                  nCount_RNA   > opt$min_counts   & nCount_RNA   < opt$max_counts   &
                  percent.mt   < opt$max_percent_mt)
  obj <- NormalizeData(obj, verbose=FALSE)
  obj <- FindVariableFeatures(obj, selection.method="mvp", nfeatures=2000, verbose=FALSE)
  obj <- ScaleData(obj, verbose=FALSE)

  # Label transfer from Allen reference
  keptDims <- seq_len(opt$anchor_dims)
  anchors  <- FindTransferAnchors(reference=allen, query=obj, dims=keptDims, reference.reduction="pca")
  preds_df <- NULL
  for (lab in label.cols) {
    ref_col <- paste0(lab, "_label")
    stopifnot(ref_col %in% colnames(allen[[]]))
    p <- TransferData(anchorset=anchors, refdata=allen[[]][, ref_col, drop=TRUE])
    p <- p[, c(1, ncol(p))]; colnames(p) <- c(lab, paste0(lab, ".prediction.score.max"))
    preds_df <- if (is.null(preds_df)) as.data.frame(p) else cbind(preds_df, p)
  }
  preds_df <- preds_df[Cells(obj), , drop=FALSE]
  obj <- AddMetaData(obj, metadata=preds_df)

  # Doublets
  sce <- as.SingleCellExperiment(obj)
  sce <- scDblFinder(sce, clusters=TRUE, dbr.sd=0, dbr.per1k=0.006)
  obj$scDblFinder_class <- sce$scDblFinder.class
  obj$scDblFinder_score <- sce$scDblFinder.score

  obj
}

gex_list <- setNames(vector("list", nrow(samples)), samples$sample_id)
for (i in seq_len(nrow(samples))) {
  sid <- samples$sample_id[i]
  sdir <- file.path(opt$aligned_dir, sid)
  message("GEX: ", sid, " <- ", sdir)
  gex_list[[sid]] <- mk_obj_from_dir(sdir, sid)
}
stopifnot(length(gex_list) >= 1)

# ---- merge + downstream dims -------------------------------------------------
# Features chosen across objects
features <- SelectIntegrationFeatures(object.list = gex_list)
if (length(gex_list) == 1) {
  gex <- gex_list[[1]]
} else {
  gex <- merge(gex_list[[1]], gex_list[2:length(gex_list)],
               add.cell.ids=names(gex_list), project="PFC")
}
# Replace second "_" with "#" in cell names
replace_second_underscore <- function(s) gsub("^(.*?_.*?)_(.*)$", "\\1#\\2", s)
gex <- RenameCells(gex, new.names = sapply(Cells(gex), replace_second_underscore, USE.NAMES=FALSE))

# PCA/UMAP/NN/clusters (keep logic)
gex <- NormalizeData(gex, verbose=FALSE)
gex <- ScaleData(gex, verbose=FALSE)
gex <- RunPCA(gex, features=features, verbose=FALSE)

keptDims <- seq_len(opt$umap_dims)
gex <- RunUMAP(gex, dims=keptDims, seed.use=42, verbose=FALSE)
gex <- FindNeighbors(gex, dims=keptDims)
gex <- FindClusters(gex, resolution=opt$cluster_res)

# ---- outputs -----------------------------------------------------------------
saveRDS(gex, opt$out_rds)

md <- gex@meta.data
qc_cols <- intersect(c("Sample","nFeature_RNA","nCount_RNA","percent.mt",
                       "scDblFinder_class","scDblFinder_score",
                       "class","subclass","neighborhood","cluster","seurat_clusters"),
                     colnames(md))
qc <- data.frame(cell = rownames(md), md[, qc_cols, drop=FALSE], check.names=FALSE)
write.table(qc, opt$out_qc, sep="\t", quote=FALSE, row.names=FALSE)

pdf(opt$out_fig, width=8, height=4)
if ("class" %in% colnames(md))   print(DimPlot(gex, reduction="umap", group.by="class", label=TRUE) + ggtitle("GEX UMAP by class"))
if ("Sample" %in% colnames(md))  print(DimPlot(gex, reduction="umap", group.by="Sample") + ggtitle("GEX UMAP by sample"))
dev.off()
