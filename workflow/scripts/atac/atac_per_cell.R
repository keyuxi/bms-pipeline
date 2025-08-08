#!/usr/bin/env Rscript

# Performs per-cell ATAC QC + filtering (and optional doublet scores)
# ATAC project and cell metadata sheets already filtered.
# Inputs:
#   --arrow_dir       : directory containing {sample}.arrow files
#   --samples         : TSV with a 'sample_id' column (optionally 'arrow_path')
#   --out_proj        : output ArchR project directory (will be created/updated)
#   --out_qc          : TSV of per-cell QC with pass/fail
#   --out_fig         : PDF with QC plots
#   --genome          : default mm10
#   --threads         : default: all-2
#   --min_tss         : default 5
#   --min_log10_nfrags: default 3.5  (i.e. nFrags >= 10^3.5)
#   --max_blacklist   : default 0.05
#   --copy_arrows     : TRUE/FALSE (default FALSE; avoid duplication)
#
# Output:
#   - ArchR project (filtered) at --out_proj
#   - QC table with columns: cell, Sample, TSSEnrichment, nFrags, BlacklistRatio, pass
#   - QC PDF at --out_fig

suppressPackageStartupMessages({
    library(optparse)
    library(ArchR)
    library(ggplot2)
    library(viridis)
})

opt_list <- list(
    make_option("--arrow_dir", type="character"),
    make_option("--samples", type="character"),
    make_option("--out_proj", type="character"),
    make_option("--out_qc", type="character"),
    make_option("--out_fig", type="character"),
    make_option("--genome", type="character", default="mm10"),
    make_option("--threads", type="integer", default=max(1, parallel::detectCores() - 2)),
    make_option("--copy_arrows", action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=opt_list))

dir.create(dirname(opt$out_qc), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(opt$out_fig), showWarnings = FALSE, recursive = TRUE)
dir.create(opt$out_proj, showWarnings = FALSE, recursive = TRUE)

addArchRGenome(opt$genome)
addArchRThreads(threads = opt$threads)

# ---- samples & arrow files ---------------------------------------------------
samps <- read.delim(opt$samples, sep=",", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
if (!"sample_id" %in% names(samps)) {
    stop("samples.csv must contain a 'sample_id' column")
}

if ("arrow_path" %in% names(samps)) {
    arrow_files <- samps$arrow_path
    names(arrow_files) <- samps$sample_id
} else {
    arrow_files <- file.path(opt$arrow_dir, paste0(samps$sample_id, ".arrow"))
    names(arrow_files) <- samps$sample_id
}
missing <- arrow_files[!file.exists(arrow_files)]
if (length(missing) > 0) {
    stop(sprintf("Missing Arrow files:\n%s", paste(missing, collapse = "\n")))
}

# ---- Load or create project --------------------------------------------------
proj <- if (file.exists(file.path(opt$out_proj, "ProjectMetadata.rds")) ||
            file.exists(file.path(opt$out_proj, "Metadata.rds"))) {
    loadArchRProject(opt$out_proj, showLogo = FALSE)
} else {
  ArchRProject(
    ArrowFiles       = unname(arrow_files),
    outputDirectory  = opt$out_proj,
    copyArrows       = isTRUE(opt$copy_arrows),
    showLogo         = FALSE
  )
}

# Filter cells & save
idx_pass <- which(
       ((proj$Sample %in% c("YK01_En7Apos", "YK04_En7Cneg", "YK05_En7Cpos")) |
        (proj$Sample == "YK02_En7Bneg" & proj$nFrags >= 10^(3.7)) |
        (proj$nFrags >= 10^(3.5))) &
           (proj$TSSEnrichment >= 5) & (proj$BlacklistRatio < 0.05)
)
cells_pass <- proj$cellNames[idx_pass]
proj <- proj[cells_pass, ]
saveArchRProject(proj)

# ---- QC plots ----------------------------------------------------------------
df <- data.frame(getCellColData(proj))
p1 <- plotGroups(ArchRProj = proj, groupBy = "Sample",
                 colorBy = "cellColData", name = "TSSEnrichment",
                 plotAs = "ridges")
p2 <- plotGroups(ArchRProj = proj, groupBy = "Sample",
                 colorBy = "cellColData", name = "nFrags",
                 plotAs = "ridges")
p3 <- ggplot(df, aes(x = log(TSSEnrichment), y = log(nFrags + 1))) +
  facet_wrap(vars(Sample)) +
  geom_density2d_filled(contour_var = "ndensity") +
  labs(x="log(TSS Enrichment)", y="log(nFrags + 1)")

pdf(opt$out_fig, width = 10, height = 7)
try(print(p1), silent = TRUE)
try(print(p2), silent = TRUE)
try(print(p3), silent = TRUE)
dev.off()

# ---- QC table ---------------------------------------------------------
qc_tab <- data.frame(
    cell = rownames(df),
    Sample = df$Sample,
    TSSEnrichment = df$TSSEnrichment,
    nFrags = df$nFrags,
    BlacklistRatio = df$BlacklistRatio,
    pass = as.integer(pass_vec),
    stringsAsFactors = FALSE
)
write.table(qc_tab, opt$out_qc, sep="\t", quote=FALSE, row.names=FALSE)