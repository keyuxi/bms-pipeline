library(here)
here::i_am(".gitignore")
library(ArchR)
addArchRGenome("mm10")
addArchRThreads(threads = parallel::detectCores() - 2)
library(dbscan)
library(dplyr)
library(tidyr)
library(ggplot2)
library("BSgenome.Mmusculus.UCSC.mm10")

library(TFBSTools)
library(JASPAR2022)

# -------- Manual Args --------
archr.proj.path <- here::here("Data/ArchRProjects/BMS_diff")
archr.proj.path.peak <- here::here("Data/ArchRProjects/BMS_diff_peak")

  # -------- Setup --------
# Make sure the package directory is attached correctly
# pkg.dir <- "/oak/stanford/groups/wjg/kyx/software/rocker-chromatin/macs2"
# if (!(pkg.dir %in% .libPaths())) {
#     .libPaths( c( .libPaths(), pkg.dir) )
# }

proj <- loadArchRProject(archr.proj.path, showLogo = FALSE)

# # -------- Make Pseudo-bulk Replicates --------
# Approximately one pseudobulk per cell_type x Sample
# Unless too few cells
proj <- addGroupCoverages(ArchRProj = proj,
                          groupBy = "cell_type",
                          minCells = 20,
                          maxCells = 300,
                          sampleRatio = 1,
                          minReplicates = 13,
                          maxReplicates = 16,
                          returnGroups = FALSE,
                          force = TRUE)
saveArchRProject(proj, outputDirectory = archr.proj.path.peak, load = TRUE)

# # -------- Call Peaks -------
# pathToMacs2 <- findMacs2()
pathToMacs2 <- "/Users/kyx/software/miniconda/envs/bms/bin/macs2"
# Reuses pseudobulks generated with the same groupBy argument
# Lower the reproducibility argument to catch de novo peaks
# Can be fixed with downstream filtering if signal is too weak
proj <- addReproduciblePeakSet(
    ArchRProj = proj,
    groupBy = "cell_type",
    reproducibility = "1",
    pathToMacs2 = pathToMacs2
)

peakSet <- getPeakSet(proj)
peakSet_noX <- peakSet[seqnames(peakSet) != "chrX"]
peakSet_noX <- dropSeqlevels(peakSet_noX, c("chrX"), pruning.mode = "coarse")
peakSet_noX$peak_name <- paste(seqnames(peakSet_noX), start(peakSet_noX), sep = "_")
names(peakSet_noX) <- peakSet_noX$peak_name

proj <- addPeakSet(
    ArchRProj = proj,
    peakSet   = peakSet_noX,
    force     = TRUE
)

arrow_dir <- proj@sampleColData@listData[["ArrowFiles"]]
proj@sampleColData@listData[["ArrowFiles"]] <- sapply(names(arrow_dir),
    function(x) paste0(here::here("Data/ArchRProjects/BMS_diff_peak/ArrowFiles/"), x, ".arrow"))

proj <- addPeakMatrix(proj, force = TRUE)
saveRDS(proj, here::here("Data/ArchRProjects/BMS_diff_peak/Save-ArchR-Project.rds"))
saveArchRProject(proj, outputDirectory = archr.proj.path.peak, load = FALSE)

# # -------- ChromVAR --------
# force <- TRUE
# if (force){
#     proj <- addMotifAnnotations(
#         ArchRProj   = proj,
#         motifSet    = "JASPAR2020",     # newest collection
#         annoName    = "JASPAR",
#         force       = TRUE
#     )
#     proj <- addBgdPeaks(proj, force = T)
# }
#
# proj <- addDeviationsMatrix(
#   ArchRProj = proj,
#   peakAnnotation = "JASPAR",
#   force = TRUE
# )
# plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
# plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
#
#
# saveArchRProject(proj, outputDirectory = new.archr.proj.path, load = FALSE)
