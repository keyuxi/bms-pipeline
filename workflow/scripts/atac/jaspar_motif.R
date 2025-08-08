library(here)
here::i_am(".gitignore")
library(ArchR)
addArchRGenome("mm10")
addArchRThreads(threads = parallel::detectCores() - 2)
library(TFBSTools)
library(JASPAR2022)
library("BSgenome.Mmusculus.UCSC.mm10")

# --------------------------------------------------------------------
archr.proj.path.peak <- here::here("Data/ArchRProjects/BMS_diff_peak")
archr.proj.path.motif <- here::here("Data/ArchRProjects/BMS_diff_motif")
outdir <- here::here("Output/Motif")
# --------------------------------------------------------------------
proj <- readRDS(here::here("Data/ArchRProjects/BMS_diff_peak/Save-ArchR-Project.rds"))
# proj.peak <- loadArchRProject(archr.proj.path.peak, showLogo = FALSE)

pfm_list <- TFBSTools::getMatrixSet(
    JASPAR2022,
    opts = list(collection = "CORE", tax_group = "vertebrates")
)

# Build a quick lookup table between IDs and TF names
ids       <- sapply(pfm_list, ID)              # e.g. "MA1123.1"
tf_names  <- sapply(pfm_list, name)            # e.g. "TWIST1"
motif_tbl <- data.frame(id = ids, tf = tf_names, stringsAsFactors = FALSE)

# Helper: check if the TFs of interest are included
# motif_tbl[ grepl("neurod|twist|tcf3|Mef2c", motif_tbl$tf, ignore.case = TRUE), ]

# Convert to PWMatrixList that ArchR takes
pwm_vec  <- lapply(pfm_list, toPWM)
pwm_list <- do.call(PWMatrixList, pwm_vec)

# Add these PWMs to ArchR as a custom motif collection
proj <- addMotifAnnotations(
    ArchRProj  = proj,
    motifPWMs  = pwm_list,
    name       = "Motif",
    force = TRUE
)

proj <- addBgdPeaks(proj, force = T)

# Build the chromVAR deviation matrix
proj <- addDeviationsMatrix(
    ArchRProj       = proj,
    peakAnnotation  = "Motif",  # must match the name you just gave
    force           = TRUE
)

saveArchRProject(proj, outputDirectory = archr.proj.path.motif, load = FALSE)
saveRDS(proj, here::here("Data/ArchRProjects/BMS_diff_peak/Save-ArchR-Project.rds"))
saveRDS(pwm_list, file.path(outdir, "pwm_list.rds"))
saveRDS(motif_tbl, file.path(outdir, "motif_tbl.rds"))
