#' Returns a Seurat object with just the RNA counts
#'
#' @param matrix.dir String, directory with `barcodes.tsv.gz`, `features.tsv.gz`,
#' and `matrix.mtx.gz`
#' @param project String, project name of the Seurat object
#'
#' @return A Seurat object

load_seurat_object <- function(matrix.dir, project = "bms"){
    # Returns a Seurat object with just the RNA counts
    barcode.path <- file.path(matrix.dir, "barcodes.tsv.gz")
    features.path <- file.path(matrix.dir, "features.tsv.gz")
    matrix.path <- file.path(matrix.dir, "matrix.mtx.gz")
    if (!file.exists(matrix.path)){
        matrix.path <- file.path(matrix.dir, "matrix.mtx")
    }
    mat <- Matrix::readMM(file = matrix.path)
    feature.names = read.delim(features.path,
                               header = FALSE,
                               stringsAsFactors = FALSE)
    barcode.names = read.delim(barcode.path,
                               header = FALSE,
                               stringsAsFactors = FALSE)

    # Remove chr M, X, Y
    gex.mask <- (feature.names$V3 == "Gene Expression") & !(feature.names$V4 %in% c("chrM", "chrX", "chrY"))
    feature.names <- feature.names[gex.mask,]
    mat <- mat[gex.mask,]

    # Collapse duplicated gene names
    # Identify the first occurrence of each gene
    dedup_indices <- !duplicated(feature.names$V2)

    # Subset the matrix and feature names to keep only the first occurrence
    mat <- mat[dedup_indices, ]
    feature.names <- feature.names[dedup_indices, ]

    # Assign gene names to the new matrix
    colnames(mat) <- barcode.names$V1
    rownames(mat) <- feature.names$V2
    mat <- as(mat, "dgCMatrix")

    # bms <- CreateSeuratObject(counts = mat, project = project, names.delim = "#")
    #' `names.delim` doesn't work as what you expected it to do!

    gex.assay <- CreateAssayObject(counts = mat, assay = "RNA")
    bms <- CreateSeuratObject(gex.assay, project = project)
    bms$RNA[["Ensembl_id"]] <- feature.names$V1
    bms$RNA[["chr"]] <- feature.names$V4
    return(bms)
}
