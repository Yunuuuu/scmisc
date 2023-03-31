#' Reads object in HDF5 files
#' @inheritParams rhdf5::h5read
#' @inheritDotParams rhdf5::h5read -file -name
#' @return A sparse matrix
#' @export 
read_h5 <- function(file, name = "matrix", ...) {
    h5 <- rhdf5::h5read(file, name = name)
    mat <- Matrix::sparseMatrix(
        dims = h5$shape,
        i = as.numeric(h5$indices),
        p = as.numeric(h5$indptr),
        x = as.numeric(h5$data),
        index1 = FALSE
    )
    colnames(mat) <- h5$barcodes
    rownames(mat) <- scuttle::uniquifyFeatureNames(
        h5$features$id, h5$features$name
    )
    mat
}
