#' Reads object in HDF5 files
#' @param file The filename (character) of the file in which the dataset is be
#'  located. It is possible to provide an object of class
#'  [H5IdComponent][rhdf5::H5IdComponent-class] representing a H5 location
#'  identifier (file or group). See [H5Fcreate][rhdf5::H5Fcreate],
#'  [H5Fopen][rhdf5::H5Fopen], [H5Gcreate][rhdf5::H5Gcreate],
#'  [H5Gopen][rhdf5::H5Gopen] to create an object of this kind. 
#' @param name The name of the dataset in the HDF5 file. 
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
