.calculate_diffusion_map <- function(x, ncomponents = 20L, ..., tips = NULL, w_width = 0.1) {
    assert_pkg("destiny")
    dm <- destiny::DiffusionMap(
        data = x, ...,
        n_eigs = ncomponents, rotate = FALSE,
        n_pcs = NA, suppress_dpt = FALSE
    )
    if (is.null(tips)) {
        dpt <- destiny::DPT(dm, w_width = w_width)
    } else {
        dpt <- destiny::DPT(dm, tips = tips, w_width = w_width)
    }
    evs <- destiny::eigenvectors(dm)
    attr(dpt, "w_width") <- w_width
    attr(evs, "DPT") <- dpt
    evs
}

#' Create a diffusion map of cells
#'
#' @param x For `calculateDiffusionMap`, a numeric matrix of log-expression
#' values where rows are features and columns are cells. Alternatively, a
#' [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment] or
#' [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] object
#' containing such a matrix.
#'
#' For `runDiffusionMap`, a
#' [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] object.
#' @inheritParams scater::runPCA
#' @inheritParams destiny::DPT
#' @param ...
#' For the `calculateDiffusionMap` generic, additional arguments to
#' pass to specific methods.
#'
#' For the [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment] or
#' [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] methods,
#' additional arguments to pass to the ANY method.
#'
#' For the ANY methods, additional arguments to pass to [destiny::DiffusionMap]
#'  - `sigma`: Diffusion scale parameter of the Gaussian kernel. One of
#'             \code{'local'}, \code{'global'}, a (\link[base]{numeric}) global
#'             sigma or a \link{Sigmas} object.  When choosing \code{'global'},
#'             a global sigma will be calculated using
#'             \code{\link{find_sigmas}}.  (Optional. default: \code{'local'}) A
#'             larger sigma might be necessary if the eigenvalues can not be
#'             found because of a singularity in the matrix.
#'  - `k`: Number of nearest neighbors to consider (default: a guess betweeen
#'         100 and \eqn{n - 1}. See \code{\link{find_dm_k}}).
#'  - `density_norm`: logical. If TRUE, use density normalisation.
#'  - `n_local`: If \code{sigma == 'local'}, the \code{n_local}th nearest
#'               neighbor(s) determine(s) the local sigma.
#'  - `censor_val`: Value regarded as uncertain. Either a single value or one
#'                  for every dimension (Optional, default: censor_val).
#'  - `censor_range`: Uncertainity range for censoring (Optional, default:
#'                   none).  A length-2-vector of certainty range start and end.
#'                   TODO: also allow \eqn{ 2 \times G } matrix.
#'  - `missing_range`: Whole data range for missing value model. Has to be
#'                    specified if NAs are in the data.
#'  - `knn_params`: Parameters passed to \code{\link{find_knn}}.
#'  - `verbose`: Show a progressbar and other progress information (default: do
#'              it if censoring is enabled).
#'
#' For `runDiffusionMap`, additional arguments to pass to
#' `calculateDiffusionMap`.
#'
#' @return
#' For `calculateDiffusionMap`, a numeric matrix of coordinates for each cell
#' (row) in each of \code{ncomponents} DiffusionMap (column).
#'
#' For `runDiffusionMap`, a
#' [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] object is
#' returned containing this matrix in \code{\link{reducedDims}(..., name)}.
#'
#' In both cases, the attributes of the DiffusionMap coordinate matrix contain
#' the `DPT`: the raw output of [DPT][destiny::DPT].
#' @export
#' @name calculateDiffusionMap
methods::setGeneric("calculateDiffusionMap", function(x, ...) {
    standardGeneric("calculateDiffusionMap")
})

#' @export
#' @rdname calculateDiffusionMap
methods::setMethod("calculateDiffusionMap", "ANY", .calculate_diffusion_map)

#' @export
#' @rdname calculateDiffusionMap
methods::setMethod(
    "calculateDiffusionMap", "SummarizedExperiment",
    function(x, ..., assay.type = "logcounts") {
        .calculate_diffusion_map(
            SummarizedExperiment::assay(x, assay.type), ...
        )
    }
)

#' @export
#' @rdname calculateDiffusionMap
methods::setMethod(
    "calculateDiffusionMap", "SingleCellExperiment",
    function(x, ..., assay.type = "logcounts", dimred = "PCA", n_dimred = NULL) {
        mat <- .get_mat_from_sce(x,
            assay.type = assay.type,
            dimred = dimred, n_dimred = n_dimred
        )
        .calculate_diffusion_map(mat, ...)
    }
)

#' @export
#' @rdname calculateDiffusionMap
runDiffusionMap <- function(x, ..., altexp = NULL, name = "DiffusionMap") {
    assert_s4_class(x, "SingleCellExperiment")
    if (!is.null(altexp)) {
        y <- SingleCellExperiment::altExp(x, altexp)
    } else {
        y <- x
    }
    SingleCellExperiment::reducedDim(x, name) <- calculateDiffusionMap(y, ...)
    x
}

.get_mat_from_sce <- function(x, assay.type, dimred, n_dimred) {
    if (!is.null(dimred)) {
        mat <- SingleCellExperiment::reducedDim(x, dimred)
        if (!is.null(n_dimred)) {
            if (length(n_dimred) == 1L) {
                n_dimred <- seq_len(n_dimred)
            }
            mat <- mat[, n_dimred, drop = FALSE]
        }
        mat
    } else {
        SummarizedExperiment::assay(x, assay.type)
    }
}
