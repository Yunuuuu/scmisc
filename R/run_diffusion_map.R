#' Create a diffusion map of cells
#'
#' @description
#' A modified version of [DiffusionMap][destiny::DiffusionMap], this function
#' will add the DiffusionMap results components into x if it is a
#' [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] object.
#' @param x A numeric matrix of counts with cells in columns and features in
#' rows.
#'
#' Alternatively, a
#' [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] object
#' containing such a matrix.
#' @param root The cell index/indices from which to calculate the DPT(s)
#'   (integer of length 1-3). If length more than `size`, we will choose the top
#'   size elements with the maximal eigenvectors CD1 values. Can be also a
#'   logical atomic value, in this way, we will sample the elements with TRUE
#'   values.
#' @param w_width Window width to use for deciding the branch cutoff.
#' @param dpt A scalar logical value indicates whether to run
#'   [DPT][destiny::DPT].
#' @param size A scalar numeric ranges from 1 to 3 specifying how many roots
#'   should be used.
#' @param ... Other arguments passed to [DiffusionMap][destiny::DiffusionMap].
#' @return A [DiffusionMap][destiny::DiffusionMap] or [DPT][destiny::DPT]
#'   object depend on `dpt`.
#' @name run_diffusion_map
NULL

#' @export
#' @rdname run_diffusion_map
setGeneric(
    "run_diffusion_map",
    function(x, ...) standardGeneric("run_diffusion_map")
)

run_diffusion_map_internal <- function(x, root = NULL, ..., w_width = 0.1, dpt = TRUE, size = 1L) {
    if (!rlang::is_scalar_logical(dpt)) {
        cli::cli_abort("{.arg dpt} must be a scalar logical value.")
    }
    if (size <= 0 || size > 3L) {
        cli::cli_abort("{.arg size} must range from 1L to 3L.")
    }
    out <- destiny::DiffusionMap(x, ...)
    if (dpt) {
        if (is.null(root)) {
            out <- destiny::DPT(out, w_width = w_width)
        } else {
            root <- handle_root(out, root, size = size)
            out <- destiny::DPT(out, tips = root, w_width = w_width)
        }
    }
    out
}

#' @export
#' @rdname run_diffusion_map
setMethod("run_diffusion_map", "ANY", run_diffusion_map_internal)

#' @export
#' @rdname run_diffusion_map
setMethod("run_diffusion_map", "DiffusionMap", function(x, root = NULL, ..., w_width = 0.1, dpt = TRUE, size = 1L) {
    if (!rlang::is_scalar_logical(dpt)) {
        cli::cli_abort("{.arg dpt} must be a scalar logical value.")
    }
    if (size <= 0 || size > 3L) {
        cli::cli_abort("{.arg size} must range from 1L to 3L.")
    }
    out <- x
    if (dpt) {
        if (is.null(root)) {
            out <- destiny::DPT(out, w_width = w_width)
        } else {
            root <- handle_root(out, root, size = size)
            out <- destiny::DPT(out, tips = root, w_width = w_width)
        }
    }
    out
})

#' @param dimred String or integer scalar specifying the existing dimensionality
#'   reduction results to use.
#' @param n_dimred Integer scalar or vector specifying the dimensions to use if
#'   dimred is specified.
#' @param assay.type A string or integer scalar indicating which `assays` in the
#'   `x` contains the count matrix. Only used when dimred is `NULL`.
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @export
#' @rdname run_diffusion_map
setMethod(
    "run_diffusion_map", "SingleCellExperiment",
    function(x, dimred = "PCA", ...,
             n_dimred = NULL,
             assay.type = "logcounts") {
        if (is.null(dimred) && is.null(assay.type)) {
            cli::cli_abort("One of {.arg dimred} or {.arg assay.type} must be specified.")
        }
        if (!is.null(dimred)) {
            diffusion_input <- SingleCellExperiment::reducedDim(x, dimred)
            if (!is.null(n_dimred)) {
                if (is_scalar_numeric(n_dimred)) {
                    n_dimred <- min(n_dimred, ncol(diffusion_input),
                        na.rm = TRUE
                    )
                    n_dimred <- seq_len(n_dimred)
                }
                diffusion_input <- diffusion_input[, n_dimred, drop = FALSE]
            }
        } else {
            diffusion_input <- SummarizedExperiment::assay(
                diffusion_input, assay.type
            )
        }
        run_diffusion_map_internal(x = diffusion_input, ...)
    }
)

handle_root <- function(dm, root, size) {
    if (is.numeric(root)) {
        if (length(root) <= size) {
            return(as.integer(root))
        } else {
            idx <- root
        }
    } else if (is.logical(root)) {
        idx <- which(root, useNames = FALSE)
    } else {
        cli::cli_abort("Unsupported type of {.arg root}")
    }
    evs <- destiny::eigenvectors(dm)[, 1L]
    seq_along(evs)[order(evs, decreasing = TRUE)][idx][
        seq_len(size)
    ]
}
