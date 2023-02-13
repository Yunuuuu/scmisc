#' Create a diffusion map of cells
#'
#' @description
#' A modified version of [DiffusionMap][destiny::DiffusionMap] and
#' [DPT][destiny::DPT].
#' @param x A numeric matrix of counts with cells in columns and features in
#'   rows.
#'
#'   Alternatively, a
#'   [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] object
#'   containing such a matrix, or a [DiffusionMap][destiny::DiffusionMap]
#'   object.
#' @param root The cell index/indices from which to calculate the DPT(s)
#'   (integer of length `1-3`). If the length of root is more than `size` or  a
#'   logical atomic vector with the compatible length of `x`(the same length of
#'   `ncol(x)` if x is a simple matrix or
#'   [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] or the
#'   same length of `x@@d` if x is a [DiffusionMap][destiny::DiffusionMap]
#'   object), we will index with root value and choose the top `size` elements
#'   with the maximal [eigenvectors][destiny::eigenvectors] CD1 values.
#' @param w_width Window width to use for deciding the branch cutoff.
#' @param dpt A scalar logical value indicates whether to run
#'   [DPT][destiny::DPT].
#' @param size A scalar integer ranges from 1 to 3 specifying how many roots
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
        out <- run_dpt_helper(out, root = root, w_width = w_width, size = size)
    }
    out
}

run_dpt_helper <- function(x, root, w_width, size) {
    if (is.null(root)) {
        destiny::DPT(x, w_width = w_width)
    } else {
        root <- handle_root(x, root, size = size)
        destiny::DPT(x, tips = root, w_width = w_width)
    }
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
    if (dpt) {
        run_dpt_helper(x, root = root, w_width = w_width, size = size)
    } else {
        x
    }
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
        root <- as.integer(root)
        if (length(root) <= size) {
            return(root)
        }
        idx <- root
    } else if (is.logical(root)) {
        if (length(root) == length(dm@d)) {
            idx <- which(root, useNames = FALSE)
        } else {
            cli::cli_abort("the length of logical {.arg root} must be compatible with {.arg x}.")
        }
    } else {
        cli::cli_abort("Unsupported type of {.arg root}")
    }
    # order the index of `evs` firstly
    evs <- destiny::eigenvectors(dm)[, 1L, drop = TRUE]
    res <- seq_along(evs)[order(evs, decreasing = TRUE, na.last = TRUE)]

    # choose the index specified in root and get the top size elments
    res[res %in% idx][seq_len(size)]
}
