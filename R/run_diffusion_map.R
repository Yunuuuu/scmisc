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
#' @param dpt A scalar logical value indicates whether to run
#'   [DPT][destiny::DPT].
#' @param ... Other arguments passed to [destiny][destiny::DiffusionMap]
#' @return A list or a modified x that contains the DiffusionMap coordinates in
#'   reducedDim(x, "DiffusionMap"), if dpt is TRUE, branch and dpt are also kept
#'   in `colData(x)`. 
#' @name run_diffusion_map
NULL

#' @export
#' @rdname run_diffusion_map
setGeneric(
    "run_diffusion_map",
    function(x, ...) standardGeneric("run_diffusion_map")
)

run_diffusion_map_internal <- function(x, dpt = TRUE, ...) {
    if (!rlang::is_scalar_logical(dpt)) {
        cli::cli_abort("{.arg dpt} must be a scalar logical value.")
    }
    dm_res <- destiny::DiffusionMap(x, ...)
    out <- list(evs = destiny::eigenvectors(dm_res))
    if (dpt) {
        dpt_res <- destiny::DPT(dm_res)

        # save branch data
        destiny_branch <- dpt_res@branch[, 1L, drop = TRUE]

        # save DPT data
        tip_cells <- destiny::tips(dpt_res) # ordered by the branch ids
        destiny_dpt <- S4Vectors::DataFrame(t(dpt_res[tip_cells]))
        S4Vectors::metadata(destiny_dpt) <- list(tip_cells = tip_cells)
        out <- c(out, list(branch = destiny_branch, dpt = destiny_dpt))
    }
    out
}

#' @export
#' @rdname run_diffusion_map
setMethod("run_diffusion_map", "ANY", run_diffusion_map_internal)

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
    function(x, ..., dpt = TRUE,
             dimred = "PCA", 
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
        res <- run_diffusion_map_internal(x = diffusion_input, dpt = dpt, ...)
        SingleCellExperiment::reducedDim(x, "DiffusionMap") <- res$evs
        if (dpt) {
            x$destiny_branch <- res$branch
            x$destiny_dpt <- res$dpt
        }
        x
    }
)
