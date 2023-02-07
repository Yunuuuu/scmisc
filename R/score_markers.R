#' Score marker genes
#' @description Compute various summary scores for potential marker genes to
#' distinguish between groups of cells.
#' @param x A matrix-like object containing log-normalized expression values,
#' with genes in rows and cells in columns. Alternatively, a
#' [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment] or
#' [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] object
#' containing such a matrix in its assays.
#' @param restricted An atomic vector coerced to characters containing the
#' subset of groups in clusters to be returned. By default, all unique groups in
#' clusters will be returned, but this can be restricted to specific groups of
#' interest.
#' @param top_n A scalar integer indicates derive how many top features from the
#' results. If cellmarker is `TRUE`, these features were then searched from
#' cellmarker database. If `NULL`, all features will be returned. Default:
#' `20L`.
#' @param order_by A character vector of column names by which to order. By
#' default, sorts over `mean.AUC`. See [scoreMarkers][scran::scoreMarkers] for
#' full summary scores names.
#' @param order An integer vector with only possible values of `1`and `-1`,
#' corresponding to ascending and descending order. The length of order must be
#' either 1 or equal to that of `order_by`. If `length(order) == 1`, it is
#' recycled to `length(order_by)`.
#' @param cellmarker A scalar logical indicates if search the CellMarker
#' database use the `top_n` genes. See [cellmarker_search()]
#' @param clusters A factor (or vector coercible into a factor) specifying the
#'   group to which each cell in `x` belongs. Alternatively, String specifying
#'   the field of `colData(x)` containing the grouping factor if `x` is a
#'   [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment], e.g.
#'   See [scoreMarkers][scran::scoreMarkers] groups argument.
#' @param lfc A numeric scalar specifying the log-fold change threshold to
#' compute effect sizes against.
#' @param features This can be a logical, integer or character vector indicating
#' the rows of x to use. See [scoreMarkers][scran::scoreMarkers] subset.row
#' argument
#' @param ... For the generic, further arguments to pass to specific methods.
#' @return A data.frame
#' @name score_markers
NULL

#' @keywords internal
score_markers_internal <- function(x, restricted = NULL, top_n = 20L, order_by = "mean.AUC", order = -1L, cellmarker = FALSE, clusters = NULL, lfc = 1L, features = NULL, ...) { # nolint
    score_markers_list <- cached_score_markers(
        x,
        clusters = clusters, lfc = lfc,
        features = features,
        ...
    )
    if (!is.null(restricted)) {
        score_markers_list <- score_markers_list[as.character(restricted)]
    }
    res <- lapply(score_markers_list, function(score_markers) {
        if (is.null(top_n)) top_n <- nrow(score_markers)
        data.table::setorderv(score_markers, order_by, order = order)
        score_markers <- score_markers[seq_len(top_n)]
        if (cellmarker) {
            cellmarker_search(score_markers$genes)
        } else {
            data.table::setDF(score_markers)[]
        }
    })
    if (length(res) == 1L) {
        return(res[[1L]])
    } else {
        res
    }
}

#' @export
#' @rdname score_markers
setGeneric(
    "score_markers",
    function(x, ...) standardGeneric("score_markers")
)

#' @export
#' @rdname score_markers
setMethod("score_markers", "ANY", score_markers_internal)

#' @export
#' @param assay.type A string or integer scalar indicating which `assays` in the
#' `x` contains the count matrix.
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @rdname score_markers
setMethod(
    "score_markers", "SingleCellExperiment",
    function(x, ..., clusters = NULL, features = NULL, assay.type = "logcounts") {
        if (is.null(clusters)) {
            clusters <- x$label
            if (is.null(clusters)) {
                cli::cli_abort("{.field label} never exist in {.arg x}.")
            }
        } else if (rlang::is_scalar_character(clusters) && ncol(x) > 1L) {
            clusters <- scater::retrieveCellInfo(
                x, clusters,
                search = "colData"
            )$value
        }
        if (is.null(features)) {
            features <- tryCatch(
                SingleCellExperiment::rowSubset(x, onAbsence = "error"),
                error = function(cnd) NULL
            )
        }
        score_markers_internal(x,
            clusters = clusters,
            features = features,
            ...,
            assay.type = assay.type
        )
    }
)

#' @export
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @rdname score_markers
setMethod(
    "score_markers", "SummarizedExperiment",
    function(x, ..., clusters = NULL, assay.type = "logcounts") {
        if (is.null(clusters)) {
            clusters <- x$label
            if (is.null(clusters)) {
                cli::cli_abort("{.field label} never exist in {.arg x}.")
            }
        }
        score_markers_internal(x,
            clusters = clusters,
            ...,
            assay.type = assay.type
        )
    }
)

#' @keywords internal
cached_score_markers <- function(x, clusters, lfc, features, ...) {
    score_markers_list <- scran::scoreMarkers(
        x,
        groups = clusters, lfc = lfc,
        subset.row = features,
        ...
    )
    lapply(score_markers_list, function(x) {
        score_markers <- as.data.frame(
            x,
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
        data.table::setDT(score_markers, keep.rownames = "genes")
        score_markers
    })
}

rlang::on_load(cached_score_markers <- memoise::memoize(
    cached_score_markers
))
