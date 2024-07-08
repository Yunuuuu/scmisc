#' Score marker genes
#' @description Compute various summary scores for potential marker genes to
#' distinguish between groups of cells.
#' @param x A matrix-like object containing log-normalized expression values,
#' with genes in rows and cells in columns. Alternatively, a
#' [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment] or
#' [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] object
#' containing such a matrix in its assays.
#' @param restricted An atomic vector coerced to characters containing the
#' subset of groups in `groups` to be returned. By default, all unique groups in
#' `groups` will be returned, but this can be restricted to specific groups of
#' interest.
#' @param top_n A scalar integer indicates derive how many top features from the
#' results. If cellmarker is `TRUE`, these features were then searched from
#' cellmarker database. If `NULL`, all features will be returned. Default:
#' `10L`.
#' @param order_by A character vector of column names by which to order. By
#' default, sorts over `mean.AUC`. See [scoreMarkers][scran::scoreMarkers] for
#' full summary scores names.
#' @param order An integer vector with only possible values of `1`and `-1`,
#' corresponding to ascending and descending order. The length of order must be
#' either 1 or equal to that of `order_by`. If `length(order) == 1`, it is
#' recycled to `length(order_by)`.
#' @param cellmarker A scalar logical indicates if search the CellMarker
#' database use the `top_n` genes. See [cellmarker_search()]
#' @param groups A factor (or vector coercible into a factor) specifying the
#'   group to which each cell in `x` belongs. Alternatively, if `x` is a
#'   [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment], e.g,
#'   String specifying the field of `colData(x)` containing the grouping factor.
#'   In this way, if `groups` is `NULL`, "label" in `colData(x)` will be
#'   extracted.
#' @param features This can be a logical, integer or character vector indicating
#' the rows of x to use. If `NULL`, and x is a
#' [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment],
#' [rowSubset][SingleCellExperiment::rowSubset] is used to derive the features.
#' See [scoreMarkers][scran::scoreMarkers] subset.row argument.
#' @param absolute_lfc A boolean value indicates whether to calculate absolute
#' logFC when filtering by `lfc_cutoff`. Only used when `lfc_cutoff` is not
#' `NULL`.
#' @param lfc_cutoff Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells.
#' @inheritDotParams scran::scoreMarkers
#' @return A data.frame
#' @name score_markers
NULL

#' @keywords internal
.score_markers <- function(x, groups = NULL,
                           restricted = NULL, top_n = 10L,
                           order_by = "mean.AUC", order = -1L,
                           ...,
                           features = NULL,
                           cellmarker = FALSE,
                           absolute_lfc = FALSE,
                           lfc_cutoff = NULL) { # nolint
    score_markers_list <- cached_score_markers(
        x,
        groups = groups,
        features = handle_row_data(object = x, features),
        ...
    )
    if (!is.null(restricted)) {
        score_markers_list <- score_markers_list[as.character(restricted)]
    }
    res <- lapply(score_markers_list, function(markers) {
        markers <- markers[, logFC := self.average - other.average]
        data.table::setcolorder(markers, "logFC", after = "other.average")
        if (!is.null(lfc_cutoff)) {
            if (absolute_lfc) {
                markers <- markers[abs(logFC) > lfc_cutoff]
            } else {
                markers <- markers[logFC > lfc_cutoff]
            }
        }
        if (is.null(top_n)) {
            top_n <- nrow(markers)
        } else {
            top_n <- min(top_n, nrow(markers))
        }
        data.table::setorderv(markers, order_by, order = order)
        markers <- markers[seq_len(top_n)]
        if (cellmarker) {
            cellmarker_search(markers$genes)
        } else {
            data.table::setDF(markers)
        }
    })
    if (length(res) == 1L) {
        return(res[[1L]])
    } else {
        res
    }
}

utils::globalVariables(c("self.average", "other.average", "logFC"))

#' @export
#' @rdname score_markers
setGeneric(
    "score_markers",
    function(x, ...) standardGeneric("score_markers")
)

#' @export
#' @rdname score_markers
setMethod("score_markers", "ANY", .score_markers)

#' @export
#' @param assay.type A string or integer scalar indicating which `assays` in the
#' `x` contains the count matrix.
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @rdname score_markers
setMethod(
    "score_markers", "SingleCellExperiment",
    function(x, ..., groups = NULL, features = NULL, assay.type = "logcounts") {
        groups <- handle_column_data(object = x, groups)
        if (is.null(features)) {
            features <- SingleCellExperiment::rowSubset(
                x,
                onAbsence = "none"
            )
        }
        .score_markers(x,
            groups = groups,
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
    function(x, ..., groups = NULL, assay.type = "logcounts") {
        .score_markers(x,
            groups = handle_column_data(object = x, groups),
            ...,
            assay.type = assay.type
        )
    }
)

#' @keywords internal
cached_score_markers <- function(x, groups, features, ...) {
    score_markers_list <- scran::scoreMarkers(
        x = x,
        groups = groups,
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
