#' Plot dots heatmap of group-level expression averages
#'
#' @description
#' `plot_grouped_dots` create a dot plot of expression values for a grouping of
#' cells, where the size and color of each dot represents the proportion of
#' detected expression values and the average expression, respectively, for each
#' feature in each group of cells.
#'
#' The order of the row or column will be the same with `unlist(marker_list)` or
#' `levels(groups)` depend on whether flip is `TRUE` or `FALSE` if
#' ComplexHeatmap clustering is turned off (if groups isn't a factor, the
#' internal will coerce it as a factor). So we can easily add annotation in
#' Heatmap as long as we provide a
#' [HeatmapAnnotation][ComplexHeatmap::HeatmapAnnotation] in the same order.
#' @param threshold Numeric value specifying the cap on the proportion of
#'   detected expression values.
#' @param ... Other arguments passed to [DotsHeatmap] and specific methods.
#' @inheritParams DotsHeatmap
#' @return A [DotsHeatmap] Object
#' @export
#' @name grouped_dots
#' @aliases plot_grouped_dots
setGeneric(
    "plot_grouped_dots",
    function(x, ...) standardGeneric("plot_grouped_dots")
)

plot_grouped_dots_internal <- function(
    x, marker_list = NULL, cluster2cell = NULL, groups = NULL, ...,
    blocks = NULL, colour = color, color = NULL, center = FALSE,
    scale = FALSE, threshold = 0L, zlim = NULL,
    flip = FALSE, row_labels = NULL, column_labels = NULL) {
    grouped_heat_internal(
        x = x, marker_list = marker_list,
        groups = groups, blocks = blocks,
        cluster2cell = cluster2cell,
        center = center, scale = scale,
        zlim = zlim, threshold = threshold,
        colour = colour, flip = flip,
        row_labels = row_labels, column_labels = column_labels,
        ..., graph_type = "dots"
    )
}

#' @inheritParams grouped_heatmap
#' @export
#' @rdname grouped_dots
setMethod("plot_grouped_dots", "ANY", plot_grouped_dots_internal)

#' @inheritParams grouped_heatmap
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
#' @rdname grouped_dots
setMethod(
    "plot_grouped_dots", "SummarizedExperiment",
    function(x, ..., groups = NULL, blocks = NULL, assay.type = "logcounts", swap_rownames = NULL) {
        x <- swap_rownames(x, swap_rownames)
        if (!is.null(blocks)) {
            blocks <- handle_column_data(object = x, blocks)
        }
        plot_grouped_dots_internal(
            x = SummarizedExperiment::assay(x, assay.type),
            groups = handle_column_data(object = x, groups),
            blocks = blocks,
            ...
        )
    }
)
