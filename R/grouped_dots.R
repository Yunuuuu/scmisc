#' Plot dots heatmap of group-level expression averages
#'
#' @description
#' `plot_grouped_dots` create a dot plot of expression values for a grouping of
#' cells, where the size and color of each dot represents the proportion of
#' detected expression values and the average expression, respectively, for each
#' feature in each group of cells.
#'
#' @inherit grouped_heatmap details
#' @param threshold Numeric value specifying the cap on the proportion of
#'   detected expression values.
#' @param ... Other arguments passed to [DotsHeatmap] and specific methods.
#' @inheritParams DotsHeatmap
#' @inherit DotsHeatmap return
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
    dots_size_legend_param = list(),
    flip = FALSE, row_labels = NULL, column_labels = NULL) {
    if (!is.null(dots_size_legend_param)) {
        dots_size_legend_param$title <- dots_size_legend_param$title %||%
            "NumDetected"
    }
    grouped_heat_internal(
        x = x, marker_list = marker_list,
        groups = groups, blocks = blocks,
        cluster2cell = cluster2cell,
        center = center, scale = scale,
        zlim = zlim, threshold = threshold,
        colour = colour, flip = flip,
        row_labels = row_labels, column_labels = column_labels,
        dots_size_legend_param = dots_size_legend_param,
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
