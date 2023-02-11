#' Plot heatmap or dots of group-level expression averages
#'
#' Create a heatmap of average expression values for each group of cells and
#' specified features in a SingleCellExperiment object.
#'
#' Create a dot plot of expression values for a grouping of cells, where the
#' size and color of each dot represents the proportion of detected expression
#' values and the average expression, respectively, for each feature in each
#' group of cells.
#' @param x A numeric matrix of counts with cells in columns and features in
#' rows.
#'
#' Alternatively, a
#' [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment] or
#' [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] object
#' containing such a matrix.
#' @inheritParams annotate_clusters
#' @param cluster2cell A named character or factor returned by
#' [`annotate_clusters()`][annotate_clusters].
#' @param groups A factor (or vector coercible into a factor) specifying the
#'   group to which each cell in `x` belongs. Alternatively, if `x` is a
#'   [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment], e.g,
#'   String specifying the field of `colData(x)` containing the grouping factor.
#'   In this way, if groups is `NULL`, "label" in `colData(x)` will be
#'   extracted.
#' @param ... Other arguments passed to [Heatmap][ComplexHeatmap::Heatmap]
#' @param colour Alias for color.
#' @param threshold Numeric value specifying the cap on the proportion of
#' detected expression values. Only used by `plot_grouped_dots`.
#' @inheritParams scater::plotDots
#' @return A [Heatmap-class][ComplexHeatmap::Heatmap-class] object
#' @name grouped_heatmap
NULL

#############################################################################
# for plot_grouped_heat
#' @export
#' @rdname grouped_heatmap
setGeneric(
    "plot_grouped_heat",
    function(x, ...) standardGeneric("plot_grouped_heat")
)

plot_grouped_heat_internal <- function(
    x, marker_list, cluster2cell = NULL, groups = NULL, ...,
    blocks = NULL, colour = color, color = NULL, center = FALSE,
    scale = FALSE, zlim = NULL) {
    grouped_heat_internal(
        x = x, marker_list = marker_list,
        groups = groups, blocks = blocks,
        cluster2cell = cluster2cell,
        center = center, scale = scale,
        zlim = zlim, threshold = 0L,
        colour = colour, ...,
        graph_type = "square"
    )
}

#' @export
#' @rdname grouped_heatmap
setMethod("plot_grouped_heat", "ANY", plot_grouped_heat_internal)

#' @export
#' @param assay.type A string or integer scalar indicating which `assays` in the
#'   `x` contains the count matrix.
#' @param swap_rownames A characters (or vector coercible into a characters)
#'    used to identify features instead of `rownames(x)` when labelling plot
#'   elements.  Alternatively, if `x` is a
#'   [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment], e.g,
#'   String specifying the field of `rowData(x)` to be used.
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @rdname grouped_heatmap
setMethod(
    "plot_grouped_heat", "SummarizedExperiment",
    function(x, ..., groups = NULL, blocks = NULL, assay.type = "logcounts", swap_rownames = NULL) {
        x <- swap_rownames(x, swap_rownames)
        if (!is.null(blocks)) {
            blocks <- handle_column_data(object = x, blocks)
        }
        plot_grouped_heat_internal(
            x = SummarizedExperiment::assay(x, assay.type),
            groups = handle_column_data(object = x, groups),
            blocks = blocks,
            ...
        )
    }
)

##############################################################################
# for plot_grouped_dots
#' @export
#' @rdname grouped_heatmap
setGeneric(
    "plot_grouped_dots",
    function(x, ...) standardGeneric("plot_grouped_dots")
)

plot_grouped_dots_internal <- function(
    x, marker_list, cluster2cell = NULL, groups = NULL, ...,
    blocks = NULL, colour = color, color = NULL, center = FALSE,
    scale = FALSE, threshold = 0L, zlim = NULL) {
    grouped_heat_internal(
        x = x, marker_list = marker_list,
        groups = groups, blocks = blocks,
        cluster2cell = cluster2cell,
        center = center, scale = scale,
        zlim = zlim, threshold = threshold,
        colour = colour, ...,
        graph_type = "dots"
    )
}

#' @export
#' @rdname grouped_heatmap
setMethod("plot_grouped_dots", "ANY", plot_grouped_dots_internal)

#' @export
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @rdname grouped_heatmap
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

##################################################################
grouped_heat_internal <- function(x, marker_list, groups = NULL,
                                  blocks = NULL, cluster2cell = NULL,
                                  center = FALSE, scale = FALSE,
                                  zlim = NULL, threshold = 0L,
                                  colour = NULL, ...,
                                  graph_type = c("dots", "square")) {
    assert_class(marker_list, is.list, "list", null_ok = FALSE)
    if (length(marker_list) > 0L) {
        if (all(!has_names(marker_list))) {
            cli::cli_abort("All elements in {.arg marker_list} must be named.")
        }
    } else {
        cli::cli_abort("Empty list is not allowed in {.arg marker_list}.")
    }
    # define the default value
    statistics <- switch(graph_type,
        dots = c("mean", "prop.detected"),
        square = "mean"
    )
    groups <- groups %||% rep_len("Group", ncol(x))

    # calculate statistics for all markers in marker_list
    markers <- unlist(marker_list, recursive = FALSE, use.names = FALSE)
    marker_groups <- rep(names(marker_list), times = lengths(marker_list))
    marker_groups <- factor(marker_groups, names(marker_list))
    stat_list <- summarize_features_by_groups(
        x = x, features = markers, groups = groups,
        statistics = statistics, blocks = blocks,
        threshold = threshold,
        check_dup = FALSE
    )$statistics

    # rows are genes
    # columns are clusters
    heat_data_list <- heatmap_scale(
        stat_list$mean,
        center = center,
        scale = scale,
        colour = colour,
        zlim = zlim
    )
    col_fn <- circlize::colorRamp2(
        heat_data_list$colour_breaks,
        colors = heat_data_list$colour
    )
    layer_fn <- switch(graph_type,
        dots = function(j, i, x, y, width, height, fill) {
            size <- ComplexHeatmap::pindex(
                stat_list$prop.detected,
                i = i, j = j
            )
            grid::grid.circle(
                x = x, y = y,
                r = abs(size) / 2 * min(grid::unit.c(width, height)),
                gp = grid::gpar(fill = fill, col = NA)
            )
        },
        square = NULL
    )

    if (is.null(cluster2cell)) {
        column_split <- NULL
    } else {
        column_split <- cluster2cell[colnames(heat_data_list$x)]
    }
    if (nlevels(marker_groups) == 1L) {
        row_split <- NULL
    } else {
        row_split <- marker_groups
    }
    # column is groups; row is genes
    ComplexHeatmap::Heatmap(
        heat_data_list$x,
        col = col_fn,
        row_labels = markers,
        row_split = row_split,
        column_labels = colnames(heat_data_list$x),
        column_split = column_split,
        layer_fun = layer_fn,
        ...
    )
}

heatmap_scale <- function(x, center, scale, colour = NULL, zlim = NULL) {
    if (center) {
        x <- x - rowMeans(x, na.rm = TRUE)
    }
    if (scale) {
        x <- x / sqrt(rowSums(x^2L, na.rm = TRUE) / (ncol(x) - 1L))
    }
    if (is.null(zlim)) {
        if (center) {
            extreme <- max(abs(x), na.rm = TRUE)
            zlim <- extreme * c(-1L, 1L)
        } else {
            zlim <- range(x, na.rm = TRUE)
        }
    }
    if (is.null(colour)) {
        if (center) {
            colour <- rev(RColorBrewer::brewer.pal(9L, "RdYlBu"))
        } else {
            colour <- viridis::viridis(9L)
        }
    }
    x[x < zlim[1L]] <- zlim[1L]
    x[x > zlim[2L]] <- zlim[2L]
    list(
        x = x, colour = colour,
        colour_breaks = seq(zlim[1L], zlim[2L], length.out = length(colour)),
        zlim = zlim
    )
}
