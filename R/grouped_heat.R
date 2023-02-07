#' Plot heatmap of group-level expression averages
#'
#' Create a heatmap of average expression values for each group of cells and
#' specified features in a SingleCellExperiment object.
#' @inheritParams facet_dots
#' @param blocks A factor (or vector coercible into a factor) specifying the
#'   blocking factor to which each cell in `x` belongs (e.g., batch of origin).
#'   Alternatively, if `x` is a
#'   [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment], e.g,
#'   String specifying the field of `colData(x)` containing the grouping factor.
#'   In this way, if clusters is `NULL`, "label" in `colData(x)` will be
#'   extracted.
#' @param colour Alias for color.
#' @param ... Other arguments passed to [Heatmap][ComplexHeatmap::Heatmap]
#' @inheritParams scater::plotDots
#' @name plot_grouped_heat
NULL

grouped_heat_internal <- function(
    x, marker_list, clusters, cluster2cell = NULL,
    blocks = NULL, center = FALSE,
    scale = FALSE, zlim = NULL,
    colour = color, color = NULL,
    ...) {
    if (is.null(colnames(x))) {
        colnames(x) <- seq_len(ncol(x))
    }
    ids <- S4Vectors::DataFrame(clusters = clusters)
    if (!is.null(blocks)) {
        ids$blocks <- blocks
    }
    gene2cell <- structure(
        factor(
            rep(names(marker_list), times = lengths(marker_list)),
            levels = names(marker_list)
        ),
        names = unlist(marker_list, recursive = FALSE, use.names = FALSE)
    )
    heat_se <- scuttle::summarizeAssayByGroup(
        x, ids,
        subset.row = names(gene2cell),
        statistic = "mean"
    )
    if (!is.null(blocks)) {
        heat_matrix <- scuttle::correctGroupSummary(
            heat_se,
            group = heat_se$group,
            block = heat_se$blocks
        )
    } else {
        heat_matrix <- SummarizedExperiment::assay(heat_se)
    }
    # rows are genes
    # columns are clusters
    heat_data_list <- heatmap_scale(
        heat_matrix,
        center = center,
        scale = scale,
        colour = colour,
        zlim = zlim
    )
    ComplexHeatmap::Heatmap(
        heat_data_list$x,
        col = circlize::colorRamp2(
            heat_data_list$colour_breaks,
            colors = heat_data_list$colour
        ),
        row_split = gene2cell[rownames(heat_data_list$x)],
        column_split = cluster2cell[colnames(heat_data_list$x)],
        ...
    )
}

#' @export
#' @rdname plot_grouped_heat
setGeneric(
    "plot_grouped_heat",
    function(x, ...) standardGeneric("plot_grouped_heat")
)

#' @export
#' @rdname plot_grouped_heat
setMethod("plot_grouped_heat", "ANY", grouped_heat_internal)

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
#' @rdname plot_grouped_heat
setMethod(
    "plot_grouped_heat", "SummarizedExperiment",
    function(x, ..., clusters = NULL, blocks = NULL, assay.type = "logcounts", swap_rownames = NULL) {
        x <- swap_rownames(x, swap_rownames)
        if (!is.null(blocks)) {
            blocks <- handle_column_data(object = x, blocks)
        }
        grouped_heat_internal(
            SummarizedExperiment::assay(x, assay.type),
            clusters = handle_column_data(object = x, clusters),
            blocks = blocks,
            ...
        )
    }
)

# Generated from function body. Editing this file has no effect.
heatmap_scale <- function(x, center, scale, colour = NULL, zlim = NULL) {
    if (center) {
        x <- x - rowMeans(x)
    }
    if (scale) {
        x <- x / sqrt(rowSums(x^2L) / (ncol(x) - 1L))
    }
    if (is.null(zlim)) {
        if (center) {
            extreme <- max(abs(x), na.rm = TRUE)
            zlim <- extreme * c(-1L, 1L)
        } else {
            zlim <- range(x)
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
        colour_breaks = seq(zlim[1L], zlim[2L],
            length.out = length(colour) + 1L
        ),
        zlim = zlim
    )
}
