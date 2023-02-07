#' Create a facet dot plot of expression values
#'
#' Create a dot plot of expression values for a grouping of cells, where the
#' size and color of each dot represents the proportion of detected expression
#' values and the average expression, respectively, for each feature in each
#' group of cells.
#' @param x Currently, only
#' [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] object is
#' supported.
#' @inheritParams annotate_clusters
#' @param clusters A factor (or vector coercible into a factor) specifying the
#'   group to which each cell in `x` belongs. Alternatively, String specifying
#'   the field of `colData(x)` containing the grouping factor if `x` is a
#'   [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment], e.g.
#' @param cluster2cell A named character or factor returned by
#' [`annotate_clusters()`][annotate_clusters].
#' @param flip A scalar logical indicates whether flipping the plot.
#' @param facet_args A named list passed to [facet_grid][ggplot2::facet_grid].
#' @param ... Other arguments passed to [plotDots][scater::plotDots].
#' @return A ggplot2 object.
#' @name facet_dots
NULL

#' @keywords internal
facet_dots_internal <- function(x, marker_list, clusters = NULL, cluster2cell = NULL, flip = TRUE, facet_args = list(scales = "free", space = "free"), ...) {
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
    gene2cell <- structure(
        factor(
            rep(names(marker_list), times = lengths(marker_list)),
            levels = names(marker_list)
        ),
        names = unlist(marker_list,
            recursive = FALSE, use.names = FALSE
        )
    )
    base_plot <- scater::plotDots(
        x,
        features = names(gene2cell),
        group = I(clusters),
        other_fields = list(S4Vectors::DataFrame(
            ..marker_celltypes.. = unname(gene2cell[rownames(x)])
        )),
        ...
    )
    facet_args_list <- list(ggplot2::vars(..marker_celltypes..)) # nolint
    if (!is.null(cluster2cell)) {
        ..cluster2cell.. <- cluster2cell
        facet_args_list <- c(
            facet_args_list,
            list(ggplot2::vars(..cluster2cell..[as.character(Group)])) # nolint
        )
    }
    if (flip) {
        names(facet_args_list) <- c("cols", "rows")[seq_along(facet_args_list)]
        base_plot <- base_plot + ggplot2::coord_flip()
    } else {
        names(facet_args_list) <- c("rows", "cols")[seq_along(facet_args_list)]
    }
    facet_args_list <- c(facet_args_list, facet_args)
    base_plot +
        rlang::exec(ggplot2::facet_grid, !!!facet_args_list) +
        ggplot2::theme(strip.clip = "off")
}

utils::globalVariables(c("..marker_celltypes..", "Group"))

#' @export
#' @rdname facet_dots
setGeneric(
    "facet_dots",
    function(x, ...) standardGeneric("facet_dots")
)

#' @export
#' @rdname facet_dots
setMethod("facet_dots", "ANY", facet_dots_internal)

#' @export
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @rdname facet_dots
setMethod(
    "facet_dots", "SingleCellExperiment",
    facet_dots_internal
)
