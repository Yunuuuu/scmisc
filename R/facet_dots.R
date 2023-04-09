#' Create a facet dot plot of expression values
#'
#' Create a dot plot of expression values for a grouping of cells, where the
#' size and color of each dot represents the proportion of detected expression
#' values and the average expression, respectively, for each feature in each
#' group of cells.
#' @param x Currently, only
#'   [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] object
#'   is supported.
#' @inheritParams annotate_clusters
#' @param clusters A factor (or vector coercible into a factor) specifying the
#'   group to which each cell in `x` belongs. Alternatively, if `x` is a
#'   [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment], e.g,
#'   String specifying the field of `colData(x)` containing the grouping factor.
#'   In this way, if clusters is `NULL`, "label" in `colData(x)` will be
#'   extracted.
#' @param cluster2cell A named character or factor returned by
#'   [`annotate_clusters()`][annotate_clusters].
#' @param flip A scalar logical indicates whether flipping the plot.
#' @param facet_args A named arguments list passed to
#'   [facet_grid][ggplot2::facet_grid].
#' @inheritDotParams scater::plotDots -object -features -group -other_fields
#' @return A ggplot2 object.
#' @name facet_dots
NULL

#' @keywords internal
facet_dots_internal <- function(x, marker_list, clusters, cluster2cell = NULL, flip = TRUE, facet_args = list(scales = "free", space = "free"), ...) {
    assert_class(marker_list, is.list,
        "{.cls list} object",
        null_ok = FALSE
    )
    if (length(marker_list) == 0L) {
        cli::cli_abort("Empty list is not allowed in {.arg marker_list}.")
    } else if (!rlang::is_named(marker_list)) {
        cli::cli_abort("All elements in {.arg marker_list} must be named.")
    } else if (anyDuplicated(names(marker_list))) {
        cli::cli_abort("Duplicated names found in {.arg marker_list}")
    }
    markers <- unlist(marker_list, recursive = FALSE, use.names = FALSE)
    dup_markers <- unique(markers[duplicated(markers)])
    if (length(dup_markers) > 0L) {
        cli::cli_abort(c(
            "Duplicated markers are provided in {.arg marker_list}",
            "!" = "Duplicated marker{?s}: {.val {dup_markers}}",
            i = "Try to use {.fn plot_grouped_dots} if you want to diplay duplicated markers."
        ))
    }
    gene2cell <- structure(
        factor(
            rep(names(marker_list), times = lengths(marker_list)),
            levels = names(marker_list)
        ),
        names = markers
    )
    base_plot <- scater::plotDots(
        x,
        features = markers,
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
    facet_obj <- do.call(ggplot2::facet_grid, facet_args_list)
    base_plot + facet_obj + ggplot2::theme(strip.clip = "off")
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
    function(x, ..., clusters = NULL) {
        facet_dots_internal(x,
            clusters = handle_column_data(object = x, clusters),
            ...
        )
    }
)
