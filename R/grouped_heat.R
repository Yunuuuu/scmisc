#' Plot heatmap of group-level expression averages
#'
#' @description
#' `plot_grouped_heat` create a heatmap of average expression values for each
#' group of cells and specified features in a SingleCellExperiment object.
#'
#' @details
#' The order of the row or column will be the same with `unlist(marker_list)` or
#' `levels(groups)` depend on whether flip is `TRUE` or `FALSE` if
#' ComplexHeatmap clustering is turned off (if groups isn't a factor, the
#' internal will coerce it as a factor). So we can easily add annotation in
#' Heatmap as long as we provide a
#' [HeatmapAnnotation][ComplexHeatmap::HeatmapAnnotation] in the same order.
#'
#' @param x A numeric matrix of counts with cells in columns and features in
#' rows.
#'
#' Alternatively, a
#' [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment] or
#' [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] object
#' containing such a matrix.
#' @inheritParams annotate_clusters
#' @param marker_list A named list contaning the markers for each cell types.
#'   names indicate the cell type label and values indicate the markers for this
#'   cell type. If NULL, defaults to all features in __x__.
#' @param cluster2cell A named character or factor returned by
#' [`annotate_clusters()`][annotate_clusters].
#' @param groups A factor (or vector coercible into a factor) specifying the
#'   group to which each cell in `x` belongs. Alternatively, if `x` is a
#'   [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment], e.g,
#'   String specifying the field of `colData(x)` containing the grouping factor.
#'   In this way, if groups is `NULL`, "label" in `colData(x)` will be
#'   extracted.
#' @param colour Alias for color.
#' @param flip A scalar logical indicates whether flipping the axis, the default
#'   `FALSE` means rows are `features` specified in `marker_list` and columns
#'   are `groups`.
#' @param row_labels,column_labels The row and column labels of the heatmap. Can
#'   be One of:
#'   - `NULL` for the default labels. If flip is `FALSE`, row_labels are the
#'     same with `unlist(marker_list)` and the column_labels are the same with
#'     `levels(groups)`. Otherwise, the reverse is also true.
#'   - A character vector of labels. If named, they will be matched by indexing
#'     with the default labels, Otherwise, used as it is.
#'   - A function that takes the default labels as input and returns a character
#'    vector of the same length. Can be also a rlang
#'     [lambda][rlang::as_function()] function notation.
#' @param ... Other arguments passed to [Heatmap][ComplexHeatmap::Heatmap] and
#'   specific methods.
#' @inheritParams DotsHeatmap
#' @inheritParams scater::plotDots
#' @return A [Heatmap-class][ComplexHeatmap::Heatmap-class] object
#' @export
#' @name grouped_heatmap
#' @aliases plot_grouped_heat
setGeneric(
    "plot_grouped_heat",
    function(x, ...) standardGeneric("plot_grouped_heat")
)

plot_grouped_heat_internal <- function(
    x, marker_list = NULL, cluster2cell = NULL, groups = NULL, ...,
    blocks = NULL, colour = color, color = NULL, center = FALSE,
    scale = FALSE, zlim = NULL, flip = FALSE,
    row_labels = NULL, column_labels = NULL) {
    grouped_heat_internal(
        x = x, marker_list = marker_list,
        groups = groups, blocks = blocks,
        cluster2cell = cluster2cell,
        center = center, scale = scale,
        zlim = zlim, threshold = 0L,
        colour = colour, flip = flip,
        row_labels = row_labels, column_labels = column_labels,
        ..., graph_type = "square"
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

##################################################################
grouped_heat_internal <- function(x, marker_list = NULL, groups = NULL,
                                  blocks = NULL, cluster2cell = NULL,
                                  center = FALSE, scale = FALSE,
                                  zlim = NULL, threshold = 0L,
                                  colour = NULL, ..., flip = FALSE,
                                  row_labels = NULL, column_labels = NULL,
                                  graph_type = c("dots", "square")) {
    assert_(marker_list, is.list, "a {.cls list}", null_ok = TRUE)
    if (!is.null(marker_list)) {
        if (length(marker_list) == 0L) {
            cli::cli_abort("Empty list is not allowed in {.arg marker_list}.")
        } else if (!rlang::is_named(marker_list)) {
            cli::cli_abort("All elements in {.arg marker_list} must be named.")
        }
        if (anyDuplicated(names(marker_list))) {
            cli::cli_abort("Duplicated names found in {.arg marker_list}")
        }
    }

    if (!rlang::is_scalar_logical(flip)) {
        cli::cli_abort("{.arg flip} must be a scalar logical value.")
    }
    # define the default value
    statistics <- switch(graph_type,
        dots = c("mean", "prop.detected"),
        square = "mean"
    )
    groups <- groups %||% rep_len("Group", ncol(x))

    if (!is.null(marker_list)) {
        # calculate statistics for all markers in marker_list
        markers <- unlist(marker_list, recursive = FALSE, use.names = FALSE)
        marker_groups <- rep(names(marker_list), times = lengths(marker_list))
        marker_groups <- factor(marker_groups, names(marker_list))
        marker_groups <- data.table::fdroplevels(marker_groups)
        markers <- as.character(markers)
        if (anyNA(markers)) {
            cli::cli_abort("{.val {NA}} is found in {.arg marker_list}")
        }
        is_existed <- markers %chin% rownames(x) # nolint
        if (!all(is_existed)) {
            cli::cli_abort(c(
                "Finding {.val {sum(!is_existed)}} non-existed feature{?s} in {.arg marker_list}",
                "!" = "Non-existed feature{?s}: {.val {markers[!is_existed]}}"
            ))
        }
    } else {
        markers <- NULL
        marker_groups <- factor(rep_len("Group", nrow(x)), "Group")
    }

    stat_list <- summarize_features_by_groups(
        x = x, features = markers, groups = groups,
        statistics = statistics, blocks = blocks,
        threshold = threshold
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

    # prepare matrix for heatmap
    heat_matrix <- heat_data_list$x
    size_matrix <- stat_list$prop.detected

    col_fn <- circlize::colorRamp2(
        heat_data_list$colour_breaks,
        colors = heat_data_list$colour
    )

    # prepare column_split and row_split
    # column is groups; row is genes if flip is FALSE
    if (flip) {
        heat_matrix <- t(heat_matrix)
        size_matrix <- t(size_matrix)
        if (nlevels(marker_groups) == 1L) {
            column_split <- NULL
        } else {
            column_split <- marker_groups
        }
        if (is.null(cluster2cell)) {
            row_split <- NULL
        } else {
            row_split <- cluster2cell[rownames(heat_matrix)]
        }
        row_label_arg <- "levels(groups)"
        column_label_arg <- "unlist(marker_list)"
    } else {
        if (is.null(cluster2cell)) {
            column_split <- NULL
        } else {
            column_split <- cluster2cell[colnames(heat_matrix)]
        }
        if (nlevels(marker_groups) == 1L) {
            row_split <- NULL
        } else {
            row_split <- marker_groups
        }
        row_label_arg <- "unlist(marker_list)"
        column_label_arg <- "levels(groups)"
    }
    # prepare row_labels
    row_labels <- label_fn_helper(row_labels,
        labels = rownames(heat_matrix),
        label_arg = row_label_arg
    )

    # prepare column labels
    column_labels <- label_fn_helper(column_labels,
        labels = colnames(heat_matrix),
        label_arg = column_label_arg
    )

    if (graph_type == "square") {
        heat_obj <- ComplexHeatmap::Heatmap(
            heat_matrix,
            col = col_fn,
            row_split = row_split,
            column_split = column_split,
            row_labels = row_labels,
            column_labels = column_labels,
            ...
        )
    } else {
        heat_obj <- DotsHeatmap(
            heat_matrix,
            matrix_size = size_matrix,
            col = col_fn,
            row_split = row_split,
            column_split = column_split,
            row_labels = row_labels,
            column_labels = column_labels,
            ...
        )
    }
    heat_obj
}

#' @keywords internal
#' @noRd
label_fn_helper <- function(x, labels, arg = rlang::caller_arg(x), label_arg = rlang::caller_arg(labels), call = rlang::caller_env()) {
    if (is.character(x)) {
        if (rlang::is_named(x)) {
            x <- x[labels]
        } else if (length(x) != length(labels)) {
            cli::cli_abort(
                "An unnamed chracter vector {.arg {arg}} must have the same length of {.arg {label_arg}}",
                call = call
            )
        }
    } else if (is.function(x) ||
        rlang::is_formula(x) ||
        rlang::is_quosure(x)) {
        x <- rlang::as_function(x)(labels)
        if (length(x) != length(labels) || !is.character(x)) {
            cli::cli_abort(
                "{.fn {arg}} must returned a {.cls character} with the same length of {.arg {label_arg}}",
                call = call
            )
        }
    } else if (is.null(x)) {
        x <- labels
    } else {
        cli::cli_abort("{.arg {arg}} must be a {.cls character} or a {.cls function} or {.val {NULL}}", call = call)
    }
    x
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
            zlim <- range(x, na.rm = TRUE, finite = TRUE)
        }
    }
    if (is.null(colour)) {
        if (center) {
            colour_fn <- scales::brewer_pal(palette = "RdYlBu", direction = -1L)
        } else {
            colour_fn <- scales::viridis_pal(option = "D")
        }
        colour <- colour_fn(9L)
    }
    x[x < zlim[1L]] <- zlim[1L]
    x[x > zlim[2L]] <- zlim[2L]
    list(
        x = x, colour = colour,
        colour_breaks = seq(zlim[1L], zlim[2L], length.out = length(colour)),
        zlim = zlim
    )
}
