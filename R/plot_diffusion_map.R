#' Plot the diffusion map
#' @param x A [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment]
#'   object with `DiffusionMap` reducedDim component and `DPT` component if
#'   add_paths is `TRUE`.
#' @param ncomponents A numeric scalar indicating the number of dimensions to
#'   plot, starting from the first dimension. Alternatively, a numeric vector
#'   specifying the dimensions to be plotted.
#' @param colour_by Specification of a column metadata field or a feature to
#'   colour by, see the by argument in
#'   [`?retrieveCellInfo`][scater::retrieveCellInfo()] for possible values. If
#'   "DPT", `x` must contain components "destiny_branch" and "destiny_dpt", and
#'   the DPT of root will be used.
#' @param add_paths A scalar logical value indicates whether to add path.
#' @param root Root branch ID. Will be used as the start of the DPT. (default:
#'   lowest branch ID).
#' @param paths_to Numeric Branch IDs. Are used as target(s) for the path(s) to
#'   draw.
#' @param path_args A list of arguments to passed to
#'   [geom_path][ggplot2::geom_path].
#' @param w_width Window width for smoothing the path (see
#'   [smth.gaussian][smoother::smth.gaussian]).
#' @param ... For `plot_diffusion_map`, Other arguments passed to
#'   [plotReducedDim][scater::plotReducedDim]; for `add_dpt_paths`, other
#'   arguments passed to [geom_path][ggplot2::geom_path]; and for
#'   `plot_dpt_paths`, other arguments passed to [add_dpt_paths]
#' @name plot_diffusion_map
#' @return 
#'   - `plot_diffusion_map`: A ggplot2 object depicting the Diffusion map, if
#'     add_paths is `TRUE`, DPT of the root is also added in this map.
#'   - `plot_dpt_paths`: A ggplot2 object depicting the DPT paths.
#'   - `add_dpt_paths`: A ggplot2 geom object containing the DPT paths.
NULL

#' @rdname plot_diffusion_map
#' @export
plot_diffusion_map <- function(x, colour_by = "DPT", add_paths = TRUE, root = NULL, paths_to = NULL, path_args = list(), w_width = 0.1, ..., ncomponents = 2L) {
    if (!methods::is(x, "SingleCellExperiment")) {
        cli::cli_abort("{.arg x} must be a {.cls SingleCellExperiment}")
    }
    if (!rlang::is_scalar_logical(add_paths)) {
        cli::cli_abort("{.arg add_paths} must be a scalar logical value.")
    }
    if (identical(x, "DPT") || add_paths) {
        assert_length(root, 1L, null_ok = TRUE)
        assert_class(root, is.numeric, "numeric", null_ok = TRUE)
        branch_data <- x$destiny_branch
        dpt_data <- x$destiny_dpt
        if (is.null(branch_data) || is.null(dpt_data)) {
            cli::cli_abort("{.arg x} must contain {.field destiny_branch} and {.field destiny_dpt}")
        }
        if (is.null(root)) {
            root <- min(branch_data, na.rm = TRUE)
        } else {
            root <- as.integer(root)
        }
        dpt_data <- dpt_data[[root]]
        if (identical(x, "DPT")) {
            colour_by <- S4Vectors::DataFrame(DPT = dpt_data)
        }
    }
    dm_plot <- scater::plotReducedDim(
        x,
        dimred = "DiffusionMap",
        ncomponents = ncomponents,
        colour_by = colour_by,
        ...
    )
    if (add_paths) {
        if (length(ncomponents) == 1L) {
            ncomponents <- seq_len(ncomponents)
        }
        evs <- SingleCellExperiment::reducedDim(
            x, "DiffusionMap"
        )[, ncomponents, drop = FALSE]
        dm_plot <- dm_plot + rlang::inject(
            geom_dpt_paths_helper(
                root = root, paths_to = paths_to,
                branch = branch_data,
                dpt = dpt_data,
                evs = evs,
                w_width = w_width,
                !!!path_args
            )
        )
    }
    dm_plot
}

#' @rdname plot_diffusion_map
#' @export
add_dpt_paths <- function(x, root = NULL, paths_to = NULL, w_width = 0.1, ..., ncomponents = 2L) {
    if (!methods::is(x, "SingleCellExperiment")) {
        cli::cli_abort("{.arg x} must be a {.cls SingleCellExperiment}")
    }
    assert_length(root, 1L, null_ok = TRUE)
    assert_class(root, is.numeric, "numeric", null_ok = TRUE)
    branch_data <- x$destiny_branch
    dpt_data <- x$destiny_dpt
    if (is.null(branch_data) || is.null(dpt_data)) {
        cli::cli_abort("{.arg x} must contain {.field destiny_branch} and {.field destiny_dpt}")
    }
    if (is.null(root)) {
        root <- min(branch_data, na.rm = TRUE)
    } else {
        root <- as.integer(root)
    }
    dpt_data <- dpt_data[[root]]
    if (!any("DiffusionMap", SingleCellExperiment::reducedDimNames(x))) {
        cli::cli_abort("Cannot find {.val DiffusionMap} in the {.field reducedDim}{.arg x}")
    }
    if (length(ncomponents) == 1L) {
        ncomponents <- seq_len(ncomponents)
    }
    evs <- SingleCellExperiment::reducedDim(
        x, "DiffusionMap"
    )[, ncomponents, drop = FALSE]
    geom_dpt_paths_helper(
        root = root, paths_to = paths_to,
        branch = branch_data,
        dpt = dpt_data,
        evs = evs,
        w_width = w_width
    )
}

#' @rdname plot_diffusion_map
#' @export 
plot_dpt_paths <- function(...) {
    ggplot2::ggplot() + add_dpt_paths(...)
}

geom_dpt_paths_helper <- function(root = integer(), paths_to, branch, dpt, evs, w_width = 0.1, ...) {
    if (is.null(paths_to)) {
        paths_to <- setdiff(branch[!is.na(branch)], root)
    }
    evs <- as.data.frame(evs)
    if (ncol(evs) >= 2L) {
        names(evs)[1:2] <- c("x", "y")
    } else {
        cli::cli_abort("Cannot plot with less than 2 {.arg ncomponents}")
    }
    paths <- lapply(paths_to, function(path_id) {
        path_idx <- which(branch %in% c(root, path_id))
        path_pt <- dpt[path_idx]
        path_evs <- evs[path_idx, , drop = FALSE]
        path_data <- path_evs[order(path_pt), ]
        path_data <- lapply(path_data, function(col) {
            smoother::smth.gaussian(col, w_width, tails = TRUE)
        })
        data.table::setDF(path_data)
        ggplot2::geom_path(
            ggplot2::aes(.data$x, .data$y), # nolint
            data = path_data,
            inherit.aes = FALSE,
            ...
        )
    })
    paths
}
