#' Plot the diffusion map and the DPT paths
#' @param x A [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment]
#'   object.
#' @param y A [DiffusionMap][destiny::DiffusionMap] or [DPT][destiny::DPT]
#'   object.
#' @param ... For `plot_diffusion_map`, other arguments passed to
#'   [plotReducedDim][scater::plotReducedDim]; for `add_dpt_paths`, other
#'   arguments passed to [geom_path][ggplot2::geom_path].
#' @name plot_diffusion_map
#' @return For `plot_diffusion_map`, A ggplot2 object depicting the Diffusion
#'   map, if add_paths is `TRUE`, DPT path of the root is also added in this
#'   map. for `add_dpt_paths`, a ggplot2 geom object is returned.
NULL

##########################################################################
#' @export
#' @rdname plot_diffusion_map
setGeneric(
    "plot_diffusion_map",
    function(x, y, ...) standardGeneric("plot_diffusion_map"),
    signature = c("x", "y")
)

#' @export
#' @rdname plot_diffusion_map
#' @importClassesFrom destiny DiffusionMap
setMethod(
    "plot_diffusion_map",
    c("SingleCellExperiment", "DiffusionMap"),
    function(x, y, ...) {
        evs <- destiny::eigenvectors(y)
        SingleCellExperiment::reducedDim(x, "DiffusionMap") <- evs
        scater::plotReducedDim(
            x,
            dimred = "DiffusionMap",
            ...
        )
    }
)

#' @param dimred A string or integer scalar indicating the reduced dimension
#'   result in `reducedDims(x)` to plot.
#' @param colour_by Specification of a column metadata field or a feature to
#'   colour by, see the by argument in
#'   [`?retrieveCellInfo`][scater::retrieveCellInfo()] for possible values. If y
#'   is a [DPT][destiny::DPT] object, this can also be "DPT" (or starting by
#'   "DPT", such as "DPT1", "DPT2" ...), "branch" or "Branch" e.g. which will be
#'   extracted from `y`. See also [plot.DPT][destiny::plot.DPT].
#' @param add_paths A scalar logical value indicates whether to add path.
#' @param root Root branch ID. Will be used as the start of the DPT. (default:
#'   lowest branch ID).
#' @param paths_to Numeric Branch IDs, which are used as target(s) for the
#'   path(s) to draw.
#' @param path_args A list of arguments to passed to
#'   [geom_path][ggplot2::geom_path]. Elements of the list can be the same
#'   length of `paths_to`. If you want to pass a complicated oject into
#'   `path_args`, try to wrap it with [list][base::list] function, since the
#'   internal use `rep_len` parallelly to get the same length of `paths_to`. for
#'   example `list(grid::arrow(ends = "last"))`.
#' @param w_width Window width for smoothing the path (see
#'   [smth.gaussian][smoother::smth.gaussian]).
#' @param ncomponents A numeric scalar indicating the number of dimensions to
#'   plot, starting from the first dimension. Alternatively, a numeric vector
#'   specifying the dimensions to be plotted.
#' @rdname plot_diffusion_map
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom destiny DPT
#' @export
setMethod(
    "plot_diffusion_map",
    c("SingleCellExperiment", "DPT"),
    function(x, y, dimred = "DiffusionMap", colour_by = "DPT",
             add_paths = TRUE, root = NULL, paths_to = NULL,
             path_args = list(),
             w_width = 0.1, ..., ncomponents = 2L) {
        if (!rlang::is_scalar_logical(add_paths)) {
            cli::cli_abort("{.arg add_paths} must be a scalar logical value.")
        }

        is_colour_by_destiny_data <- grepl("^DPT", colour_by, perl = TRUE) ||
            identical(colour_by, "Branch") ||
            identical(colour_by, "branch")

        if (is_colour_by_destiny_data || add_paths) {
            assert_length(root, 1L, null_ok = TRUE)
            assert_class(root, is.numeric, "numeric", null_ok = TRUE)

            # extract branch data
            branch_data <- y@branch[, 1L, drop = TRUE]

            # extract DPT data
            dpt_data <- as.data.frame(t(y[destiny::tips(y)]),
                make.names = FALSE
            )
            if (is.null(root)) {
                root <- min(branch_data, na.rm = TRUE)
            } else {
                root <- as.integer(root)
            }
            dpt_data <- dpt_data[[root]]
        }

        if (is_colour_by_destiny_data) {
            colour_by <- rlang::exec(
                quote(S4Vectors::DataFrame),
                !!colour_by := switch(colour_by,
                    DPT = dpt_data,
                    branch = ,
                    Branch = branch_data,
                    y[[colour_by]]
                ),
                check.names = FALSE
            )
        }
        if (add_paths || identical(dimred, "DiffusionMap")) {
            dm <- y@dm
            if (!methods::is(dm, "DiffusionMap")) {
                cli::cli_abort("{.field dm} slot in {.arg y} must be a {.cls DiffusionMap} object.")
            }
        }

        if (identical(dimred, "DiffusionMap")) {
            plot_res <- plot_diffusion_map(
                x, dm,
                colour_by = colour_by,
                ncomponents = ncomponents,
                ...
            )
        } else {
            plot_res <- scater::plotReducedDim(
                x,
                dimred = dimred,
                colour_by = colour_by,
                ncomponents = ncomponents,
                ...
            )
        }

        if (add_paths) {
            plot_res <- plot_res + rlang::exec(
                add_dpt_paths,
                x = x, y = y, root = root, paths_to = paths_to,
                w_width = w_width, ncomponents = ncomponents,
                !!!path_args
            )
        }
        plot_res
    }
)
######################################################################


######################################################################
#' @rdname plot_diffusion_map
#' @export
setGeneric(
    "add_dpt_paths",
    function(x, y, ...) standardGeneric("add_dpt_paths"),
    signature = c("x", "y")
)

#' @param colour The color of the paths.
#' @param color The color of the paths.
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom destiny DPT
#' @rdname plot_diffusion_map
#' @export
setMethod(
    "add_dpt_paths",
    c("SingleCellExperiment", "DPT"),
    function(x, y, root = NULL, paths_to = NULL,
             colour = c(
                 "gray62", "#F5C710", "#CD0BBC", "#28E2E5",
                 "#2297E6", "#61D04F", "#DF536B", "black"
             ),
             color = colour,
             w_width = 0.1, ..., ncomponents = 2L) {

        assert_length(root, 1L, null_ok = TRUE)
        assert_class(root, is.numeric, "numeric", null_ok = TRUE)

        # extract branch data
        branch_data <- y@branch[, 1L, drop = TRUE]

        # extract DPT data
        dpt_data <- as.data.frame(t(y[destiny::tips(y)]),
            make.names = FALSE
        )
        if (is.null(root)) {
            root <- min(branch_data, na.rm = TRUE)
        } else {
            root <- as.integer(root)
        }
        dpt_data <- dpt_data[[root]]

        dm <- y@dm
        if (!methods::is(dm, "DiffusionMap")) {
            cli::cli_abort("{.field dm} slot in {.arg y} must be a {.cls DiffusionMap} object.")
        }
        if (length(ncomponents) == 1L) {
            ncomponents <- seq_len(ncomponents)
        }
        evs <- destiny::eigenvectors(dm)[, ncomponents, drop = FALSE]
        geom_dpt_paths_helper(
            root = root, paths_to = paths_to,
            branch = branch_data,
            dpt = dpt_data, evs = evs,
            w_width = w_width,
            path_args = rlang::list2(color = color, ...)
        )
    }
)
#####################################################################

geom_dpt_paths_helper <- function(root = integer(), paths_to, branch, dpt, evs, w_width = 0.1, path_args = list()) {
    if (is.null(paths_to)) {
        paths_to <- setdiff(branch[!is.na(branch)], root)
    }
    if (length(paths_to) > 0L) {
        args_list <- list(path_id = paths_to)
    } else {
        return(NULL)
    }
    evs <- as.data.frame(evs)
    if (ncol(evs) >= 2L) {
        evs <- evs[1:2]
        names(evs) <- c("x", "y")
    } else {
        cli::cli_abort("Cannot plot with less than 2 {.arg ncomponents}")
    }
    if (length(path_args) > 0L) {
        path_args <- lapply(path_args, rep_len, length.out = length(paths_to))
        args_list <- c(args_list, path_args)
    }
    .mapply(function(path_id, ...) {
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
    }, args_list, MoreArgs = NULL)
}
