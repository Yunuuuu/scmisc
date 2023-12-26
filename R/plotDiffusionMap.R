#' Plot the diffusion map and the DPT paths
#'
#' @param object A
#'   [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] object.
#' @param root Root branch ID. Will be used as the start of the DPT. (default:
#'   lowest branch ID)
#' @param paths_to Numeric Branch IDs, which are used as target(s) for the
#'   path(s) to draw.
#' @param add_paths A scalar logical value indicates whether to add path.
#' @param path_args A list of arguments to passed to
#'   [geom_path][ggplot2::geom_path].
#' @param w_width Window width for smoothing the path (see
#'   [smth.gaussian][smoother::smth.gaussian]). If `NULL`, will use the same
#'   `w_width` used when running [calculateDiffusionMap].
#' @param dpt A string, specify the `DPT` to use, usually starts with "DPT". If
#' `NULL`, will use `root` to extract the `DPT`. `colour_by`, `shape_by`,
#' `size_by`, `order_by`, `text_by` and `color_by` can use "DPT" to mapping from
#' the extracted `DPT`.
#' @inheritParams scater::plotReducedDim
#' @inheritDotParams scater::plotReducedDim -object -dimred -colour_by -shape_by -size_by -order_by -text_by -color_by -percentVar
#' @return A ggplot object
#' @export
plotDiffusionMap <- function(
    object, root = NULL, paths_to = NULL,
    add_paths = FALSE, path_args = list(), w_width = NULL, ...,
    dpt = NULL, colour_by = color_by, shape_by = NULL,
    size_by = NULL, order_by = NULL, text_by = NULL,
    color_by = NULL, dimred = "DiffusionMap", ncomponents = 2L) {
    assert_s4_class(object, "SingleCellExperiment")
    assert_bool(add_paths)
    assert_string(dpt, null_ok = TRUE)
    evs <- SingleCellExperiment::reducedDim(object, dimred)
    dpt_obj <- attr(evs, "DPT")
    branch <- dpt_obj@branch[, 1L, drop = TRUE]
    branch_values <- unique(branch[!is.na(branch)])
    assert_inclusive(root, branch_values, null_ok = TRUE)
    assert_inclusive(paths_to, branch_values, null_ok = TRUE)
    if (any(ncomponents > ncol(evs))) {
        cli::cli_abort("{.arg ncomponents} is larger than {.code ncols(reducedDim(x,{dimred}))}")
    }
    if (!methods::is(dpt_obj, "DPT")) {
        cli::cli_abort(c(
            "Cannot use {.fn plotDiffusionMap} with {.val {dimred}} DiffusionMap",
            "you must run {.fn runDM} with {.code name = {dimred}}"
        ))
    }
    root <- root %||% destiny::tips(dpt_obj)[[1L]]
    if (is.null(dpt)) {
        dpt_nm <- "DPT"
        dpt <- dpt_obj[root, , drop = TRUE]
    } else {
        dpt_nm <- dpt # nolint
        dpt <- dpt_obj[[dpt]]
        if (is.null(dpt)) {
            cli::cli_abort("Cannot find {.val {dpt_nm}}")
        }
    }
    out <- scater::plotReducedDim(
        object = object, dimred = dimred, percentVar = NULL,
        colour_by = allow_diffusion_map_cols(colour_by, dpt, branch, dpt_nm),
        shape_by = allow_diffusion_map_cols(shape_by, dpt, branch, dpt_nm),
        size_by = allow_diffusion_map_cols(size_by, dpt, branch, dpt_nm),
        order_by = allow_diffusion_map_cols(order_by, dpt, branch, dpt_nm),
        text_by = allow_diffusion_map_cols(text_by, dpt, branch, dpt_nm),
        ...
    )
    if (add_paths) {
        w_width <- w_width %||% attr(dpt_obj, "w_width")
        if (length(ncomponents) == 1L) {
            ncomponents <- seq_len(ncomponents)
        }
        data <- as.data.frame(
            evs[, ncomponents, drop = FALSE],
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
        if (ncol(data) >= 2L) {
            names(data)[1:2] <- c("x", "y")
        } else {
            cli::cli_abort("Cannot plot with less than 2 {.arg ncomponents}")
        }
        paths_data <- dpt_paths_data(root, paths_to,
            dpt = dpt, branch = branch, data = data, w_width = w_width
        )
        out <- out + lapply(paths_data, function(x) {
            rlang::inject(ggplot2::geom_path(
                mapping = ggplot2::aes(.data$x, .data$y),
                data = x, inherit.aes = FALSE, !!!path_args
            ))
        })
    }
    out
}

allow_diffusion_map_cols <- function(by, dpt, branch, dpt_nm) {
    if (is.null(by)) {
        return(by)
    }
    if (identical(by, "DPT") || identical(by, "dpt")) {
        by <- dpt_nm
        value <- dpt
    } else if (identical(by, "Branch") || identical(by, "branch")) {
        value <- branch
    } else {
        return(by)
    }
    out <- S4Vectors::DataFrame(nm = value)
    names(out) <- by
    out
}

dpt_paths_data <- function(root, paths_to, dpt, branch, data, w_width) {
    if (is.null(paths_to)) {
        paths_to <- setdiff(branch[!is.na(branch)], branch[root])
    }
    if (!length(paths_to)) {
        return(NULL)
    }
    .mapply(function(path_id) {
        path_idx <- which(branch %in% c(branch[root], path_id))
        path_pt <- dpt[path_idx]
        path_data <- data[path_idx, , drop = FALSE]
        path_data <- path_data[order(path_pt), ]
        path_data <- lapply(path_data, function(col) {
            smoother::smth.gaussian(col, w_width, tails = TRUE)
        })
        data.table::setDF(path_data)
    }, list(path_id = paths_to), MoreArgs = NULL)
}
