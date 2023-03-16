#' Constructor method for DotsHeatmap class
#' @param matrix A matrix. Either numeric or character. If it is a simple
#'  vector, it will be converted to a one-column matrix. Only used for dot
#'  color.
#' @param matrix_size A numeric matrix or vector. If it is a simple vector, it
#'  will be converted to a one-column matrix. Only used for dot size.
#' @param circle_gp Graphic parameters for drawing point circle (for heatmap
#'   body). The value should be specified by gpar and fill parameter is ignored.
#' @param radius_range A numeric vector of length one or two specifying the
#'  minimal and maximal dots radius (measured in "mm"). If length one, all dots
#'  will have the same size. If length two with different values, the
#'  `matrix_size` will be scaled into `radius_range`. Default value: `c(1, 6)`.
#' @param dots_size_legend_param Other arguments passed to
#'  [Legend][ComplexHeatmap::Legend]. This will be used to define the dots size
#'  legend. Two extra arguments (`radius_range` and `radius`) are added to
#'  enhance the circle legend creation. `radius_range` define the minimal and
#'  maximal legend dots radius (measured in "mm") and an [unit][grid::unit] or
#'  numeric `radius` directly specify the dots radius. If both `radius_range`
#'  and `radius` are specified, only `radius_range` will be used. `radius` must
#'  have the same length of `at`, the default `at` is created with
#'  [breaks_extended][scales::breaks_extended].  `dots_size_legend_param` can
#'  also be `NULL`, in this way, no legend for dots size will be created and a
#'  [Heatmap][ComplexHeatmap::Heatmap] object will be returned.
#' @param slice_border_gp Graphic parameters for drawing slice rectangles. The
#'   value should be specified by [gpar][grid::gpar] and fill parameter is
#'   always setted to "transparent". If `NULL`, no slice border will be drawn.
#' @inheritDotParams ComplexHeatmap::Heatmap -matrix -rect_gp -layer_fun -cell_fun
#' @return If dots_size_legend_param is a list, a [DotsHeatmap] object will be
#'   returned, if dots_size_legend_param is NULL, a
#'   [Heatmap][ComplexHeatmap::Heatmap] object.
#' @name DotsHeatmap
#' @export
DotsHeatmap <- function(matrix, matrix_size = NULL, circle_gp = gpar(col = NA), radius_range = c(1, 5), dots_size_legend_param = list(), slice_border_gp = NULL, ...) {
    assert_class(radius_range, "numeric", is.numeric)
    if (!data.table::between(length(radius_range), 1L, 2L)) {
        cli::cli_abort("{.arg radius_range} must be a numeric with a length 1 or 2")
    }
    if (anyNA(radius_range) || any(radius_range < 0L)) {
        cli::cli_abort("{.arg radius_range} must not be negative or {.val NA}")
    }
    assert_class(dots_size_legend_param, "list", is.list, null_ok = TRUE)
    assert_class(slice_border_gp, "gpar", null_ok = TRUE)
    matrix_size <- matrix_size %||% matrix
    assert_class(matrix_size, "numeric", is.numeric)
    if (!is.matrix(matrix)) {
        if (is.atomic(matrix)) {
            cli::cli_alert_info("convert simple vector to one-column matrix")
            matrix <- matrix(matrix, ncol = 1L)
        } else {
            cli::cli_abort("{.arg matrix} must be a matrix or a simple vector.")
        }
    }
    if (!is.matrix(matrix_size)) {
        if (is.atomic(matrix_size)) {
            cli::cli_alert_info("convert simple vector to one-column matrix")
            matrix_size <- matrix(matrix_size, ncol = 1L)
        } else {
            cli::cli_abort("{.arg matrix_size} must be a matrix or a simple vector.")
        }
    }
    if (!all(dim(matrix) == dim(matrix_size))) {
        cli::cli_abort("{.arg matrix} and {.arg matrix_size} must have the same dimensions")
    }
    scale_size <- scale_range(unclass(radius_range), matrix_size)
    if (!is.null(slice_border_gp)) {
        slice_border_gp$fill <- "transparent"
    }
    heatmap <- ComplexHeatmap::Heatmap(
        matrix = matrix,
        rect_gp = gpar(type = "none"),
        layer_fun = function(j, i, x, y, width, height, fill) {
            size_values <- ComplexHeatmap::pindex(matrix_size, i = i, j = j)
            circle_gp$fill <- fill
            grid::grid.circle(
                x = x, y = y,
                r = scale_size(size_values),
                gp = circle_gp
            )
            if (!is.null(slice_border_gp)) {
                grid::grid.rect(gp = slice_border_gp)
            }
        },
        ...
    )
    if (is.null(dots_size_legend_param)) {
        return(heatmap)
    } else {
        dots_size_legend_param <- do.call(
            prepare_legend_param_for_dots_size,
            c(dots_size_legend_param, list(
                matrix = matrix_size, scale_size = scale_size
            ))
        )
    }

    methods::new("DotsHeatmap",
        heatmap = heatmap,
        dots_legend = do.call(ComplexHeatmap::Legend, dots_size_legend_param)
    )
}

prepare_legend_param_for_dots_size <- function(matrix, at = NULL, radius_range = NULL, radius = NULL, ..., scale_size) {
    at <- unique(at) %||% scales::breaks_extended()(matrix)
    if (!is.null(radius_range)) {
        radius <- scale_range(unclass(radius_range), matrix)(at)
    } else if (!is.null(radius)) {
        breaks_len <- length(at)
        radius_len <- length(radius)
        if (radius_len != breaks_len) {
            cli::cli_abort(c(
                "{.arg radius} and {.arg at} must have the same length",
                i = "the length of {.arg at}: {.val {breaks_len}}",
                i = "the length of {.arg radius}: {.val {radius_len}}"
            ))
        }
        # unit function will not omit the existed units
        radius <- grid::unit(radius, "mm")
    } else {
        radius <- scale_size(at)
    }
    graphics <- vector("list", length = length(radius))
    # lapply will omit unit attributes, so we just use for loop
    for (i in seq_along(graphics)) {
        graphics[[i]] <- circle_legend(radius[[i]])
    }
    list(at = at, graphics = graphics)
}

#' @return A scale function taking a numeric vector and return a scaled unit
#'   values.
#' @noRd
scale_range <- function(range, matrix) {
    force(range)
    from <- range(matrix, na.rm = TRUE, finite = TRUE)
    function(x) {
        x <- scales::rescale(x, to = range, from = from)
        grid::unit(x, "mm")
    }
}

circle_legend <- function(size) {
    force(size)
    function(x, y, w, h) {
        grid::grid.circle(x, y,
            r = size,
            gp = gpar(fill = "black", col = NA)
        )
    }
}

#' @importFrom grid gpar
#' @export
grid::gpar

#' @importClassesFrom ComplexHeatmap Heatmap
#' @importClassesFrom ComplexHeatmap Legends
methods::setClass(
    "DotsHeatmap",
    slots = list(
        heatmap = "Heatmap",
        dots_legend = "Legends"
    )
)

#' @importFrom ComplexHeatmap draw
#' @export
#' @noRd
ComplexHeatmap::draw

#' Draw a Dots Heatmap
#' @description
#'  These objects are imported from other packages. Follow the links below to
#'  see their documentation.
#' \describe{
#'   \item{draw}{\code{\link[ComplexHeatmap]{draw}}}
#' }
#' @param object A [DotsHeatmap] object.
#' @param ... Arguments passed to
#'   [Heatmap][ComplexHeatmap::draw,HeatmapList-method].
#' @export
#' @rdname draw
methods::setMethod("draw", signature = "DotsHeatmap", function(object, ...) {
    draw(object@heatmap, heatmap_legend_list = object@dots_legend, ...)
})

#' Draw a DotsHeatmap with default parameters
#' @details Actually it calls [draw], but only with default parameters. If users
#'   want to customize the DotsHeatmap, they can pass parameters directly to
#'   [draw].
#' @param object A [DotsHeatmap] object.
#' @importFrom methods show
#' @export
methods::setMethod("show", signature = "DotsHeatmap", function(object) {
    draw(object)
})
