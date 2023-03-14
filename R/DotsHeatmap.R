#' Constructor method for DotsHeatmap class
#' @param matrix A matrix. Either numeric or character. If it is a simple
#'   vector, it will be converted to a one-column matrix. Only used for dot
#'   color.
#' @param matrix_size A numeric matrix or vector. If it is a simple vector, it
#'   will be converted to a one-column matrix. Only used for dot size.
#' @param dots_size A numeric vector of length one or two specifying the the
#'   minimal and maximal relative dots radius compared with the heatmap
#'   `min(grid::unit.c(width, height)) / 2L`. If length one, all dots will have
#'   the same size. If length two with different values, the `matrix_size` will
#'   be scaled into `dots_size`. Default value: `c(0.1, 1)`, means the maximal
#'   value in matrix_size will have a radius equal to the heatmap
#'   `min(grid::unit.c(width, height)) / 2L` and the minimal value will equal to
#'   `0.1 * min(grid::unit.c(width, height)) / 2L`.
#' @param dots_size_legend_param Other arguments passed to
#'   [Legend][ComplexHeatmap::Legend]. This will be used to define the dots size
#'   legend.
#' @param slice_border_gp Graphic parameters for drawing slice rectangles. The
#'   value should be specified by [gpar][grid::gpar] and fill parameter is
#'   always setted to "transparent". If `NULL`, no slice border will be drawn.
#' @inheritDotParams ComplexHeatmap::Heatmap -matrix -rect_gp -layer_fun
#' @return A [DotsHeatmap] object
#' @name DotsHeatmap
#' @export
DotsHeatmap <- function(matrix, matrix_size = NULL, dots_size = c(0.1, 1), dots_size_legend_param = list(), slice_border_gp = NULL, ...) {
    assert_class(dots_size, is.numeric, "numeric")
    if (!data.table::between(length(dots_size), 1L, 2L)) {
        cli::cli_abort("{.arg dots_size} must be a numeric with a length 1 or 2")
    }
    if (anyNA(dots_size) || any(dots_size < 0L)) {
        cli::cli_abort("{.arg dots_size} must not be negative or {.val NA}")
    }
    dots_size <- rep_len(dots_size, 2L)
    matrix_size <- matrix_size %||% matrix
    assert_class(matrix_size, is.numeric, "numeric")
    if (is.atomic(matrix_size)) {
        cli::cli_alert_info("convert simple vector to one-column matrix")
        matrix_size <- matrix(matrix_size, ncol = 1L)
    } else if (!is.matrix(matrix)) {
        cli::cli_abort("{.arg matrix_size} must be a matrix or a simple vector.")
    }

    if (!all(dim(matrix) == dim(matrix_size))) {
        cli::cli_abort("{.arg matrix} and {.arg matrix_size} must have the same dimensions")
    }

    scale_size <- radius_fn(dots_size, matrix_size)
    dots_size_legend_param$at <- dots_size_legend_param$at %||%
        scales::breaks_extended()(matrix_size)
    dots_size_legend_param$legend_gp <- dots_size_legend_param$legend_gp %||%
        gpar()
    if (is.null(dots_size_legend_param$legend_gp$size)) {
        dots_size_legend_param$legend_gp$size <- scale_size(
            dots_size_legend_param$at
        )
    }
    dots_size_legend_param$type <- "point"
    methods::new(
        "DotsHeatmap",
        heatmap = ComplexHeatmap::Heatmap(
            matrix,
            rect_gp = gpar(type = "none"),
            layer_fun = function(j, i, x, y, width, height, fill) {
                size_values <- ComplexHeatmap::pindex(matrix_size, i = i, j = j)
                grid::grid.circle(
                    x = x, y = y,
                    r = scale_size(size_values),
                    gp = gpar(fill = fill, col = NA)
                )
                if (!is.null(slice_border_gp)) {
                    grid::grid.rect(gp = slice_border_gp)
                }
            },
            ...
        ),
        dots_legend = do.call(ComplexHeatmap::Legend, dots_size_legend_param)
    )
}

radius_fn <- function(dots_size, matrix) {
    force(dots_size)
    min_unit <- min(grid::unit.c(
        grid::unit(1 / ncol(matrix), "npc"),
        grid::unit(1 / nrow(matrix), "npc")
    ))
    from <- range(matrix, na.rm = TRUE, finite = TRUE)
    function(x) {
        scales::rescale(x, to = dots_size, from = from) * min_unit / 2
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
