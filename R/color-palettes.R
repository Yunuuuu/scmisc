#' SCANPY Color Palettes
#'
#' Color palette inspired by scanpy in
#' <https://github.com/scverse/scanpy/blob/master/scanpy/plotting/palettes.py>.
#'
#' @param palette Palette type.
#'   Currently there is two available option: `"zeileis"` (28-color palette) and
#'   `"godsnot"` (102-color palette).
#' @param alpha Transparency level, a real number in `(0, 1]`.
#'   See `alpha` in [rgb][grDevices::rgb] for details.
#' @export
pal_scanpy <- function(palette = c("zeileis", "godsnot"), alpha = 1) {
    palette <- match.arg(palette)

    if (alpha > 1L || alpha <= 0L) {
        cli::cli_abort("alpha must be in {.filed (0, 1]}")
    }
    raw_cols <- scanpy_pal[[palette]]
    raw_cols_rgb <- grDevices::col2rgb(raw_cols)
    alpha_cols <- grDevices::rgb(
        raw_cols_rgb[1L, , drop = TRUE], 
        raw_cols_rgb[2L, , drop = TRUE], 
        raw_cols_rgb[3L, , drop = TRUE],
        alpha = alpha * 255L,
        names = names(raw_cols),
        maxColorValue = 255L
    )
    scales::manual_pal(unname(alpha_cols))
}

#' SCANPY Color Scales
#'
#' See \code{\link{pal_scanpy}} for details.
#'
#' @inheritParams pal_scanpy
#' @param ... additional parameters for
#' [discrete_scale][ggplot2::discrete_scale]
#'
#' @rdname scale_scanpy
#' @export 
scale_color_scanpy <- function(palette = c("zeileis", "godsnot"), alpha = 1, ...) {
    palette <- match.arg(palette)
    ggplot2::discrete_scale(
        "colour", "scanpy", pal_scanpy(palette, alpha), ...
    )
}

#' @rdname scale_scanpy
#' @export 
scale_colour_scanpy <- scale_color_scanpy

#' @rdname scale_scanpy
#' @export 
scale_fill_scanpy <- function(palette = c("zeileis", "godsnot"), alpha = 1, ...) {
    palette <- match.arg(palette)
    ggplot2::discrete_scale(
        "fill", "scanpy", pal_scanpy(palette, alpha), ...
    )
}

scanpy_pal <- list(
    godsnot = c(
        "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941",
        "#006FA6", "#A30059", "#FFDBE5", "#7A4900", "#0000A6",
        "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
        "#5A0007", "#809693", "#6A3A4C", "#1B4400", "#4FC601",
        "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900",
        "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA",
        "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
        "#300018", "#0AA6D8", "#013349", "#00846F", "#372101",
        "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2",
        "#C2FF99", "#001E09", "#00489C", "#6F0062", "#0CBD66",
        "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
        "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459",
        "#456648", "#0086ED", "#886F4C", "#34362D", "#B4A8BD",
        "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F",
        "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF",
        "#3B9700", "#04F757", "#C8A1A1", "#1E6E00", "#7900D7",
        "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600",
        "#D790FF", "#9B9700", "#549E79", "#FFF69F", "#201625",
        "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
        "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98",
        "#A4E804", "#324E72"
    ),
    zeileis = c(
        "#023fa5", "#7d87b9", "#bec1d4", "#d6bcc0", "#bb7784",
        "#8e063b", "#4a6fe3", "#8595e1", "#b5bbe3", "#e6afb9",
        "#e07b91", "#d33f6a", "#11c638", "#8dd593", "#c6dec7",
        "#ead3c6", "#f0b98d", "#ef9708", "#0fcfc0", "#9cded6",
        "#d5eae7", "#f3e1eb", "#f6c4e1", "#f79cd4", "#7f7f7f",
        "#c7c7c7", "#1CE6FF", "#336600"
    )
)
