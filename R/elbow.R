#' Find the elbow point in the curve of variance explained by each successive PC
#'
#' Find the elbow point in the curve of variance explained by each successive
#' PC. This can be used to determine the number of PCs to retain.
#'
#' @param x A Numeric vector containing the variance explained by each PC.
#' @param dimred A string or integer scalar indicating the reduced dimension
#'   containing the `PCA` with a `percentVar` attribute.
#' @param ... Other arguments passed to specific methods.
#' @return An integer scalar specifying the number of PCs at the elbow point.
#' @export
#' @name pca_elbow
setGeneric(
    "pca_elbow",
    function(x, ...) standardGeneric("pca_elbow")
)

#' @rdname pca_elbow
#' @export
setMethod(
    "pca_elbow", "numeric",
    function(x) {
        # Percentage of variance explained is tucked away in the attributes.
        PCAtools::findElbowPoint(sort(x, decreasing = TRUE))
    }
)

#' @rdname pca_elbow
#' @export
setMethod(
    "pca_elbow", "SingleCellExperiment",
    function(x, dimred = "PCA") {
        var <- attr(SingleCellExperiment::reducedDim(x, dimred), "percentVar")
        pca_elbow(var)
    }
)
