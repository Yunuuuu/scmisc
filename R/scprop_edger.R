#' Changes in cluster abundance
#'
#' @param data A data.frame.
#' @param sample A string of the column in `data` containing sample id.
#' @param celltype A string of the column in `data` containing cell types.
#' @param formula A formula object to build design matrix from `data`.
#' @param filter A bool, if `TRUE`, will use [filterByExpr][edgeR::filterByExpr]
#' to filter.
#' @param normalize A bool, if `TRUE`, will use
#' [calcNormFactors][edgeR::calcNormFactors] to compute normalization factors.
#' @param coef An integer or character index vector indicating which
#' coefficients of the linear model are to be tested equal to zero. If `NULL`,
#' the last coef will be used.
#' @param p.adjust Correction method, a string. Can be abbreviated. See
#' [p.adjust][stats::p.adjust].
#' @param ... For the generic, further arguments to pass to specific methods.
#' @return A [data.table][data.table::data.table].
#' @references 
#' <https://bioconductor.org/books/3.18/OSCA.multisample/differential-abundance.html#>
#' @name scprop_edger
NULL

scprop_edger_internal <- function(data, sample, celltype, formula, filter = TRUE, normalize = FALSE, coef = NULL, p.adjust = "BH") {
    assert_pkg("edgeR")
    assert_string(sample)
    assert_string(celltype)
    assert_data_frame_columns(data, c(sample, celltype))
    assert_(formula, rlang::is_formula, "a {.cls formula}")
    data <- data.table::as.data.table(data)
    abundances <- unclass(table(data[[celltype]], data[[sample]]))
    idx <- match(colnames(abundances), data[[sample]])
    data <- data[idx]
    data.table::setDF(data)
    ## edgeR ----------------------------------------------------
    edger_data <- edgeR::DGEList(abundances, samples = data)
    design <- stats::model.matrix(formula, data = edger_data$samples)
    if (isTRUE(filter)) {
        keep <- edgeR::filterByExpr(edger_data, design = design)
        edger_data <- edger_data[keep, ]
    }
    if (isTRUE(normalize)) {
        edger_data <- edgeR::calcNormFactors(edger_data)
    }
    edger_data <- edgeR::estimateDisp(edger_data, design, trend.method = "none")
    fit_res <- edgeR::glmQLFit(edger_data, design,
        robust = TRUE, abundance.trend = FALSE
    )
    res <- edgeR::glmQLFTest(fit_res, coef = coef %||% ncol(design))
    out <- edgeR::topTags(res, n = Inf, adjust.method = p.adjust)
    data.table::as.data.table(as.data.frame(out), keep.rownames = "cells")
}

#' @export
#' @rdname scprop_edger
setGeneric(
    "scprop_edger",
    function(data, ...) standardGeneric("scprop_edger")
)

#' @export
#' @rdname scprop_edger
setMethod("scprop_edger", "data.frame", scprop_edger_internal)

#' @export
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @rdname scprop_edger
setMethod("scprop_edger", "SummarizedExperiment", function(data, ...) {
    scprop_edger_internal(
        data = as.data.frame(
            SummarizedExperiment::colData(data),
            check.names = FALSE, make.names = FALSE
        ),
        ...
    )
})
