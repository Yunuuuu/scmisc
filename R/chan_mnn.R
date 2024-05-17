#' Fast mutual nearest neighbors correction with `scran.chan`
#'
#' @param x A single
#'   [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] object.
#' @inheritParams scran.chan::logNormCounts.chan
#' @param assay.type A string specifying which assay values contains the counts.
#' @param subset.row Integer, logical or character vector specifying which
#' features to use in the PCA (e.g., highly variable genes). If NULL, all
#' features in x are used.
#' @param d Integer scalar specifying the number of top PCs to obtain.
#' @param k Integer scalar specifying the number of neighbors to use when
#' identifying MNN pairs.
#' @inheritParams scran.chan::runPCA.chan
#' @inheritDotParams scran.chan::mnnCorrect.chan -x -batch -k -approximate -num.threads
#' @param approximate Logical scalar specifying whether to perform an
#' approximate neighbor search.
#' @param threads Integer scalar specifying the number of threads to use.
#' @seealso
#' - [logNormCounts.chan][scran.chan::logNormCounts.chan]
#' - [runPCA.chan][scran.chan::runPCA.chan]
#' - [mnnCorrect.chan][scran.chan::mnnCorrect.chan]
#' @export
chan_mnn <- function(
    x, batch = NULL,
    assay.type = "counts",
    subset.row = NULL,
    d = 50L, k = 20L,
    scale = FALSE,
    batch.mode = NULL,
    batch.method = NULL,
    ...,
    approximate = TRUE,
    threads = 1L) {
    assert_pkg("scran.chan")
    assert_s4_class(x, "SingleCellExperiment")
    assert_bool(approximate)
    threads <- as.integer(threads)

    # nromalization, adjust for differences in sequencing depth --------
    assay <- SummarizedExperiment::assay(x, assay.type)
    norm <- scran.chan::logNormCounts.chan(
        x = scran.chan::initializeSparseMatrix(assay),
        size.factors = SingleCellExperiment::sizeFactors(x),
        batch = batch, batch.mode = batch.mode,
        num.threads = threads
    )

    # dimensionality reduction -----------------------------------------
    batch_pcs <- scran.chan::runPCA.chan(
        x = norm,
        num.comp = as.integer(d),
        subset = subset.row,
        scale = scale,
        num.threads = threads,
        batch = batch,
        batch.method = batch.method,
        rotation = TRUE
    )

    # run MNN --------------------------------------------------------
    mnn_res <- scran.chan::mnnCorrect.chan(
        x = batch_pcs$components,
        batch = batch, k = as.integer(k), ...,
        approximate = approximate,
        num.threads = threads
    )

    # Saving the results
    # SingleCellExperiment::sizeFactors(x) <- norm$size.factors
    # SummarizedExperiment::assay(x, "ChanBatchNorm") <-
    #     log1p(t(t(DelayedArray::DelayedArray(assay)) / norm$size.factors)) /
    #         log(2L)
    pca <- t(batch_pcs$components)
    attr(pca, "percentVar") <- batch_pcs$prop.variance
    attr(pca, "rotation") <- batch_pcs$rotation
    SingleCellExperiment::reducedDim(x, "PCA") <- pca
    corrected <- t(mnn_res$corrected)
    attr(corrected, "merge.order") <- mnn_res$merge.order
    attr(corrected, "num.pairs") <- mnn_res$num.pairs
    SingleCellExperiment::reducedDim(x, "corrected") <- corrected
    x
}
