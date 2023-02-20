#' Fast mutual nearest neighbors correction
#' @description 
#' Correct for batch effects in single-cell expression data using a fast version
#' of the mutual nearest neighbors (MNN) method. This is modified from
#' [fastMNN][batchelor::fastMNN], but use
#' [multiBatchNorm][batchelor::multiBatchNorm] to normalize differences in
#' sequencing depth rather than [cosineNorm][batchelor::cosineNorm]. `cos.norm`
#' argument of [fastMNN][batchelor::fastMNN] is always set to `FALSE`.
#' @param subset.row A vector specifying which features used to run
#'   [multiBatchPCA][batchelor::multiBatchPCA]. 
#' @param BSPARAM A [BiocSingularParam][BiocSingular::IrlbaParam] object
#'   specifying the algorithm to use for PCA in
#'   [multiBatchPCA][batchelor::multiBatchPCA].
#' @param BNPARAM A [BiocNeighborParam][BiocNeighbors::BiocNeighborParam] object
#'   specifying the nearest neighbor algorithm.
#' @param BPPARAM A [BiocParallelParam][BiocParallel::BiocParallelParam] object
#'   specifying whether the PCA and nearest-neighbor searches should be
#'   parallelized.
#' @inheritParams batchelor::multiBatchNorm
#' @inherit batchelor::fastMNN
#' @examples 
#'  d1 <- matrix(rnbinom(50000, mu=10, size=1), ncol=100)
#'  sce1 <- SingleCellExperiment::SingleCellExperiment(list(counts=d1))
#'  d2 <- matrix(rnbinom(20000, mu=50, size=1), ncol=40)
#'  sce2 <- SingleCellExperiment::SingleCellExperiment(list(counts=d2))
#'  out <- fast_mnn(sce1, sce2)
#'  # Corrected values for use in clustering, etc.
#'  str(SingleCellExperiment::reducedDim(out)) 
#'  # Extracting corrected expression values for gene 10.
#'  summary(SummarizedExperiment::assay(out)[10,])
#' @importClassesFrom S4Vectors List
#' @importClassesFrom S4Vectors DataFrame
#' @export 
fast_mnn <- function(
    ...,
    batch = NULL,
    norm.args = list(),
    min.mean = 1L,
    assay.type = "counts",
    k = 20L, prop.k = NULL,
    restrict = NULL, ndist = 3,
    d = 50L, deferred = TRUE,
    weights = NULL, get.variance = FALSE,
    merge.order = NULL, auto.merge = FALSE,
    min.batch.skip = 0L, subset.row = NULL,
    correct.all = FALSE,
    BSPARAM = BiocSingular::IrlbaParam(),
    BNPARAM = BiocNeighbors::AnnoyParam(),
    BPPARAM = BiocParallel::SerialParam()) {

    if (is.null(norm.args$name)) {
        norm.args$name <- "multiBatchNorm"
    }
    # nromalization, adjust for differences in sequencing depth --------
    # This can be done before variance-modelling or after.
    # - http://bioconductor.org/books/3.15/OSCA.multisample/integrating-datasets.html#slower-setup;
    #   in this book, variance-modelling is performed before `multiBatchNorm`
    # - https://github.com/LTLA/batchelor/blob/master/R/quickCorrect.R; in the
    #   source code of `quickCorrect`, variance-modelling is performed after
    #   multiBatchNorm.
    # 
    # so both method should be okay. but I prefered to use multiBatchNorm for
    # the batch corrected logcounts. As batchelor::quickCorrect use this firstly
    # and then run modelGeneVar on the returned `logcounts` assay, we just do
    # this firstly. The class of the returned `logcounts` will determined by the
    # class of "counts" assay.  If the computer has a large memory, it's better
    # to use `dgCMatrix` since it's much faster and accurater.
    # 
    # - multiBatchNorm need size factor, if it's NULL, the internal will
    #   calculate it with `librarySizeFactors` function
    # - normalizetion shouldn't use the `subset.row` since library size is the
    #   total sum of counts across all genes for each cell
    batch_list <- batchelor::multiBatchNorm(
        ...,
        batch = batch,
        norm.args = norm.args,
        min.mean = min.mean,
        subset.row = NULL,
        normalize.all = TRUE,
        preserve.single = FALSE,
        assay.type = assay.type,
        BPPARAM = BPPARAM
    )
    batchelor::fastMNN(
        batch_list,
        batch = NULL, k = k, prop.k = prop.k,
        restrict = restrict, cos.norm = FALSE, ndist = ndist,
        d = d, deferred = deferred,
        weights = weights, get.variance = get.variance,
        merge.order = merge.order, auto.merge = auto.merge,
        min.batch.skip = min.batch.skip, subset.row = subset.row,
        correct.all = correct.all, 
        assay.type = norm.args$name,
        BSPARAM = BSPARAM,
        BNPARAM = BNPARAM,
        BPPARAM = BPPARAM
    )
}
