#' Fast mutual nearest neighbors correction
#' @description
#' Correct for batch effects in single-cell expression data using a fast version
#' of the mutual nearest neighbors (MNN) method. This is similar with
#' [fastMNN][batchelor::fastMNN], but use
#' [multiBatchNorm][batchelor::multiBatchNorm] to normalize differences in
#' sequencing depth rather than [cosineNorm][batchelor::cosineNorm].
#' @param x A single
#'   [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] object.
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
#' @inheritParams batchelor::multiBatchPCA
#' @inheritParams batchelor::reducedMNN
#' @details
#' This function provides a variant of the [fastMNN][batchelor::fastMNN]
#' function, modified for simplicity
#' In particular:
#' \itemize{
#' \item It performs a batch normalization via
#' [multiBatchNorm][batchelor::multiBatchNorm].
#'
#' \item It performs a multi-sample PCA via
#' [multiBatchPCA][batchelor::multiBatchPCA] and subsequently performs all
#' calculations in the PC space.  This reduces computational work and provides
#' some denoising for improved neighbour detection.  As a result, though, the
#' corrected output cannot be interpreted on a gene level and is useful only for
#' cell-level comparisons, e.g., clustering and visualization.
#'
#' \item The correction vector for each cell is directly computed from its
#' \code{k} nearest neighbours in the same batch.  Specifically, only the
#' \code{k} nearest neighbouring cells that \emph{also} participate in MNN pairs
#' are used.  Each MNN-participating neighbour is weighted by distance from the
#' current cell, using a tricube scheme with bandwidth equal to the median
#' distance multiplied by \code{ndist}.  This ensures that the correction vector
#' only uses information from the closest cells, improving the fidelity of local
#' correction.
#'
#' \item Issues with \dQuote{kissing} are avoided with a two-step procedure that
#' removes variation along the batch effect vector.  First, the average
#' correction vector across all MNN pairs is computed.  Cell coordinates are
#' adjusted such that all cells in a single batch have the same position along
#' this vector.  The correction vectors are then recalculated with the adjusted
#' coordinates (but the same MNN pairs).
#' }
#'
#' The \code{batch} argument allows users to easily perform batch correction
#' when all cells have already been combined into a single object.  This avoids
#' the need to manually split the matrix or SingleCellExperiment object into
#' separate objects for input into \code{fastMNN}.  In this situation, the order
#' of input batches is defined by the order of levels in \code{batch}.
#'
#' @return A `SingleCellExperiment` objects with normalized log-expression
#' values in the "multiBatchNorm" assay (depending on values in norm.args). A
#' "PCA" matrix in the `reducedDims` slot, containing the first d PCs returned
#' by `multiBatchPCA` for each cell. A "corrected" matrix in the `reducedDims`
#' slot, containing corrected low-dimensional coordinates for each cell.
#' @seealso 
#' - [multiBatchNorm][batchelor::multiBatchNorm]
#' - [multiBatchPCA][batchelor::multiBatchPCA]
#' - [reducedMNN][batchelor::reducedMNN]
#' @examples
#' d1 <- matrix(rnbinom(50000, mu = 10, size = 1), ncol = 100)
#' d2 <- matrix(rnbinom(20000, mu = 50, size = 1), ncol = 40)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'     list(counts = cbind(d1, d2))
#' )
#' sce$batch <- rep(c("a", "b"), c(100L, 40L))
#' out <- fast_mnn(sce, batch = sce$batch)
#' # Corrected values for use in clustering, etc.
#' str(SingleCellExperiment::reducedDim(out))
#' # Extracting corrected expression values for gene 10.
#' summary(SummarizedExperiment::assay(out)[10, ])
#' @export
fast_mnn <- function(
    x,
    batch = NULL,
    norm.args = list(),
    assay.type = "counts",
    subset.row = NULL,
    d = 50L, k = 20L, prop.k = NULL,
    min.mean = 1L,
    weights = NULL,
    restrict = NULL, ndist = 3L,
    merge.order = NULL, auto.merge = FALSE,
    min.batch.skip = 0L,
    deferred = TRUE,
    BSPARAM = BiocSingular::IrlbaParam(),
    BNPARAM = BiocNeighbors::AnnoyParam(),
    BPPARAM = BiocParallel::SerialParam()) {
    assert_class(
        x,
        function(x) methods::is(x, "SingleCellExperiment"), "SingleCellExperiment"
    )
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
    if (is.null(norm.args$name)) {
        norm.args$name <- "multiBatchNorm"
    }
    x <- batchelor::multiBatchNorm(
        x,
        batch = batch,
        norm.args = norm.args,
        min.mean = min.mean,
        subset.row = NULL,
        normalize.all = TRUE,
        preserve.single = TRUE,
        assay.type = assay.type,
        BPPARAM = BPPARAM
    )
    batch_pcs <- batchelor::multiBatchPCA(x,
        batch = batch, d = d,
        subset.row = subset.row,
        weights = weights,
        get.all.genes = FALSE,
        get.variance = FALSE,
        preserve.single = TRUE,
        assay.type = norm.args$name,
        BSPARAM = BSPARAM,
        deferred = deferred,
        BPPARAM = BPPARAM
    )
    SingleCellExperiment::reducedDim(x, "PCA") <- batch_pcs[[1L]]

    mnn_res <- batchelor::reducedMNN(
        batch_pcs[[1L]],
        batch = batch,
        k = k, prop.k = prop.k,
        restrict = restrict, ndist = ndist,
        merge.order = merge.order, auto.merge = auto.merge,
        min.batch.skip = min.batch.skip,
        BNPARAM = BNPARAM, BPPARAM = BPPARAM
    )

    ## keep reducedDim data -------------------------------------------
    SingleCellExperiment::reducedDim(x, "corrected") <- mnn_res$corrected

    x
}
