#' Quick and dirty subclustering
#'
#' Performs a quick subclustering for all cells within each group.
#'
#' @param ... Other arguments passed to
#' [quickSubCluster][scran::quickSubCluster] 
#'   - `x` A matrix of counts or log-normalized expression values (if
#'     \code{normalize=FALSE}), where each row corresponds to a gene and each
#'     column corresponds to a cell. Alternatively, a
#'     [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment] or
#'     [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] object
#'     containing such a matrix. 
#'   - `groups`: A vector of group assignments for all cells, usually
#'     corresponding to cluster identities.
#'   - `normalize`: Logical scalar indicating whether each subset of \code{x}
#'     should be log-transformed prior to further analysis.
#'   - `restricted`: Character vector containing the subset of groups in
#'     \code{groups} to be subclustered. By default, all unique groups in
#'     \code{groups} are used for subclustering, but this can be restricted to
#'     specific groups of interest to save compute time.
#'   - `prepFUN`: A function that accepts a single
#'     [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] object
#'     and returns another
#'     [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment]
#'     containing any additional elements required for clustering (e.g., PCA
#'     results).
#'   - `min.ncells`: An integer scalar specifying the minimum number of cells in
#'     a group to be considered for subclustering.
#'   - clusterFUN A function that accepts a single
#'     [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] object
#'     and returns a vector of cluster assignments for each cell in that object.
#'   - `BLUSPARAM`: A [BlusterParam][bluster::BlusterParam-class] object that is
#'     used to specify the clustering via [clusterRows][bluster::clusterRows].
#'     Only used when \code{clusterFUN=NULL}.
#'   - `assay.type`: String or integer scalar specifying the relevant assay.
#'
#' @return A factor
#' @seealso [quickSubCluster][scran::quickSubCluster]
sub_cluster <- function(...) {
    assert_pkg("scran")
    labels <- scran::quickSubCluster(..., format = "%s.%s", simplify = TRUE)
    order_labels(labels)
}

order_labels <- function(x) {
    lvls <- unique(x)
    main_lvls <- as.integer(lvls)
    int_lvls <- main_lvls[main_lvls == lvls] # main cluster labels
    dbl_lvls <- as.numeric(setdiff(lvls, int_lvls)) # sub-clustering labels
    dbl_lvls_list <- split(dbl_lvls, as.integer(dbl_lvls))
    mapping <- structure(int_lvls, names = int_lvls)
    max_lvls <- max(mapping)
    for (i in seq_along(dbl_lvls_list)) {
        cur_lvls <- sort(dbl_lvls_list[[i]])
        mapping <- c(mapping, structure(
            c(cur_lvls[1L], seq_len(length(cur_lvls) - 1L) + max_lvls),
            names = cur_lvls
        ))
        max_lvls <- max(mapping)
    }
    factor(x, levels = names(mapping), labels = mapping)
}
