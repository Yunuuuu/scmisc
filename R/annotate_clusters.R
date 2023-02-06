#' Annotate cluster cell types with provided markers
#'
#' @description This function calculate the average expression levels of all
#'   markers in a group of cells. Then annotate this cluster with the cell type
#'   whose markers have the maximal average expression levels in this cluster.
#'   Finally, the function will also adjust to the results according to manual
#'   if it is not `NULL`.
#' @param x A numeric matrix of counts with cells in columns and features in
#' rows.
#'
#' Alternatively, a
#' [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment] or
#' [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment] object
#' containing such a matrix.
#' @param clusters A factor (or vector coercible into a factor) specifying the
#'   group to which each cell in `x` belongs. Alternatively, String specifying
#'   the field of `colData(x)` containing the grouping factor if `x` is a
#'   [SingleCellExperiment][SingleCellExperiment::SingleCellExperiment], e.g.,
#'   usually clusters.
#' @param marker_list A named list contaning the markers for each cell types.
#'   names indicate the cell type label and values indicate the markers for this
#'   cell type.
#' @param manual A named list contaning the manually annotated cell types for
#'   the clusters. The names indicate the cell type label and the values
#'   indicate the cluster labels (will be coerced to character).
#' @param ... For the generic, further arguments to pass to specific methods.
#' @return A named character, names indicates the cluster labels and values
#'   indicates the cell type names.
#' @name annotate_clusters
NULL

annotate_clusters_internal <- function(x, clusters, marker_list, manual = NULL) {
    assert_class(marker_list, is.list, "list", null_ok = FALSE)
    assert_class(manual, is.list, "list", null_ok = TRUE)
    if (length(marker_list) > 0L) {
        if (all(!has_names(marker_list))) {
            cli::cli_abort("All elements in {.arg marker_list} must be named.")
        }
    }
    if (length(manual) > 0L) {
        if (all(!has_names(manual))) {
            cli::cli_abort("All elements in {.arg manual} must be named.")
        }
    }
    cluster2cell <- lapply(marker_list, function(markers) {
        markers <- unique(markers)
        # we firstly calculate the sum expression values of all markers in a
        # each cluster.
        sum_array <- scuttle::summarizeAssayByGroup(
            x,
            subset.row = markers,
            ids = clusters,
            statistics = "sum"
        )
        # a named numeric vector
        cluster_sums <- colSums(
            SummarizedExperiment::assay(sum_array, "sum"),
            na.rm = TRUE
        )
        i <- names(cluster_sums)

        # we then get the number of items: cells in this cluster * number of
        # markers
        # use `c` function to removing attributes except names
        numbers <- c(table(clusters))
        numbers <- numbers[i] * length(markers)
        data.table::data.table(
            clusters = i,
            means = cluster_sums / numbers
        )
    })
    # we annotate the cluster as the celltype whose
    cluster2cell <- data.table::rbindlist(
        cluster2cell,
        use.names = TRUE,
        idcol = "celltype"
    )[, list(celltype = celltype[[which.max(means)]]), by = "clusters"][ # nolint
        , structure(celltype, names = clusters) # nolint
    ]

    # adjust results by manual.
    if (length(manual) == 0L) {
        cluster2cell_levels <- intersect(names(marker_list), cluster2cell)
    } else {
        cluster2cell[as.character(
            unlist(manual, recursive = FALSE, use.names = FALSE)
        )] <- rep(names(manual), times = lengths(manual))
        cluster2cell_levels <- c(
            intersect(names(marker_list), cluster2cell),
            setdiff(names(manual), names(marker_list))
        )
    }
    factor(cluster2cell, levels = cluster2cell_levels)
}

utils::globalVariables(c("means", "celltype"))

#' @export
#' @rdname annotate_clusters
setGeneric(
    "annotate_clusters",
    function(x, ...) standardGeneric("annotate_clusters")
)

#' @export
#' @rdname annotate_clusters
setMethod("annotate_clusters", "ANY", annotate_clusters_internal)

#' @export
#' @param assay.type A string or integer scalar indicating which `assays` in the
#' `x` contains the count matrix.
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @rdname annotate_clusters
setMethod(
    "annotate_clusters", "SingleCellExperiment",
    function(x, clusters = NULL, ..., assay.type = "logcounts") {
        if (is.null(clusters)) {
            clusters <- x$label
            if (is.null(clusters)) {
                cli::cli_abort("{.field label} never exist in {.arg x}.")
            }
        } else if (rlang::is_scalar_character(clusters) && ncol(x) > 1L) {
            clusters <- scater::retrieveCellInfo(
                x, clusters,
                search = "colData"
            )$value
        }
        annotate_clusters_internal(
            SummarizedExperiment::assay(x, assay.type),
            clusters = clusters,
            ...
        )
    }
)

#' @export
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @rdname annotate_clusters
setMethod(
    "annotate_clusters", "SummarizedExperiment",
    function(x, clusters = NULL, ..., assay.type = "logcounts") {
        if (is.null(clusters)) {
            clusters <- x$label
            if (is.null(clusters)) {
                cli::cli_abort("{.field label} never exist in {.arg x}.")
            }
        }
        annotate_clusters_internal(
            SummarizedExperiment::assay(x, assay.type),
            clusters = clusters,
            ...
        )
    }
)

#' @export
#' @importClassesFrom SeuratObject Seurat
#' @rdname annotate_clusters
setMethod(
    "annotate_clusters", "Seurat",
    function(x, clusters = NULL, ..., assay.use = NULL) { # nolint
        assay.use <- assay.use %||% Seurat::DefaultAssay(x) 
        if (is.null(clusters)) {
            clusters <- SeuratObject::Idents(x)
            if (is.null(clusters)) {
                cli::cli_abort("{.field Idents} never exist in {.arg x}.")
            }
        } else if (rlang::is_scalar_character(clusters) && ncol(x) > 1L) {
            x <- SeuratObject::SetIdent(x, value = clusters)
            clusters <- SeuratObject::Idents(x)
        }
        annotate_clusters_internal(
            SeuratObject::GetAssayData(
                object = x, slot = "data",
                assay = assay.use
            ),
            clusters = clusters,
            ...
        )
    }
)
