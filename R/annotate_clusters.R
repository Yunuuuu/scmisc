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
#'   group to which each cell in `x` belongs. Alternatively, if `x` is a
#'   [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment], e.g,
#'   String specifying the field of `colData(x)` containing the grouping factor.
#'   In this way, if clusters is `NULL`, "label" in `colData(x)` will be
#'   extracted.
#' @param marker_list A named list contaning the markers for each cell types.
#'   names indicate the cell type label and values indicate the markers for this
#'   cell type.
#' @param manual A named list contaning the manually annotated cell types for
#'   the clusters. The names indicate the cell type label and the values
#'   indicate the cluster labels (will be coerced to character).
#' @param blocks A factor (or vector coercible into a factor) specifying the
#'   blocking factor to which each cell in `x` belongs (e.g., batch of origin).
#'   Alternatively, if `x` is a
#'   [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment], e.g,
#'   String specifying the field of `colData(x)` containing the grouping factor.
#'   In this way, if clusters is `NULL`, "label" in `colData(x)` will be
#'   extracted.
#' @param ... For the generic, further arguments to pass to specific methods.
#' @return A named character, names indicates the cluster labels and values
#'   indicates the cell type names.
#' @name annotate_clusters
NULL

annotate_clusters_internal <- function(x, clusters, marker_list, manual = NULL, blocks = NULL) {
    assert_class(marker_list, "list", is.list, null_ok = FALSE)
    assert_class(manual, "list", is.list, null_ok = TRUE)
    if (length(marker_list) == 0L) {
        cli::cli_abort("Empty list is not allowed in {.arg marker_list}.")
    } else if (!rlang::is_named(marker_list)) {
        cli::cli_abort("All elements in {.arg marker_list} must be named.")
    }

    # prepare manual annotation
    if (length(manual) > 0L) {
        if (!rlang::is_named(manual)) {
            cli::cli_abort("All elements in {.arg manual} must be named.")
        }
        manual_clusters <- as.character(
            unlist(manual, recursive = FALSE, use.names = FALSE)
        )
        manual_labels <- rep(names(manual), times = lengths(manual))
        if (anyNA(manual_clusters)) {
            cli::cli_warn(c(
                "{.val {NA}} is found in {.arg manual}",
                "i" = "will omit {.val {NA}}"
            ))
            manual_clusters <- manual_clusters[!is.na(manual_clusters)]
            manual_labels <- manual_labels[!is.na(manual_clusters)]
        }
        missed_clusters <- setdiff(manual_clusters, as.character(clusters))
        if (length(missed_clusters) > 0L) {
            cli::cli_abort(c(
                "Cannot find all items provided by {.arg manual} in {.arg clusters}",
                "!" = "Missed items: {.val {missed_clusters}}"
            ))
        }
        is_dup_clusters <- duplicated(manual_clusters, fromLast = TRUE)
        dup_clusters <- unique(manual_clusters[is_dup_clusters])
        if (length(dup_clusters) > 0L) {
            cli::cli_warn(c(
                "Duplicated clusters are provided in {.arg manual}",
                "!" = "Duplicated items: {.val {dup_clusters}}",
                "i" = "will use the later one"
            ))
            manual_clusters <- manual_clusters[!is_dup_clusters]
            manual_labels <- manual_labels[!is_dup_clusters]
        }
    } else {
        manual_clusters <- NULL
    }

    # duplicated markers shouldn't be calculated twice
    cluster2cell <- imap(marker_list, function(markers, i) {
        msg <- sprintf("{.field %s} of {.arg %s}", i, "marker_list")
        if (anyNA(markers)) {
            cli::cli_warn(c(
                sprintf("{.val {NA}} is found in %s", msg),
                "!" = "Will omit {.val NA} value"
            ))
            markers <- markers[!is.na(markers)]
        }
        dup_markers <- unique(markers[duplicated(markers)])
        if (length(dup_markers) > 0L) {
            cli::cli_warn(c(
                sprintf("Duplicated markers are provided in %s", msg),
                "!" = "Duplicated feature{?s}: {.val {dup_markers}}",
                "i" = "Only one will be used"
            ))
            markers <- unique(markers)
        }
        is_existed <- markers %chin% rownames(x) # nolint
        if (!all(is_existed)) {
            cli::cli_warn(c(
                sprintf(
                    "Finding {.val {sum(!is_existed)}} non-existed feature{?s} in %s", msg
                ),
                "!" = "Non-existed feature{?s}: {.val {markers[!is_existed]}}",
                "i" = "Will be omitted directly"
            ))
            markers <- markers[is_existed]
        }
        if (length(markers) == 0L) {
            return(NULL)
        }
        sum_stats <- summarize_features_by_groups(
            x,
            features = markers, groups = clusters,
            statistics = "sum", blocks = blocks
        )
        sums <- sum_stats$statistics$sum
        numbers <- sum_stats$numbers
        data.table::data.table(
            clusters = colnames(sums),
            means = colSums(sums, na.rm = TRUE) / (numbers * nrow(sums))
        )
    })

    # we annotate the cluster as the celltype whose
    cluster2cell <- data.table::rbindlist(
        cluster2cell,
        use.names = FALSE,
        idcol = "celltype"
    )[, list(celltype = celltype[[which.max(means)]]), by = "clusters"][
        , structure(celltype, names = clusters) # nolint
    ]

    # adjust results by manual.
    # we put the annotated value of marker_list in the first and put manual
    # label next to them.
    if (length(manual_clusters) == 0L) {
        cluster2cell_levels <- intersect(names(marker_list), cluster2cell)
    } else {
        cluster2cell[manual_clusters] <- manual_labels
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
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @rdname annotate_clusters
setMethod(
    "annotate_clusters", "SummarizedExperiment",
    function(x, clusters = NULL, ..., assay.type = "logcounts") {
        annotate_clusters_internal(
            SummarizedExperiment::assay(x, assay.type),
            clusters = handle_column_data(object = x, clusters),
            ...
        )
    }
)

#' @export
#' @importClassesFrom SeuratObject Seurat
#' @rdname annotate_clusters
setMethod(
    "annotate_clusters", "Seurat",
    function(x, clusters = NULL, ..., assay.type = NULL) { # nolint
        assay.type <- assay.type %||% SeuratObject::DefaultAssay(x)
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
                assay = assay.type
            ),
            clusters = clusters,
            ...
        )
    }
)
