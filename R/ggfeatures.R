#' Create a ggplot from a SingleCellExperiment
#'
#' `ggcells` treats cells as the data values so users can reference row names of
#' `x` (if provided in `features`), column metadata variables and dimensionality
#' reduction results. They can also reference row names and metadata variables
#' for alternative Experiments.
#'
#' `ggfeatures` treats features as the data values so users can reference column
#' names of `x` (if provided in `cells`) and row metadata variables.
#'
#' @inheritParams scuttle::makePerCellDF
#' @inheritParams ggplot2::ggplot
#' @param cells Character vector specifying the features for which to extract
#' expression profiles across cells. Alternative, A named list contaning
#' the cell groups, where names indicate the cell group and values indicate
#' the cells for this group. In this way, data will be transformed into
#' long-format and a column `.features_group` will be added into final data.
#' @inheritDotParams scuttle::makePerCellDF -x -features -assay.type -exprs_values
#' @inheritParams scater::ggfeatures
#' @importFrom ggplot2 aes
#' @export
#' @name make_gg
ggfeatures <- function(x, mapping = aes(), cells = NULL, ..., assay.type = "logcounts", extract_mapping = TRUE) {
    make_gg_data(
        x = x, type = "features", mapping = mapping,
        items = cells, ..., assay.type = assay.type,
        extract_mapping = extract_mapping
    )
}

#' @param features Character vector specifying the features for which to extract
#' expression profiles across cells. May also include features in alternative
#' Experiments if permitted by use.altexps. Alternative, A named list contaning
#' the marker groups, where names indicate the marker group and values indicate
#' the markers for this group. In this way, data will be transformed into
#' long-format and a column `.cells_group` will be added into final data.
#' @rdname make_gg
#' @export
ggcells <- function(x, mapping = aes(), features = NULL, ..., assay.type = "logcounts", extract_mapping = TRUE) {
    make_gg_data(
        x = x, type = "cells", mapping = mapping,
        items = features, ..., assay.type = assay.type,
        extract_mapping = extract_mapping
    )
}

make_gg_data <- function(x, type, mapping = aes(), items = NULL, ..., assay.type = "logcounts", extract_mapping = TRUE, arg = rlang::caller_arg(items)) {
    assert_s4_class(x, "SingleCellExperiment")
    if (is.null(items) || is.character(items)) {
        extracted <- c(items, .aes_in_use(mapping, extract_mapping))
    } else if (is.list(items)) {
        nms <- names(items)
        if (is.null(nms) || anyNA(nms)) {
            cli::cli_abort("All elements in {.arg {arg}} must be named.")
        } else if (anyDuplicated(nms)) {
            cli::cli_abort("Duplicated names found in {.arg {arg}}")
        }
        extracted <- unlist(items, recursive = FALSE, use.names = FALSE)
    } else {
        cli::cli_abort("{.arg arg} must be a {.cls list} or {.cls character}")
    }
    data <- switch(type,
        features = scuttle::makePerFeatureDF(
            x = x, cells = extracted, ..., assay.type = assay.type
        ),
        cells = scuttle::makePerCellDF(
            x = x, features = extracted, ..., assay.type = assay.type
        )
    )
    missed_items <- setdiff(extracted, names(data))
    if (length(missed_items)) {
        cli::cli_abort(c(
            "all {.arg {arg}} must in {.arg x}",
            i = "Missed {arg}: {missed_items}"
        ))
    }
    data <- data.table::as.data.table(data)
    if (is.list(items)) {
        var_nm <- setdiff(c("features", "cells"), type)
        assay.type <- paste0(".", assay.type)
        data <- data.table::melt(data,
            measure.vars = extracted,
            variable.name = paste0(".", var_nm),
            value.name = assay.type
        )
        group_nm <- paste0(".", var_nm, "_class")
        n_group <- switch(type,
            cells = ncol(x),
            features = nrow(x)
        )
        data[, (group_nm) := rep(rep(names(items), times = lengths(items)),
            each = n_group
        )]
    }
    ggplot2::ggplot(data[], mapping)
}

.aes_in_use <- function(mapping, extract_mapping) {
    collected <- NULL
    if (extract_mapping) {
        for (i in seq_along(mapping)) {
            collected <- c(collected, rlang::as_label(mapping[[i]]))
        }
    }
    collected
}
