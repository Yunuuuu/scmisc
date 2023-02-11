#' @param groups A factor or vector specifying the group identity for each
#' column of x, usually clusters or cell types.
#' @param statistics A Character vector specifying the type of statistics to be
#' computed, see Details. `c("mean", "sum", "num.detected", "prop.detected",
#' "median")`
#' @param blocks A factor or vector specifying the blocking level for each
#' column of x, e.g., batch of origin.
#' @param threshold A numeric scalar specifying the threshold above which a gene
#' is considered to be detected.
#' @keywords internal
#' @noRd 
summarize_features_by_groups <- function(x, features, groups, statistics, blocks = NULL, threshold = 0L, id = NULL, check_dup = TRUE) {
    if (is.null(colnames(x))) {
        colnames(x) <- seq_len(ncol(x))
    }
    ids <- S4Vectors::DataFrame(groups = groups)
    if (!is.null(blocks)) {
        ids$blocks <- blocks
    }
    msg <- id %||% ""
    if (anyNA(features)) {
        cli::cli_warn(c(
            sprintf("{.val {NA}} is found in %s.", msg),
            "i" = "will omit {.val {NA}}"
        ))
        features <- features[!is.na(features)]
    }

    if (check_dup) {
        dup_features <- unique(features[duplicated(features)])
        if (length(dup_features) > 0L) {
            cli::cli_warn(c(
                sprintf("Duplicated features are provided in %s.", msg),
                "x" = "Duplicated items: {.val {dup_features}}",
                "i" = "will use only once"
            ))
            features <- unique(features)
        }
    }

    # calculate the statistics then corrected by blocks
    stat_se <- scuttle::summarizeAssayByGroup(
        x,
        ids = ids,
        subset.row = features,
        statistics = statistics,
        store.number = "ncells",
        threshold = threshold
    )
    numbers <- stat_se$ncells
    if (!is.null(blocks)) {
        stat_matrix_list <- lapply(
            SummarizedExperiment::assayNames(stat_se),
            function(assay) {
                transform <- switch(assay,
                    mean = ,
                    sum = "raw",
                    num.detected = , 
                    prop.detected = "logit"
                )
                scuttle::correctGroupSummary(
                    SummarizedExperiment::assay(stat_se, assay),
                    group = stat_se$groups,
                    block = stat_se$blocks,
                    transform = transform
                )
            }
        )
        numbers <- data.table::data.table(
            .numbers. = numbers,
            .groups. = stat_se$groups
        )[, list(.sums. = sum(.numbers., na.rm = TRUE)), by = ".groups."][
            , structure(.sums., names = as.character(.groups.))
        ]
        numbers <- numbers[as.character(stat_se$groups)]
    } else {
        colnames(stat_se) <- stat_se$groups
        stat_matrix_list <- SummarizedExperiment::assays(stat_se)
    }
    # rows are genes; columns are groups
    list(statistics = stat_matrix_list, numbers = numbers)
}

utils::globalVariables(c(".sums.", ".groups.", ".numbers."))