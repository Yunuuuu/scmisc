#' @param statistic Character vector specifying the type of statistics to be
#' computed, see Details. `c("mean", "sum", "num.detected", "prop.detected",
#' "median")`
#' @param ... Other arguments passed to `correctGroupSummary`.
#' @keywords internal
summarize_features_by_groups <- function(x, features, groups, statistic, blocks = NULL, id = NULL, check_dup = TRUE, threshold = 0L, ...) {
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
        statistic = statistic,
        store.number = "ncells",
        threshold = threshold
    )
    numbers <- stat_se$ncells
    if (!is.null(blocks)) {
        stat_matrix <- scuttle::correctGroupSummary(
            SummarizedExperiment::assay(stat_se, statistic),
            group = stat_se$groups,
            block = stat_se$blocks,
            ...
        )
    } else {
        colnames(stat_se) <- stat_se$groups
        stat_matrix <- SummarizedExperiment::assay(stat_se, statistic)
        numbers <- aggregate(numbers,
            by = list(groups = stat_se$groups),
            "sum", simplify = TRUE, drop = FALSE,
            na.rm = TRUE
        )
        numbers <- numbers$x[
            match(colnames(stat_matrix), numbers$groups)
        ]
    }
    # rows are genes; columns are groups
    list(statistic = stat_matrix, numbers = numbers)
}
