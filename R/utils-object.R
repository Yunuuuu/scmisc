# for the feature data, we usually want to all or some features (rownames)
#' @keywords internal
handle_row_data <- function(object, idx) {
    if (is.null(idx)) {
        rownames(object)
    } else if (rlang::is_scalar_character(idx) && ncol(object) > 1L) {
        get_rowData_column(object, idx)
    } else {
        index_features(object, idx)
    }
}

# for the column data, we usually want to a characters or factors with the
# length equal to `ncol(object)`
#' @keywords internal
handle_column_data <- function(object, idx, arg = rlang::caller_arg(idx)) {
    if (is.null(idx)) {
        idx <- get_colData_column(object, "label")
    } else if (rlang::is_scalar_character(idx) && ncol(object) > 1L) {
        idx <- get_colData_column(object, idx)
    } else {
        if (length(idx) != nrow(object)) {
            cli::cli_abort("{.arg {arg}} should be of length {.val {nrow(object)}}.")
        }
    }
    idx
}

#' @keywords internal
swap_rownames <- function(object, swap_rownames = NULL) {
    if (is.null(swap_rownames)) {
        return(object)
    } else if (rlang::is_scalar_character(swap_rownames) && nrow(object) > 1L) {
        rownames(object) <- get_rowData_column(object, swap_rownames)
    } else {
        if (length(swap_rownames) != nrow(object)) {
            cli::cli_abort("{.arg swap_rownames} should be {.val NULL} or a string or of length {.val {nrow(object)}} characters.")
        }
        rownames(object) <- as.character(swap_rownames)
    }
    object
}

get_rowData_column <- function(object, column) {
    if (!any(column == colnames(SummarizedExperiment::rowData(object)))) {
        cli::cli_abort("Cannot find column {.val {column}} in rowData")
    }
    SummarizedExperiment::rowData(object)[[column]]
}

get_colData_column <- function(object, column) {
    if (!any(column == colnames(SummarizedExperiment::colData(object)))) {
        cli::cli_abort("Cannot find column {.val {column}} in colData")
    }
    SummarizedExperiment::colData(object)[[column]]
}

index_features <- function(object, idx) {
    if (is.logical(idx)) {
        if (length(idx) != nrow(object)) {
            cli::cli_abort("logical features should be of length {.val {nrow(object)}}.")
        }
        idx <- rownames(object)[idx]
    } else if (is.numeric(idx)) {
        if (any(idx > nrow(object)) || any(idx <= 0) ||
            any(idx != round(idx))) {
            cli::cli_abort("All features should be round numbers; > 0 and <= {.val {nrow(object)}}.")
        }
        idx <- rownames(object)[idx]
    } else if (is.character(idx) || is.factor(idx)) {
        idx <- as.character(idx)
        if (!all(idx %in% rownames(object))) {
            cli::cli_abort("Some features are not in the input object.")
        }
    } else {
        cli::cli_abort("Unsupported features type: {.cls {typeof(idx)}}")
    }
    idx
}
