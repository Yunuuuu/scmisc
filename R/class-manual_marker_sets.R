#' Cell markers sets collected from research articles
#'
#' @description
#' We manullay collected some cell marker sets from research articles, these
#' funciton juse provide some methods to operate with the manual cell marker
#' sets (manual_cms).
#'
#' `manual_cms` function just returns the database with a class of
#' `manual_cell_marker_datasets` which has a specific `print` method.
#'
#' `get_markers` function extract the specific marker set from the datasets,
#' which return a `marker_set` class object.
#' @param name The name of the marker set to extract. `name` usually is the
#' surname of the first author followed by the years of the published article
#' connected with "_".
#' @param x An `manual_cms` or `marker_set` object.
#' @param ... Not used currently
#' @name manual_cms
NULL

#' @export
#' @rdname manual_cms
get_markers <- function(name) {
    get(name, pos = manual_cell_marker_datasets, inherits = FALSE)
}

to_marker_set <- function(name, ..., reference) {
    x <- marker_set(..., reference = reference)
    add_marker_set(name = name, marker_set = x)
}

#' @export
#' @rdname manual_cms
manual_cms <- function() {
    manual_cell_marker_datasets
}

#' @export
#' @rdname manual_cms
print.manual_cms <- function(x, ...) {
    cli::cli_text("A total of {rlang::env_length(x)} cell markers set{?s}")
}

manual_cell_marker_datasets <- structure(
    new.env(parent = emptyenv()),
    class = "manual_cms"
)

marker_set <- function(..., reference) {
    x <- rlang::dots_list(..., .named = NULL, .homonyms = "error")
    x <- new_marker_set(x, reference = reference)
    validate_marker_set(x)
    x
}

new_marker_set <- function(list, reference) {
    stopifnot(is.list(list))
    structure(
        list,
        class = "marker_set",
        reference = reference
    )
}

# convention: name should be the surname of the first author followed by the
# years of the published article connected with "_"
add_marker_set <- function(name, marker_set) {
    if (exists(name, where = manual_cell_marker_datasets, inherits = FALSE)) {
        cli::cli_abort("Existed cell marker set found in datasets")
    }
    manual_cell_marker_datasets[[name]] <- marker_set
}

#' @export
#' @rdname manual_cms
print.marker_set <- function(x, ...) {
    cli::cli_text("Cell marker set derived from: {.url {attr(x, \"reference\")}}")
    cli::cli_rule("marker set details")
    # cat main markers if it exits
    if (any("main" == names(x))) {
        main_par_id <- cli::cli_par()
        cli::cli_text("Main cell types")
        main_item_lid <- cli::cli_ul()
        .mapply(
            cli_list,
            list(
                label = names(x$main),
                items = x$main
            ),
            MoreArgs = list(add_num = FALSE)
        )
        cli::cli_end(id = main_item_lid)
        cli::cli_end(id = main_par_id)
    }
    # cat other markers
    other_items <- x[setdiff(names(x), "main")]
    if (length(other_items)) {
        cli::cli_ul()
        cli::cli_text("Including markers for {length(other_items)} group{?s}")
        .mapply(
            cli_list,
            list(
                label = names(other_items),
                items = lapply(other_items, names)
            ),
            MoreArgs = NULL
        )
        cli::cli_end()
    }
    invisible(x)
}

validate_marker_set <- function(x) {
    x_type <- typeof(x)
    if (x_type == "character") {
        # No need names for character
        is_right <- TRUE
    } else if (x_type == "list") {
        # for a list, we should check all elements have names,
        # and then recall this function to check every elments
        is_right <- rlang::is_named(x) &&
            all(vapply(x, validate_marker_set, logical(1L), USE.NAMES = FALSE))
    } else {
        cli::cli_abort("all elements must be a {.cls list} or a {.cls chracter}")
    }
    if (!is_right) {
        cli::cli_abort("all elements or sub-elements should be named")
    }
    is_right
}

wrap_cat <- function(label, items, sep = ": ", collapse = ", ", indent = 0L, exdent = 2L) {
    total <- length(items)

    ext <- if (total == 0L) {
        "none"
    } else if (total == 1L) {
        items
    } else if (total <= 6L) {
        paste(
            paste(items[-length(items)], collapse = collapse),
            "and", items[[length(items)]],
            sep = " "
        )
    } else {
        paste(
            paste(
                paste(items[1:3], collapse = collapse),
                "...",
                paste(items[(total - 1L):total], collapse = collapse),
                sep = ", "
            ),
            sprintf("(%d total)", total),
            sep = " "
        )
    }
    cat(strwrap(
        paste(label, ext, sep = sep),
        indent = indent, exdent = exdent
    ), sep = "\n")
}
