`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

pindex <- function(array, ...) {
    if (length(dim(array)) != ...length()) {
        stop("Indexing must have as many as the number of dimentions of array")
    }
    dots <- list(...)
    # all index must be atomic
    is_right <- vapply(dots, function(x) {
        is.atomic(x) && !is.null(x)
    }, logical(1L))
    if (!all(is_right)) {
        stop("All elements in ... must be atomic (`NULL` is not allowed)")
    }
    dots_len <- lengths(dots)
    if (any(dots_len == 0L)) {
        stop("Empty index is not allowed")
    }
    common_len <- max(dots_len)
    if (any(dots_len > 1L & dots_len < common_len)) {
        stop("Only index of length one are recycled")
    }
    if (common_len != 1L) {
        dots[dots_len == 1L] <- lapply(dots[dots_len == 1L], function(x) {
            rep_len(x, common_len)
        })
    }
    array[do.call("cbind", dots)]
}

is_scalar <- function(x) {
    length(x) == 1L
}

is_scalar_numeric <- function(x) {
    length(x) == 1L && is.numeric(x)
}

#' Implement similar purrr::imap function
#' @note this function won't keep the names of .x in the final result
#' @keywords internal
#' @noRd
imap <- function(.x, .f, ...) {
    names_idx <- names(.x) %||% as.character(seq_along(.x))
    out <- .mapply(
        rlang::as_function(.f),
        dots = list(.x, names_idx),
        MoreArgs = list(...)
    )
    names(out) <- names_idx
    out
}

#' Rename elements in a list, data.frame or vector
#'
#' This is akin to `dplyr::rename` and `plyr::rename`. It renames elements given
#' as names in the `replace` vector to the values in the `replace` vector
#' without touching elements not referenced.
#'
#' @param x A data.frame or a named vector or list
#' @param replace A named character vector. The names identifies the elements in
#' `x` that should be renamed and the values gives the new names.
#'
#' @return `x`, with new names according to `replace`
#'
#' @keywords internal
#' @noRd
rename <- function(x, replace) {
    current_names <- names(x)
    old_names <- names(replace)
    missing_names <- setdiff(old_names, current_names)
    if (length(missing_names)) {
        replace <- replace[!old_names %in% missing_names]
        old_names <- names(replace)
    }
    names(x)[match(old_names, current_names)] <- as.vector(replace)
    x
}

# `NULL` passed to ... will also be kept
#' @param x A list to modify
#' @param replace A named list with components to replace (or add) corresponding
#' components in `x`.
#' @noRd
modify_list <- function(x, replace) {
    for (name in names(replace)) {
        value <- replace[[name]]
        if (rlang::is_zap(value)) {
            x[[name]] <- NULL
        } else {
            x[name] <- list(value)
        }
    }
    x
}

check_dots_named <- function(..., call = parent.frame()) {
    dots <- rlang::dots_list(..., .named = NULL, .homonyms = "error")
    if (!rlang::is_named2(dots)) {
        cli::cli_abort("All elements in {.arg ...} must be named", call = call)
    }
    dots
}

# call object utils ----------------------------------------
call_standardise <- function(call, env = parent.frame()) {
    expr <- rlang::get_expr(call)
    env <- rlang::get_env(call, env)
    fn <- rlang::eval_bare(rlang::node_car(expr), env)
    if (rlang::is_primitive(fn)) {
        call
    } else {
        matched <- rlang::call_match(expr, fn, defaults = TRUE)
        rlang::set_expr(call, matched)
    }
}

# find_expr_deps <- function(expr) {
#     # substitute "~" with "base::identity"
#     expr <- rlang::call2(rlang::expr(substitute), expr,
#         env = rlang::expr(list(`~` = base::identity)),
#         .ns = "base"
#     )
#     expr <- rlang::eval_bare(expr, env = rlang::base_env())
#     codetools::findGlobals(rlang::new_function(args = NULL, body = expr))
# }

change_expr <- function(expr, from, to) {
    switch(expr_type(expr),

        # special cases
        # for missing argument in pairlist
        missing = rlang::missing_arg(),

        # seems like this will always live in the end of anonymous function call
        # we just return as it is
        integer = expr,

        # Base cases
        constant = ,
        symbol = if (identical(expr, from)) to else expr,

        # Recursive cases
        call = as.call(lapply(expr, change_expr, from = from, to = to)),
        pairlist = as.pairlist(lapply(expr, change_expr, from = from, to = to)),
        cli::cli_abort("Don't know how to handle type {.cls {expr_type(expr)}}")
    )
}

expr_type <- function(x) {
    if (rlang::is_missing(x)) {
        "missing"
    } else if (rlang::is_syntactic_literal(x)) {
        "constant"
    } else if (is.symbol(x)) {
        "symbol"
    } else if (is.call(x)) {
        "call"
    } else if (is.pairlist(x)) {
        "pairlist"
    } else {
        typeof(x)
    }
}

# environment helper
# Wrapper around list2env with a NULL check. In R <3.2.0, if an empty unnamed
# list is passed to list2env(), it errors. But an empty named list is OK. For
# R >=3.2.0, this wrapper is not necessary.
# @param empty_to_null Controls what to do when x is NULL or empty list.
#   If TRUE, return NULL. If FALSE, return an empty list.
list2env2 <- function(x, envir = NULL, parent = emptyenv(),
                      hash = (length(x) > 100),
                      size = max(29L, length(x))) {
    if (is.null(envir)) {
        envir <- new.env(hash = hash, parent = parent, size = size)
    }
    if (length(x) == 0L) {
        return(envir)
    }
    list2env(x, envir)
}

# cli output helper
cli_list <- function(label, items, sep = ": ", add_num = TRUE) {
    items <- cli::cli_vec(items, list("vec-trunc" = 3L))
    message <- "{.field {label}}{sep}{.val {items}}"
    if (add_num) {
        message <- paste(message, "({length(items)} item{?s})", sep = " ")
    }
    cli::cli_li(message)
}

conf_level_to_prob <- function(x) {
    x * c(-1L, 1L) / 2L + 0.5
}

set_seed <- function(seed = NULL, envir = rlang::caller_env()) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        oseed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
        oseed <- NULL
    }
    run <- substitute(on.exit(restore_rng(oseed)), list(oseed = oseed))
    eval(run, envir = envir)
    seed <- seed %||% random_seed(1L)
    set.seed(seed)
}

restore_rng <- function(oseed) {
    if (is.null(oseed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        }
    } else {
        assign(".Random.seed", oseed, envir = .GlobalEnv, inherits = FALSE)
    }
}

random_seed <- function(n) {
    sample.int(1e6L, n)
}
