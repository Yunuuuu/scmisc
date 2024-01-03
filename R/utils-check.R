# S4 object -------------------------------------
assert_s4_class <- function(
    x, is_class, what, ...,
    arg = rlang::caller_arg(x),
    call = rlang::caller_env()) {
    if (rlang::is_string(is_class)) {
        class <- is_class
        is_class <- function(x) {
            methods::is(x, class)
        }
        if (missing(what)) {
            what <- class
        }
    }
    assert_(
        x = x, assert_fn = is_class, what = what,
        ..., arg = arg, call = call
    )
}

assert_data_frame_columns <- function(x, columns, ..., args = rlang::caller_arg(x), call = rlang::caller_env()) {
    missing_cols <- setdiff(columns, names(x))
    if (length(missing_cols)) {
        args <- style_arg(args)
        if (length(args) == 1L) {
            msg <- args
        } else {
            msg <- sprintf("One of %s", oxford_comma(args, final = "or"))
        }
        rlang::abort(
            c(
                sprintf(
                    "%s must contain columns: %s", msg,
                    oxford_comma(columns)
                ),
                x = sprintf("missing columns: %s", oxford_comma(missing_cols))
            ),
            ...,
            call = call
        )
    }
}
