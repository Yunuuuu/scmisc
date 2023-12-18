#' Find roots in a DiffusionMap object
#'
#' @description
#' Find roots in a DiffusionMap object. If we know the `start` and `end` of a
#' root, we then check if the tips of the root are in the start and ends node we
#' already known (see [find_tips][destiny::find_tips]). This function in this
#' way check `n_root` nodes with the top `DPT` in the start cluster and return
#' roots meet this criteria.
#' @param dm A [DiffusionMap][destiny::DiffusionMap] object.
#' @param start,ends The start and ends cluster identity. start must have a
#'   length 1L, while the length of ends ranges from 1L to 2L. All start and
#'   ends must exist in `ref`. If ends is `NULL`, this only check if one tip is
#'   in the start cluster.
#' @param ref All cell identity of `dm`, this must have the same length of
#'   `dm@@d`.
#' @param n_root The number of nodes we should test.
#' @return An integer index .
#' @export
find_roots <- function(dm, start, ends = NULL, ref, n_root = 100L) {
    assert_s4_class(dm, "DiffusionMap")
    assert_(start, is_scalar, "scalar", null_ok = FALSE)
    assert_(n_root, is_scalar_numeric, "a number", null_ok = FALSE)

    if (length(dm@d) != length(ref)) {
        cli::cli_abort("{.arg ref} must have the same length {.code length(dm@d)} ({length(dm@d)})")
    }
    if (!all(c(start, ends) %in% ref)) {
        cli::cli_abort("{.arg start} and {.arg ends} must exist in {.arg ref}.")
    }
    if (!is.null(ends)) {
        if (length(ends) > 2L || length(ends) < 1L) {
            cli::cli_abort("{.arg ends} must be {.code NULL} or the length of {.arg ends} must ranges from 1L to 2L.")
        }
    }

    # we extract the top n_root with maximal DPT from start notes.
    roots <- which(ref == start, useNames = FALSE)
    dpt <- methods::new("DPT",
        branch = matrix(), tips = matrix(),
        dm = dm
    )
    dpt <- dpt[sample(roots, size = 1L)][roots]
    n_root <- min(n_root, length(roots), na.rm = TRUE)
    roots <- roots[order(dpt, decreasing = TRUE)][seq_len(n_root)]

    roots[vapply(
        roots, is_matched_root, logical(1L),
        dm = dm, start = start,
        ends = ends, ref = ref
    )]
}

# root is the index
# start, end1 and end2 are all a single string.
# ref is the reference of start, end1 and end2.
is_matched_root <- function(dm, root, start, ends, ref) {
    tips <- destiny::find_tips(dm, root = root)
    if (is.null(ends)) {
        any(tips %in% start, na.rm = TRUE)
    } else {
        setequal(ref[tips], c(start, ends))
    }
}
