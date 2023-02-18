#' Find roots in a DiffusionMap object
#'
#' @description
#' Find roots in a DiffusionMap object. If we know the `start` and `end` of a
#' root, we then check if the tips of the root are in the start and ends node we
#' already known. This function in this way check `n_root` nodes with the top
#' `GPT` in the start cluster and return roots meet this criteria.
#' @param dm A [DiffusionMap][destiny::DiffusionMap] object.
#' @param start,ends The start and ends cluster identity. start must have a
#'   length 1L, while the length of ends ranges from 1L to 2L. All start and
#'   ends must exist in `ref`.
#' @param ref All cell identity of `dm`, this must have the same length of
#'   `dm@@d`.
#' @param n_root The number of nodes we should test.
#' @return An integer index .
#' @export
find_roots <- function(dm, start, ends, ref, n_root = 100L) {
    if (!methods::is(dm, "DiffusionMap")) {
        cli::cli_abort("{.arg dm} must be a {.cls DiffusionMap}.")
    }
    assert_length(start, 1L, null_ok = FALSE)
    assert_length(n_root, 1L, null_ok = FALSE)
    assert_class(n_root, is.numeric, "numeric", null_ok = FALSE)
    if (length(ends) > 2L || length(ends) < 1L) {
        cli::cli_abort("the length of {.arg ends} must ranges from 1L to 2L.")
    }
    if (length(dm@d) != length(ref)) {
        cli::cli_abort("{.arg ref} must have the same length {.code length(dm@d)} ({length(dm@d)})")
    }
    if (!all(c(start, ends) %in% ref)) {
        cli::cli_abort("{.arg start} and {.arg ends} must exist in {.arg ref}.")
    }
    # we extract the top n_root with maximal DPT from start notes.
    roots <- which(ref == start, useNames = FALSE)
    dpt <- methods::new("DPT",
        branch = matrix(), tips = matrix(), 
        dm = dm
    )
    dpt <- dpt[sample(roots, size = 1L)][roots]
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
    setequal(ref[tips], c(start, ends))
}
