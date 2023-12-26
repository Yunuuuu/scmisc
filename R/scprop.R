#' Permutation Test For Proportions
#'
#' @param identity An atomic vector define the cell types
#' @param compare An atomic vector coerced into factor to define the comparison
#' groups, must have two unique values.
#' @param stat A string giving the type of the statistic to calculate. Current
#' options include "odds_ratio", "risk_ratio", "risk_diff".
#' @param times Conventient way to set both `n_permutation` and `n_bootstrap`.
#' @param n_permutation Number of permutation to calcualte `p.value`.
#' @param n_bootstrap Number of bootstrap to calcualte `ci.low` and `ci.high`.
#' @param p.adjust Correction method, a string. Can be abbreviated. See
#' [p.adjust][stats::p.adjust].
#' @param conf.int Logical indicating whether or not to include a confidence
#' interval in the tidied output. Defaults to `TRUE`.
#' @param conf.level The confidence level to use for the confidence interval if
#' conf.int = `TRUE`. Must be strictly greater than 0 and less than 1. Defaults
#' to 0.95, which corresponds to a 95 percent confidence interval.
#' @param BPPARAM A [BiocParallelParam][BiocParallel::BiocParallelParam] object
#' specifying whether Permutation (and Bootstraping) should be parallelized.
#' @return A [data.table][data.table::data.table]
#' @seealso
#' <https://github.com/rpolicastro/scProportionTest>
#' @export
scprop_compare <- function(
    identity, compare, stat = "risk_ratio", times = 2000L,
    n_permutation = times, n_bootstrap = times,
    p.adjust = "BH", conf.int = TRUE, conf.level = 0.95,
    BPPARAM = BiocParallel::SerialParam()) {
    assert_pkg("matrixStats")
    assert_bool(conf.int)
    assert_(conf.level, function(x) {
        is_scalar_numeric(x) && data.table::between(x, 0, 1L, incbounds = TRUE)
    }, "a scalar numeric in [0, 1]")
    assert_(identity, rlang::is_atomic, "an atomic")
    assert_s4_class(BPPARAM, "BiocParallelParam")
    compare <- as.factor(compare)
    if (nlevels(compare) != 2L) {
        cli::cli_abort(
            "Provided {.arg compare} must have two unique groups to compare"
        )
    }
    if (length(identity) != length(compare)) {
        cli::cli_abort(
            "{.arg identity} and {.arg compare} must have the same length"
        )
    }
    stat <- match.arg(stat, c("risk_ratio", "odds_ratio", "risk_diff"))
    stat <- eval(rlang::sym(sprintf("stat_%s", stat)))

    # Get observed differences in fraction --------
    obs_diff <- scprop_diff(compare, identity, stat)

    # allocate seed -----------------------------------
    seed <- BiocParallel::bpRNGseed(BPPARAM)
    if (conf.int && length(seed) < 2L) {
        set_seed(seed[1L])
        seed <- random_seed(2L)
    }

    # Permutation test ----------------------------
    permuted <- BiocParallel::bplapply(seq_len(n_permutation), function(i) {
        scprop_diff(sample(compare), identity, stat,
            id = paste0("Permutation", i)
        )
    }, BPPARAM = BPPARAM)
    BiocParallel::bpRNGseed(BPPARAM) <- seed[1L]
    permuted <- Reduce(function(x, y) {
        merge(x, y, all = TRUE, by = "identity")
    }, permuted)
    out <- merge(obs_diff, permuted, by = "identity")
    out[, increased := rowSums(sapply(.SD, function(x) x >= estimate)), # nolint
        .SDcols = patterns("^Permutation")
    ]
    out[, decreased := rowSums(sapply(.SD, function(x) x <= estimate)), # nolint
        .SDcols = patterns("^Permutation")
    ]
    out[, c("increased", "decreased") := lapply(.SD, function(x) {
        (x + 1) / (n_permutation + 1)
    }), .SDcols = c("increased", "decreased")]
    out <- out[, .SD, .SDcols = !patterns("^Permutation")]
    if (conf.int) {
        # Boostrap for confidence interval ----------
        data <- data.table::data.table(.compare = compare, .identity = identity)
        BiocParallel::bpRNGseed(BPPARAM) <- seed[2L]
        booted <- BiocParallel::bplapply(seq_len(n_bootstrap), function(i) {
            dd <- data[, list(.identity = sample(.identity, replace = TRUE)), # nolint
                by = ".compare"
            ]
            scprop_diff(dd$.compare, dd$.identity, stat,
                id = paste0("Bootstrap", i)
            )
        }, BPPARAM = BPPARAM)
        booted <- Reduce(function(x, y) {
            merge(x, y, all = TRUE, by = "identity")
        }, booted)
        bootstrap_cols <- setdiff(names(booted), "identity")
        booted[, (bootstrap_cols) := lapply(.SD, function(x) {
            data.table::fifelse(!is.finite(x), NA_real_, x)
        }), .SDcols = bootstrap_cols]
        out <- merge(out, booted, by = "identity")
        out[, booted_mean := rowMeans(.SD, na.rm = TRUE), # nolint
            .SDcols = bootstrap_cols
        ]
        out[, c("ci.low", "ci.high") := data.table::as.data.table(
            matrixStats::rowQuantiles(as.matrix(.SD),
                probs = conf_level_to_prob(conf.level),
                na.rm = TRUE
            )
        ), .SDcols = bootstrap_cols]
        out <- out[, .SD, .SDcols = !bootstrap_cols]
    }
    out[, p.value := data.table::fifelse(estimate > 0, increased, decreased)]
    out[, p.adjust := stats::p.adjust(p.value, p.adjust)] # nolint
    out[, term := levels(compare)[2L]] # nolint
    data.table::setcolorder(out, "term")
    out[, .SD, .SDcols = !c("increased", "decreased")]
}

utils::globalVariables(c(
    "estimate", "increased", "decreased", "p.value", "p.adj",
    ".identity", ".fraction", "booted_mean", "patterns",
    "count", "a", "b", "c", "d", "term"
))

scprop_diff <- function(compare, identity, stat, id = "estimate") {
    data <- data.table::data.table(.compare = compare, .identity = identity)
    data <- data[, list(count = .N), by = c(".compare", ".identity")]
    data <- data.table::dcast(data, .identity ~ .compare, value.var = "count")
    data.table::setnames(data, levels(compare), c("c", "a"))
    data[, c("b", "d") := lapply(.SD, function(x) {
        sum(x) - x
    }), .SDcols = c("a", "c")]
    data[, estimate := stat(a, b, c, d)] # nolint
    data <- data[, c(".identity", "estimate")]
    data.table::setnames(data, c("identity", id))
}

stat_risk <- function(a, b, c, d) {
    log(a / (a + b))
}

stat_odds <- function(a, b, c, d) {
    log(a / b)
}

stat_risk_ratio <- function(a, b, c, d) {
    log((a / (a + b)) / (c / (c + d)))
}

stat_odds_ratio <- function(a, b, c, d) {
    log((a / b) * (d / c))
}

stat_risk_diff <- function(a, b, c, d) {
    (a / (a + b)) - (c / (c + d))
}
