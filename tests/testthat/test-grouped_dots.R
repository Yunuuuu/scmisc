test_that("grouped_dots works well", {
    sce <- scuttle::logNormCounts(sce)
    set.seed(1L)
    testthat::expect_s4_class(plot_grouped_dots(sce,
        list(a = rownames(sce)[sample(nrow(sce), 10L)]),
        groups = "batch"
    ), "DotsHeatmap")
    testthat::expect_s4_class(plot_grouped_dots(sce,
        list(a = rownames(sce)[sample(nrow(sce), 10L)]),
        groups = "batch", dots_size_legend_param = NULL
    ), "Heatmap")
    testthat::expect_s4_class(plot_grouped_dots(sce,
        list(a = rownames(sce)[sample(nrow(sce), 10L)]),
        groups = "batch", blocks = "block"
    ), "DotsHeatmap")
    testthat::expect_s4_class(plot_grouped_dots(sce,
        list(a = rownames(sce)[sample(nrow(sce), 10L)]),
        groups = "batch", , blocks = "block",
        dots_size_legend_param = NULL
    ), "Heatmap")
})
