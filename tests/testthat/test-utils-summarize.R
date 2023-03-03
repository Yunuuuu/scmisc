test_that("summarize_features_by_groups keep the order", {
    set.seed(1L)
    features <- sample(rownames(sce), 100L)
    sce$batch <- factor(sce$batch)
    levels(sce$batch) <- unique(sce$batch)[sample(seq_along(unique(sce$batch)))]
    stats_summary <- summarize_features_by_groups(sce,
        features = features,
        groups = sce$batch,
        statistics = "sum"
    )
    expect_identical(rownames(stats_summary$statistics$sum), features)
    expect_identical(
        colnames(stats_summary$statistics$sum),
        levels(sce$batch)
    )

    stats_summary <- summarize_features_by_groups(sce,
        features = features,
        groups = sce$batch,
        statistics = "sum",
        blocks = sce$block
    )
    expect_identical(rownames(stats_summary$statistics$sum), features)
    expect_identical(
        colnames(stats_summary$statistics$sum),
        levels(sce$batch)
    )
})
