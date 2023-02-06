step_sce_sizefactor <- function(sce, batch = NULL, type = NULL, ...,
                                cluster_param = list(
                                    method = "igraph",
                                    graph.fun = "leiden",
                                    cluster.args = list(
                                        objective_function = "modularity",
                                        resolution_parameter = 1L,
                                        n_iterations = -1L
                                    )
                                ),
                                step_param = list()) {
    type <- match.arg(type, c("deconvolution", "librarysize"))
    sce <- quo_or_symbol(sce)
    batch <- rlang::enquo(batch)
    call <- switch(type,
        deconvolution = rlang::expr(
            scuttle::pooledSizeFactors(
                x = !!sce,
                ...,
                clusters = scran::quickCluster(
                    !!sce,
                    block = !!batch,
                    !!!cluster_param
                ),
                positive = FALSE
            )
        ),
        librarysize = rlang::expr(
            scuttle::librarySizeFactors(
                x = !!sce, ...
            )
        )
    )
    build_step(
        id = "sizefactor", expr = call, step_param,
        list(deps = NULL, finished = FALSE, bind = TRUE, seed = FALSE)
    )
}

utils::globalVariables("sizefactor")

step_sce_normalize <- function(sce, type = NULL, ..., step_param = list()) {
    sce <- quo_or_symbol(sce)
    call <- rlang::expr(scuttle::logNormCounts(
        x = !!sce, size.factors = sizefactor,
        ...
    ))
    build_step(
        id = "normalize", expr = call, step_param,
        list(deps = "sizefactor", finished = FALSE, bind = TRUE, seed = FALSE)
    )
}

step_sce_batchnorm <- function(sce, batch, ..., step_param = list()) {
    sce <- quo_or_symbol(sce)
    batch <- rlang::enquo(batch)
    call <- rlang::expr(batchelor::multiBatchNorm(
        x = !!sce,
        batch = !!batch,
        preserve.single = TRUE,
        normalize.all = TRUE,
        ...
    ))
    build_step(
        id = "batchnorm", expr = call, step_param,
        list(deps = "sizefactor", finished = FALSE, bind = TRUE, seed = FALSE)
    )
}

step_sce_cluster <- function() {

}

step_sce_anno <- function() {

}
