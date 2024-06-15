#' Run infercnv workflow
#'
#' @param counts Single cell gene expression counts matrix.
#' @param granges A [GRanges][GenomicRanges::GRanges-class] object indicates the
#' positions of each gene. It's also possible to provide a
#' [EnsDb][ensembldb::EnsDb] or [TxDb][GenomicFeatures::TxDb-class] object, in
#' this way, `keytype` will be used to define the gene identifiers.
#' @param groups A character indicates the cell type of each column.
#' @param ref_groups Groups used as the reference.
#' @param keytype The gene type in the counts matrix, either "gene_id" or
#' "symbol". This parameter is only applicable when the `granges` object is an
#' [EnsDb][ensembldb::EnsDb] with the default value set to "symbol", or a
#' [TxDb][GenomicFeatures::TxDb-class] with the default value set to "gene_id".
#' When "gene_id" is used, it represents the ensemble id for
#' [EnsDb][ensembldb::EnsDb] or entrez id for
#' [TxDb][GenomicFeatures::TxDb-class].
# ' Only "gene_id" can be used for [TxDb][GenomicFeatures::TxDb-class].
#' @param standard_chr A boolean value indicates whether keep only 'standard'
#' chromosomes. See
#' [keepStandardChromosomes][GenomeInfoDb::keepStandardChromosomes].
#' @param species A string indicates the genus and species of the organism.
#' Supported species can be seen with `names(GenomeInfoDb::genomeStyles())`.
#' Only used when standard_chr is `TRUE`.
#' @param seqstyle A character string that sets the seqlevels style for
#' `granges`.
#' @inheritParams infercnv::CreateInfercnvObject
#' @inheritParams infercnv::run
#' @inheritDotParams infercnv::run -infercnv_obj -png_res
#' @seealso
#' - [infercnv::CreateInfercnvObject()]
#' - [infercnv::run()]
#' @references
#' <https://github.com/broadinstitute/infercnv/wiki>
run_infercnv <- function(
    counts, granges, groups, ref_groups = character(),
    keytype = NULL, standard_chr = TRUE, species = NULL,
    seqstyle = "UCSC", chr_exclude = c("chrX", "chrY", "chrM"),
    max_cells_per_group = NULL,
    min_max_counts_per_cell = c(100, Inf),
    png_res = 600L, ...) {
    assert_pkg("infercnv")
    assert_pkg("GenomicRanges")
    assert_pkg("GenomicFeatures")
    assert_pkg("GenomeInfoDb")
    assert_pkg("BiocGenerics")
    assert_(counts, function(x) inherits(x, c("dgCMatrix", "matrix")))
    assert_bool(standard_chr)
    if (ncol(counts) != length(groups)) {
        cli::cli_abort(paste(
            "the length of {.arg groups}",
            "must be compatible with {.code ncol(x)}"
        ))
    }
    # The first column is the cell name,
    # and the 2nd column indicates the known cell type.
    anno <- data.frame(
        celltypes = as.character(groups),
        row.names = colnames(counts)
    )

    # prepare gene order file -------------------------
    if (methods::is(granges, "EnsDb")) {
        keytype <- keytype %||% "symbol"
        if (!identical(keytype, "gene_id") && !identical(keytype, "symbol")) {
            cli::cli_abort(paste(
                "{.arg keytype} must be {.val gene_id} (ENSEMBL ID)",
                "or {.val symbol} (GENE SYMBOL) for {.cls EnsDb}"
            ))
        }
        granges <- GenomicFeatures::genes(granges, columns = keytype)
        granges <- granges[granges[[keytype]] != ""]
        if (standard_chr) {
            granges <- GenomeInfoDb::keepStandardChromosomes(
                granges,
                species = species,
                pruning.mode = "tidy"
            )
        }
        granges <- GenomicRanges::reduce(
            split(granges, granges[[keytype]]),
            ignore.strand = TRUE
        )
        granges <- unlist(granges[lengths(granges) == 1])
    } else if (methods::is(granges, "TxDb")) {
        keytype <- keytype %||% "gene_id"
        if (!identical(keytype, "gene_id")) {
            cli::cli_abort(paste(
                "{.arg keytype} must be",
                "{.val gene_id} (ENTREZID) for {.cls TxDb}"
            ))
        }
        granges <- GenomicFeatures::genes(
            granges,
            columns = keytype,
            single.strand.genes.only = FALSE
        )
        if (standard_chr) {
            granges <- GenomeInfoDb::keepStandardChromosomes(
                granges,
                species = species,
                pruning.mode = "tidy"
            )
        }
        granges <- GenomicRanges::reduce(granges, ignore.strand = TRUE)
        granges <- unlist(granges[lengths(granges) == 1])
    } else if (methods::is(granges, "GRanges")) {
        if (standard_chr) {
            granges <- GenomeInfoDb::keepStandardChromosomes(
                granges,
                species = species,
                pruning.mode = "tidy"
            )
        }
    } else {
        cli::cli_abort(
            "{.arg granges} must be a",
            "{.cls GRanges} or {.cls EnsDb} or {.cls TxDb}"
        )
    }
    GenomeInfoDb::seqlevelsStyle(granges) <- seqstyle
    gene_order <- data.frame(
        seqnames = GenomeInfoDb::seqnames(granges),
        start = BiocGenerics::start(granges),
        end = BiocGenerics::end(granges),
        row.names = names(granges)
    )
    gene_order <- gene_order[
        order(gene_order$seqnames, gene_order$start, gene_order$end), ,
        drop = FALSE
    ]
    gene_order$seqnames <- as.character(gene_order$seqnames)
    infercnv_obj <- infercnv::CreateInfercnvObject(
        raw_counts_matrix = counts,
        gene_order_file = gene_order,
        annotations_file = anno,
        ref_group_names = as.character(ref_groups),
        max_cells_per_group = max_cells_per_group,
        min_max_counts_per_cell = min_max_counts_per_cell,
        chr_exclude = as.character(chr_exclude)
    )
    infercnv::run(infercnv_obj = infercnv_obj, ..., png_res = png_res)
}
