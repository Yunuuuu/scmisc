#' Seach CellMarker database
#'
#' Search http://xteam.xbio.top/CellMarker.
#' @param markers an atomic character, the markers to search in the CellMarker
#' database, can be the Gene Symbol, Gene ID, Protein Symbol or Protein ID
#' (usually starts with "P", "Q" or "O").
#' @param species a scalar string, "human" or "mouse".
#' @param internal logical value, indicates whether to use internal CellMarker
#' data. If `NULL`, this will be determined automatically; if the CellMarker
#' data has been downloaded once, namely, we have already used this function
#' once with a `internal` value `FALSE`, then the `NULL` will indicate `FALSE`.
#' Otherwise `TRUE`. The internal data was downloaded from CellMarker on
#' 2023-01-10.
#' @return a data.frame of the searching results. A column named `targeted`
#' containing the matched markers from CellMarker data, the row will be sorted
#' descendingly by the number of matched markers.
#' @export
cellmarker_search <- function(markers, species = "human", internal = NULL) {
    data <- data.table::copy(cellmarker_get(species, internal))
    # Since gene_list may contain duplicated values with different case
    # it's better to use the input markers as the reference
    data[, targeted := lapply(gene_list, function(.genes, .markers) {
        .markers[match(
            tolower(.markers), tolower(.genes),
            nomatch = 0L, incomparables = NA_character_
        ) > 0L]
    }, .markers = markers)]
    data[, targeted_size := lengths(targeted)]
    data[, targeted_prop := targeted_size / length(markers)]
    data.table::setcolorder(
        data, c("targeted", "targeted_size", "targeted_prop"),
        after = "CellOntologyID"
    )
    data.table::setcolorder(
        data, intersect(cellmarker_gene_cols, names(data)),
        after = "gene_list"
    )
    data <- data[targeted_size > 0L][
        order(-targeted_size, na.last = TRUE)
    ]
    data.table::setDF(data)[]
}

utils::globalVariables(
    c("targeted", "gene_list", "targeted_size", "targeted_prop")
)

cellmarker_get <- function(species = "human", internal = NULL) {
    species <- match.arg(species, c("human", "mouse"))
    if (is.null(internal)) {
        if (exists(species, where = cellmarker_database_external, inherits = FALSE)) {
            internal <- FALSE
        } else {
            internal <- TRUE
        }
    }
    if (internal) {
        cellmarker_database[[species]] # nolint
    } else {
        cellmarker_download(species)
    }
}

cellmarker_prepare <- function(data) { # nolint styler: off
    data[, gene_list := .mapply( # nolint styler: off
        function(...) {
            Reduce(union, list(...))
        },
        unname(lapply(.SD, function(markers) { # nolint styler: off
            markers_list <- strsplit(
                gsub("\\s*\\[\\s*|\\s*\\]\\s*", "", markers, perl = TRUE),
                "\\s*,\\s*",
                perl = TRUE
            )
            lapply(markers_list, function(markers_trim) {
                markers_trim <- trimws(markers_trim, "both", "[\\h\\v]")
                markers_trim[!is.na(markers_trim) & markers_trim != "NA"]
            })
        })),
        MoreArgs = NULL
    ), .SDcols = intersect(cellmarker_gene_cols, names(data))]
}

utils::globalVariables("gene_list")

cellmarker_download <- function(species) {
    if (!exists(species, where = cellmarker_database_external, inherits = FALSE)) {
        data_link <- switch(species,
            human = "http://xteam.xbio.top/CellMarker/download/Human_cell_markers.txt",
            mouse = "http://xteam.xbio.top/CellMarker/download/Mouse_cell_markers.txt"
        )
        cli::cli_alert_info("Reading data from {.url {data_link}}")
        cellmarker_database_external[[species]] <- cellmarker_prepare(
            data.table::fread(data_link)
        )
        # envir <- topenv(environment(NULL))
        # unlockBinding("cellmarker_database_external", envir)
        # utils::assignInMyNamespace(
        #     "cellmarker_database_external",
        #     cellmarker_database_external
        # )
        # lockBinding("cellmarker_database_external", envir)
    }
    get(species, pos = cellmarker_database_external, inherits = FALSE)
}

# cellmarker_database_external <- list(
#     human = NULL, mouse = NULL
# )
cellmarker_database_external <- new.env(parent = emptyenv())
cellmarker_gene_cols <- c(
    "cellMarker", "geneSymbol", "geneID", "proteinName", "proteinID"
)
