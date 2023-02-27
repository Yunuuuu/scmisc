#' Run cellranger-count
#' @description Count gene expression (targeted or whole-transcriptome) and/or feature barcode reads from a single sample and GEM well
#' @param id A unique run id and output folder name.
#' @param fastqs  Path to input FASTQ data.
#' @param transcriptome Path of folder containing 10x-compatible transcriptome
#'   reference. 
#' @param sample Prefix of the filenames of FASTQs to select.
#' @param cores Set max cores the pipeline may request at one time. Only applies
#'   to local jobs. 
#' @param other_args Other arguments passed to cellranger count.
#' @inheritParams run_cellranger
#' @export 
run_cellranger_count <- function(id, fastqs, transcriptome = NULL, sample = NULL, cores = NULL, other_args = character("--nosecondary"), cellranger_cmd = NULL, sys_args = list()) {
    args <- c(
        handle_arg(id, "--id", sep = "="),
        handle_arg(fastqs, "--fastqs", sep = "="),
        handle_arg(transcriptome, "--transcriptome", sep = "="),
        handle_arg(sample, "--sample", sep = "="),
        handle_arg(cores, "--localcores", sep = "="),
        other_args
    )
    run_cellranger(
        subcmd = "count", cellranger_args = args,
        cellranger_cmd = cellranger_cmd, sys_args = sys_args
    )
}

#' Run cellranger
#' @param subcmd Sub-Command of cellranger.
#' @param cellranger_args arguments passed to cellranger sub-command.
#' @param cellranger_cmd Path to cellranger cmd.
#' @param sys_args Other arguments passed to [system2].
#' @export 
run_cellranger <- function(subcmd, cellranger_args = character(), cellranger_cmd = NULL, sys_args = list()) {
    run_command(
        args = c(subcmd, cellranger_args),
        cmd = cellranger_cmd,
        name = "cellranger",
        sys_args = sys_args
    )
}

run_command <- function(args = character(), cmd, name = NULL, sys_args = list(), verbose = TRUE) {
    if (!is.null(cmd)) {
        if (!file.exists(cmd) &&
            nzchar(Sys.which(cmd)) > 0L) {
            cmd <- Sys.which(cmd)
        }
        cmd <- normalizePath(cmd, mustWork = TRUE)
    } else if (!is.null(name)) {
        cmd <- Sys.which(name)
        if (nzchar(cmd) == 0L) {
            cli::cli_abort("Cannot find {.field {name}} command")
        }
    }
    sys_args <- c(list(command = cmd, args = args), sys_args)
    if (verbose) {
        cli_args <- cli::cli_vec( # nolint
            args,
            list("vec-sep" = " ", "vec-last" = " ")
        )
        cli::cli_alert("Running command {.field {cmd} {cli_args}}")
    }
    do.call(system2, sys_args)
}

handle_arg <- function(arg, name, format = "%s", sep = " ") {
    if (is.null(arg) || isFALSE(arg)) {
        return(NULL)
    } else if (isTRUE(arg)) {
        return(name)
    } else {
        return(sprintf(paste(name, format, sep = sep), arg))
    }
}
