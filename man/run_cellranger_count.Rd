% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellranger.R
\name{run_cellranger_count}
\alias{run_cellranger_count}
\title{Run cellranger-count}
\usage{
run_cellranger_count(
  id,
  fastqs,
  transcriptome = NULL,
  sample = NULL,
  cores = NULL,
  other_args = c("--nosecondary"),
  cellranger_cmd = NULL,
  sys_args = list()
)
}
\arguments{
\item{id}{A unique run id and output folder name.}

\item{fastqs}{Path to input FASTQ data.}

\item{transcriptome}{Path of folder containing 10x-compatible transcriptome
reference.}

\item{sample}{Prefix of the filenames of FASTQs to select.}

\item{cores}{Set max cores the pipeline may request at one time. Only applies
to local jobs.}

\item{other_args}{Other arguments passed to cellranger count.}

\item{cellranger_cmd}{Path to cellranger cmd.}

\item{sys_args}{Other arguments passed to \link{system2}.}
}
\description{
Count gene expression (targeted or whole-transcriptome) and/or feature barcode reads from a single sample and GEM well
}
