% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellranger.R
\name{run_cellranger}
\alias{run_cellranger}
\title{Run cellranger}
\usage{
run_cellranger(
  subcmd,
  cellranger_args = character(),
  cellranger_cmd = NULL,
  sys_args = list()
)
}
\arguments{
\item{subcmd}{Sub-Command of cellranger.}

\item{cellranger_args}{arguments passed to cellranger sub-command.}

\item{cellranger_cmd}{Path to cellranger cmd.}

\item{sys_args}{Other arguments passed to \link{system2}.}
}
\description{
Run cellranger
}
