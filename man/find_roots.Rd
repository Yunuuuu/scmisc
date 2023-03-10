% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_roots.R
\name{find_roots}
\alias{find_roots}
\title{Find roots in a DiffusionMap object}
\usage{
find_roots(dm, start, ends = NULL, ref, n_root = 100L)
}
\arguments{
\item{dm}{A \link[destiny:DiffusionMap-class]{DiffusionMap} object.}

\item{start, ends}{The start and ends cluster identity. start must have a
length 1L, while the length of ends ranges from 1L to 2L. All start and
ends must exist in \code{ref}. If ends is \code{NULL}, this only check if one tip is
in the start cluster.}

\item{ref}{All cell identity of \code{dm}, this must have the same length of
\code{dm@d}.}

\item{n_root}{The number of nodes we should test.}
}
\value{
An integer index .
}
\description{
Find roots in a DiffusionMap object. If we know the \code{start} and \code{end} of a
root, we then check if the tips of the root are in the start and ends node we
already known (see \link[destiny:find_tips]{find_tips}). This function in this
way check \code{n_root} nodes with the top \code{DPT} in the start cluster and return
roots meet this criteria.
}
