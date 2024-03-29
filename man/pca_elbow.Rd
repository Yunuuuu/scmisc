% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elbow.R
\name{pca_elbow}
\alias{pca_elbow}
\alias{pca_elbow,numeric-method}
\alias{pca_elbow,SingleCellExperiment-method}
\title{Find the elbow point in the curve of variance explained by each successive PC}
\usage{
pca_elbow(x, ...)

\S4method{pca_elbow}{numeric}(x)

\S4method{pca_elbow}{SingleCellExperiment}(x, dimred = "PCA")
}
\arguments{
\item{x}{A Numeric vector containing the variance explained by each PC.}

\item{...}{Other arguments passed to specific methods.}

\item{dimred}{A string or integer scalar indicating the reduced dimension
containing the \code{PCA} with a \code{percentVar} attribute.}
}
\value{
An integer scalar specifying the number of PCs at the elbow point.
}
\description{
Find the elbow point in the curve of variance explained by each successive
PC. This can be used to determine the number of PCs to retain.
}
