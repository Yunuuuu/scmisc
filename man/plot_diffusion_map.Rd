% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_diffusion_map.R
\name{plot_diffusion_map}
\alias{plot_diffusion_map}
\alias{plot_diffusion_map,SingleCellExperiment,DiffusionMap-method}
\alias{plot_diffusion_map,SingleCellExperiment,DPT-method}
\alias{add_dpt_paths}
\alias{add_dpt_paths,SingleCellExperiment,DPT-method}
\title{Plot the diffusion map and the DPT paths}
\usage{
plot_diffusion_map(x, y, ...)

\S4method{plot_diffusion_map}{SingleCellExperiment,DiffusionMap}(
  x,
  y,
  dimred = "DiffusionMap",
  colour_by = "DPT",
  root = NULL,
  ...
)

\S4method{plot_diffusion_map}{SingleCellExperiment,DPT}(
  x,
  y,
  dimred = "DiffusionMap",
  colour_by = "DPT",
  root = NULL,
  add_paths = TRUE,
  paths_to = NULL,
  path_args = list(),
  w_width = 0.1,
  ...,
  ncomponents = 2L
)

add_dpt_paths(x, y, ...)

\S4method{add_dpt_paths}{SingleCellExperiment,DPT}(
  x,
  y,
  root = NULL,
  paths_to = NULL,
  colour = c("gray62", "#F5C710", "#CD0BBC", "#28E2E5", "#2297E6", "#61D04F", "#DF536B",
    "black"),
  color = colour,
  w_width = 0.1,
  ...,
  ncomponents = 2L
)
}
\arguments{
\item{x}{A \link[SingleCellExperiment:SingleCellExperiment]{SingleCellExperiment}
object.}

\item{y}{A \link[destiny:DiffusionMap-class]{DiffusionMap} or \link[destiny:DPT]{DPT}
object.}

\item{...}{For \code{plot_diffusion_map}, other arguments passed to
\link[scater:plotReducedDim]{plotReducedDim}; for \code{add_dpt_paths}, other
arguments passed to \link[ggplot2:geom_path]{geom_path}.}

\item{dimred}{A string or integer scalar indicating the reduced dimension
result in \code{reducedDims(x)} to plot. If "DiffusionMap", the \code{eigenvectors}
of \code{y} will be ploted.}

\item{colour_by}{Specification of a column metadata field or a feature to
colour by, see the by argument in
\code{\link[scater:retrieveCellInfo]{?retrieveCellInfo}} for possible values. If y
is a \link[destiny:DPT]{DPT} object, this can also be "DPT" (or starting by
"DPT", such as "DPT1", "DPT2" ...), "branch" or "Branch" (only for \code{DPT} of
y method) e.g. which will be extracted from \code{y}. If See also
\link[destiny:plot.DPT]{plot.DPT}.}

\item{root}{Root cell index. Will be used as the start of the \code{DPT} (default:
the first tip index for \code{DPT} method or \link[destiny:random_root]{random_root}
for \code{DiffusionMap} method).}

\item{add_paths}{A scalar logical value indicates whether to add path.}

\item{paths_to}{Numeric Branch IDs, which are used as target(s) for the
path(s) to draw.}

\item{path_args}{A list of arguments to passed to
\link[ggplot2:geom_path]{geom_path}. Elements of the list can be the same
length of \code{paths_to}. If you want to pass a complicated oject into
\code{path_args}, try to wrap it with \link[base:list]{list} function, since the
internal use \code{rep_len} parallelly to get the same length of \code{paths_to}. for
example \code{list(grid::arrow(ends = "last"))}.}

\item{w_width}{Window width for smoothing the path (see
\link[smoother:smth.gaussian]{smth.gaussian}).}

\item{ncomponents}{A numeric scalar indicating the number of dimensions to
plot, starting from the first dimension. Alternatively, a numeric vector
specifying the dimensions to be plotted.}

\item{colour}{The color of the paths.}

\item{color}{The color of the paths.}
}
\value{
For \code{plot_diffusion_map}, A ggplot2 object depicting the Diffusion
map, if add_paths is \code{TRUE}, DPT path of the root is also added in this
map. for \code{add_dpt_paths}, a ggplot2 geom object is returned.
}
\description{
Plot the diffusion map and the DPT paths
}
