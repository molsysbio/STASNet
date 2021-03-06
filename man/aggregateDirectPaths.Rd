% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameters_comparison.R
\name{aggregateDirectPaths}
\alias{aggregateDirectPaths}
\alias{subsetDirectPathsAggregation}
\title{Aggregate the paths values for plotting}
\usage{
aggregateDirectPaths(direct_paths)

subsetDirectPathsAggregation(aggregated_paths, regex, not = FALSE)
}
\arguments{
\item{direct_paths}{A list of direct paths (as outputed by getDirectPaths). The names of the elements are considered to be the names of the models.}

\item{aggregated_paths}{Aggregated paths as returned by aggregateDirectPaths}

\item{regex}{A vector of regular expression}

\item{not}{Whether the subpaths should not match the regex}
}
\value{
A list with two entries. 'paths' is a matrix with the path and the models names as rownames, and columns 'value', 'hv' and 'lv' representing the value with confidence interval of the path. 'model_names' is the name of the models.
}
\description{
Extract the subset of paths aggregated from multiple models that match a gived regex
}
