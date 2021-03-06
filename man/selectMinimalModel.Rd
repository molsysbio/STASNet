% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_manipulation.R
\name{reduceModel}
\alias{reduceModel}
\alias{selectMinimalModel}
\title{Selection of a minimal model by iteratively removing insignificant links using the likelihood ratio test.}
\usage{
reduceModel(original_model, accuracy = 0.95)

selectMinimalModel(original_model, accuracy = 0.95, verbose = F)
}
\arguments{
\item{original_model}{An MRAmodel object, as the one produced by createModel or importModel}

\item{accuracy}{Probability threshold, the type I error for each link will be 1-accuracy. Multiple testing is not taken into account.}
}
\value{
An MRAmodel object of the reduced model with the data
}
\description{
Selection of a minimal model by iteratively removing insignificant links using the likelihood ratio test.
}
\author{
Mathurin Dorel \email{dorel@horus.ens.fr}

Bertram Klinger \email{bertram.klinger@charite.de}
}
