% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot.disag_model}
\alias{plot.disag_model}
\title{Plot results of fitted model}
\usage{
\method{plot}{disag_model}(x, ...)
}
\arguments{
\item{x}{Object of class \emph{disag_model} to be plotted.}

\item{...}{Further arguments to \emph{plot} function.}
}
\value{
A list of two ggplot plots: results of the fixed effects and an in-sample observed vs predicted plot
}
\description{
Plotting function for class \emph{disag_model} (the result of the disaggregation fitting).
}
\details{
Produces two plots: results of the fixed effects and in-sample observed vs predicted plot.
}
