% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot.SANvb}
\alias{plot.SANvb}
\title{Plotting the variational inference output}
\usage{
\method{plot}{SANvb}(x, ...)
}
\arguments{
\item{x}{object of class \code{SANvb} (the result of a call to \code{\link{variational_CAM}}, 
\code{\link{variational_fiSAN}}, \code{\link{variational_fSAN}}.}

\item{...}{additional graphical parameters to be passed.}
}
\value{
The function plots a summary of the fitted model.
}
\description{
Plot method for objects of class \code{SANvb}. 
The function displays two graphs, meant to analyze the estimated distributional and observational clusters.
}
\seealso{
\code{\link{print.SANvb}}, \code{\link{variational_CAM}}, 
\code{\link{variational_fiSAN}}, \code{\link{variational_fSAN}}.
}
