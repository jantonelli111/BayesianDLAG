% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SummaryFunctions.R
\name{PlotMainDLAG}
\alias{PlotMainDLAG}
\title{Plot the posterior mean and 95 percent credible interval for a main effect distributed lag}
\usage{
PlotMainDLAG(DLAGmodel, whichExposure, main = NULL, ylim = NULL)
}
\arguments{
\item{DLAGmodel}{A fitted model from the DLAGbayes function}

\item{whichExposure}{an integer indicating which exposure is to be plotted.
Must be between 1 and the number of exposures}

\item{main}{The title of the plot}

\item{ylim}{The limits of the y-axis of the plot}
}
\description{
Plot the posterior mean and 95 percent credible interval for a main effect distributed lag
}
