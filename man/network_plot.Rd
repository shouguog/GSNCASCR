% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_2.R
\name{network_plot}
\alias{network_plot}
\title{Do permutation for calculate score}
\usage{
network_plot(
  networklist,
  name1 = "status1",
  name2 = "status2",
  fileName = NULL
)
}
\arguments{
\item{networklist}{list with two networks}

\item{name1}{the first status for analysis}

\item{name2}{the second status for analysis}

\item{filename}{png file name to save the networks, no figure generated when fileName is NULL}
}
\value{
gglist with three ggplot2 objects
}
\description{
Do permutation for calculate score
}
