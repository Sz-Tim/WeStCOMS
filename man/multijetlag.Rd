% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{multijetlag}
\alias{multijetlag}
\title{Lag multiple variables at once}
\usage{
multijetlag(data, ..., n = 10)
}
\arguments{
\item{data}{Dataframe}

\item{...}{Unquoted variable names to lag}

\item{n}{Number of lags}
}
\description{
From https://stackoverflow.com/questions/55814028/multiple-lags-with-dplyr
}
