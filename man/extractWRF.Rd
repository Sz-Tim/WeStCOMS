% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadAndExtract.R
\name{extractWRF}
\alias{extractWRF}
\title{Extract variables from WeStCOMS-WRF (Weather Research Forecasting) model}
\usage{
extractWRF(sampling.df, wrf.dir, wrf_i, returnFullDf = FALSE)
}
\arguments{
\item{sampling.df}{Dataframe with a row for each sample. Columns must include
`site.id`, `obs.id`, `date` (YYYYMMDD), `wrf_i`, `lon`, `lat`, where
`wrf_i` indicates the corresponding rownumber in the dataframe wrf_i}

\item{wrf.dir}{Character vector with directory for WRF files}

\item{wrf_i}{Dataframe with WRF filenames and date ranges in yyyy-mm-dd
(required column names: fname, date_0, date_1)}

\item{returnFullDf}{Logical: return full dataframe or only obs.id +
hydrovars?}
}
\value{
dataframe with site.id, date, and WRF variables
}
\description{
Given an input dataframe with site locations and dates, this function
extracts the corresponding WeStCOMS-WRF data. Output is summarised for the
day as given in `daySummaryFn`.
}
