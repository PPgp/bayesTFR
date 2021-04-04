\name{tfr_raw_data}
\alias{rawTFR}
\docType{data}
\title{
Raw TFR Data
}
\description{
Data set containing the raw TFR estimates for all countries and the data quality covariates.
}
\usage{data("rawTFR")}
\format{
  A data frame with 12709 observations on the following 5 variables.
  \describe{
    \item{\code{country_code}}{Three-digit UN ISO-3166 code for the country of that observation is for.}
    \item{\code{year}}{a numeric vector for the year of the observation data.}
    \item{\code{tfr}}{TFR value.}
    \item{\code{method}}{Estimation method to obtain this value. One of the categorical data quality indicator.}
    \item{\code{source}}{Source of the data. One of the categorical data quality indicator.}
  }
}
\details{
It is used as the default raw TFR data in a \code{\link{run.tfr.mcmc}} simulation. It can be used as a template for a user-defined data which can be provided via the \code{my.tfr.raw.file} argument of \code{\link{run.tfr.mcmc}}. The \dQuote{method} and \dQuote{source} columns are used as the default data quality covariates.
}
\source{
Data provided by the United Nations Population Division.
}
\examples{
data(rawTFR)
head(rawTFR)
}
\keyword{datasets}
