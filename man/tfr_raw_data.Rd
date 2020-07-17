\name{tfr_raw_data}
\alias{rawTFR}
\docType{data}
\title{
Raw TFR Data
}
\description{
Data set containing the raw TFR estimation of all past TFR estimation and the data quality covariates.
}
\usage{data("rawTFR")}
\format{
  A data frame with 12709 observations on the following 5 variables.
  \describe{
    \item{\code{country_code}}{Three-digit UN ISO-3166 code for the country of that observation is for.}
    \item{\code{year}}{a numeric vector for the year of the observation data.}
    \item{\code{tfr}}{TFR value.}
    \item{\code{method}}{Estimating methods to obtain this value. One of the categorical data quality indicator.}
    \item{\code{source}}{Source of the data. One of the categorical data quality indicator.}
  }
}
\details{
It works as the default raw TFR data in a simulation when \code{my.tfr.raw.file} is \code{NULL}. If then, \code{DataProcess} and \code{Estimating.Methods} will work as default data quality covariates.
}
\source{
Data provided by the United Nations Population Division.
}
\examples{
data(rawTFR)
head(rawTFR)
}
\keyword{datasets}
