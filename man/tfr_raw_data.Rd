\name{tfr_raw_data}
\alias{TFR_cleaned_2019}
\docType{data}
\title{
Raw TFR Data
}
\description{
Data set containing the raw TFR estimation of all past TFR estimation and the data quality covariates.
}
\usage{data("TFR_cleaned_2019")}
\format{
  A data frame with 12709 observations on the following 5 variables.
  \describe{
    \item{\code{ISO.code}}{Three-digit UN ISO-3166 code for the country of that observation is for.}
    \item{\code{Year}}{a numeric vector for the year of the observation data.}
    \item{\code{DataValue}}{TFR value.}
    \item{\code{Estimating.Methods}}{Estimating methods to obtain this value. One of the categorical data quality indicator.}
    \item{\code{DataProcess}}{Source of the data. One of the categorical data quality indicator.}
  }
}
\details{
It works as the default raw TFR data in a simulation when \code{my.tfr.raw.file} is \code{NULL}. If then, \code{DataProcess} and \code{Estimating.Methods} will work as default data quality covariates.
}
\source{
Data provided by the United Nations Population Division.
}
\examples{
data(TFR_cleaned_2019)
head(TFR_cleaned_2019)
}
\keyword{datasets}
