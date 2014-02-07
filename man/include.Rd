\name{include}
\alias{include}
\alias{include_2008}
\alias{include_2010}
\alias{include_2012}

\docType{data}
\title{
Inclusion Codes
}
\description{
Data sets containing codes that determine which countries are to be included into a simulation or/and projections. 
}
\usage{
data(include_2008)
data(include_2010)
data(include_2012)
}

\format{
  Data frames containing one record per country or region. It has the following variables:
  \describe{
    \item{\code{country}}{Name of country or region. Not used.}
    \item{\code{country_code}}{Numerical Location Code (3-digit codes following ISO 3166-1 numeric standard) - see \url{http://en.wikipedia.org/wiki/ISO_3166-1_numeric}.}
    \item{include_code}{Entries for which \code{include_code=2} are included in MCMC simulations (i.e. estimation of the model parameters). Entries for which \code{include_code} is 1 or 2 are included in the prediction.}}
}
\details{
	In a simulation, an \code{include_*} dataset is selected that corresponds to the given \code{wpp.year} passed to the function \code{\link{run.tfr.mcmc}}. It is merged with a \code{\link[wpp2012:tfr]{tfr}} dataset from the corresponding wpp package using the \code{country_code} column. Thus, the country entries in this dataset should correspond to entries in the \code{tfr} dataset.

The package contains also a dataset called \file{my_tfr_template} (in \file{extdata} directory) which is a template for user-specified TFR time series. It has the same structure as the \code{\link[wpp2012:tfr]{tfr}} dataset, except that most of the columns are optional. The only required column is \code{country_code} (see description of the argument \code{my.tfr.file} in \code{\link{run.tfr.mcmc}}).
}
\source{
Data provided by the United Nations Population Division.
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(include_2012)
head(include_2012)
}
\keyword{datasets}