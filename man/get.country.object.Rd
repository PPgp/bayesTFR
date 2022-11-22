\name{get.country.object}
\alias{get.country.object}
\alias{get.countries.table}
\alias{get.countries.table.bayesTFR.mcmc.set}
\alias{get.countries.table.bayesTFR.prediction}
\alias{get.countries.phase}
\alias{get.countries.phase.bayesTFR.mcmc.set}
\alias{get.countries.phase.bayesTFR.prediction}


\title{
Accessing Country Information
}
\description{
Function \code{get.country.object} returns an object containing country name, code and index. Functions \code{get.countries.table} return a data frame containing codes, names and optionally ISO character codes of all countries. Functions \code{get.countries.phase} return countries table with the TFR phase they are currently in (1, 2, or 3).
}
\usage{
get.country.object(country, meta = NULL, country.table = NULL, index = FALSE)

\method{get.countries.table}{bayesTFR.mcmc.set}(object, iso = FALSE, \dots)
\method{get.countries.table}{bayesTFR.prediction}(object, iso = FALSE, \dots)

\method{get.countries.phase}{bayesTFR.mcmc.set}(object, \dots)
\method{get.countries.phase}{bayesTFR.prediction}(object, \dots)
}
\arguments{
  \item{country}{Country name, code or index. If it is an index, the argument \code{index} must be set to \code{TRUE}. The code can be either numeric or ISO-2 or ISO-3 characters.}
  \item{meta}{Object of class \code{\link{bayesTFR.mcmc.meta}}. If it is not given, the argument \code{country.table} must be given.}
  \item{country.table}{A table containing columns \dQuote{name} and \dQuote{code} from which the country info can be extracted. Only relevant, if \code{meta} is \code{NULL}. }
  \item{index}{Logical determining if the argument \code{country} is an index.}
  \item{object}{Object of class \code{\link{bayesTFR.mcmc.set}} or \code{\link{bayesTFR.prediction}}.}
  \item{iso}{Logical. If \code{TRUE}, two extra columns are added to the table, namely 2- and 3-characters ISO codes.}
  \item{\dots}{Not used.}
}
\details{
Given partial information about a country (i.e. having either name or code or index), the function \code{get.country.object} returns  an object containing all three pieces of information. Only countries are considered that are included in the simulations and predictions. Country index is an internal index used in various components of a \code{\link{bayesTFR.mcmc.meta}} object.
}
\value{
Function \code{get.country.object} returns a list with components:
\item{name}{Country name}
\item{code}{Country numeric code}
\item{index}{Country index}

Function \code{get.countries.table} returns a data frame with columns \code{code}, \code{name}, and optionally (if \code{iso} is \code{TRUE}) \code{iso2} and \code{iso3}.

Function \code{get.countries.phase} returns a data frame with columns \code{code}, \code{name} and \code{phase}.

}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{country.names}}
}
\examples{
sim.dir <- file.path(find.package("bayesTFR"), "ex-data", "bayesTFR.output")
m <- get.tfr.mcmc(sim.dir)
# all five calls should give the same answer
get.country.object('China', m$meta)
get.country.object('CN', m$meta)
get.country.object(156, m$meta)
get.country.object(56, m$meta, index=TRUE)
get.country.object(156, NULL, country.table=get.countries.table(m))

# phase 3 countries
subset(get.countries.phase(m), phase == 3)
}

\keyword{ attribute }

