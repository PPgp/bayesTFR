\name{tfr.median.set.all}
\alias{tfr.median.set.all}
\alias{tfr.median.reset.estimation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Editing median for estimation and projections.
}
\description{
These functions are to be used by expert analysts. They allow to change the estimation and projection medians to specific values.
}
\usage{
tfr.median.set.all(sim.dir, country, values, years=NULL)

tfr.median.reset.estimation(sim.dir, countries)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sim.dir}{Directory containing the prediction object.}
  \item{country}{Name or numerical code of a country.}
  \item{countries}{Vector of country names or codes.}
  \item{values}{Array of the new median values.}
  \item{years}{Numeric vector giving years which \code{values} correspond to. Ideally it should be of the same length as \code{values}.}
}
\details{
These functions are set for experts to adjust both prediction and estimation medians.
}
\value{
Output is a list. If there are years in estimation, then there will be an object of class \code{\link{bayesTFR.mcmc.meta}}. If there are years in prediction, then there will be an object of class \code{\link{bayesTFR.prediction}}
}
\author{
Peiran Liu
}
