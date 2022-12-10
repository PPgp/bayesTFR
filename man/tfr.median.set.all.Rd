\name{tfr.median.set.all}
\alias{tfr.median.set.all}
\alias{tfr.median.reset.estimation}
\alias{tfr.shift.estimation.to.wpp}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Editing median for estimation and projections.
}
\description{
These functions are to be used by expert analysts. They allow to change the estimation and projection medians to specific values or to the WPP values.
}
\usage{
tfr.median.set.all(sim.dir, country, values, years = NULL, 
    burnin = 0, thin = 1)

tfr.median.reset.estimation(sim.dir, countries)

tfr.shift.estimation.to.wpp(sim.dir, \dots, verbose = TRUE)
}

\arguments{
  \item{sim.dir}{Directory containing the prediction object.}
  \item{country}{Name or numerical code of a country.}
  \item{countries}{Vector of country names or codes.}
  \item{values}{Array of the new median values.}
  \item{years}{Numeric vector giving years which \code{values} correspond to. Ideally it should be of the same length as \code{values}.}
  \item{burnin}{Burnin to use when computing the estimation median.}
  \item{thin}{Thinning interval to use when computing the estimation median.}
  \item{\dots}{Can be used to pass \code{burnin} \code{thin} to the underlying funcions.}
  \item{verbose}{Logical. If \code{TRUE} a progress of the adjustment is shown.}
}
\details{
Expert analysts can use these functions to adjust both prediction and estimation medians. Estimation medians can only be adjusted if the simulation was performed with \code{uncertainty = TRUE}. In such a case \code{years} can include past time periods. By default a union of estimation and projection time periods is considered when matched to \code{values}. 

Function \code{tfr.shift.estimation.to.wpp} shifts the median estimation of all countries so that they match the values in the \code{tfr} dataset of the corresponding WPP package. Argument \code{burnin} and \code{thin} should be passed to compute the estimation medians.

Function \code{tfr.median.reset.estimation} resets previous adjustments obtained using \code{tfr.median.set.all}.
}
\value{
Output is a list. If there are time periods matched to estimation, an object of class \code{\link{bayesTFR.mcmc.meta}} is included in the element \code{meta}. If there are time periods matched to years in prediction, then an object of class \code{\link{bayesTFR.prediction}} is included in the element \code{pred}.

Function \code{tfr.shift.estimation.to.wpp} returns the updated \code{mcmc.set} object.
}
\author{
Peiran Liu
}

\seealso{\code{\link{tfr.shift.prediction.to.wpp}} for shifting prediction medians to WPP values.}
