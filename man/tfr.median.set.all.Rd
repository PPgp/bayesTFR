\name{tfr.median.set.all}
\alias{tfr.median.set.all}
\alias{tfr.mean.set.all}
\alias{tfr.median.reset.estimation}
\alias{tfr.shift.reset.estimation}
\alias{tfr.shift.estimation.to.wpp}

\title{
Shifting the Distribution of the Estimation and Projection.
}
\description{
These functions are to be used by expert analysts. They allow to change the estimation and projection medians or means to specific values or to the WPP values.
}
\usage{
tfr.median.set.all(sim.dir, country, values, years = NULL, 
    burnin = 0, thin = 1, subdir = "predictions")

tfr.mean.set.all(sim.dir, country, values, years = NULL, 
    burnin = 0, thin = 1, subdir = "predictions")

tfr.shift.reset.estimation(sim.dir, countries = NULL)

tfr.median.reset.estimation(...) # obsolete

tfr.shift.estimation.to.wpp(sim.dir, \dots, stat = "median", 
    verbose = TRUE)
}

\arguments{
  \item{sim.dir}{Directory containing the prediction object.}
  \item{country}{Name or numerical code of a country.}
  \item{countries}{Vector of country names or codes. If \code{NULL}, the reset is done for all countries.}
  \item{values}{Array of the new median values.}
  \item{years}{Numeric vector giving years which \code{values} correspond to. Ideally it should be of the same length as \code{values}.}
  \item{burnin}{Burnin to use when computing the estimation median.}
  \item{thin}{Thinning interval to use when computing the estimation median.}
  \item{subdir}{Subdirectory of \code{sim.dir} containing the predictions.}
  \item{\dots}{Can be used to pass \code{burnin} \code{thin} to the underlying funcions.}
  \item{stat}{Statistics to adjust to. Possible values are \dQuote{median} (default) and \dQuote{mean}.}
  \item{verbose}{Logical. If \code{TRUE} a progress of the adjustment is shown.}
}
\details{
Expert analysts can use these functions to adjust both prediction and estimation medians or means. The estimation distribution can only be adjusted if the simulation was performed with \code{uncertainty = TRUE}. In such a case \code{years} can include past time periods. By default a union of estimation and projection time periods is considered when matched to \code{values}. 

Functions \code{tfr.median.set.all} and \code{tfr.mean.set.all} move the estimation and projection distribution of given country so that their medians/means match the given values. Arguments \code{burnin} and \code{thin} should be passed to compute the estimation statistics.

Function \code{tfr.shift.estimation.to.wpp} shifts the medians or means (if \code{stat = "mean"}) of the estimation of all countries so that they match the values in the \code{tfr} dataset of the corresponding WPP package. 

Function \code{tfr.shift.reset.estimation} resets previous adjustments obtained using the functions above. By default it resets adjustments for all countries. Note that it resets both, estimation as well as projection. Function \code{tfr.median.reset.estimation} is obsolete, replaced by \code{tfr.shift.reset.estimation}.

The offset of the estimation distribution is stored in the list item \code{traj.shift.estimation} of the meta object stored in the simulation directory. When these functions are used, the meta object is updated and written to disk. Similarly, if projection is adjusted via these functions, the updated item \code{traj.shift} of the prediction object is written to the prediction directory.
}
\value{
List with two elements:
\item{meta}{If there are time periods matched to estimation, this is an object of class \code{\link{bayesTFR.mcmc.meta}} containing the updated element \code{traj.shift.estimation}.} 
\item{pred}{If there are time periods matched to years in prediction, this is an object of class \code{\link{bayesTFR.prediction}} with the updated element \code{traj.shift}.}

Function \code{tfr.shift.estimation.to.wpp} returns the updated \code{mcmc.set} object.
}
\author{
Peiran Liu, Hana Sevcikova
}

\seealso{\code{\link{tfr.shift.prediction.to.wpp}} for shifting prediction medians/means to WPP values.}
