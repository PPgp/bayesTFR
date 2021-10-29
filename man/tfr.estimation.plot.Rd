\name{tfr.estimation.plot}
\alias{tfr.estimation.plot}
\title{
Plot TFR Estimation
}
\description{
Plot past TFR estimation results from a simulation that accounted for past TFR uncertainty.
}
\usage{
tfr.estimation.plot(mcmc.list = NULL, country = NULL, sim.dir = NULL, 
    burnin = 0, thin = 1, pis = c(80, 95), plot.raw = TRUE, 
    grouping = "source", save.image = TRUE, plot.dir = "Estimation.plot", 
    adjust = TRUE, country.code = deprecated(), ISO.code = deprecated())
}
\arguments{
  \item{mcmc.list}{Object of class \code{\link{bayesTFR.mcmc.set}} corresponding Phase II MCMCs. If it is \code{NULL}, the object is loaded from the directory given by \code{sim.dir}.}
  \item{country}{Name or numerical code of a country. It can also be given as ISO-2 or ISO-3 characters.}
  \item{sim.dir}{Directory with the MCMC simulation results.}
  \item{burnin}{Burn-in for getting trajectories and quantiles. A positive burn-in \eqn{x} will remove first \eqn{x} iterations from each chain.}
  \item{thin}{Thin for getting trajectories and quantiles. Thinning level \eqn{x} greater than 1 will store one iteration per \eqn{x} samples}
  \item{pis}{Probability interval. It can be a single number or an array of two numbers.}
  \item{plot.raw}{Whether raw data used for the estimation should be plotted.}
  \item{grouping}{If raw data is plotted, then grouping should be one of the categorical feature in the data, so that the color and shape of the raw data will differ for different groups.}
  \item{save.image}{Logical. Whether the resulting plot will be saved.}
  \item{plot.dir}{If \code{save.image=TRUE}, specify the directory for saving the plot.}
  \item{adjust}{Logical. By default, if the estimation median is adjusted using e.g. \code{\link{tfr.median.set.all}}, the function plots the adjusted median. If \code{adjust=FALSE} the original (non-adjusted) median is plotted.}
  \item{country.code, ISO.code}{Deprecated arguments. Use argument \code{country} instead.}
}

\details{
\code{tfr.estimation.plot} plots posterior distribution of past TFR estimations for a given country. It only works if uncertainty is considered in the MCMC process.
}
\author{
Peiran Liu, Hana Sevcikova
}

\examples{
\dontrun{
sim.dir <- tempfile()
mcmc.set <- run.tfr.mcmc(nr.chains = 1, iter = 10, output.dir = sim.dir, 
    replace.output = TRUE, uncertainty = TRUE)
tfr.estimation.plot(mcmc.set, "Nigeria", save.image = FALSE)
unlink(sim.dir, recursive = TRUE)}
}
