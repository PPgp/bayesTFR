\name{get.tfr.estimation}
\alias{get.tfr.estimation}
\title{
Get Past TFR Estimation
}
\description{
Get past TFR estimation, including trajectories and quantiles if required.
}
\usage{
get.tfr.estimation(mcmc.list = NULL, country.code = NULL, ISO.code = NULL, 
    sim.dir = NULL, burnin = 0, thin = 1, probs = NULL, adjust = TRUE)
}
\arguments{
  \item{mcmc.list}{Object of class \code{\link{bayesTFR.mcmc.set}} corresponding Phase II MCMCs. If it is \code{NULL}, the object is loaded from the directory given by \code{sim.dir}.}
  \item{country.code}{Name or numerical code of a country. If it is \code{NULL}, \code{ISO.code} needs to be provided as ISO-3 character.}
  \item{ISO.code}{Country code as ISO-3 character. Only used if \code{country.code} is set to be \code{NULL}.}
  \item{sim.dir}{Directory with the MCMC simulation results. Only used if \code{mcmc.list} is \code{NULL}.}
  \item{burnin}{Burn-in for getting trajectories and quantiles. A positive burn-in \eqn{x} will remove first \eqn{x} iterations from each chain.}
  \item{thin}{Thin for getting trajectories and quantiles. Thinning level \eqn{x} greater than 1 will store one iteration per \eqn{x} samples.}
  \item{probs}{A vector of numbers between \eqn{[0,1]} specifying which estimation quantiles should be outputted. If it is set to \code{NULL} no quantiles are returned.}
  \item{adjust}{Logical indicating whether the adjusted median and trajectories should be returned.}
}
\details{
This function is used to obtain the TFR estimation trajectories as well as corresponding quantiles if the \code{mcmc.list} has been obtained while taking account for uncertainty about the past, i.e. \code{uncertainty=TRUE} in \code{\link{run.tfr.mcmc}}. Quantiles are included in the results if \code{probs} is not \code{NULL}.
}
\value{
  \item{tfr_table}{Table storing the trajectories. It is a matrix with rows equal to number of trajectories, and column equal to number of time periods.}
  \item{country.obj}{A list storing information about the country which the trajectories and quantiles correspond to. It corresponds to the output of \code{\link{get.country.object}}.}
  \item{tfr_quantile}{Optional. A data.table object, storing the quantiles of estimates for each time period as specified by the \code{probs} argument. The time periods are contained in the \code{year} column.}
}
\author{
Peiran Liu, Hana Sevcikova
}
\examples{
\dontrun{
sim.dir <- tempfile()
m <- run.tfr.mcmc(nr.chains = 1, iter = 10, output.dir = sim.dir, uncertainty = TRUE)
get.tfr.estimation(m, "Nigeria", probs = c(0.1, 0.5, 0.9))
unlink(sim.dir, recursive = TRUE)}
}
