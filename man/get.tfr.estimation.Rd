\name{get.tfr.estimation}
\alias{get.tfr.estimation}
\title{
Get Past TFR Estimation
}
\description{
Get past TFR estimation, including trajectories and quantiles if required.
}
\usage{
get.tfr.estimation(mcmc.list=NULL, country.code=NULL, ISO.code=NULL, sim.dir=NULL, burnin=0, thin = 1, probs=NULL)
}
\arguments{
  \item{mcmc.list}{Object of class \code{\link{bayesTFR.mcmc.set}} corresponding Phase II MCMCs. If it is \code{NULL}, the object is loaded from the directory given by \code{sim.dir}.}
  \item{country.code}{Name or numerical code of a country. If it is \code{NULL}, \code{ISO.code} need to be provided with ISO-3 character.}
  \item{ISO.code}{Country code as ISO-3 character. Only used if \code{country.code} is set to be \code{NULL}.}
  \item{sim.dir}{Directory with the MCMC simulation results.}
  \item{burnin}{Burn-in for getting trajectories and quantiles. A positive burn-in \eqn{x} will throw first \eqn{x} iterations of estimations for each chain.}
  \item{thin}{Thin for getting trajectories and quantiles. Thinning level \eqn{x} greater than 1 will store one iteration per \eqn{x} samples}
  \item{probs}{A vector of numbers between \eqn{[0,1]}. If it is set to be \code{NULL} then this step will be skipped.  Otherwise, quantiles of estimation for each year will be outputted.}
}
\details{
This function is used to obtain the TFR trajectories as well as corresponding quantiles if the \code{mcmc.list} is obtained when \code{uncertainty=TRUE}. This is used to visualize past estimates as well as prediction. Quantiles could also be obtained if \code{probs} is not \code{NULL}.
}
\value{
  \item{tfr_table}{Table storing the trajectories. The table will be a matrix with rows equal to number of trajectories, and column equal to number of periods.}
  \item{country.obj}{A list storing the information of the country for the trajectory table is for. Include \code{name} for the name of the country, \code{code} for the ISO-3166 standard 3-digit code, and \code{index} for the number of column this country stays in the meta file of the \code{mcmc.set}.}
  \item{tfr_quantile}{Optional. A data.table object, storing the quantiles of estimates of each period required by \code{probs}. Another column \code{year} for marking the period.}
}
\author{
Peiran Liu, Hana Sevcikova
}
\examples{
\dontrun{
mcmc.set <- run.tfr.mcmc(nr.chains=1, iter=10, output.dir='test.dir', replace.output=TRUE, uncertainty=TRUE)
# add the median of the hierarchical DL curves
res <- get.tfr.estimation(mcmc.set, 'Nigeria', probs=0.5)}

}
