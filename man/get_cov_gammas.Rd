\name{get.cov.gammas}
\Rdversion{1.1}
\alias{get.cov.gammas}
\title{Covariance Matrices of Gamma Parameters
}
\description{
From a given MCMC, obtain a covariance matrix of the \eqn{\gamma_{ci}}{gamma_ci} (\eqn{i=1,2,3}) parameters for each country \eqn{c}.
}
\usage{
get.cov.gammas(mcmc.set = NULL, sim.dir = NULL, burnin = 200, chain.id = 1)
}
\arguments{
  \item{mcmc.set}{Object of class \code{bayesTFR.mcmc.set}. If it is \code{NULL}, the \code{sim.dir} is used to load existing simulation results.}
  \item{sim.dir}{Directory with existing MCMC simulation results. It is only used if \code{mcmc.set} is \code{NULL}.}
  \item{burnin}{Number of burn-in iterations to be discarded from the begining of the chain.}
  \item{chain.id}{Identifier of the MCMC to be used. By default the first chain is used.}
}
\details{
In order to speed-up MCMC convergence, the package contains default values of gamma covariance that were obtained from a previous run (they can be loaded using \code{data(proposal_cov_gammas_cii)}). However, this function allows to obtain new values and overwrite the default values by 
passing the resulting object to the \code{\link{run.tfr.mcmc}} function as the \code{proposal_cov_gammas} argument.
}
\value{
	A list with components:
	  \item{values}{An array of size nr_countries \eqn{\times}{x} 3 \eqn{\times}{x} 3 containing values of the covariance matrix of \eqn{\gamma_{ci}}{gamma_ci} (\eqn{i=1,2,3}) for each country \eqn{c}.}
	  \item{country_codes}{A vector of size nr_countries. A covariance matrix \code{values[i,,]} corresponds to a country with the code \code{country_codes[i]}.}
}

\author{
Leontine Alkema, Hana Sevcikova
}

\seealso{
\code{\link{run.tfr.mcmc}}
}
\examples{
\dontrun{
sim.dir <- file.path(find.package("bayesTFR"), "ex-data", "bayesTFR.output")
cov.gammas <- get.cov.gammas(sim.dir=sim.dir, burnin=0)
m <- run.tfr.mcmc(nr.chains=1, iter=10, proposal_cov_gammas=cov.gammas, verbose=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
