\name{run.tfr3.mcmc}
\alias{run.tfr3.mcmc}
\alias{continue.tfr3.mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Running Markov Chain Monte Carlo for Parameters of Total Fertility Rate in Phase III
}
\description{
Runs (or continues running) MCMCs for simulating phase III total fertility rate, using a Bayesian hierarchical version of an AR(1) model.
}
\usage{
run.tfr3.mcmc(sim.dir, nr.chains = 3, iter = 50000, thin = 10, 
    replace.output = FALSE, my.tfr.file = NULL, buffer.size = 100, 
    use.extra.countries = FALSE, 
    mu.prior.range = c(0, 2.1), rho.prior.range = c(0, 1 - .Machine$double.xmin), 
    sigma.mu.prior.range = c(1e-05, 0.318), sigma.rho.prior.range = c(1e-05, 0.289), 
    sigma.eps.prior.range = c(1e-05, 0.5), 
    mu.ini = NULL, mu.ini.range = mu.prior.range, 
    rho.ini = NULL, rho.ini.range = rho.prior.range, 
    sigma.mu.ini = NULL, sigma.mu.ini.range = sigma.mu.prior.range, 
    sigma.rho.ini = NULL, sigma.rho.ini.range = sigma.rho.prior.range, 
    sigma.eps.ini = NULL, sigma.eps.ini.range = sigma.eps.prior.range, 
    seed = NULL, parallel = FALSE, nr.nodes = nr.chains, compression.type = "None", 
    auto.conf = list(max.loops = 5, iter = 50000, iter.incr = 20000, nr.chains = 3, 
                    thin = 60, burnin = 10000), 
    verbose = FALSE, verbose.iter = 1000, \dots)
        
continue.tfr3.mcmc(sim.dir, iter, chain.ids=NULL, 
    parallel = FALSE, nr.nodes = NULL, auto.conf = NULL,
    verbose=FALSE, verbose.iter = 1000, \dots)
}

\arguments{
  \item{sim.dir}{Directory with an existing simulation of phase II TFR (see \code{\link{run.tfr.mcmc}}).}
  \item{nr.chains}{Number of MCMC chains to run.}
  \item{iter}{Number of iterations to run in each chain. In addition to a single value, it can have the value \sQuote{auto} in which case the function runs for the number of iterations given in the \code{auto.conf} list (see below), then checks if the MCMCs converged (using the \code{auto.conf} settings). If it did not converge, the procedure is repeated until convergence is reached or the number of repetition exceeded \code{auto.conf$max.loops}.}
  \item{thin}{Thinning interval between consecutive observations to be stored on disk.}
  \item{replace.output}{If \code{TRUE}, previously stored results of a phase III simulation will be overwritten.}
  \item{my.tfr.file}{File name containing user-specified TFR time series for one or more countries. See description of this argument in \code{\link{run.tfr.mcmc}}.}
  \item{buffer.size}{Buffer size (in number of iterations) for keeping data in the memory.}
  \item{use.extra.countries}{By default, only countries are used in the MCMCs that were assigned for estimation (i.e. their \sQuote{include_code} is 2 in the \link{include}) dataset and are in phase III at present time (argument \code{present.year} in \code{\link{run.tfr.mcmc}}). If this argument is \code{TRUE}, countries that were added using \code{\link{run.tfr.mcmc.extra}} and are in phase III are also included.}
  \item{mu.prior.range, rho.prior.range, sigma.mu.prior.range, sigma.rho.prior.range, sigma.eps.prior.range}{Min and max for the prior (uniform) distribution of these paraemters.}
  \item{mu.ini, rho.ini, sigma.mu.ini, sigma.rho.ini, sigma.eps.ini}{Initial value(s) of the parameters. It can be a single value or an array of the size \code{nr.chains}.  By default, if \code{nr.chains} is 1, it is the middle point of the corresponding range. Otherwise, it is uniformly randomly distributed within the range.
}
  \item{mu.ini.range, rho.ini.range, sigma.mu.ini.range, sigma.rho.ini.range, sigma.eps.ini.range}{Min and max for the initial values.}
  \item{seed}{Seed of the random number generator. If \code{NULL} no seed is set.}
  \item{parallel}{Logical determining if the simulation should run multiple chains in parallel. If it is \code{TRUE}, the package \pkg{\link[snowFT]{snowFT}} is required.}
  \item{nr.nodes}{Relevant only if \code{parallel} is \code{TRUE}. It gives the number of nodes for running the simulation in parallel.}
  \item{compression.type}{One of \sQuote{None}, \sQuote{gz}, \sQuote{xz}, \sQuote{bz}, determining type of a compression of the MCMC files.}
  \item{auto.conf}{List containing a configuration for an \sQuote{automatic} run (see description of argument \code{iter}). Item \code{iter} gives the number of iterations in the first chunk of the MCMC simulation; item \code{iter.incr} gives the number of iterations in the following chunks; \code{nr.chains} gives the number of chains in all chunks of the MCMC simulation; items \code{thin} and \code{burnin} are used in the convergence diagnostics following each chunk; \code{max.loops} controls the maximum number of chunks. All items must be integer values. This argument is only used if the function argument \code{iter} is set to \sQuote{auto}.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{verbose.iter}{Integer determining how often (in number of iterations) messages are outputted during the estimation.}
  \item{\dots}{Additional arguments to be passed to the function \code{\link[snowFT]{performParallel}}, if \code{parallel} is \code{TRUE}.}
  \item{chain.ids}{Array of chain identifiers that should be resumed. If it is \code{NULL}, all existing chains are resumed.}
}
\details{
The MCMCs are stored in \code{sim.dir} in a subdirectory called \dQuote{phaseIII}. It has exactly the same structure as phase II MCMCs described in \code{\link{run.tfr.mcmc}}.
}
\value{
	An object of class \code{bayesTFR.mcmc.set} which is a list with two components:
	\item{meta}{An object of class \code{\link{bayesTFR.mcmc.meta}}.}
	\item{mcmc.list}{A list of objects of class \code{\link{bayesTFR.mcmc}}, one for each MCMC.}
}
\references{
Raftery, A.E., Alkema, L. and Gerland, P. (2013). Bayesian Population Projections for the United Nations. 
Statistical Science, in press.
}
\author{
Hana Sevcikova
}

\seealso{
\code{\link{run.tfr.mcmc}}, \code{\link{get.tfr3.mcmc}}
}

\examples{
\dontrun{
sim.dir <- tempfile()
# Runs Phase II MCMCs (must be run before Phase III)
m <- run.tfr.mcmc(nr.chains=1, iter=5, output.dir=sim.dir, verbose=TRUE)
# Runs Phase III MCMCs
m3 <- run.tfr3.mcmc(sim.dir=sim.dir, nr.chains=2, iter=50, thin=1, verbose=TRUE)
m3 <- continue.tfr3.mcmc(sim.dir=sim.dir, iter=10, verbose=TRUE)
summary(m3, burnin=10)
unlink(sim.dir, recursive=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{multivariate}
