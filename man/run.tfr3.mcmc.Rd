\name{run.tfr3.mcmc}
\alias{run.tfr3.mcmc}
\alias{continue.tfr3.mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Running Markov Chain Monte Carlo for Parameters of Total Fertility Rate in Phase III
}
\description{
Runs (or continues running) MCMCs for simulating Phase III total fertility rate, using a Bayesian hierarchical version of an AR(1) or AR(2) model.
}
\usage{
run.tfr3.mcmc(sim.dir, nr.chains = 3, iter = 50000, thin = 10, 
    replace.output = FALSE, my.tfr.file = NULL, buffer.size = 500, 
    use.extra.countries = FALSE, ar.order = 1, momentum = FALSE, momentum.lag = 1,
    mu.prior.range = c(0, 2.1), mu.prior.mean = NULL, mu.prior.sd = NULL, mu.fixed = NULL,
    rho.prior.range = c(0, 1 - .Machine$double.xmin), 
    sigma.mu.prior.range = c(1e-05, 0.318), sigma.rho.prior.range = c(1e-05, 0.289), 
    sigma.eps.prior.range = c(1e-05, 0.5), 
    alpha1.prior.range = c(0, 1), alpha2.prior.range = c(-1, 1),
    sigma.alpha1.prior.range = c(1e-05, 0.289), 
    sigma.alpha2.prior.range = c(1e-05, 2*0.289),
    beta.prior.range = c(0, 1), sigma.beta.prior.range = c(1e-05, 0.289),
    mu.ini = NULL, mu.ini.range = mu.prior.range, 
    rho.ini = NULL, rho.ini.range = rho.prior.range, 
    sigma.mu.ini = NULL, sigma.mu.ini.range = sigma.mu.prior.range, 
    sigma.rho.ini = NULL, sigma.rho.ini.range = sigma.rho.prior.range, 
    sigma.eps.ini = NULL, sigma.eps.ini.range = sigma.eps.prior.range, 
    alpha1.ini = NULL, alpha1.ini.range = alpha1.prior.range,
    alpha2.ini = NULL, alpha2.ini.range = alpha2.prior.range,
    sigma.alpha1.ini = NULL, sigma.alpha1.ini.range = sigma.alpha1.prior.range,
    sigma.alpha2.ini = NULL, sigma.alpha2.ini.range = sigma.alpha2.prior.range,
    beta.ini = NULL, beta.ini.range = beta.prior.range,
    sigma.beta.ini = NULL, sigma.beta.ini.range = sigma.beta.prior.range,
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
  \item{ar.order}{Order of the autoregressive process, either 1 (default) or 2. See Details.}
  \item{momentum}{Logical. If \code{TRUE}, the AR(2) process is parametrized by the persistence \eqn{\rho_c}{rho.c} and the momentum \eqn{\beta_c}{beta.c} instead of the partial autocorrelations, and \code{ar.order} is set to \code{momentum.lag} + 1. See Details.}
  \item{momentum.lag}{Number of years (or time periods) \eqn{L} over which the velocity in the momentum model is computed. Only used if \code{momentum} is \code{TRUE}. See Details.}
  \item{mu.prior.range, rho.prior.range, sigma.mu.prior.range, sigma.rho.prior.range, sigma.eps.prior.range}{Min and max for the prior (uniform) distribution of these paraemters. \code{rho.prior.range} and \code{sigma.rho.prior.range} are only used if \code{ar.order} is 1 or \code{momentum} is \code{TRUE}.}
  \item{mu.prior.mean, mu.prior.sd}{If given, the prior distribution of the world mean \eqn{\bar{\mu}}{mu} is normal with this mean and standard deviation, truncated to \code{mu.prior.range}. By default, the prior is uniform on \code{mu.prior.range}. Both arguments must be given together.}
  \item{mu.fixed}{If given, the world mean \eqn{\bar{\mu}}{mu} is fixed to this value and not estimated; the country-specific \eqn{\mu_c}{mu.c} and the other parameters are estimated conditionally on it. This can be used for example to keep \eqn{\bar{\mu}}{mu} at the value estimated by an AR(1) model while estimating the AR(2) or momentum parameters, as \eqn{\bar{\mu}}{mu} is only weakly identified if the persistence is close to 1. It cannot be combined with \code{mu.prior.mean}.}
  \item{alpha1.prior.range, alpha2.prior.range, sigma.alpha1.prior.range, sigma.alpha2.prior.range}{Min and max for the prior (uniform) distribution of the world means and standard deviations of the partial autocorrelations. Only used if \code{ar.order} is 2 and \code{momentum} is \code{FALSE}.}
  \item{beta.prior.range, sigma.beta.prior.range}{Min and max for the prior (uniform) distribution of the world mean and standard deviation of the momentum parameter. Only used if \code{momentum} is \code{TRUE}.}
  \item{mu.ini, rho.ini, sigma.mu.ini, sigma.rho.ini, sigma.eps.ini, alpha1.ini, alpha2.ini, sigma.alpha1.ini, sigma.alpha2.ini, beta.ini, sigma.beta.ini}{Initial value(s) of the parameters. It can be a single value or an array of the size \code{nr.chains}.  By default, if \code{nr.chains} is 1, it is the middle point of the corresponding range. Otherwise, it is uniformly randomly distributed within the range.
}
  \item{mu.ini.range, rho.ini.range, sigma.mu.ini.range, sigma.rho.ini.range, sigma.eps.ini.range, alpha1.ini.range, alpha2.ini.range, sigma.alpha1.ini.range, sigma.alpha2.ini.range, beta.ini.range, sigma.beta.ini.range}{Min and max for the initial values.}
  \item{seed}{Seed of the random number generator. If \code{NULL} no seed is set.}
  \item{parallel}{Logical determining if the simulation should run multiple chains in parallel. If it is \code{TRUE}, the package \pkg{\link[snowFT]{snowFT}} is required.}
  \item{nr.nodes}{Relevant only if \code{parallel} is \code{TRUE}. It gives the number of nodes for running the simulation in parallel.}
  \item{compression.type}{One of \sQuote{None}, \sQuote{gz}, \sQuote{xz}, \sQuote{bz}, determining type of a compression of the MCMC files. Important: Do not use this option for a long MCMC simulation as this tends to cause very long run times due to slow reading!}
  \item{auto.conf}{List containing a configuration for an \sQuote{automatic} run (see description of argument \code{iter}). Item \code{iter} gives the number of iterations in the first chunk of the MCMC simulation; item \code{iter.incr} gives the number of iterations in the following chunks; \code{nr.chains} gives the number of chains in all chunks of the MCMC simulation; items \code{thin} and \code{burnin} are used in the convergence diagnostics following each chunk; \code{max.loops} controls the maximum number of chunks. All items must be integer values. This argument is only used if the function argument \code{iter} is set to \sQuote{auto}.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{verbose.iter}{Integer determining how often (in number of iterations) messages are outputted during the estimation.}
  \item{\dots}{Additional parameters to be passed to the function \code{\link[snowFT]{performParallel}}, if \code{parallel} is \code{TRUE}.}
  \item{chain.ids}{Array of chain identifiers that should be resumed. If it is \code{NULL}, all existing chains are resumed.}
}
\details{
By default (\code{ar.order = 1}), the model is given by equations (12)-(14) in Raftery et al. (2014), i.e. \eqn{f_{c,t+1} - \mu_c = \rho_c (f_{c,t} - \mu_c) + \varepsilon_{c,t}}.

If \code{ar.order = 2}, the model is extended to
\deqn{f_{c,t+1} - \mu_c = \rho_{1,c} (f_{c,t} - \mu_c) + \rho_{2,c} (f_{c,t-1} - \mu_c) + \varepsilon_{c,t}.}
To ensure stationarity, the AR coefficients are parametrized by the partial autocorrelations \eqn{\alpha_{1,c}}{alpha1.c} and \eqn{\alpha_{2,c}}{alpha2.c} as \eqn{\rho_{1,c} = \alpha_{1,c} (1 - \alpha_{2,c})}, \eqn{\rho_{2,c} = \alpha_{2,c}}, with
\eqn{\alpha_{1,c} \sim TN_{[0,1)}(\bar{\alpha}_1, \sigma^2_{1,\alpha})} and \eqn{\alpha_{2,c} \sim TN_{[-1,1)}(\bar{\alpha}_2, \sigma^2_{2,\alpha})}. The world-level parameters are called \sQuote{alpha1}, \sQuote{alpha2}, \sQuote{sigma.alpha1} and \sQuote{sigma.alpha2}, the country-specific parameters \sQuote{alpha1.c} and \sQuote{alpha2.c} (see \code{\link{tfr3.parameter.names}}). The process is conditioned on the first two Phase III observations, and thus only countries with at least three Phase III observations are included in the estimation. Other Phase III countries use the world-level parameters in projections.

If \code{momentum = TRUE}, the process is written as
\deqn{f_{c,t+1} - \mu_c = \rho_c (f_{c,t} - \mu_c) + \beta_c (f_{c,t} - f_{c,t-L})/L + \varepsilon_{c,t},}
where \eqn{L} is given by \code{momentum.lag}, i.e. the velocity is the average change of TFR over the last \eqn{L} periods. This is an AR(\eqn{L+1}) process with coefficients \eqn{\rho_c + \beta_c/L} at lag 1 and \eqn{-\beta_c/L} at lag \eqn{L+1}; for \eqn{L=1}, \eqn{\rho_{1,c} = \rho_c + \beta_c} and \eqn{\rho_{2,c} = -\beta_c}. The process is conditioned on the first \eqn{L+1} Phase III observations, and thus only countries with at least \eqn{L+2} Phase III observations are included in the estimation. The priors are
\eqn{\rho_c \sim TN_{[0,1)}(\bar{\rho}, \sigma^2_{\rho})} and \eqn{\beta_c \sim TN_{[0,1)}(\bar{\beta}, \sigma^2_{\beta})}. The persistence \eqn{\rho_c} equals the sum of the AR coefficients and thus plays the same role as \eqn{\rho_c} in AR(1), including the amount of information the data carry about \eqn{\mu_c}. The momentum \eqn{\beta_c} carries a part of the recent velocity of TFR into the next period; in projections, this results in a drift that decays over time. The model reduces to AR(1) for \eqn{\beta_c = 0}. The world-level parameters are called \sQuote{rho}, \sQuote{beta}, \sQuote{sigma.rho} and \sQuote{sigma.beta}, the country-specific parameters \sQuote{rho.c} and \sQuote{beta.c}.

The MCMCs are stored in \code{sim.dir} in a subdirectory called \dQuote{phaseIII}. It has exactly the same structure as phase II MCMCs described in \code{\link{run.tfr.mcmc}}.
}
\value{
	An object of class \code{bayesTFR.mcmc.set} which is a list with two components:
	\item{meta}{An object of class \code{\link{bayesTFR.mcmc.meta}}.}
	\item{mcmc.list}{A list of objects of class \code{\link{bayesTFR.mcmc}}, one for each MCMC.}
}
\references{
Raftery, A.E., Alkema, L. and Gerland, P. (2014). Bayesian Population Projections for the United Nations.
Statistical Science, Vol. 29, 58-68. \doi{10.1214/13-STS419}. 
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
# AR(2) model for Phase III
m3 <- run.tfr3.mcmc(sim.dir=sim.dir, nr.chains=2, iter=50, thin=1, ar.order=2, replace.output=TRUE)
summary(m3, burnin=10)
unlink(sim.dir, recursive=TRUE)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{multivariate}
