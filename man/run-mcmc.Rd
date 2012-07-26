\name{run.tfr.mcmc}
\Rdversion{1.1}
\alias{run.tfr.mcmc}
\alias{continue.tfr.mcmc}
\alias{bayesTFR.mcmc.set}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Running Markov Chain Monte Carlo for Parameters of Total Fertility Rate 
}
\description{
Runs (or continues running) MCMCs for simulating the total fertility rate of all countries of the world, using a Bayesian hierarchical model.
}
\usage{
run.tfr.mcmc(nr.chains = 3, iter = 62000, 
    output.dir = file.path(getwd(), "bayesTFR.output"), 
    thin = 1, replace.output = FALSE, 
    start.year = 1950, present.year = 2010, wpp.year = 2010, 
    my.tfr.file = NULL, buffer.size = 100, 
    U.c.low = 5.5, U.up = 8.8, U.width = 3,
    mean.eps.tau0 = -0.25, sd.eps.tau0 = 0.4, nu.tau0 = 2, 
    Triangle_c4.low = 1, Triangle_c4.up = 2.5, 
    Triangle_c4.trans.width = 2,
    Triangle4.0 = 0.3, delta4.0 = 0.8, nu4 = 2,
    S.low = 3.5, S.up = 6.5, S.width = 0.5, 
    a.low = 0, a.up = 0.2, a.width = 0.02, 
    b.low = a.low, b.up = a.up, b.width = 0.02, 
    sigma0.low = 0.01, sigma0.up = 0.6, sigma0.width = 0.1, 
    sigma0.min = 0.001, 
    const.low = 0.8, const.up = 2, const.width = 0.3, 
    d.low = 0.05, d.up = 0.5, d.trans.width = 1, 
    chi0 = -1.5, psi0 = 0.6, nu.psi0 = 2, 
    alpha0.p = c(-1, 0.5, 1.5), delta0 = 1, nu.delta0 = 2, 
    dl.p1 = 9, dl.p2 = 9, 
    S.ini = NULL, a.ini = NULL, b.ini = NULL, sigma0.ini = NULL, 
    Triangle_c4.ini = NULL, const.ini = NULL, gamma.ini = 1, 
    proposal_cov_gammas = NULL,
    seed = NULL, parallel = FALSE, nr.nodes = nr.chains, 
    save.all.parameters = FALSE, 
    auto.conf = list(max.loops=5, iter=62000, iter.incr=10000, 
        nr.chains=3, thin=80, burnin=2000),
    verbose = FALSE, verbose.iter = 10, ...)
		
continue.tfr.mcmc(iter, chain.ids=NULL, 
    output.dir=file.path(getwd(), "bayesTFR.output"), 
    parallel = FALSE, nr.nodes = NULL, auto.conf = NULL,
    verbose=FALSE, verbose.iter = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nr.chains}{Number of MCMC chains to run.}
  \item{iter}{Number of iterations to run in each chain. In addition to a single value, it can have the value \sQuote{auto} in which case the function runs for the number of iterations given in the \code{auto.conf} list (see below), then checks if the MCMCs converged (using the \code{auto.conf} settings). If it did not converge, the procedure is repeated until convergence is reached or the number of repetition exceeded \code{auto.conf$max.loops}.}
  \item{output.dir}{Directory which the simulation output should be written into.}
  \item{thin}{Thinning interval between consecutive observations to be stored on disk.}
  \item{replace.output}{If \code{TRUE}, existing outputs in \code{output.dir} will be replaced by results of this simulation.}
  \item{start.year}{Start year for using historical data.}
  \item{present.year}{End year for using historical data.}
  \item{wpp.year}{Year for which WPP data is available (see also \link{UN2010}).}
  \item{my.tfr.file}{File name containing user-specified TFR time series for one or more countries. It can include a subset of countries contained in the UN WPP file, as well as a set of new countries. In the former case,
  	the function replaces the corresponding country data from the UN WPP file by values in this file. Only columns are replaced that match column names of the WPP file, and in addition, columns \sQuote{last.observed} and \sQuote{include_code} are used, if present. Countries are matched using the column \sQuote{country_code} (see also \link{UN2010}). In addition, in order the countries to be included in the simulation, in both cases (whether they are included in the UN TFR file or not), they must be contained in the UN table of locations and the \sQuote{include_code} in that table must be set to 2 (see \link{WPP2010_LOCATIONS}). If the column \sQuote{include_code} is present in \code{my.tfr.file}, its value overwrites the corresponding include code in the location file, unless is -1.  See also Details below.}
  \item{buffer.size}{Buffer size (in number of iterations) for keeping data in the memory. The smaller the \code{buffer.size} the more often will the process access the hard disk and thus, the slower the run. On the other hand, the smaller the \code{buffer.size} the less data will be lost in case of failure.}
  \item{U.c.low, U.up}{Lower and upper bound of the parameter \eqn{U_c}, the start level of the fertility transition. The lower bound is set for each country as the maximum of \code{U.c.low} and the minimum of historical TFR for that country.}
  \item{U.width}{Width for slice sampling used when updating the \eqn{U_c} parameter.}
  \item{mean.eps.tau0, sd.eps.tau0}{Mean and standard deviation of the prior distribution of \eqn{m_{\tau}}{m_tau} which is the mean of the distortion terms \eqn{\epsilon_{c,\tau}}{eps_{c,tau}} in start periods \eqn{\tau_c}{tau_c}.}
  \item{nu.tau0}{The shape parameter of the prior Gamma distribution of \eqn{1/s_{\tau}^2}{1/s_tau^2} is \code{nu.tau0}/2. \eqn{s_{\tau}}{s_tau} is standard deviation of the distortion terms in start periods \eqn{\tau_c}{tau_c}.}
  \item{Triangle_c4.low, Triangle_c4.up}{Lower and upper bound of the \eqn{\Delta_{c4}}{Triangle_c4} parameter.}
  \item{Triangle_c4.trans.width}{Width for slice sampling used when updating the logit-transformed \eqn{\Delta_{c4}}{Triangle_c4} parameter.}
   \item{Triangle4.0, delta4.0}{Mean and standard deviation of the prior distribution of the \eqn{\Delta_4}{Triangle_4} parameter which is the hierarchical mean of the logit-transformed \eqn{\Delta_{c4}}{Triangle_c4}.}
   \item{nu4}{The shape parameter of the prior Gamma distribution of \eqn{1/\delta_4^2}{1/delta_4^2} is \code{nu4}/2. \eqn{\delta_4}{delta_4} is standard deviation of the logit-transformed \eqn{\Delta_{c4}}{Triangle_c4}.}
  \item{S.low, S.up}{Lower and upper bound of the uniform prior distribution of the \eqn{S} parameter which is the TFR at which the distortion has maximum variance.}
  \item{S.width}{Width for slice sampling used when updating the \eqn{S} parameter.}
  \item{a.low, a.up}{Lower and upper bound of the uniform prior distribution of the \eqn{a} parameter which is a coefficient for linear decrease of the TFR for TFR larger than \eqn{S}.}
  \item{a.width}{Width for slice sampling used when updating the \eqn{a} parameter.}  
  \item{b.low, b.up}{Lower and upper bound of the uniform prior distribution of the \eqn{b} parameter which is a  coefficient for linear decrease of the TFR for TFR smaller than \eqn{S}.}
  \item{b.width}{Width for slice sampling used when updating  the \eqn{b} parameter.}
  \item{sigma0.low, sigma0.up}{Lower and upper bound of the uniform prior distribution of the \eqn{\sigma_0}{sigma_0} parameter. \eqn{\sigma_0^2}{sigma_0^2} is the maximum variance of the distortions at TFR equals \eqn{S}.}
  \item{sigma0.width}{Width for slice sampling used when updating the \eqn{\sigma_0}{sigma_0} parameter.}
  \item{sigma0.min}{Minimum value that \eqn{\sigma_0}{sigma_0} can take.}
  \item{const.low, const.up}{Lower and upper bound of the uniform prior distribution of the \eqn{c} parameter which is the multiplier of standard deviation of the distortions before 1975.}
  \item{const.width}{Width for slice sampling used when updating the \eqn{c} parameter.}
  \item{d.low, d.up}{Lower and upper bound of the parameter \eqn{d_c}, the maximum annual decrement for country \eqn{c}. (Note that in Alkema et al. this parameter is a five-years decrement.)}
  \item{d.trans.width}{Width for slice sampling used when updating the logit-transformed \eqn{d_c} parameter.}
  \item{chi0, psi0}{Prior mean and standard deviation of the \eqn{\chi}{chi} parameter which is the hierarchical mean of logit-transformed maximum decline parameter \eqn{d_c}.} 
  \item{nu.psi0}{The shape parameter of the prior Gamma distribution of \eqn{1/\psi^2}{1/psi^2} is \code{nu.psi0}/2. \eqn{\psi}{psi} is the standard devation of logit-transformed maximum decline parameter \eqn{d_c}.}
  \item{alpha0.p}{Vector of prior means of the \eqn{\alpha_i}{alpha_i} parameters, \eqn{i=1,2,3}. \eqn{\alpha_i}{alpha_i} is the hierarchical mean of \eqn{\gamma_{ci}}{gamma_{ci}}.}
  \item{delta0}{Prior standard deviation of the \eqn{\alpha_i}{alpha_i} parameters. It is a single value, i.e. the same standard deviation is used for all \eqn{i}.}
  \item{nu.delta0}{The shape parameter of the prior Gamma distribution of \eqn{1/\delta_i^2}{1/delta_i^2} is \code{nu.delta0}/2. \eqn{\delta_i}{1/delta_i} is the standard deviation of \eqn{\gamma_{ci}}{gamma_{ci}}.}
  \item{dl.p1, dl.p2}{Values of the parameters \eqn{p_1} and \eqn{p_2} of the double logistic function.}
  \item{S.ini}{Initial value for the \eqn{S} parameter. It can be a single value or an array of the size \code{nr.chains}.  By default, if \code{nr.chains} is 1, it is the middle point of the interval [\code{S.low, S.up}]. Otherwise, it is equally spaced distributed between \code{S.low} and \code{S.up}.}
  \item{a.ini}{Initial value for the \eqn{a} parameter. It can be a single value or an array of the size \code{nr.chains}.  By default, if \code{nr.chains} is 1, it is the middle point of the interval [\code{a.low, a.up}]. Otherwise, it is equally spaced distributed between \code{a.low} and \code{a.up}.}
  \item{b.ini}{Initial value for the \eqn{b} parameter. It can be a single value or an array of the size \code{nr.chains}.  By default, if \code{nr.chains} is 1, it is the middle point of the interval [\code{b.low, b.up}]. Otherwise, it is equally spaced distributed between \code{b.low} and \code{b.up}.}
  \item{sigma0.ini}{Initial value for the \eqn{\sigma_0}{sigma_0} parameter. It can be a single value or an array of the size \code{nr.chains}.  By default, if \code{nr.chains} is 1, it is the middle point of the interval [\code{sigma0.low, sigma0.up}]. Otherwise, it is equally spaced distributed between \code{sigma0.low} and \code{sigma0.up}.}
  \item{Triangle_c4.ini}{Initial value for the \eqn{\Delta_{c4}}{Triangle_c4} parameter. It can be a single value or an array of the size \code{nr.chains}.  By default, if \code{nr.chains} is 1, it is the middle point of the interval [\code{Triangle_c4.low, Triangle_c4.up}]. Otherwise, it is equally spaced distributed between \code{Triangle_c4.low} and \code{Triangle_c4.up}.}
  \item{const.ini}{Initial value for the \eqn{c} parameter. It can be a single value or an array of the size \code{nr.chains}.  By default, if \code{nr.chains} is 1, it is the middle point of the interval [\code{const.low, const.up}]. Otherwise, it is equally spaced distributed between \code{const.low} and \code{const.up}.}
  \item{gamma.ini}{Initial value for the \eqn{\gamma_c}{gamma_c} parameter. The same value is used for all countries.}
  \item{proposal_cov_gammas}{Proposal for the gamma covariance matrices for each country. It should be a list with two values: \code{values} and \code{country_codes}. The structure corresponds to the object returned by the function \code{\link{get.cov.gammas}}. By default the covariance matrices are obtained using \code{data(proposal_cov_gammas_cii)}. This argument overwrite the defaults for countries contained the argument.}
  \item{seed}{Seed of the random number generator. If \code{NULL} no seed is set. It can be used to generate reproducible results.}
  \item{parallel}{Logical determining if the simulation should run multiple chains in parallel. If it is \code{TRUE}, the package \pkg{\link[snowFT]{snowFT}} is required.}
  \item{nr.nodes}{Relevant only if \code{parallel} is \code{TRUE}. It gives the number of nodes for running the simulation in parallel. By default it equals to the number of chains.}
  \item{save.all.parameters}{If \code{TRUE}, the distortion terms \eqn{\epsilon_{c,t}}{eps_ct} for all \eqn{t} are stored on disk, otherwise not.}
  \item{auto.conf}{List containing a configuration for an \sQuote{automatic} run (see description of argument \code{iter}). Item \code{iter} gives the number of iterations in the first chunk of the MCMC simulation; item \code{iter.incr} gives the number of iterations in the following chunks; \code{nr.chains} gives the number of chains in all chunks of the MCMC simulation; items \code{thin} and \code{burnin} are used in the convergence diagnostics following each chunk; \code{max.loops} controls the maximum number of chunks. All items must be integer values. This argument is only used if the function argument \code{iter} is set to \sQuote{auto}.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{verbose.iter}{Integer determining how often (in number of iterations) log messages are outputted during the estimation.}
  \item{\dots}{Additional parameters to be passed to the function \code{\link[snowFT]{performParallel}}, if \code{parallel} is \code{TRUE}.}
  \item{chain.ids}{Array of chain identifiers that should be resumed. If it is \code{NULL}, all existing chains in \code{output.dir} are resumed.}
  }
\details{
	The function \code{run.tfr.mcmc} creates an object of class \code{\link{bayesTFR.mcmc.meta}} and stores it in \code{output.dir}.  It launches \code{nr.chains} MCMCs, either sequentially or in parallel.  Parameter traces of each chain are stored (as compressed ASCII files) in a subdirectory of \code{output.dir}, called \code{mc}\emph{x} where \emph{x} is the identifier of that chain.  There is one file per parameter, named after the parameter with the suffix \dQuote{.txt.gz}.  Country-specific parameters (\eqn{U, d, \gamma}{U, d, gamma}) have the suffix \code{_c}\emph{y} where \emph{y} is the country code.  In addition to the trace files, each \code{mc}\emph{x} directory contains the object \code{\link{bayesTFR.mcmc}} in binary format.  All chain-specific files  are written into disk after the first, last and each \code{buffer.size}-th iteration.
	
Using the function \code{continue.tfr.mcmc} one can continue simulating an existing MCMCs by \code{iter} iterations for either all or selected chains.

For simulation of the hyperparameters of the Bayesian hierarchical model, all countries are used that are included in the UN WPP file (\link{UN2010}), possibly complemented by the \code{my.tfr.file}, that have \code{include_code} (in the UN location file) equal to 2. The hyperparameters are used to simulate country-specific parameters, which is done for all countries with \code{include_code} equal 1 or 2. Note that users can overwrite the UN include code by including the column \code{include_code} into the the \code{my.tfr.file}. The following values are recognized: -1 (do not overwrite the UN WPP include code), 0 (ignore), 1 (include in prediction but not estimation), 2 (include in both, estimation and prediction). Thus, the set of countries included in the estimation and prediction can be fully user-specific.

}
\value{
	An object of class \code{bayesTFR.mcmc.set} which is a list with two components:
	\item{meta}{An object of class \code{\link{bayesTFR.mcmc.meta}}.}
	\item{mcmc.list}{A list of objects of class \code{\link{bayesTFR.mcmc}}, one for each MCMC.}
}
\references{
L. Alkema, A. E. Raftery, P. Gerland, S. J. Clark, F. Pelletier, Buettner, T., Heilig, G.K. (2011). Probabilistic Projections of the Total Fertility Rate for All Countries. Demography, Vol. 48, 815-839. (A Working Paper version at \url{http://www.csss.washington.edu/Papers} [nr. 97])
}
\author{
	Hana Sevcikova, Leontine Alkema
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{get.tfr.mcmc}}, \code{\link{summary.bayesTFR.mcmc.set}}.
}
\examples{
\dontrun{
m <- run.tfr.mcmc(nr.chains=1, iter=5, verbose=TRUE)
summary(m)
m <- continue.tfr.mcmc(iter=5, verbose=TRUE)
summary(m)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{multivariate}
%