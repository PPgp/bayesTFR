\name{bayesTFR-package}
\alias{bayesTFR-package}
\alias{bayesTFR}
\docType{package}
\title{
Bayesian Fertility Projection
}

\description{
Collection of functions for making probabilistic projections of total fertility rate (TFR) for all countries of the world, 
using a Bayesian hierarchical model (BHM) and the United Nations demographic time series. Functions for subnational projections are also available.
}

\details{
\tabular{ll}{
Package: \tab bayesTFR\cr
Version: \tab 6.1-2\cr
Date: \tab 2018-01-18\cr
License: \tab GPL (>= 2)\cr
URL: \tab \url{https://bayespop.csss.washington.edu}\cr
}

The projection follows a method developed by Alkema et al. (2011) and Raftery et al (2013). It uses historical data provided by the United Nations to simulate a posterior distribution of total fertility rates for all countries in the world simultaneously. 

The estimation is split into two parts: 
\enumerate{
\item BHM for fertility in a transition phase (Phase II), as described in Alkema et al. (2011),
\item BHM for fertility in a post-transition phase (Phase III), as described in Raftery et al (2013).
}
The second part is optional and can be replaced by a simple AR(1) process.

The main functions of the package are:
\itemize{
\item \link{run.tfr.mcmc}: Evokes running a Markov Chain Monte Carlo (MCMC) simulation for TFR in Phase II using one or more chains, possibly in parallel.  It results in a posterior sample of the mcmc parameters.  Existing simulation runs can be resumed using \link{continue.tfr.mcmc}.
\item \link{run.tfr3.mcmc}: Starts MCMCs for TFR in Phase III. Existing simulation runs can be resumed using \link{continue.tfr3.mcmc}.
\item \link{tfr.predict}: Using the posterior parameter samples it derives posterior trajectories of the total fertility rate for all countries.
\item \link{run.tfr.mcmc.extra}: Runs MCMC for extra countries or regions, i.e. for countries not included in the Bayesian hierarchical model. It can be also used for aggregations.
\item \link{tfr.predict.extra}: Generates predictions for extra countries or aggregated regions.
}

The order of the functions above roughly corresponds to a typical workflow when using the package: 1. run a Phase II MCMC simulation, 2. run a Phase III MCMC simulation (optional but recommended), 3. generate predictions, 4. analyze results (using the functions below). If there are countries that were not included in steps 1.-3., or if there are aggregated regions for which a prediction is desired, one proceeds with the two functions at the bottom of the list above, followed by the analyzing functions below.

A number of functions analyzing results are included in the package:
\itemize{
\item \link{tfr.trajectories.plot}: Shows the posterior trajectories for a given country, including their median and given probability intervals.
\item \link{tfr.trajectories.table}: Shows the posterior trajectories for a given country in a tabular form.
\item \link{tfr.map}: Shows a TFR world map for a given projection period.
\item \link{DLcurve.plot}: Shows the posterior curves of the double logistic function used in the simulation of PhaseII, including their median and given probability intervals.
\item \link{tfr.partraces.plot} and \link{tfr.partraces.cs.plot}: Plot the Phase II MCMC traces of country-independent parameters and country-specific parameters, respectively. \link{tfr3.partraces.plot} and \link{tfr3.partraces.cs.plot} do the same for Phase III MCMCs.
\item \link{tfr.pardensity.plot} and \link{tfr.pardensity.cs.plot}: Plot the posterior density of the Phase II MCMCs for country-independent parameters and country-specific parameters, respectively. \link{tfr3.pardensity.plot} and \link{tfr3.pardensity.cs.plot} do the same for Phase III MCMCs.
\item \link{summary.bayesTFR.mcmc.set}: Summary function for the MCMC results.
\item \link{summary.bayesTFR.prediction}: Summary function for the prediction results.
}

For MCMC diagnostics, functions \link{coda.list.mcmc} and \link{coda.list.mcmc3} create an object of type \dQuote{mcmc.list} that can be used with the \pkg{coda} package.  Furthermore, function \link{tfr.diagnose} and \link{tfr3.diagnose} analyze the MCMCs using the Raftery diagnostics implemented in the \pkg{coda} package and gives information about parameters that did not converge. 
%In case of having multiple chains, function \link{tfr.gelman.diag} computes the Gelman diagnostics for different burnins which can by plotted using \link{tfr.plot.gelman}.  

Existing simulation results can be accessed using the \link{get.tfr.mcmc} (Phase II) and \link{get.tfr3.mcmc} (Phase III) functions.  An existing prediction can be accessed via \link{get.tfr.prediction}. Existing convergence diagnostics can be accessed using the \link{get.tfr.convergence}, \link{get.tfr.convergence.all}, \link{get.tfr3.convergence} and \link{get.tfr3.convergence.all} functions.

The historical national TFR data are taken from one of the packages \pkg{wpp2017} (default), \pkg{wpp2015}, \pkg{wpp2012} or \pkg{wpp2010}, depending on users settings.

Subnational TFR projections can be generated using \link{tfr.predict.subnat}. In this case, historical data must be provided by the user. Existing projections can be accessed from disk via \link{get.regtfr.prediction}. 
}

\note{
There is a directory \code{ex-data} shipped with the package which contains results from an example simulation, containing one chain with 60 iterations.  The Example section below shows how these results were created.  These data are used in Example sections throughout the manual.  The user can either reproduce the data in her/his local directory, or use the ones from the package.   
}

\author{
Hana Sevcikova <hanas@uw.edu>, Leontine Alkema <alkema@nus.edu.sg>, Adrian Raftery <raftery@uw.edu>, Bailey Fosdick <bfosdick@uw.edu>, Patrick Gerland (gerland@un.org)

Maintainer: Hana Sevcikova <hanas@uw.edu>
}
\references{
Hana Sevcikova, Leontine Alkema, Adrian E. Raftery (2011). bayesTFR: An
  R Package for Probabilistic Projections of the Total Fertility Rate.
  Journal of Statistical Software, 43(1), 1-29.
  \url{http://www.jstatsoft.org/v43/i01/}.

L. Alkema, A. E. Raftery, P. Gerland, S. J. Clark, F. Pelletier, Buettner, T., Heilig, G.K. (2011). \href{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3367999}{Probabilistic Projections of the Total Fertility Rate for All Countries.} Demography, Vol. 48, 815-839.

Raftery, A.E., Alkema, L. and Gerland, P. (2014). \href{http://www.ncbi.nlm.nih.gov/pubmed/25324591}{Bayesian Population Projections for the United Nations.}
Statistical Science, Vol. 29, 58-68. 
}

\keyword{ package }

\examples{
\dontrun{
# This command produces output data such as in the directory ex-data
sim.dir <- tempfile()
# Phase II MCMCs
m <- run.tfr.mcmc(nr.chains=1, iter=60, output.dir=sim.dir, seed=1, verbose=TRUE)
# Phase III MCMCs (not included in the package)
m3 <- run.tfr3.mcmc(sim.dir=sim.dir, nr.chains=2, iter=100, thin=1, seed=1, verbose=TRUE)
# Prediction
pred <- tfr.predict(m, burnin=30, burnin3=50, verbose=TRUE)
summary(pred, country='Ghana')
unlink(sim.dir, recursive=TRUE)
}
}
% to reproduce the example included in the package exactly:
% m <- run.tfr.mcmc(nr.chains=1, iter=60, seed=1, verbose=TRUE, compression.type='bz', start.year=1950)
% pred <- tfr.predict(m, burnin=30, verbose=TRUE, save.as.ascii=0, use.tfr3=FALSE)
% # Directory with thinned chains and summary files in "predictions" removed to save disk space.  