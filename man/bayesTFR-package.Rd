\name{bayesTFR-package}
\Rdversion{1.1}
\alias{bayesTFR-package}
\alias{bayesTFR}
\docType{package}
\title{
Bayesian Fertility Projection
}

\description{
Collection of functions for making probabilistic projections of total fertility rate (TFR) for all countries of the world, 
using a Bayesian hierarchical model and the United Nations demographic time series.
}

\details{
\tabular{ll}{
Package: \tab bayesTFR\cr
Version: \tab 1.4-2\cr
Date: \tab 2011-10-13\cr
Depends: \tab R (>= 2.9.0), mvtnorm, coda\cr
Suggests: \tab rworldmap, fields, snowFT, googleVis \cr
License: \tab GPL (>= 2)\cr
URL: \tab \url{http://bayespop.csss.washington.edu}\cr
}

The projection follows a method developed by Alkema et al.  It uses historical data provided by the United Nations (included in the package) to simulate a posterior distribution of total fertility rates for all countries in the world simultaneously. 

The main functions of the package are:
\itemize{
\item \link{run.tfr.mcmc}: Evokes running a Markov Chain Monte Carlo (MCMC) simulation for one or more chains, possibly in parallel.  It results in a posterior sample of the mcmc parameters.  Existing simulation runs can be resumed using \link{continue.tfr.mcmc}.
\item \link{tfr.predict}: Using the posterior parameter samples it derives posterior trajectories of the total fertility rate for all countries.
\item \link{run.tfr.mcmc.extra}: Runs MCMC for extra countries or regions, i.e. for countries not included in the Bayesian hierarchical model. It can be also used for aggregations.
\item \link{tfr.predict.extra}: Generates predictions for extra countries or aggregated regions.
}

The order of the functions above roughly corresponds to a typical workflow when using the package: 1. run an MCMC simulation, 2. generate predictions, 3. analyze results (using the functions below). If there are countries that were not included in steps 1. and 2., or if there are aggregated regions for which a prediction is desired, one proceeds with the two functions at the bottom of the list above, followed by the analyzing functions below.

A number of functions analyzing results are included in the package:
\itemize{
\item \link{tfr.trajectories.plot}: Shows the posterior trajectories for a given country, including their median and given probability intervals.
\item \link{tfr.trajectories.table}: Shows the posterior trajectories for a given country in a tabular form.
\item \link{tfr.map}: Shows a TFR world map for a given projection period.
\item \link{DLcurve.plot}: Shows the posterior curves of the double logistic function used in the simulation, including their median and given probability intervals.
\item \link{tfr.partraces.plot} and \link{tfr.partraces.cs.plot}: Plot the MCMC traces of country-independent parameters and country-specific parameters, respectively.
\item \link{tfr.pardensity.plot} and \link{tfr.pardensity.cs.plot}: Plot the posterior density of the MCMCs for country-independent parameters and country-specific parameters, respectively.
\item \link{summary.bayesTFR.mcmc.set}: Summary function for the MCMC results.
\item \link{summary.bayesTFR.prediction}: Summary function for the prediction results.
}

For MCMC diagnostics, function \link{coda.mcmc.list} creates an object of type \dQuote{mcmc.list} that can be used with the \pkg{coda} package.  Furthermore, function \link{tfr.diagnose} analyzes the MCMCs using the Raftery diagnostics implemented in the \pkg{coda} package and gives information about parameters that did not converge. 
#In case of having multiple chains, function \link{tfr.gelman.diag} computes the Gelman diagnostics for different burnins which can by plotted using \link{tfr.plot.gelman}.  

Existing simulation results can be accessed using the \link{get.tfr.mcmc} function.  An existing prediction can be accessed via \link{get.tfr.prediction}. Existing convergence diagnostics can be accessed using the \link{get.tfr.convergence} and \link{get.tfr.convergence.all} functions.

There are two important datasets used as an input into the MCMC simulations, provided by the United Nations Population Division:
\itemize{
\item \link{UN2008}: Contains the historical data of total fertility rates as available for WPP 2008.
\item \link{WPP2008_LOCATIONS}: Location dataset (for WPP 2008) used as a filter and look-up table.
}
}

\note{
There is a directory \code{ex-data} shipped with the package which contains results from an example simulation, containing three chains with 50 iterations each.  The Example section below shows how these results were created.  These data are used in Example sections throughout the manual.  The user can either reproduce the data in her/his local directory, or use the ones from the package.   
}

\author{
Hana Sevcikova <hanas@uw.edu>, Leontine Alkema <alkema@nus.edu.sg>, Adrian Raftery <raftery@uw.edu>

Maintainer: Hana Sevcikova <hanas@uw.edu>
}
\references{
Hana Sevcikova, Leontine Alkema, Adrian E. Raftery (2011). bayesTFR: An
  R Package for Probabilistic Projections of the Total Fertility Rate.
  Journal of Statistical Software, 43(1), 1-29.
  \url{http://www.jstatsoft.org/v43/i01/}.

L. Alkema, A. E. Raftery, P. Gerland, S. J. Clark, F. Pelletier, Buettner, T.: 
Probabilistic projections of the total fertility rate for all countries. Working Paper nr. 97, Center for Statistics and the Social Sciences, University of Washington. \url{http://www.csss.washington.edu/Papers}
}

\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
\dontrun{
# This command produces output data such as in the directory ex-data
m <- run.tfr.mcmc(nr.chains=3, iter=50, seed=1, verbose=TRUE)
pred <- tfr.predict(m, burnin=30, verbose=TRUE)
summary(pred, country='Ghana')
}
}
