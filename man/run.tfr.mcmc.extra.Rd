\name{run.tfr.mcmc.extra}
\Rdversion{1.1}
\alias{run.tfr.mcmc.extra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run MCMC for Extra Countries, Areas or Regions
}
\description{
Run MCMC for extra countries, areas or regions. It uses the posterior distribution of model hyperparameters from an existing simulation to generate country-specific parameters. 
}
\usage{
run.tfr.mcmc.extra(sim.dir = file.path(getwd(), "bayesTFR.output"), 
    countries = NULL, my.tfr.file = NULL, 
    iter = NULL, thin = 1, burnin = 2000,
    parallel = FALSE, nr.nodes = NULL,  verbose = FALSE, ...)
}
\arguments{
  \item{sim.dir}{
Directory with an existing simulation.
}
  \item{countries}{Vector of country codes. These include codes of areas and regions (see column \code{country_code} in \code{\link{UN2010}}).
}
  \item{my.tfr.file}{
File name containing user-specified TFR time series for countries for which the simulation should run (see Details below).
}
  \item{iter}{
Number of iterations to be used for sampling from the posterior distribution of the hyperparameters. By default, the number of iterations used in the existing simulation is taken.
}
  \item{thin}{
Thinning interval for sampling from the posterior distribution of the hyperparameters.
}
  \item{burnin}{
Number of iterations discarded before sampling from the posterior distribution of the hyperparameters. It is also used when computing proposal of gamma covariance matrices (see \code{\link{get.cov.gammas}}).
}
  \item{parallel}{
Logical determining if the simulation should run multiple chains in parallel.
}
  \item{nr.nodes}{
Relevant only if \code{parallel} is \code{TRUE}. It gives the number of nodes for running the simulation in parallel. By default it equals to the number of chains contained in the existing simulation.
}
  \item{verbose}{
Logical switching log messages on and off.
}
\item{\dots}{Additional parameters to be passed to the function \code{\link[snowFT]{performParallel}}, if \code{parallel} is \code{TRUE}.}
}
\details{
The function can be used to make predictions for countries, areas or regions (further denoted as \sQuote{countries}) that were not included in the MCMC estimation (invoked by \code{\link{run.tfr.mcmc}}). It creates MCMC traces for country-specific parameters. The purpose of this function is to have country-specific parameters available in order to be able to generate projections for additional countries or their aggregations, without having to re-run the often time-expensive MCMC simulation.

The set of countries to be considered by this function can be given either by their codes, using the argument \code{countries}, in which case the countries must be included in the UN WPP file (see \code{\link{UN2010}}). Or, it can be given by a user-specific TFR file, using the argument \code{my.tfr.file}. The function considers a union of both arguments. The function will ignore all countries that were used in the existing MCMC simulation for estimating the hyperparameters. Countries that already own country-specific parameters (e.g. because they were included in \code{my.tfr.file} passed to \code{\link{run.tfr.mcmc}}) get their parameters recomputed. Note that all countries must be included in the UN location file (see \link{WPP2010_LOCATIONS}), but unlike in \code{\link{run.tfr.mcmc}}, their  \code{include_code} is ignored. 
}
\note{If there is an existing projection for the directory \code{sim.dir}, use \code{\link{tfr.predict.extra}} to obtain projections for the extra countries used in this function.}
\value{
An object of class \code{\link{bayesTFR.mcmc.set}}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova, Leontine Alkema
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{run.tfr.mcmc}}, \code{\link{tfr.predict.extra}}
}
\examples{
\dontrun{
m <- run.tfr.mcmc(nr.chains=1, iter=20, verbose=TRUE)
m <- run.tfr.mcmc.extra(countries=c(908,924), burnin=10, verbose=TRUE)
summary(m, country=924)
pred <- tfr.predict(burnin=10, verbose=TRUE)
summary(pred, country=908)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{multivariate}

