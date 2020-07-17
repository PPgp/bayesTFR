\name{run.tfr.mcmc.extra}
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
    iter = NULL, thin = 1, thin.extra = 1, burnin = 2000,
    parallel = FALSE, nr.nodes = NULL,  my.locations.file = NULL,
    verbose = FALSE, verbose.iter = 100, uncertainty=FALSE, 
    my.tfr.raw.file=ifelse(uncertainty, file.path(find.package("bayesTFR"), "data", "rawTFR.csv"), NULL),
    iso.unbiased=NULL, covariates=c('source', 'method'), cont_covariates=NULL, ...)
}
\arguments{
  \item{sim.dir}{
Directory with an existing simulation.
}
  \item{countries}{Vector of country codes. These include codes of areas and regions (see column \code{country_code} in \code{\link[wpp2012]{UNlocations}}).
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
  \item{thin.extra}{
Thinning interval for the MCMC run for extra countries.
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
\item{my.locations.file}{File name containing user-specified locations. See Details below.}
  \item{verbose}{
Logical switching log messages on and off.
}
\item{verbose.iter}{Integer determining how often (in number of iterations) log messages are outputted during the estimation.}
\item{uncertainty}{Whether past TFR uncertainty is considered. If past TFR uncertainty is considered, countries listed in \code{countries} will be re-simulated with the model account for past TFR estimation. It will take the data from my.tfr.raw.file as the observations, estimate the distribution of these observations with respect to the true TFR. Then instead of treating the observed data as true data, assume true TFR are unknown and include an extra step for estimating past TFR.}
\item{my.tfr.raw.file}{File name of the raw TFR, used when uncertainty is \code{TRUE}. See details in \code{run.tfr.mcmc}.}
\item{iso.unbiased}{Countries whose Vital registration TFR estimates are considered unbiased with 0.0161 standard deviation. See details in \code{run.tfr.mcmc}.}
\item{covariates, cont_covariates}{Categorical and continuous features used in estimating bias and standard deviation. See details in \code{run.tfr.mcmc}.}
\item{\dots}{Additional parameters to be passed to the function \code{\link[snowFT]{performParallel}}, if \code{parallel} is \code{TRUE}.}
}
\details{
The function can be used to make predictions for countries, areas or regions (further denoted as \sQuote{countries}) that were not included in the MCMC estimation (invoked by \code{\link{run.tfr.mcmc}}). It creates MCMC traces for country-specific parameters. The purpose of this function is to have country-specific parameters available in order to be able to generate projections for additional countries or their aggregations, without having to re-run the often time-expensive MCMC simulation.

The set of countries to be considered by this function can be given either by their codes, using the argument \code{countries}, in which case the countries must be included in the UN WPP \code{\link[wpp2012]{tfr}} dataset. Or, it can be given by a user-specific TFR file, using the argument \code{my.tfr.file}. The function considers a union of both arguments. The function will ignore all countries that were used in the existing MCMC simulation for estimating the hyperparameters. Countries that already own country-specific parameters (e.g. because they were included in \code{my.tfr.file} passed to \code{\link{run.tfr.mcmc}}) get their parameters recomputed. Note that all countries should be included in the \code{\link[wpp2012]{UNlocations}} dataset, but unlike in \code{\link{run.tfr.mcmc}}, their  \code{include_code} is ignored. As in the case of \code{\link{run.tfr.mcmc}}, the default dataset of locations \code{\link[wpp2012]{UNlocations}} can be overwritten using a file of the same structure as \code{\link[wpp2012]{UNlocations}} passed via the \code{my.locations.file} argument. This file should be especially used, if TFR is simulated for new locations that are not included in \code{\link[wpp2012]{UNlocations}}.
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
sim.dir <- tempfile()
m <- run.tfr.mcmc(nr.chains=1, iter=20, output.dir=sim.dir, verbose=TRUE)
m <- run.tfr.mcmc.extra(sim.dir=sim.dir, countries=c(908,924), burnin=10, verbose=TRUE)
summary(m, country=924)
pred <- tfr.predict(m, burnin=10, use.tfr3=FALSE, verbose=TRUE)
summary(pred, country=908)
unlink(sim.dir, recursive=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{multivariate}

