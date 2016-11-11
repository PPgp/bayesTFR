\name{DLcurve.plot}
\Rdversion{1.1}
\alias{DLcurve.plot}
\alias{DLcurve.plot.all}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting  Posterior Distribution of the Double Logistic Function
}
\description{
The functions for plotting and retrieving the posterior distribution of the double logistic function used in the simulation of Phase II. Plots include the median and given probability intervals of the distribution.
}
\usage{
DLcurve.plot(mcmc.list, country, burnin = NULL, pi = 80, tfr.max = 10, 
    nr.curves = NULL, predictive.distr = FALSE, ylim = NULL, 
    xlab = 'TFR (reversed)', ylab = 'TFR decrement', main = NULL, 
    show.legend = TRUE, col=c('black', 'red', "#00000020"), \dots)
	
DLcurve.plot.all(mcmc.list = NULL, sim.dir = NULL, 
    output.dir = file.path(getwd(), 'DLcurves'),
    output.type = "png", burnin = NULL, verbose = FALSE, \dots)
    
tfr.world.dlcurves(x, mcmc.list, burnin=NULL, countryUc=NULL, \dots)

tfr.country.dlcurves(x, mcmc.list, country, burnin=NULL, \dots)
}
\arguments{
  \item{mcmc.list}{List of \code{\link{bayesTFR.mcmc}} objects, an object of class \code{\link{bayesTFR.mcmc.set}} or of class \code{\link{bayesTFR.prediction}}. In case of \code{DLcurve.plot.all} if it si \code{NULL}, it is loaded from \code{sim.dir}.}
  \item{country}{Name or numerical code of a country.}
  \item{burnin}{Number of iterations to be discarded from the beginning of parameter traces.}
  \item{pi}{Probability interval. It can be a single number or an array.}
  \item{tfr.max}{Maximum TFR to be shown in the plot.}
  \item{nr.curves}{Number of curves to be plotted. If \code{NULL}, all curves are plotted.}
  \item{predictive.distr}{Logical. If \code{TRUE}, an error term is added to each trajectory.}
  \item{ylim, xlab, ylab, main}{Graphical parameters passed to the \code{plot} function.}
  \item{show.legend}{Logical determining if the legend should be shown.}
  \item{col}{Vector of colors in this order: 1. observed data points, 2. quantiles, 3. trajectories}
  \item{\dots}{For the plotting functions, there are additional graphical parameters. For \code{DLcurve.plot.all}, \code{\dots} contains also arguments
  		\code{pi}, \code{tfr.max} and \code{nr.curves}. For the \code{tfr.*.dlcurves} functions, 
  		these are arguments passed to the underlying functions (\code{predictive.distr} and \code{return.sigma} for obtaining a sample of the standard deviation of the error term ).}
  \item{sim.dir}{Directory with the simulation results. Only relevant, if \code{mcmc.list} is \code{NULL}.}
  \item{output.dir}{Directory into which resulting graphs are stored.}
  \item{output.type}{Type of the resulting files. It can be \dQuote{png}, \dQuote{pdf}, \dQuote{jpeg}, \dQuote{bmp}, \dQuote{tiff}, or \dQuote{postscript}.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{x}{TFR values for which the double logistic should be computed.}
  \item{countryUc}{Country to use the parameter U_c from (start of the fertility transition). If it is not given, the middle point of the prior distribution is used.}
}
\details{
\code{DLcurve.plot} plots double logistic curves for the given country.  \code{DLcurve.plot.all} creates such plots for all countries and stores them in \code{output.dir}. Parameters inputting the double logistic function are either thinned traces created by the \code{\link{tfr.predict}} function (if \code{mcmc.list} is an object of class \code{\link{bayesTFR.prediction}}), or they are selected by equal spacing from the MCMC traces. In the former case, \code{burnin} is set automatically; in the latter case, \code{burnin} defaults to 0 since such object has already been \dQuote{burned}. If \code{nr.curves} is smaller than 2000, the median and probability intervals are computed on a sample of 2000 equally spaced data points, otherwise on all plotted curves.

Function \code{tfr.world.dlcurves} returns the DL curves of the hierarchical distribution, conditioned on the starting point of the fertility transition in a given country (given by the \code{countryUc} argument). Function \code{tfr.country.dlcurves} returns DL curves for the given country. If \code{mcmc.list} is a prediction object, \code{burnin} should not be given, as such object has already been \dQuote{burned}.  
}

\value{
	\code{tfr.world.dlcurves} and \code{tfr.country.dlcurves} return a matrix of size \eqn{N \times M} where \eqn{N} is the number of trajectories and \eqn{M} is the number of values of \eqn{x}. If the argument \code{return.sigma} is set to \code{TRUE}, the return value is a list with the first element being the DL values and the second element being a matrix of the standard deviation of the DL error term \code{sigma_eps}.
}
\author{
Hana Sevcikova, Leontine Alkema
}
\examples{
\dontrun{
sim.dir <- file.path(find.package("bayesTFR"), "ex-data", "bayesTFR.output")
mcmc.set <- get.tfr.mcmc(sim.dir=sim.dir)
DLcurve.plot(country="Burkina Faso", mcmc.set, burnin=15)

# add the median of the hierarchical DL curves
x <- seq(0, 10, length=100)
world <- tfr.world.dlcurves(x, mcmc.set, burnin=15, countryUc="Burkina Faso")
qw <- apply(world, 2, median) 
lines(x, qw, col='blue')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

