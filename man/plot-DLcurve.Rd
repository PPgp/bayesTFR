\name{DLcurve.plot}
\Rdversion{1.1}
\alias{DLcurve.plot}
\alias{DLcurve.plot.all}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting  Posterior Distribution of the Double Logistic Function
}
\description{
The functions plot the posterior distribution of the double logistic function used in the simulation of Phase II, including their median and given probability intervals.
}
\usage{
DLcurve.plot(mcmc.list, country, burnin = NULL, pi = 80, tfr.max = 10, 
    nr.curves = NULL, predictive.distr = FALSE, ylim = NULL, 
    xlab = 'TFR (reversed)', ylab = 'TFR decrement', main = NULL, 
    show.legend = TRUE, \dots)
	
DLcurve.plot.all(mcmc.list = NULL, sim.dir = NULL, 
    output.dir = file.path(getwd(), 'DLcurves'),
    output.type = "png", burnin = NULL, verbose = FALSE, \dots)
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
  \item{\dots}{Additional graphical parameters. For \code{DLcurve.plot.all}, \code{\dots} contains also arguments
  		\code{pi}, \code{tfr.max} and \code{nr.curves}.}
  \item{sim.dir}{Directory with the simulation results. Only relevant, if \code{mcmc.list} is \code{NULL}.}
  \item{output.dir}{Directory into which resulting graphs are stored.}
  \item{output.type}{Type of the resulting files. It can be \dQuote{png}, \dQuote{pdf}, \dQuote{jpeg}, \dQuote{bmp}, \dQuote{tiff}, or \dQuote{postscript}.}
  \item{verbose}{Logical switching log messages on and off.}
}
\details{
\code{DLcurve.plot} plots double logistic curves for the given country.  \code{DLcurve.plot.all} creates such plots for all countries and stores them in \code{output.dir}. Parameters inputting the double logistic function are either thinned traces created by the \code{\link{tfr.predict}} function (if \code{mcmc.list} is an object of class \code{\link{bayesTFR.prediction}}), or they are selected by equal spacing from the MCMC traces. In the former case, \code{burnin} is set automatically; in the latter case, \code{burnin} defaults to 0. If \code{nr.curves} is smaller than 2000, the median and probability intervals are computed on a sample of 2000 equally spaced data points, otherwise on all plotted curves.
}

\author{
Hana Sevcikova, Leontine Alkema
}
\examples{
\dontrun{
sim.dir <- file.path(find.package("bayesTFR"), "ex-data", "bayesTFR.output")
mcmc.set <- get.tfr.mcmc(sim.dir=sim.dir)
DLcurve.plot(country="Burkina Faso", mcmc.set, burnin=15)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }

