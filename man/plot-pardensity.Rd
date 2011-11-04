\name{tfr.pardensity.plot}
\Rdversion{1.1}
\alias{tfr.pardensity.plot}
\alias{tfr.pardensity.cs.plot}
\title{
Plotting MCMC Parameter Density
}
\description{
Functions for plotting density of the posterior distribution of the MCMC parameters.
}
\usage{
tfr.pardensity.plot(mcmc.list = NULL, 
    sim.dir = file.path(getwd(), "bayesTFR.output"), 
    chain.ids = NULL, par.names = tfr.parameter.names(trans = TRUE), 
    burnin = NULL, dev.ncol=5, low.memory = TRUE, ...)
    
tfr.pardensity.cs.plot(country, mcmc.list=NULL, 
    sim.dir=file.path(getwd(), 'bayesTFR.output'), 
    chain.ids=NULL, par.names=tfr.parameter.names.cs(trans=TRUE), 
    burnin=NULL, dev.ncol=3, low.memory=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{country}{Name or numerical code of a country.}
  \item{mcmc.list}{List of \code{\link{bayesTFR.mcmc}} objects, or an object of class \code{\link{bayesTFR.mcmc.set}} or of class \code{\link{bayesTFR.prediction}}. If it is \code{NULL}, the parameter values are loaded from \code{sim.dir}.}
  \item{sim.dir}{Directory with the MCMC simulation results. It is only used if \code{mcmc.list} is \code{NULL}.}
  \item{chain.ids}{List of MCMC identifiers to be plotted. If it is \code{NULL}, all chains found in \code{mcmc.list} or \code{sim.dir} are plotted.}
  \item{par.names}{Names of parameters for which density should be plotted. By default all (possibly transformed) country-independent parameters are plotted if used within \code{tfr.pardensity.plot}, or country-specific parameters are plotted if used within \code{tfr.pardensity.cs.plot}.}
  \item{burnin}{Number of iterations to be discarded from the beginning of each chain.}
  \item{dev.ncol}{Number of column for the graphics device. If the number of parameters is smaller than \code{dev.ncol}, the number of columns is automatically decreased.}
  \item{low.memory}{Logical indicating if the processing should run in a low-memory mode. If it is \code{FALSE}, traces of all available parameters are loaded into memory. Otherwise, parameters are loaded as they are needed and are not kept in the memory.}
  \item{\dots}{Further arguments passed to the \code{density} function.
}
}
\details{
The functions plot the density of the posterior distribution either for country-independent parameters (\code{tfr.pardensity.plot}) or for country-specific parameters (\code{tfr.pardensity.cs.plot}), one graph per parameter.  One can restrict it to specific chains by setting the \code{chain.ids} argument and to specific parameters by setting the \code{par.names} argument. 

If \code{mcmc.list} is an object of class \code{\link{bayesTFR.prediction}} and if this object contains thinned traces, they are used instead of the full chains. In such a case, \code{burnin} and \code{chain.ids} cannot be modified - their value is set to the one used when the thinned traces were created, namely when running \code{\link{tfr.predict}}. In a situation with long MCMC chains, this approach can  significantly speed-up creation of the density plots.
}


\author{
Hana Sevcikova
}

\seealso{
\code{\link{tfr.partraces.plot}}
}
\examples{
\dontrun{
sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
tfr.pardensity.plot(sim.dir=sim.dir)
tfr.pardensity.cs.plot(country='Ireland', sim.dir=sim.dir, bw=0.2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
