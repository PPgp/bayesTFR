\name{tfr.partraces.plot}
\Rdversion{1.1}
\alias{tfr.partraces.plot}
\alias{tfr.partraces.cs.plot}
\alias{tfr3.partraces.plot}
\alias{tfr3.partraces.cs.plot}
\title{
Plotting MCMC Parameter Traces
}
\description{
Functions for plotting the MCMC parameter traces.
}
\usage{
tfr.partraces.plot(mcmc.list = NULL, 
    sim.dir = file.path(getwd(), "bayesTFR.output"), chain.ids = NULL, 
    par.names = tfr.parameter.names(trans = TRUE), 
    nr.points = NULL, dev.ncol=5, low.memory = TRUE, \dots)
	
tfr.partraces.cs.plot(country, mcmc.list = NULL, 
    sim.dir = file.path(getwd(), "bayesTFR.output"), chain.ids = NULL, 
    par.names = tfr.parameter.names.cs(trans = TRUE), 
    nr.points = NULL, dev.ncol=3, low.memory = TRUE, \dots)
    
tfr3.partraces.plot(mcmc.list = NULL, 
    sim.dir = file.path(getwd(), "bayesTFR.output"), chain.ids = NULL, 
    par.names = tfr3.parameter.names(), 
    nr.points = NULL, dev.ncol=3, low.memory = TRUE, \dots)
	
tfr3.partraces.cs.plot(country, mcmc.list = NULL, 
    sim.dir = file.path(getwd(), "bayesTFR.output"), chain.ids = NULL, 
    par.names = tfr3.parameter.names.cs(), 
    nr.points = NULL, dev.ncol=2, low.memory = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{country}{Name or numerical code of a country. The code can be either numeric or ISO-2 or ISO-3 characters.}
  \item{mcmc.list}{List of \code{bayesTFR.mcmc} objects, or an object of class \code{bayesTFR.mcmc.set} or of class \code{bayesTFR.prediction} (allowed only for Phase II MCMCs). If it is \code{NULL}, the traces are loaded from \code{sim.dir}.}
  \item{sim.dir}{Directory with the MCMC simulation results. It is only used if \code{mcmc.list} is \code{NULL}.}
  \item{chain.ids}{List of MCMC identifiers to be plotted. If it is \code{NULL}, all chains found in \code{mcmc.list} or \code{sim.dir} are plotted.}
  \item{par.names}{Names of parameters for which traces should be plotted. By default all (possibly transformed) country-independent parameters are plotted if used within \code{tfr.partraces.plot} and \code{tfr3.partraces.plot}, or country-specific parameters are plotted if used within \code{tfr.partraces.cs.plot} and \code{tfr3.partraces.cs.plot}.}
  \item{nr.points}{Number of points to be plotted. If \code{NULL}, all points are plotted, otherwise the traces are thinned evenly.}
  \item{dev.ncol}{Number of column for the graphics device. If the number of parameters is smaller than \code{dev.ncol}, the number of columns is automatically decreased.}
  \item{low.memory}{Logical indicating if the processing should run in a low-memory mode. If it is \code{FALSE}, traces of all available parameters are loaded into memory. Otherwise, parameters are loaded as they are needed and are not kept in the memory.}
  \item{\dots}{Additional graphical arguments.}
}
\details{
The functions plot MCMC traces either for country-independent parameters (\code{tfr.partraces.plot} for phase II MCMCs and \code{tfr3.partraces.plot} for phase III MCMCs) or for country-specific parameters (\code{tfr.partraces.cs.plot} for phase II MCMCs and \code{tfr3.partraces.cs.plot} for phase III MCMCs), one graph per parameter.  One can restrict it to specific chains by setting the \code{chain.ids} argument, and to specific parameters by setting the \code{par.names} argument.
}

\author{
Hana Sevcikova
}

\seealso{
\code{\link{coda.list.mcmc}} for retrieving raw values of the traces.
}
\examples{
\dontrun{
sim.dir <- file.path(find.package("bayesTFR"), "ex-data", "bayesTFR.output")
tfr.partraces.plot(sim.dir=sim.dir)
tfr.partraces.cs.plot(country="Netherlands", sim.dir=sim.dir)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
