\name{coda.list.mcmc}
\Rdversion{1.1}
\alias{coda.list.mcmc}
\alias{coda.list.mcmc3}
\alias{coda.mcmc}
\alias{coda.mcmc.bayesTFR.mcmc}

\title{
Convertion to coda's Objects
}
\description{
The functions convert MCMC traces (simulated using \code{\link{run.tfr.mcmc}} and \code{\link{run.tfr3.mcmc}}) into objects that can be used with the \pkg{coda} package.
}
\usage{
coda.list.mcmc(mcmc = NULL, country = NULL, chain.ids = NULL, 
    sim.dir = file.path(getwd(), "bayesTFR.output"), 
    par.names = tfr.parameter.names(), 
    par.names.cs = tfr.parameter.names.cs(), 
    rm.const.pars = FALSE, burnin = 0, 
    low.memory = FALSE, \dots)
    
coda.list.mcmc3(mcmc = NULL, country = NULL, chain.ids = NULL, 
    sim.dir = file.path(getwd(), "bayesTFR.output"), 
    par.names = tfr3.parameter.names(), 
    par.names.cs = tfr3.parameter.names.cs(), 
    burnin = 0, low.memory = FALSE, \dots)
	
\method{coda.mcmc}{bayesTFR.mcmc}(mcmc, country = NULL, 
    par.names = NULL, par.names.cs = NULL, 
    burnin = 0, thin = 1, \dots)
}
\arguments{
\item{mcmc}{In \code{coda.mcmc}, it is an object of class \code{\link{bayesTFR.mcmc}}. In \code{coda.list.mcmc} and \code{coda.list.mcmc3}, it is either a list of \code{\link{bayesTFR.mcmc}} objects, or an object of class \code{\link{bayesTFR.mcmc.set}} or in case of \code{coda.list.mcmc} it can be \code{\link{bayesTFR.prediction}}. If it is \code{NULL}, the MCMCs are loaded from \code{sim.dir}. Either \code{mcmc} or \code{sim.dir} must be given.}
\item{country}{Country name or code. It is used in connection with the \code{par.names.cs} argument (see below). 
}
  \item{chain.ids}{
Vector of chain identifiers. By default, all chains available in the \code{mcmc.list} object are included.
}
  \item{sim.dir}{
Directory with the MCMC simulation results. Only used if \code{mcmc.list} is \code{NULL}.
}
  \item{par.names}{
Names of country-independent parameters to be included. In \code{coda.mcmc} the default names are \code{\link{tfr.parameter.names}()} if the \code{mcmc} object is an MCMC of phase II or  \code{\link{tfr3.parameter.names}()} if the MCMC is of phase III.
}
  \item{par.names.cs}{
Names of country-specific parameters to be included. The argument \code{country} is used to filter out traces that correspond to a specific country. If \code{country} is not given, for each parameter, traces for all countries are included. In \code{coda.mcmc} the default names are \code{\link{tfr.parameter.names.cs}()} if the \code{mcmc} object is an MCMC of phase II or  \code{\link{tfr3.parameter.names.cs}()} if the MCMC is of phase III.
}
  \item{rm.const.pars}{
Logical indicating if parameters with constant values should be removed.
}
  \item{burnin}{
Burnin indicating how many iterations should be removed from the beginning of each chain.
}
  \item{low.memory}{
Logical indicating if the function should run in a memory-efficient mode.
}
\item{thin}{Thinning interval.}
  \item{\dots}{
Additional arguments passed to the \pkg{coda}'s \code{\link[coda]{mcmc}} function.
}
}
\details{
Function \code{coda.list.mcmc} is for accessing all chains of phase II MCMCs; Function \code{coda.list.mcmc3} is for accessing all chains of phase III MCMCs.
}
\value{
The function \code{coda.list.mcmc} and \code{coda.list.mcmc3} return an object of class \dQuote{mcmc.list}. The function \code{coda.mcmc} returns an object of class \dQuote{mcmc}, both defined in the \pkg{coda} package.
}

\author{
Hana Sevcikova
}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
\keyword{ ts }
