\name{get.tfr.parameter.traces}
\Rdversion{1.1}
\alias{get.tfr.parameter.traces}
\alias{get.tfr.parameter.traces.cs}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Accessing MCMC Parameter Traces
}
\description{
Functions for accessing traces of the MCMC parameters, either country-independent or country-specific.
}
\usage{
get.tfr.parameter.traces(mcmc.list, par.names = tfr.parameter.names(), 
    burnin = 0, thinning.index = NULL, thin = NULL)
    
get.tfr.parameter.traces.cs(mcmc.list, country.obj, 
    par.names=tfr.parameter.names.cs(), 
    burnin=0, thinning.index=NULL, thin=NULL)
}
\arguments{
  \item{mcmc.list}{List of \code{\link{bayesTFR.mcmc}} objects.}
  \item{country.obj}{Country object list (see \code{\link{get.country.object}}).}
  \item{par.names}{Names of country-independent parameters (in case of \code{get.tfr.parameter.traces}) or country-specific parameters (in case of \code{get.tfr.parameter.traces.cs}) to be included.}
  \item{burnin}{Burnin indicating how many iterations should be removed from the beginning of each chain.}
  \item{thinning.index}{Index of the traces for thinning. If it is \code{NULL}, \code{thin} is used. \code{thinning.index} does not include \code{burnin}. For example, if there are two MCMC chains of length 1000, \code{burnin=200} and we want a sample of length 400, then the value should be \code{thinning.index=seq(1,1600, length=400)}.}
  \item{thin}{Alternative to \code{thinning.index}. In the above example it would be \code{thin=4}.}
}

\value{
Both functions return a matrix with columns being the parameters and rows being the MCMC values, attached to one another in case of multiple chains. \code{get.tfr.parameter.traces} contains country-independent parameters, \code{get.tfr.parameter.traces.cs} contains country-specific parameters.
}


\author{
Hana Sevcikova
}


\seealso{
\code{\link{coda.mcmc.list}} for another way of retrieving parameter traces. 
}

\examples{
sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
m <- get.tfr.mcmc(sim.dir)
tfr.values <- get.tfr.parameter.traces(m$mcmc.list, burnin=10, par.names="sigma0")
\dontrun{
hist(tfr.values, main=colnames(tfr.values))
}
tfr.values.cs <- get.tfr.parameter.traces.cs(m$mcmc.list, 
                    get.country.object('Canada', meta=m$meta),
                    burnin=10, par.names="Triangle_c4")
\dontrun{
hist(tfr.values.cs, main=colnames(tfr.values.cs))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

