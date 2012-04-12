\name{get.tfr.mcmc}
\Rdversion{1.1}
\alias{get.tfr.mcmc}
\alias{has.tfr.mcmc}
\alias{tfr.mcmc}
\alias{tfr.mcmc.list}
\title{
Accessing MCMC Results
}
\description{
The function \code{get.tfr.mcmc} retrieves results of an MCMC simulation and creates an object of class \code{\link{bayesTFR.mcmc.set}}. Function \code{has.tfr.mcmc} checks the existence of such results. Function \code{tfr.mcmc} extracts a single chain and \code{tfr.mcmc.list} extracts several or all chains from the simulation results. 
}
\usage{
get.tfr.mcmc(sim.dir = file.path(getwd(), "bayesTFR.output"), 
    chain.ids = NULL, low.memory = TRUE, burnin = 0, verbose = FALSE)

has.tfr.mcmc(sim.dir)

tfr.mcmc(mcmc.set, chain.id)

tfr.mcmc.list(mcmc.set, chain.ids=NULL)
}
\arguments{
  \item{sim.dir}{Directory where the simulation results are stored.}
  \item{chain.ids}{Chain identifiers in case only specific chains should be included in the resulting object. By default, all available chains are included.}
  \item{low.memory}{If \code{FALSE} full MCMC traces are loaded into memory.}
  \item{burnin}{Burnin used for loading traces. Only relevant, if \code{low.memory=FALSE}.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{chain.id}{Chain identifier.}
  \item{mcmc.set}{Object of class \code{\link{bayesTFR.mcmc.set}}.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{get.tfr.mcmc} returns an object of class \code{\link{bayesTFR.mcmc.set}}. \code{has.tfr.mcmc} returns a logical value. \code{tfr.mcmc} returns an object of class \code{\link{bayesTFR.mcmc}}, and \code{tfr.mcmc.list} returns a list of \code{\link{bayesTFR.mcmc}} objects.
}

\author{
Hana Sevcikova
}

\seealso{
\code{\link{bayesTFR.mcmc.set}}
}
\examples{
sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", "bayesTFR.output")
m <- get.tfr.mcmc(sim.dir)
summary(m)

# summary of the single chains
for(mc in tfr.mcmc.list(m)) print(summary(mc))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
