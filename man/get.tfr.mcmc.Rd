\name{get.tfr.mcmc}
\Rdversion{1.1}
\alias{get.tfr.mcmc}
\alias{has.tfr.mcmc}
\alias{tfr.mcmc}
\alias{tfr.mcmc.list}
\alias{get.tfr3.mcmc}
\alias{has.tfr3.mcmc}
\title{
Accessing MCMC Results
}
\description{
The function \code{get.tfr.mcmc} retrieves results of an MCMC simulation of Phase II and creates an object of class \code{\link{bayesTFR.mcmc.set}}. Function \code{has.tfr.mcmc} checks the existence of such results. Functions \code{get.tfr3.mcmc} and \code{has.tfr3.mcmc} do the same for Phase III MCMCs. Function \code{tfr.mcmc} extracts a single chain and \code{tfr.mcmc.list} extracts several or all chains from the simulation results. 
}
\usage{
get.tfr.mcmc(sim.dir = file.path(getwd(), "bayesTFR.output"), 
    chain.ids = NULL, low.memory = TRUE, burnin = 0, verbose = FALSE)

has.tfr.mcmc(sim.dir)

get.tfr3.mcmc(sim.dir = file.path(getwd(), "bayesTFR.output"), \dots)

has.tfr3.mcmc(sim.dir)

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
  \item{\dots}{Arguments passed to \code{get.tfr.mcmc}.}
}
\details{
Function \code{get.tfr.mcmc} is an accessor of results generated using \code{\link{run.tfr.mcmc}} and \code{\link{continue.tfr.mcmc}}. Function \code{get.tfr3.mcmc} can be used to access results generated using \code{\link{run.tfr3.mcmc}} and \code{\link{continue.tfr3.mcmc}}. 
}
\value{
\code{get.tfr.mcmc} and \code{get.tfr3.mcmc} return an object of class \code{\link{bayesTFR.mcmc.set}}. \code{has.tfr.mcmc} and \code{has.tfr3.mcmc} return a logical value. \code{tfr.mcmc} returns an object of class \code{\link{bayesTFR.mcmc}}, and \code{tfr.mcmc.list} returns a list of \code{\link{bayesTFR.mcmc}} objects.
}

\author{
Hana Sevcikova
}

\seealso{
\code{\link{bayesTFR.mcmc.set}}
}
\examples{
sim.dir <- file.path(find.package("bayesTFR"), "ex-data", "bayesTFR.output")
m <- get.tfr.mcmc(sim.dir)
summary(m)

# summary of the single chains
for(mc in tfr.mcmc.list(m)) print(summary(mc))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

