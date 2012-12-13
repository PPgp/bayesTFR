\name{get.total.iterations}
\Rdversion{1.1}
\alias{get.total.iterations}
\alias{get.stored.mcmc.length}
\title{
Total Number of Iterations
}
\description{
Function \code{get.total.iterations} gives the total number of iterations of MCMCs summed over chains whith burnin being subtracted from each chain. Function \code{get.stored.mcmc.length} gives the total length of the MCMCs stored on disk minus those iterations that correspond to burnin. Result of the latter will be different from the former only if the MCMCs were run with value of \code{thin} larger than one. 
}
\usage{
get.total.iterations(mcmc.list, burnin = 0)

get.stored.mcmc.length(mcmc.list, burnin = 0)
}
\arguments{
  \item{mcmc.list}{List of \code{\link{bayesTFR.mcmc}} objects.}
  \item{burnin}{Number of iterations to be subtracted from each chain.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A single number.
}

\author{
Hana Sevcikova
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{bayesTFR.mcmc}}
}
\examples{
sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", "bayesTFR.output")
mcmc.set <- get.tfr.mcmc(sim.dir=sim.dir)
get.total.iterations(mcmc.set$mcmc.list) # 60=2x30
get.total.iterations(mcmc.set$mcmc.list, burnin=10) # 40=2x(30-10)

\dontrun{
sim.dir <- tempfile()
m <- run.tfr.mcmc(iter=10, nr.chains=2, output.dir=sim.dir, thin=5, verbose=TRUE)
get.total.iterations(m$mcmc.list) # 20=2x10
get.stored.mcmc.length(m$mcmc.list) # 6=2x3
unlink(sim.dir, recursive=TRUE)
}
}

\keyword{ manip }

