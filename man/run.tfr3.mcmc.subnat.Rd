\name{run.tfr3.mcmc.subnat}
\alias{run.tfr3.mcmc.subnat}

\title{
Running Phase III MCMCs for Subnational TFR Series
}
\description{
Runs MCMCs to estimate Phase III parameters of TFR for sub-national data, using a Bayesian hierarchical version of an AR(1) model.
}
\usage{
run.tfr3.mcmc.subnat(countries, 
    sim.dir = file.path(getwd(), "bayesTFR.output"), iter = 30000,
    use.world.posterior = FALSE, world.sim.dir = sim.dir, post.burnin = 10000, 
    auto.conf = list(max.loops = 5, iter = 30000, iter.incr = 10000, 
        nr.chains = 3, thin = 40, burnin = 2000), verbose = FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{countries}{Vector of numerical codes of countries the sub-national simulation to run for. They must correspond to country codes from \code{\link[wpp2012]{UNlocations}} and sub-national Phase II MCMC (generated using \code{\link{run.tfr.mcmc.subnat}}) have to exist for those countries.}
  \item{sim.dir}{Directory with a world simulation where the sub-national Phase II MCMCs are stored in its subdirectory \sQuote{subnat}.}
  \item{iter}{Number of iterations to run in each chain. In addition to a single value, it can have the value \sQuote{auto} in which case the function runs for the number of iterations given in the \code{auto.conf} list (see below), then checks if the MCMCs converged (using the \code{auto.conf} settings). If it did not converge, the procedure is repeated until convergence is reached or the number of repetitions exceeded \code{auto.conf$max.loops}.}
  \item{use.world.posterior}{Logical determining of the posterior distribution from a world simulation should be used as prior in this simulation. If \code{TRUE}, the posterior is extracted from \code{world.sim.dir}. If \code{FALSE}, the priors are set to reasonable values.}
  \item{world.sim.dir}{Directory with a world simulation from which the posterior and other parameters are extracted and used for the sub-national simulation. Nothing is written into this directory.}
  \item{post.burnin}{Burnin used for extracting posterior from a world simulation. It is only used if the \code{world.sim.dir} does not contain projections. If it does, the same burnin is used as for generating world projections.}
  \item{auto.conf}{List containing a configuration for an \sQuote{automatic} run (see description of argument \code{iter}). Item \code{iter} gives the number of iterations in the first chunk of the MCMC simulation; item \code{iter.incr} gives the number of iterations in the following chunks; \code{nr.chains} gives the number of chains in all chunks of the MCMC simulation; items \code{thin} and \code{burnin} are used in the convergence diagnostics following each chunk; \code{max.loops} controls the maximum number of chunks. All items must be integer values. This argument is only used if the function argument \code{iter} is set to \sQuote{auto}.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{\dots}{Additional arguments passed to \code{\link{run.tfr3.mcmc}}.}
}
\details{
The function operates on country-specific directories created by \code{\link{run.tfr.mcmc.subnat}}. It calls \code{\link{run.tfr3.mcmc}} for each such directory. As in teh case of \code{\link{run.tfr.mcmc.subnat}}, the prior distributions are set to the posterior of a world simulation, either extracted dynamically (if \code{use.world.posterior=TRUE}) or hard-coded otherwise.
}
\value{
List of \code{\link{bayesTFR.mcmc.set}} objects, one per country.}

%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova and Adrian Raftery
}

\seealso{
\code{\link{run.tfr3.mcmc}}, \code{\link{run.tfr.mcmc.subnat}}, \code{\link{tfr.predict.subnat}}
}

\examples{
\dontrun{
sim.dir <- tempfile()
# run a parent world simulation
m <- run.tfr.mcmc(nr.chains=1, iter=5, output.dir=sim.dir)
# run sub-national Phase II simulation
tfrdata <- file.path(find.package("bayesTFR"), "extdata", "UN2012SubnationalExample.txt")
msub <- run.tfr.mcmc.subnat(countries=124, my.tfr.file=tfrdata,
    nr.chains=1, iter=10, thin=1, sim.dir=sim.dir)
# run sub-national Phase III simulation
msub3 <- run.tfr3.mcmc.subnat(countries=124, iter=10, thin=1, sim.dir=sim.dir, verbose=TRUE)
# continue Phase II MCMC for country 124
m124 <- continue.tfr3.mcmc(iter=10, sim.dir=file.path(sim.dir, 'subnat/c124'))
summary(m124)

# get the m124 object at later time 
m124 <- get.tfr3.mcmc(file.path(sim.dir, 'subnat/c124'))
summary(m124)
unlink(sim.dir, recursive=TRUE)
}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{multivariate}

