\name{run.tfr.mcmc.subnat}
\alias{run.tfr.mcmc.subnat}

\title{
Running Phase II MCMCs for Subnational TFR Series
}
\description{
Runs MCMCs to estimate Phase II parameters of TFR for sub-national data, using a Bayesian hierarchical model.  
}
\usage{
run.tfr.mcmc.subnat(countries, my.tfr.file, 
    sim.dir = file.path(getwd(), "bayesTFR.output"), 
    nr.chains = 3, iter = 30000, thin = 10, 
    start.year = 1860, present.year = 2010,     
    use.world.posterior = FALSE, world.sim.dir = NULL, 
    post.burnin = 2000, buffer.size = 100, seed = NULL,
    parallel = FALSE, nr.nodes = nr.chains, compression.type = 'None',
    auto.conf = list(max.loops = 5, iter = 30000, iter.incr = 10000, 
        nr.chains = 3, thin = 40, burnin = 2000), 
    verbose = FALSE, verbose.iter = 100, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{countries}{Vector of numerical codes of countries the sub-national simulation to run for. They must correspond to country codes from \code{\link[wpp2012]{UNlocations}}.}
	\item{my.tfr.file}{File name containing TFR data for the sub-regions of \code{countries}. It is a tab delimited ASCII file and must have columns \sQuote{country_code} that corresponds to codes given in \code{countries} (i.e. it has the same value for all rows of the same country); \sQuote{reg_code} giving the numerical code of the sub-regions (must be unique within country); \sQuote{name} giving the name of the sub-regions; multiple columns with the TFR data each named as \eqn{xxxx-yyyy} where \eqn{xxxx} and \eqn{yyyy} are the start and end year of the corresponding period, e.g. \sQuote{1990-1995}. Optionally, the file can have columns \sQuote{last_observed} containing the year of the last observation, and \sQuote{include_code} determining if the record shold be used for estimation (see \code{\link{include}}). The package contains an example of such file, called \dQuote{UN2012SubnationalExample.txt} located in the \dQuote{extdata} directory.}
  \item{sim.dir}{Directory with a world simulation. Results of the subnational simulation are stored into a subdirectory \sQuote{subnat} of this directory. If \code{world.sim.dir} is not given, various meta parameters, e.g. \code{wpp.year}, are taken from this simulation.}
  \item{nr.chains}{Number of MCMC chains to run.}
  \item{iter}{Number of iterations to run in each chain. In addition to a single value, it can have the value \sQuote{auto} in which case the function runs for the number of iterations given in the \code{auto.conf} list (see below), then checks if the MCMCs converged (using the \code{auto.conf} settings). If it did not converge, the procedure is repeated until convergence is reached or the number of repetitions exceeded \code{auto.conf$max.loops}.}
  \item{thin}{Thinning interval between consecutive observations to be stored on disk.}
  \item{start.year}{Start year for using historical data.}
  \item{present.year}{End year for using historical data.}
  \item{use.world.posterior}{Logical determining of the posterior distribution from a world simulation should be used as prior in this simulation. If \code{TRUE}, the posterior is extracted from \code{world.sim.dir}. If \code{FALSE}, the priors are set to reasonable values.}
  \item{world.sim.dir}{Directory with a world simulation from which the posterior and other parameters are extracted and used for the sub-national simulation. Nothing is written into this directory.}
  \item{post.burnin}{Burnin used for extracting posterior from a world simulation. It is only used if the \code{world.sim.dir} does not contain projections. If it does, the same burnin is used as for generating world projections.}
  \item{seed}{Seed of the random number generator. If \code{NULL} no seed is set. It can be used to generate reproducible results.}
  \item{parallel}{Logical determining if the simulation should run multiple chains in parallel. If it is \code{TRUE}, the package \pkg{\link[snowFT]{snowFT}} is required.}
  \item{nr.nodes}{Relevant only if \code{parallel} is \code{TRUE}. It gives the number of nodes for running the simulation in parallel. By default it equals to the number of chains.}
  \item{compression.type}{One of \sQuote{None}, \sQuote{gz}, \sQuote{xz}, \sQuote{bz}, determining type of a compression of the MCMC files.}
  \item{auto.conf}{List containing a configuration for an \sQuote{automatic} run (see description of argument \code{iter}). Item \code{iter} gives the number of iterations in the first chunk of the MCMC simulation; item \code{iter.incr} gives the number of iterations in the following chunks; \code{nr.chains} gives the number of chains in all chunks of the MCMC simulation; items \code{thin} and \code{burnin} are used in the convergence diagnostics following each chunk; \code{max.loops} controls the maximum number of chunks. All items must be integer values. This argument is only used if the function argument \code{iter} is set to \sQuote{auto}.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{verbose.iter}{Integer determining how often (in number of iterations) log messages are outputted during the estimation.}
  \item{\dots}{Additional parameters to be passed to the function \code{\link[snowFT]{performParallel}}, if \code{parallel} is \code{TRUE}.}
}
\details{
The sub-national simulation runs one MCMC simulation per country, analogous to a world simulation with sub-regions being countries and country being the world. The prior distribution is set to be a posterior from a world simulation. If \code{use.world.posterior} is \code{FALSE}, the values of the priors are hard-coded to values from the latest UN wpp 2012 simulation. Results of each such country simulation are stored in a directory called \dQuote{c\eqn{xxx}} where \eqn{xxx} is the country's numerical code. The structure of such directory is the same as in the case of a world simulation as described in \code{\link{run.tfr.mcmc}}. Thus, an object of class \code{\link{bayesTFR.mcmc.set}} can be extracted from each country directory and used with various functions of the package.
}
\value{
	List of \code{\link{bayesTFR.mcmc.set}} objects, one per country.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova and Adrian Raftery
}
%\note{
%%  ~~further notes~~
%}


\seealso{
\code{\link{run.tfr.mcmc}}, \code{\link{run.tfr3.mcmc.subnat}}, \code{\link{tfr.predict.subnat}}
}
\examples{
\dontrun{
sim.dir <- tempfile()
# run a parent world simulation
m <- run.tfr.mcmc(nr.chains=1, iter=5, output.dir=sim.dir)
# run sub-national simulation
tfrdata <- file.path(find.package("bayesTFR"), "extdata", "UN2012SubnationalExample.txt")
msub <- run.tfr.mcmc.subnat(countries=c(36, 124), my.tfr.file=tfrdata,
    nr.chains=1, iter=10, thin=1, sim.dir=sim.dir, verbose=TRUE)
# continue MCMC for country 36
m36 <- continue.tfr.mcmc(iter=10, output.dir=file.path(sim.dir, 'subnat/c36'))
summary(m36)
# summary for 124
# either
summary(msub[['124']])
# or
m124 <- get.tfr.mcmc(file.path(sim.dir, 'subnat/c124'))
summary(m124)
unlink(sim.dir, recursive=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
\keyword{ multivariate }
