\name{tfr.predict.subnat}
\alias{tfr.predict.subnat}

\title{
Generating Posterior Trajectories of Subnational TFR
}
\description{
Generates posterior trajectories of the total fertility rate for subregions of given countries, using the Scale-AR(1) method.
}
\usage{
tfr.predict.subnat(countries, my.tfr.file, 
                sim.dir = file.path(getwd(), "bayesTFR.output"), 
                end.year = 2100, start.year = NULL, subdir = "predictions", 
                output.dir = NULL, annual = NULL, nr.traj = NULL, seed = NULL, 
                min.tfr = 0.5, ar.pars = NULL, save.as.ascii = 0, verbose = TRUE)
}

\arguments{
  \item{countries}{Vector of numerical country codes or country names.}
  \item{my.tfr.file}{Tab-separated ASCII file containing the subnational TFR data. See Details for more information on its format.}
  \item{sim.dir}{Simulation directory with the national projections generated using \code{\link{tfr.predict}}.}
 \item{end.year}{End year of the projections.}
 \item{start.year}{Start year of the projections. By default, projections start at the same time point as the national projections.} 
   \item{subdir}{Subdirectory of \code{sim.dir} containing the national predictions.}
   \item{output.dir}{Directory into which the resulting prediction objects and the trajectories are stored. See below for details.}
\item{annual}{Logical indicating if the subnational projection should be on an annual scale or a 5-year scale. By default, 
the scale is matched to the national simulation. If the subnational and national scales are not the same, the national trajectories are either interpolated (if \code{annual = TRUE} and the national simulation is not annual) or averaged to 5-year values (if \code{annual = FALSE} and the national simulation is annual).}
\item{nr.traj}{Number of trajectories to be generated. If \code{NULL}, the number of trajectories in the national projections is used.}
\item{seed}{Seed of the random number generator. If NULL no seed is set. It can be used to generate reproducible projections.}
\item{min.tfr}{Lower bound on TFR.}
\item{ar.pars}{Named vector containing the parameter estimates of the AR(1) process. If given, it must have elements called \code{mu}, \code{rho} and \code{sigma}. By default for a 5-year simulation, \code{c(mu = 1, rho = 0.92464, sigma = 0.04522)} is used. For an annual simulation these default parameters are scaled to \code{c(mu = 1, rho = 0.98445, sigma = 0.02086)}, see details below.}
\item{save.as.ascii}{Either a number determining how many trajectories should be converted into an ASCII file, or \dQuote{all} in which case all trajectories are converted. By default no conversion is performed.}
\item{verbose}{Logical switching log messages on and off.}
}

\details{
The function implements the methodology described in Sevcikova et al (2017). Given a set of national bayesTFR projections, it applies the Scale-AR(1) model to each national trajectory and each subregion of given countries which yield subnational TFR projections. 

The file on subnational data passed in \code{my.tfr.file} has to have a column \dQuote{country_code} with numerical values corresponding to countries given in the argument \code{countries}, and column \dQuote{reg_code} giving the numerical identifier of each subregion. Column \dQuote{name} should be used for subregion name, and column \dQuote{country} for country name. An optional column  \dQuote{include_code} can be used to eliminate entries from processing. Entries with values of 1 or 2 will be included, all others will be ignored. Column \dQuote{last.observed} can be used to define which time period contains the last observed data point (given as integer, e.g. year in the middle of the time period). Remaining columns define the time periods, e.g. \dQuote{2000-2005}, \dQuote{2005-2010} for a 5-year simulation, or \dQuote{2020}, \dQuote{2021} for an annual simulation. The package contains an example of such dataset, see Example below.

The default AR(1) parameters were designed for a 5-year simulation, see Sevcikova et al (2017) for more detail. These are \eqn{\mu = 1, \rho = 0.92464, \sigma = 0.04522}. We use the following conversion for the autoregressive parameter \eqn{\rho} and the standard deviation \eqn{\sigma} if an annual AR(1) process is desired: \eqn{\rho^* = \rho^{(1/5)}, \sigma^* = \sigma * \sqrt{((1 - \rho^{(2/5)})/(1 - \rho^2))}}. The long-term mean \eqn{\mu} stays the same for both processes. Thus, the annual parameters are  \code{c(mu = 1, rho = 0.98445, sigma = 0.02086)}. Note that if the \code{ar.pars} argument is specified by the user, it is assumed that the parameters have been scaled appropriately and thus, no conversion takes place.

Argument \code{output.dir} gives a location on disk where results of the function should be stored. If it is \code{NULL} (default),  results are stored in the same directory as the national projections. In both cases a subdirectory called \dQuote{subnat} is created in which each country has its own subfolder with the country code in its name. Each such subfolder contains the same type of outputs as in the national case generated using \code{\link{tfr.predict}}, most importantly a directory \dQuote{predictions} with trajectories for each region.
}
\value{
A list of objects of class \code{\link{bayesTFR.prediction}}. The name of each element includes its country code. Not all elements of the class \code{\link{bayesTFR.prediction}} are available. For example, no \code{mcmc.set} is attached to these objects. Thus, not all functions that work with \code{\link{bayesTFR.prediction}} can be applied to these results.
}

\note{
Even though the resulting object contains subnational results, the names  of its elements are the same as in the national case. This allows to apply the same functions on both objects (subnational and national). However, it means that sometimes the meaning of the elements or function arguments does not match the subnational context. For example, various functions expect the argument \code{country}. When a subnational object is passed to such function, \code{country} means a subregion. 
}

\references{
Hana Sevcikova, Adrian E. Raftery, Patrick Gerland (2018). Probabilistic Projection of Subnational Total Fertility Rates. Demographic Research, Vol. 38(60), 1843-1884. \doi{10.4054/DemRes.2018.38.60}.
}

\author{
Hana Sevcikova
}


\seealso{
\code{\link{get.regtfr.prediction}}, \code{\link{tfr.predict}}
}

\examples{
# View the example data
my.subtfr.file <- file.path(find.package("bayesTFR"), 'extdata', 'subnational_tfr_template.txt')
subtfr <- read.delim(my.subtfr.file, check.names=FALSE)
head(subtfr)

# Directory with national projections (contains 30 trajectories for each country)
nat.dir <- file.path(find.package("bayesTFR"), "ex-data", "bayesTFR.output")

# Subnational projections for Australia and Canada
subnat.dir <- tempfile()
preds <- tfr.predict.subnat(c(36, 124), my.tfr.file=my.subtfr.file,
    sim.dir=nat.dir, output.dir=subnat.dir, start.year=2013)
names(preds)
get.countries.table(preds[["36"]])
summary(preds[["36"]], "Queensland")
tfr.trajectories.plot(preds[["36"]], "Queensland")

# plot subnational and national TFR in one plot
nat.pred <- get.tfr.prediction(nat.dir)
tfr.trajectories.plot(preds[["36"]], 186, pi=80, half.child.variant=FALSE)
tfr.trajectories.plot(nat.pred, "Australia", half.child.variant=FALSE,
      add=TRUE, col=rep("darkgreen", 5), nr.traj=0, show.legend=FALSE)
legend("topright", c("regional TFR", "national TFR"), col=c("red", "darkgreen"), 
  lty=1, bty='n')

# Retrieve trajectories
trajs.Alberta <- get.tfr.trajectories(preds[["124"]], "Alberta")
summary(t(trajs.Alberta))

# cleanup
unlink(subnat.dir)

# See more examples in ?get.regtfr.prediction
}

\keyword{models}
\keyword{ts}

