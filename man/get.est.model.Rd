\name{get.est.model}
\alias{get.bias.model}
\alias{get.std.model}
\alias{tfr.bias.sd}

\title{
Accessing estimated bias and standard deviations
}
\description{
Functions for obtaining bias and standard deviation of the estimated models as well as the model fits.
}
\usage{
tfr.bias.sd(mcmc.list = NULL, country = NULL, sim.dir = NULL, \dots)
get.bias.model(mcmc.list = NULL, country = NULL, sim.dir = NULL, \dots)
get.std.model(mcmc.list = NULL, country = NULL, sim.dir = NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmc.list}{Object of class \code{\link{bayesTFR.mcmc.set}} corresponding to Phase II MCMCs. If it is \code{NULL}, the object is loaded from the directory given by \code{sim.dir}.}
  \item{country}{Name or numerical code of a country. It can also be given as ISO-2 or ISO-3 characters.}
  \item{sim.dir}{Directory with the MCMC simulation results. Only used if \code{mcmc.list} is not given.}
  \item{\dots}{Not used.}
}
\details{
Functions \code{get.bias.model} and \code{get.std.model} are used to obtain the model fit for estimated bias and standard deviation, respectively, when uncertainty about input data is taken into account. These are used in the MCMC steps stored in \code{mcmc.list}. Function \code{tfr.bias.sd} combines both infos into one object.
}
\value{
  Functions \code{get.bias.model} and \code{get.std.model} return a list with
  \item{model}{\code{lm} object corresponding to the linear model used to estimate the bias (in case of \code{get.bias.model}) and standard deviation (in case of \code{get.std.model}).}
  \item{table}{\code{data.frame} object storing the bias/standard deviation of all possible combinations in the raw data sets for the given country.}
  Function \code{tfr.bias.sd} consolidates these items into a single list where the elements are \code{model_bias}, \code{model_sd} and \code{table}.
}
\author{
Peiran Liu, Hana Sevcikova
}
\examples{
\dontrun{
sim.dir <- tempfile()
mcmc.set <- run.tfr.mcmc(nr.chains = 1, iter = 10, 
    output.dir = sim.dir, uncertainty = TRUE)
tfr.bias.sd(mcmc.set, "Nigeria")
unlink(sim.dir, recursive = TRUE)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ get model }
