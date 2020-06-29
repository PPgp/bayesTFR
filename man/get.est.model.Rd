\name{get.est.model}
\Rdversion{1.1}
\alias{get.bias.model}
\alias{get.std.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Accessing estimated bias and standard deviations
}
\description{
Functions for obtaining bias and standard deviation estimated models as well as their fits.
}
\usage{
get.bias.model(mcmc.list=NULL, country.code=NULL, ISO.code=NULL, sim.dir=NULL)
get.std.model(mcmc.list=NULL, country.code=NULL, ISO.code=NULL, sim.dir=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmc.list}{Object of class \code{\link{bayesTFR.mcmc.set}} corresponding Phase II MCMCs. If it is \code{NULL}, the object is loaded from the directory given by \code{sim.dir}.}
  \item{country.code}{Name or numerical code of a country. If it is \code{NULL}, \code{ISO.code} need to be provided with ISO-3 character.}
  \item{ISO.code}{Country code as ISO-3 character. Only used if \code{country.code} is set to be \code{NULL}.}
  \item{sim.dir}{Directory with the MCMC simulation results.}
}
\details{
Functions \code{get.bias.model} and \code{get.std.model} are used to obtain bias and standard deviation estimated models as well as their fits. They are used in the MCMC steps stored in \code{mcmc.list}.
}
\value{
  \item{model}{\code{lm} object describe the linear model estimated for bias or standard deviation model.}
  \item{table}{\code{data.frame} object storing the bias or standard deviation of all possible combinations in the raw data sets for a certain country.}
}
\author{
Peiran Liu, Hana Sevcikova
}
\examples{
\dontrun{
mcmc.set <- run.tfr.mcmc(nr.chains=1, iter=10, output.dir='test.dir', replace.output=TRUE, uncertainty=TRUE)
# add the median of the hierarchical DL curves
get.bias.model(mcmc.set, 'Nigeria')}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ get model }
