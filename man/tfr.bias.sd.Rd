\name{tfr.bias.sd}
\alias{tfr.bias.sd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Accessing estimated bias, standard deviation estimates together.
}
\description{
Function for obtaining bias and standard deviation estimated models simultaneously.
}
\usage{
tfr.bias.sd(mcmc.list=NULL, country.code=NULL, ISO.code=NULL, sim.dir=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmc.list}{Object of class \code{\link{bayesTFR.mcmc.set}} corresponding Phase II MCMCs. If it is \code{NULL}, the object is loaded from the directory given by \code{sim.dir}.}
  \item{country.code}{Name or numerical code of a country. If it is \code{NULL}, \code{ISO.code} need to be provided with ISO-3 character.}
  \item{ISO.code}{Country code as ISO-3 character. Only used if \code{country.code} is set to be \code{NULL}.}
  \item{sim.dir}{Directory with the MCMC simulation results.}
}
\details{
Function \code{tfr.bias.sd} is used to obtain bias and standard deviation estimated models together. They are used in the MCMC steps stored in \code{mcmc.list}.
}
\value{
  \item{model_bias}{\code{lm} object describe the linear model estimated for bias model.}
  \item{model_sd}{\code{lm} object describe the linear model estimated for standard deviation model.}
  \item{table}{\code{data.frame} object storing the bias and standard deviation of all possible combinations in the raw data sets for a certain country.}
}
\author{
Peiran Liu, Hana Sevcikova
}
\examples{
\dontrun{
mcmc.set <- run.tfr.mcmc(nr.chains = 1, iter = 10, output.dir = "test.dir", 
    replace.output = TRUE, uncertainty = TRUE)
# add the median of the hierarchical DL curves
get.bias.model(mcmc.set, 'Nigeria')}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ get model }
