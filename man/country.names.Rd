\name{country.names}
\alias{country.names}

\title{
Accessing Country Names
}
\description{
The function returns country names for countries given either by their codes or by index.
}
\usage{
country.names(meta, countries = NULL, index = FALSE)
}
\arguments{
  \item{meta}{Object of class \code{\link{bayesTFR.mcmc.meta}}, \code{\link{bayesTFR.mcmc.set}}, \code{\link{bayesTFR.mcmc}}, or \code{\link{bayesTFR.prediction}}.}
  \item{countries}{Vector of country codes or indices. If it is not given, names of all countries are returned.}
  \item{index}{Logical indicating if the argument \code{countries} is an index.}
}
\details{
The function considers countries that are included in the simulations and predictions. If the argument \code{countries} is not given, all countries are returned in the same order as they are stored in the meta object. 
}
\value{
Vector of country names.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{get.countries.table}}, \code{\link{get.country.object}}
}
\examples{
sim.dir <- file.path(find.package("bayesTFR"), "ex-data", "bayesTFR.output")
m <- get.tfr.mcmc(sim.dir)
country.names(m)
# these two calls should give the same answer
country.names(m, c(800, 120))
country.names(m, c(15, 20), index=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ attribute }

