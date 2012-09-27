\name{tfr.dl.coverage}
\alias{tfr.dl.coverage}
\title{
Goodness of Fit of the Double Logistic Function
}
\description{
The function computes coverage, i.e. the ratio of observed data fitted within the given probability intervals of the predictive posterior distribution of the double logistic function, as well as the root mean square error of the simulation.
}
\usage{
tfr.dl.coverage(sim.dir, pi = c(80, 90, 95), burnin = 2000, verbose = TRUE)
}

\arguments{
  \item{sim.dir}{Directory with the MCMC simulation results. If a prediction and its corresponding thinned mcmcs are available in the simulation directory, those are taken for assessing the goodness of fit.}
  \item{pi}{Probability interval. It can be a single number or an array.}
  \item{burnin}{Burnin. Only relevant if \code{sim.dir} does not contain thinned chains.}
  \item{verbose}{Logical switching log messages on and off.}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{List with the following components:
  \item{total.coverage}{Vector of the coverage, one element per probability interval. For each \code{pi}, it is the ratio of the number of observed data points that fall within the probability interval of the posterior distribution over the total number of data points, i.e. TFR for all countries and historical time periods.}
  \item{time.coverage}{Matrix corresponding to the coverage computed per time period. (Rows correspond to probability intervals, columns correspond to time.) It is derived as \code{total.coverage} but both, the nominator and denominator contain only data points belonging to the corresponding time period.}
  \item{country.coverage}{Matrix corresponding to the coverage computed per country. (Rows correspond to probability intervals, columns correspond to countries.) It is derived as \code{total.coverage} but both, the nominator and denominator contain only data points belonging to the corresponding country.}
  \item{total.rmse}{Root mean square error as \eqn{\sqrt(1/n\sum(x-m)^2)} where \eqn{x} are observed data points, \eqn{m} is the median of the posterior distribution and \eqn{n} is the number of data points. Here the sum is taken over all countries and historical time periods.}
  \item{time.rmse}{As \code{total.rmse} but each time period is considered separately.}
  \item{country.rmse}{As \code{total.rmse} but each country is considered separately.} 
  \item{n}{0-1 matrix of time periods times countries, which indicates if the corresponding data point was included in the goodness of fit computation. Zeros indicate missing historical values.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova
}
\note{
To see the fit visually per country, use \code{\link{DLcurve.plot}(\dots, predictive.distr=TRUE,\dots)}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{DLcurve.plot}}
}
\examples{
sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", "bayesTFR.output")
tfr <- get.tfr.mcmc(sim.dir)
# Note that this simulation is a toy example and thus has not converged.
gof <- tfr.dl.coverage(sim.dir)
gof$time.coverage
DLcurve.plot(tfr, country=608, predictive.distr=TRUE, pi=c(80, 90, 95))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }

