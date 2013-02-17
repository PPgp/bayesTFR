\name{tfr.dl.coverage}
\alias{tfr.dl.coverage}
\title{
Goodness of Fit of the Double Logistic Function
}
\description{
The function computes coverage, i.e. the ratio of observed data fitted within the given probability intervals of the predictive posterior distribution of the double logistic function, as well as the root mean square error and mean absolute error of the simulation.
}
\usage{
tfr.dl.coverage(sim.dir, pi = c(80, 90, 95), burnin = 2000, verbose = TRUE)
}

\arguments{
  \item{sim.dir}{Directory with the MCMC simulation results. If a prediction and its corresponding thinned MCMCs are available in the simulation directory, those are taken for assessing the goodness of fit.}
  \item{pi}{Probability interval. It can be a single number or an array.}
  \item{burnin}{Burnin. Only relevant if \code{sim.dir} does not contain thinned chains.}
  \item{verbose}{Logical switching log messages on and off.}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{List with the following components:
  \item{total.coverage}{Vector of the coverage, one element per probability interval. For each \code{pi}, it is the ratio of the number of observed data points that fall within the probability interval of the posterior distribution over the total number of data points, i.e. TFR for all countries and historical time periods.}
  \item{time.coverage}{Matrix corresponding to the coverage computed per time period. (Rows correspond to probability intervals, columns correspond to time.) It is derived like \code{total.coverage} except that both, the nominator and denominator, contain only data points belonging to the corresponding time period.}
  \item{country.coverage}{Matrix corresponding to the coverage computed per country. (Rows correspond to probability intervals, columns correspond to countries.) It is derived like \code{total.coverage} except that both, the nominator and denominator, contain only data points belonging to the corresponding country.}
  \item{total.rmse}{Root mean square error as \eqn{\sqrt{(1/n\sum(x-m)^2)}} where \eqn{x} are observed data points, \eqn{m} is the mean of the posterior distribution and \eqn{n} is the number of data points. Here the sum is taken over all countries and historical time periods.}
  \item{time.rmse}{Like \code{total.rmse} except that each time period is considered separately.}
  \item{country.rmse}{Like \code{total.rmse} except that each country is considered separately.}
  \item{total.mae}{Mean absolute error as \eqn{1/n\sum|x-m|} where \eqn{x} are observed data points, \eqn{m} is the median of the posterior distribution and \eqn{n} is the number of data points. Here the sum is taken over all countries and historical time periods.}
    \item{time.mae}{Like \code{total.mae} except that each time period is considered separately.}
  \item{country.mae}{Like \code{total.mae} except that each country is considered separately.}
  \item{pred.cdf}{\eqn{T \times C} matrix (with \eqn{T} being the number of time periods and \eqn{C} being the number of countries), containing the predictive CDF of the observation, i.e. the quantile of each data point within the predictive posterior distribution.}
  \item{n}{0-1 \eqn{T \times C} matrix indicating if the corresponding data point was included in the goodness of fit computation. Zeros indicate missing historical values.}
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
\dontrun{
sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", "bayesTFR.output")
tfr <- get.tfr.mcmc(sim.dir)
# Note that this simulation is a toy example and thus has not converged.
gof <- tfr.dl.coverage(sim.dir)
gof$time.coverage
DLcurve.plot(tfr, country=608, predictive.distr=TRUE, pi=c(80, 90, 95))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }

