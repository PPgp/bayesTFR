\name{summary.bayesTFR.mcmc.set}
\Rdversion{1.1}
\alias{summary.bayesTFR.mcmc.set}
\alias{print.summary.bayesTFR.mcmc.set}
\alias{print.summary.bayesTFR.mcmc.set.meta}
\alias{summary.bayesTFR.mcmc}
\title{
Summary Statistics for TFR Markov Chain Monte Carlo Chains
}
\description{
Summary of an object \code{\link{bayesTFR.mcmc.set}} or \code{\link{bayesTFR.mcmc}}, computed via \code{\link{run.tfr.mcmc}}.  It can be obtained either for all countries or for a specific country, and either for all parameters or for specific parameters.  The function uses the \code{\link[coda]{summary.mcmc}} function of the \pkg{coda} package.
}
\usage{
\method{summary}{bayesTFR.mcmc.set}(object, country = NULL, chain.id = NULL, 
    par.names = tfr.parameter.names(trans = TRUE), 
    par.names.cs = tfr.parameter.names.cs(trans = TRUE), 
    meta.only = FALSE, thin = 1, burnin = 0, \dots)
	
\method{summary}{bayesTFR.mcmc}(object, country = NULL,
    par.names = tfr.parameter.names(trans = TRUE), 
    par.names.cs = tfr.parameter.names.cs(trans = TRUE), 
    thin = 1, burnin = 0, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{\link{bayesTFR.mcmc.set}} or \code{\link{bayesTFR.mcmc}}.}
  \item{country}{Country name or code if a country-specific summary is desired. By default, summary for all countries is generated.}
  \item{chain.id}{Identifiers of MCMC chains. By default, all chains are considered.}
  \item{par.names}{Country independent parameters to be included in the summary.}
  \item{par.names.cs}{Country-specific parameters to be included in the summary.}
  \item{meta.only}{If it is \code{TRUE}, only meta information of the simulation is included.}
  \item{thin}{Thinning interval. Only used if larger than the \code{thin} argument used in \code{\link{run.tfr.mcmc}}.}
  \item{burnin}{Number of iterations to be discarded from the beginning of each chain before computing the summary.}
  \item{\dots}{Additional arguments passed to the \code{\link[coda]{summary.mcmc}} function of the \pkg{coda} package.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%}

\author{
Hana Sevcikova
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{bayesTFR.mcmc.set}}, \code{\link[coda]{summary.mcmc}}
}
\examples{
sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
m <- get.tfr.mcmc(sim.dir)
summary(m, country='Czech Republic', burnin=20)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
