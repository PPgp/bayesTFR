\name{tfr.raftery.diag}
\Rdversion{1.1}
\alias{tfr.raftery.diag}

\title{
Raftery Diagnostics for Parameters of the Total Fertility Rate
}
\description{
The function computes the Raftery diagnostics for each parameter, taking median over chains.}
\usage{
tfr.raftery.diag(mcmc = NULL, 
    sim.dir = file.path(getwd(), "bayesTFR.output"), 
    burnin = 0, country = NULL,
    par.names = tfr.parameter.names(trans = TRUE), 
    par.names.cs = tfr.parameter.names.cs(trans=TRUE),
    country.sampling.prop = 1, verbose=TRUE, \dots)
}

\arguments{
  \item{mcmc}{
A \code{bayesTFR.mcmc} or \code{bayesTFR.mcmc.set} object.
}
  \item{sim.dir}{
Directory with the MCMC simulation results. Only used if \code{mcmc} is \code{NULL}.
}
  \item{burnin}{
Burnin.
}
  \item{country}{Name or code of a country. If it is given, country-specific parameters are reduced to parameters of that country.}
  \item{par.names}{
Names of country-independent parameters for which the Raftery diagnostics should be computed.
}
  \item{par.names.cs}{
Names of country-specific parameters for which the Raftery diagnostics should be computed.
}
\item{country.sampling.prop}{Proportion of countries that are included in the diagnostics. It should be between 0 and 1. If it is smaller than 1, the countries are randomly sampled. It is only relevant if \code{par.names.cs} is not \code{NULL}. }
\item{verbose}{Logical switching log messages on and off.}
\item{\dots}{Additional arguments passed to the \code{\link{coda.list.mcmc}} function.}
}
\details{
The Raftery diagnostics is computed for each parameter, using \pkg{coda}'s \code{raftery.diag} with r=0.0125, q=0.025 and q=0.975. Values of \eqn{N} and burnin are taken as the median over chains. For each country-specific parameter, the maximum over all included countries of such medians is taken. 
}
\value{
	List with the components:
	\item{Nmedian}{2-d array of \eqn{N}  values (processed as described in Details) with two rows: first corresponding to q=0.025, second corresponding to q=0.975. Each column corresponds to one parameter.}
	\item{burnin}{2-d array of the same structure as \code{Nmedian}, containing the burnin values (processed as described in Details).}
	\item{not.converged.parameters}{List with two elements, each of which is a data frame containing columns \dQuote{parameter.name}, \dQuote{chain.id}, and \dQuote{N}. These are parameters for which the computed value of Raftery diagnostics \eqn{N} is larger than the total number of finished iterations summed over all chains. The first element of the list corresponds to q=0.025, second corresponds to q=0.975.}
	\item{not.converged.inchain.parameters}{List of the same structure as \code{not.converged.parameters}. The parameters included are those for which the computed value of Raftery diagnostics \eqn{N} is larger than the number of finished iterations in the corresponding chain.}
	\item{N.country.indep}{Data frame containing columns \dQuote{parameter.name}, \dQuote{chain.id}, \dQuote{N0.025}, and \dQuote{N0.975}. Each row gives \eqn{N} computed with the two different \eqn{q} for each country-independent parameter and chain.}
	\item{N.country.spec}{The same as \code{N.country.indep}, but here the country-specific parameters are considered.}
	\item{Nmedian.country.spec}{2-d array of \eqn{N} values for country-specific parameters containing medians over chains.}
	\item{thin.ind}{List with elements '0.025', '0.975' and \code{median}. The first two elements are matrices with one row per chain and one column per parameter. They contain values of \code{thin} that makes the MCMC independent, for  q=0.025 and q=0.975, respectively. The \code{median} element is of the same structure as \code{Nmedian}, containing medians ove rows in the two matrices in this list.}
	\item{nr.countries}{Vector with elements \code{used} (number of countries used in this diagnostics) and \code{total} (number of countries that this \code{mcmc} object was estimated on).}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova, Adrian Raftery
}

\seealso{
\code{\link[coda]{raftery.diag}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
