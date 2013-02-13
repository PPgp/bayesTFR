\name{tfr.diagnose}
\Rdversion{1.1}
\alias{tfr.diagnose}
\alias{tfr3.diagnose}
\alias{bayesTFR.convergence}
\alias{has.mcmc.converged}

\title{
Convergence Diagnostics of TFR Markov Chain Monte Carlo
}
\description{
Functions \code{tfr.diagnose} and \code{tfr3.diagnose} run convergence diagnostics of existing TFR MCMCs for phase II and phase III, respectively, using the \code{raftery.diag} function from the \pkg{coda} package. \code{has.mcmc.converged} checks if the existing diagnostics converged.
}
\usage{
tfr.diagnose(sim.dir, thin = 80, burnin = 2000, express = FALSE, 
    country.sampling.prop = NULL, keep.thin.mcmc=FALSE, verbose = TRUE)
    
tfr3.diagnose(sim.dir, thin = 60, burnin = 10000, express = TRUE, 
    country.sampling.prop = NULL, verbose = TRUE)
    
has.mcmc.converged(diag)
}
\arguments{
  \item{sim.dir}{Directory with the MCMC simulation results.}
  \item{thin}{Thinning interval.}
  \item{burnin}{Number of iterations to be discarded from the beginning of the parameter traces.}
  \item{express}{Logical. If \code{TRUE}, the convergence diagnostics is run only on the country-independent parameters. If \code{FALSE}, the country-specific parameters are included in the diagnostics. The number of countries can be controlled by \code{country.sampling.prop}.}
  \item{country.sampling.prop}{Proportion of countries that are included in the diagnostics. If it is \code{NULL} and \code{express=FALSE}, all countries are included. Setting here a number between 0 and 1, one can limit the number of countries which are then randomly sampled. Note that for long MCMCs, this argument may significantly influence the run-time of this function.}
  \item{keep.thin.mcmc}{Logical. If \code{TRUE} the thinned traces used for computing the diagnostics are stored on disk (see \code{\link{create.thinned.tfr.mcmc}}). It is only available for phase II MCMCs.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{diag}{Object of class \code{bayesTFR.convergence}.}
}
\details{
The diagnose functions invoke the \code{\link{tfr.raftery.diag}} (or \code{\link{tfr3.raftery.diag}}) function separately for country-independent parameters and for country-specific parameters. It results in two possible states: red, i.e. it did not converge, and green, i.e. it converged.
The resulting object from \code{tfr.diagnose} is stored in \cr \file{\{sim.dir\}/diagnostics/bayesTFR.convergence_\{thin\}_\{burnin\}.rda} and can be accessed using the function \code{\link{get.tfr.convergence}}. Function \code{tfr3.diagnose} stores its result into \cr \file{\{sim.dir\}/phaseIII/diagnostics/bayesTFR.convergence_\{thin\}_\{burnin\}.rda} which can be accessed via \code{\link{get.tfr3.convergence}}. 
}
\value{
	\code{has.mcmc.converged} returns a logical value determining if there is convergence or not.
	
	\code{tfr.diagnose} and \code{tfr3.diagnose} return an object of class \code{bayesTFR.convergence} with components:
	\item{result}{Table containing all not-converged parameters. Its columns include \sQuote{Total iterations needed} and \sQuote{Remaining iterations}.}
	\item{lresult.country.independent}{Number of rows in \code{result} that correspond to country-independent paramters. These rows are groupped at the beginning of the table.}
	\item{country.independent}{Result of \code{\link{tfr.raftery.diag}} processed on country-independent parameters.}
	\item{country.specific}{Result of \code{\link{tfr.raftery.diag}} processed on country-specific parameters.}
	\item{iter.needed}{Number of additional iterations suggested in order to achieve convergence.}
	\item{iter.total}{Total number of iterations of the original unthinned set of chains.}
	\item{use.nr.traj}{Suggestion for number of trajectories in generating predictions.}
	\item{burnin}{Burnin used.}
	\item{thin}{Thinning interval used.}
	\item{status}{Vector of character strings containing the result status. Possible values: \sQuote{green}, \sQuote{red}.}
	\item{mcmc.set}{Object of class \code{\link{bayesTFR.mcmc.set}} that corresponds to the original set of MCMCs on which the diagnostics was run.}
	\item{thin.mcmc}{If \code{keep.thin.mcmc} is \code{TRUE}, it is an object of class \code{\link{bayesTFR.mcmc.set}} that corresponds to the thinned mcmc set on which the diagnostics was run, otherwise \code{NULL}.}
	\item{express}{Value of the input argument \code{express}.}
	\item{nr.countries}{Vector with elements \code{used} - number of countries used in this diagnostics, and \code{total} - number of countries that this \code{mcmc.set} object was estimated on.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova, Leontine Alkema, Adrian Raftery
}
%\note{
%%  ~~further notes~~
%}


\seealso{
\code{\link{tfr.raftery.diag}}, \code{\link[coda]{raftery.diag}}, \code{\link{summary.bayesTFR.convergence}}, \code{\link{get.tfr.convergence}}, \code{\link{create.thinned.tfr.mcmc}}
}

\keyword{ htest }

