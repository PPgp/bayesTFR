\name{bayesTFR.mcmc}
\Rdversion{1.1}
\alias{bayesTFR.mcmc}

\title{
MCMC Simulation Object
}
\description{
MCMC simulation object \code{bayesTFR.mcmc} containing information about one MCMC chain, either from Phase II or Phase III simulation. A set of such objects belonging to the same simulation together with a \code{\link{bayesTFR.mcmc.meta}} object constitute a \code{\link{bayesTFR.mcmc.set}} object. 
}

\value{
A \code{bayesTFR.mcmc} object contains parameters of the Bayesian hierarchical model, more specifically, their values from the last iteration.  If it is a \emph{Phase II} object these parameters are: \cr
\code{psi, chi, a_sd, b_sd, const_sd, S_sd, sigma0, mean_eps_tau, sd_eps_tau, Triangle4} - non-country specific parameters, containing one value each.\cr
\code{alpha, delta} -  non-country specific parameters, containing three values each.\cr
\code{U_c, d_c, Triangle_c4} - country-specific parameters (1d array).\cr
\code{gamma_ci} - country-specific parameter with three values for each country, i.e. an \eqn{n \times 3}{n x 3} matrix where \eqn{n} is the number of countries.\cr
\emph{Phase III} MCMC objects contain single-value parameters \code{mu}, \code{rho}, \code{sigma.mu}, \code{sigma.rho}, \code{sigma.eps} and \eqn{n}-sized vectors \code{mu.c} and \code{rho.c}.\cr
Furthermore, the object (independent of Phase) contains components:
\item{iter}{Total number of iterations the simulation was started with.}
\item{finished.iter}{Number of iterations that were finished.  Results from the last finished iteration are stored in the parameters above.}
\item{length}{Length of the MCMC stored on disk. It differs from \code{finished.iter} only if \code{thin} is larger than one.}
\item{thin}{Thinning interval used when simulating the MCMCs.}
\item{id}{Identifier of this chain.}
\item{output.dir}{Subdirectory (relative to \code{output.dir} in the \code{\link{bayesTFR.mcmc.meta}} object) where results of this chain are stored.}
\item{traces}{This is a placeholder for keeping whole parameter traces in the memory.  If the processing operates in a low memory mode, it will be 0.  It can be filled in using the function \code{\link{get.tfr.mcmc}(\dots, low.memory=FALSE)}. In such a case, \code{traces} is a \eqn{I \times J}{I x J} array where \eqn{I} is the MCMC \code{length} and \eqn{J} is the number of parameters.}
\item{traces.burnin}{Burnin used to retrieve the traces, i.e. how many stored iterations are missing from the beginning in the \code{traces} array comparing to the \sQuote{raw} traces on the disk.}
\item{rng.state}{State of the random number generator at the end of the last finished interation.}
\item{compression.type}{Type of compression of the underlying files.}
\item{meta}{Object of class \code{\link{bayesTFR.mcmc.meta}} used for simulation of this chain.}
}
\details{
An object \code{bayesTFR.mcmc} points to a place on disk (element \code{output.dir}) where MCMC results from all iterations are stored. They can be retrieved to the memory using \code{\link{get.tfr.mcmc}(\dots)} (Phase II) or \code{\link{get.tfr3.mcmc}(\dots)} (Phase III), and \code{\link{tfr.mcmc}(\dots)}. 

The object is in standard cases not to be manipulated by itself, but rather as part of a \code{\link{bayesTFR.mcmc.set}} object.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{run.tfr.mcmc}}, \code{\link{get.tfr.mcmc}}, \code{\link{run.tfr3.mcmc}}, \code{\link{get.tfr3.mcmc}}, \code{\link{bayesTFR.mcmc.set}}, \code{\link{bayesTFR.mcmc.meta}}
}
\examples{
sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", "bayesTFR.output")
# loads traces from one chain
m <- get.tfr.mcmc(sim.dir, low.memory=FALSE, burnin=35, chain.ids=1)
# should have 25 rows, since 60 iterations in total minus 35 burnin
dim(tfr.mcmc(m, 1)$traces)
summary(m, chain.id=1)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classes}

