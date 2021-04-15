\name{get.thinned.tfr.mcmc}
\Rdversion{1.1}
\alias{get.thinned.tfr.mcmc}
\alias{create.thinned.tfr.mcmc}

\title{
Creating and Accessing Thinned MCMCs
}
\description{
The function \code{get.thinned.tfr.mcmc} accesses 
a thinned and burned version of the given Phase II MCMC set. \code{create.thinned.tfr.mcmc} creates or updates such a set.
}

\usage{
get.thinned.tfr.mcmc(mcmc.set, thin = 1, burnin = 0)

create.thinned.tfr.mcmc(mcmc.set, thin = 1, burnin = 0, 
    output.dir = NULL, verbose = TRUE, uncertainty = FALSE,
    update.with.countries = NULL)
}
\arguments{
  \item{mcmc.set}{Object of class \code{\link{bayesTFR.mcmc.set}} of Phase II.}
  \item{thin, burnin}{Thinning interval and burnin used for creating or identifying the thinned object.}
  \item{output.dir}{Output directory. It is only used if the output goes to a non-standard directory.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{uncertainty}{If users want to save the thinned estimated TFR in the new mcmc object, this parameter should be set \code{TRUE}.}
  \item{update.with.countries}{If an existing set is to be updated, this should be a vector of country indices for the update.}
}
\details{
The function \code{create.thinned.tfr.mcmc} is called from \code{\link{tfr.predict}} and thus, the resulting object contains exactly the same MCMCs used for generating projections. In addition, it can be also called from \code{\link{tfr.diagnose}} if desired, so that the projection process can re-use such a set that leads to a convergence. 

The thinning is done as follows: The given \code{burnin} is removed from the beginning of each chain in the original MCMC set. Then each chain is thinned by \code{thin} using equal spacing and all chains are collapsed into one single chain per parameter. They are stored in the main simulation directory under the name \file{thinned_mcmc_\emph{t}_\emph{b}} where \emph{t} is the value  of \code{thin} and \emph{b} the value of \code{burnin}.   

If \code{uncertainty=TRUE}, the estimated TFR is thinned and saved as well.
}
\value{
Both functions return an object of class \code{\link{bayesTFR.mcmc.set}}. \code{get.thinned.tfr.mcmc} returns \code{NULL} if such object does not exist.
}

\author{
Hana Sevcikova
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bayesTFR.mcmc.set}}, \code{\link{tfr.predict}}, \code{\link{tfr.diagnose}}
}
\examples{
\dontrun{
sim.dir <- tempfile()
m <- run.tfr.mcmc(nr.chains=2, iter=30, seed=1, output.dir=sim.dir, verbose=TRUE)
tfr.predict(m, burnin=15, use.tfr3=FALSE) # creates thinned MCMCs
mb <- get.thinned.tfr.mcmc(m, thin=1, burnin=15)
summary(mb, meta.only=TRUE) # length 30 = 2chains x (30-15)iters.
unlink(sim.dir, recursive=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ programming }
\keyword{ manip }

