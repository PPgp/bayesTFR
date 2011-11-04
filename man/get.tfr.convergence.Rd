\name{get.tfr.convergence}
\Rdversion{1.1}
\alias{get.tfr.convergence}
\alias{get.tfr.convergence.all}

\title{
Accessing a Convergence Object
}
\description{
The function loads objects of class \code{\link{bayesTFR.convergence}} from disk.
}
\usage{
get.tfr.convergence(sim.dir = file.path(getwd(), "bayesTFR.output"), 
    thin=80, burnin = 2000)
	
get.tfr.convergence.all(sim.dir = file.path(getwd(), "bayesTFR.output"))
}
\arguments{
  \item{sim.dir}{Simulation directory used for computing the diagnostics.}
  \item{thin}{Thinning interval used with this diagnostics.}
  \item{burnin}{Burnin used for computing the diagnostics.}
}
\details{
Function \code{get.tfr.convergence} loads an object of class \code{\link{bayesTFR.convergence}} for the specific \code{thin} and \code{burnin}. Function \code{get.tfr.convergence.all} loads all \code{\link{bayesTFR.convergence}} objects available for \code{sim.dir}. 
}
\value{
\code{get.tfr.convergence} returns an object of class \code{\link{bayesTFR.convergence}}; \cr
\code{get.tfr.convergence.all} returns a list of objects of class \code{\link{bayesTFR.convergence}}.
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
\code{\link{bayesTFR.convergence}}, \code{\link{summary.bayesTFR.convergence}}.
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

