\name{tfr.gelman.diag}
\Rdversion{1.1}
\alias{tfr.gelman.diag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gelman Diagnostics for Multiple Burnins
}
\description{
The function computes the Gelman diagnostics for given parameters in the interval \eqn{[t_0, t_1]} where either \eqn{t_0} or \eqn{t_1} can be moved by a given step. 
}
\usage{
tfr.gelman.diag(sim.dir = file.path(getwd(), "bayesTFR.output"), 
    t.start = 0, t.end = NULL, step = NULL, move.t0 = TRUE, 
    par.names = tfr.parameter.names(trans = TRUE), 
    par.names.cs = tfr.parameter.names.cs(trans = TRUE), 
    file.name = paste("gelman.diag.", 
        ifelse(move.t0, "t0", "t1"), ".txt", sep = ""), 
    append = FALSE, verbose = TRUE)
}
\arguments{
  \item{sim.dir}{
Directory with the MCMC simulation results.
}
  \item{t.start}{
Minimum value for the moving parameter \eqn{t}.
}
  \item{t.end}{
Maximum value for the moving parameter \eqn{t}.  By default it is \eqn{N}, the number of iterations of the MCMC.
}
  \item{step}{
Step size for the moving parameter \eqn{t}.  By default the step size is set to a value that corresponds to having 10 values of the Gelman diagnostics equally spaced between \code{t.start} and \code{t.end}.
}
  \item{move.t0}{
Logical indicating if \eqn{t_0} should be moved forward and \eqn{t_1} held fixed. If \code{FALSE}, \eqn{t_1} is moved backward and \eqn{t_0} is kept fixed (see details below).
}
  \item{par.names}{
Names of country-independent parameters for which the Gelman diagnostics should be computed.
}
  \item{par.names.cs}{
Names of country-specific parameters for which the Gelman diagnostics should be computed.
}
  \item{file.name}{
Name of file into which results are written.
}
  \item{append}{
Logical indicating if the results should be appended to existing data in \code{file.name}. 
}
  \item{verbose}{
Logical switching log messages on and off.}
}
\details{
If \code{move.t0} is \code{TRUE}, the Gelman diagnostics is computed in the interval \code{[t_0(i), N]}
where \code{t_0(i)=t.start+(i-1)*step, i=1,2,..., t_0(i) <= t.end}. 
If \code{move.t0} is \code{FALSE}, the Gelman diagnostics is computed in the interval \code{[0, t_1(i)]}
where \code{t_1(i)=i*step, i=1,2,..., t.start <= t_1(i) <= t.end}.

The first column of the resulting file is called \dQuote{t0}, if  \code{move.t0} is \code{TRUE}, or \dQuote{t1}, if \code{move.t0} is \code{FALSE}, and it contains values of the moving \eqn{t}.  The remaining columns contain the Gelman diagnostics for each parameter. 

Results can be plotted using \code{\link{tfr.plot.gelman}}.
}
\value{
Table read from from \code{file.name}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova, Adrian Raftery
}

\seealso{
\code{\link{tfr.plot.gelman}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
