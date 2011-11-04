\name{tfr.plot.gelman}
\Rdversion{1.1}
\alias{tfr.plot.gelman}
\title{
Plotting Results of tfr.gelman.diag
}
\description{
Plotting results of the \code{\link{tfr.gelman.diag}} function.
}
\usage{
tfr.plot.gelman(data = NULL, 
    file = file.path(getwd(), "bayesTFR.output", 
        "gelman.diag.t0.txt"), 
    show.pars = "slow", nr.pars = 10, quant = 0.95, 
    par.names = tfr.parameter.names(trans = TRUE), 
    min.burnin = 0, max.burnin = NULL, optimalR = 1.1)
}
\arguments{
  \item{data}{
Results of the function \code{\link{tfr.gelman.diag}}. If it is \code{NULL}, it is read from \code{file}. 
}
  \item{file}{
Name of the file written by \code{\link{tfr.gelman.diag}}. It is only used if \code{data} is \code{NULL}.
}
  \item{show.pars}{
If \code{show.pars} is \dQuote{slow}, the function plots \code{nr.pars} parameters with the highest Gelman diagnostics. If the value is \dQuote{quantile}, it plots quantile curves given by \code{quant}. Otherwise, it plots all parameters. If \code{par.names} is given, this argument is ignored.
}
  \item{nr.pars}{
Number of parameters to be plotted. Only used if \code{show.pars} is \dQuote{slow}.
}
  \item{quant}{
Quantile value to be plotted.
}
  \item{par.names}{
If \code{par.names} is given, only those parameters are plotted and \code{show.pars} is ignored.
}
  \item{min.burnin}{
Minimum value of the \dQuote{t0} (\dQuote{t1}) column to be plotted.
}
  \item{max.burnin}{
Maximum value of the \dQuote{t0} (\dQuote{t1}) column to be plotted.
}
  \item{optimalR}{
Value of the Gelman dignostics considered as the optimum.  The function draws a horizontal line at this value.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
