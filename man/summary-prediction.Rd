\name{summary.bayesTFR.prediction}
\Rdversion{1.1}
\alias{summary.bayesTFR.prediction}
\alias{print.summary.bayesTFR.prediction}

\title{
Summary of a Prediction of the Total Fertility Rate
}
\description{
Country-specific summary of an object of class \code{\link{bayesTFR.prediction}}, created using the function \code{\link{tfr.predict}}. The summary contains the mean, standard deviation and several commonly used quantiles of the simulated trajectories.
}
\usage{
\method{summary}{bayesTFR.prediction}(object, country=NULL, compact = TRUE, \dots)
}
\arguments{
  \item{object}{Object of class \code{\link{bayesTFR.prediction}}.}
  \item{country}{Country name or code. If it is \code{NULL}, only prediction parameters are included.}
  \item{compact}{Logical switching between a smaller and larger number of displayed quantiles.}
  \item{\dots}{A list of further arguments.}
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
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{bayesTFR.prediction}}
}
\examples{
\dontrun{
sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
pred <- tfr.predict(sim.dir=sim.dir, 
                    output.dir=file.path(getwd(), 'exampleTFRpred'), 
                    burnin=25, verbose=TRUE)
# If the above function was run previously, do
# pred <- get.tfr.prediction(sim.dir=file.path(getwd(), 'exampleTFRpred'))
                                                        
summary(pred, country='Ireland')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
\keyword{ print }
