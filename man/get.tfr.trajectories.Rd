\name{get.tfr.trajectories}
\alias{get.tfr.trajectories}

\title{
Accessing TFR Trajectories
}
\description{
Function for accessing TFR trajectories.
}
\usage{
get.tfr.trajectories(tfr.pred, country)
}

\arguments{
  \item{tfr.pred}{Object of class \code{\link{bayesTFR.prediction}}.}
  \item{country}{Name or numerical code of a country.}
}
\details{
The function loads TFR trajectories for the given country from disk, offsets it if needed (see \code{\link{tfr.median.shift}}) and returns it as a matrix.
}
\value{
Array of size number of projection periods (including the present year) times the number of trajectories. The row names correspond to the mid-years of the prediction periods.
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
\code{\link{bayesTFR.prediction}}, \code{\link{get.tfr.prediction}}, \code{\link{tfr.trajectories.table}}, \code{\link{tfr.median.shift}}
}
\examples{
sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", "bayesTFR.output") 
pred <- get.tfr.prediction(sim.dir=sim.dir)
get.tfr.trajectories(pred, "Germany")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
