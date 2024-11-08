\name{get.tfr.prediction}
\Rdversion{1.1}
\alias{get.tfr.prediction}
\alias{has.tfr.prediction}
\alias{available.tfr.predictions}

\title{
Accessing a Prediction Object
}
\description{
Function \code{get.tfr.prediction} retrieves results of a prediction and creates an object of class \code{\link{bayesTFR.prediction}}. Function \code{has.tfr.prediction} checks an existence of such results. Function \code{available.tfr.predictions} lists predictions available in the given simulation directory.
}
\usage{
get.tfr.prediction(mcmc = NULL, sim.dir = NULL, mcmc.dir = NULL,
                    subdir = "predictions")

has.tfr.prediction(mcmc = NULL, sim.dir = NULL, subdir = "predictions")

available.tfr.predictions(mcmc = NULL, sim.dir = NULL, full.names = FALSE)
}

\arguments{
  \item{mcmc}{Object of class \code{\link{bayesTFR.mcmc.set}} used to make the prediction. It must correspond to a Phase II MCMC. If it is \code{NULL}, the prediction is loaded from directory given by \code{sim.dir}.}
  \item{sim.dir}{Directory where the prediction is stored.  It should correspond to the value of the \code{output.dir} argument used in the \code{\link{tfr.predict}} function.  Only relevant if \code{mcmc} is \code{NULL}.}
  \item{mcmc.dir}{Optional argument to be used only in a special case when the mcmc object contained in the prediction object was estimated in different directory than in the one to which it points to (for example due to moving or renaming the original directory). The argument causes that the mcmc is redirected to the given directory. It can be set to \code{NA} if no loading of the mcmc object is desired.}
  \item{subdir}{Subdirectory of \code{sim.dir} for this particular prediction.}
  \item{full.names}{Logical. If \code{TRUE}, the directory names are given as full paths, otherwise (default) only the base names.}
}
\details{
If \code{mcmc} is not \code{NULL}, the search directory is set to \code{mcmc$meta$output.dir}. This approach assumes that the prediction was stored in the same directory as the MCMC simulation, i.e. the \code{output.dir} argument of the \code{\link{tfr.predict}} function was set to \code{NULL}. If it is not the case, the argument \code{mcmc.dir} should be used.

Usually, all predictions are stored in the subdirectory \dQuote{predictions} of the simulation directory. If the subdirectory has a different name, the argument \code{subdir} should be used. This allows to keep multiple predictions in one (MCMC) simulation directory.

The function \code{available.tfr.predictions} can be used to view all available predictions in the simulation directory.  
}

\value{
	Function \code{has.tfr.prediction} returns a logical indicating if a prediction exists for the given \code{mcmc}. 
	
	Function \code{available.tfr.predictions} returns a vector of directory names containing TFR predictions. 
	
	Function \code{get.tfr.prediction} returns an
object of class \code{\link{bayesTFR.prediction}}.
}

\author{
Hana Sevcikova
}
\seealso{
\code{\link{bayesTFR.prediction}}, \code{\link{tfr.predict}}, \code{\link{summary.bayesTFR.prediction}}
}
\examples{
sim.dir <- file.path(find.package("bayesTFR"), "ex-data", "bayesTFR.output")
pred <- get.tfr.prediction(sim.dir=sim.dir)
summary(pred, country="Canada")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
