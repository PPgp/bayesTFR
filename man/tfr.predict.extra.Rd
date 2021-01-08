\name{tfr.predict.extra}
\Rdversion{1.1}
\alias{tfr.predict.extra}

\title{
Generating Posterior Trajectories of the Total Fertility Rate for Specific Countries or Regions
}
\description{
Using the posterior parameter samples the function  generates posterior trajectories of the total fertility rate for given countries or regions. It is intended to be used after running \code{\link{run.tfr.mcmc.extra}}, but it can be also used for purposes of testing specific settings on one or a few countries. 
}
\usage{
tfr.predict.extra(sim.dir = file.path(getwd(), 'bayesTFR.output'), 
    prediction.dir = sim.dir, countries = NULL, 
    save.as.ascii = 0, verbose = TRUE, uncertainty=FALSE)
	}
\arguments{
  \item{sim.dir}{Directory with the MCMC simulation results.}
  \item{prediction.dir}{Directory where the prediction object and the trajectories are stored.}
  \item{countries}{Vector of country codes for which the prediction should be made. If it is \code{NULL}, the prediction is run for all countries that are included in the MCMC object but for which no prediction was generated.}
  \item{save.as.ascii}{Either a number determining how many trajectories should be converted into an ascii file, or \dQuote{all} in which case all trajectories are converted. It should be set to 0, if no converions is desired. Note that the convertion is done on all countries.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{uncertainty}{Logical. If the MCMC steps has considered uncertainty of past TFR and \code{uncertainty=TRUE}, starting point of prediction trajectories will be the last estimated trajectories of TFR. Otherwise, it will use last observed TFR as starting point of prediction.}
  }
\details{
In order to use this function, a prediction object must exist, i.e. the function \code{\link{tfr.predict}} must have been processed prior to using this function.

Trajectories for given countries or regions are generated and stored in binary format along with other countries (in \code{prediction_dir}). The existing prediction object is updated and stored in the same directory. If \code{save.as.ascii} is larger than zero, trajectories of ALL countries are converted to an ascii format. 
}
\value{
	Updated object of class \code{\link{bayesTFR.prediction}}.
	}
	
\author{
Hana Sevcikova
}
\seealso{
\code{\link{tfr.predict}}
}
%\examples{

%}
\keyword{distribution}
\keyword{multivariate}
