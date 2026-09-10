\name{tfr.traj.shift}
\alias{tfr.median.set}
\alias{tfr.mean.set}
\alias{tfr.median.shift}
\alias{tfr.traj.shift}
\alias{tfr.median.adjust}
\alias{tfr.mean.adjust}
\alias{tfr.median.reset}
\alias{tfr.shift.reset}
\alias{tfr.shift.prediction.to.wpp}

\title{
Shifting the Distribution of the Projection
}
\description{
These functions are to be used by expert analysts. They allow to change the projection medians/means either to specific values, including the WPP values, or shift the distribution by a given constant, or by a specific adjusting procedure.
}
\usage{
tfr.median.set(sim.dir, country, values, years = NULL, \dots)

tfr.mean.set(sim.dir, country, values, years = NULL, \dots)

tfr.traj.shift(sim.dir, country, reset = FALSE, shift = 0, 
    from = NULL, to = NULL, \dots)
    
tfr.median.shift(\dots) # obsolete

tfr.median.adjust(sim.dir, countries, factor1 = 2/3, factor2 = 1/3, forceAR1 = FALSE,
    subdir = "predictions")
    
tfr.mean.adjust(sim.dir, countries, factor1 = 2/3, factor2 = 1/3, forceAR1 = FALSE,
    subdir = "predictions")

tfr.shift.reset(sim.dir, countries = NULL, \dots)

tfr.median.reset(\dots) # obsolete

tfr.shift.prediction.to.wpp(sim.dir, subdir = "predictions", \dots)
}
\arguments{
  \item{sim.dir}{Directory containing the prediction object.}
  \item{country}{Name or numerical code of a country.}
  \item{countries}{Vector of country names or codes. If \code{NULL} in the \code{tfr.median.reset} function, the reset is done for all countries.}
  \item{values}{Array of the new median values.}
  \item{years}{Numeric vector giving years which \code{values} correspond to. Ideally it should be of the same length as \code{values}. If it is \code{NULL},
  	\code{values} are set starting from the first prediction period. If \code{values} correspond to consecutive years, only the first year might be given here. A year \eqn{t} represents a prediction period \eqn{[t_i, t_{i+1}]} if \eqn{t_i < t \leq t_{i+1}}.}
  \item{reset}{Logical. If \code{TRUE} medians in a range of \code{from} and \code{to} are reset to their original values.}
  \item{shift}{Constant by which the medians should be offset. It is not used if \code{reset} is \code{TRUE}.}
  \item{from}{Year from which the offset/reset should start. By default, it starts at the first prediction period.}
  \item{to}{Year until which the offset/reset should be done. By default, it is set to the last prediction period.}
  \item{factor1, factor2}{Adjusting factors for the first and second projection period, respectively (see below).}
  \item{forceAR1}{Logical. If \code{TRUE}, the given countries are forced to enter Phase III (i.e. the AR(1) process) in the first projection period.}
  \item{subdir}{Subdirectory of \code{sim.dir} containing the predictions.}
  \item{\dots}{Additional arguments passed to the underlying adjustment function. For \code{tfr.shift.prediction.to.wpp} it can be \code{stat} with values \dQuote{median} (default) or \dQuote{mean} to specify which statistics should be adjusted; \code{verbose} to show/hide the progress of the adjustment and \code{wpp.year} to adjust it to if it differs from the wpp year of the simulation. For the other functions it can be \code{subdir} to specify the location of the prediction.}
}
\details{
The functions \code{tfr.median.set} and \code{tfr.mean.set} can be used to set the medians/means of the given country to specific values. Function \code{tfr.traj.shift} can be used to offset the trajectories by a specific constant, or to reset the distributions to their original BHM values. The function \code{tfr.median.shift} is obsolete, replaced by \code{tfr.traj.shift}.
Functions \code{tfr.median.adjust} and \code{tfr.mean.adjust} run the prediction procedure for the given countries with an additional decrement in the model in the first two projection periods. In the first projection period it is computed as \code{factor1*S} where \code{S} is a difference between observed decrement and the expected decrement (by the double logistic function) in the last observed period. In the second projection period, in the formula \code{factor1} is replaced by \code{factor2}. If \code{forceAR1} is set to \code{TRUE}, we recommend to set \code{factor1} and \code{factor2} to 0. The function then calls \code{tfr.median.set} or \code{tfr.mean.set} in order to store the new median/mean for each country.

Function\code{tfr.shift.prediction.to.wpp} shifts the projected medians or means (if \code{stat} is \dQuote{mean}), so that they correspond to the values found in the \code{tfrprojMed} datasets of the \pkg{wpp} package that either corresponds to the package used for the simulation itself or is given by the \code{wpp.year} argument. If using \pkg{wpp2022} or higher, the dataset name is automatically adjusted depending if it is an annual or a 5-year simulation. 

Function \code{tfr.shift.reset} resets the distribution of the given countries to the original values. By default it deletes adjustments for all countries. Function \code{tfr.median.reset} is obsolete, replaced by \code{tfr.shift.reset}.

In all these functions, if a median or mean is modified, the corresponding offset is stored in the prediction object (element \code{traj.shift}) and the updated prediction object is written back to disk. All functions in the package that use trajectories and trajectory statistics use the \code{traj.shift} values to offset the results correspondingly by shifting the trajectories.
}
\value{
All functions return an updated object of class \code{\link{bayesTFR.prediction}}.
}

\author{Hana Sevcikova, Leontine Alkema}

\seealso{\code{\link{tfr.median.set.all}} and \code{\link{tfr.mean.set.all}} for shifting estimation medians and means.}

%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

