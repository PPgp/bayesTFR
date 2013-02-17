\name{tfr.median.set}
\alias{tfr.median.set}
\alias{tfr.median.shift}
\alias{tfr.median.adjust}

\title{
Editing Medians of the Projection
}
\description{
These functions are to be used by expert analysts. They allow to change the projection medians either to specific values or shift the medians by a given constant, or by a specific adjusting procedure.
}
\usage{
tfr.median.set(sim.dir, country, values, years = NULL)

tfr.median.shift(sim.dir, country, reset = FALSE, shift = 0, 
    from = NULL, to = NULL)
    
tfr.median.adjust(sim.dir, countries, factor1 = 2/3, factor2 = 1/3, forceAR1 = FALSE)
}
\arguments{
  \item{sim.dir}{Directory containing the prediction object.}
  \item{country}{Name or numerical code of a country.}
  \item{countries}{Vector of country names or codes.}
  \item{values}{Array of the new median values.}
  \item{years}{Numeric vector giving years which \code{values} correspond to. Ideally it should be of the same length as \code{values}. If it is \code{NULL},
  	\code{values} are set starting from the first prediction period. If \code{values} correspond to consecutive years, only the first year might be given here. A year \eqn{t} represents a prediction period \eqn{[t_i, t_{i+1}]} if \eqn{t_i < t \leq t_{i+1}}.}
  \item{reset}{Logical. If \code{TRUE} medians in a range of \code{from} and \code{to} are reset to their original values.}
  \item{shift}{Constant by which the medians should be offset. It is not used if \code{reset} is \code{TRUE}.}
  \item{from}{Year from which the offset/reset should start. By default, it starts at the first prediction period.}
  \item{to}{Year until which the offset/reset should be done. By default, it is set to the last prediction period.}
  \item{factor1, factor2}{Adjusting factors for the first and second projection period, respectively (see below).}
  \item{forceAR1}{Logical. If \code{TRUE}, the given countries are forced to enter Phase III (i.e. the AR(1) process) in the first projection period.}
}
\details{
The function \code{tfr.median.set} can be used to set the medians of the given country to specific values. Function \code{tfr.median.shift} can be used to offset the medians by a specific constant, or to reset the medians to their original BHM values.
Function \code{tfr.median.adjust} runs the prediction procedure for the given countries with an additional decrement in the model in the first two projection periods. In the first projection period it is computed as \code{factor1*S} where \code{S} is a difference between observed decrement and the expected decrement (by the double logistic function) in the last observed period. In the second projection period, in the formula \code{factor1} is replaced by \code{factor2}. If \code{forceAR1} is set to \code{TRUE}, we recommend to set \code{factor1} and \code{factor2} to 0. The function then calls \code{tfr.median.set} in order to store the new median for each country.

In all three functions, if a median is modified, the corresponding offset is stored in the prediction object (element \code{median.shift}) and the updated prediction object is written back to disk. All functions in the package that use trajectories and trajectory statistics use the \code{median.shift} values to offset the results correspondingly.  
}
\value{
All three functions return an updated object of class \code{\link{bayesTFR.prediction}}.
}

\author{Hana Sevcikova, Leontine Alkema}

%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

