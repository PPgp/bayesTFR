\name{UN_time}
\Rdversion{1.1}
\alias{UN_time}
\docType{data}
\title{
Dataset with UN-specific Time Coding
}
\description{
Dataset used by the UN for coding time. It is an TAB-separated ASCII file called \dQuote{UN_time.txt}.
}
\usage{data(UN_time)}
\format{
  A data frame with 1034 observations on the following 4 variables.
  \describe{
    \item{\code{TimeID}}{Time identifier.}
    \item{\code{TLabel}}{Label of the time, with minimum values of \code{1950} and \code{1950-1955}, and maximum values of   \code{2399}, \code{2400} and \code{2400-2405}.}
    \item{\code{TDate}}{Equal to \code{TLabel} if it is a single year, or the starting year of \code{TLabel}, if it is an interval.}
    \item{\code{Tinterval}}{Length of the time interval, or zero, if it is a single year.}
  }
}
\details{
For 5-year period data, fertility rates are defined from 1 July of year (t) to 1 July of year (t+5), with 1 January of year (t+3) as exact mid-date. This means for example that data for 2000-2005, refer to the period between 2000.5 and 2005.5, with 2003.0 as exact mid-point.
}
\source{
Data provided by the United Nations Population Division
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(UN_time)
str(UN_time)
}
\keyword{datasets}
