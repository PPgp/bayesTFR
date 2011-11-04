\name{UN2010}
\Rdversion{1.1}
\alias{UN2008}
\alias{UN2010}
\alias{my_tfr_template}
\docType{data}
\title{
United Nations Time Series of the Total Fertility Rate
}
\description{
Data set containing the United Nations time series of the total fertility rate for all countries of the world as available in 2008 and 2010.  It is stored in a TAB-separated ASCII files \dQuote{UN2008.txt} and \dQuote{UN2010.txt}.
}
\usage{
data(UN2010)
data(UN2008)
}

\format{
  A data frame with 229 observations, one per country or region. It contains the following variables:
  \describe{
    \item{\code{Index}}{Numeric vector enumerating the observations. Not used by the package.}
    \item{\code{country}}{Name of country or region (following ISO 3166 official short names in English - see \url{http://www.iso.org/iso/country_codes/iso_3166_code_lists/english_country_names_and_code_elements.htm} and United Nations Multilingual Terminology Database - see \url{http://unterm.un.org}).}
    \item{\code{country_code}}{Numerical Location Code (3-digit codes following ISO 3166-1 numeric standard) - see \url{http://en.wikipedia.org/wiki/ISO_3166-1_numeric}.}
    \item{\code{1950-1955}}{TFR in time interval from 1950 to 1955.}
    \item{\code{1955-1960}}{TFR in time interval from 1955 to 1960.}
    \item{\code{1960-1965}}{TFR in time interval from 1960 to 1965.}
    \item{\code{1965-1970}}{TFR in the given time interval }
    \item{\code{1970-1975}}{TFR in the given time interval  }
    \item{\code{1975-1980}}{TFR in the given time interval  }
    \item{\code{1980-1985}}{TFR in the given time interval  }
    \item{\code{1985-1990}}{TFR in the given time interval  }
    \item{\code{1990-1995}}{TFR in the given time interval  }
    \item{\code{1995-2000}}{TFR in the given time interval  }
    \item{\code{2000-2005}}{TFR in the given time interval  }
    \item{\code{2005-2010}}{TFR in the given time interval  }
    \item{\code{2010-2015}}{TFR in the given time interval  }
    \item{\code{2015-2020}}{TFR in the given time interval  }
    \item{\code{2020-2025}}{TFR in the given time interval  }
    \item{\code{2025-2030}}{TFR in the given time interval  }
    \item{\code{2030-2035}}{TFR in the given time interval  }
    \item{\code{2035-2040}}{TFR in the given time interval  }
    \item{\code{2040-2045}}{TFR in the given time interval  }
    \item{\code{2045-2050}}{TFR in the given time interval  }
  }
}
\details{
All entries in \code{UN2008} and \code{UN2010}, respectively, must be included in the data set of locations, i.e. \code{\link{WPP2008_LOCATIONS}} and \code{\link{WPP2010_LOCATIONS}}, respectively. The the corresponding data set pair is matched using the column \code{country_code}.

Regarding the time columns, only those are necessary that correspond to an earlier or equal time point than the value of the argument \code{present.year} passed to the function \code{\link{run.tfr.mcmc}}. Thus, if \code{present.year=2009}, all columns after \code{2005-2010} are not used.

Optionally, the data set can contain a column called \code{last.observed} containing the year of the last observation for each country. In such a case, the code would ignore any data after that time point. Furthermore, the function \code{\link{tfr.predict}} fills in the missing values using the median of the BHM procedure (stored in \code{tfr_matrix_reconstructed} of the \code{\link{bayesTFR.prediction}} object). For \code{last.observed} values that are below a middle year of a time interval \eqn{[t_i, t_{i+1}]} (computed as \eqn{t_i+3}) the last valid data point is the time interval \eqn{[t_{i-1}, t_i]}, whereas for values larger equal a middle year, the data point in \eqn{[t_i, t_{i+1}]} is valid.

The package contains also a data set called \sQuote{my_tfr_template} which is a template for user-specified TFR time series. It has the same structure, except that it can contain any subset of the set of the columns above. The only required column is \code{country_code} (see description of the argument \code{my.tfr.file} in \code{\link{run.tfr.mcmc}}).
}
\source{
Data provided by the United Nations Population Division; see also \url{http://esa.un.org/unpd/wpp}.
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(UN2008)
head(UN2008)
}
\keyword{datasets}
