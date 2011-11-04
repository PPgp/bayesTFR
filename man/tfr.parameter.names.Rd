\name{tfr.parameter.names}
\Rdversion{1.1}
\alias{tfr.parameter.names}
\alias{tfr.parameter.names.cs}
\alias{tfr.parameter.names.extended}
\alias{tfr.parameter.names.cs.extended}

\title{
Accessing Parameter Names
}
\description{
Functions for accessing names of the MCMC parameters, either country-independent or country-specific.
}
\usage{
tfr.parameter.names(trans = NULL)
tfr.parameter.names.cs(country.code = NULL, trans = NULL)
tfr.parameter.names.extended()
tfr.parameter.names.cs.extended(country.code = NULL)
}
\arguments{
  \item{trans}{It can be \code{NULL} or logical. If \code{TRUE},
	 names of the transformable parameters (i.e. \sQuote{alpha} and \sQuote{delta} in case of country-independent parameters, or \sQuote{gamma} in case of country-specific parameters) are replaced by the names of the transformed parameters (i.e. \sQuote{alphat}, \sQuote{deltat}, or \sQuote{gammat}). If \code{trans=FALSE}, there is no such replacement. If \code{trans=NULL}, all parameter names, 
including the transformable parameters are returned.}
\item{country.code}{Country code. If it is given, the country-specific parameter names contain the postfix \sQuote{_c\eqn{x}} where \eqn{x} is the \code{country.code}.}
}

\value{
	\code{tfr.parameter.names} returns names of the country-independent parameters.\cr
	\code{tfr.parameter.names.cs} returns names of the country-specific parameters.\cr
	\code{tfr.parameter.names.extended} returns names of all country-independent parameters, including the transformed parameters. Parameters \sQuote{alpha}, \sQuote{delta}, \sQuote{alphat}, and\sQuote{deltat} are in their extended format with the postfix \sQuote{_1}, \sQuote{_2} and \sQuote{_3}. \cr
	\code{tfr.parameter.names.cs.extended} returns names of all country-specific parameters, including the transformed parameters. Parameters \sQuote{gamma} and\sQuote{gammat} are in their extended format with the postfix \sQuote{_1}, \sQuote{_2} and \sQuote{_3}.
}

\author{
Hana Sevcikova
}

\examples{
tfr.parameter.names()
tfr.parameter.names.extended()
tfr.parameter.names.cs()
tfr.parameter.names.cs.extended()
}

\keyword{ manip }

