\name{tfr.parameter.names}
\Rdversion{1.1}
\alias{tfr.parameter.names}
\alias{tfr.parameter.names.cs}
\alias{tfr.parameter.names.extended}
\alias{tfr.parameter.names.cs.extended}
\alias{tfr3.parameter.names}
\alias{tfr3.parameter.names.cs}

\title{
Accessing Parameter Names
}
\description{
Functions for accessing names of the MCMC parameters, either country-independent or country-specific.
}
\usage{
tfr.parameter.names(trans = NULL, meta = NULL)
tfr.parameter.names.cs(country.code = NULL, trans = NULL, back.trans = TRUE)
tfr.parameter.names.extended()
tfr.parameter.names.cs.extended(country.code = NULL)

tfr3.parameter.names()
tfr3.parameter.names.cs(country.code = NULL)
}
\arguments{
  \item{trans}{It can be \code{NULL} or logical. If \code{TRUE},
	 names of the transformable parameters (i.e. \sQuote{alpha} in case of country-independent parameters, or \sQuote{gamma} in case of country-specific parameters) are replaced by the names of the transformed parameters (i.e. \sQuote{alphat}, or \sQuote{gammat}). If \code{trans=FALSE}, there is no such replacement. If \code{trans=NULL}, all parameter names, 
including the transformable parameters are returned.}
\item{meta}{It can be \code{NULL} or bayesTFR.mcmc.meta object. If \code{NULL}, names of parameters \eqn{phi} of annual model in Phase II will not be considered. If the meta object is provided, then \eqn{phi} will be included if \eqn{rho\_phase2} is considered and estimated in the MCMC process.}
\item{country.code}{Country code. If it is given, the country-specific parameter names contain the postfix \sQuote{_c\eqn{x}} where \eqn{x} is the \code{country.code}.}
\item{back.trans}{Logical indicating if back-transformable parameter names (i.e. \sQuote{Triangle_c1}, \dots, \sQuote{Triangle_c3}) should be returned.}
}


\value{
	\code{tfr.parameter.names} returns names of the country-independent Phase II parameters.\cr
	\code{tfr.parameter.names.cs} returns names of the country-specific Phase II parameters.\cr
	\code{tfr.parameter.names.extended} returns names of all country-independent Phase II parameters, including the transformed parameters. Parameters \sQuote{alpha}, \sQuote{delta}, \sQuote{alphat}, and\sQuote{deltat} are in their extended format with the postfix \sQuote{_1}, \sQuote{_2} and \sQuote{_3}. \cr
	\code{tfr.parameter.names.cs.extended} returns names of all country-specific Phase II parameters, including the transformed parameters. Parameters \sQuote{gamma} and\sQuote{gammat} are in their extended format with the postfix \sQuote{_1}, \sQuote{_2} and \sQuote{_3}.\cr
	\code{tfr3.parameter.names} returns names of the country-independent Phase III parameters.\cr
	\code{tfr3.parameter.names.cs} returns names of the country-specific Phase III parameters.\cr
}

\author{
Hana Sevcikova
}

\examples{
tfr.parameter.names()
tfr.parameter.names.extended()
tfr.parameter.names.cs()
tfr.parameter.names.cs.extended()
tfr3.parameter.names()
tfr3.parameter.names.cs()
}

\keyword{ manip }

