\name{bayesTFR.mcmc.meta}
\Rdversion{1.1}
\alias{bayesTFR.mcmc.meta}

\title{
MCMC Simulation Meta Object
}
\description{
Simulation meta object \code{bayesTFR.mcmc.meta} used by all chains of the same MCMC simulation. 
It contains information that is common to all chains. It is part of a \code{\link{bayesTFR.mcmc.set}} object.
}


\value{
	A \code{bayesTFR.mcmc.meta} object contains components \code{U.width, U.c.low, U.up,
			a.width, a.low, a.up, 
			b.width, b.low, b.up,
			const.width, const.low, const.up,
			S.width, S.low, S.up,
			sigma0.width, sigma0.low, sigma0.up,
			sigma0.min, 
			mean.eps.tau0, sd.eps.tau0, nu.tau0,
			d.low, d.up, 
			Triangle4.0, delta4.0, nu4,
            Triangle_c4.low, Triangle_c4.up,
            Triangle_c4.trans.width,
			chi0, nu.psi0, psi0, d.trans.width, 
			dl.p1, dl.p2, 
			nu.delta0, delta0, alpha0.p, start.year, present.year, wpp.year, buffer.size, auto.conf}. Their meaning and value are the same as the input arguments 
			to the \code{\link{run.tfr.mcmc}} function. Furthermore, it contains components:
	\item{nr.chains}{Number of MCMC chains.}
	\item{tfr_matrix_all}{A \eqn{q \times n}{q x n} matrix with the United Nations TFR estimates. \eqn{q} is number of years (see \code{T_end} below), \eqn{n} is number of countries (see \code{nr_countries} below). The first \eqn{n_e} columns correspond to countries included in the MCMC estimation (see \code{nr_countries_estimation} below), where \eqn{n_e <= n}.}
	\item{tfr_matrix}{Like \code{tfr_matrix_all}, but it has \code{NA} values before and after country's fertility transition.}
	\item{tfr_matrix_observed}{Like \code{tfr_matrix}, but it has \code{NA} values for years where no historical data is available.}
	\item{nr_countries}{Number of countries included in the tfr matrices.}
	\item{nr_countries_estimation}{Number of countries included in the MCMC estimation. It must be smaller or equal to \code{nr_countries}.}
	\item{tau_c}{Estimated start year of the fertility decline for each country (as a row index within the tfr matrices). -1 means that the decline started before \code{start.year}.}
	\item{id_Tistau}{Index of countries for which \code{present.year} is equal to \code{tau_c}.}
	\item{id_DL}{Index of countries for which the projection is made using the double logistic function, i.e. high fertility countries.}
	\item{id_early}{Index of countries with early decline, i.e. countries for which \code{tau_c=-1}.}
	\item{id_notearly}{Index of countries with not early decline.}
	\item{lambda_c}{Start period of the recovery phase for each country (as an index of the \code{tfr_matrix}).}
	\item{start_c}{Maximum of \code{tau_c} and 1 for each country. Thus, it is a row index of the \code{tfr_matrix} where the fertility decline starts.}
	\item{proposal_cov_gammas_cii}{Proposal covariance matrices of \eqn{\gamma_{ci}}{gamma_ci} for each country.}
	\item{T_end}{Number of years for which United Nations historical data are available (i.e. number of rows of \code{tfr_matrix}).}
	\item{T_end_c}{Like \code{T_end} but country specific.}
	\item{regions}{List of arrays of length \code{nr_countries}. These are: \cr
		\code{name} - Region name for each country.\cr
		\code{code} - Region code for each country. \cr
		\code{area_name} - Area name for each country. \cr
		\code{area_code} - Area code for each country. \cr
		\code{country_name} - Array of country names. \cr
		\code{country_code} - Array of country codes. \cr
		Any country indices in the \code{bayesTFR.mcmc.meta} object are derived from this component.}
	\item{output.dir}{Directory for storing simulation output.}
}

\details{
The object is in standard cases not to be manipulated by itself, but rather as part of a \code{\link{bayesTFR.mcmc.set}} object.
}

%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova, Leontine Alkema
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{run.tfr.mcmc}}, \code{\link{get.tfr.mcmc}}
}
\examples{
sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
m <- get.tfr.mcmc(sim.dir)
summary(m, meta.only = TRUE)
names(m$meta)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classes }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
