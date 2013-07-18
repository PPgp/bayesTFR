\name{tfr.predict}
\Rdversion{1.1}
\alias{tfr.predict}
\alias{bayesTFR.prediction}

\title{
Generating Posterior Trajectories of the Total Fertility Rate
}
\description{
Using the posterior parameter samples simulated by \code{\link{run.tfr.mcmc}} (and possibly \code{\link{run.tfr3.mcmc}}) the function  generates posterior trajectories for the total fertility rate for all countries of the world.
}
\usage{
tfr.predict(mcmc.set = NULL, end.year = 2100, 
    sim.dir = file.path(getwd(), "bayesTFR.output"), 
    replace.output = FALSE, start.year = NULL, 
    nr.traj = NULL, thin = NULL, burnin = 2000,
    use.diagnostics = FALSE, use.tfr3 = TRUE, burnin3 = 10000,
    mu = 2.1, rho = 0.8859, sigmaAR1 = 0.1016,
    use.correlation = FALSE, save.as.ascii = 1000, output.dir = NULL, 
    low.memory = TRUE, seed = NULL, verbose = TRUE, \dots)
}

\arguments{
  \item{mcmc.set}{Object of class \code{\link{bayesTFR.mcmc.set}} corresponding Phase II MCMCs. If it is \code{NULL}, the object is loaded from the directory given by \code{sim.dir}.}
  \item{end.year}{End year of the prediction.}
  \item{sim.dir}{Directory with the MCMC simulation results. It should equal to the \code{output.dir} argument in \code{\link{run.tfr.mcmc}}.}
  \item{replace.output}{Logical. If \code{TRUE}, existing predictions in \code{output.dir} will be replaced by results of this run.}
  \item{start.year}{Start year of the prediction. By default the prediction is started at the next time period after \code{present.year} set in the estimation step. If \code{start.year} is smaller than the default, projections for countries and time periods that have data available after \code{start.year} are set to those data.}
  \item{nr.traj}{Number of trajectories to be generated. If \code{NULL}, the argument \code{thin} is taken to determine the number of trajectories. If both are \code{NULL}, the number of trajectories corresponds to the size of the parameter sample.}
  \item{thin}{Thinning interval used for determining the number of trajectories. Only relevant, if \code{nr.traj} is \code{NULL}.}
  \item{burnin}{Number of iterations to be discarded from the beginning of the parameter traces.}
  \item{use.diagnostics}{Logical determining if an existing convergence diagnostics for phase II MCMCs should be used for choosing the values of \code{thin} and \code{burnin}. In such a case, arguments \code{nr.traj}, \code{thin} and \code{burnin} are ignored. The \sQuote{best} values are chosen from results of running the \code{\link{tfr.diagnose}} function. Only diagnostics can be used that suggest a convergence of the underlying MCMCs. If there are more than one such objects, the one is chosen whose recommendation for the number of trajectories is larger and closest to 2000.}
  \item{use.tfr3}{Logical determining if phase III should be predicted via MCMCs (simulated via \code{\link{run.tfr3.mcmc}}) or a classic AR(1) process. If \code{TRUE} but no phase III MCMCs were simulated, a warning is given and the prediction switches automatically to a classic AR(1) process.}
  \item{burnin3}{Burnin used for phase III MCMCs. Only relevant if \code{use.tfr3} is \code{TRUE}.}
  \item{save.as.ascii}{Either a number determining how many trajectories should be converted into an ASCII file, or \dQuote{all} in which case all trajectories are converted. It should be set to 0, if no conversion is desired.}
  \item{output.dir}{Directory into which the resulting prediction object and the trajectories are stored. If it is \code{NULL}, it is set to either \code{sim.dir}, or to \code{output.dir} of \code{mcmc.set$meta} if \code{mcmc.set} is given.}
  \item{low.memory}{Logical indicating if the prediction should run in a low-memory mode. If it is \code{FALSE}, the whole traces of all parameters, including the burnin, are loaded into memory. Otherwise, burnins are discarded and parameters are loaded as they are needed and are not kept in the memory.}
  \item{mu}{Long-term mean \eqn{\mu}{mu} in the AR(1) projection model. Only used if \code{use.tfr3} is \code{FALSE}.}
  \item{rho}{Autoregressive parameter \eqn{\rho}{rho} in AR(1) projection model. If it is \code{NULL} it is estimated from the data. Only used if \code{use.tfr3} is \code{FALSE}.}
  \item{sigmaAR1}{Standard deviation \eqn{s} of AR(1) distortion terms in short-term projections. If it is \code{NULL} it is estimated from the data. It can be a single value or a vector giving the standard deviations for single projections. If the vector is shorter than the number of projections simulated via the AR(1) process, the last value is repeated for the remaining projections. In case of a single value (default), the standard deviation is kept constant over all AR(1) projections. Only used if \code{use.tfr3} is \code{FALSE}.}
  \item{use.correlation}{Logical. If \code{TRUE} the model errors are sampled jointly for all countries (Fosdick and Raftery, 2012).}
  \item{seed}{Seed of the random number generator. If \code{NULL} no seed is set. It can be used to generate reproducible projections.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{\dots}{Further arguments passed to the underlying functions.}
}
\details{
The trajectories are generated using a distribution of country-specific decline curves (Alkema et al 2011) and either a classic AR(1) process or a country-specific AR(1) process (Raftery et al 2013).  Phase II parameter samples  simulated using \code{\link{run.tfr.mcmc}} are used from all chains, from which the given burnin was discarded. They are evenly thinned to match \code{nr.traj} or using the \code{thin} argument. Such thinned parameter traces, collapsed into one chain, if they do not already exist, are stored on disk into the sub-directory \file{\{thinned_mcmc_\emph{t}_\emph{b}} where \emph{t} is the value  of \code{thin} and \emph{b} the value of \code{burnin} (see \code{\link{create.thinned.tfr.mcmc}}). 
	
If Phase III is projected using a BHM (i.e. if \code{use.tfr3} is \code{TRUE}), parameter samples simulated via \code{\link{run.tfr3.mcmc}} are used from which burnin (given by \code{burnin3}) is discarded and the chains are evenly thinned in a way that the total size corresponds to the final size of the Phase II parameter samples. Countries for which there are no simulated country-specific Phase III parameters (e.g. because their TFR is still in Phase II or it is an aggregated region) use samples of the  \dQuote{world} AR(1) parameters.

The projection is run for all missing values before the present year, if any. Medians over the trajectories are used as  imputed values and the trajectories are discarded. The process then continues by projecting the future values where all generated trajectories are kept.

The resulting prediction object is saved into \file{\{output.dir\}/predictions}. Trajectories for all countries are saved into the same directory in a binary format, one file per country. At the end of the projection, if \code{save.as.ascii} is larger than 0, the function converts the given number of trajectories into a CSV file of a UN-specific format. They are selected by equal spacing (see function \code{\link{convert.tfr.trajectories}} for more details on the conversion). In addition, two summary files are created: one in a user-friendly format, the other using a UN-specific coding of the variants and time (see \code{\link{write.projection.summary}} for more details).
}

\value{
	Object of class \code{bayesTFR.prediction} which is a list containing components:
	\item{quantiles}{A \eqn{n \times q \times p}{n x q x p} array of quantile values computed on all trajectories. \eqn{n} is the number of countries,
		\eqn{q} is the number of quantile bounds and \eqn{p} is the number of projections.}
	\item{traj.mean.sd}{A \eqn{n \times 2 \times p}{n x 2 x p} array holding the mean of all trajectories in the first column and the standard deviation in the second column.
		\eqn{n} and \eqn{p} are the number of countries and number of projections, respectively.}
	\item{nr.traj}{Number of trajectories.}
	\item{trf_matrix_reconstructed}{Matrix containing imputed TFR values on spots where the original TFR matrix has missing values, i.e. between the last observed data point and the present year. }
	\item{output.directory}{Directory where trajectories corresponding to this prediction are stored.}
	\item{nr.projections}{Number of projections.}
	\item{burnin, thin, burnin3, thin3}{Burnin and thin used for this prediction for Phase II and Phase III, respectively.}
	\item{end.year}{The end year of this prediction.}
	\item{mu, rho, sigma_t, sigmaAR1}{Parameters of the AR(1) process. \code{sigma_t} is a vector of actual values of the standard deviation \eqn{s} used for each projection.}
	\item{na.index}{Index of trajectories for which at least one country got \code{NA} values.}
	\item{mcmc.set}{Object of class \code{\link{bayesTFR.mcmc.set}} used for this prediction, i.e. the burned, thinned, and collapsed MCMC chain.}
	}
	
\references{
L. Alkema, A. E. Raftery, P. Gerland, S. J. Clark, F. Pelletier, Buettner, T., Heilig, G.K. (2011). Probabilistic Projections of the Total Fertility Rate for All Countries. Demography, Vol. 48, 815-839. (A Working Paper version at \url{http://www.csss.washington.edu/Papers} [nr. 97])

Raftery, A.E., Alkema, L. and Gerland, P. (2013). Bayesian Population Projections for the United Nations. 
Statistical Science, in press.

Fosdick, B., Raftery, A.E. (2012). Regional Probabilistic Fertility Forecasting by Modeling Between-Country Correlations. Working paper, CSSS, University of Washington. \url{http://www.csss.washington.edu/Papers} [nr. 126]

}
\author{
Hana Sevcikova, Leontine Alkema, Bailey Fosdick
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{run.tfr.mcmc}}, \code{\link{run.tfr3.mcmc}}, \code{\link{create.thinned.tfr.mcmc}}, \code{\link{convert.tfr.trajectories}}, \code{\link{write.projection.summary}},    
\code{\link{get.tfr.prediction}}, \code{\link{summary.bayesTFR.prediction}}
}
\examples{
\dontrun{
sim.dir <- tempfile()
m <- run.tfr.mcmc(nr.chains=1, iter=10, output.dir=sim.dir, verbose=TRUE)
m3 <- run.tfr3.mcmc(sim.dir=sim.dir, nr.chains=2, iter=40, thin=1, verbose=TRUE)
pred <- tfr.predict(m, burnin=0, burnin3=10, verbose=TRUE)
summary(pred, country="Iceland")
unlink(sim.dir, recursive=TRUE)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{multivariate}
