\name{tfr.predict.subnat}
\alias{tfr.predict.subnat}

\title{
Generating Posterior Trajectories of Subnational TFR
}
\description{
Using the posterior parameter samples simulated by \code{\link{run.tfr.mcmc.subnat}} (and possibly \code{\link{run.tfr3.mcmc.subnat}}) the function  generates posterior trajectories of the total fertility rate for sub-regions of given countries.
}
\usage{
tfr.predict.subnat(mcmc.set.list = NULL, sim.dir = file.path(getwd(), "bayesTFR.output"), 
    countries = NULL, nr.traj = NULL, thin = NULL, burnin = 2000, use.diagnostics = FALSE, 
    use.tfr3 = TRUE, burnin3 = 2000, mu = 2.1, rho = 0.8859, use.phase3.from = NULL, 
    use.correlation = FALSE, cor.method = c("bayes", "median", "mean"), 
    cor.start.year = NULL, output.dir = NULL, \dots, verbose = TRUE)
}

\arguments{
  \item{mcmc.set.list}{List of \code{\link{bayesTFR.mcmc.set}} objects, as returned by \code{\link{run.tfr.mcmc.subnat}}. They correspond to sub-national Phase II MCMCs, one object per country. Names of the list items must correspond to country codes. If this argument is \code{NULL}, the objects are loaded from the directory given by \code{sim.dir}.}
  \item{sim.dir}{Simulation directory. It should equal to the \code{output.dir} argument in \code{\link{run.tfr.mcmc.subnat}}.}
  \item{countries}{Numerical country codes. If it is \code{NULL} projections are generated for all countries for which MCMCs exist in the \sQuote{subnat} sub-directory.}
  \item{nr.traj}{Number of trajectories to be generated. If \code{NULL}, the argument \code{thin} is taken to determine the number of trajectories. If both are \code{NULL}, the number of trajectories corresponds to the size of the parameter sample.}
  \item{thin}{Thinning interval used for determining the number of trajectories. Only relevant, if \code{nr.traj} is \code{NULL}.}
  \item{burnin}{Number of iterations to be discarded from the beginning of the parameter traces.}
  \item{use.diagnostics}{Logical determining if an existing convergence diagnostics for phase II MCMCs should be used for choosing the values of \code{thin} and \code{burnin}. It has the same meaning as in \code{\link{tfr.predict}}.}
  \item{use.tfr3}{Logical determining if phase III should be predicted via MCMCs (simulated via \code{\link{run.tfr3.mcmc.subnat}}) or a classic AR(1) process. If \code{TRUE} but no phase III MCMCs were simulated, a warning is given and the prediction switches automatically to a classic AR(1) process.}
  \item{burnin3}{
%%     ~~Describe \code{burnin3} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{rho}{
%%     ~~Describe \code{rho} here~~
}
  \item{use.phase3.from}{
%%     ~~Describe \code{use.phase3.from} here~~
}
  \item{use.correlation}{
%%     ~~Describe \code{use.correlation} here~~
}
  \item{cor.method}{
%%     ~~Describe \code{cor.method} here~~
}
  \item{cor.start.year}{
%%     ~~Describe \code{cor.start.year} here~~
}
  \item{output.dir}{
%%     ~~Describe \code{output.dir} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (mcmc.set.list = NULL, sim.dir = file.path(getwd(), 
    "bayesTFR.output"), countries = NULL, nr.traj = NULL, thin = NULL, 
    burnin = 2000, use.diagnostics = FALSE, use.tfr3 = TRUE, 
    burnin3 = 10000, mu = 2.1, rho = 0.8859, use.phase3.from = NULL, 
    use.correlation = FALSE, cor.method = c("bayes", "median", 
        "mean"), cor.start.year = NULL, output.dir = NULL, ..., 
    verbose = TRUE) 
{
    if (is.null(mcmc.set.list)) {
        cdirs <- list.files(file.path(sim.dir, "subnat"), pattern = "^c[0-9]+", 
            full.names = TRUE)
        for (i in 1:length(cdirs)) {
            cdir.short <- strsplit(cdirs[i], split = .Platform$file.sep)
            cdir.short <- cdir.short[[1]][length(cdir.short[[1]])]
            country.code <- substring(cdir.short, 2)
            if (is.null(countries) || is.element(country.code, 
                countries)) 
                mcmc.set.list[[country.code]] <- get.tfr.mcmc(cdirs[i])
        }
    }
    if (is.null(nr.traj) && is.null(thin) && !use.diagnostics) 
        stop("Either nr.traj or thin must be given or use.diagnostics must be TRUE.")
    result <- cor.mat <- NULL
    outdir <- NULL
    for (country in names(mcmc.set.list)) {
        if (use.correlation) {
            if (verbose) 
                cat("\nCompute correlation for country", country, 
                  "...")
            snr.traj <- nr.traj
            sburnin <- burnin
            if (use.diagnostics) {
                nrtraj.burnin <- get.burnin.nrtraj.from.diagnostics(mcmc.set.list[[country]])
                sburnin <- nrtraj.burnin[2]
            }
            has.phase3 <- FALSE
            ar1pars <- NULL
            use.external.phase3 <- FALSE
            if (use.tfr3) {
                has.phase3 <- has.tfr3.mcmc(mcmc.set.list[[country]]$meta$output.dir)
                if (!has.phase3 && !is.null(use.phase3.from) && 
                  has.tfr3.mcmc(use.phase3.from)) {
                  has.phase3 <- TRUE
                  use.external.phase3 <- TRUE
                }
            }
            if (!has.phase3) {
                if (is.null(rho) || is.na(rho)) 
                  ar1pars <- get.ar1.parameters(mu = mu, mcmc.set.list[[country]]$meta)
            }
            cor.mat <- tfr.correlation.subnat(mcmc.set.list[[country]], 
                method = cor.method, burnin = sburnin, burnin3 = burnin3, 
                use.external.phase3 = use.external.phase3, use.phase3.from = use.phase3.from, 
                ar1pars = ar1pars, cor.start.year = cor.start.year, 
                verbose = verbose)
            if (verbose) 
                cat(" done.\n")
        }
        if (!is.null(output.dir)) {
            outdir <- file.path(output.dir, "subnat", paste0("c", 
                country))
        }
        result[[country]] <- tfr.predict(mcmc.set = mcmc.set.list[[country]], 
            nr.traj = nr.traj, thin = thin, burnin = burnin, 
            use.tfr3 = use.tfr3, burnin3 = burnin3, use.diagnostics = use.diagnostics, 
            mu = mu, rho = rho, use.phase3.from = use.phase3.from, 
            use.correlation = use.correlation, correlation.matrices = cor.mat, 
            output.dir = outdir, verbose = verbose, ...)
    }
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
