\name{tfr.trajectories.plot}
\Rdversion{1.1}
\alias{tfr.trajectories.plot}
\alias{tfr.trajectories.plot.all}
\alias{tfr.trajectories.table}

\title{
Output of Posterior Distribution of TFR Trajectories}
\description{
The functions plot/tabulate the posterior distribution of TFR trajectories for a given country, or for all countries, including their median and given probability intervals. 
}
\usage{
tfr.trajectories.plot(tfr.pred, country, pi = c(80, 95), 
    half.child.variant = TRUE, nr.traj = NULL, 
    adjusted.only = TRUE, typical.trajectory = FALSE, 
    mark.estimation.points = FALSE,
    xlim = NULL, ylim = NULL, type = 'b', xlab = 'Year', ylab = 'TFR', 
    main = NULL, lwd = c(2, 2, 2, 2, 2, 1), 
    col=c('black', 'green', 'red', 'red', 'blue', '#00000020'),
    show.legend = TRUE, add = FALSE, uncertainty = FALSE, 
    col_unc = "purple", \dots)
	
tfr.trajectories.plot.all(tfr.pred, 
    output.dir = file.path(getwd(), 'TFRtrajectories'),
    output.type = "png", main = NULL, verbose = FALSE, \dots)
	
tfr.trajectories.table(tfr.pred, country, pi = c(80, 95), 
    half.child.variant = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tfr.pred}{Object of class \code{\link{bayesTFR.prediction}}.}
  \item{country}{Name or numerical code of a country.}
  \item{pi}{Probability interval. It can be a single number or an array.}
  \item{half.child.variant}{If \code{TRUE} the United Nations variant of \dQuote{\code{+/-0.5} child} (relative to the median) is shown.}
  \item{nr.traj}{Number of trajectories to be plotted. If \code{NULL}, all trajectories are plotted, otherwise they are thinned evenly.}
  \item{adjusted.only}{Logical. By default, if the projection median is adjusted using e.g. \code{\link{tfr.median.set}}, the function plots the adjusted median. If \code{adjusted.only=FALSE} the original (non-adjusted) median is plotted as well.}
  \item{typical.trajectory}{Logical. If \code{TRUE} one trajectory is shown that has the smallest distance to the median.}
  \item{mark.estimation.points}{Logical. If \code{TRUE}, points that were not used in the estimation (phase I) are shown in lighter color than points in phase II and III.}
  \item{xlim, ylim, type, xlab, ylab}{Graphical parameters passed to the \code{plot} function.}
    \item{lwd, col}{Vector of six elements giving the line width and color for: 1. observed data, 2. imputed missing data, 3. median, 4. quantiles, 5. half-child variant, 6. trajectories.}
  \item{show.legend}{Logical controlling whether the legend should be drawn.}
  \item{add}{Logical controlling whether the trajectories should be plotted into a new graphic device (\code{FALSE}) or into an existing device (\code{TRUE}). One can use this argument to plot trajectories from multiple countries into one graphics.}
  \item{\dots}{Additional graphical parameters. In addition, for \code{tfr.trajectories.plot.all}, \code{\dots} contains any of the arguments of \code{tfr.trajectories.plot}.}
  \item{output.dir}{Directory into which resulting graphs are stored.}
  \item{output.type}{Type of the resulting files. It can be \dQuote{png}, \dQuote{pdf}, \dQuote{jpeg}, \dQuote{bmp}, \dQuote{tiff}, or \dQuote{postscript}.}
  \item{main}{Main title for the plot(s). In \code{tfr.trajectories.plot.all} any occurence of the string \dQuote{XXX} is replaced by the country name.}
  \item{verbose}{Logical switching log messages on and off.}
  \item{uncertainty}{Logical: \code{TRUE} means uncertainty of past TFR should be plotted with the same level of uncertainty interval.}
  \item{col_unc}{Color of past TFR estimation uncertainty plot.}
}
\details{
\code{tfr.trajectories.plot} plots posterior distribution of TFR trajectories for a given country. \code{tfr.trajectories.table} gives the same output as a table. \code{tfr.trajectories.plot.all} creates a set of graphs (one per country) that are stored in \code{output.dir}. 

The median and given probability intervals are computed using all available trajectories. Thus, \code{nr.traj} does not influence those values - it is used only to control the number of trajectories plotted.
}

\author{
Hana Sevcikova, Leontine Alkema
}

\seealso{
\code{\link{bayesTFR.prediction}}
}
\examples{
\dontrun{
sim.dir <- file.path(find.package("bayesTFR"), "ex-data", "bayesTFR.output")
pred <- get.tfr.prediction(sim.dir)
tfr.trajectories.plot(pred, country="Burkina Faso", pi=c(80, 95))
tfr.trajectories.table(pred, country="Burkina Faso", pi=c(80, 95))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
