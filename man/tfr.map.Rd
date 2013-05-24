\name{tfr.map}
\Rdversion{1.1}
\alias{tfr.map}
\alias{tfr.map.all}
\alias{get.tfr.map.parameters}
\alias{tfr.map.gvis}
\title{
TFR World Map 
}
\description{
Generates a world map of the total fertility rate for given projection period and quantile. In addition, country specific Phase II MCMC parameters can be projected into the world map. 
}
\usage{
tfr.map(pred, quantile = 0.5, 
    year = NULL, par.name = NULL, adjusted = FALSE,
    projection.index = 1, device = "dev.new", main = NULL, 
    resolution=c("coarse","low","less islands","li","high"),
    device.args = NULL, data.args = NULL, \dots)
	
tfr.map.all(pred, output.dir, output.type = "png", 
    tfr.range = NULL, nr.cats = 50, same.scale = TRUE, 
    quantile = 0.5, file.prefix='TFRwrldmap_', \dots)
			
get.tfr.map.parameters(pred, tfr.range = NULL, 
    nr.cats = 50, same.scale = TRUE, quantile = 0.5, \dots)
    
tfr.map.gvis(pred, year = NULL, quantile = 0.5, pi = 80, 
    par.name = NULL, adjusted = FALSE, \dots)
}

\arguments{
  \item{pred}{Object of class \code{\link{bayesTFR.prediction}}.}
  \item{quantile}{Quantile for which the map should be generated. It must be equal to one of the values in \code{dimnames(pred$quantiles[[2]])}, i.e. 0, 0.025, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.975, 1. Value 0.5 corresponds to the median.}
  \item{year}{Year to be plotted. It can be a year within a projection period or a year within an estimation period. In the latter case, the observed data are plotted. If not given, \code{projection.index} determines the projection year.}
  \item{par.name}{Name of a country-specific parameter to be plotted. If \code{NULL}, the TFR is plotted. Allowed values are any of those returned by \code{tfr.parameter.names.cs.extended()} and \sQuote{lambda} (see Details).}
  \item{adjusted}{Logical indicating if the measure to be plotted is based on adjusted TFRs.}
  \item{projection.index}{Index of the projection to be displayed. It is only relevant if \code{year} is \code{NULL}. \code{projection.index=1} means the present year, \code{projection.index=2} means the first projection period after present year, etc..
}
  \item{device}{Device for displaying the map. It is passed to the \code{\link[rworldmap]{mapDevice}} function of the \pkg{rworldmap} package. If it is equal to \sQuote{dev.cur}, the current device is used. Otherwise, it can be \sQuote{dev.new}, \sQuote{png}, \sQuote{pdf} etc.
}
  \item{main}{
Title for the map. If it is \code{NULL}, a default title is constructed from the projection year and quantile.
}
\item{resolution}{Map resolution as implemented in \code{\link[rworldmap]{getMap}}. High resolution requires the \pkg{rworldxtra} package.}
  \item{device.args}{List of 
additional arguments to be passed to the \code{\link[rworldmap]{mapDevice}} function of the \pkg{rworldmap} package.
}
 \item{data.args}{List of additional arguments to be passed to the underlying data retrieving function.} 
 \item{output.dir}{Directory into which resulting maps are stored.}
 \item{output.type}{Type of the resulting files. It can be \dQuote{png}, \dQuote{pdf}, \dQuote{jpeg}, \dQuote{bmp}, \dQuote{tiff}, or \dQuote{postscript}.}
 \item{tfr.range}{Range of the total fertility rate to be displayed. It is of the form \code{c(}\var{tfr.min}, \var{tfr.max}\code{)}.
 By default, the whole range is considered. Note that countries with values outside of the given range will appear white.}
  \item{nr.cats}{Number of color categories.}
  \item{same.scale}{Logical controlling if maps for all projection years of this prediction object should be on the same color scale.}
    \item{file.prefix}{Prefix for file names.}
  \item{\dots}{
Arguments passed to the \code{\link[rworldmap]{mapCountryData}} function of the \pkg{rworldmap} package. In case of \code{tfr.map.gvis} these are passed to the underlying data retrieving function (the same as \code{data.args}).
}
\item{pi}{Probability interval to be shown when a country is selected in an interactive map. The corresponding quantiles must be available (see argument \code{quantile} above).}
}
\details{
\code{tfr.map} creates a single map for a given projection period and quantile using the \pkg{rworldmap} package. \code{tfr.map.all} generates a sequence of such maps, namely one for each projection period. If the package \pkg{fields} is installed, a color bar legend at the botom of the map is created.

Function \code{get.tfr.map.parameters} can be used in combination with \code{tfr.map}. (Note that \code{get.tfr.map.parameters} is called from inside of \code{tfr.map.all}.) It sets breakpoints for the color scheme using quantiles of a fitted gamma distribution.

Function \code{tfr.map.gvis} creates an interactive map using the \pkg{googleVis} package and opens it in an internet browser. It also generates a table of TFRs that can be sorted by columns interactively in the browser. 

By default, both \code{tfr.map} and \code{tfr.map.gvis} produce maps of TFRs. Alternatively, the functions can be used to plot country-specific Phase II MCMC parameters into a world map. They are given by the argument \code{par.name}. In addition to the MCMC parameters, if \code{par.name='lambda'}, the period of the end of TFR decline (i.e. start of Phase III) is computed for each country and projected into the map. In such a case, we recommend to adjust the color scale in \code{tfr.map} e.g. using the arguments \code{catMethod='pretty'} and \code{numCats=20} (see \code{\link[rworldmap]{mapCountryData}}). 
}
\value{
\code{get.tfr.map.parameters} returns a list with elements:
  \item{pred}{The object of class \code{\link{bayesTFR.prediction}} used in the function.}
  \item{quantile}{Value of the argument \code{quantile}.}
  \item{catMethod}{If the argument \code{same.scale} is \code{TRUE}, this element contains breakpoints for categorization. It is generated from a fitted gamma distribution. Otherwise, it is \code{NULL}.}
  \item{numCats}{Number of categories.}
  \item{coulourPalette}{Subset of the rainbow palette, starting from dark blue and ending at red.}
  \item{\dots}{Additional arguments passed to the function.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova, Patrick Gerland, Adrian Raftery
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
sim.dir <- file.path(find.package("bayesTFR"), "ex-data", "bayesTFR.output")
pred <- get.tfr.prediction(sim.dir=sim.dir)
# Uses heat colors and seven categories by default
tfr.map(pred)
# Uses more colors with more suitable categorization
params <- get.tfr.map.parameters(pred)
do.call("tfr.map", params)
# Another projection year on the same scale
do.call("tfr.map", c(list(year=2043), params))

# Using Google Vizualization tool
tfr.map.gvis(pred)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
