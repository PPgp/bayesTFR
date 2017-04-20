\name{get.regtfr.prediction}
\alias{get.regtfr.prediction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Accessing Subnational Prediction Objects
}
\description{
Retrieve subnational (regional) prediction results produced by \code{\link{tfr.predict.subnat}}, either for one country or for all available countries.
}

\usage{
get.regtfr.prediction(sim.dir, country = NULL)
}

\arguments{
  \item{sim.dir}{Simulation directory of the subnational predictions. It corresponds to the argument \code{output.dir} in \code{\link{tfr.predict.subnat}}.
}
  \item{country}{Numerical country code. If it is not given, all available subnational predictions are retrieved.}
}
\details{
Predictions for country \eqn{x} are assumed to be stored in \dQuote{\code{sim.dir}/subnat/c\eqn{x}}.
}

\value{
If argument \code{country} is given, the function returns an object of class \code{\link{bayesTFR.prediction}}. If it is \code{NULL}, it returns a list of such objects. Names of the list items are the country codes.
}




\seealso{
\code{\link{tfr.predict.subnat}}
}
\examples{
# Subnational example data
my.subtfr.file <- file.path(find.package("bayesTFR"), 'extdata', 'subnational_tfr_template.txt')
subtfr <- read.delim(my.subtfr.file, check.names=FALSE, stringsAsFactor=FALSE)
countries <- unique(subtfr[, c("country_code", "country")])

# Directory with national projections (contains 30 trajectories for each country)
nat.dir <- file.path(find.package("bayesTFR"), "ex-data", "bayesTFR.output")

# Subnational projections for all three countries ()
subnat.dir <- tempfile()
tfr.predict.subnat(countries$country_code, my.tfr.file=my.subtfr.file,
    sim.dir=nat.dir, output.dir=subnat.dir, start.year=2013)
    
# Retrieve results for all countries
preds <- get.regtfr.prediction(subnat.dir)
names(preds)

# View tables of subregions for each country
for(i in 1:nrow(countries)) {
  cat("\n\n", countries$country[i], "\n")
  print(get.countries.table(preds[[as.character(countries$country_code[i])]]))
}
# Quantiles for individual subregions
tfr.trajectories.table(preds[["218"]], "Bolivar")

# Retrieve results for one country
pred <- get.regtfr.prediction(subnat.dir, 218)
tfr.trajectories.plot(pred, "Loja")

# cleanup
unlink(subnat.dir)

# See more examples in ?tfr.predict.subnat
}

\keyword{ manip }% use one of  RShowDoc("KEYWORDS")

