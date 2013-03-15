\name{write.projection.summary}
\Rdversion{1.1}
\alias{write.projection.summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Writing Projection Summary Files
}
\description{
The function creates two files containing projection summaries, such as the median, the lower and upper bound of the 80 and 90\% probability intervals, respectively, the +/- 0.5 child variant and the constant variant. One file is in a user-friendly format, whereas the other is in a UN-specific format with internal coding of the time and the variants. In addition, a file containing  some of the model parameters is created.
}
\usage{
write.projection.summary(dir = file.path(getwd(), "bayesTFR.output"), 
    output.dir = NULL, revision = NULL, adjusted = FALSE)
}
\arguments{
  \item{dir}{Directory containing the prediction object. It should correspond to the \code{output.dir} argument of the \code{\link{tfr.predict}} function.}
  \item{output.dir}{Directory in which the resulting file will be stored. If \code{NULL} the same directory is used as for the prediction.}
  \item{revision}{UN WPP revision number. If \code{NULL} it is determined from the corresponding WPP year: WPP 2008 corresponds to revision 13, every subsequent WPP increases the revision number by one. Used as a constant in the second file only.}
  \item{adjusted}{Logical. By default the function writes summary using the original BHM projections. If the projection medians are adjusted (using e.g. \code{\link{tfr.median.set}}), setting this argument to \code{TRUE} causes using the adjusted projections.}
}
\details{
The first file that the function creates is called \file{projection_summary_user_friendly.csv} (or \file{projection_summary_user_friendly_adjusted.csv} if \code{adjusted=TRUE}), it is a comma-separated table with the following columns: 
\itemize{\item{\dQuote{country_name}: }{country name}
		 \item{\dQuote{country_code}: }{country code}
		 \item{\dQuote{variant}: }{name of the variant, such as \dQuote{median}, \dQuote{lower 80}, \dQuote{upper 80}, \dQuote{lower 95}, \dQuote{upper 95}, \dQuote{-0.5child}, \dQuote{+0.5child}, \dQuote{constant}}
		\item{period1: }{e.g. \dQuote{2005-2010}: TFR for the first time period}
		\item{period2: }{e.g. \dQuote{2010-2015}: TFR for the second time period}
		\item{\dots }{further columns with TFR projections}
	}

The second file, called \file{projection_summary.csv} (or \file{projection_summary_adjusted.csv} if \code{adjusted=TRUE}), also comma-separated table, contains the same information as above in a UN-specific format:
\itemize{\item{\dQuote{RevID}: }{revision number, passed to the function as an argument}
		  \item{\dQuote{VarID}: }{variant identifier, extracted from the \code{\link{UN_variants}} dataset}
		  \item{\dQuote{LocID}: }{country code}
		  \item{\dQuote{TimeID}: }{time identifier, extracted from the \code{\link{UN_time}} dataset}
		  \item{\dQuote{TFR}: }{the total fertility rate for this variant, location and time period}
}

The third comma-separated file, called \file{projection_summary_parameters.csv} contains the following columns:
\itemize{\item{\dQuote{country_name}: }{country name}
		 \item{\dQuote{country_code}: }{country code}
		 \item{\dQuote{TF_time_start_decline}: }{start period of TFR decline}
		 \item{\dQuote{TF_max}: }{TFR at the onset of the fertitility transition (median of the \eqn{U_c} parameter)}
		 \item{\dQuote{TF_max_decrement}: }{maximum decrement of TFR decline (median of the \eqn{d_c} parameter)}
		 \item{\dQuote{TF_end_level}: }{median of the end level of the fertility transition (\eqn{\Delta_{c4}}{Triangle_c4} parameter)}
		 \item{\dQuote{TF_end_level_low}: }{2.5 percentile of the \eqn{\Delta_{c4}}{Triangle_c4} distribution}
		 \item{\dQuote{TF_end_level_high}: }{97.5 percentile of the \eqn{\Delta_{c4}}{Triangle_c4} distribution}
		 \item{\dQuote{TF_time_end_decline}: }{period of the end decline, measured using the prediction median}
	}
Note that this file is not created if \code{adjusted=TRUE}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Hana Sevcikova
}
\note{
This function is automatically called from the \code{\link{tfr.predict}} function, therefore in standard cases it will not be needed to call it directly.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{convert.tfr.trajectories}}, \code{\link{tfr.predict}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
