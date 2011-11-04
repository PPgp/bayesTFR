\name{UN_variants}
\Rdversion{1.1}
\alias{UN_variants}
\docType{data}
\title{
Dataset with UN-specific Coding of Variants
}
\description{
Dataset used by the UN for coding variants. It also includes variants for the lower and upper bounds of the 80 and 95\% probability intervals, respectively, resulting from the Bayesian hierarchical model. The dataset is stored in a TAB-separated ASCII file called \dQuote{UN_variants.txt}.
}
\usage{data(UN_variants)}
\format{
  A data frame with 23 observations on the following 5 variables.
  \describe{
    \item{\code{RevID}}{Revision identifier.}
    \item{\code{VarID}}{Identifier of the variant.}
    \item{\code{Vshort}}{Short name of the variant.}
    \item{\code{VName}}{Full name of the variant.}
    \item{\code{VariantDomain}}{Domain of the variant.}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
Data provided by the United Nations Population Division
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(UN_variants)
str(UN_variants)
}
\keyword{datasets}
