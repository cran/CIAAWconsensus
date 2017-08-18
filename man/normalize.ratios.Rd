\name{normalize.ratios}
\alias{normalize.ratios}
\title{Normalize isotope amount ratios to a common reference isotope}
\description{
This function converts the isotope amount ratios of an element from various studies to a single common reference isotope so that all isotope ratios can be directly compared to one another.
The conversion involves a direct application of the law of propagation of uncertainty and this function discards the possible covariances between the isotope ratios.
}
\usage{normalize.ratios(dat, element, ref.isotope, expand = FALSE)}
\arguments{
\item{dat}{
A data frame of results from each study where each study reports one or more isotope ratios (outcomes).
The data frame must include the following named columns: Study, Year, Author, Outcome, Value, Unc, k_extra (see Details).
}
\item{element}{Lowercase english name of the element, e.g., "antimony"}
\item{ref.isotope}{Desired reference isotope, e.g., "121Sb"}
\item{expand}{(Logical) Specification of whether or not to expand the isotope ratio uncertainties using the values of \code{dat$k_extra}}
}
\details{
The isotope ratio vector \strong{\eqn{R}} is transformed to the reference isotope \eqn{R_{ref}} by dividing each element of the set to the chosen reference isotope.
The covariances of the transformed isotope ratios are obtained using the Law of Propagation of Uncertainty. This function assumes all isotope ratios reported by a given study as uncorrelated.
While this is not strictly true in practice, such assumption is made largely because of the lack of reported correlations in the literature.

The format of \code{dat} data frame for a simple dataset is as follows:
\tabular{lllllll}{
    \strong{Study} \tab \strong{Year} \tab \strong{Author} \tab \strong{Outcome} \tab \strong{Value} \tab \strong{Unc} \tab \strong{k_extra} \cr
    1 \tab 1954 \tab Howard  \tab 191Ir/193Ir \tab 0.5949  \tab 0.0025  \tab 9 \cr
    2 \tab 1991 \tab Creaser \tab 191Ir/193Ir \tab 0.5948  \tab 0.0001  \tab 9 \cr
    3 \tab 1992 \tab Chang   \tab 191Ir/193Ir \tab 0.59399 \tab 0.00103 \tab 6 \cr
	4 \tab 1993 \tab Walczyk \tab 191Ir/193Ir \tab 0.59418 \tab 0.00037 \tab 9
  }
}
\value{
\item{R}{A list of the normalized isotope amount ratios}
\item{u.R}{A list of standard uncertainties for \code{R}}
\item{cov.R}{A list of covariance matrices for \code{R} (forced to zero covariances)}
}
\references{
J. Meija and Z. Mester (2008) \href{http://dx.doi.org/10.1088/0026-1394/45/1/008}{Uncertainty propagation of atomic weight measurement results}. \emph{Metrologia}, 45, 53-62

J. Meija and A. Possolo (2017) \href{https://doi.org/10.1088/1681-7575/aa634d}{Data reduction framework for standard atomic weights and isotopic compositions of the elements}. \emph{Metrologia}, 54, 229-238
}
\author{Juris Meija <juris.meija@nrc-cnrc.gc.ca> and Antonio Possolo}

\examples{
## Normalize all platinum isotope data to platinum-195
normalize.ratios(platinum.data, "platinum", "195Pt")
}