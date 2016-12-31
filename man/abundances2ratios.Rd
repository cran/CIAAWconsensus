\name{abundances2ratios}
\alias{abundances2ratios}
\title{Isotope ratios of a chemical element from isotopic abundances}
\description{
This function calculates the isotope ratios of a chemical element from the given isotopic abundances and their uncertainties.
The uncertainty evaluation is done using the propagation of uncertainty and the missing correlations between the isotopic abundances are reconstructed using Monte Carlo methods.
}
\usage{
abundances2ratios(x, ux, ref=1, iterations=1e4)
}

\arguments{
  \item{x}{A vector of isotopic abundances of an element}
  \item{ux}{Standard uncertainties of \code{x}}
  \item{ref}{Index to specify the desired reference isotope for isotope amount ratios}
  \item{iterations}{Number of iterations for isotopic abundance correlation mapping}
}
\details{
Situations are often encountered where isotopic abundances are reported but not the isotope ratios.
In such cases we reconstruct the isotope ratios that are consistent with the abundances and their uncertainties. 
Given only the abundances and their uncertainties, for elements with four or more isotopes one cannot unambiguously infer the uncertainties of the ratios due to the unknown correlations between isotopic abundances.
This missing information can be reconstructed by mapping all possible correlations between isotopic abundances.
}
\value{
\item{R}{Isotope ratio vector, excluding the trivial ratio}
\item{R.u}{Standard uncertainties of the isotope ratios}
\item{R.cov}{Isotope ratio covariance matrix}
\item{N}{Number of sucessful Monte Carlo iterations}
}
\references{
J. Meija and Z. Mester (2008) \href{http://dx.doi.org/10.1088/0026-1394/45/4/012}{Atomic weight uncertainty calculation from isotopic composition of the elements}. \emph{Metrologia}, 45, 459

J. Meija, A. Possolo (2017) Data reduction framework for standard atomic weights and isotopic compositions of the elements. \emph{Metrologia}, submitted

JCGM 101:2008 Evaluation of measurement data - Supplement 1 to the "Guide to the expression of uncertainty in measurement" - Propagation of distributions using a Monte Carlo method
}
\author{Juris Meija <juris.meija@nrc-cnrc.gc.ca> and Antonio Possolo}
\examples{
## Isotope ratios of zinc from the isotopic abundances
x = c(0.48630, 0.27900, 0.04100, 0.18750, 0.00620)
ux = c(0.00091, 0.00076, 0.00031, 0.00135, 0.00010)
abundances2ratios(x,ux,ref=2)

## The corresponding atomic weight can be obtained using at.weight(z$R,z$R.cov,"zinc","66Zn")
}