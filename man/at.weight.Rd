\name{at.weight}
\alias{at.weight}
\title{Atomic weight and isotopic abundances of a chemical element from isotope ratios}
\description{
This function calculates the isotopic abundances and the atomic weight of a chemical element from the given isotope amount ratios and their uncertainties.
The uncertainty evaluation is done using the Monte Carlo method and the relevant masses of the isotopes are extracted from the www.ciaaw.org.
}
\usage{
at.weight(ratio, ratio.cov, element, ref.isotope)
}

\arguments{
  \item{ratio}{A vector of nontrivial isotope amount ratios of an element}
  \item{ratio.cov}{A covariance matrix of \code{ratio}}
  \item{element}{A string consisting of the lowercase English name of the element. For example, "zinc"}
  \item{ref.isotope}{A string which specifies the reference isotope. For example, "64Zn"}
}
\details{
The isotopic composition of an element with \code{N} isotopes is characterized using a set of \code{N-1} nontrivial isotope amount ratios.
As an example, silicon has three stable isotopes (silicon-28, silicon-29, and silicon-30) and its isotope ratios can be reported against
either of its stable isotopes in three distinct ways: (1) 29Si/28Si and 29Si/28Si or (2) 28Si/29Si and 30Si/29Si, or (3) 28Si/30Si and 29Si/30Si.
}
\value{
\item{aw}{Atomic weight}
\item{aw.u}{Standard uncertainty of the atomic weight}
\item{aw.U95}{Expanded uncertainty of the atomic weight corresponding to 95\% confidence}
\item{abundances}{Isotopic abundances}
\item{abundances.u}{Standard uncertainty of the isotopic abundances}
\item{abundances.U95}{Expanded uncertainty of the isotopic abundances corresponding to 95\% confidence}
\item{abundances.cov}{Covariance matrix of the isotopic abundances}
}
\references{
J.Meija and Z. Mester (2008) \href{http://dx.doi.org/10.1088/0026-1394/45/1/008}{Uncertainty propagation of atomic weight measurement results}. \emph{Metrologia}, 45, 53-62

J. Meija and A. Possolo (2017) \href{https://doi.org/10.1088/1681-7575/aa634d}{Data reduction framework for standard atomic weights and isotopic compositions of the elements}. \emph{Metrologia}, 54, 229-238

JCGM 101:2008 Evaluation of measurement data - Supplement 1 to the "Guide to the expression of uncertainty in measurement" - Propagation of distributions using a Monte Carlo method
}
\author{Juris Meija <juris.meija@nrc-cnrc.gc.ca> and Antonio Possolo}
\examples{
## Atomic weight and isotopic abundances of iridium which correspond
## to the isotope ratio 191Ir/193Ir = 0.59471(13)
at.weight(0.59471, matrix(0.00013^2), "iridium", "193Ir")

## Atomic weight and isotopic abundances of silicon which correspond
## to isotope ratios 28Si/29Si = 1.074(69) and 30Si/29Si = 260(11)
## with a correlation of 0.80 between the two isotope ratios
ratios = c(1.074,260)
r.cov = matrix(c(0.069^2,0.80*0.069*11,0.80*0.069*11,11^2),ncol=2,byrow=TRUE)
at.weight(ratios, r.cov, "silicon", "29Si")
}