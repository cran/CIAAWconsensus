\name{mmm}
\alias{mmm}
\title{Multivariate meta-analysis of correlated effects}
\description{
This function provides meta-analysis of multivariate correlated data using the marginal method of moments with working independence assumption as described by Chen et al (2016).
As such, the meta-analysis does not require correlations between the outcomes within each dataset.
}
\usage{
mmm(y, uy, knha = TRUE, verbose = TRUE)
}
\arguments{
\item{y}{A matrix of results from each of the \code{n} laboratories (rows) where each study reports \code{m} isotope ratios (columns)
}
\item{uy}{A matrix with uncertainties of the results given in \code{y}
}
\item{knha}{
(Logical) Allows for the adjustment of consensus uncertainties using the Birge ratio (Knapp-Hartung adjustment)
}
\item{verbose}{
(Logical) Requests annotated summary output of the results
}
}
\details{
The marginal method of moments delivers the inference for correlated effect sizes using multiple univariate meta-analyses.
}
\value{
\item{studies}{The number of independent studies}
\item{beta}{The consensus estimates for all outcomes}
\item{beta.u}{Standard uncertainties of the consensus estimates}
\item{beta.U95}{Expanded uncertainties of the consensus estimates corresponding to 95\% confidence}
\item{beta.cov}{Covariance matrix of the consensus estimates}
\item{beta.cor}{Correlation matrix of the consensus estimates}
\item{H}{Birge ratios (Knapp-Hartung adjustment) which were applied to adjust the standard uncertainties of each consensus outcome}
\item{I2}{Relative total variability due to heterogeneity (in percent) for each outcome}
}
\references{
Y. Chen, Y. Cai, C. Hong, and D. Jackson (2016) \href{http://dx.doi.org/10.1002/sim.6789}{Inference for correlated effect sizes using multiple univariate meta-analyses}. \emph{Statistics in Medicine}, 35, 1405-1422

J. Meija and A. Possolo (2017) \href{https://doi.org/10.1088/1681-7575/aa634d}{Data reduction framework for standard atomic weights and isotopic compositions of the elements}. \emph{Metrologia}, 54, 229-238
}
\author{Juris Meija <juris.meija@nrc-cnrc.gc.ca> and Antonio Possolo}

\examples{
## Consensus isotope amount ratios for platinum
df=normalize.ratios(platinum.data, "platinum", "195Pt")
mmm(df$R, df$u.R)
}