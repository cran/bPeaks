\name{baseLineCalc}
\alias{baseLineCalc}
\title{
Function to calculate the average number of reads mapped on each nucleotide in a genome 
}
\description{
This function calculates the mean genome-wide read depth. 
}
\usage{
baseLineCalc(covData)
}
\arguments{
  \item{covData}{
A vector with the numbers of sequences aligned at each genomic position to be considered in the analysis 
}
}
\details{
This function adds the numbers of sequences observed at each position and divides this number by the genome size (total number of nucleotides).
}
\value{
The average number of reads mapped on each nucleotide in the genome.
}
\references{
\code{http://bpeaks.gene-networks.net/}
}
\author{
Gaelle LELANDAIS
}
\note{
Detailed information and tutorials can be found online \code{http://bpeaks.gene-networks.net/}.
}
\seealso{
\code{\link{bPeaksAnalysis}}
\code{\link{peakLocation}}
}
\examples{
# get library
library(bPeaks)

# get PDR1 data
data(dataPDR1)

# mean genome-wide read depth in IP data
meanIPcov = baseLineCalc(dataPDR1$IPdata[,3])
print(meanIPcov)

# mean genome-wide read depth in control data
meanContCov = baseLineCalc(dataPDR1$controlData[,3])
print(meanContCov) 
}
\keyword{genome coverage}
\keyword{read depth}


