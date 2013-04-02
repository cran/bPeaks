\name{peakDrawing}
\alias{peakDrawing}
\title{
Function to draw graphical representations of genomic regions detected using bPeaks methodology
}
\description{
This function allows to create PDF files, with graphical representations of the detected basic peaks (bPeaks). Genomic regions are shown together with the values of the parameters used to detect the region. 
}
\usage{
peakDrawing(vecIP, vecControl, lineIP, lineControl, lineFC, lineAverage, 
            posInf = 1, posSup = NULL, add = 10, title = "")
}
\arguments{
  \item{vecIP}{
Vector with sequencing depth at each nucleotide (from start pos = 1 to end)
}
  \item{vecControl}{
Vector with sequencing depth at each nucleotide (from start pos = 1 to end)
}
  \item{lineIP}{
Threshold values used for peak detection (IP signal)
}
  \item{lineControl}{
Threshold values used for peak detection (control signal)
}
  \item{lineFC}{
Threshold values used for peak detection (log2(IP/control) values)
}
  \item{lineAverage}{
Threshold values used for peak detection (average log2(IP) and log2(control) values)
}
  \item{posInf}{
Genomic position to start the representation
}
  \item{posSup}{
Genomic position to end the representation
}
  \item{add}{
Number of bases before and after posInf and posSup to add
}
  \item{title}{
Graphic main title
}
}
\details{
More information can be found online: \code{http://bpeaks.gene-networks.net/}.
}
\value{
Image in x11() terminal
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
}
\examples{
# get library 
library(bPeaks)

# get PDR1 data
data(dataPDR1)

# IP signal (smoothed) - Chromosome IV
IPsignal = dataSmoothing(dataPDR1$IPdata[dataPDR1$IPdata[,1] == "chrIV",3], 20)
# control signal (smoothed)
controlSignal = dataSmoothing(dataPDR1$controlData[dataPDR1$controlData[,1] == "chrIV",3], 20)

# draw all chromosome
peakDrawing(vecIP = IPsignal, vecControl = controlSignal, 
            lineIP = 0, lineControl = 0, lineFC = 0, lineAverage = 0, 
            posInf = 465000, posSup = 465550, 
            add = 10, title = "PDR1 data - chromosome #4")

}
\keyword{bPeaks}

