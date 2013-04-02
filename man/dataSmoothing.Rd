\name{dataSmoothing}
\alias{dataSmoothing}
\title{
Function to smooth sequencing coverage along a chromosome
}
\description{
This function allows to obtain a smoothed signal of the genome-wide read depth. Simple moving average (SMA) procedure is used. At each genomic position, the sequencing coverage is replaced by the unweighted mean of the n surrounding positions (n/2 before and n/2 after).  
}
\usage{
dataSmoothing(vecData, widthValue = 20)
}
\arguments{
  \item{vecData}{
A vector with the numbers of sequences aligned at each genomic position to be considered in the analysis
}
  \item{widthValue}{
The number (n/2) of surrounding positions to use for mean calculation
}
}
\details{
Detailed information and tutorials can be found online \code{http://bpeaks.gene-networks.net/}.
}
\value{
A vector with the smoothed signal. Note that the SMA procedure creates missing values at the beginning and at the end of the vector with the smoothed signal.
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
\code{\link{baseLineCalc}}
}
\examples{
# get data
data(dataPDR1)

# inital IP signal
iniIPsignal = dataPDR1$IPdata[,3]

par(mfrow = c(2,2))

# plot initial IP signal
plot(iniIPsignal[416900:417400], type = "h", 
    xlab = "genomic position", ylab = "sequencing coverage",
    main = "IP sample (PDR1 data)\nno smoothing", 
    col = "red")

# calculate and plot smoothed signal (widthValue = 5)
smoothedIPsignal = dataSmoothing(vecData = iniIPsignal, widthValue = 5)
plot(smoothedIPsignal[416900:417400], type = "h", 
    xlab = "genomic position", ylab = "sequencing coverage",
    main = "IP sample (PDR1 data)\nsmoothing (widthValue = 5)", 
    col = "pink")

# calculate and plot smoothed signal (widthValue = 10)
smoothedIPsignal = dataSmoothing(vecData = iniIPsignal, widthValue = 10)
plot(smoothedIPsignal[416900:417400], type = "h", 
    xlab = "genomic position", ylab = "sequencing coverage",
    main = "IP sample (PDR1 data)\nsmoothing (widthValue = 10)", 
    col = "pink")

# calculate and plot smoothed signal (widthValue = 20)
smoothedIPsignal = dataSmoothing(vecData = iniIPsignal, widthValue = 20)
plot(smoothedIPsignal[416900:417400], type = "h", 
    xlab = "genomic position", ylab = "sequencing coverage",
    main = "IP sample (PDR1 data)\nsmoothing (widthValue = 20)", 
    col = "pink")
}
\keyword{moving average}
\keyword{signal smoothing}

