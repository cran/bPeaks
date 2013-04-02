\name{peakDetection}
\alias{peakDetection}
\title{
Peak calling method, i. e. identification of genomic regions with a high density of sequences (reads)
}
\description{
bPeaks uses a sliding window to scan the genomic sequence. Four criterion define interesting regions: 1) a high number of reads in the IP sample (T1 = IPthreshold); 2) a low number of reads in the control sample (T2 = controlThreshold); 3) a high value of log(IP/control) (T3 = ratioThreshold) and 4) a good sequencing coverage (IP and control samples) (T4 = averageThreshold). Parameters for peak detection, therefore, rely on four threshold values (T1, T2, T3 and T4). Graphical outputs and BED files are provided allowing the user to rapidly assess the relevance of the chosen parameters.
}
\usage{
peakDetection(IPdata, controlData, chrName, windowSize = 150, windowOverlap = 50,
              outputName = "bPeaks_results", 
              baseLineIP = NULL, baseLineControl = NULL, 
              IPthreshold = 6, controlThreshold = 4, ratioThreshold = 2,
              averageThreshold = 0.7, peakDrawing = TRUE)
}
\arguments{
  \item{IPdata}{
A dataframe with sequencing results of the IP sample. This dataframe has three columns (chromosome, position, number of sequences) and should have been created with the \code{\link{dataReading}} function
}
  \item{controlData}{
A dataframe with sequencing results of the control sample. This dataframe has three columns (chromosome, position, number of sequences) and should have been created with the \code{\link{dataReading}} function
}
  \item{chrName}{
Name of the chromosome to be scanned with the sliding windows (to compare IP and control signals and detect interesting regions)
}
  \item{windowSize}{
Size of the sliding window to scan chromosomes
}
  \item{windowOverlap}{
Size of the overlap between two successive windows
}
  \item{outputName}{
Name for output files created during bPeaks procedure
}
  \item{baseLineIP}{
Value of the mean genome-wide read depth (IP sample). This value is calculated using the \code{\link{baseLineCalc}} function
}
  \item{baseLineControl}{
Value of the mean genome-wide read depth (control sample). This value is calculated using the \code{\link{baseLineCalc}} function
}
  \item{IPthreshold}{
Threshold T1. Threshold to consider IP signal as sufficiently important to be interesting. Note that these threshold is a multiplicative parameter that will be combined with the calculated "baseLineIP" (see before) value. As an illustration, if the IPthreshold = 6, it means that to be selected, the IP signal should be GREATER than 6 * baseLineIP
}
  \item{controlThreshold}{
Threshold T2. Threshold to consider control signal as sufficiently low to be interesting. Note that these threshold is a multiplicative parameter that will be combined with the calculated "baseLineControl" (see before) value. As an illustration, if the controlThreshold = 2, it means that to be selected, the control signal should remain LOWER than 2 * baseLineControl
}
  \item{ratioThreshold}{
Threshold T3. Threshold to consider log2(IP/control) values as sufficiently important to be interesting
}
  \item{averageThreshold}{
Threshold T4. Threshold to consider (log2(IP) + log2(control)) / 2 as sufficiently important to be interesting. These parameter is important to ensure that the analyzed genomic region has enough sequencing coverage to be reliable. These threshold should be between [0, 1] and refers to the quantile value of the global distribution observed with the analyzed chromosome  
}
  \item{peakDrawing}{
TRUE or FLASE. If TRUE, the function \code{\link{peakDrawing}} is called and PDF files with graphical representations of detected peaks are created.
}
}
\details{
Detailed description of the bPeaks procedure together with tutorials can be found online:
 
\code{http://bpeaks.gene-networks.net/}.
}
\value{
A matrix with genomic positions of the detected peaks. Summaries of parameter calculations and peak detection criteria are shown in PDF files (saved in the working directory). 
}
\references{
\code{http://bpeaks.gene-networks.net/}
}
\author{
Gaelle LELANDAIS
}
\note{
Detailed information and tutorials can be found online \code{http://bpeaks.gene-networks.net/}. Don't hesitate to contact us for further discussions.
}
\seealso{
\code{\link{dataReading}}
\code{\link{baseLineCalc}}
\code{\link{peakDrawing}}
\code{\link{bPeaksAnalysis}}
}
\examples{
# get library
library(bPeaks)

# get PDR1 data
data(dataPDR1)

# combine IP and control data
allData = cbind(dataPDR1$IPdata, dataPDR1$controlData)
colnames(allData) = c("chr", "pos", "IPsignal", "chr", "pos", "controlSignal")

print("**********************************************")
# calculate baseline IP and control values
lineIP    = baseLineCalc(allData$IPsignal)
print(paste("Baseline coverage value in IP sample : ", round(lineIP, 3)))
lineControl = baseLineCalc(allData$controlSignal)
print(paste("Baseline coverage value in control sample : ", round(lineControl, 3)))
print("**********************************************")
print("")

# get list of chromosomes
chromNames = unique(allData[,1])

# start peak detection on the first chromosome
print("**********************************************")
print(paste("Starting analysis of chromosome ", chromNames[1]))

# information for one chromosome
subData = allData[allData[,1] == chromNames[1],]

# only 10 kb are analyzed here (as an illustration)
vecIP      = subData[40000:50000,3]
vecControl = subData[40000:50000,6]

# smooth of the data
smoothedIP    = dataSmoothing(vecData = vecIP, widthValue = 20) 
smoothedControl = dataSmoothing(vecData = vecControl, widthValue = 20) 

# peak detection
detectedPeaks = peakDetection(IPdata = smoothedIP, controlData = smoothedControl, 
                chrName = as.character(chromNames[1]), 
                windowSize = 150, windowOverlap = 50,
                outputName = paste("bPeaks_example_", chromNames[1], sep = ""), 
                baseLineIP = lineIP, baseLineControl = lineControl,
                IPthreshold = 4, controlThreshold = 2,
                ratioThreshold = 1, averageThreshold = 0.5,
                peakDrawing = TRUE)

# print detected genomic positions
print(detectedPeaks)
}
\keyword{peak calling}
\keyword{ChIP-seq}
\keyword{protein binding sites}
\keyword{protein-DNA interactions}
\keyword{deep sequencing}
