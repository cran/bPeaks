\name{bPeaksAnalysis}
\alias{bPeaksAnalysis}
\title{
Function to run the entire bPeaks procedure
}
\description{
This function allows to detect basic peaks (bPeaks) using the procedure described in the function \code{\link{peakDetection}}. Chromosomes are analyzed successively. Several values(regarding thresholds T1, T2, T3 and T4 and other parameters) can be specified simultaneously in order to rapidly compare the obtained results and evaluate parameter relevance.
}
\usage{
bPeaksAnalysis(IPdata, controlData, cdsPositions = NULL, 
		smoothingValue = 20,
               	windowSize = 150, windowOverlap = 50, 
		IPcoeff = 2, controlCoeff = 2,
               	log2FC = 2, averageQuantiles = 0.9, 
		resultName = "bPeaks", 
		peakDrawing = TRUE, promSize = 800, withoutOverlap = FALSE)
}
\arguments{
  \item{IPdata}{
A dataframe with sequencing results of the IP sample. This dataframe has three columns (chromosome, position, number of sequences) and should have been created with the \code{\link{dataReading}} function
}
  \item{controlData}{
A dataframe with sequencing results of the control sample. This dataframe has three columns (chromosome, position, number of sequences) and should have been created with the \code{\link{dataReading}} function
}
  \item{cdsPositions}{
Not mandatory. A table (matrix) with positions of CDS (genes). Four columns are required (chromosome, starting position, ending position, strand (W or C), description). CDS positions for several yeast species are stored in bPeaks package 
(see the dataset \code{\link{yeastCDS}} and also \code{\link{peakLocation}} function) 
}
  \item{smoothingValue}{
The number (n/2) of surrounding positions to use for mean calculation in the \code{\link{dataSmoothing}} function
}
  \item{windowSize}{
Size of the sliding windows to scan chromosomes
}
  \item{windowOverlap}{
Size of the overlap between two successive windows
}
  \item{IPcoeff}{
Threshold T1. Value for the multiplicative parameter that will be combined with the value of the mean genome-wide read depth (see \code{\link{baseLineCalc}}). As an illustration, if the IPcoeff = 6, it means that to be selected, the IP signal should be GREATER than 6 * (the mean genome-wide read depth). Note that a vector with different values can be specified, the bPeaks analysis will be therefore repeated using successively each value for peak detection 
}
  \item{controlCoeff}{
 Threshold T2. Value for the multiplicative parameter that will be combined with the value of the mean genome-wide read depth (see \code{\link{baseLineCalc}}). As an illustration, if the controlCoeff = 2, it means that to be selected, the control signal should be LOWER than 2 * (the mean genome-wide read depth). Note that a vector with different values can be specified, the bPeaks analysis will be therefore repeated using successively each value for peak detection
}
  \item{log2FC}{
Threshold T3. Threshold to consider log2(IP/control) values as sufficiently important to be interesting. Note that a vector with different values can be specified, the bPeaks analysis will be therefore repeated using successively each value for peak detection
}
  \item{averageQuantiles}{
Threshold T4. Threshold to consider (log2(IP) + log2(control)) / 2 as sufficiently important to be interesting. This parameter ensures that the analyzed genomic region has enough sequencing coverage to be reliable. These threshold should be between [0, 1] and refers to the quantile value of the global distribution observed with the analyzed chromosome 
}
  \item{resultName}{
Name for output files created during bPeaks procedure
}
  \item{peakDrawing}{
TRUE or FLASE. If TRUE, the function \code{\link{peakDrawing}} is called and PDF files with graphical representations of detected peaks are created.
}
  \item{promSize}{
Size of the genomic regions to be considered as "upstream" to the annotated genomic features (see documentation of the function \code{\link{peakLocation}} for more information).
}
  \item{withoutOverlap}{
If TRUE, this option allows to filter peak that are located in a promoter AND a CDS.
}
}
\details{
More information together with tutorials can be found online \code{http://bpeaks.gene-networks.net/}.
}
\value{
BED files for each chromosomes and a final BED file combining all the results with information regarding detected peaks (genomic positions, mean IP signal, etc.). These files are all saved in the R working directory. Summaries of parameter calculations and peak detection criteria are shown in PDF files (saved in the working directory). 
}
\references{
\code{http://bpeaks.gene-networks.net/}
}
\author{
Gaelle LELANDAIS
}
\note{
Detailed information and tutorials can be found online \code{http://bpeaks.gene-networks.net/}
}
\seealso{
\code{\link{peakDetection}}
\code{\link{dataReading}}
\code{\link{dataSmoothing}}
\code{\link{baseLineCalc}}
\code{\link{peakDrawing}}
\code{\link{peakLocation}}
}
\examples{
# get library
library(bPeaks)

# STEP 1: get PDR1 data
data(dataPDR1)

# STEP 2 : bPeaks analysis (only 10 kb of chrIV are analyzed here, 
#          as an illustration)
bPeaksAnalysis(IPdata = dataPDR1$IPdata[40000:50000,], 
               controlData = dataPDR1$controlData[40000:50000,], 
               windowSize = 150, windowOverlap = 50, 
               IPcoeff = 4, controlCoeff = 2, log2FC = 1, 
               averageQuantiles = 0.5,
               resultName = "bPeaks_example", 
               peakDrawing = TRUE, promSize = 800)

\dontrun{
# STEP 2 : bPeaks analysis (all chromosome)
bPeaksAnalysis(IPdata = dataPDR1$IPdata, controlData = dataPDR1$controlData, 
                cdsPositions = dataPDR1$cdsPositions, 
                smoothingValue = c(20), 
                windowSize = c(150), windowOverlap = 50, 
                IPcoeff = c(2), controlCoeff = c(2), log2FC = c(2), 
                averageQuantiles = c(0.9),
                resultName = "bPeaks_PDR1_chr4", 
                peakDrawing = TRUE, promSize = 800)

# To repeat the bPeaks analysis with different parameters
bPeaksAnalysis(IPdata = dataPDR1$IPdata, controlData = dataPDR1$controlData, 
                cdsPositions = dataPDR1$cdsPositions, 
                smoothingValue = c(20), 
                windowSize = c(150), windowOverlap = 50, 
                IPcoeff = c(2, 4, 6), controlCoeff = c(2, 4, 6), log2FC = c(2, 3), 
                averageQuantiles = c(0.7, 0.9),
                resultName = "bPeaks_PDR1_chr4_paremeterEval", 
                peakDrawing = FALSE, promSize = 800)

# -> Summary table is created and saved as "peakStats.Robject" in the working directory
# as well as a text file named "_bPeaks_parameterSummary.txt"...
load("peakStats.Robject")
# This table comprises different information regarding peak detection (number of peaks,
# mean size of peaks, mean IP signal, mean control signal, etc.)
peakStats[1:2,]

#     smoothingValue windowSize windowOverlap IPcoeff controlCoeff log2FC
# [1,]             20        150            50       1            1      1
# [2,]             20        150            50       1            1      1
#     averageQuantiles bPeakNumber meanSize meanIPsignal meanControlSignal
# [1,]              0.5         308  209.091      276.047            71.534
# [2,]              0.7         294  205.782      287.808            74.002
#     meanLog2FC bPeakNumber_beforeFeatures bPeakNumber_afterFeatures
# [1,]      1.571                         99                        80
# [2,]      1.589                         94                        77
#     bPeakNumber_inFeatures
# [1,]                     52
# [2,]                     53

}
}
\keyword{peak calling}
\keyword{ChIP-seq}
\keyword{protein binding sites}
\keyword{protein-DNA interactions}
\keyword{deep sequencing}
