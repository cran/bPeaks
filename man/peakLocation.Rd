\name{peakLocation}
\alias{peakLocation}
\title{
Function to locate detected basic peaks (bPeaks) according to predifined chromosomal features
}
\description{
Starting from a BED file with positions of detected peaks and a file with positions of particular chromosomal features (genes for instance), this function allows to dentify the peaks that are located "before", "after" or "in" the chromosomal features.
}
\usage{
peakLocation(bedFile, chromosomalFeatureFile = "", genomicInfo = NULL, 
            outputName = "bPeakLocation", promSize = 800)
}
\arguments{
  \item{bedFile}{
Name of a BED file with positions of detected peaks (using bPeaks or another program)
}
  \item{chromosomalFeatureFile}{
Name of a file with positions of particular chromosomal features. This file must comprise four columns (chromosome, starting position, ending position, description)
}
  \item{genomicInfo}{
A table with four columns (chromosome, starting position, ending position, description). This table should have  been created with the \code{\link{dataReading}} function  
}
  \item{outputName}{
Name for output files
}
  \item{promSize}{
Genomic size to be considered "before" and "after" the chromosomal features 
}
}
\details{
More information can be found online \code{http://bpeaks.gene-networks.net/}.
}
\value{
Graphics and text files (saved in the R working directory).
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
\code{\link{dataReading}}
}
\examples{
\dontrun{
# -> bPeaks analysis with (all chromosome and default parameters optimized for yeasts)

# STEP 1: get PDR1 data (ChIP-seq experiments - IP and control samples - 
# with transcription factor Pdr1 in yeast Saccharomyces cerevisiae) 
data(dataPDR1)

# STEP 2: bPeaks analysis
bPeaksAnalysis(IPdata = dataPDR1$IPdata, 
               controlData = dataPDR1$controlData, 
               windowSize = 150, windowOverlap = 50, 
               IPcoeff = 6, controlCoeff = 4, 
               log2FC = 2, averageQuantiles = 0.9,
               resultName = "bPeaks_PDR1", 
               peakDrawing = TRUE)

# STEP 3 : procedure to locate peaks according to 
# predefined chromosomal features
peakLocation(bedFile = "bPeaks_PDR1_bPeaks_allGenome.bed", 
            genomicInfo = dataPDR1$chromosomalFeatures,
            outputName = "bPeakLocation_finalPDR1", promSize = 800)
}
}
