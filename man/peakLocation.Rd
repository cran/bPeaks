\name{peakLocation}
\alias{peakLocation}
\title{
Function to locate detected basic peaks (bPeaks) according to predifined chromosomal features
}
\description{
Starting from a BED file with positions of detected peaks and a table with positions of CDS (genes), this function allows to dentify the peaks that are located "upstream" or "in" annotated CDS. Annotations of CDS for different yeast species are available in bPeaks package (see data \code{\link{yeastCDS}}).
}
\usage{
peakLocation(bedFile, cdsPositions, withoutOverlap = FALSE,
            outputName = "bPeaksLocation", promSize = 800)
}
\arguments{
  \item{bedFile}{
Name of a BED file with positions of detected peaks (using bPeaks or another program)
}
  \item{cdsPositions}{
A table (matrix) with positions of CDS (genes). Four columns are required (chromosome, starting position, ending position, strand (W or C), description)
}
  \item{withoutOverlap}{
If TRUE, this option allows to filter peak that are located in a promoter AND a CDS.
}
  \item{outputName}{
Name for output files
}
  \item{promSize}{
Genomic size to be considered as promoter (upstream to CDS) 
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
\code{\link{yeastCDS}}
}
\examples{
\dontrun{
# -> bPeaks analysis with (all chromosome and default parameters optimized for yeasts)

# STEP 1: get PDR1 data and annotations in yeasts
data(dataPDR1)
data(yeastCDS)

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
            cdsPositions = yeastCDS$Saccharomyces.cerevisiae,
	    withoutOverlap = FALSE,
            outputName = "bPeakLocation_finalPDR1", promSize = 800)
}
}
