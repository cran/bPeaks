\name{bPeaks-package}
\alias{bPeaks-package}
\alias{bPeaks}
\docType{package}
\title{
bPeaks: a simple and intuitive approach for detection of basic peaks (bPeaks) from ChIP-seq data
}
\description{
bPeaks is a simple approach to identify transcription factor binding sites from ChIP-seq data. Our general philosophy is to provide an easy-to-use tool. Instead of complex statistical models, bPeaks uses a combination of 4 cutoffs to mimic "good peak" properties as described by biologists who visually inspect the ChIP-seq data on a genome browser.
}
\details{
\tabular{ll}{
Package: \tab bPeaks\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-03-22\cr
License: \tab GPL \cr
}
}
\author{
Jawad MERHEJ and Gaelle LELANDAIS
Maintainer: Gaelle LELANDAIS <gaelle.lelandais@univ-paris-diderot.fr>
}
\references{
More information can be found online: \code{http://bpeaks.gene-networks.net/}
}
\keyword{peak calling}
\keyword{ChIP-seq}
\keyword{protein binding sites}
\keyword{protein-DNA interactions}
\keyword{deep sequencing}
\seealso{
\code{http://bpeaks.gene-networks.net/}
}
\examples{
# get library
library(bPeaks)

# STEP 1: get PDR1 data (ChIP-seq experiments - IP and control samples - 
# with transcription factor Pdr1 in yeast Saccharomyces cerevisiae) 
data(dataPDR1)

# STEP 2 : bPeaks analysis (only 10 kb of chrIV are analyzed here, 
#          as an illustration)
bPeaksAnalysis(IPdata = dataPDR1$IPdata[40000:50000,], 
                controlData = dataPDR1$controlData[40000:50000,], 
                windowSize = 150, windowOverlap = 50, 
                IPcoeff = 4, controlCoeff = 2, log2FC = 1, 
                averageQuantiles = 0.5,
                resultName = "bPeaks_example")

# --> Result files (PDF and BED) are written in the working directory.

\dontrun{
# -> bPeaks analysis (all chromosome 4 and default parameters optimized for yeasts)

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
