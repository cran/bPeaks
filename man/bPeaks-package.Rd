\name{bPeaks-package}
\alias{bPeaks-package}
\alias{bPeaks}
\docType{package}
\title{
bPeaks: an intuitive peak-calling strategy to detect transcription factor binding sites from ChIP-seq data in small eukaryotic genomes
}
\description{
bPeaks is a simple approach to identify transcription factor binding sites from ChIP-seq data. Our general philosophy is to provide an easy-to-use tool, well-adapted for small eukaryotic genomes (< 20 Mb). bPeaks uses a combination of 4 cutoffs (T1, T2, T3 and T4) to mimic "good peak" properties as described by biologists who visually inspect the ChIP-seq data on a genome browser. For yeast genomes, bPeaks calculates the proportion of peaks that fall in promoter sequences. These peaks are good candidates as transcription factor binding sites. 
}
\details{
\tabular{ll}{
Package: \tab bPeaks\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2014-02-28\cr
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
\keyword{small eukaryotic genomes}
\seealso{
\code{http://bpeaks.gene-networks.net/}
}
\examples{
# get library
library(bPeaks)

# STEP 1: get PDR1 data (ChIP-seq experiments, IP and control samples, 
# 	related to the transcription factor Pdr1 in yeast Saccharomyces 
# 	cerevisiae) 
data(dataPDR1)

# STEP 2 : bPeaks analysis (only 10 kb of chrIV are analyzed here, 
#          as an illustration)
bPeaksAnalysis(IPdata = dataPDR1$IPdata[40000:50000,], 
               controlData  = dataPDR1$controlData[40000:50000,], 
	       cdsPositions = dataPDR1$cdsPositions,
               windowSize = 150, windowOverlap = 50, 
               IPcoeff = 4, controlCoeff = 2, 
	       log2FC  = 1, averageQuantiles = 0.5,
               resultName = "bPeaks_example")

# --> Result files (PDF and BED) are written in the working directory.

\dontrun{
# -> bPeaks analysis, all chromosome IV and default parameters (optimized for yeasts)

# STEP 1: get PDR1 data (ChIP-seq experiments, IP and control samples, 
# 	related to the transcription factor Pdr1 in yeast Saccharomyces 
# 	cerevisiae) 
data(dataPDR1)

# STEP 2: bPeaks analysis
bPeaksAnalysis(IPdata = dataPDR1$IPdata, 
               controlData  = dataPDR1$controlData,
	       cdsPositions = dataPDR1$cdsPositions, 
               windowSize = 150, windowOverlap = 50, 
               IPcoeff = 2, controlCoeff = 2, 
               log2FC = 2, averageQuantiles = 0.9,
               resultName = "bPeaks_PDR1", 
               peakDrawing = TRUE)

# STEP 3 : procedure to locate peaks according to 
# 	   gene positions
peakLocation(bedFile = "bPeaks_PDR1_bPeaks_allGenome.bed", 
            cdsPositions = yeastCDS$Saccharomyces.cerevisiae,
            outputName = "bPeakLocation_finalPDR1", promSize = 800)

# -> Note that cds (genes) positions are stored in bPeaks package for several yeast
# species
data(yeastCDS)

summary(yeastCDS)
#                         Length Class      Mode     
# Debaryomyces.hansenii    31370  -none-     character
# Eremothecium.gossypii    23615  -none-     character
# Kluyveromyces.lactis     25380  -none-     character
# Pichia.sorbitophila      55875  -none-     character
# Saccharomyces.kluyveri   27790  -none-     character
# Yarrowia.lipolytica      32235  -none-     character
# Zygosaccharomyces.rouxii 24955  -none-     character
# Saccharomyces.cerevisiae     5  data.frame list     
# Candida.albicans             5  data.frame list     
# Candida.glabrata             5  data.frame list
}

}
