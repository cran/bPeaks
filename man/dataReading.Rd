\name{dataReading}
\alias{dataReading}
\title{
Function to import sequencing results in R
}
\description{
Sequencing results should be converted to datafiles with numbers of sequences mapped on each nucleotide in the reference genome. These files can be generated from indexed BAM files (mapping results) using the "genomeCoverageBed" tool available from the BEDTOOLS suite [1]. 
More information concerning file descriptions can be found online: 

\code{http://bpeaks.gene-networks.net/}.

File format should be as follows (chromosome, position, number of sequences):
\tabular{lll}{
chrI\tab	1\tab4\cr
chrI\tab	2\tab4\cr
chrI\tab	3\tab4\cr
chrI\tab	4\tab4\cr
chrI\tab	5\tab7\cr
chrI\tab	6\tab7\cr
chrI\tab7\tab9\cr
chrI\tab	8\tab9\cr
chrI\tab	9\tab10\cr
chrI\tab	10\tab13\cr
}
}
\usage{
dataReading(IPfile, controlFile, yeastSpecies = NULL)
}
\arguments{
  \item{IPfile}{
Name of the file with sequencing results related to IP sample
}
  \item{controlFile}{
Name of the file with sequencing results related to control sample
}
  \item{yeastSpecies}{
Not mandatory. Annotations to be used for locations of peaks in promoters.
Annotations of CDS are available in bPeaks for yeasts: Debaryomyces.hansenii, 
Eremothecium.gossypii, Kluyveromyces.lactis, Pichia.sorbitophila, 
Saccharomyces.kluyveri, Yarrowia.lipolytica, Zygosaccharomyces.rouxii,
Saccharomyces.cerevisiae, Candida.albicans, Candida.glabrata (see data {\link{yeastCDS}})
}
}
\details{
To obtain a required file from a BAM file (resultFile.bam), the command line is (SHELL): 

genomeCoverageBed -ibam resultFile.bam -d > resultFile.txt 

More information concerning file conversions can be found online: 

\code{http://bpeaks.gene-networks.net/}.
}
\value{
A list with three elements ($IPdata, $controlData, $cdsPositions): IP data, control data and (if specified by user) CDS positions for locations of peaks in promoters.
}
\references{
[1] Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 26, 6, pp. 841 842.
}
\author{
Gaelle LELANDAIS
}
\note{
Conversion of file formats regarding sequencing results can be a tricky task. Detailed information can be found online \code{http://bpeaks.gene-networks.net/}. Don't hesitate to contact us for further discussions.
}
\seealso{
\code{\link{peakLocation}}
\code{\link{dataPDR1}}
\code{\link{yeastCDS}}
}
\examples{
# get library
library(bPeaks)

\dontrun{
# Sequencing result files associated to PDR1 datasets (IP and control files) 
# can be downloaded from our website http://bpeaks.gene-networks.net/. 
# They are respectively named "IP_genomeCoverage.txt" (IP sample) and 
# "INPUT_genomeCoverage.txt" (control sample).

# Import in R the sequencing results using S. cerevisiae CDS annotations.
data(yeastCDS)
seqResult = dataReading("IP_genomeCoverage.txt", "INPUT_genomeCoverage.txt",
			yeastSpecies = yeastCDS$Saccharomyces.cerevisiae)

# IP data
seqResult$IPdata

# control data
seqResult$controlData

# run peak detection from IP and control data (with default parameters)
bPeaksAnalysis(IPdata = seqResult$IPdata, 
	       controlData  = seqResult$controlData,
	       cdsPositions = seqResult$cdsPositions)
}
}
\keyword{data reading}

