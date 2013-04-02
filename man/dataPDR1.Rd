\name{dataPDR1}
\alias{dataPDR1}
\docType{data}
\title{
ChIP-seq results (IP and control samples) obtained with the transcription factor Pdr1 in yeast Saccharomyces cerevisiae
}
\description{
ChIP-seq experiments were performed in order to identify the genomic regions that interact with the transcription factor Pdr1, in yeast Saccharomyces cerevisiae. Two samples (IP and control) were sequenced simultaneously using the Illumina technology (ENS transcriptome platform). Only the data for chr14 are available here, but complete datasets can be found online: \code{http://bpeaks.gene-networks.net}
}
\usage{data(dataPDR1)}
\format{
\code{dataPDR1$IPdata}:
IPdata and controlData are dataframes with three columns. The first column comprises chromosome names, the second column comprises base positions and the third column comprises the numbers of sequences mapped at the considered position.
\code{dataPDR1$controlData}:
IPdata and controlData are dataframes with three columns. The first column comprises chromosome names, the second column comprises base positions and the third column comprises the numbers of sequences mapped at the considered position.
\code{dataPDR1$chromosomalFeatures}:
chromosomalFeatures is a table with annotated positions of genes in yeast S. cerevisiae. The first column indicates chromosome names, the second and third columns indicate respectively "start" and "end" positions of genes, and the fourth column indicates the gene annotation (according to the Saccharomyces Genome Database (SGD \code{http://www.yeastgenome.org/})).
}
\details{
Complete procedure to analyze sequencing data (intial FASTQ files) can be found 
online:
\code{http://bpeaks.gene-networks.net}. 
Initial read length was 50 bases. After quality controls and filtering of low quality bases, around 30.000.000 of sequences (IP sample) and around 88.000.000 of sequences (control sample) were independantly mapped on the genome using the bowtie algorithm [1]. Output files (SAM format) were converted into BAM files and indexed using the SAMTOOLS suite [2]. Numbers of sequences mapped on each nucleotide in the reference genome were finally calculated using the "genomeCoverageBed" tool available from the BEDTOOLS suite [3].    
}
\source{
Sample sequencing was performed at the "transcriptome platform", ENS institute in Paris (France), \code{http://www.transcriptome.ens.fr}.
}
\references{
More information concerning this dataset can found online : \code{http://bpeaks.gene-networks.net}.

[1] Langmead B, Trapnell C, Pop M, Salzberg SL. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biol 10:R25.

[2] Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078 9. [PMID: 19505943]

[3] Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 26, 6, pp. 841 842.
}
\examples{
# get library
library(bPeaks)

# get data
data(dataPDR1)

summary(dataPDR1)
#                    Length Class      Mode     
#IPdata                  3  data.frame list     
#controlData             3  data.frame list     
#chromosomalFeatures 26412  -none-     character

# run peak calling, comparing IP and control samples
# (only 10 kb of chrIV are analyzed here, as an illustration)
bPeaksAnalysis(IPdata = dataPDR1$IPdata[40000:50000,], 
                controlData = dataPDR1$controlData[40000:50000,], 
                windowSize = 150, windowOverlap = 50, 
                IPcoeff = 4, controlCoeff = 2, log2FC = 1, 
                averageQuantiles = 0.5,
                resultName = "bPeaks_example", 
                peakDrawing = TRUE, promSize = 800)

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
resLocation = peakLocation(bedFile = "bPeaks_PDR1_bPeaks_allGenome.bed", 
                           genomicInfo = dataPDR1$chromosomalFeatures,
                           outputName = "bPeakLocation_finalPDR1", promSize = 800)

# number of detected peaks
print(resLocation$numPeaks)

# number of peaks "before" annotated genes
print(resLocation$beforeFeatures)

# number of peaks "in" annotated genes
print(resLocation$inFeatures)

# number of peaks "after" annotated genes
print(resLocation$afterFeatures)
}
}
\keyword{Pdr1 transcription factor}
\keyword{Saccharomyces cerevisiae}
\keyword{ChIP-seq results}
