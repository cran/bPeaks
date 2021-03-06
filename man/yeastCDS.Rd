\name{yeastCDS}
\alias{yeastCDS}
\docType{data}
\title{
Annotations of CDS for different yeast species
}
\description{
Annotations of CDS for different yeast species. Data were collected (July 2013) from: 
Genolevures \code{http://genolevures.org/}, 
SGD \code{http://www.yeastgenome.org/},
CGD \code{http://www.candidagenome.org/} databases.
}
\usage{data(yeastCDS)}
\format{
\code{yeastCDS$Saccharomyces.cerevisiae}:
A table (matrix) with positions of CDS (genes) in yeast S. cerevisiae. Four columns are required (chromosome, starting position, ending position, strand (W or C), description).
}
\details{
This R object is a list composed of multiple tables (matrices) - one for each species - with positions 
of CDS (genes). 
}
\source{
Data were collected (July 2013) from 
Genolevures \code{http://genolevures.org/}, SGD \code{http://www.yeastgenome.org/} and 
CGD \code{http://www.candidagenome.org/} databases.
}
\references{
Genolevures database \code{http://genolevures.org/}

SGD database \code{http://www.yeastgenome.org/}

CGD database \code{http://www.candidagenome.org/}
}
\examples{
# get library
library(bPeaks)

# get data
data(yeastCDS)

# different species for wich information is available
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

# CDS annotations for yeast Debaryomyces hansenii
yeastCDS$Debaryomyces.hansenii[1:10,]
#      chrm     start   end     strand geneName      
# [1,] "Deha2A" "2023"  "6370"  "C"    "DEHA2A00110g"
# [2,] "Deha2A" "6587"  "7810"  "C"    "DEHA2A00132g"
# [3,] "Deha2A" "8314"  "9354"  "W"    "DEHA2A00154g"
# [4,] "Deha2A" "9632"  "9844"  "C"    "DEHA2A00176g"
# [5,] "Deha2A" "13806" "14132" "W"    "DEHA2A00198g"
# [6,] "Deha2A" "14558" "16519" "C"    "DEHA2A00220g"
# [7,] "Deha2A" "17520" "19442" "W"    "DEHA2A00242g"
# [8,] "Deha2A" "22619" "23977" "C"    "DEHA2A00264g"
# [9,] "Deha2A" "24949" "25434" "W"    "DEHA2A00286g"
#[10,] "Deha2A" "26440" "26640" "C"    "DEHA2A00308g"
}
\keyword{Pdr1 transcription factor}
\keyword{Saccharomyces cerevisiae}
\keyword{ChIP-seq results}
