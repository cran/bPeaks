peakLocation <-
function(bedFile, chromosomalFeatureFile = "", genomicInfo = NULL, 
                        outputName = "bPeakLocation", promSize = 800){

#######
# written by G. LELANDAIS <gaelle.lelandais@univ-paris-diderot.fr>    
#######


    # argument names were changed for user convenience
    if(!is.null(genomicInfo)){
        ORFinfo = genomicInfo
    }else if(chromosomalFeatureFile != ""){
        print("Reading chromosomal features for peak location... ")
        ORFinfo = as.matrix(read.table(chromosomalFeatureFile, sep = "\t"))
        print("... done")
    }else{
        print("ERROR : No file for chromosomalFeatures is specified")
    }

    print("**********************************************")
    # data reading
    print("Opening BED file with peak information:")
    print(bedFile) 
    peakRes = as.matrix(read.table(bedFile, sep = "\t")) 
    print("**********************************************")

    # to store all the results    
    locRes = NULL

    # statitics
    numPeaks = nrow(peakRes)
    numPromD  = 0
    numPromI  = 0
    numORF    = 0

    print("")
    print("Starting peak location regarding chromosomal features...")
    print("")
 
    # for each detected peak....
    for(i in 1:nrow(peakRes)){

        currentPeakRes = NULL

        # chromosome
        chrmInfo = peakRes[i,1]
        peakStart = as.numeric(peakRes[i,2])
        peakEnd   = as.numeric(peakRes[i,3])        
        # the middle position will be used to map the peak on the genome
        middlePos = mean(c(peakEnd, peakStart))

        # get ORF on the analysed chromosome
        subORFinfo = ORFinfo[ORFinfo[,1] == chrmInfo,]
    
        if(nrow(subORFinfo) > 0){
            
            promDirect = (middlePos < as.numeric(subORFinfo[,2])) & 
                  (middlePos > (as.numeric(subORFinfo[,2]) - promSize))

            promIndirect = (middlePos > as.numeric(subORFinfo[,3])) & 
                  (middlePos < (as.numeric(subORFinfo[,3]) + promSize))
                
            ORFloc = (middlePos > as.numeric(subORFinfo[,2])) &
                     (middlePos < as.numeric(subORFinfo[,3]))  

            numElements = sum(promDirect) + sum(promIndirect) + sum(ORFloc)

       #     print(paste("number of detected elements for", 
                       # peakRes[i,4], "=", numElements))    

            if(sum(promDirect) > 0){
                locRes = c(locRes, paste(peakRes[i,4], subORFinfo[promDirect,4], sep = "\tbefore\t"))
            }
    
            if(sum(promIndirect) > 0){
                locRes = c(locRes, paste(peakRes[i,4], subORFinfo[promIndirect,4], sep = "\tafter\t"))
            }

            if(sum(ORFloc) > 0){
                locRes = c(locRes, paste(peakRes[i,4], subORFinfo[ORFloc,4], sep = "\tin\t"))
            }
                           
            # number of detected elements
             numPromD  = numPromD + sum(promDirect)
             numPromI  = numPromI + sum(promIndirect)
             numORF    = numORF + sum(ORFloc)

        }else{
            print(paste("WARNING : chromosome", chrmInfo, "is unknown..."))
        }

    # end of for()
    }
   
    # statistics printing
    print(paste("# of analyzed peaks: ", numPeaks))
    print(paste("# of peaks BEFORE annotated elements : ", numPromD))
    print(paste("# of peaks AFTER annotated elements : ", numPromI))
    print(paste("# of peaks IN annotated elements : ", numORF))
    print("")

    pdf(paste(outputName, "_peakLocation.pdf", sep = ""))
    barplot(c(numPromD, numPromI, numORF), names = c("before elements", "after elements", "in elements"),
            col = c("yellow", "orange", "green"), horiz = T,
            main = paste("Peak location regarding chromosomal features\n# of analyzed peaks:", numPeaks))    
    dev.off()

    print("Saving the results in:")
    print(paste(outputName, "_peakLocation.txt", sep = ""))
 
    # writing the results
    write.table(locRes, file = paste(outputName, "_peakLocation.txt", sep = ""), quote = F, 
                row.names = F, col.names = F)

    statsRes = list(numPeaks =  numPeaks, beforeFeatures = numPromD, 
                inFeatures = numORF, afterFeatures =  numPromI)   

   print("**********************************************")

   return(statsRes)


# end of function peakLocation()
}
