bPeaksAnalysis <-
function(IPdata, controlData, chromosomalFeatures = NULL, 
                  smoothingValue = 20, windowSize = 150, windowOverlap = 50, 
                  IPcoeff = 6, controlCoeff = 4, log2FC = 2, averageQuantiles = 0.9,
                  resultName = "bPeaks", peakDrawing = TRUE, promSize = 800){


#######
# written by G. LELANDAIS <gaelle.lelandais@univ-paris-diderot.fr>    
#######


    # argument names were changed for user convenience 
    INPUTdata = controlData

    # combine IP and INPUT data
    allData = cbind(IPdata, INPUTdata)
    colnames(allData) = c("chr", "pos", "IPsignal", "chr", "pos", "controlSignal")

    print("**********************************************")
    print("Starting analysis of genome-wide read depth...")
    # calculate baseline IP and INPUT values
    lineIP    = baseLineCalc(allData$IPsignal)
    print(paste("Baseline IP signal (average value): ", round(lineIP, 3)))
    lineINPUT = baseLineCalc(allData$controlSignal)
    print(paste("Baseline ontrol signal (average value): ", round(lineINPUT, 3)))
    print("**********************************************")
    print("")

    # get list of chromosomes
    chromNames = unique(allData[,1])

    print(paste(length(chromNames), "different chromosomes will be analyzed successively"))
    print("")
    
    # to store results
    peakStats = NULL
    #load("peakStats.Robject")

    for(smoo in smoothingValue){

    for(w in windowSize){
    
    for(t in windowOverlap){

    for(j in IPcoeff){

    for(k in controlCoeff){

    for(r in log2FC){

    for(s in averageQuantiles){

    allPeaks = NULL

    print("**********************************************")
    print("Parameter values for basic peaks (bPeaks) detection :")
    print("")
    print(paste("Smoothing:", smoo))
    print(paste("Window size:", w))
    print(paste("Window overlap:", t))
    print(paste("IP threshold (T1): > x", j, sep = ""))
    print(paste("Control threshold (T2): < x", k, sep = ""))
    print(paste("Log2FC threshlod (T3): >", r, sep = ""))
    print(paste("Average log2 signals (quantiles, between 0 and 1) (T4): > ", s, sep = ""))
    print("**********************************************")
    print("")
    print("")
    print("")

    for(i in 1:length(chromNames)){

        print("**********************************************")
        print(paste("Starting analysis of chromosome ", chromNames[i]))

        # information for one chromosome
        subData = allData[allData[,1] == chromNames[i],]

        # test of the smoothing function
        vecIP = subData[,3]
        vecINPUT = subData[,6]

        # smooth of the data
        print("Starting data smoothing...")
        smoothedIP    = dataSmoothing(vecData = vecIP, widthValue = smoo) 
        print("...done (IP signal)")
        smoothedINPUT = dataSmoothing(vecData = vecINPUT, widthValue = smoo) 
        print("...done (control signal)")

        # seed detection
        res = peakDetection(IPdata = smoothedIP, controlData = smoothedINPUT, 
                           chrName = as.character(chromNames[i]), 
                            windowSize = w, windowOverlap = t,
                            outputName = paste(resultName, "-", 
                                            as.character(chromNames[i]), sep = ""), 
                            baseLineIP = lineIP, baseLineControl = lineINPUT,
                            IPthreshold = j, controlThreshold = k,
                            ratioThreshold = r, averageThreshold = s,
                            peakDrawing = peakDrawing)

    # add chromosom information
    if(is.null(res) == F){
        res2 = cbind(as.character(chromNames[i]), res)
        allPeaks = rbind(allPeaks, res2)
    }

    print("**********************************************")
    print("")


    }

    # bug correction in case of no peak is detected
    if(is.null(allPeaks) == F){

        print("Saving BED file with all bPeak information:")
        print(paste(resultName, "_bPeaks_allGenome.bed"))
        # write final BED final 
        allPeaks = allPeaks[order(as.numeric(allPeaks[,colnames(allPeaks) == "IP"]), decreasing = T),]
        bedFinalData = cbind(allPeaks[,1:3], paste("allGenome_bPeak_", 1:nrow(allPeaks), sep = ""), 
                                allPeaks[,4:6])

        write.table(bedFinalData, file = paste(resultName, "_bPeaks_allGenome.bed", sep = ""), 
                quote = F, sep = "\t", eol = "\n", 
                col.names = F, row.names = F)
        print("**********************************************")
        print("")
 
        if(is.null(chromosomalFeatures) == F){

           ORFinfo = chromosomalFeatures
           bedFile = paste(resultName, "_bPeaks_allGenome.bed", sep = "")

           resPeakLoc = peakLocation(bedFile = bedFile, genomicInfo = ORFinfo, outputName = resultName, promSize = promSize)

            NumIn = resPeakLoc$inFeatures
            NumBefore = resPeakLoc$beforeFeatures
            NumAfter  = resPeakLoc$afterFeatures

         }else{

              NumIn = NA
              NumBefore = NA
              NumAfter  = NA
          }

      # statistics related to the detected peaks
      peakStats = rbind(peakStats, c(smoo, w, t, j, k, r, s, nrow(allPeaks), round(mean(as.numeric(allPeaks[,3]) - as.numeric(allPeaks[,2])),3), round(mean(as.numeric(allPeaks[,4])),3), round(mean(as.numeric(allPeaks[,5])), 3), round(mean(as.numeric(allPeaks[,6])), 3), NumBefore, NumAfter, NumIn))

      save(peakStats, file = "peakStats.Robject")

    # end of if(is.null(allPeaks) == F)
    }

    # end of x
    }

    # end of s
    }

    # end of r
    }

    # end of k
    }

    # end of j
    }

    # end of t
    }

    # end of w
    }

    # at the end only
    if(is.null(nrow(peakStats) == F)){

    colnames(peakStats) = c("smoothingValue", "windowSize", "windowOverlap", "IPcoeff", "controlCoeff", "log2FC", "averageQuantiles", "bPeakNumber", "meanSize", "meanIPsignal", "meanControlSignal", "meanLog2FC", "bPeakNumber_beforeFeatures", "bPeakNumber_afterFeatures", "bPeakNumber_inFeatures")

    write.table(peakStats, file = paste(resultName, "_bPeaks_parameterSummary.txt", sep = ""), quote = F,
                sep = "\t", row.names = F, col.names = T)

    save(peakStats, file = "peakStats.Robject")

    }

    print("")
    print("**********************************************************************")
    print("END OF THE PROCEDURE")
    print("Thanks for using bPeaks library. Feedbacks are welcomed :)")
    print("**********************************************************************")

# end of function bPeaks()
}
