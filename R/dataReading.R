dataReading <-
function(IPfile, controlFile, chromosomalFeatures = ""){

#######
# written by G. LELANDAIS <gaelle.lelandais@univ-paris-diderot.fr>    
#######

    print("**********************************************")
    print("Reading IP and control datasets... ")

    IPdata = read.table(IPfile)
    controlData = read.table(controlFile)

    print("... done")

    if(chromosomalFeatures != ""){
    
        print("Reading chromosomal features for peak location... ")

        ORFinfo = as.matrix(read.table(chromosomalFeatures, sep = "\t"))

        print("... done")
    }else{
        ORFinfo = NULL
    }
    
    print("**********************************************")
    print("")

    return(list(IPdata = IPdata, controlData = controlData, chromosomalFeatures = ORFinfo))

# end of function dataReading()
}
