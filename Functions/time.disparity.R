##########################
#time.disparity
##########################
#Calculates the disparity for binned pco.data and output a bin.disparity table object
#v0.2
##########################
#SYNTAX :
#<time_pco> time intervals or slices from a pco
#<...> disparity arguments (see ?disparity for information)
##########################
#----
#guillert(at)tcd.ie 16/03/2014
##########################

time.disparity<-function(time_pco, method=c("centroid", "sum.range", "product.range", "sum.variance", "product.variance"), CI=c(50, 95), bootstraps=1000, central_tendency=median, rarefaction=FALSE, verbose=FALSE, rm.last.axis=FALSE) {
    #SANITIZING
    #time_pco
    check.class(time_pco, "list", " must be a list of binned pco data.")
    if(length(names(time_pco))!=length(time_pco)) {
        stop("time_pco data must have bins names.")
    }

    #rarefaction
    if(rarefaction == TRUE) {
        message("Rarefaction is calculated and slows down the disparity calculation.\nUse Rarefaction=FALSE to speed up the calculations.")
    }

    #CALCULATING THE DISPARITY FOR EACH BIN
    disparity_interval<-lapply(time_pco, disparity, method=method, CI=CI, bootstraps=bootstraps, central_tendency=central_tendency, rarefaction=rarefaction, verbose=verbose, rm.last.axis=rm.last.axis)

    #Sorting the data as a table
    if(rarefaction == FALSE) {
        #create the table's first row
        disparity_intervals_table<-disparity_interval[[1]]
        #Loop through the other elements of the table
        for(bin in 2:length(disparity_interval)) {
            disparity_intervals_table<-rbind(disparity_intervals_table, disparity_interval[[bin]])
        }
        #Renaming the rarefaction column bin
        colnames(disparity_intervals_table)[1]<-"time section"
        #Saving the bin names
        disparity_intervals_table[,1]<-names(time_pco)
    } else {
        #If rarefaction has been calculated, only get the last element of each rarefaction table
        #create the table's first row
        disparity_intervals_table<-disparity_interval[[1]][nrow(disparity_interval[[1]]),]
        #Loop through the other elements of the table
        for(bin in 2:length(disparity_interval)) {
            disparity_intervals_table<-rbind(disparity_intervals_table, disparity_interval[[bin]][nrow(disparity_interval[[bin]]),])
        }
        #Renaming the rarefaction column bin
        colnames(disparity_intervals_table)[1]<-"time section"
        #Saving the bin names
        disparity_intervals_table[,1]<-names(time_pco)
    }
    return(disparity_intervals_table)
}
