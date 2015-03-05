##########################
#bin.disparity
##########################
#Calculates the disparity for binned pco.data and output a bin.disparity table object
#v0.1
##########################
#SYNTAX :
#<pco_binned> binned pco data
#<...> disparity arguments (see ?disparity for information)
##########################
#----
#guillert(at)tcd.ie 06/03/2014
##########################

bin.disparity<-function(pco_binned, method=c("centroid", "sum.range", "product.range", "sum.variance", "product.variance"), CI=c(50, 95), bootstraps=1000, central_tendency=median, rarefaction=FALSE, verbose=FALSE) {
    #SANITIZING
    #pco_binned
    check.class(pco_binned, "list", " must be a list of binned pco data.")
    if(length(names(pco_binned))!=length(pco_binned)) {
        stop("pco_binned data must have bins names.")
    }

    #rarefaction
    if(rarefaction == TRUE) {
        message("Rarefaction is calculated and slows down the disparity calculation.\nUse Rarefaction=FALSE to speed up the calculations.")
    }

    #CALCULATING THE DISPARITY FOR EACH BIN
    disparity_binned<-lapply(pco_binned, disparity, method=method, CI=CI, bootstraps=bootstraps, central_tendency=central_tendency, rarefaction=rarefaction, verbose=verbose)

    #Sorting the data as a table
    if(rarefaction == FALSE) {
        #create the table's first row
        disparity_binned_table<-disparity_binned[[1]]
        #Loop through the other elements of the table
        for(bin in 2:length(disparity_binned)) {
            disparity_binned_table<-rbind(disparity_binned_table, disparity_binned[[bin]])
        }
        #Renaming the rarefaction column bin
        colnames(disparity_binned_table)[1]<-"bin"
        #Saving the bin names
        disparity_binned_table[,1]<-names(pco_binned)
    } else {
        #If rarefaction has been calculated, only get the last element of each rarefaction table
        #create the table's first row
        disparity_binned_table<-disparity_binned[[1]][nrow(disparity_binned[[1]]),]
        #Loop through the other elements of the table
        for(bin in 2:length(disparity_binned)) {
            disparity_binned_table<-rbind(disparity_binned_table, disparity_binned[[bin]][nrow(disparity_binned[[bin]]),])
        }
        #Renaming the rarefaction column bin
        colnames(disparity_binned_table)[1]<-"bin"
        #Saving the bin names
        disparity_binned_table[,1]<-names(pco_binned)
    }
    return(disparity_binned_table)
}
