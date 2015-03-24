##########################
#time.disparity
##########################
#Calculates the disparity for interval pco.data and output a interval.disparity table object
#v0.2.3
##########################
#SYNTAX :
#<time_pco> time intervals or slices from a pco
#<...> disparity arguments (see ?disparity for information)
##########################
#----
#guillert(at)tcd.ie 16/03/2014
##########################

time.disparity<-function(time_pco, method=c("centroid", "sum.range", "product.range", "sum.variance", "product.variance"), CI=c(50, 95), bootstraps=1000, central_tendency=median, rarefaction=FALSE, verbose=FALSE, rm.last.axis=FALSE, save.all=FALSE) {
    #SANITIZING
    #time_pco
    check.class(time_pco, "list", " must be a list of time sections of pco data.")
    if(length(names(time_pco))!=length(time_pco)) {
        stop("time_pco data must have time sections names.")
    }

    #rarefaction
    if(rarefaction == TRUE) {
        message("Rarefaction is calculated and slows down the disparity calculation.\nUse Rarefaction=FALSE to speed up the calculations.")
    }

    #Managing bins with only one data point
    while(any(unlist(lapply(time_pco, nrow)) < 2)) {
        wrong_intervals<-which(unlist(lapply(time_pco, nrow)) < 2)
        #Moving the first wrong interval to the next interval in time (unless the wrong interval is the last one)
        if(wrong_intervals[1] != length(time_pco)) {
            host_interval<-wrong_intervals[1]+1
            message("Intervals ", names(time_pco)[wrong_intervals[1]], " and ", names(time_pco)[host_interval], " are combined due to insufficient data.")
        } else {
            #Moving the wrong interval in the preceding one
            host_interval<-wrong_intervals[1]-1
            message("Intervals ", names(time_pco)[host_interval], " and ", names(time_pco)[wrong_intervals[1]], " are combined due to insufficient data.")
        }
        #Creating the new interval
        new_interval<-rbind(time_pco[[wrong_intervals[1]]], time_pco[[host_interval]])
        #Creating the new time_pco data
        new_time_pco<-time_pco ; names(new_time_pco)<-names(time_pco)
        #replacing the wrong interval
        new_time_pco[[host_interval]]<-new_interval
        #renaming the interval
        if(wrong_intervals[1] != length(time_pco)) {
            names(new_time_pco)[host_interval]<-paste(strsplit(names(new_time_pco)[wrong_intervals[1]], split="-")[[1]][1],strsplit(names(new_time_pco)[host_interval], split="-")[[1]][2],sep="-")
        } else {
            names(new_time_pco)[host_interval]<-paste(strsplit(names(new_time_pco)[host_interval], split="-")[[1]][1],strsplit(names(new_time_pco)[wrong_intervals[1]], split="-")[[1]][2],sep="-")
        }
        #removing empty interval
        new_time_pco[[wrong_intervals[1]]]<-NULL
        time_pco<-new_time_pco
    }

    #CALCULATING THE DISPARITY FOR EACH BIN
    disparity_interval<-lapply(time_pco, disparity, method=method, CI=CI, bootstraps=bootstraps, central_tendency=central_tendency, rarefaction=rarefaction, verbose=verbose, rm.last.axis=rm.last.axis, save.all=save.all)

    #Return the table only
    if(save.all == FALSE) {
        #Sorting the data as a table
        if(rarefaction == FALSE) {
            #create the table's first row
            disparity_intervals_table<-disparity_interval[[1]]
            #Loop through the other elements of the table
            for(interval in 2:length(disparity_interval)) {
                disparity_intervals_table<-rbind(disparity_intervals_table, disparity_interval[[interval]])
            }
            #Renaming the rarefaction column interval
            colnames(disparity_intervals_table)[1]<-"time"
            #Saving the interval names
            disparity_intervals_table[,1]<-names(time_pco)

        } else {

            #If rarefaction has been calculated, only get the last element of each rarefaction table

            #create an interval row
            interval_row<-matrix(nrow=(nrow(disparity_interval[[1]])), data=rep(names(time_pco)[[1]]))
            #add the disparity results (with rarefaction)
            interval_tab<-cbind(interval_row, disparity_interval[[1]])
            #binding the interval table
            disparity_intervals_table<-rbind(interval_tab)

            #Loop through the other intervals
            for(interval in 2:length(time_pco)) {
                #create an interval row
                interval_row<-matrix(nrow=(nrow(disparity_interval[[interval]])), data=rep(names(time_pco)[[interval]]))
                #add the disparity results (with rarefaction)
                interval_tab<-cbind(interval_row, disparity_interval[[interval]])
                #binding the interval table
                disparity_intervals_table<-rbind(disparity_intervals_table, interval_tab)
            }

            #Renaming the rarefaction column interval
            colnames(disparity_intervals_table)[1]<-"time"
        }

        return(disparity_intervals_table)
    
    } else {

        #Sorting the data as a table
        if(rarefaction == FALSE) {
            #Creating the quantile table
            #create the table's first row
            disparity_intervals_table<-disparity_interval[[1]][[1]]
            #Loop through the other elements of the table
            for(interval in 2:length(disparity_interval)) {
                disparity_intervals_table<-rbind(disparity_intervals_table, disparity_interval[[interval]][[1]])
            }
            #Renaming the rarefaction column interval
            colnames(disparity_intervals_table)[1]<-"time"
            #Saving the interval names
            disparity_intervals_table[,1]<-names(time_pco)

        } else {

            #If rarefaction has been calculated, only get the last element of each rarefaction table

            #create an interval row
            interval_row<-matrix(nrow=(nrow(disparity_interval[[1]][[1]])), data=rep(names(time_pco)[[1]]))
            #add the disparity results (with rarefaction)
            interval_tab<-cbind(interval_row, disparity_interval[[1]][[1]])
            #binding the interval table
            disparity_intervals_table<-rbind(interval_tab)

            #Loop through the other intervals
            for(interval in 2:length(time_pco)) {
                #create an interval row
                interval_row<-matrix(nrow=(nrow(disparity_interval[[interval]][[1]])), data=rep(names(time_pco)[[interval]]))
                #add the disparity results (with rarefaction)
                interval_tab<-cbind(interval_row, disparity_interval[[interval]][[1]])
                #binding the interval table
                disparity_intervals_table<-rbind(disparity_intervals_table, interval_tab)
            }

            #Renaming the rarefaction column interval
            colnames(disparity_intervals_table)[1]<-"time"
        }

        #saving the results per time section
        disparity_intervals_values<-list()
        #Loop through all the elements of the list to extract the values
        for(interval in 1:length(disparity_interval)) {
            disparity_intervals_values[[interval]]<-disparity_interval[[interval]][[2]]
        }
        #Renaming the list elements
        names(disparity_intervals_values)<-names(time_pco)

        #output
        output<-list("quantiles"=disparity_intervals_table, "values"=disparity_intervals_values)
        return(output)
    }
}
