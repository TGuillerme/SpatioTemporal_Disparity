##########################
#Correct a diversity vector
##########################
#Corrects a diversity vector the same way as in time.disparity does (when diversity is equal to one in an interval, this interval is combined with the following one.)
#v0.1
##########################
#SYNTAX :
#<diversity> a vector of taxa count as produced by the disparity function
##########################
#----
#guillert(at)tcd.ie 23/03/2015
##########################

cor.diversity<-function(diversity) {
    #SANITIZING
    check.class(diversity, "integer", " must be a vector of taxa counts.")
    if(is.null(names(diversity))) {
        stop("diversity must be a vector of taxa counts with intervals names.")
    }
    
    #CORRECTING THE DIVERSITY VECTOR
    while(any(diversity < 2)) {
        
        #Selecting the wrong interval
        wrong_intervals<-which(diversity < 2)
        names(wrong_intervals)<-NULL
        
        #Moving the first wrong interval to the next interval in time (unless the wrong interval is the last one)
        if(wrong_intervals[1] != length(diversity)) {
            host_interval<-wrong_intervals[1]+1
            names(host_interval)<-NULL
            #message("Intervals ", names(diversity)[wrong_intervals[1]], " and ", names(diversity)[host_interval], " are combined due to insufficient data.")
        } else {
            #Moving the wrong interval in the preceding one
            host_interval<-wrong_intervals[1]-1
            names(host_interval)<-NULL
            #message("Intervals ", names(diversity)[host_interval], " and ", names(diversity)[wrong_intervals[1]], " are combined due to insufficient data.")
        }

        #Creating the new interval
        new_interval<-diversity[wrong_intervals[1]]+diversity[host_interval] ; names(new_interval)<-NULL
        #Creating the new diversity data
        new_diversity<-diversity ; names(new_diversity)<-names(diversity)
        #replacing the wrong interval
        new_diversity[host_interval]<-new_interval
        #renaming the interval
        if(wrong_intervals[1] != length(diversity)) {
            names(new_diversity)[host_interval]<-paste(strsplit(names(new_diversity)[wrong_intervals[1]], split="-")[[1]][1],strsplit(names(new_diversity)[host_interval], split="-")[[1]][2],sep="-")
        } else {
            names(new_diversity)[host_interval]<-paste(strsplit(names(new_diversity)[host_interval], split="-")[[1]][1],strsplit(names(new_diversity)[wrong_intervals[1]], split="-")[[1]][2],sep="-")
        }
        #removing empty interval
        new_diversity<-new_diversity[-wrong_intervals[1]]
        diversity<-new_diversity
    }

    return(diversity)
}