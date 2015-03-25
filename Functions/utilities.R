###########################
#Utility functions
##########################

##########################
#make.nexus
##########################
#Generates a list following Claddis::ReadMorphNexus format
#----
#SYNTAX :
#<matrix> a character matrix
#<header> optional. a header name
#<ordering> optional. a vector of "unord" / "ord" of the length of the matrix columns
#<weights> optional. a numeric vector of the length of the matrix columns
#----
make.nexus<-function(matrix, header, ordering, weights) {
    #SANITIZING
    #matrix
    check.class(matrix, "matrix")

    #header
    if(missing(header)) {
        header<-NA
    } else {
        check.class(header, "character")
        check.length(header, 1, " must be a single character string.", errorif=FALSE)
    }

    #ordering
    if(missing(ordering)) {
        ordering<-rep("unord", ncol(matrix))
    } else {
        check.class(ordering, "character")
        check.length(ordering, ncol(matrix), " must be the same length as the matrix.", errorif=FALSE)
        options(warn=-1)
        if(any(ordering != c("unord", "ord"))) {
            stop("Ordering vector must contain only 'unord' or/and 'ord' values.")
        }
        options(warn=0)
    }

    #weights
    if(missing(weights)) {
        weights<-rep(1, ncol(matrix))
    } else {
        check.class(weights, "integer")
        check.length(weights, ncol(matrix), " must be the same length as the matrix.", errorif=FALSE)
    }

    #BUILD THE NEXUS OBJECT
    nexus<-list()
    nexus$header<-header
    nexus$matrix<-matrix
    nexus$ordering<-ordering
    nexus$weights<-weights
    nexus$max.vals<-apply(matrix, 2, max, na.rm=TRUE)
    nexus$min.vals<-apply(matrix, 2, min, na.rm=TRUE)

    return(nexus)
}

##########################
#cor.diversity
##########################
#Correct a diversity vector
#----
#SYNTAX :
#<diversity> a vector of taxa count as produced by the disparity function
#----
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

##########################
#states.count
##########################
#Count the number of states per characters
#----
#SYNTAX :
#<character> a vector of character states
#----
states.count<-function(character) {
    #Isolate the states
    states<-levels(as.factor(character))
    #Check if multi states
    if(length(grep("&", states)) > 0) {
        #Isolating the multi states
        multi_states<-states[grep("&", states)]
        multi_states<-as.factor(unlist(strsplit(multi_states, split="&")))
        #Removing the multi states from the states list
        states<-states[-grep("&", states)]
        #Check if any of the multi states is not yet present in the states list
        if(any(is.na(match(levels(multi_states), states)))) {
            states<-c(states, multi_states[which(is.na(match(levels(multi_states), states)))])
        }
    }
    #Count the number of states
    return(length(states))
}

##########################
#extract.dist
##########################
#extract a distance from a distance matrix build using Claddis::MorphDistMatrix
#----
#SYNTAX :
#<dist.list> the distances list output from Claddis::MorphDistMatrix
#<distance> the selected distance
#----
extract.dist<-function(dist.list, distance) {
    #SANITIZING
    #dist.list
    check.class(dist.list, "list")

    #Extracting the distance element from the distance list
    if(distance == "numeric") {
        output<-dist.list[[distance]]
    } else {
        output<-dist.list[[grep(distance, names(dist.list))]]
    }
    return(output)
}