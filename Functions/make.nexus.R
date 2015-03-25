###########################
#Utility function for making a dummy nexus format matrix (empty)
##########################
#Generates a list following Claddis::ReadMorphNexus format
#v0.1
##########################
#SYNTAX :
#<matrix> a character matrix
#<header> optional. a header name
#<ordering> optional. a vector of "unord" / "ord" of the length of the matrix columns
#<weights> optional. a numeric vector of the length of the matrix columns
##########################
#----
#guillert(at)tcd.ie 25/03/2015
##########################

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