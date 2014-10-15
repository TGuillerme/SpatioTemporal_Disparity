##########################
#PCO for STD
##########################
#Pipeline for running a pco on the ancestral matrix from anc.state
#v0.1.1
#Update:
#This version is methodologically wrong:
#If NA's are introduced in the distance matrix they are replaced by the mean distance.
##########################
#SYNTAX :
#<anc.matrix> the ancestral matrix (a list of states and probabilities returned from anc.state).
#<distance> the method for calculating the distance in vegdist(){vegan}.
#<scale> whether to scale the data (default=FALSE).
#<center> whether to center the data (default=FALSE).
#<na.rm> Pairwise deletion of missing observations when computing dissimilarities (default = FALSE).
#<correction> Correction methods for negative eigenvalues (details in pcoa(){ape}) (default = "none").
#<...> optional arguments to be passed to vegdist(){vegan}.
##########################
#----
#guillert(at)tcd.ie 29/09/2014
##########################

pco.std<-function(anc.matrix, distance="euclidean", scale=FALSE, center=FALSE, na.rm=FALSE, correction="none", ...){

    #SANITIZING
    #packages
    require(vegan)
    require(ape)

    #DISCLAIMER
    warning("This function is in development.\nIf NA's are introduced in the distance matrix they are replaced by the mean distance.\nThis is methodological wrong!")

    #anc.matrix
    check.class(anc.matrix, "list", " must be a list from anc.state containing two elements: \'state\' and \'prob\'.")
    check.length(anc.matrix, 2, " must be a list from anc.state containing two elements: \'state\' and \'prob\'.")
    if(names(anc.matrix)[1] != "state") {
        stop(as.character(substitute(anc.matrix)), " must be a list from anc.state containing two elements: \'state\' and \'prob\'.", call.=FALSE)
    } 
    if(names(anc.matrix)[2] != "prob") {
        stop(as.character(substitute(anc.matrix)), " must be a list from anc.state containing two elements: \'state\' and \'prob\'.", call.=FALSE)
    }
    matrix<-anc.matrix$state

    #distance
    METHODS<-c("manhattan", "euclidean", "canberra", "bray", "kulczynski", "gower", "morisita", "horn", "mountford", "jaccard", "raup", "binomial", "chao", "altGower", "cao")
    if(any(distance == METHODS)) {
        method<-distance
    } else {
        stop("Wrong distance method: \n see ?vegdist for available distances.", call.=FALSE)
    }

    #scale
    check.class(scale, 'logical', " must be logical.")

    #center
    check.class(center, 'logical', " must be logical.")

    #na.rm
    check.class(na.rm, 'logical', " must be logical.")

    #correction
    CORRECTIONS <- c("none", "lingoes", "cailliez")
    if(any(correction == CORRECTIONS)) {
        correction<-correction
    } else {
        stop("Wrong correction method: \n see ?pcoa for available corrections.", call.=FALSE)
    }

    #PCO

    #Transforming "?" in NAs in the matrix
    for (character in 1:ncol(matrix)) {
        for (taxa in 1:nrow(matrix)) {
            if(as.character(matrix[taxa, character]) == "?") {
                matrix[taxa, character] <- NA
            }
        }
    }

    #Transforming the matrix to numeric
    matrix<-as.data.frame(matrix)
    for (character in 1:ncol(matrix)) {
        matrix[,character]<-as.numeric(matrix[,character])
    }

    #Scaling
    scaled.matrix<-scale(matrix, scale, center)

    #Calculating the distance matrix
    if(na.rm == TRUE) {
        distance.matrix<-vegdist(scaled.matrix, method, na.rm=TRUE, ...)
    } else {
        distance.matrix<-vegdist(scaled.matrix, method, na.rm=FALSE, ...)
    }


    #DEVELOPEMENT PART (FIX!)
    if(length(which(is.na(distance.matrix))) > 0) {
        nas<-length(which(is.na(distance.matrix)))
        distance.matrix[which(is.na(distance.matrix))]<-mean(distance.matrix, na.rm=TRUE)
        warning("Version in development:\n", nas, " NAs have been replaced by the mean distance.\nThis is methodologically wrong!")
    }





    #Calculating the pco
    pco.std<-pcoa(distance.matrix, correction)

    #OUTPUT
    return(pco.std)
#End
}