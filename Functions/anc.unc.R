##########################
#Ancestral uncertainty
##########################
#Transforms the ancestral state matrix by excluding all infered states below a certain threshold.
#v0.1
##########################
#SYNTAX :
#<anc.matrix> the ancestral matrix (a list of states and probabilities returned from anc.state)
#<threshold> the probability threshold (a value between 0.5 and 1)
##########################
#----
#guillert(at)tcd.ie 29/09/2014
##########################

anc.unc<-function(anc.matrix, threshold=0.5) {
    #SANITIZING
    #anc.matrix
    check.class(anc.matrix, 'list', ' must the ancestral matrix with a list of states and a list of probabilities.')
    check.length(anc.matrix, 2, ' must the ancestral matrix with a list of states and a list of probabilities.')
    if(names(anc.matrix)[1] != 'state') {
        stop('anc.matrix must the ancestral matrix with a list of states and a list of probabilities.', call.=FALSE)
    }
    if(names(anc.matrix)[2] != 'prob') {
        stop('anc.matrix must the ancestral matrix with a list of states and a list of probabilities.', call.=FALSE)
    }

    #threshold
    check.class(threshold, 'numeric', ' must be a numerical value between 0.5 and 1.')
    check.length(threshold, 1, ' must be a numerical value between 0.5 and 1.')
    if(threshold < 0.5) {
        stop('threshold must be a numerical value between 0.5 and 1.', call.=FALSE) 
    }
    if(threshold > 1) {
        stop('threshold must be a numerical value between 0.5 and 1.', call.=FALSE) 
    }

    #REPLACING CHARACTER STATES BY "?" IN THE MATRIX IF STATE PROBABILITY < threshold
    for (character in 1:ncol(anc.matrix$prob)) {
        for (taxa in 1:nrow(anc.matrix$prob)) {
            if(anc.matrix$prob[taxa, character] < threshold) {
                anc.matrix$state[taxa, character] <- "?"
            }
        }
    }

    #Return
    return(anc.matrix)
}