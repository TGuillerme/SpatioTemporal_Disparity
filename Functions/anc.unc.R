##########################
#Ancestral uncertainty
##########################
#Transforms the ancestral state matrix by excluding all infered states below a certain threshold.
#v0.1.1
#Update:
#Fixed sanitizing to match with the new anc.matrix type (+rate)
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
    check.class(anc.matrix, "list", " must be a list from anc.state containing three elements: \'state\', \'prob\' and \'rate\'.")
    check.length(anc.matrix, 3, " must be a list from anc.state containing three elements: \'state\', \'prob\' and \'rate\'.")
    if(names(anc.matrix)[1] != "state") {
        stop(as.character(substitute(anc.matrix)), " must be a list from anc.state containing three elements: \'state\', \'prob\' and \'rate\'.", call.=FALSE)
    } 
    if(names(anc.matrix)[2] != "prob") {
        stop(as.character(substitute(anc.matrix)), " must be a list from anc.state containing three elements: \'state\', \'prob\' and \'rate\'.", call.=FALSE)
    }
    if(names(anc.matrix)[3] != "rate") {
        stop(as.character(substitute(anc.matrix)), " must be a list from anc.state containing three elements: \'state\', \'prob\' and \'rate\'.", call.=FALSE)
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