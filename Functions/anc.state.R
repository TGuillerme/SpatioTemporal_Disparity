##########################
#Ancestral state matrix
##########################
#Recreate the ancestral matrix for each node and each character
#v0.1
#To do: allow model to be not only ML
##########################
#SYNTAX :
#<tree> a 'phylo' object
#<matrix> the character matrix
#<model> the model to use for ancestral state reconstruction
##########################
#----
#guillert(at)tcd.ie 25/09/2014
##########################


anc.state<-function(tree, matrix, model='ML', verbose=TRUE){

#SANITYZING

    #tree
    check.class(tree, 'phylo', ' must be a phylo object.')
    #Is binary?
    tree<-bin.tree(tree)

    #matrix
    if(class(matrix) == 'data.frame') {
        matrix<-as.matrix(matrix)
    }
    check.class(matrix, 'matrix', ' must be a matrix.')

    #model
    check.class(model, 'character', ' must be \'ML\'.')
    if(model !='ML') {
        stop('type must be \'ML\'.')
    }

    #verbose
    check.class(verbose, 'logical', ' must be logical.')

#ESTIMATING THE ANCESTRAL MATRIX FOR EACH CHARATER AND NODE

    #Ancestral states estimations from a matrix
    anc.list<-anc.state_ace(tree, matrix, model, verbose)

    #Creating the state probability matrix for the nodes and the tips
    anc.prob<-anc.state_prob(tree, matrix, anc.list)

    #Creating the state matrix for the nodes and the tips
    anc.state<-anc.state_state(tree, matrix, anc.list)

#OUTPUT

    anc.matrix<-list("state"=anc.state, "prob"=anc.prob)
    return(anc.matrix)

}
