##########################
#Ancestral state matrix
##########################
#Recreate the ancestral matrix for each node and each character
#v0.2
#Update: allow to use ML Bayesian or Threshold method
#Update: added a saving option
##########################
#SYNTAX :
#<tree> a 'phylo' object
#<matrix> the character matrix
#<method> the method to use for ancestral state reconstruction ('ML' or 'Bayesian' or 'Threshold')
#<verbose> whether to be verbose or not
#<save> whether to save the matrix or not. If TRUE, the default saving name is "anc.state.save", a name can be provided as save="name".
#<...> any optional arguments to be passed to anc.Bayes{phytools} or threshBayes{phytools} functions.
##########################
#----
#guillert(at)tcd.ie 29/09/2014
##########################


anc.state<-function(tree, matrix, method='ML', verbose=TRUE, save=FALSE, ...){

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

    #method
    check.class(method, 'character', ' must be \'ML\', \'Bayesian\' or \'Threshold\'.')
    if(method !='ML') {
        if(method != 'Bayesian') {
            if(method != 'Threshold') {
                stop('Method must be \'ML\', \'Bayesian\' or \'Threshold\'.')
            }
        }
    }

    #verbose
    check.class(verbose, 'logical', ' must be logical.')

    #save
    if(class(save) == 'logical') {
        if(save == TRUE) {
            save.name="anc.state.save"
        }
    } else {
        check.class(save, 'character', ' must be \'logical\' or \'character\'')
        save.name<-save
        save<-TRUE
    }

#ESTIMATING THE ANCESTRAL MATRIX FOR EACH CHARATER AND NODE

    #Ancestral states estimations from a matrix
    anc.list<-anc.state_ace(tree, matrix, method, verbose, ...)

    #Creating the state probability matrix for the nodes and the tips
    anc.prob<-anc.state_prob(tree, matrix, anc.list)

    #Creating the state matrix for the nodes and the tips
    anc.state<-anc.state_state(tree, matrix, anc.list)

    #Creating the rate matrix
    anc.rate<-anc.state_rate(tree, matrix, anc.list)

#OUTPUT

    anc.matrix<-list("state"=anc.state, "prob"=anc.prob, "rate"=anc.rate)
    if (save == TRUE) {
        if(method == 'ML') {
            save(anc.matrix, files=paste(save.name, ".rda", sep=""))
        } else {
            save(anc.matrix, files=paste(save.name, ".rda", sep=""))
            save(anc.matrix.trace, files=paste(save.name, ".trace", sep=""))
        }
    }
    return(anc.matrix)

}
