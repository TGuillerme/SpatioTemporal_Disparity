##########################
#Ancestral state matrix
##########################
#Recreate the ancestral matrix for each node and each character
#v1.0
#Update: allow to use ML Bayesian or Threshold method
#Update: added a saving option
#Update: treats multistates characters
#Update: methods revisited
##########################
#SYNTAX :
#<tree> a 'phylo' object
#<nexus> a nexus file list containing the matrix and a list of ordering (use ReadMorphNexus{Claddis} for proper format)
#<method> the method to use for ancestral state reconstruction ('ML-ape' or 'ML-phytools')
#<verbose> whether to be verbose or not
###############################################<...> any optional arguments to be passed to anc.Bayes{phytools} or threshBayes{phytools} functions.
##########################
#----
#guillert(at)tcd.ie 06/03/2015
##########################


anc.state<-function(tree, nexus, method='ML', verbose=TRUE, ...){

#SANITYZING

    #tree
    check.class(tree, 'phylo', ' must be a phylo object.')
    #Is binary?
    tree<-bin.tree(tree)

    #nexus
    check.class(nexus, 'list', ' must be a nexus list.\n Use Claddis::ReadMorphNexus() for generating the proper formatted object.')
    #matrix element present?
    #$matrix
    if(!any(names(nexus) == "matrix")) {
        stop('nexus must be a nexus list.\n Use Claddis::ReadMorphNexus() for generating the proper formatted object.')
    }
    #$ordering
    if(!any(names(nexus) == "ordering")) {
        message('There was no character ordering list available in the nexus object:\n characters are now all considered as unordered.\n Use Claddis::ReadMorphNexus() for generating the proper formatted object.')
        #Generate default ordering (none)
        nexus$ordering<-c(rep("unord", ncol(nexus$matrix)))
    }

    #method
    check.class(method, 'character', ' must be \'ML-ape\' or \'ML-Claddis\'.')
    if(method !='ML-ape') {
        if(method != 'ML-Claddis') {
            stop('Method must be must be \'ML-ape\' or \'ML-Claddis\'.')
        }
    }

    #nexus (again)
    #$max and min values (if method = 'ML-Claddis')
    if(method == 'ML-Claddis' & !any(names(nexus) == "max.vals")) {
        stop('Nexus object needs to contain a \'max.vals\' vector if chosen method is \'ML-Claddis\'.\n Use Claddis::ReadMorphNexus() for generating the proper formatted object.')
    }
    if(method == 'ML-Claddis' & !any(names(nexus) == "min.vals")) {
        stop('Nexus object needs to contain a \'min.vals\' vector if chosen method is \'ML-Claddis\'.\n Use Claddis::ReadMorphNexus() for generating the proper formatted object.')
    }

    #verbose
    check.class(verbose, 'logical', ' must be logical.')

#ESTIMATING THE ANCESTRAL MATRIX FOR EACH CHARATER AND NODE

    #Ancestral states estimations from a matrix
    anc.list<-anc.state_ace(tree, nexus, method, verbose, ...)

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
