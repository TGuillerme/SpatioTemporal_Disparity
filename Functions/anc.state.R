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


anc.state<-function(tree, nexus, method='ML-ape', verbose=TRUE, ...){

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
    check.class(method, 'character', ' must be \'ML-ape\' or \'ML-claddis\'.')
    if(method !='ML-ape') {
        if(method != 'ML-claddis') {
            stop('Method must be must be \'ML-ape\' or \'ML-claddis\'.')
        }
    }

    #nexus (again)
    #$max and min values (if method = 'ML-claddis')
    if(method == 'ML-claddis' & !any(names(nexus) == "max.vals")) {
        stop('Nexus object needs to contain a \'max.vals\' vector if chosen method is \'ML-claddis\'.\n Use Claddis::ReadMorphNexus() for generating the proper formatted object.')
    }
    if(method == 'ML-claddis' & !any(names(nexus) == "min.vals")) {
        stop('Nexus object needs to contain a \'min.vals\' vector if chosen method is \'ML-claddis\'.\n Use Claddis::ReadMorphNexus() for generating the proper formatted object.')
    }

    #verbose
    check.class(verbose, 'logical', ' must be logical.')

#ESTIMATING THE ANCESTRAL MATRIX FOR EACH CHARATER AND NODE

    #Ancestral states estimations from a matrix
    anc.list<-anc.state_ace(tree, nexus, method, verbose, ...)

    #Creating the state probability matrix for the nodes and the tips
    anc.prob<-anc.state_prob(tree, nexus$matrix, anc.list)

    #Creating the state matrix for the nodes and the tips
    anc.state<-anc.state_state(tree, nexus$matrix, anc.list)

    #Creating the rate matrix
    if(method == 'ML-ape') {
        anc.rate<-anc.state_rate(tree, nexus$matrix, anc.list)
    }

#OUTPUT

    if(method == 'ML-ape') {
        anc.matrix<-list("state"=anc.state, "prob"=anc.prob, "rate"=anc.rate)
    } else {
        anc.matrix<-list("state"=anc.state, "prob"=anc.prob)
    }
    return(anc.matrix)

}
