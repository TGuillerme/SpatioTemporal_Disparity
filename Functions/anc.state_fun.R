#FUNCTIONS FOR anc.state

#Ancestral states estimations from a matrix
anc.state_ace<-function(tree, nexus, method, verbose, ...) {
    #Empty return list
    anc.list<-list()
    
    #If method is ML-ace and all characters are not unordered print a warning
    if(!all(nexus$ordering == "unord") & method == 'ML-ape') {
        message("Some ordered characters where found but the chosen method (\'ML-ape\') will count them as unordered.")
    }

    #If method is ML-claddis and some characters are multistate print a warning
    #if(grep("&", nexus$matrix) > 0 & method == 'ML-claddis') {
    #   message("Some multistate characters have been detected and chosen method is \'ML-claddis\'.\n
    #       The ancestral states of these characters will be estimated based on continuous ancestral character estimation (ape::ace)\n
    #       and the displayed likelihood score corresponds ")
    #}

    #Be verbose
    if(verbose == TRUE) {
        message('Estimating the ancestral states for ', ncol(nexus$matrix), ' characters:', appendLF=FALSE)
    }

    if(method == 'ML-ape') {
        for (character in 1:ncol(nexus$matrix)) {
            #ML
            anc.list[[character]]<-ace(as.factor(nexus$matrix[,character]), tree, type="d", model="ER", ...)
            #Be verbose
            if(verbose == TRUE) {
                message('.', appendLF=FALSE)
            }
        }
    }

    if(method == 'ML-claddis') {
        for (character in 1:ncol(nexus$matrix)) {
            #The following code is modified from claddis::AncStateEstMatrix.
            #Get minimum value for character:
            minval<-nexus$min.vals[character]
            #Get maximum value for character:
            maxval<-nexus$max.vals[character]

            #If all character is only NA
            if(is.na(minval) && is.na(maxval)) {
                anc.lik<-matrix(c(rep(1,Nnode(tree)), rep(0, Nnode(tree))), ncol=2)
                #column names
                colnames(anc.lik)<-c(NA, 0)
                #rownames
                if(length(grep(tree$node.label)) > 0) {
                    rownames(anc.lik)<-tree$node.label
                }
                #Save the result
                anc.list[[character]]<-list(anc.lik=anc.lik)

            } else {
                
                #Check if character is variable (non-constant)
                if(maxval == minval) {
                    anc.lik<-matrix(c(rep(1,Nnode(tree)), rep(0, Nnode(tree))), ncol=2)
                    #column names
                    colnames(anc.lik)<-rep(minval,2)
                    #rownames
                    if(length(grep(tree$node.label)) > 0) {
                        rownames(anc.lik)<-tree$node.label
                    }
                    #Save the result
                    anc.list[[character]]<-list(anc.lik=anc.lik)

                } else {

                    #Estimate the ancestral state based on Claddis method
                    #If estimating states for all taxa then treat missing values as all possible states:
                    nexus$matrix[grep(TRUE, is.na(nexus$matrix[, character])), character]<-paste(minval:maxval, collapse="&")
                    #Find tips which cannot be used due to missing data:
                    tipstogo<-rownames(nexus$matrix)[grep(TRUE, is.na(nexus$matrix[, character]))]
                        
                    #Only continue if at least three tips in pruned tree:
                    if(length(tipstogo) > (Ntip(tree) - 2)) {
                        stop("Input matrix contains not enough data.")
                    }
                    #Create the character tree
                    if(length(tipstogo) > 0) {
                        chartree<-drop.tip(tree, tipstogo)
                    } else {
                        chartree<-tree
                    }
                                        #Get tip values for the pruned tree:
                    tipvals<-nexus$matrix[chartree$tip.label, character]
                                
                    #Set discrete character estimation model if unordered and or binary character:
                    if(maxval - minval == 1 || maxval - minval > 1 && nexus$ordering[character] == "unord") {
                        mymodel<-"ER"
                    } else {
                        #Set discrete character estimation model if ordered multistate character:
                        #Create all zero matrix:
                        mymodel<-matrix(0, nrow=(maxval - minval) + 1, ncol=(maxval - minval) + 1)
                        #Name rows and columns as states:
                        rownames(mymodel)<-colnames(mymodel)<-minval:maxval    
                        #Enter one for all the off-diagonal diagonals (an ordered change model):
                        for(j in 1:(length(mymodel[1, ]) - 1)) mymodel[j + 1, j]<-mymodel[j, j + 1]<-1
                    }
                                
                    #Create matrix to store probabilities of tip values:
                    tipvals.mat<-matrix(0, nrow=length(tipvals), ncol=maxval - minval + 1)
                    #Add rownames (tip labels):
                    rownames(tipvals.mat)<-names(tipvals)
                    #Add colunames (state values):
                    colnames(tipvals.mat)<-minval:maxval
                    #Fill all probabilities equal to one (non-polymorphisms):
                    for(j in colnames(tipvals.mat)) {
                        tipvals.mat[grep(TRUE, tipvals == j), j]<-1
                    }

                    #Set polymorphisms
                    #If there are polymorphisms make all observed states equally probable:
                    if(any(apply(tipvals.mat, 1, sum) == 0)) {
                    #Get list of tip values with polymorphisms:
                        polymorphism.values<-grep(TRUE, apply(tipvals.mat, 1, sum) == 0)
                        #Go through each polymorphism:
                        for(j in polymorphism.values) {
                        #Get list of each state:
                        states<-strsplit(tipvals[j], "&")[[1]]
                            #Make each state equally probable in tip values matrix:
                            tipvals.mat[j, states]<-1 / length(states)
                        }
                    }

                    #Estimating the ancestral states.
                    #Remove any potential node labels on the character tree to avoid an error from rerootingMethod():
                    chartree$node.label<-NULL 
                    #Get likelihoods for each state in taxa and ancestors:
                    state_likelihoods<-rerootingMethod(chartree, tipvals.mat, model=mymodel, ...)$marginal.anc

                    #Selecting only the nodes
                    lik.anc<-state_likelihoods[-c(1:Ntip(chartree)),]
                    #adding the node names (if available)
                    if(length(tree$node.label) != 0) {
                        rownames(lik.anc)<-tree$node.label
                    }
                                  
                    #Set in the same format as ML-ape output
                    anc.list[[character]]<-list(lik.anc=lik.anc)
                    #Be verbose
                    if(verbose == TRUE) {
                        message('.', appendLF=FALSE)
                    }

                #End condition (if not constant)
                }        
            #End condition (if not NA)
            }
        #End loop (character)
        }
    #End condition (method == ML-Claddis)
    }

    if(method == 'Bayesian') {
        stop("Bayesian method in development.")
        #Bayesian
        require(phangorn)

        #ancestral.pml(type="bayes")


        characters<-as.factor(matrix[,character])
        names(characters)<-row.names(matrix)
        anc.list[[character]]<-anc.Bayes(tree, characters, ...)
        #Generate manageable output (transform it + allow saving the trace).
        if(verbose == TRUE) {
            message('.', appendLF=FALSE)
        }
        
        }
    if(method == 'Threshold') {
        stop("Threshol method in development.")
        #Revell's threshold method
        threshBayes()
    }

    if(method == 'MrBayes') {
        stop("MrBayes method in development.")
        #Creates a MrBayes script
    }
    

    if(verbose == TRUE) {
        message('Done.', appendLF=FALSE)
    }
    return(anc.list)

}

#Creating the state probability matrix for the nodes and the tips
anc.state_prob<-function(tree, matrix, anc.state_ace) {
    #Creating the empty matrix
    prob.matrix<-matrix(data=1, ncol=ncol(matrix), nrow=(nrow(matrix)+Nnode(tree)))
    if(!is.null(tree$node.label)) {
        row.names(prob.matrix)<-c(row.names(matrix), tree$node.label)
    } else {
        row.names(prob.matrix)<-c(row.names(matrix), paste("n",seq(1:Nnode(tree)), sep=""))
    }
    #Filling the matrix
    for (character in 1:ncol(matrix)) {
        for (edge in 1:Nnode(tree)) {
            prob.matrix[(nrow(matrix)+edge),character]<-max(anc.state_ace[[character]]$lik.anc[edge,])
        }
    }
    return(prob.matrix)
}

#Creating the state matrix for the nodes and the tips
anc.state_state<-function(tree, matrix, anc.state_ace) {
    #Creating the empty matrix
    state.matrix<-matrix(NA, ncol=ncol(matrix), nrow=(nrow(matrix)+Nnode(tree)))
    if(!is.null(tree$node.label)) {
        row.names(state.matrix)<-c(row.names(matrix), tree$node.label)
    } else {
        row.names(state.matrix)<-c(row.names(matrix), paste("n",seq(1:Nnode(tree)), sep=""))
    }
    #Filling the matrix
    #Tips states (observed)
    state.matrix[1:nrow(matrix), 1:ncol(matrix)]<-as.matrix(matrix)
    #Node states(estimated)
    for (character in 1:ncol(matrix)) {
        for (edge in 1:Nnode(tree)) {
            #Extracting the column name for each edge and character
            max.prob<-which(anc.state_ace[[character]]$lik.anc[edge,]==max(anc.state_ace[[character]]$lik.anc[edge,]))[[1]]
            state.matrix[(nrow(matrix)+edge),character]<-colnames(anc.state_ace[[character]]$lik.anc)[max.prob] #-1 is to make the character start at 0 (which greps 1,2,3, etc... instead of 0,1,2, ...)
        }
    }
    return(state.matrix)
}


#Creating the rate matrix for the nodes and the tips
anc.state_rate<-function(tree, matrix, anc.state_ace) {
    #Creating the empty matrix
    rate.matrix<-as.data.frame(matrix(NA, ncol=2, nrow=(ncol(matrix))))
    names(rate.matrix)<-c("rate.estimate", "std.err")
    #Adding the character names if available:
    if(!is.null(colnames(matrix))) {
        rownames(rate.matrix)<-colnames(matrix)
    }

    #Filling the matrix
    for(character in 1:ncol(matrix)) {
        #rate
        rate.matrix[character,1]<-anc.state_ace[[character]]$rates
        #std.err
        rate.matrix[character,2]<-anc.state_ace[[character]]$se
    }

    return(rate.matrix)
}
