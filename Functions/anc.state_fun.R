#FUNCTIONS FOR anc.state

#Ancestral states estimations from a matrix
anc.state_ace<-function(tree, matrix, model, verbose) {
    anc.list<-list()
    if(verbose == TRUE) {
        message('Estimating the ancestral states for ', ncol(matrix), ' characters:', appendLF=FALSE)
    }
    for (character in 1:ncol(matrix)) {
        if(model == 'ML') {
            anc.list[[character]]<-ace(as.factor(matrix[,character]), tree, type="d")
            if(verbose == TRUE) {
                message('.', appendLF=FALSE)
            }
        }
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
            state.matrix[(nrow(matrix)+edge),character]<-colnames(anc.state_ace[[character]]$lik.anc)[max.prob] # -1 is to make the character start at 0 (which greps 1,2,3, etc... instead of 0,1,2, ...)
        }
    }
    return(state.matrix)
}
