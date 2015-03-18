#Counting the number of states per characters
states.count<-function(character) {
    #Isolate the states
    states<-levels(as.factor(character))
    #Check if multi states
    if(length(grep("&", states)) > 0) {
        #Isolating the multi states
        multi_states<-states[grep("&", states)]
        multi_states<-as.factor(unlist(strsplit(multi_states, split="&")))
        #Removing the multi states from the states list
        states<-states[-grep("&", states)]
        #Check if any of the multi states is not yet present in the states list
        if(any(is.na(match(levels(multi_states), states)))) {
            states<-c(states, multi_states[which(is.na(match(levels(multi_states), states)))])
        }
    }
    #Count the number of states
    return(length(states))
}

#Generating random tree parameters (birth death)
gen.param.tree<-function(x=1) {
    lambda<-runif(x)
    mu<-runif(x,0,lambda)
    return(cbind(lambda, mu))
}

#Adding a root time and node labels to a tree (for lapply loops)
lapply.root<-function(tree, root) {
    tree$root.time<-root
    tree$node.label<-paste("n",seq(1:Nnode(tree)), sep="")
    return(tree)
}

#Generates a Q matrix
Q.matrix<-function(states, rate) {
    #Creating the empty matrix
    Q_mat<-matrix(nrow=states, ncol=states, data=0)
    colnames(Q_mat)<-rownames(Q_mat)<-0:(states-1)
    #Adding the rates to the matrix (rates are equal)
    Q_mat[upper.tri(Q_mat)]<-Q_mat[lower.tri(Q_mat)]<-rate/(states-1)
    diag(Q_mat)<- -(states-1)*(rate/(states-1))
    return(Q_mat)
}


#Generates a fully random matrix
random.mat<-function(n.tips, characters, states) {
    #creates and empty matrix
    rand_matrix<-matrix(nrow=n.tips, ncol=characters, data=NA)
    #loop through the characters
    for (character in 1:characters) {
        rand_matrix[,character]<-sample(0:(states[character]-1), size=n.tips, replace=TRUE)
    }
    return(rand_matrix)
}

#Renaming the column names for all the matrices
renaming.rows<-function(matrix, names) {
    rownames(matrix)<-names
    return(matrix)
}

#Generates a simulated matrix
simulate.mat<-function(phy, n.tips, matrix_characters, matrix_states, max.mat.rate) {
    #Creating the empty matrix
    rand_matrix<-matrix(nrow=n.tips, ncol=matrix_characters, data=NA)
    #character loop
    for(character in 1:matrix_characters) {
        if(matrix_states[character] == 2) {
            #binary
            rand_matrix[,character]<-sim.character(phy, rep(runif(1, 0, max.mat.rate),2), x0=sample(0:1, 1), model="mk2")
        } else {
            #Create a Q matrix
            Q_mat<-Q.matrix(matrix_states[[character]], runif(1, 0, max.mat.rate))
            #Change the states from 0:k to 1:k+1
            rownames(Q_mat)<-colnames(Q_mat)<-as.numeric(colnames(Q_mat))+1
            #not binary
            rand_matrix[,character]<-sim.character(phy, Q_mat, x0=sample(1:matrix_states[character], 1), model="mkn")-1
        }
    }
    return(rand_matrix)
}