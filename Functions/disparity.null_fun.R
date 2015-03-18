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
Q.matrix<-function(states, max.rate=0.5) {
    #Creating the empty matrix
    Q_mat<-matrix(nrow=states, ncol=states, data=0)
    colnames(Q_mat)<-rownames(Q_mat)<-0:(states-1)
    #Random rate generation (rates are equal)
    rate<-runif(1, 0, max.rate)
    #Adding the rates to the matrix (rates are equal)
    Q_mat[upper.tri(Q_mat)]<-Q_mat[lower.tri(Q_mat)]<-rate
    diag(Q_mat)<-1-rate-1
    return(Q_mat)
}
set.seed(1)
bla<-sim.character(phy, rep(rate,2), x0=0, model="mkn")
set.seed(1)
bli<-sim.character(phy, Q_mat, x0=0, model="mkn")

#Generates a fully random matrix
full.random.mat<-function(tree, states) {}