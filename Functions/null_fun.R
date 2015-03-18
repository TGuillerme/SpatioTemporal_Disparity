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