#Probability of a point in the branch to be in the ancestral state
P.anc<-function(state.prob, rate, brlen) {
    #exp(rate*brlen) is the expectation of changing state.
    P.anc<-state.prob*(exp(-rate*brlen))
    return(P.anc)
}

#Randomly generates a state with a likelihood probability for any time point along a branch
branch.state<-function(anc.state, off.state, state.prob, rate, brlen) {
    branch.state<-sample(c(anc.state, off.state), size=1, prob=c(P.anc(state.prob, rate, brlen), (1-P.anc(state.prob, rate, brlen))))
    return(branch.state)
}

