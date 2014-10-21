##########################
#slice.rate
##########################
#RATES slicing method for slice.tree
#Choosing between the state of the ancestor or the state of the offspring depending on the rate of the character and then length of the branch.
#v0.1
##########################
#SYNTAX :
#<tree> a 'phylo' object
#<pco.scores> a pco.scores object containing PC axis data and optional taxonomic data.
#<slices> slices ages. Can be either a single value (for a single slice) or a series of values (for multiple slices)
#<method> the slicing method (what becomes of the sliced branches): can be ACCTRAN, DELTRAN or RATE.
##########################
#----
#guillert(at)tcd.ie 20/10/2014
##########################

slice.rate<-function(tree, slice, rate, prob){
    #SANITIZING
    #embedded within std.slice

    #SLICE RATE

    #DEBUG
    warning("DEBUG MODE")
    slice=0

    #Creating the ACCTRAN/DELTRAN table
    #subtrees
    sub_tree_acc<-slice.tree(tree, slice, method="ACCTRAN")
    sub_tree_del<-slice.tree(tree, slice, method="DELTRAN") 
    #subtaxa table
    del_acc.table<-data.frame("DELTRAN"=sub_tree_del$tip.label, "ACCTRAN"=sub_tree_acc$tip.label)

    #From this table, selecting the right tip depending on the rate and the branch length.
    terms <- sub_tree_acc$edge[, 2] <= Ntip(sub_tree_acc)
    terminal.edges <- sub_tree_acc$edge.length[terms]
    names(terminal.edges) <- sub_tree_acc$tip.label[sub_tree_acc$edge[terms, 2]]
    del_acc.table$"brlen"<-terminal.edges

    #Calculate the branch state
    global.state<-vector()
    for (tip in 1:nrow(del_acc.table)) {
        brstate<-vector()
        for (character in 1:ncol(prob)) {
            #state.prob=prob[which(rownames(prob) == del_acc.table[tip,1]) ,character]
            brstate[[character]]<-branch.state(0,1, 1, rate=rate[tip, 1], brlen=del_acc.table[tip, 3])
        }
        global.state[[tip]]<-sum(brstate)/ncol(prob)
    }
    del_acc.table$"global.state"<-global.state

    #Determinate the branch label
    del_acc.table$"tip.label"<-NA
    for (tip in 1:nrow(del_acc.table)) {
        if(del_acc.table[tip,4] > 0.5) {
            del_acc.table[tip,5]<-as.character(del_acc.table[tip,2])
        } else {
            del_acc.table[tip,5]<-as.character(del_acc.table[tip,1])
        }
    }



}


#RATES METHOD:
#Choosing between the state of the ancestor or the state of the offspring depending on the rate of the character and then length of the branch.

#branch.state<-function(anc.stat, anc.stat.prob, off.stat, br.length, trait.rate) {
#    #Randomly generates a state with a likelihood probability for any time point along a branch
#    branch.prob<-anc.sta.prob*exp(-trait.rate*br.length)
#    branch.stat<-sample(c(anc.trait,des.trait), 1, prob=c(branch.prob, (1-branch.prob)))
#}

#RATES<-function(...) {
#    stop("to do")
#    #use branch.state
#}