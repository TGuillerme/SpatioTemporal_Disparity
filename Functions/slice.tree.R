##########################
#slice.tree
##########################
#Slices a tree given a specific age
#Modyfied timeSliceTree function from Davd Bapst (paleotrree)
#v0.1
##########################
#SYNTAX :
#<tree> a 'phylo' object
#<age> where to slice the tree
#<method> the slicing method (what becomes of the sliced branches): can be ACCTRAN, DELTRAN or RATE
##########################
#----
#guillert(at)tcd.ie 25/09/2014
##########################

slice.tree<-function(tree, age, method) {

    #SANITIZING
    #tree
    check.class(tree, 'phylo', ' must be a phylo object.')
    #must have node labels
    if(is.null(tree$node.label)) {
        stop('The tree must have node label names.')
    }

    #age
    check.class(age, 'numeric', " must be numeric.")
    #age must be at least higher than the root age
    if(age > max(tree.age(tree)[,1])) {
        stop('Age is lower than the root age!')
    }

    #method
    check.class(method, 'character', " must be \'ACCTRAN\', \'DELTRAN\' or \'RATES\'.")
    if(method != "ACCTRAN") {
        if(method != "DELTRAN") {
            if(method != "RATES") {
                        stop(as.character(substitute(method)), " must be \'ACCTRAN\', \'DELTRAN\' or \'RATES\'." , call.=FALSE)
            }
        }
    }

    #SLICING A TREE
    #Creating the tree.age matrix
    tree_age<-tree.age(tree)

    #Running the timeSliceTree function (remove warning, called as a message in the original function)
    suppressMessages(
        tree_slice<-timeSliceTree(tree, age, drop.extinct=TRUE, plot=FALSE)
    )
    #Error with trees with two taxa
    if(Ntip(tree_slice) < 3) {
        stop('To few taxa for the tree slice at age ', age, '!', call.=FALSE)
    }

    #Selecting the tips
    tips<-tree_slice$tip.label

    #renaming the tree_slice
    tree_sliced<-tree_slice

    #Correcting the sliced tree
    for (tip in 1:Ntip(tree_slice)) {

        #Check if the tree is sliced at the exact age of a tip (e.g. time=0)
        if(tree_age[which(tree_age[,2]==tips[tip]),1] != age) {
            if(method == "DELTRAN") {
                tree_sliced$tip.label[tip]<-slice.tree_DELTRAN(tree, tips[tip], tree_slice)
            }
            if(method == "ACCTRAN") {
                tree_sliced$tip.label[tip]<-slice.tree_ACCTRAN(tree, tips[tip], tree_slice)
            }
            if(method == "RATES") {
                stop("TO DO!")
            }
        }
    }

    return(tree_sliced)

}