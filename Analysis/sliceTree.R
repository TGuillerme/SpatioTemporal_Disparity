#Modyfied timeSliceTree function from Davd Bapst (paleotrree)


#setwd('~/PhD/Projects/SpatioTemporal_Disparity/Analysis')
#library(ape)
#tree<-read.tree('BDtree.tre')
#source('PCO_test.R')
#source('treeAge.R')
#tree.full<-tree
#tree.full$node.label<-row.names(state.matrix[(nrow(table)+1):nrow(state.matrix),])

#plot(tree.full, cex=0.7)
#nodelabels(cex=0.7)
#axisPhylo()
#library(paleotree)


#age<-4
#treeAge.4<-timeSliceTree(tree.full, age, drop.extinct=TRUE, plot=FALSE)
#plot(treeAge.4)
#nodelabels()
#axisPhylo()

parent.node<-function(tree, tip) {
    #Selecting parent edge in the full tree
    parent.edge<-tree$edge[which(tree$edge[,2] == grep(tip, c(tree$tip.label, tree$node.label))[1]), 1]
    #Selecting parent node in the full tree
    parent_node<-tree$node.label[parent.edge-Ntip(tree)]
    return(parent_node)
}

offspring.node<-function(tree, parent_node, tip) {
    #Extracting the subtrees connected to the parent node
    offsprings<-tree$edge[which(tree$edge[,1] == (grep(parent_node, tree$node.label)+Ntip(tree))), 2]
    #Testing which subtree contains tip
    for (node in 1:length(offsprings)) {
        #Check if the "node" is a node or a tip
        if(offsprings[node] > Ntip(tree)) {
            subtree<-extract.clade(tree, offsprings[node])
            if(length(grep(tip, subtree$tip.label))==1) {
                offspring.edge<-offsprings[node]
            }
        } else {
            subtree<-tree$tip.label[offsprings[node]]
            if(length(grep(tip, subtree))==1) {
                offspring.edge<-offsprings[node]
            }           
        }
    }
    #Returning the name of the offspring node
    if(offspring.edge > Ntip(tree)) {
        offspring_node<-tree$node.label[offspring.edge-Ntip(tree)]
    } else {
        offspring_node<-tree$tip.label[offspring.edge]
    }
    return(offspring_node)
}


DELTRAN<-function(tree, tip, tree_slice) {
    parent_node<-parent.node(tree, tip)
    while(length(grep(parent_node, tree_slice$node.label)) != 1) {
        #Repeat if parent.node is not present in the sliced tree
        parent_node<-parent.node(tree, parent_node)
    }
    return(parent_node)
}

ACCTRAN<-function(tree, tip, tree_slice) {
    parent_node<-DELTRAN(tree, tip, tree_slice)
    offspring_node<-offspring.node(tree, parent_node, tip)
    return(offspring_node)
}


check.slice.tree.method<-function(method, msg){
    if(method != "ACCTRAN") {
        if(method != "DELTRAN") {
            if(method != "RATES") {
                        stop(as.character(substitute(method)), msg , call.=FALSE)
            }
        }
    }
}

branch.state<-function(anc.stat, anc.stat.prob, off.stat, br.length, trait.rate) {
    #Randomly generates a state with a likelihood probability for any time point along a branch
    branch.prob<-anc.sta.prob*exp(-trait.rate*br.length)
    branch.stat<-sample(c(anc.trait,des.trait), 1, prob=c(branch.prob, (1-branch.prob)))
}

RATES<-function(...) {
    stop("to do")
    #use branch.state
}


#tree<-tree.full
#age<-4
#method="DELTRAN"

slice.tree<-function(tree, age, method) {
    #function for doing sharp time slices (using Dave Bapst's timSliceTree function)
    #method must be ACCTRAN, DELTRAN or RATES
    require(paleotree)

    #SANITIZING
    #needs tree with node labels
    #ages must be numeric
    source('~/Packaging/mulTree/R/check.class.R')
    check.class(age, 'numeric', " must be numeric.")
    check.slice.tree.method(method, " must be 'ACCTRAN', 'DELTRAN' or 'RATES'.")

    #Creating the treeAge matrix
    tree.age<-treeAge(tree)

    #Running the timeSliceTree function (remove warning, called as a message in the original function)
    suppressMessages(tree_slice<-timeSliceTree(tree, age, drop.extinct=TRUE, plot=FALSE))

    #Selecting the tips
    tips<-tree_slice$tip.label

    #renaming the tips on the tree according the the method"
    #ACCTRAN transforms the tip into the offspring node/tip in the tree
    #DELTRAN transforms the tip into the parent node/tip in the tree
    #RATES recalculates the tip with a trait state likelihood probability for any time point along a branch - NEED TO RECALCULATE THE PCO MATRIX


    tree_sliced<-tree_slice

    for (tip in 1:Ntip(tree_slice)) {
        #Check if the tree is sliced at the exact age of a tip (e.g. time=0)
        if(tree.age[grep(tips[tip], tree.age[,2]),1] != age) {
            if(method == "DELTRAN") {
                tree_sliced$tip.label[tip]<-DELTRAN(tree, tips[tip], tree_slice)
            }
            if(method == "ACCTRAN") {
                tree_sliced$tip.label[tip]<-ACCTRAN(tree, tips[tip], tree_slice)
            }
            if(method == "RATES") {
                stop("TO DO!")
            }
        }
    }

    return(tree_sliced)


}







#TEST parent/offspring.node
tr<-rtree(10)
tr$node.label<-paste("n",seq(1:9), sep="")
tip<-tr$tip.label[5]
#Should be NULL (no node.label)
test.p1<-parent.node(rtree(5), 1)
#Should be NULL (no node.label)
test.p2<-parent.node(rtree(5), tip)
#Should be character
test.p3<-parent.node(tr, 1)
#Should be a character
test.p4<-parent.node(tr, tip)
#Should be 0 (parent of the root)
test.p5<-parent.node(tr, "n1")

#Should be a character
test.o1<-offspring.node(tr, "n1", tip)
#Should be a character
test.o2<-offspring.node(tr, test.p4, tip)

library(testthat)
test_that("Testing parent.node", {
    expect_null(test.p1)
    expect_null(test.p2)
    expect_is(test.p3, "character")
    expect_that(length(test.p3), equals(1))
    expect_is(test.p4, "character")
    expect_that(length(test.p4), equals(1))
    expect_that(test.p5, equals(character(0)))
})
test_that("Testing offspring.node", {
    expect_that(length(test.o1), equals(1))
    expect_is(test.o1, "character")
    expect_that(length(test.o2), equals(1))   
    expect_is(test.o2, "character")
})

