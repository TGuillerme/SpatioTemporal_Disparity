#FUNCTIONS FOR slice.tree

#Select the parent node of a tip
slice.tree_parent.node<-function(tree, tip) {
    #Selecting parent edge in the full tree
    parent.edge<-tree$edge[which(tree$edge[,2] == grep(tip, c(tree$tip.label, tree$node.label))[1]), 1]
    #Selecting parent node in the full tree
    parent_node<-tree$node.label[parent.edge-Ntip(tree)]
    #error if not working
    if (length(parent_node) != 1) {
        stop('No parent node found.')
    }
    return(parent_node)
}

#Select the offspring node/tip of a node towards a tip
slice.tree_offspring.node<-function(tree, parent_node, tip) {
    #Stop if parent node is the same as tip
    if(parent_node == tip) {
        stop('Parent node is the tip!')
    }
    #Extracting the subtrees connected to the parent node
    offsprings<-tree$edge[which(tree$edge[,1] == (match(parent_node, tree$node.label)+Ntip(tree))), 2]
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


#Modify the tree slicing by replacing the tips that are not at the cut by the parent node
slice.tree_DELTRAN<-function(tree, tip, tree_slice) {
    parent_node<-slice.tree_parent.node(tree, tip)
    while(match(parent_node, tree_slice$node.label, nomatch=FALSE) == 0) {
        #Repeat if slice.tree_parent.node is not present in the sliced tree
        parent_node<-slice.tree_parent.node(tree, parent_node)
    }
    return(parent_node)
}

#Modify the tree slicing by replacing the tips that are not at the cut by the offspring node towards the tip
slice.tree_ACCTRAN<-function(tree, tip, tree_slice) {
    parent_node<-slice.tree_DELTRAN(tree, tip, tree_slice)
    offspring_node<-slice.tree_offspring.node(tree, parent_node, tip)
    return(offspring_node)
}
