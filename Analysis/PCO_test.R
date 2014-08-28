#Header
setwd('~/PhD/Projects/SpatioTemporal_Disparity/Analysis')
library(ape)


#Plotting the polygon outlines per type

plot.pco<-function(pco.scores ,clade.1 , clade.2, ...) {
    library(grDevices)
    plot(1,1, col="white", xlab="PC1", ylab="PC2", xlim=c(min(pco.scores[,1]), max(pco.scores[,1])), ylim=c(min(pco.scores[,2]), max(pco.scores[,2])), ...)
    points(pco.scores[,1], pco.scores[,2], col=c("red", "blue")[as.factor(pco.scores[,clade.col])])
    polygon(clade.1[chull(clade.1),], border="red")
    polygon(clade.2[chull(clade.2),], border="blue")
}

source('treeAge.R')

#data
tree<-read.tree('BDtree.tre')
table<-read.dna('MorphoMat.phylip', as.character=TRUE, as.matrix=TRUE)
tree.age<-treeAge(tree, max(tree$edge.length))
plot(tree) ; axisPhylo()

#Euclidian matrix
eucl.table<-dist(table, method = "euclidean")

#pcoa
pco<-pcoa(eucl.table)

#Extracting clades
clade1<-extract.clade(tree, 88)$tip.label
clade2<-extract.clade(tree, 54)$tip.label
pco.scores<-as.data.frame(pco$vectors)
clade.col<-ncol(pco.scores)+1
pco.scores[, clade.col]<-"NA"
pco.scores[match(clade1, row.names(pco.scores)),clade.col]<-"clade1"
pco.scores[match(clade2, row.names(pco.scores)),clade.col]<-"clade2"
clade.1<-pco.scores[match(clade1, row.names(pco.scores)),1:2]
clade.2<-pco.scores[match(clade2, row.names(pco.scores)),1:2]

dev.new()
plot.pco(pco.scores, clade.1, clade.2, main="Plot all data")


#Extracting only the living taxa
age.0<-tree$tip.label[which(tree.age == 0)]
pco.scores.0<-pco.scores[match(age.0, row.names(pco.scores)),]
clade.1.0<-clade.1[match(age.0, rownames(clade.1))[-which(is.na(match(age.0, rownames(clade.1))))], ]
clade.2.0<-clade.2[match(age.0, rownames(clade.2))[-which(is.na(match(age.0, rownames(clade.2))))], ]

plot.pco(pco.scores.0, clade.1.0, clade.2.0, main="Plot age=0")

#Extracting the taxa at any age
#Now it gets more tricky:
#-must take into account the fossil occurrence time?
#-reconstruct ancestral state per character?
#-redo PCO on that? or do PCO on the whole tree (including ancestors, etc...)



#------------------------
#Ancestral states
#------------------------
verbose=TRUE
ancestral.list<-list()
for (character in 1:ncol(table)) {
    ancestral.list[[character]]<-ace(as.factor(table[,character]), tree, type = "d")
    if(verbose==TRUE) {
        message('.', appendLF=FALSE)
    }
}

#Plotting (visual)
for (character in 1:ncol(table)) {
    plot(tree, label.offset = 0.1, cex=0.5)
    co <- c("blue", "yellow" , "red")
    tiplabels(pch = 22, bg = co[as.numeric(as.factor(table[,character]))], cex = 1)
    nodelabels(thermo = ancestral.list[[character]]$lik.anc, piecol = co, cex = 0.4)
    Sys.sleep(1)
}
#Should test character reconstruction by a 'up an down' method (simulate/estimate)





#------------------------
#Global PCO
#------------------------
#Create a matrix containing the characters for each taxa (real) and the nodes (estimated) and redo the PCO from there.
#Create the probability matrix
#The probability is the probability of element (tip/node/branch) of being in the right (observed) state for a given discrete character
#For all the species present in the matrix with data, this probability = 1.

#Creating the probability matrix for all tips and nodes with probability = 1 for every cell
nodes=tree$Nnode
prob.matrix<-matrix(data=1, ncol=ncol(table), nrow=(nrow(table)+nodes))
row.names(prob.matrix)<-c(row.names(table), paste("n",seq(1:nodes), sep=""))
#TEST
library(testthat)
expect_true(all(prob.matrix[1:nrow(table),] == 1))


#Replacing the probability of each node (was 1) into the highest1 probability calculated in ancestral.list
for (character in 1:ncol(table)) {
    for (edge in 1:nodes) {
        prob.matrix[(nrow(table)+edge),character]<-max(ancestral.list[[character]]$lik.anc[edge,])
    }
}

#Creating the new characters state matrix for every tips and nodes
state.matrix<-matrix(NA, ncol=ncol(table), nrow=(nrow(table)+nodes))
row.names(state.matrix)<-c(row.names(table), paste("n",seq(1:nodes), sep=""))

#Adding the character states for the tips (observed)
state.matrix[1:nrow(table), 1:ncol(table)]<-as.matrix(table)
#TEST
library(testthat)
expect_true(all(state.matrix[1:nrow(table),] == table))

#Adding the character states for the nodes (estimated)
for (character in 1:ncol(table)) {
    for (edge in 1:nodes) {
        #Extracting the column name for each edge and character
        state.matrix[(nrow(table)+edge),character]<-which(ancestral.list[[character]]$lik.anc[edge,]==max(ancestral.list[[character]]$lik.anc[edge,]))[[1]]-1 # -1 is to make the character start at 0 (which greps 1,2,3, etc... instead of 0,1,2, ...)
    }
}
#Test (visual)
for (character in 1:ncol(table)) {
    co<-c("blue", "yellow" , "red")
    op<-par(mfrow=c(1, 2))
    #Plot from the tree (observed)
    plot(tree, label.offset=0.1, cex=0.5, main="From tree")
    tiplabels(pch=22, bg=co[as.numeric(as.factor(table[,character]))], cex=1)
    nodelabels(bg=co[as.numeric(as.factor(round(ancestral.list[[character]]$lik.anc[,2])))], cex=0.3)
    #Plot from the state.matrix (estimated)
    plot(tree, label.offset=0.1, cex=0.5, main="From matrix")
    tiplabels(pch=22, bg=co[as.numeric(as.factor(state.matrix[1:nrow(table),character]))], cex=1)
    nodelabels(bg=co[as.numeric(as.factor(state.matrix[(nrow(table)+1):nrow(state.matrix),character]))], cex=0.3)
    par(op)
    Sys.sleep(1)
}


#Redoing PCO on the state.matrix
pco<-pcoa(dist(state.matrix, method="euclidian"))
#biplot(pco)


#Adding edge names to the trees
tree.full<-tree
tree.full$node.label<-row.names(state.matrix[(nrow(table)+1):nrow(state.matrix),])

clade1<-c(extract.clade(tree.full, 88)$tip.label, extract.clade(tree.full, 88)$node.label)
clade2<-c(extract.clade(tree.full, 54)$tip.label, extract.clade(tree.full, 54)$node.label)
pco.scores<-as.data.frame(pco$vectors)
clade.col<-ncol(pco.scores)+1
pco.scores[, clade.col]<-"NA"
pco.scores[match(clade1, row.names(pco.scores)),clade.col]<-"clade1"
pco.scores[match(clade2, row.names(pco.scores)),clade.col]<-"clade2"
clade.1<-pco.scores[match(clade1, row.names(pco.scores)),1:2]
clade.2<-pco.scores[match(clade2, row.names(pco.scores)),1:2]

dev.new()
plot.pco(pco.scores, clade.1, clade.2, main="with nodes")



#Extracting age:
#getting all nodes or tips that are the given age
#For tips, keep only if age == given.age
#For nodes, take node value if node < edge < tip

#Different methods for extracting branch character state
#ACCTRAN where changes are assigned along branches of a phylogenetic tree as close to the root as possible.
#DELTRAN where changes are assigned along branches as close to the tips as possible.
#RATE where the changes are randomly assigned along the branches according to the ancestral state and the overall of the character.

#When ACCTRAN, the state of the branch is equal to the state of it's upper node or tip
#When DELTRAN, the state of the branch is equal to the state of it's lower node
#When RATE, the state of the branch is randomly sampled with a probability of being the ancestral state being equal to:
#The probability value of the ancestral state multiplied by the probability of not changing this state
#Where the probability of not changing the state of a character is proportional to the exponential of the inverse of the rate (i.e. the rate of not changing) multiplied by the branch length.
#P(not changing the state)=exp(-rate*branch length)
#Therefore, the probability of staying in the ancestral state is equal to:
#P(ancestral state)*P(not changing the state)
#WARNING: this assumes the rates are constant!
#Maybe put that before the slicing on the final pco?

#Slicing plots
source('sliceTree.R')

#Six slices with ACCTRAN
tree0.acc<-slice.tree(tree.full, age=0, method="ACCTRAN")
tree1.acc<-slice.tree(tree.full, age=1, method="ACCTRAN")
tree2.acc<-slice.tree(tree.full, age=2, method="ACCTRAN")
tree3.acc<-slice.tree(tree.full, age=3, method="ACCTRAN")
tree4.acc<-slice.tree(tree.full, age=4, method="ACCTRAN")
tree5.acc<-slice.tree(tree.full, age=5, method="ACCTRAN")

#Six slices with DELTRAN
tree0.del<-slice.tree(tree.full, age=0, method="DELTRAN")
tree1.del<-slice.tree(tree.full, age=1, method="DELTRAN")
tree2.del<-slice.tree(tree.full, age=2, method="DELTRAN")
tree3.del<-slice.tree(tree.full, age=3, method="DELTRAN")
tree4.del<-slice.tree(tree.full, age=4, method="DELTRAN")
tree5.del<-slice.tree(tree.full, age=5, method="DELTRAN")

#make plot pco
pcoSlice<-function(tree, main="slice") {
    pcoSlice<-pco.scores[tree$tip.label,]
    clade.1<-pcoSlice[match(clade1, row.names(pcoSlice)),1:2]
    clade.2<-pcoSlice[match(clade2, row.names(pcoSlice)),1:2]
    plot.pco(pcoSlice, clade.1[-which(is.na(clade.1)),], clade.2[-which(is.na(clade.2)),], main=main)
}


op<-par(mfrow=c(3, 2))
pcoSlice(tree0.del, main="slice 0")
pcoSlice(tree1.del, main="slice 1")
pcoSlice(tree2.del, main="slice 2")
pcoSlice(tree3.del, main="slice 3")
pcoSlice(tree4.del, main="slice 4")
pcoSlice(tree5.del, main="slice 5")
par(op)