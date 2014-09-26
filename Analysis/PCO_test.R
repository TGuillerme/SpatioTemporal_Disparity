#Header
setwd('~/PhD/Projects/SpatioTemporal_Disparity/Analysis')
library(ape)
setwd('~/PhD/Projects/SpatioTemporal_Disparity/Functions')
source('~/PhD/Projects/SpatioTemporal_Disparity/Functions/refresh.R')


#Plotting the polygon outlines per type

source('~/PhD/Projects/SpatioTemporal_Disparity/Functions/treeAge.R')
source('~/PhD/Projects/SpatioTemporal_Disparity/Functions/sliceTree.R')
source('~/PhD/Projects/SpatioTemporal_Disparity/Functions/plot.fun.R')

#data
tree<-read.tree('BDtree.tre')
table<-read.dna('MorphoMat.phylip', as.character=TRUE, as.matrix=TRUE)
#table<-read.dna('GappedMat.phylip', as.character=TRUE, as.matrix=TRUE)
#table<-table[26:51,] #for a 75% missing data table
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

op<-par(mfrow=c(3, 2))
pcoSlice(tree0.acc, main="slice 0")
pcoSlice(tree1.acc, main="slice 1")
pcoSlice(tree2.acc, main="slice 2")
pcoSlice(tree3.acc, main="slice 3")
pcoSlice(tree4.acc, main="slice 4")
pcoSlice(tree5.acc, main="slice 5")
par(op)




#------------------------
#Empirical data - Slater
#------------------------
#And now with real data
Slater.table<-read.table("../Data/2013-Slater-MEE-morpho.phylip", header=F, sep=" ", row.names=1) 
Slater.tree<-read.nexus('../Data/2013-Slater-MEE-TEM.tre')
plot(Slater.tree, show.tip.label=FALSE) ; axisPhylo()

#Correct the table
#Remove species only present in the table
library(caper)
missing.species<-comparative.data(Slater.tree, data.frame("species"=row.names(Slater.table), "dummy"=rnorm(nrow(Slater.table)), "dumb"=rnorm(nrow(Slater.table))), "species")$dropped$unmatched.rows
suppressWarnings(table.tmp<-Slater.table[-c(which(row.names(Slater.table) == missing.species)),])
bugged<-c("Abrocomidae","Aegialodon","Murtoilestes") #not sure why

#Remove species with only missing data from the table
missing<-which(duplicated(table.tmp)==TRUE)
missing.na.species<-rownames(table.tmp[missing,])
table.tmp<-table.tmp[-missing,]
table<-as.matrix(table.tmp)

#list of missing species (to drop in the tree)
missing.data.species<-c(missing.species, missing.na.species, bugged) #not sure why?

#TEST
library(testthat)
expect_equal((nrow(table.tmp)+length(missing.species)+ length(missing)), nrow(Slater.table))
#expect_that(Ntip(Slater.tree), equals(nrow(table.tmp)))
expect_is(table, "matrix")


table<-table[-which(rownames(table) == bugged),] #for some reason Abrocomidae is not removed

#Replacing "?" by NA's
#for (character in 1:ncol(table)) {
#    character.corrected<-as.factor(table[,character])
#    wrong.levels<-levels(character.corrected)
#    levels(character.corrected)<-c(NA, wrong.levels[-1])
#    table[,character]<-character.corrected
#}
#if(any(grep("?", table.tmp))) {
#    expect_true(anyNA(table))
#} else {
#    expect_false(anyNA(table))
#}

#68, 71
stable<-table[-69,]
table<-stable
#stable<-table[-c(1, 68,71), ]
eucl.table<-dist(table, method = "euclidean")
pco<-pcoa(eucl.table)
expect_is(pco, "pcoa")


#------------------------
#Empirical data - Beck
#------------------------

Beck.tree<-read.nexus('~/PhD/Projects/SpatioTemporal_Disparity/Empirical_mammals/2014-Beck-ProcB-TEM.tre')
Beck.table<-read.table('~/PhD/Projects/SpatioTemporal_Disparity/Empirical_mammals/2014-Beck-Morpho.table', header=F, sep=" ", row.names=1)
plot(Beck.tree, cex=0.5) ; axisPhylo()


library(caper)
missing.in.data<-comparative.data(Beck.tree, data.frame("species"=row.names(Beck.table), "dummy"=rnorm(nrow(Beck.table)), "dumb"=rnorm(nrow(Beck.table))), "species")$dropped$unmatched.rows
missing.in.tree<-comparative.data(Beck.tree, data.frame("species"=row.names(Beck.table), "dummy"=rnorm(nrow(Beck.table)), "dumb"=rnorm(nrow(Beck.table))), "species")$dropped$tips

remove<-NULL
for (n in 1:length(missing.in.data)) {
    remove<-c(remove, which(row.names(Beck.table) == missing.in.data[n]))
}
table.tmp<-Beck.table[-remove,]

#missing<-which(duplicated(table.tmp)==TRUE)
#missing.na.species<-rownames(table.tmp[missing,])
#table.tmp<-table.tmp[-missing,]
table<-as.matrix(table.tmp)

#missing.data.species<-c(missing.species, missing.na.species)


eucl.table<-dist(table, method = "euclidean")
pco<-pcoa(eucl.table)
expect_is(pco, "pcoa")

#Remove the species with only missing data from the tree
tree.tmp<-drop.tip(Beck.tree, missing.in.tree)
#TEST
expect_that(Ntip(tree.tmp), equals(nrow(table)))

#Making the tree fully dichotomous
tree<-multi2di(tree.tmp)
expect_true(is.binary.tree(tree))
#Adding a small branch length to the newly generated dichotomies (0.1*min(branch length))
tree$edge.length[which(tree$edge.length == 0)]<-0.1*min(tree$edge.length[-which(tree$edge.length == 0)])
expect_that(length(which(tree$edge.length == 0)), equals(0))
#Adding tree nodes names

#Calculate the tree age
tree.age<-treeAge(tree)

#plot the cleaned tree
plot(tree, cex=0.5) ; axisPhylo()

#Testing pco
#Euclidian matrix
eucl.table<-dist(table, method = "euclidean")

#pcoa
pco<-pcoa(eucl.table)













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

#------------------------
#Creating the prob and state
#------------------------
nodes=tree$Nnode
prob.matrix<-matrix(data=1, ncol=ncol(table), nrow=(nrow(table)+nodes))
row.names(prob.matrix)<-c(row.names(table), paste("n",seq(1:nodes), sep=""))

#Replacing the probability of each node (was 1) into the highest probability calculated in ancestral.list
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

#Adding the character states for the nodes (estimated)
for (character in 1:ncol(table)) {
    for (edge in 1:nodes) {
        #Extracting the column name for each edge and character
        max.prob<-which(ancestral.list[[character]]$lik.anc[edge,]==max(ancestral.list[[character]]$lik.anc[edge,]))[[1]]
        state.matrix[(nrow(table)+edge),character]<-colnames(ancestral.list[[character]]$lik.anc)[max.prob] # -1 is to make the character start at 0 (which greps 1,2,3, etc... instead of 0,1,2, ...)
    }
}

#Test all columns (no extra levels)
for(character in 1:ncol(table)) {
    expect_that(levels(as.factor(table[,character])) %in% as.vector(state.matrix[,character]), equals(rep(TRUE, length(levels(as.factor(table[,character]))))))
}

n.tips=1:Ntip(tree)
n.nodes=(Ntip(tree)+1):(Ntip(tree)+Nnode(tree))

expect_equal(Ntip(tree), nrow(table))
expect_equal((Ntip(tree)+Nnode(tree)), nrow(state.matrix))
expect_equal((nrow(state.matrix)-Ntip(tree)), Nnode(tree))
expect_equal(length(n.tips), Ntip(tree))
expect_equal(length(n.nodes), Nnode(tree))


#Visual test
for (character in 1:ncol(table)) {
    co<-palette()
    op<-par(mfrow=c(1, 2))
    #Plot from the tree (observed)
    plot(tree, label.offset=0.1, cex=0.5, main=paste("From tree ", character, sep="") , show.tip.label=FALSE)
    tiplabels(pch=22, bg=co[as.numeric(as.factor(table[,character]))], cex=1)    
    #Extract the state
    anc.states<-vector()
    for (row in 1:nrow(ancestral.list[[character]]$lik.anc)) {
        anc.states[[row]]<-which(ancestral.list[[character]]$lik.anc[row,] == max(ancestral.list[[character]]$lik.anc[row,]))[[1]]
    }
    nodelabels(bg=co[as.numeric(anc.states)], cex=0.3)
    #Plot from the state.matrix (estimated)
    plot(tree, label.offset=0.1, cex=0.5, main=paste("From matrix ", character, sep=""), show.tip.label=FALSE)
    character.levels<-levels(as.factor(state.matrix[,character]))
    tiplabels(pch=22, bg=co[factor(state.matrix[n.tips,character], levels=character.levels)], cex=1)  #1:nrow(table)
    nodelabels(bg=co[factor(state.matrix[n.nodes,character], levels=character.levels)], cex=0.3) #(nrow(table))+1):nrow(state.matrix)
    par(op)
    Sys.sleep(1)
}



tree.full<-tree
tree.full$node.label<-row.names(state.matrix[(nrow(table)+1):nrow(state.matrix),])
plot(tree.full, cex=0.5) ; axisPhylo() ; nodelabels(text=tree.full$node.label, cex=0.4)

#------------------------
#Global PCO
#------------------------

dist.matrix<-dist(state.matrix, method="euclidian")
pco<-pcoa(dist.matrix) #Error in eigen(delta1) : infinite or missing values in 'x'
#biplot(pco)


#------------------------
#PCO plot
#------------------------
#Adding edge names to the trees
tree.full<-tree
tree.full$node.label<-row.names(state.matrix[(nrow(table)+1):nrow(state.matrix),])

clade1<-extract.clade(tree.full, 153)
clade2<-drop.tip(tree.full, clade1$tip.label)
pco.scores<-as.data.frame(pco$vectors)
clade.col<-ncol(pco.scores)+1
pco.scores[, clade.col]<-"NA"
pco.scores[match(clade1$tip.label, row.names(pco.scores)), clade.col]<-"eutheria"
pco.scores[match(clade2$tip.label, row.names(pco.scores)),clade.col]<-"non-eutheria"
pco.scores[match(clade1$node.label, row.names(pco.scores)), clade.col]<-"eutheria"
pco.scores[match(clade2$node.label, row.names(pco.scores)),clade.col]<-"non-eutheria"
clade.1<-pco.scores[which(pco.scores[,clade.col] == "eutheria"),1:2]
clade.2<-pco.scores[which(pco.scores[,clade.col] == "non-eutheria"),1:2]
expect_that(sort(row.names(clade.1)), equals(sort(c(clade1$tip.label, clade1$node.label))))
expect_that(sort(row.names(clade.2)), equals(sort(c(clade2$tip.label, clade2$node.label))))

dev.new()

plot.pco(pco.scores, clade.1, clade.2, main="with nodes")

#------------------------
#time slices
#------------------------
#Slicing plots
source('sliceTree.R')

#Six slices with ACCTRAN
tree0.acc<-slice.tree(tree.full, age=0, method="ACCTRAN")
tree1.acc<-slice.tree(tree.full, age=18, method="ACCTRAN")
tree2.acc<-slice.tree(tree.full, age=37, method="ACCTRAN")
tree3.acc<-slice.tree(tree.full, age=56, method="ACCTRAN")
tree4.acc<-slice.tree(tree.full, age=75, method="ACCTRAN")
tree5.acc<-slice.tree(tree.full, age=93, method="ACCTRAN")
tree6.acc<-slice.tree(tree.full, age=112, method="ACCTRAN")
tree7.acc<-slice.tree(tree.full, age=131, method="ACCTRAN")
tree8.acc<-slice.tree(tree.full, age=150, method="ACCTRAN")

#Six slices with DELTRAN
tree0.del<-slice.tree(tree.full, age=0, method="DELTRAN")
tree1.del<-slice.tree(tree.full, age=18, method="DELTRAN")
tree2.del<-slice.tree(tree.full, age=37, method="DELTRAN")
tree3.del<-slice.tree(tree.full, age=56, method="DELTRAN")
tree4.del<-slice.tree(tree.full, age=75, method="DELTRAN")
tree5.del<-slice.tree(tree.full, age=93, method="DELTRAN")
tree6.del<-slice.tree(tree.full, age=112, method="DELTRAN")
tree7.del<-slice.tree(tree.full, age=131, method="DELTRAN")
tree8.del<-slice.tree(tree.full, age=150, method="DELTRAN")

#make plot pco


op<-par(mfrow=c(3, 3))
pcoSlice(tree0.del, main="0 Mya")
pcoSlice(tree1.del, main="18 Mya")
pcoSlice(tree2.del, main="37 Mya")
pcoSlice(tree3.del, main="56 Mya")
pcoSlice(tree4.del, main="75 Mya")
pcoSlice(tree5.del, main="93 Mya")
pcoSlice(tree6.del, main="112 Mya")
pcoSlice(tree7.del, main="131 Mya")
pcoSlice(tree8.del, main="150 Mya")
par(op)


op<-par(mfrow=c(3, 3))
pcoSlice(tree0.acc, main="0 Mya")
pcoSlice(tree1.acc, main="18 Mya")
pcoSlice(tree2.acc, main="37 Mya")
pcoSlice(tree3.acc, main="56 Mya")
pcoSlice(tree4.acc, main="75 Mya")
pcoSlice(tree5.acc, main="93 Mya")
pcoSlice(tree6.acc, main="112 Mya")
pcoSlice(tree7.acc, main="131 Mya")
pcoSlice(tree8.acc, main="150 Mya")
par(op)