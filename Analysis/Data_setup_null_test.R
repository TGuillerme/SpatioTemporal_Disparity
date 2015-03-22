#Script for testing the properties of the cladistic space

#Setwd
if(length(grep("TGuillerme", getwd()))) {
    setwd('~/PhD/Projects/SpatioTemporal_Disparity/Analysis')
} else {
    warning("You might have to change the directory!")
}
if(length(grep("SpatioTemporal_Disparity/Analysis", getwd()))==0) {
    if(length(grep("SpatioTemporal_Disparity-master/Analysis", getwd()))==0) {
        stop("Wrong directory!\nThe current directory must be:\nSpatioTemporal_Disparity/Analysis/ OR SpatioTemporal_Disparity-master/Analysis/\nYou can clone the whole repository from:\nhttps://github.com/TGuillerme/SpatioTemporal_Disparity")
    }
}

#Load the functions and the packages
source("functions.R")

#Load Beck data (up to line 159)
#BECK 2014 ProcB
chain_name<-"null_test"
data_path<-"../Data/"
file_matrix<-"../Data/2014-Beck-ProcB-matrix-morpho.nex"
file_tree<-"../Data/2014-Beck-ProcB-TEM.tre"
int_breaks<-rev(seq(from=0, to=150, by=20))+5
int_breaks[length(int_breaks)]<-0
slices<-rev(seq(from=0, to=150, by=10))
KT_bin=4.5
KT_sli=9.5

######################
#Tree and matrix
######################

#matrix
Nexus_data<-ReadMorphNexus(file_matrix)
Nexus_matrix<-Nexus_data$matrix
#tree
Tree_data<-read.nexus(file_tree)

#Cleaning the matrices and the trees
#Remove species with only missing data before hand
if (any(apply(as.matrix(Nexus_matrix), 1, function(x) levels(as.factor((x)))) == "?")) {
    Nexus_matrix<-Nexus_matrix[-c(as.vector(which(apply(as.matrix(Nexus_matrix), 1, function(x) levels(as.factor(x))) == "?"))),]
}

#Cleaning the tree and the table
#making the saving folder
tree<-clean.tree(Tree_data, Nexus_matrix)
table<-clean.table(Nexus_matrix, Tree_data)
Nexus_data$matrix<-table

#Forcing the tree to be binary
tree<-bin.tree(tree)

#Adding node labels to the tree
tree$node.label<-paste("n",seq(1:Nnode(tree)), sep="")
#Setting the tree root age
ages_data<-tree.age(tree)
tree$root.time<-max(ages_data[,1])

######################
#FADLAD file
######################

#Load the F/LAD for Beck
FADLAD<-read.csv(paste(data_path, "Beck2014_FADLAD.csv", sep=""), row.names=1)

######################
#Selecting only stems
######################

#subtree
tree<-extract.clade(tree, node=150)
ages_data<-tree.age(tree)
tree$root.time<-max(ages_data[,1])

#submatrix
match(tree$tip.label, rownames(Nexus_data$matrix))
Nexus_data$matrix<-Nexus_data$matrix[match(tree$tip.label, rownames(Nexus_data$matrix)) ,]

#Isolating the states list
states_list<-apply(Nexus_data$matrix, 2, states.count) #states.count function is available in the sanitizing functions

######################
#Generating all the models
######################

observed_mat<-Nexus_data
observed_tree<-tree
ran.mat_obs.tre_init<-null.data(tree=observed_tree, matrix=states_list, matrix.model="random", replicates=20, verbose=TRUE)
save(ran.mat_obs.tre_init, file=paste("../Data/",chain_name,"/ran.mat_obs.tre_init.Rda", sep=""))

sim.mat_obs.tre_init<-null.data(tree=observed_tree, matrix=states_list, matrix.model="sim.char", replicates=20, verbose=TRUE)
save(sim.mat_obs.tre_init, file=paste("../Data/",chain_name,"/sim.mat_obs.tre_init.Rda", sep=""))

#obs.mat_yul.tre_init<-null.data(tree="yule", matrix=observed_mat$matrix, replicates=10, verbose=TRUE, root.time=tree$root.time)
#save(obs.mat_yul.tre_init, file=paste("../Data/",chain_name,"/obs.mat_yul.tre_init.Rda", sep=""))

#obs.mat_bde.tre_init<-null.data(tree="bd", matrix=observed_mat$matrix, replicates=10, verbose=TRUE, root.time=tree$root.time)
#save(obs.mat_bde.tre_init, file=paste("../Data/",chain_name,"/obs.mat_bde.tre_init.Rda", sep=""))

#ran.mat_yul.tre_init<-null.data(tree="yule", matrix=states_list, matrix.model="random", replicates=1, verbose=TRUE, root.time=tree$root.time, n.tips=Ntip(tree))
#ran.mat_bde.tre_init<-null.data(tree="bd", matrix=states_list, matrix.model="random", replicates=1, verbose=TRUE, root.time=tree$root.time, n.tips=Ntip(tree))
#sim.mat_yul.tre_init<-null.data(tree="yule", matrix=states_list, matrix.model="sim.char", replicates=1, verbose=TRUE, root.time=tree$root.time, n.tips=Ntip(tree))
#sim.mat_bde.tre_init<-null.data(tree="bd", matrix=states_list, matrix.model="sim.char", replicates=1, verbose=TRUE, root.time=tree$root.time, n.tips=Ntip(tree))

#Recreating the proper nexus format file with the new matrices
ran.mat_obs.tre<-sim.mat_obs.tre<-obs.mat_yul.tre<-obs.mat_bde.tre<-ran.mat_yul.tre<-ran.mat_bde.tre<-sim.mat_yul.tre<-sim.mat_bde.tre<-list()
for (replicate in 1:length(ran.mat_obs.tre_init)) {
    ran.mat_obs.tre[[replicate]]<-observed_mat
    ran.mat_obs.tre[[replicate]]$matrix<-ran.mat_obs.tre_init[[replicate]]
    sim.mat_obs.tre[[replicate]]<-observed_mat
    sim.mat_obs.tre[[replicate]]$matrix<-sim.mat_obs.tre_init[[replicate]]
    #obs.mat_yul.tre[[replicate]]<-observed_mat
    #obs.mat_yul.tre[[replicate]]$matrix<-obs.mat_yul.tre_init[[replicate]]
    #obs.mat_bde.tre[[replicate]]<-observed_mat
    #obs.mat_bde.tre[[replicate]]$matrix<-obs.mat_bde.tre_init[[replicate]]
    #ran.mat_yul.tre[[replicate]]<-observed_mat
    #ran.mat_yul.tre[[replicate]]$matrix<-ran.mat_yul.tre_init[[replicate]]
    #ran.mat_bde.tre[[replicate]]<-observed_mat
    #ran.mat_bde.tre[[replicate]]$matrix<-ran.mat_bde.tre_init[[replicate]]
    #sim.mat_yul.tre[[replicate]]<-observed_mat
    #sim.mat_yul.tre[[replicate]]$matrix<-sim.mat_yul.tre_init[[replicate]]
    #sim.mat_bde.tre[[replicate]]<-observed_mat
    #sim.mat_bde.tre[[replicate]]$matrix<-sim.mat_bde.tre_init[[replicate]]
}

####################################
#Ancestral states reconstruction - Fast version (ACE)
####################################

ace_obs.mat_obs.tre<-anc.state(observed_tree, observed_mat, method='ML-ape', verbose=TRUE)
ace_obs.mat_obs.tre$state<-apply(ace_obs.mat_obs.tre$state, 2, replace.na) #replace.na comes from the sanitizing functions
save(ace_obs.mat_obs.tre, file=paste("../Data/",chain_name,"/ace_obs.mat_obs.tre.Rda", sep=""))


ace_ran.mat_obs.tre<-ace_sim.mat_obs.tre<-ace_obs.mat_yul.tre<-ace_obs.mat_bde.tre<-list()
for (replicate in 1:length(ran.mat_obs.tre_init)) {
    ace_ran.mat_obs.tre[[replicate]]<-anc.state(observed_tree, ran.mat_obs.tre[[replicate]], method='ML-ape', verbose=TRUE)
    ace_ran.mat_obs.tre[[replicate]]$state<-apply(ace_ran.mat_obs.tre[[replicate]]$state, 2, replace.na)

    ace_sim.mat_obs.tre[[replicate]]<-anc.state(observed_tree, sim.mat_obs.tre[[replicate]], method='ML-ape', verbose=TRUE)
    ace_sim.mat_obs.tre[[replicate]]$state<-apply(ace_sim.mat_obs.tre[[replicate]]$state, 2, replace.na)

    #ace_obs.mat_yul.tre[[replicate]]<-anc.state(obs.mat_yul.tre_init[[replicate]], observed_mat, method='ML-ape', verbose=TRUE)
    #ace_obs.mat_yul.tre[[replicate]]$state<-apply(ace_obs.mat_yul.tre[[replicate]]$state, 2, replace.na)

    #ace_obs.mat_bde.tre[[replicate]]<-anc.state(obs.mat_bde.tre_init[[replicate]], observed_mat, method='ML-ape', verbose=TRUE)
    #ace_obs.mat_bde.tre[[replicate]]$state<-apply(ace_obs.mat_bde.tre[[replicate]]$state, 2, replace.na)
}
save(ace_ran.mat_obs.tre, file=paste("../Data/",chain_name,"/ace_ran.mat_obs.tre.Rda", sep=""))
save(ace_sim.mat_obs.tre, file=paste("../Data/",chain_name,"/ace_sim.mat_obs.tre.Rda", sep=""))
#save(ace_obs.mat_yul.tre, file=paste("../Data/",chain_name,"/ace_obs.mat_yul.tre.Rda", sep=""))
#save(ace_obs.mat_bde.tre, file=paste("../Data/",chain_name,"/ace_obs.mat_bde.tre.Rda", sep=""))

#Adding nodes to the nexus matrices
observed_mat95<-observed_mat
observed_mat95$matrix<-anc.unc(ace_obs.mat_obs.tre, 0.95, missing=NA)$state
observed_mat$matrix<-ace_obs.mat_obs.tre$state

ran.mat_obs.tre95<-ran.mat_obs.tre
sim.mat_obs.tre95<-sim.mat_obs.tre
#obs.mat_yul.tre95<-obs.mat_yul.tre
#obs.mat_bde.tre95<-obs.mat_bde.tre
for (replicate in 1:length(ran.mat_obs.tre_init)) {
    ran.mat_obs.tre95[[replicate]]$matrix<-anc.unc(ace_ran.mat_obs.tre[[replicate]], 0.95, missing=NA)$state
    ran.mat_obs.tre[[replicate]]$matrix<-ace_ran.mat_obs.tre[[replicate]]$state

    sim.mat_obs.tre95[[replicate]]$matrix<-anc.unc(ace_sim.mat_obs.tre[[replicate]], 0.95, missing=NA)$state
    sim.mat_obs.tre[[replicate]]$matrix<-ace_sim.mat_obs.tre[[replicate]]$state

    #obs.mat_yul.tre95[[replicate]]$matrix<-anc.unc(ace_obs.mat_yul.tre[[replicate]], 0.95, missing=NA)$state
    #obs.mat_yul.tre[[replicate]]$matrix<-ace_obs.mat_yul.tre[[replicate]]$state

    #obs.mat_bde.tre95[[replicate]]$matrix<-anc.unc(ace_obs.mat_bde.tre[[replicate]], 0.95, missing=NA)$state
    #obs.mat_bde.tre[[replicate]]$matrix<-ace_obs.mat_bde.tre[[replicate]]$state
}

####################################
#Distance matrix
####################################

#Distance matrix using also nodes
dist_obs.mat_obs.tre<-MorphDistMatrix.verbose(observed_mat, verbose=TRUE)
save(dist_obs.mat_obs.tre, file=paste("../Data/",chain_name,"/dist_obs.mat_obs.tre.Rda", sep=""))

dist_ran.mat_obs.tre<-lapply(ran.mat_obs.tre, MorphDistMatrix.verbose, verbose=TRUE)
save(dist_ran.mat_obs.tre, file=paste("../Data/",chain_name,"/dist_ran.mat_obs.tre.Rda", sep=""))

dist_sim.mat_obs.tre<-lapply(sim.mat_obs.tre, MorphDistMatrix.verbose, verbose=TRUE) 
save(dist_sim.mat_obs.tre, file=paste("../Data/",chain_name,"/dist_sim.mat_obs.tre.Rda", sep=""))

#dist_obs.mat_yul.tre<-lapply(obs.mat_yul.tre, MorphDistMatrix.verbose, verbose=TRUE) 
#save(dist_obs.mat_yul.tre, file=paste("../Data/",chain_name,"/dist_obs.mat_yul.tre.Rda", sep=""))

#dist_obs.mat_bde.tre<-lapply(obs.mat_bde.tre, MorphDistMatrix.verbose, verbose=TRUE) 
#save(dist_obs.mat_bde.tre, file=paste("../Data/",chain_name,"/dist_obs.mat_bde.tre.Rda", sep=""))

#Distance matrix using also nodes95
#Distance matrix using also nodes
dist_obs.mat_obs.tre95<-MorphDistMatrix.verbose(observed_mat, verbose=TRUE)
save(dist_obs.mat_obs.tree95, file=paste("../Data/",chain_name,"/dist_obs.mat_obs.tree95.Rda", sep=""))

dist_ran.mat_obs.tre95<-lapply(ran.mat_obs.tre95, MorphDistMatrix.verbose, verbose=TRUE)
save(dist_ran.mat_obs.tre95, file=paste("../Data/",chain_name,"/dist_ran.mat_obs.tree95.Rda", sep=""))

dist_sim.mat_obs.tre95<-lapply(sim.mat_obs.tre95, MorphDistMatrix.verbose, verbose=TRUE) 
save(dist_sim.mat_obs.tre95, file=paste("../Data/",chain_name,"/dist_sim.mat_obs.tree95.Rda", sep=""))

#dist_obs.mat_yul.tre95<-lapply(obs.mat_yul.tre95, MorphDistMatrix.verbose, verbose=TRUE) 
#save(dist_obs.mat_yul.tre95, file=paste("../Data/",chain_name,"/dist_obs.mat_yul.tree95.Rda", sep=""))

#dist_obs.mat_bde.tre95<-lapply(obs.mat_bde.tre95, MorphDistMatrix.verbose, verbose=TRUE) 
#save(dist_obs.mat_bde.tre95, file=paste("../Data/",chain_name,"/dist_obs.mat_bde.tree95.Rda", sep=""))