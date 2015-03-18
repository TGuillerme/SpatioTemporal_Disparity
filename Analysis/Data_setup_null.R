#Script for doing ancestral states reconstruction and distance matrix from an input nexus file

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

###################
#Reading the files
###################

#Selecting the file
chain_name<-"Beck2014"
data_path<-"../Data/"
file_matrix<-"../Data/2014-Beck-ProcB-matrix-morpho.nex"
file_tree<-"../Data/2014-Beck-ProcB-TEM.tre"

#matrix
Nexus_data<-ReadMorphNexus(file_matrix)
Nexus_matrix<-Nexus_data$matrix
#tree
Tree_data<-read.nexus(file_tree)

######################################
#Cleaning the matrices and the trees
######################################

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

#Isolating the states list
states_list<-apply(Nexus_data$matrix, 2, states.count) #states.count function is available in the sanitizing functions

#Generating null matrices
null_matrices_random<-null.data(tree=tree, matrix=states_list, matrix.model="random", replicates=10, verbose=TRUE)
null_matrices_simcha<-null.data(tree=tree, matrix=states_list, matrix.model="sim.char", replicates=10, verbose=TRUE)

#Recreating the proper nexus format file with the new matrices
null_mat_random<-null_mat_simcha<-list()
for (replicate in 1:10) {
    null_mat_random[[replicate]]<-Nexus_data
    null_mat_random[[replicate]]$matrix<-null_matrices_random[[replicate]]
    null_mat_simcha[[replicate]]<-Nexus_data
    null_mat_simcha[[replicate]]$matrix<-null_matrices_simcha[[replicate]]
}

####################################
#Ancestral states reconstruction - Skipped for testing
####################################

#anc_states<-anc.state(tree, Nexus_data, method='ML-ape', verbose=TRUE)
#save(anc_states, file=paste("../Data/",chain_name,"_ancestral_states-ape.Rda", sep=""))

#anc_states<-anc.state(tree, Nexus_data, method='ML-claddis', verbose=TRUE)
#save(anc_states, file=paste(data_path, chain_name, "/",chain_name,"_ancestral_states-claddis.Rda", sep=""))

####################################
#Distance matrix
####################################

#Distance matrix using tips only
message("\nCalculating the distance matrix for the tips only...", appendLF=FALSE)
null_mat_rand_dist<-lapply(null_mat_random, MorphDistMatrix.verbose, verbose=TRUE) 
null_mat_simchar_dist<-lapply(null_mat_simcha, MorphDistMatrix.verbose, verbose=TRUE) 
message("Done.\n", appendLF=FALSE)
save(null_mat_rand_dist, file=paste(data_path, chain_name, "/",chain_name,"_null_mat_rand_dist.Rda", sep="")) #null_mat_rand_dist
save(null_mat_simchar_dist, file=paste(data_path, chain_name, "/",chain_name,"_null_mat_simchar_dist.Rda", sep="")) #null_mat_simchar_dist

save(dist_tips, file=paste(data_path, chain_name, "/",chain_name,"_distance-tips.Rda", sep="")) #dist_tips
