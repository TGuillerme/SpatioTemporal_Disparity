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
chain_name<-"Beck2014"
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
FADLAD<-read.csv(paste(data_path, chain_name, "_FADLAD.csv", sep=""), row.names=1)

######################
#Ancestral states reconstruction files
######################

#load(paste(data_path, chain_name, "/", chain_name, "_ancestral_states-claddis.Rda", sep="")) #anc_states

######################
#Distance matrices
######################

load(paste(data_path, chain_name, "/", chain_name, "_distance-tips.Rda", sep="")) #dist_tips
load(paste(data_path, chain_name, "/", chain_name, "_null_mat_rand_dist.Rda", sep="")) #null_mat_rand_dist
load(paste(data_path, chain_name, "/", chain_name, "_null_mat_simchar_dist.Rda", sep="")) #null_mat_simchar_dist

#Extracting a list of max.dist.matrix
max_null_rand<-list()
max_null_simc<-list()
for (replicate in 1:length(null_mat_rand_dist)) {
    max_null_rand[[replicate]]<-null_mat_rand_dist[[replicate]]$max.dist.matrix
    max_null_simc[[replicate]]<-null_mat_simchar_dist[[replicate]]$max.dist.matrix
}
#Remove the inapplicable characters
trimmed_max_data_tips<-TrimMorphDistMatrix(dist_tips$max.dist.matrix)
trimmed_max_null_rand<-lapply(max_null_rand, TrimMorphDistMatrix)
trimmed_max_null_simc<-lapply(max_null_simc, TrimMorphDistMatrix)

#Extracting list of removed taxa
max_null_rand_rm<-list()
max_null_simc_rm<-list()
for (replicate in 1:length(null_mat_rand_dist)) {
    max_null_rand_rm[[replicate]]<-trimmed_max_null_rand[[replicate]]$removed.taxa
    max_null_simc_rm[[replicate]]<-trimmed_max_null_simc[[replicate]]$removed.taxa
}

#Remove the dropped taxa from the tree
tree_tips<-drop.tip(tree, trimmed_max_data_tips$removed.taxa)
tree_rand<-lapply(max_nul_rand_rm, function(x) drop.tip(tree, x))
tree_simc<-lapply(max_nul_simc_rm, function(x) drop.tip(tree, x))

#List of trees
#trees<-list("tips"=tree_tips, "nodes"=tree_nodes, "nodes95"=tree_nodes95)

######################
#PCO
######################

pco_data_tips<-cmdscale(trimmed_max_data_tips$dist.matrix, k=nrow(trimmed_max_data_tips$dist.matrix) - 1, add=T)$points
pco_null_rand<-lapply(trimmed_max_null_rand, cmdscale, k=nrow(trimmed_max_null_rand) - 1, add=T)
pco_null_simc<-lapply(trimmed_max_null_simc, cmdscale, k=nrow(trimmed_max_null_simc) - 1, add=T)

#Isolate only the pco scores for the lists
for (replicate in 1:length(null_mat_rand_dist)) {
    pco_null_rand[[replicate]]<-pco_null_rand[[replicate]]$points
    pco_null_simc[[replicate]]<-pco_null_simc[[replicate]]$points
}

#Storing as a list
#pco_data<-list("tips"=pco_data_tips, "nodes"=pco_data_nodes, "nodes95"=pco_data_nodes95)

######################
#Disparity - per interval only for now - centroid only for now
######################

#Calculating the rarefaction
#rarefaction_median<-disparity(pco_data, rarefaction=TRUE, verbose=TRUE, central_tendency=median)

#Generating the different intervals PCOs
pco_int_tips<-int.pco(pco_data_tips, tree_tips, int_breaks, include.nodes=FALSE, FAD_LAD=FADLAD)

pco_int_null_rand<-list()
pco_int_null_simc<-list()
for (replicate in 1:length(null_mat_rand_dist)) {
    pco_int_null_rand[[replicate]]<-int.pco(pco_int_null_rand[[replicate]], tree_rand[[replicate]], int_breaks, include.nodes=FALSE, FAD_LAD=FADLAD)
    pco_int_null_simc[[replicate]]<-int.pco(pco_int_null_simc[[replicate]], tree_simc[[replicate]], int_breaks, include.nodes=FALSE, FAD_LAD=FADLAD)
}

#Calculating the disparity per intervals
disp_int_tips<-time.disparity(pco_int_tips, method="centroid", verbose=TRUE)
disp_int_null_rand<-lapply(pco_int_null_rand, time.disparity, , method="centroid", verbose=TRUE)
disp_int_null_simc<-lapply(pco_int_null_simc, time.disparity, , method="centroid", verbose=TRUE)

#Combining the list data together

#Plotting the different models