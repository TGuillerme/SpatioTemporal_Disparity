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
tree_rand<-list()
tree_simc<-list()
if(length(max_null_rand_rm) == 0) {
    tree_rand<-rep(tree, length(null_mat_rand_dist))
} else {
    for (replicate in 1:length(null_mat_rand_dist)) {
        tree_rand[[replicate]]<-drop.tip(tree, max_null_rand_rm[[replicate]])
    }
}

if(length(max_null_simc_rm) == 0) {
    tree_simc<-rep(tree, length(null_mat_rand_dist))
} else {
    for (replicate in 1:length(null_mat_rand_dist)) {
        tree_simc[[replicate]]<-drop.tip(tree, max_null_simc_rm[[replicate]])
    }
}

#Extracting a list of trimmed matrices
trimmed_mat_max_null_rand<-list()
trimmed_mat_max_null_simc<-list()
for (replicate in 1:length(null_mat_rand_dist)) {
    trimmed_mat_max_null_rand[[replicate]]<-trimmed_max_null_rand[[replicate]]$dist.matrix
    trimmed_mat_max_null_simc[[replicate]]<-trimmed_max_null_simc[[replicate]]$dist.matrix
}

#List of trees
#trees<-list("tips"=tree_tips, "nodes"=tree_nodes, "nodes95"=tree_nodes95)

######################
#PCO
######################

pco_data_tips<-cmdscale(trimmed_max_data_tips$dist.matrix, k=nrow(trimmed_max_data_tips$dist.matrix) - 1, add=T)$points
krows<-nrow(trimmed_mat_max_null_rand[[1]])
pco_null_rand<-lapply(trimmed_mat_max_null_rand, cmdscale, k=krows - 1, add=T)
krows<-nrow(trimmed_mat_max_null_simc[[1]])
pco_null_simc<-lapply(trimmed_mat_max_null_simc, cmdscale, k=krows - 1, add=T)

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
    pco_int_null_rand[[replicate]]<-int.pco(pco_null_rand[[replicate]], tree_rand[[replicate]], int_breaks, include.nodes=FALSE, FAD_LAD=FADLAD)
    pco_int_null_simc[[replicate]]<-int.pco(pco_null_simc[[replicate]], tree_simc[[replicate]], int_breaks, include.nodes=FALSE, FAD_LAD=FADLAD)
}


#Calculating the disparity per intervals
disp_int_tips<-time.disparity(pco_int_tips, method="centroid", verbose=TRUE, save.all=TRUE)
disp_int_null_rand<-lapply(pco_int_null_rand, time.disparity, method="centroid", verbose=TRUE, save.all=TRUE)
disp_int_null_simc<-lapply(pco_int_null_simc, time.disparity, method="centroid", verbose=TRUE, save.all=TRUE)

#Isolating the values and the quantiles
disp_tips_values<-disp_int_tips[[2]] ; disp_tips_quantiles<-disp_int_tips[[1]]
rand_quantiles<-list()
rand_values<-list()
simc_quantiles<-list()
simc_values<-list()
for (replicate in 1:length(null_mat_rand_dist)) {
    rand_quantiles[[replicate]]<-disp_int_null_rand[[replicate]][[1]]
    rand_values[[replicate]]<-disp_int_null_rand[[replicate]][[2]]
    simc_quantiles[[replicate]]<-disp_int_null_simc[[replicate]][[1]]
    simc_values[[replicate]]<-disp_int_null_simc[[replicate]][[2]]
}

#Combining the list data together (average among all the replicates)
disp_int_rand<-rand_quantiles[[1]]
disp_int_simc<-simc_quantiles[[1]]
for (replicate in 2:length(null_mat_rand_dist)) {
    disp_int_rand[,2:ncol(disp_int_rand)]<-( disp_int_rand[,2:ncol(disp_int_rand)]+rand_quantiles[[replicate]][,2:ncol(disp_int_rand)] ) /2
    disp_int_simc[,2:ncol(disp_int_simc)]<-( disp_int_simc[,2:ncol(disp_int_simc)]+simc_quantiles[[replicate]][,2:ncol(disp_int_simc)] ) /2
}

#Comparing the distribution for each slice using Bhattacharya
#Creating vector lists
for(interval in 1:length(disp_tips_values)) {
    disp_tips_values[[interval]]<-as.vector(disp_tips_values[[interval]])
}
rand_values_list<-list()
simc_values_list<-list()
for(interval in 1:length(rand_values[[1]])) {
    for(replicate in 1:length(rand_values)) {
        rand_values_list[[interval]]<-as.vector(rand_values[[replicate]][[interval]])
        simc_values_list[[interval]]<-as.vector(simc_values[[replicate]][[interval]])
    }
}

#Calculating the Bhattacharrya coefficients
rand_vs_simc<-pair.bhatt.coeff(rand_values_list, simc_values_list)
rand_vs_obs<-pair.bhatt.coeff(rand_values_list, disp_tips_values)
simc_vs_obs<-pair.bhatt.coeff(simc_values_list, disp_tips_values)

#Plotting the different models
op<-par(mfrow=c(2, 2), bty="l")# oma=c(bottom, left, top, right)
#Centroid
plot.disparity(disp_tips_quantiles, rarefaction=FALSE, xlab="", ylab="Distance from centroid", measure="Cent.dist", main="Observed")
plot.disparity(disp_int_rand, rarefaction=FALSE, xlab="", ylab="Distance from centroid", measure="Cent.dist", main="Random")
plot.disparity(disp_int_simc, rarefaction=FALSE, xlab="", ylab="Distance from centroid", measure="Cent.dist", main="Simulated")
#Bhattacharrya
plot(rand_vs_simc, type="l", ylim=c(0,1), col="black", ylab="Bhattacharrya coefficient", xaxt='n', xlab="")
abline(h=0.975, col="grey", lty=3)
abline(h=0.025, col="grey", lty=3)
axis(side = 1, at=1:length(disp_tips_values), labels=names(disp_tips_values), las=2)
points(rand_vs_obs, type="l", col="red")
points(simc_vs_obs, type="l", col="blue")
xaxis()
par(op)