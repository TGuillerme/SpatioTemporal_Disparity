
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
int_breaks<-rev(seq(from=5, to=85, by=10))
int_breaks[length(int_breaks)]<-0
slices<-rev(seq(from=5, to=85, by=10))
slices[length(slices)]<-0
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
Nexus_data$matrix<-Nexus_data$matrix[match(rownames(Nexus_data$matrix), tree$tip.label, nomatch=0),]

#Isolating the states list
states_list<-apply(Nexus_data$matrix, 2, states.count) #states.count function is available in the sanitizing functions

######################
#Generating all the models
######################

observed_mat<-Nexus_data
observed_tree<-tree

load(file=paste("../Data/",chain_name,"/ran.mat_obs.tre_init.Rda", sep=""))
load(file=paste("../Data/",chain_name,"/sim.mat_obs.tre_init.Rda", sep=""))
load(file=paste("../Data/",chain_name,"/obs.mat_yul.tre_init.Rda", sep=""))
load(file=paste("../Data/",chain_name,"/obs.mat_bde.tre_init.Rda", sep=""))

#ran.mat_yul.tre_init<-null.data(tree="yule", matrix=states_list, matrix.model="random", replicates=1, verbose=TRUE, root.time=tree$root.time, n.tips=Ntip(tree))
#ran.mat_bde.tre_init<-null.data(tree="bd", matrix=states_list, matrix.model="random", replicates=1, verbose=TRUE, root.time=tree$root.time, n.tips=Ntip(tree))
#sim.mat_yul.tre_init<-null.data(tree="yule", matrix=states_list, matrix.model="sim.char", replicates=1, verbose=TRUE, root.time=tree$root.time, n.tips=Ntip(tree))
#sim.mat_bde.tre_init<-null.data(tree="bd", matrix=states_list, matrix.model="sim.char", replicates=1, verbose=TRUE, root.time=tree$root.time, n.tips=Ntip(tree))

####################################
#Ancestral states reconstruction - Fast version (ACE)
####################################

load(file=paste("../Data/",chain_name,"/ace_obs.mat_obs.tre.Rda", sep=""))
load(file=paste("../Data/",chain_name,"/ace_ran.mat_obs.tre.Rda", sep=""))
load(file=paste("../Data/",chain_name,"/ace_sim.mat_obs.tre.Rda", sep=""))
load(file=paste("../Data/",chain_name,"/ace_obs.mat_yul.tre.Rda", sep=""))
load(file=paste("../Data/",chain_name,"/ace_obs.mat_bde.tre.Rda", sep=""))

####################################
#Distance matrix
####################################

#Distance matrix using also nodes
load(file=paste("../Data/",chain_name,"/dist_obs.mat_obs.tre.Rda", sep=""))
load(file=paste("../Data/",chain_name,"/dist_ran.mat_obs.tre.Rda", sep=""))
load(file=paste("../Data/",chain_name,"/dist_sim.mat_obs.tre.Rda", sep=""))
load(file=paste("../Data/",chain_name,"/dist_obs.mat_yul.tre.Rda", sep=""))
load(file=paste("../Data/",chain_name,"/dist_obs.mat_bde.tre.Rda", sep=""))

#Distance matrix using also nodes95
load(file=paste("../Data/",chain_name,"/dist_obs.mat_obs.tre95.Rda", sep=""))
load(file=paste("../Data/",chain_name,"/dist_ran.mat_obs.tre95.Rda", sep=""))
load(file=paste("../Data/",chain_name,"/dist_sim.mat_obs.tre95.Rda", sep=""))
load(file=paste("../Data/",chain_name,"/dist_obs.mat_yul.tre95.Rda", sep=""))
load(file=paste("../Data/",chain_name,"/dist_obs.mat_bde.tre95.Rda", sep=""))

######################
#Distance matrices
######################

#Extracting a list of max.dist.matrix
max_obs.mat_obs.tre<-dist_obs.mat_obs.tre$max.dist.matrix
max_ran.mat_obs.tre<-max_sim.mat_obs.tre<-max_obs.mat_yul.tre<-max_obs.mat_bde.tre<-list()

for (replicate in 1:length(dist_ran.mat_obs.tre)) {
    max_ran.mat_obs.tre[[replicate]]<-dist_ran.mat_obs.tre[[replicate]]$max.dist.matrix
    max_sim.mat_obs.tre[[replicate]]<-dist_sim.mat_obs.tre[[replicate]]$max.dist.matrix
    max_obs.mat_yul.tre[[replicate]]<-dist_obs.mat_yul.tre[[replicate]]$max.dist.matrix
    max_obs.mat_bde.tre[[replicate]]<-dist_obs.mat_bde.tre[[replicate]]$max.dist.matrix
}

#Remove the inapplicable characters
trimmed_obs.mat_obs.tre<-TrimMorphDistMatrix(max_obs.mat_obs.tre)
trimmed_ran.mat_obs.tre<-lapply(max_ran.mat_obs.tre, TrimMorphDistMatrix)
trimmed_sim.mat_obs.tre<-lapply(max_sim.mat_obs.tre, TrimMorphDistMatrix)
trimmed_obs.mat_yul.tre<-lapply(max_obs.mat_yul.tre, TrimMorphDistMatrix)
trimmed_obs.mat_bde.tre<-lapply(max_obs.mat_bde.tre, TrimMorphDistMatrix)

#Extracting list of removed taxa
#max_null_rand_rm<-list()
#max_null_simc_rm<-list()
#for (replicate in 1:length(null_mat_rand_dist)) {
#    max_null_rand_rm[[replicate]]<-trimmed_max_null_rand[[replicate]]$removed.taxa
#    max_null_simc_rm[[replicate]]<-trimmed_max_null_simc[[replicate]]$removed.taxa
#}

#Remove the dropped taxa from the tree
#tree_tips<-drop.tip(tree, trimmed_max_data_tips$removed.taxa)
#tree_rand<-list()
#tree_simc<-list()
#if(length(max_null_rand_rm) == 0) {
#    tree_rand<-rep(tree, length(null_mat_rand_dist))
#} else {
#    for (replicate in 1:length(null_mat_rand_dist)) {
#        tree_rand[[replicate]]<-drop.tip(tree, max_null_rand_rm[[replicate]])
#    }
#}

#if(length(max_null_simc_rm) == 0) {
#    tree_simc<-rep(tree, length(null_mat_rand_dist))
#} else {
#    for (replicate in 1:length(null_mat_rand_dist)) {
#        tree_simc[[replicate]]<-drop.tip(tree, max_null_simc_rm[[replicate]])
#    }
#}

#Extracting a list of trimmed matrices
matrix_obs.mat_obs.tre<-trimmed_obs.mat_obs.tre$dist.matrix
matrix_ran.mat_obs.tre<-matrix_sim.mat_obs.tre<-matrix_obs.mat_yul.tre<-matrix_obs.mat_bde.tre<-list()
for (replicate in 1:length(dist_ran.mat_obs.tre)) {
    matrix_ran.mat_obs.tre[[replicate]]<-trimmed_ran.mat_obs.tre[[replicate]]$dist.matrix
    matrix_sim.mat_obs.tre[[replicate]]<-trimmed_sim.mat_obs.tre[[replicate]]$dist.matrix
    matrix_obs.mat_yul.tre[[replicate]]<-trimmed_obs.mat_yul.tre[[replicate]]$dist.matrix
    matrix_obs.mat_bde.tre[[replicate]]<-trimmed_obs.mat_bde.tre[[replicate]]$dist.matrix
}

#List of trees
#trees<-list("tips"=tree_tips, "nodes"=tree_nodes, "nodes95"=tree_nodes95)

######################
#PCO
######################

pco_obs.mat_obs.tre<-cmdscale(matrix_obs.mat_obs.tre, k=nrow(matrix_obs.mat_obs.tre) - 1, add=TRUE)$points
pco_ran.mat_obs.tre<-lapply(matrix_ran.mat_obs.tre, function(X) cmdscale(X, k=nrow(X)-1, add=TRUE))
pco_sim.mat_obs.tre<-lapply(matrix_sim.mat_obs.tre, function(X) cmdscale(X, k=nrow(X)-1, add=TRUE))
pco_obs.mat_yul.tre<-lapply(matrix_obs.mat_yul.tre, function(X) cmdscale(X, k=nrow(X)-1, add=TRUE))
pco_obs.mat_bde.tre<-lapply(matrix_obs.mat_bde.tre, function(X) cmdscale(X, k=nrow(X)-1, add=TRUE))

#Isolate only the pco scores for the lists
for (replicate in 1:length(dist_ran.mat_obs.tre)) {
    pco_ran.mat_obs.tre[[replicate]]<-pco_ran.mat_obs.tre[[replicate]]$points
    pco_sim.mat_obs.tre[[replicate]]<-pco_sim.mat_obs.tre[[replicate]]$points
    pco_obs.mat_yul.tre[[replicate]]<-pco_obs.mat_yul.tre[[replicate]]$points
    pco_obs.mat_bde.tre[[replicate]]<-pco_obs.mat_bde.tre[[replicate]]$points
}
######################
#Disparity - per interval only for now - centroid only for now #PROBLEM WITH INTERVALS!
######################

#Generating the different intervals PCOs
intpco_obs.mat_obs.tre<-int.pco(pco_obs.mat_obs.tre, observed_tree, int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD, diversity=TRUE)
div_obs.mat_obs.tre<-intpco_obs.mat_obs.tre[[2]] ; intpco_obs.mat_obs.tre<-intpco_obs.mat_obs.tre[[1]]

intpco_ran.mat_obs.tre<-intpco_sim.mat_obs.tre<-intpco_obs.mat_yul.tre<-intpco_obs.mat_bde.tre<-list()
div_ran.mat_obs.tre<-div_sim.mat_obs.tre<-div_obs.mat_yul.tre<-div_obs.mat_bde.tre<-list()
for (replicate in 1:length(dist_ran.mat_obs.tre)) {
    intpco_ran.mat_obs.tre<-int.pco(pco_ran.mat_obs.tre[[replicate]], observed_tree, int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD, diversity=TRUE)
    div_ran.mat_obs.tre[[replicate]]<-intpco_obs.mat_obs.tre[[2]] ; intpco_obs.mat_obs.tre[[replicate]]<-intpco_obs.mat_obs.tre[[1]]

    intpco_sim.mat_obs.tre<-int.pco(pco_sim.mat_obs.tre[[replicate]], observed_tree, int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD, diversity=TRUE)
    div_sim.mat_obs.tre[[replicate]]<-intpco_sim.mat_obs.tre[[2]] ; intpco_sim.mat_obs.tre[[replicate]]<-intpco_sim.mat_obs.tre[[1]]

    intpco_obs.mat_yul.tre<-int.pco(pco_obs.mat_yul.tre[[replicate]], obs.mat_yul.tre_init[[replicate]], int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD, diversity=TRUE)
    div_obs.mat_yul.tre[[replicate]]<-intpco_obs.mat_yul.tre[[2]] ; intpco_obs.mat_yul.tre[[replicate]]<-intpco_obs.mat_yul.tre[[1]]

    intpco_obs.mat_bde.tre<-int.pco(pco_obs.mat_bde.tre[[replicate]], obs.mat_bde.tre_init[[replicate]], int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD, diversity=TRUE)
    div_obs.mat_bde.tre[[replicate]]<-intpco_obs.mat_bde.tre[[2]] ; intpco_obs.mat_bde.tre[[replicate]]<-intpco_obs.mat_bde.tre[[1]]
}

#Calculating the disparity per intervals
disp_obs.mat_obs.tre<-time.disparity(intpco_obs.mat_obs.tre, method="centroid", verbose=TRUE, save.all=TRUE)
disp_ran.mat_obs.tre<-lapply(intpco_ran.mat_obs.tre, time.disparity, method="centroid", verbose=TRUE, save.all=TRUE)
disp_sim.mat_obs.tre<-lapply(intpco_sim.mat_obs.tre, time.disparity, method="centroid", verbose=TRUE, save.all=TRUE)
disp_obs.mat_yul.tre<-lapply(intpco_obs.mat_yul.tre, time.disparity, method="centroid", verbose=TRUE, save.all=TRUE)
disp_obs.mat_bde.tre<-lapply(intpco_obs.mat_bde.tre, time.disparity, method="centroid", verbose=TRUE, save.all=TRUE)

#Isolating the values and the quantiles
values_obs.mat_obs.tre<-disp_obs.mat_obs.tre[[2]] ; quantiles_obs.mat_obs.tre<-disp_obs.mat_obs.tre[[1]]

values_ran.mat_obs.tre<-values_sim.mat_obs.tre<-values_obs.mat_yul.tre<-values_obs.mat_bde.tre<-list()
quantiles_ran.mat_obs.tre<-quantiles_sim.mat_obs.tre<-quantiles_obs.mat_yul.tre<-quantiles_obs.mat_bde.tre<-list()

for (replicate in 1:length(dist_ran.mat_obs.tre)) {
    values_ran.mat_obs.tre[[replicate]]<-disp_ran.mat_obs.tre[[replicate][[2]] ; quantiles_ran.mat_obs.tre[[replicate]]<-disp_ran.mat_obs.tre[[replicate][[1]]
    values_sim.mat_obs.tre[[replicate]]<-disp_sim.mat_obs.tre[[replicate][[2]] ; quantiles_sim.mat_obs.tre[[replicate]]<-disp_sim.mat_obs.tre[[replicate][[1]]
    values_obs.mat_yul.tre[[replicate]]<-disp_obs.mat_yul.tre[[replicate][[2]] ; quantiles_obs.mat_yul.tre[[replicate]]<-disp_obs.mat_yul.tre[[replicate][[1]]
    values_obs.mat_bde.tre[[replicate]]<-disp_obs.mat_bde.tre[[replicate][[2]] ; quantiles_obs.mat_bde.tre[[replicate]]<-disp_obs.mat_bde.tre[[replicate][[1]]
}

#Combining the list data together (average among all the replicates)
quant_ran.mat_obs.tre<-quantiles_ran.mat_obs.tre[[1]]
quant_sim.mat_obs.tre<-quantiles_sim.mat_obs.tre[[1]]
quant_obs.mat_yul.tre<-quantiles_obs.mat_yul.tre[[1]]
quant_obs.mat_bde.tre<-quantiles_obs.mat_bde.tre[[1]]

for (replicate in 2:length(dist_ran.mat_obs.tre)) {
    quant_ran.mat_obs.tre[, 2:ncol(quant_ran.mat_obs.tre)]<-( quant_ran.mat_obs.tre[, 2:ncol(quant_ran.mat_obs.tre)] + quantiles_ran.mat_obs.tre[[replicate]][, 2:ncol(quant_ran.mat_obs.tre)] ) /2
    quant_sim.mat_obs.tre[, 2:ncol(quant_sim.mat_obs.tre)]<-( quant_sim.mat_obs.tre[, 2:ncol(quant_sim.mat_obs.tre)] + quantiles_sim.mat_obs.tre[[replicate]][, 2:ncol(quant_sim.mat_obs.tre)] ) /2
    quant_obs.mat_yul.tre[, 2:ncol(quant_obs.mat_yul.tre)]<-( quant_obs.mat_yul.tre[, 2:ncol(quant_obs.mat_yul.tre)] + quantiles_obs.mat_yul.tre[[replicate]][, 2:ncol(quant_obs.mat_yul.tre)] ) /2
    quant_obs.mat_bde.tre[, 2:ncol(quant_obs.mat_bde.tre)]<-( quant_obs.mat_bde.tre[, 2:ncol(quant_obs.mat_bde.tre)] + quantiles_obs.mat_bde.tre[[replicate]][, 2:ncol(quant_obs.mat_bde.tre)] ) /2
}

#Comparing the distribution for each slice using Bhattacharya
#Creating vector lists
val_obs.mat_obs.tre<-list()
for(interval in 1:length(values_obs.mat_obs.tre)) {
    val_obs.mat_obs.tre[[interval]]<-as.vector(values_obs.mat_obs.tre[[interval]])
}

val_ran.mat_obs.tre<-val_sim.mat_obs.tre<-val_obs.mat_yul.tre<-val_obs.mat_bde.tre<-list()
for(replicate in 1:length(val_ran.mat_obs.tre)) {
    for(interval in 1:length(val_ran.mat_obs.tre[[replicate]])) {
        val_ran.mat_obs.tre[[interval]]<-as.vector(val_ran.mat_obs.tre[[replicate]][[interval]])
    }
}
for(replicate in 1:length(val_sim.mat_obs.tre)) {
    for(interval in 1:length(val_sim.mat_obs.tre[[replicate]])) {
        val_sim.mat_obs.tre[[interval]]<-as.vector(val_sim.mat_obs.tre[[replicate]][[interval]])
    }
}
for(replicate in 1:length(val_obs.mat_yul.tre)) {
    for(interval in 1:length(val_obs.mat_yul.tre[[replicate]])) {
        val_obs.mat_yul.tre[[interval]]<-as.vector(val_obs.mat_yul.tre[[replicate]][[interval]])
    }
}
for(replicate in 1:length(val_obs.mat_bde.tre)) {
    for(interval in 1:length(val_obs.mat_bde.tre[[replicate]])) {
        val_obs.mat_bde.tre[[interval]]<-as.vector(val_obs.mat_bde.tre[[replicate]][[interval]])
    }
}

#Calculating the Bhattacharrya coefficients
obs_vs_ran.mat<-pair.bhatt.coeff(val_obs.mat_obs.tre, val_ran.mat_obs.tre)
obs_vs_sim.mat<-pair.bhatt.coeff(val_obs.mat_obs.tre, val_sim.mat_obs.tre)
obs_vs_yul.tre<-pair.bhatt.coeff(val_obs.mat_obs.tre, val_obs.mat_yul.tre)
obs_vs_bde.tre<-pair.bhatt.coeff(val_obs.mat_obs.tre, val_obs.mat_bde.tre)
#
sim_vs_ran.mat<-pair.bhatt.coeff(val_sim.mat_obs.tre, val_ran.mat_obs.tre)
bde_vs_yul.tre<-pair.bhatt.coeff(val_obs.mat_bde.tre, val_obs.mat_yul.tre)



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