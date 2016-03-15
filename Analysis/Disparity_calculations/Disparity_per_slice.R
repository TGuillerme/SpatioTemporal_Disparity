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

######################
#Testing with Beck data
######################

#After running the Data_setup script, load the different results

######################
#Tree and matrix
######################

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
tree$root.time<-max(ages.data[,1])

######################
#FADLAD file
######################

#Load the F/LAD for Beck
FADLAD<-read.csv("../Data/Beck_FADLAD.csv", row.names=1)

######################
#Ancestral states reconstruction files
######################

load(paste(data_path, chain_name, "/Beck2014_ancestral_states-claddis.Rda", sep="")) #anc_states

######################
#Distance matrices
######################

#load(paste(data_path, chain_name, "/Beck2014_distance-tips.Rda", sep="")) #dist_tips
load(paste(data_path, chain_name, "/Beck2014_distance-nodes.Rda", sep="")) #dist_nodes
load(paste(data_path, chain_name, "/Beck2014_distance-nodes95.Rda", sep="")) #dist_nodes95


dist.data<-dist_nodes95
include_nodes<-TRUE

#Remove the unaplicable characters
trimmed.max.data<-TrimMorphDistMatrix(dist.data$max.dist.matrix)
#We can see what taxa have been removed by typing:
trimmed.max.data$removed.taxa
#Remove the droped taxa from the tree
tree<-drop.tip(tree, trimmed.max.data$removed.taxa)
#Check gaps in the matrix
any(is.na(trimmed.max.data$dist.matrix))

#PCO

# We can maximise our axes by upping the value "k" (an option in the function) to N - 1 (the maximum number of axes for N objects, i.e., N taxa).
# In addition we want to use another option in the function (add) which gets around the negative eigenvalue problem that can cause downstream problems (e.g., a scree plot with negative values).
# We can specify these options fairly easily and store our answer in a new variable (pco.data) and this time we will just part of the output ($points) which are the values for our taxa on every ordination axis:
pco.data<-cmdscale(trimmed.max.data$dist.matrix, k=nrow(trimmed.max.data$dist.matrix) - 1, add=T)$points

######################
#Disparity
######################

slices<-rev(seq(from=0, to=150, by=10))

pco_sliced_ran<-slice.pco(pco.data, tree, slices, method="random", FAD_LAD=FADLAD, verbose=TRUE)
pco_sliced_acc<-slice.pco(pco.data, tree, slices, method="acctran", FAD_LAD=FADLAD, verbose=TRUE)
pco_sliced_del<-slice.pco(pco.data, tree, slices, method="deltran", FAD_LAD=FADLAD, verbose=TRUE)
pco_sliced_pro<-slice.pco(pco.data, tree, slices, method="proximity", FAD_LAD=FADLAD, verbose=TRUE)
#Calculating the disparity per bins
disparity_ran_table<-time.disparity(pco_sliced_ran, verbose=TRUE)
disparity_acc_table<-time.disparity(pco_sliced_acc, verbose=TRUE)
disparity_del_table<-time.disparity(pco_sliced_del, verbose=TRUE)
disparity_pro_table<-time.disparity(pco_sliced_pro, verbose=TRUE)

#Plot
disparity_data<-disparity_acc_table
#disparity_data[,1]<-as.numeric(disparity_data[,1])

dev.new()
op<-par(mfrow=c(3,2))
plot.disparity(disparity_data, rarefaction=FALSE, xlab="Mya", ylab="Distance from centroid", measure="Cent.dist", las=1)
abline(v=c(10.5), col="red")
plot.disparity(disparity_data, rarefaction=FALSE, xlab="Mya", ylab="Sum of ranges", measure="Sum.range", las=1)
abline(v=c(10.5), col="red")
plot.disparity(disparity_data, rarefaction=FALSE, xlab="Mya", ylab="Sum of variance", measure="Sum.var", las=1)
abline(v=c(10.5), col="red")
plot.disparity(disparity_data, rarefaction=FALSE, xlab="Mya", ylab="Product of ranges", measure="Prod.range", las=1)
abline(v=c(10.5), col="red")
plot.disparity(disparity_data, rarefaction=FALSE, xlab="Mya", ylab="Product of variance", measure="Prod.var", las=1)
abline(v=c(10.5), col="red")
par(op)