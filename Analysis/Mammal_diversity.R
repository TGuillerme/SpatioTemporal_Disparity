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

#Living mammals diversity through time
mammal_tree<-read.nexus("../Data/FritzTree.rs200k.1tree.tre")

#Add node labels and root time
tree<-lapply.root(mammal_tree, max(tree.age(mammal_tree)$age))

#Slice tree per million years to get the diversity
slices<-rev(seq(from=0, to=100, by=1))

#Slicing the tree and counting the tips
slice_count<-NULL
slice_count<-list()

for(slice in 1:length(slices)) {
    sub_tree<-timeSliceTree(tree, slices[slice], drop.extinct=TRUE, plot=FALSE)
    slice_count[slice]<-Ntip(sub_tree)
    message(".", appendLF=FALSE)
}
