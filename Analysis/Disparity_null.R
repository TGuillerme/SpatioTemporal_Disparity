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

load(paste(data_path, chain_name, "/", chain_name, "_ancestral_states-claddis.Rda", sep="")) #anc_states

######################
#Distance matrices
######################

load(paste(data_path, chain_name, "/", chain_name, "_distance-tips.Rda", sep="")) #dist_tips
load(paste(data_path, chain_name, "/", chain_name, "_distance-nodes.Rda", sep="")) #dist_nodes
load(paste(data_path, chain_name, "/", chain_name, "_distance-nodes95.Rda", sep="")) #dist_nodes95

#Remove the inapplicable characters
trimmed_max_data_tips<-TrimMorphDistMatrix(dist_tips$max.dist.matrix)
trimmed_max_data_nodes<-TrimMorphDistMatrix(dist_nodes$max.dist.matrix)
trimmed_max_data_nodes95<-TrimMorphDistMatrix(dist_nodes95$max.dist.matrix)
#Remove the dropped taxa from the tree
tree_tips<-drop.tip(tree, trimmed_max_data_tips$removed.taxa)
tree_nodes<-drop.tip(tree, trimmed_max_data_nodes$removed.taxa)
tree_nodes95<-drop.tip(tree, trimmed_max_data_nodes95$removed.taxa)
#Remove the eventual inapplicable nodes
trimmed_max_data_nodes$dist.matrix<-trimmed_max_data_nodes$dist.matrix[c(tree_nodes$tip.label, tree_nodes$node.label),c(tree_nodes$tip.label, tree_nodes$node.label)]
trimmed_max_data_nodes95$dist.matrix<-trimmed_max_data_nodes95$dist.matrix[c(tree_nodes$tip.label, tree_nodes$node.label),c(tree_nodes95$tip.label, tree_nodes95$node.label)]

#List of trees
trees<-list("tips"=tree_tips, "nodes"=tree_nodes, "nodes95"=tree_nodes95)

######################
#PCO
######################

pco_data_tips<-cmdscale(trimmed_max_data_tips$dist.matrix, k=nrow(trimmed_max_data_tips$dist.matrix) - 1, add=T)$points
pco_data_nodes<-cmdscale(trimmed_max_data_nodes$dist.matrix, k=nrow(trimmed_max_data_nodes$dist.matrix) - 1, add=T)$points
pco_data_nodes95<-cmdscale(trimmed_max_data_nodes95$dist.matrix, k=nrow(trimmed_max_data_nodes95$dist.matrix) - 1, add=T)$points

#Storing as a list
pco_data<-list("tips"=pco_data_tips, "nodes"=pco_data_nodes, "nodes95"=pco_data_nodes95)

######################
#Disparity
######################

#Calculating the rarefaction
#rarefaction_median<-disparity(pco_data, rarefaction=TRUE, verbose=TRUE, central_tendency=median)

#Generating the different intervals PCOs
pco_int_tips<-int.pco(pco_data_tips, tree_tips, int_breaks, include.nodes=FALSE, FAD_LAD=FADLAD)
pco_int_nodes<-int.pco(pco_data_nodes, tree_nodes, int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD)
pco_int_nodes95<-int.pco(pco_data_nodes95, tree_nodes95, int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD)

#Calculating the disparity per intervals
disp_int_tips<-time.disparity(pco_int_tips, verbose=TRUE)
#disp_int_tips_95axis<-time.disparity(pco_int_tips, verbose=TRUE, rm.last.axis=TRUE)
disp_int_nodes<-time.disparity(pco_int_nodes, verbose=TRUE)
#disp_int_nodes_95axis<-time.disparity(pco_int_nodes, verbose=TRUE, rm.last.axis=TRUE)
disp_int_nodes95<-time.disparity(pco_int_nodes95, verbose=TRUE)
#disp_int_nodes95_95axis<-time.disparity(pco_int_nodes95, verbose=TRUE, rm.last.axis=TRUE)

#Generating the different PCO slices
#methods list
methods=c("random", "acctran", "deltran", "proximity")
#nodes
pco_slices_nodes<-list()
for (type in 1:length(methods)) {
    pco_slices_nodes[[type]]<-slice.pco(pco_data_nodes, tree_nodes, slices, method=methods[[type]], FAD_LAD=FADLAD, verbose=TRUE) 
}
names(pco_slices_nodes)<-methods

#nodes95
pco_slices_nodes95<-list()
for (type in 1:length(methods)) {
    pco_slices_nodes95[[type]]<-slice.pco(pco_data_nodes95, tree_nodes95, slices, method=methods[[type]], FAD_LAD=FADLAD, verbose=TRUE)
}
names(pco_slices_nodes95)<-paste(methods, "95", sep="")

#Calculating the disparity per interval per list
#nodes
disp_slices_nodes<-list()
for (type in 1:length(methods)) {
    disp_slices_nodes[[type]]<-time.disparity(pco_slices_nodes[[type]], verbose=TRUE)
}
names(disp_slices_nodes)<-names(pco_slices_nodes)

#nodes95
disp_slices_nodes95<-list()
for (type in 1:length(methods)) {
    disp_slices_nodes95[[type]]<-time.disparity(pco_slices_nodes95[[type]], verbose=TRUE)
}
names(disp_slices_nodes95)<-names(pco_slices_nodes95)



############################
#Maximum disparity
############################
#Extract the number of states per characters
states.count<-function(character) {
    #Isolate the states
    states<-levels(as.factor(character))
    #Check if multi states
    if(length(grep("&", states)) > 0) {
        #Isolating the multi states
        multi_states<-states[grep("&", states)]
        multi_states<-as.factor(unlist(strsplit(multi_states, split="&")))
        #Removing the multi states from the states list
        states<-states[-grep("&", states)]
        #Check if any of the multi states is not yet present in the states list
        if(any(is.na(match(levels(multi_states), states)))) {
            states<-c(states, multi_states[which(is.na(match(levels(multi_states), states)))])
        }
    }
    #Count the number of states
    return(length(states))
}

states_list<-apply(Nexus_data$matrix, 2, states.count)
max_sp<-prod(states_list) #really high number!

