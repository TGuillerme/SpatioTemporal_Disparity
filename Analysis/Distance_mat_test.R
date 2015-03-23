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

#After running the Data_setup script, load the different results

######################
#Select the data
######################

#Selecting the file
#SLATER 2013 MEE
#chain_name<-"Slater2013"
#data_path<-"../Data/"
#file_matrix<-"../Data/2013-Slater-MEE-matrix-morpho.nex"
#file_tree<-"../Data/2013-Slater-MEE-TEM.tre"
#int_breaks<-rev(seq(from=0, to=250, by=32.5))
#slices<-rev(seq(from=0, to=250, by=20))+5
#slices[length(slices)]<-0
#KT_bin=5.5
#KT_sli=10

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
tree$root.time<-max(tree.age(tree)[,1])

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

tips_raw_dist<-dist_tips$raw.dist.matrix
tips_GED_dist<-dist_tips$GED.dist.matrix
tips_gow_dist<-dist_tips$gower.dist.matrix
tips_max_dist<-dist_tips$max.dist.matrix
tips_com_dist<-dist_tips$comp.char.matrix

tips_raw_dist<-TrimMorphDistMatrix(tips_raw_dist)$dist.matrix
tips_GED_dist<-TrimMorphDistMatrix(tips_GED_dist)$dist.matrix
tips_gow_dist<-TrimMorphDistMatrix(tips_gow_dist)$dist.matrix
tips_max_dist<-TrimMorphDistMatrix(tips_max_dist)$dist.matrix
tips_com_dist<-TrimMorphDistMatrix(tips_com_dist)$dist.matrix

nodes_raw_dist<-dist_nodes$raw.dist.matrix
nodes_GED_dist<-dist_nodes$GED.dist.matrix
nodes_gow_dist<-dist_nodes$gower.dist.matrix
nodes_max_dist<-dist_nodes$max.dist.matrix
nodes_com_dist<-dist_nodes$comp.char.matrix

nodes_raw_dist<-TrimMorphDistMatrix(nodes_raw_dist)$dist.matrix
nodes_GED_dist<-TrimMorphDistMatrix(nodes_GED_dist)$dist.matrix
nodes_gow_dist<-TrimMorphDistMatrix(nodes_gow_dist)$dist.matrix
nodes_max_dist<-TrimMorphDistMatrix(nodes_max_dist)$dist.matrix
nodes_com_dist<-TrimMorphDistMatrix(nodes_com_dist)$dist.matrix

nodes95_raw_dist<-dist_nodes95$raw.dist.matrix
nodes95_GED_dist<-dist_nodes95$GED.dist.matrix
nodes95_gow_dist<-dist_nodes95$gower.dist.matrix
nodes95_max_dist<-dist_nodes95$max.dist.matrix
nodes95_com_dist<-dist_nodes95$comp.char.matrix

nodes95_raw_dist<-TrimMorphDistMatrix(nodes95_raw_dist)$dist.matrix
nodes95_GED_dist<-TrimMorphDistMatrix(nodes95_GED_dist)$dist.matrix
nodes95_gow_dist<-TrimMorphDistMatrix(nodes95_gow_dist)$dist.matrix
nodes95_max_dist<-TrimMorphDistMatrix(nodes95_max_dist)$dist.matrix
nodes95_com_dist<-TrimMorphDistMatrix(nodes95_com_dist)$dist.matrix

######################
#PCOs
######################

pco_tips_raw_dist<-cmdscale(tips_raw_dist, k=nrow(tips_raw_dist) - 1, add=T)$points
pco_tips_GED_dist<-cmdscale(tips_GED_dist, k=nrow(tips_GED_dist) - 1, add=T)$points
pco_tips_gow_dist<-cmdscale(tips_gow_dist, k=nrow(tips_gow_dist) - 1, add=T)$points
pco_tips_max_dist<-cmdscale(tips_max_dist, k=nrow(tips_max_dist) - 1, add=T)$points
pco_tips_com_dist<-cmdscale(tips_com_dist, k=nrow(tips_com_dist) - 1, add=T)$points

pco_nodes_raw_dist<-cmdscale(nodes_raw_dist, k=nrow(nodes_raw_dist) - 1, add=T)$points
pco_nodes_GED_dist<-cmdscale(nodes_GED_dist, k=nrow(nodes_GED_dist) - 1, add=T)$points
pco_nodes_gow_dist<-cmdscale(nodes_gow_dist, k=nrow(nodes_gow_dist) - 1, add=T)$points
pco_nodes_max_dist<-cmdscale(nodes_max_dist, k=nrow(nodes_max_dist) - 1, add=T)$points
pco_nodes_com_dist<-cmdscale(nodes_com_dist, k=nrow(nodes_com_dist) - 1, add=T)$points

pco_nodes95_raw_dist<-cmdscale(nodes95_raw_dist, k=nrow(nodes95_raw_dist) - 1, add=T)$points
pco_nodes95_GED_dist<-cmdscale(nodes95_GED_dist, k=nrow(nodes95_GED_dist) - 1, add=T)$points
pco_nodes95_gow_dist<-cmdscale(nodes95_gow_dist, k=nrow(nodes95_gow_dist) - 1, add=T)$points
pco_nodes95_max_dist<-cmdscale(nodes95_max_dist, k=nrow(nodes95_max_dist) - 1, add=T)$points
pco_nodes95_com_dist<-cmdscale(nodes95_com_dist, k=nrow(nodes95_com_dist) - 1, add=T)$points

######################
#Disparity
######################

#Calculating the rarefaction
#rarefaction_median<-disparity(pco_data, rarefaction=TRUE, verbose=TRUE, central_tendency=median)

#Generating the different intervals PCOs
int_pco_tips_raw_dist<-int.pco(pco_tips_raw_dist, tree, int_breaks, include.nodes=FALSE, FAD_LAD=FADLAD, diversity=TRUE)
int_div_tips_raw_dist<-int_pco_tips_raw_dist[[2]] ; int_pco_tips_raw_dist<-int_pco_tips_raw_dist[[1]]
int_pco_tips_GED_dist<-int.pco(pco_tips_GED_dist, tree, int_breaks, include.nodes=FALSE, FAD_LAD=FADLAD, diversity=TRUE)
int_div_tips_GED_dist<-int_pco_tips_GED_dist[[2]] ; int_pco_tips_GED_dist<-int_pco_tips_GED_dist[[1]]
int_pco_tips_gow_dist<-int.pco(pco_tips_gow_dist, tree, int_breaks, include.nodes=FALSE, FAD_LAD=FADLAD, diversity=TRUE)
int_div_tips_gow_dist<-int_pco_tips_gow_dist[[2]] ; int_pco_tips_gow_dist<-int_pco_tips_gow_dist[[1]]
int_pco_tips_max_dist<-int.pco(pco_tips_max_dist, tree, int_breaks, include.nodes=FALSE, FAD_LAD=FADLAD, diversity=TRUE)
int_div_tips_max_dist<-int_pco_tips_max_dist[[2]] ; int_pco_tips_max_dist<-int_pco_tips_max_dist[[1]]
int_pco_tips_com_dist<-int.pco(pco_tips_com_dist, tree, int_breaks, include.nodes=FALSE, FAD_LAD=FADLAD, diversity=TRUE)
int_div_tips_com_dist<-int_pco_tips_com_dist[[2]] ; int_pco_tips_com_dist<-int_pco_tips_com_dist[[1]]

int_pco_nodes_raw_dist<-int.pco(pco_nodes_raw_dist, tree, int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD, diversity=TRUE)
int_div_nodes_raw_dist<-int_pco_nodes_raw_dist[[2]] ; int_pco_nodes_raw_dist<-int_pco_nodes_raw_dist[[1]]
int_pco_nodes_GED_dist<-int.pco(pco_nodes_GED_dist, tree, int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD, diversity=TRUE)
int_div_nodes_GED_dist<-int_pco_nodes_GED_dist[[2]] ; int_pco_nodes_GED_dist<-int_pco_nodes_GED_dist[[1]]
int_pco_nodes_gow_dist<-int.pco(pco_nodes_gow_dist, tree, int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD, diversity=TRUE)
int_div_nodes_gow_dist<-int_pco_nodes_gow_dist[[2]] ; int_pco_nodes_gow_dist<-int_pco_nodes_gow_dist[[1]]
int_pco_nodes_max_dist<-int.pco(pco_nodes_max_dist, tree, int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD, diversity=TRUE)
int_div_nodes_max_dist<-int_pco_nodes_max_dist[[2]] ; int_pco_nodes_max_dist<-int_pco_nodes_max_dist[[1]]
int_pco_nodes_com_dist<-int.pco(pco_nodes_com_dist, tree, int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD, diversity=TRUE)
int_div_nodes_com_dist<-int_pco_nodes_com_dist[[2]] ; int_pco_nodes_com_dist<-int_pco_nodes_com_dist[[1]]

int_pco_nodes95_raw_dist<-int.pco(pco_nodes95_raw_dist, tree, int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD, diversity=TRUE)
int_div_nodes95_raw_dist<-int_pco_nodes95_raw_dist[[2]] ; int_pco_nodes95_raw_dist<-int_pco_nodes95_raw_dist[[1]]
int_pco_nodes95_GED_dist<-int.pco(pco_nodes95_GED_dist, tree, int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD, diversity=TRUE)
int_div_nodes95_GED_dist<-int_pco_nodes95_GED_dist[[2]] ; int_pco_nodes95_GED_dist<-int_pco_nodes95_GED_dist[[1]]
int_pco_nodes95_gow_dist<-int.pco(pco_nodes95_gow_dist, tree, int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD, diversity=TRUE)
int_div_nodes95_gow_dist<-int_pco_nodes95_gow_dist[[2]] ; int_pco_nodes95_gow_dist<-int_pco_nodes95_gow_dist[[1]]
int_pco_nodes95_max_dist<-int.pco(pco_nodes95_max_dist, tree, int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD, diversity=TRUE)
int_div_nodes95_max_dist<-int_pco_nodes95_max_dist[[2]] ; int_pco_nodes95_max_dist<-int_pco_nodes95_max_dist[[1]]
int_pco_nodes95_com_dist<-int.pco(pco_nodes95_com_dist, tree, int_breaks, include.nodes=TRUE, FAD_LAD=FADLAD, diversity=TRUE)
int_div_nodes95_com_dist<-int_pco_nodes95_com_dist[[2]] ; int_pco_nodes95_com_dist<-int_pco_nodes95_com_dist[[1]]


#Calculating the disparity per intervals
disp_int_tips_raw<-time.disparity(int_pco_tips_raw_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_int_tips_GED<-time.disparity(int_pco_tips_GED_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_int_tips_gow<-time.disparity(int_pco_tips_gow_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_int_tips_max<-time.disparity(int_pco_tips_max_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_int_tips_com<-time.disparity(int_pco_tips_max_dist, method="centroid",verbose=TRUE, bootstraps=100)

disp_int_nodes_raw<-time.disparity(int_pco_nodes_raw_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_int_nodes_GED<-time.disparity(int_pco_nodes_GED_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_int_nodes_gow<-time.disparity(int_pco_nodes_gow_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_int_nodes_max<-time.disparity(int_pco_nodes_max_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_int_nodes_com<-time.disparity(int_pco_nodes_max_dist, method="centroid",verbose=TRUE, bootstraps=100)

disp_int_nodes95_raw<-time.disparity(int_pco_nodes95_raw_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_int_nodes95_GED<-time.disparity(int_pco_nodes95_GED_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_int_nodes95_gow<-time.disparity(int_pco_nodes95_gow_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_int_nodes95_max<-time.disparity(int_pco_nodes95_max_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_int_nodes95_com<-time.disparity(int_pco_nodes95_max_dist, method="centroid",verbose=TRUE, bootstraps=100)


#Generating the different PCO slices
#methods list
sli_pco_nodes_raw_dist<-slice.pco(pco_nodes_raw_dist, tree, slices, method="random", FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
sli_div_nodes_raw_dist<-sli_pco_nodes_raw_dist[[2]] ; sli_pco_nodes_raw_dist<-sli_pco_nodes_raw_dist[[1]]
sli_pco_nodes_GED_dist<-slice.pco(pco_nodes_GED_dist, tree, slices, method="random", FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
sli_div_nodes_GED_dist<-sli_pco_nodes_GED_dist[[2]] ; sli_pco_nodes_GED_dist<-sli_pco_nodes_GED_dist[[1]]
sli_pco_nodes_gow_dist<-slice.pco(pco_nodes_gow_dist, tree, slices, method="random", FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
sli_div_nodes_gow_dist<-sli_pco_nodes_gow_dist[[2]] ; sli_pco_nodes_gow_dist<-sli_pco_nodes_gow_dist[[1]]
sli_pco_nodes_max_dist<-slice.pco(pco_nodes_max_dist, tree, slices, method="random", FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
sli_div_nodes_max_dist<-sli_pco_nodes_max_dist[[2]] ; sli_pco_nodes_max_dist<-sli_pco_nodes_max_dist[[1]]
sli_pco_nodes_com_dist<-slice.pco(pco_nodes_com_dist, tree, slices, method="random", FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
sli_div_nodes_com_dist<-sli_pco_nodes_com_dist[[2]] ; sli_pco_nodes_com_dist<-sli_pco_nodes_com_dist[[1]]

sli_pco_nodes95_raw_dist<-slice.pco(pco_nodes95_raw_dist, tree, slices, method="random", FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
sli_div_nodes95_raw_dist<-sli_pco_nodes95_raw_dist[[2]] ; sli_pco_nodes95_raw_dist<-sli_pco_nodes95_raw_dist[[1]]
sli_pco_nodes95_GED_dist<-slice.pco(pco_nodes95_GED_dist, tree, slices, method="random", FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
sli_div_nodes95_GED_dist<-sli_pco_nodes95_GED_dist[[2]] ; sli_pco_nodes95_GED_dist<-sli_pco_nodes95_GED_dist[[1]]
sli_pco_nodes95_gow_dist<-slice.pco(pco_nodes95_gow_dist, tree, slices, method="random", FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
sli_div_nodes95_gow_dist<-sli_pco_nodes95_gow_dist[[2]] ; sli_pco_nodes95_gow_dist<-sli_pco_nodes95_gow_dist[[1]]
sli_pco_nodes95_max_dist<-slice.pco(pco_nodes95_max_dist, tree, slices, method="random", FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
sli_div_nodes95_max_dist<-sli_pco_nodes95_max_dist[[2]] ; sli_pco_nodes95_max_dist<-sli_pco_nodes95_max_dist[[1]]
sli_pco_nodes95_com_dist<-slice.pco(pco_nodes95_com_dist, tree, slices, method="random", FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)
sli_div_nodes95_com_dist<-sli_pco_nodes95_com_dist[[2]] ; sli_pco_nodes95_com_dist<-sli_pco_nodes95_com_dist[[1]]

disp_sli_nodes_raw<-time.disparity(sli_pco_nodes_raw_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_sli_nodes_GED<-time.disparity(sli_pco_nodes_GED_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_sli_nodes_gow<-time.disparity(sli_pco_nodes_gow_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_sli_nodes_max<-time.disparity(sli_pco_nodes_max_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_sli_nodes_com<-time.disparity(sli_pco_nodes_max_dist, method="centroid",verbose=TRUE, bootstraps=100)

disp_sli_nodes95_raw<-time.disparity(sli_pco_nodes95_raw_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_sli_nodes95_GED<-time.disparity(sli_pco_nodes95_GED_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_sli_nodes95_gow<-time.disparity(sli_pco_nodes95_gow_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_sli_nodes95_max<-time.disparity(sli_pco_nodes95_max_dist, method="centroid",verbose=TRUE, bootstraps=100)
disp_sli_nodes95_com<-time.disparity(sli_pco_nodes95_max_dist, method="centroid",verbose=TRUE, bootstraps=100)


######################
#Visualisation
######################

quartz(width = 22.4, height = 15.6) #A5 landscape
#Windows dimensions
op<-par(mfrow=c(5, 5), bty="l")# oma=c(bottom, left, top, right)
#Centroid
plot.disparity(disp_int_tips_raw, rarefaction=FALSE, xlab="", ylab="Raw Euclidean Distance (cent)", measure="Cent.dist", main="Intervals: tips", diversity=int_div_tips_raw_dist)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes_raw, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Intervals: nodes", diversity=int_div_nodes_raw_dist)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes95_raw, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Intervals: nodes(95)", diversity=int_div_nodes95_raw_dist)
abline(v= KT_bin, col="red")
plot.disparity(disp_sli_nodes_raw, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Slices: random", diversity=sli_div_nodes_raw_dist)
abline(v= KT_sli, col="red")
plot.disparity(disp_sli_nodes95_raw, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Slices: random (95)", diversity=sli_div_nodes95_raw_dist)
abline(v= KT_sli, col="red")

plot.disparity(disp_int_tips_GED, rarefaction=FALSE, xlab="", ylab="General Euclidean Distance (cent)", measure="Cent.dist", main="Intervals: tips", diversity=int_div_tips_GED_dist)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes_GED, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Intervals: nodes", diversity=int_div_nodes_GED_dist)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes95_GED, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Intervals: nodes(95)", diversity=int_div_nodes95_GED_dist)
abline(v= KT_bin, col="red")
plot.disparity(disp_sli_nodes_GED, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Slices: random", diversity=sli_div_nodes_GED_dist)
abline(v= KT_sli, col="red")
plot.disparity(disp_sli_nodes95_GED, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Slices: random (95)", diversity=sli_div_nodes95_GED_dist)
abline(v= KT_sli, col="red")

plot.disparity(disp_int_tips_gow, rarefaction=FALSE, xlab="", ylab="Corrected raw distance (Gower) (cent)", measure="Cent.dist", main="Intervals: tips", diversity=int_div_tips_gow_dist)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes_gow, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Intervals: nodes", diversity=int_div_nodes_gow_dist)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes95_gow, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Intervals: nodes(95)", diversity=int_div_nodes95_gow_dist)
abline(v= KT_bin, col="red")
plot.disparity(disp_sli_nodes_gow, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Slices: random", diversity=sli_div_nodes_gow_dist)
abline(v= KT_sli, col="red")
plot.disparity(disp_sli_nodes95_gow, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Slices: random (95)", diversity=sli_div_nodes95_gow_dist)
abline(v= KT_sli, col="red")

plot.disparity(disp_int_tips_max, rarefaction=FALSE, xlab="", ylab="Correct raw distance (MOD) (cent)", measure="Cent.dist", main="Intervals: tips", diversity=int_div_tips_max_dist)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes_max, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Intervals: nodes", diversity=int_div_nodes_max_dist)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes95_max, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Intervals: nodes(95)", diversity=int_div_nodes95_max_dist)
abline(v= KT_bin, col="red")
plot.disparity(disp_sli_nodes_max, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Slices: random", diversity=sli_div_nodes_max_dist)
abline(v= KT_sli, col="red")
plot.disparity(disp_sli_nodes95_max, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Slices: random (95)", diversity=sli_div_nodes95_max_dist)
abline(v= KT_sli, col="red")

plot.disparity(disp_int_tips_com, rarefaction=FALSE, xlab="", ylab="Number of common characters (cent)", measure="Cent.dist", main="Intervals: tips", diversity=int_div_tips_com_dist)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes_com, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Intervals: nodes", diversity=int_div_nodes_com_dist)
abline(v= KT_bin, col="red")
plot.disparity(disp_int_nodes95_com, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Intervals: nodes(95)", diversity=int_div_nodes95_com_dist)
abline(v= KT_bin, col="red")
plot.disparity(disp_sli_nodes_com, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Slices: random", diversity=sli_div_nodes_com_dist)
abline(v= KT_sli, col="red")
plot.disparity(disp_sli_nodes95_com, rarefaction=FALSE, xlab="", ylab="", measure="Cent.dist", main="Slices: random (95)", diversity=sli_div_nodes95_com_dist)
abline(v= KT_sli, col="red")

par(op)
