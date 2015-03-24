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

library(testthat)

#Disparity null dummy

#Constant evolution tree (step_tree)
step_newick<-"(((((((a:1,b:1):1,c:1):1,d:1):1,e:1):1,f:1):1,g:1):1,h:1);"
step_tree<-read.tree(text=step_newick)
step_tree$node.label<-c("n1","n2","n3","n4","n5","n6","n7")
step_tree$root.time<-7
#plot(step_tree)
#axisPhylo()
#nodelabels(text=step_tree$node.label)

#Character matrix
matrix<-matrix(ncol=80, nrow=15)
matrix[1 ,]<-rep(c(0,0,0,0,0,0,0,0),10)
matrix[2 ,]<-rep(c(0,0,0,0,0,0,0,1),10)
matrix[3 ,]<-rep(c(0,0,0,0,0,0,1,0),10)
matrix[4 ,]<-rep(c(0,0,0,0,0,0,1,1),10)
matrix[5 ,]<-rep(c(0,0,0,0,0,1,1,0),10)
matrix[6 ,]<-rep(c(0,0,0,0,0,1,1,1),10)
matrix[7 ,]<-rep(c(0,0,0,0,1,1,1,0),10)
matrix[8 ,]<-rep(c(0,0,0,0,1,1,1,1),10)
matrix[9 ,]<-rep(c(0,0,0,1,1,1,1,0),10)
matrix[10,]<-rep(c(0,0,0,1,1,1,1,1),10)
matrix[11,]<-rep(c(0,0,1,1,1,1,1,0),10)
matrix[12,]<-rep(c(0,0,1,1,1,1,1,1),10)
matrix[13,]<-rep(c(0,1,1,1,1,1,1,0),10)
matrix[14,]<-rep(c(0,1,1,1,1,1,1,1),10)
matrix[15,]<-rep(c(1,1,1,1,1,1,1,0),10)
rownames(matrix)<-c("n1", "h", "n2", "g", "n3", "f", "n4", "e", "n5", "d", "n6", "c", "n7", "b", "a")

#Proper matrix format
Nexus_data<-list()
Nexus_data$header<-"dummy"
Nexus_data$matrix<-matrix
Nexus_data$ordering<-rep("unord", 80)
Nexus_data$weights<-rep(1, 80)
Nexus_data$max.vals<-rep(1, 80)
Nexus_data$min.vals<-rep(0, 80)

#Distance matrix
dist_data<-MorphDistMatrix.verbose(Nexus_data, verbose=TRUE)

#Distance between n1 and n1 is 0
expect_true(dist_data$max.dist.matrix["n1", "n1"] == 0)
#Distance between n1 and h is equal to distance between n2 and n2
expect_true(dist_data$max.dist.matrix["n1", "h"] == dist_data$max.dist.matrix["n1", "n2"])
#Distance between nodes is increasing
expect_true(dist_data$max.dist.matrix["n1", "n1"] < dist_data$max.dist.matrix["n1", "n2"])
expect_true(dist_data$max.dist.matrix["n1", "n2"] < dist_data$max.dist.matrix["n1", "n3"])
expect_true(dist_data$max.dist.matrix["n1", "n3"] < dist_data$max.dist.matrix["n1", "n4"])
expect_true(dist_data$max.dist.matrix["n1", "n4"] < dist_data$max.dist.matrix["n1", "n5"])
expect_true(dist_data$max.dist.matrix["n1", "n5"] < dist_data$max.dist.matrix["n1", "n6"])
expect_true(dist_data$max.dist.matrix["n1", "n6"] < dist_data$max.dist.matrix["n1", "n7"])

#pco
expect_warning(pco_data<-cmdscale(dist_data$max.dist.matrix, k=nrow(dist_data$max.dist.matrix) - 1, add=T)$points)
#barplot(scree_data<-cumsum(apply(pco_data, 2, var) / sum(apply(pco_data, 2, var))))


#disparity_interval
intervals<-c(6,5,4,3,2,1,0)
pco_intervals<-int.pco(pco_data, step_tree, intervals, include.nodes=TRUE, diversity=TRUE)
div_intervals<-pco_intervals[[2]] ; pco_intervals<-pco_intervals[[1]]
set.seed(0)
disp_intervals<-time.disparity(pco_intervals, verbose=TRUE, method="centroid")

#######################
#compare now to random
#######################

#Random character matrix
rand_mat<-matrix(ncol=80, nrow=15, data=sample(0:1, 80*15, replace=TRUE))
rownames(rand_mat)<-c("n1", "h", "n2", "g", "n3", "f", "n4", "e", "n5", "d", "n6", "c", "n7", "b", "a")

#Proper matrix format
Nexus_rand<-Nexus_data
Nexus_rand$matrix<-rand_mat

#Distance matrix
dist_rand<-MorphDistMatrix.verbose(Nexus_rand, verbose=TRUE)

#pco
pco_rand<-cmdscale(dist_rand$max.dist.matrix, k=nrow(dist_rand$max.dist.matrix) - 1, add=T)$points

#disparity_intervals
intervals<-c(6,5,4,3,2,1,0)
pco_rand_int<-int.pco(pco_rand, step_tree, intervals, include.nodes=TRUE, diversity=TRUE)
div_rand_int<-pco_rand_int[[2]] ; pco_rand_int<-pco_rand_int[[1]]
set.seed(0)
disp_rand_int<-time.disparity(pco_rand_int, verbose=TRUE, method="centroid")


#######################
#compare to simulated
#######################
sim_mat<-matrix(ncol=80, nrow=15, data=NA)
node_state<-attr(tip_state<-sim.character(step_tree, rep(runif(1, 0, 1),2), x0=sample(0:1, 1), model="mk2"), "node.state")
states<-c(node_state, tip_state)
rownames(sim_mat)<-names(states)

for(character in 1:80) {
    node_state<-attr(tip_state<-sim.character(step_tree, rep(runif(1, 0, 1),2), x0=sample(0:1, 1), model="mk2"), "node.state")
    sim_mat[,character]<-c(node_state, tip_state)
}

#Proper matrix format
Nexus_sim<-Nexus_data
Nexus_sim$matrix<-sim_mat

#Distance matrix
dist_sim<-MorphDistMatrix.verbose(Nexus_sim, verbose=TRUE)

#pco
pco_sim<-cmdscale(dist_sim$max.dist.matrix, k=nrow(dist_sim$max.dist.matrix) - 1, add=T)$points

#disparity_intervals
intervals<-c(6,5,4,3,2,1,0)
pco_sim_int<-int.pco(pco_sim, step_tree, intervals, include.nodes=TRUE, diversity=TRUE)
div_pco_sim<-pco_sim_int[[2]] ; pco_sim_int<-pco_sim_int[[1]]
set.seed(0)
disp_sim_int<-time.disparity(pco_sim_int, verbose=TRUE, method="centroid")

dev.new()
#Plot comparison
op<-par(mfrow=c(2, 2), bty="l")# oma=c(bottom, left, top, right)
plot(step_tree)
axisPhylo()
nodelabels(text=step_tree$node.label)
plot.disparity(disp_intervals,diversity=div_intervals, ylim=c(0,1), main="constant data")
plot.disparity(disp_rand_int,diversity=div_intervals, ylim=c(0,1), main="random data")
plot.disparity(disp_sim_int,diversity=div_intervals, ylim=c(0,1), main="simulated data")
par(op)