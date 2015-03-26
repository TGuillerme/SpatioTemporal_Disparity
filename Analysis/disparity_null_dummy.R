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

##########################
#Disparity null dummy
##########################

#Constant evolution tree (step_tree)
step_newick<-"(((((((a:1,b:1):1,c:1):1,d:1):1,e:1):1,f:1):1,g:1):1,h:1);"
step_tree<-read.tree(text=step_newick)
step_tree$node.label<-c("n1","n2","n3","n4","n5","n6","n7")
step_tree$root.time<-7


#Character matrix
matrix<-matrix(ncol=80, nrow=15)
matrix[1 ,]<-rep(c(1,1,1,1,1,1,1,0),10) #a
matrix[2 ,]<-rep(c(0,1,1,1,1,1,1,1),10) #b
matrix[3 ,]<-rep(c(0,0,1,1,1,1,1,1),10) #c
matrix[4 ,]<-rep(c(0,0,0,1,1,1,1,1),10) #d
matrix[5 ,]<-rep(c(0,0,0,0,1,1,1,1),10) #e
matrix[6 ,]<-rep(c(0,0,0,0,0,1,1,1),10) #f
matrix[7 ,]<-rep(c(0,0,0,0,0,1,1,1),10) #g
matrix[8 ,]<-rep(c(0,0,0,0,0,0,0,1),10) #h
matrix[9 ,]<-rep(c(0,0,0,0,0,0,0,0),10) #n1
matrix[10,]<-rep(c(0,0,0,0,0,0,1,0),10) #n2
matrix[11,]<-rep(c(0,0,0,0,0,1,1,0),10) #n3
matrix[12,]<-rep(c(0,0,0,0,1,1,1,0),10) #n4
matrix[13,]<-rep(c(0,0,0,1,1,1,1,0),10) #n5
matrix[14,]<-rep(c(0,0,1,1,1,1,1,0),10) #n6
matrix[15,]<-rep(c(0,1,1,1,1,1,1,0),10) #n7
rownames(matrix)<-c(step_tree$tip.label, step_tree$node.label)
#Proper matrix format
Nexus_data<-make.nexus(matrix)
#Distance matrix
dist_data<-MorphDistMatrix.verbose(Nexus_data, verbose=TRUE)
#pco
pco_data<-cmdscale(dist_data$max.dist.matrix, k=nrow(dist_data$max.dist.matrix) - 1, add=T)$points
#disparity_interval
intervals<-c(6,5,4,3,2,1,0)
pco_intervals<-int.pco(pco_data, step_tree, intervals, include.nodes=TRUE, diversity=TRUE)
div_intervals<-pco_intervals[[2]] ; pco_intervals<-pco_intervals[[1]]
disp_intervals<-time.disparity(pco_intervals, verbose=TRUE, method="centroid", save.all=TRUE, bootstraps=10000)


#Random character matrix
rand_data<-null.data(tree=step_tree, matrix=rep(2,80), matrix.model="random", replicates=100, verbose=TRUE, include.nodes=TRUE)
ran_matrix<-lapply(rand_data, make.nexus)
#distance
dist_ran<-lapply(ran_matrix, MorphDistMatrix.verbose, verbose=TRUE)
dist_ran<-lapply(dist_ran, extract.dist, distance="max.dist.matrix")
#pco
pco_ran<-lapply(dist_ran, function(X) cmdscale(X, k=nrow(X)-1, add=TRUE)$points)
#intervals
pco_ran_int<-lapply(pco_ran, function(X) int.pco(X, step_tree, intervals, include.nodes=TRUE))
#diversity
div_ran_int<-int.pco(pco_ran[[1]], step_tree, intervals, include.nodes=TRUE, diversity=TRUE)[[2]]
#disparity
disp_ran_int<-lapply(pco_ran_int, time.disparity, verbose=TRUE, method="centroid", save.all=TRUE)
#Combine results
disp_ran_int<-combine.disp(disp_ran_int)


#Simulated data
sim_data<-null.data(tree=step_tree, matrix=rep(2,80), matrix.model="sim.char", replicates=100, verbose=TRUE, include.nodes=TRUE)
sim_matrix<-lapply(sim_data, make.nexus)
#distance
dist_sim<-lapply(sim_matrix, MorphDistMatrix.verbose, verbose=TRUE)
dist_sim<-lapply(dist_sim, extract.dist, distance="max.dist.matrix")
#pco
pco_sim<-lapply(dist_sim, function(X) cmdscale(X, k=nrow(X)-1, add=TRUE)$points)
#intervals
pco_sim_int<-lapply(pco_sim, function(X) int.pco(X, step_tree, intervals, include.nodes=TRUE))
#diversity
div_sim_int<-int.pco(pco_sim[[1]], step_tree, intervals, include.nodes=TRUE, diversity=TRUE)[[2]]
#disparity
disp_sim_int<-lapply(pco_sim_int, time.disparity, verbose=TRUE, method="centroid", save.all=TRUE)
#Combine results
disp_sim_int<-combine.disp(disp_sim_int)

#differences
ran_sim<-pair.bhatt.coeff(lapply(disp_ran_int[[2]], as.vector), lapply(disp_sim_int[[2]], as.vector))
ran_obs<-pair.bhatt.coeff(lapply(disp_ran_int[[2]], as.vector), lapply(disp_intervals[[2]], as.vector))
obs_sim<-pair.bhatt.coeff(lapply(disp_intervals[[2]], as.vector), lapply(disp_sim_int[[2]], as.vector))

#Plot comparison
op<-par(mfrow=c(2, 3), bty="l")
plot(step_tree, main="Step tree 80 characters")
axisPhylo()
nodelabels(text=step_tree$node.label, cex=0.7)
plot.disparity(disp_intervals[[1]],diversity=div_intervals, ylim=c(0,1), main="constant data")
plot.disparity(disp_ran_int[[1]] ,diversity=div_intervals, ylim=c(0,1), main="random data")
plot.disparity(disp_sim_int[[1]]  ,diversity=div_intervals, ylim=c(0,1), main="simulated data")
plot(ran_sim, type="l", ylim=c(0,1), col="black", ylab="Similarity(BC)", xlab="") #Random vs simulated
legend(3, 0.3, legend=c("ran vs sim", "obs vs ran", "obs vs sim"), col=c("black", "red", "blue"), lty=1, cex=0.8)
abline(h=0.95, col="grey", lty=3)
abline(h=0.05, col="grey", lty=3)
#axis(side = 1, at=1:length(names(div_intervals)), labels=names(div_intervals), las=2)
points(ran_obs, type="l", col="red") #Observed vs random
points(obs_sim, type="l", col="blue") #Observed vs random
par(op)


#############################
#Constant disparity increase
#############################

#Fully balanced tree with constant branch length
bal_tree<-stree(16, type="balanced")
bal_tree$edge.length<-rep(1, Ntip(bal_tree)+Nnode(bal_tree))
bal_tree<-lapply.root(bal_tree, max(tree.age(bal_tree)$ages))

#Max parsimony matrix
matrix<-matrix(ncol=110, nrow=31)
matrix[1 ,]<-rep(c(0,0,0,0,0,0,0,0,0,0,0),10) #n1
matrix[2 ,]<-rep(c(0,0,0,0,1,0,0,0,0,0,0),10) #n2
matrix[3 ,]<-rep(c(0,0,0,0,1,0,0,1,0,0,0),10) #n3
matrix[4 ,]<-rep(c(0,0,0,0,1,0,1,1,0,0,0),10) #n4
matrix[5 ,]<-rep(c(0,1,0,0,1,0,0,1,0,0,0),10) #n5
matrix[6 ,]<-rep(c(0,0,1,0,1,0,0,0,0,0,0),10) #n6
matrix[7 ,]<-rep(c(0,0,1,0,1,1,0,0,0,0,0),10) #n7
matrix[8 ,]<-rep(c(0,1,1,0,1,0,0,0,0,0,0),10) #n8
matrix[9 ,]<-rep(c(0,0,0,1,0,0,0,0,0,0,0),10) #n9
matrix[10,]<-rep(c(0,0,0,1,0,0,0,0,0,1,0),10) #n10
matrix[11,]<-rep(c(0,0,0,1,0,0,0,0,1,1,0),10) #n11
matrix[12,]<-rep(c(0,1,0,1,0,0,0,0,0,1,0),10) #n12
matrix[13,]<-rep(c(0,0,1,1,0,0,0,0,0,0,0),10) #n13
matrix[14,]<-rep(c(0,0,1,1,0,0,0,0,0,0,1),10) #n14
matrix[15,]<-rep(c(0,1,1,1,0,0,0,0,0,0,0),10) #n15

matrix[16,]<-rep(c(0,0,0,0,1,0,1,1,0,0,0),10) #t1
matrix[17,]<-rep(c(1,0,0,0,1,0,1,1,0,0,0),10) #t2
matrix[18,]<-rep(c(0,1,0,0,1,0,0,1,0,0,0),10) #t3
matrix[19,]<-rep(c(1,1,0,0,1,0,0,1,0,0,0),10) #t4
matrix[20,]<-rep(c(0,0,1,0,1,1,0,0,0,0,0),10) #t5
matrix[21,]<-rep(c(1,0,1,0,1,1,0,0,0,0,0),10) #t6
matrix[22,]<-rep(c(0,0,1,0,1,0,0,0,0,0,0),10) #t7
matrix[23,]<-rep(c(1,0,1,0,1,0,0,0,0,0,0),10) #t8
matrix[24,]<-rep(c(0,0,0,1,0,0,0,0,1,1,0),10) #t9
matrix[25,]<-rep(c(1,0,0,1,0,0,0,0,1,1,0),10) #t10
matrix[26,]<-rep(c(0,1,0,1,0,0,0,0,0,1,0),10) #t11
matrix[27,]<-rep(c(1,1,0,1,0,0,0,0,0,1,0),10) #t12
matrix[28,]<-rep(c(0,0,1,1,0,0,0,0,0,0,1),10) #t13
matrix[29,]<-rep(c(1,0,1,1,0,0,0,0,0,0,1),10) #t14
matrix[30,]<-rep(c(0,1,1,1,0,0,0,0,0,0,0),10) #t15
matrix[31,]<-rep(c(1,1,1,1,0,0,0,0,0,0,0),10) #t16
rownames(matrix)<-c(bal_tree$node.label, bal_tree$tip.label)

Nexus_data<-make.nexus(matrix)

#Distance matrix
dist_data<-MorphDistMatrix.verbose(Nexus_data, verbose=TRUE)

#pco
pco_data<-cmdscale(dist_data$max.dist.matrix, k=nrow(dist_data$max.dist.matrix) - 1, add=T)$points

#disparity_interval
intervals<-c(4,3.5,3,2.5,2,1.5,1,0.5,0)
pco_intervals<-int.pco(pco_data, bal_tree, intervals, diversity=TRUE, include.nodes=TRUE)
div_intervals<-pco_intervals[[2]] ; pco_intervals<-pco_intervals[[1]]
disp_intervals<-time.disparity(pco_intervals, verbose=TRUE, method="centroid", save.all=TRUE)

#Random character matrix
ran_data<-null.data(tree=bal_tree, matrix=rep(2,110), matrix.model="random", replicates=100, verbose=TRUE, include.nodes=TRUE)
ran_matrix<-lapply(ran_data, make.nexus)
#distance
dist_ran<-lapply(ran_matrix, MorphDistMatrix.verbose, verbose=TRUE)
dist_ran<-lapply(dist_ran, extract.dist, distance="max.dist.matrix")
#pco
pco_ran<-lapply(dist_ran, function(X) cmdscale(X, k=nrow(X)-1, add=TRUE)$points)
#intervals
pco_ran_int<-lapply(pco_ran, function(X) int.pco(X, bal_tree, intervals, include.nodes=TRUE))
#diversity
div_ran_int<-int.pco(pco_ran[[1]], bal_tree, intervals, include.nodes=TRUE, diversity=TRUE)[[2]]
#disparity
disp_ran_int<-lapply(pco_ran_int, time.disparity, verbose=TRUE, method="centroid", save.all=TRUE)
#Combine results
disp_ran_int<-combine.disp(disp_ran_int)

#Random character matrix
sim_data<-null.data(tree=bal_tree, matrix=rep(2,110), matrix.model="sim.char", replicates=100, verbose=TRUE, include.nodes=TRUE)
sim_matrix<-lapply(sim_data, make.nexus)
#distance
dist_sim<-lapply(sim_matrix, MorphDistMatrix.verbose, verbose=TRUE)
dist_sim<-lapply(dist_sim, extract.dist, distance="max.dist.matrix")
#pco
pco_sim<-lapply(dist_sim, function(X) cmdscale(X, k=nrow(X)-1, add=TRUE)$points)
#intervals
pco_sim_int<-lapply(pco_sim, function(X) int.pco(X, bal_tree, intervals, include.nodes=TRUE))
#diversity
div_sim_int<-int.pco(pco_sim[[1]], bal_tree, intervals, include.nodes=TRUE, diversity=TRUE)[[2]]
#disparity
disp_sim_int<-lapply(pco_sim_int, time.disparity, verbose=TRUE, method="centroid", save.all=TRUE)
#Combine results
disp_sim_int<-combine.disp(disp_sim_int)

#differences
ran_sim<-pair.bhatt.coeff(lapply(disp_ran_int[[2]], as.vector), lapply(disp_sim_int[[2]], as.vector))
ran_obs<-pair.bhatt.coeff(lapply(disp_ran_int[[2]], as.vector), lapply(disp_intervals[[2]], as.vector))
obs_sim<-pair.bhatt.coeff(lapply(disp_intervals[[2]], as.vector), lapply(disp_sim_int[[2]], as.vector))


#Plot comparison
op<-par(mfrow=c(2, 3), bty="l")
plot(bal_tree, main="Balanced tree 110 characters")
axisPhylo()
nodelabels(text=bal_tree$node.label, cex=0.7)
plot.disparity(disp_intervals[[1]],diversity=log(cor.diversity(div_intervals)), ylim=c(0,1), main="constant data")
plot.disparity(disp_ran_int[[1]]  ,diversity=log(cor.diversity(div_ran_int)), ylim=c(0,1), main="random data")
plot.disparity(disp_sim_int[[1]]  ,diversity=log(cor.diversity(div_sim_int)), ylim=c(0,1), main="simulated data")
plot(ran_sim, type="l", ylim=c(0,1), col="black", ylab="Similarity(BC)", xlab="") #Random vs simulated
legend(4, 0.4, legend=c("ran vs sim", "obs vs ran", "obs vs sim"), col=c("black", "red", "blue"), lty=1, cex=0.8)
abline(h=0.95, col="grey", lty=3)
abline(h=0.05, col="grey", lty=3)
axis(side = 1, at=1:length(names(cor.diversity(div_intervals))), labels=names(cor.diversity(div_intervals)), las=2)
points(ran_obs, type="l", col="red") #Observed vs random
points(obs_sim, type="l", col="blue") #Observed vs random
par(op)