setwd('~/PhD/Projects/SpatioTemporal_Disparity/Analysis')

#Loading the functions
source("functions.R")
#Loading the data
source("Beck.data.R")
#Loading the first and last apparition datum list
BeckFADLAD<-read.csv("../Data/Beck_FADLAD.csv", row.names=1)
#loading the distance matrix
load("../Data/2014-Beck-dist_matrices2.Rda")

#PCO

#Clean the distance matrix
#Remove the unaplicable characters
trimmed.max.data <-TrimMorphDistMatrix(dist.data$max.dist.matrix)
#Remove the droped taxa from the tree
tree<-drop.tip(tree, trimmed.max.data$removed.taxa)
#Is there any gap in the matrix?
any(is.na(trimmed.max.data$dist.matrix))

#Run the PCO
pco.data <- cmdscale(trimmed.max.data$dist.matrix,
                     k=nrow(trimmed.max.data$dist.matrix) - 1,
                     add=T)$points

#Phylogeny

#Renaming the tree
tree.data<-tree
#Tree ages 
ages.data<-tree.age(tree.data)
#Adding the root age to the tree
tree.data$root.time<-max(ages.data[,1])
#Plot the tree
geoscalePhylo(ladderize(tree.data), cex.age=0.6, cex.ts=0.8, cex.tip=0.5)

#Disparity through time

#Selecting the number of bins
bins_breaks<-rev(hist(ages.data[,1], plot=FALSE)$breaks)+5 ; bins_breaks[10]<-0
#Showing our selected bins
bins_breaks
#Bin the pco data
pco_binned<-bin.pco(pco.data, tree.data, bins_breaks, include.nodes=TRUE, FAD_LAD=BeckFADLAD)
#Calculate disparity
disparity_binned_table<-bin.disparity(pco_binned, method="centroid", verbose=FALSE)
#Plotting the results
plot.disparity(disparity_binned_table, rarefaction=FALSE,
               xlab="bins (Mya)", ylab="Distance from centroid", measure="Cent.dist", las=2)
#Adding the KT boundary
abline(v=c(6.5), col="red")
