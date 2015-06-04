library(devtools)
install_github("TGuillerme/SpatioTemporal_Disparity/Functions/disparity")
library(disparity)

#############################
# Data loading
#############################
# This is just exactly your workshop script up until the distance matrix calculations (without the "install.packages")

# Load the morphological matrix
nexus.data <- ReadMorphNexus("http://www.graemetlloyd.com/nexus/Cullen_etal_2013a.nex")

# Load the tree
tree.data <- read.tree("http://www.graemetlloyd.com/firstmpt/Cullen_etal_2013a.tre")

# Fix polytomies.
tree.data <- multi2di(tree.data)

# Read the age data
ages.data <- read.table("http://www.graemetlloyd.com/teaching/RE2014/Cullenages.txt", row.names=1, sep="\t", header=T)

# Building the tree
tree.data <- timePaleoPhy(tree.data, ages.data, type="mbl", vartime=2)

# Safe taxonomic reduction
safe.data <- SafeTaxonomicReduction(nexus.data)

#############################
# Distance matrix
#############################

# Distance matrix
dist.data <- MorphDistMatrix(nexus.data)

# For small datasets this is instantly but I'm playing around with reasonably big ones (100 taxa * 400 characters) and some can take a day or so to compute.
# To make sure the script is not frozen, I've just added a verbose function.
# I clumsily renamed it MorphDistMatrix.verbose to avoid any confusion.
#browseURL("https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Functions/disparity/R/MorphDistMatrix.verbose.R")

# The function is exactly the same but can take an optional additional argument: verbose (TRUE (default) or FALSE)
# With verbose = FALSE, the function is exactly the same
dist.data.verbose <- MorphDistMatrix.verbose(nexus.data, verbose=FALSE)

#The verbose option just adds some dots and some text (probably worth modifying but that's the general idea)
dist.data.verbose <- MorphDistMatrix.verbose(nexus.data, verbose=TRUE)

# Just for the sake of testing the differences
for (different_distances in 1:length(dist.data)) {
    print( all( (dist.data.verbose[[different_distances]] == dist.data[[different_distances]]), na.rm=TRUE) )
}

# Let's just continue through your workshop to get to the disparity part

# Trim inapplicable data
trimmed.max.data <-TrimMorphDistMatrix(dist.data$max.dist.matrix)

# Remove the trimmed taxa from the tree for later
tree.data<-drop.tip(tree.data, trimmed.max.data$removed.taxa)

# PCO on MOD distance
pco <- cmdscale(trimmed.max.data$dist.matrix, k=nrow(trimmed.max.data$dist.matrix) - 1, add=T, eig=TRUE)

#########################
# ABOUT THE PC
#########################

# Difference between pcoa and cmdscale??
cmd <- cmdscale(trimmed.max.data$dist.matrix, k=9, add=F, eig=TRUE)
pcoa<- pcoa(trimmed.max.data$dist.matrix)

round(cmd$eig, digit=10) == round(pcoa$values$Eigenvalues, digit=10)
as.matrix(round(cmd$points, digit=10)) == as.matrix(round(pcoa$vectors, digit=10)) #weirdly the 7th dimension is reversed???

#In cmdscale the add=T uses a Cailliez (1983) correction correcting for negative eigenvalues. 
#To maximise the variance in our analysis, we create a space of n dimensions where n cam be represented by exactly in at most k-1 dimensions. Here we use k-2 dimensions because the two last eigenvalues are null (following Cailliez correction) 

cmd <- cmdscale(trimmed.max.data$dist.matrix, k=nrow(trimmed.max.data$dist.matrix) - 2, add=T, eig=TRUE)
pcoa<- pcoa(trimmed.max.data$dist.matrix, correction="cailliez")

X<-cmd$points

#centroid.apply<-function(X) {
    #Euclidean distances to the centroid
    #This function is based on euc.dist.cent() from Finlay & Cooper 2015 - PeerJ (https://github.com/SiveFinlay/Diversity_Paper/blob/master/functions/Morpho_diversity_functions.r) 
    #Centroid (mean score of each PC axis)
    centroid<-apply(X, 2, mean)
        #Outputs length(X) values that are the mean value of each eigenvector (X)

    #Euclidean distances to the centroid
    cent.dist<-NULL
    for (j in 1:nrow(X)){
        cent.dist[j] <- dist(rbind(X[j,], centroid), method="euclidean")
            #Outputs 1 value that is the euclidean distance between the scores of taxa j for each eigenvector and the mean of each eigenvector? 
    }
    #return(cent.dist)
#}

#centroid.writen.apply<-function(X) {
    #Setting K (number of taxa) and N (number of eigenvectors)
    K<-nrow(X)
    N<-ncol(X)
    #Calculate the centroid for each eigen vector as sum(eigenvector)/number of taxa (k)
    centroid.w<-NULL
    for(n in 1:N) {
        centroid.w[n]<-sum(X[,n])/K
    }
    # Is approximatively the same than centroid
    round(centroid.w, digit=15) == round(centroid, digit=15) #note that it's equal to 0!!!!!

    #Eucliden distances to centroid
    cent.dist.w<-X
    for (n in 1:N) {
        cent.dist.w[,n]<-sqrt((X[,n]-centroid[n])^2)
    }
        #This outputs the distance from each species to the centroid of each eigenvector
#}

#centroid.n.apply<-function(X) {
    #Setting K (number of taxa) and N (number of eigenvectors)
    centroid<-apply(X, 2, mean)

    #Euclidean distance to centroid
    cent.dist.n<-NULL
    for(i in 1:K) {
        cent.dist.n[i]<-sqrt(sum((X[i,]-centroid)^2))
    }

    #And this is exactly the same
    round(cent.dist, digit=15)==round(cent.dist.n, digit=15)
#}

#VISUAL (in 3D)
install.packages("scatterplot3d")
library(scatterplot3d)

cmd <- cmdscale(trimmed.max.data$dist.matrix, k=3, add=T, eig=TRUE)

X<-cmd$points
centroid<-apply(X, 2, mean)


op<-par(mfrow=c(2,1), mar=c(0,0,0,0))
nf<-layout(matrix(c(1,2),2,1,byrow = TRUE), c(1,2),c(2,1), FALSE)
#layout.show(nf)
#Plot each observation (ordinated)
s3d<-scatterplot3d(X)
#Add the centroid
s3d$points3d(centroid[1],centroid[2],centroid[3], col="red",pch=16)
#Add the distances
s3d$points3d(c(rbind(centroid[1],X[,1])), c(rbind(centroid[2],X[,2])), c(rbind(centroid[3],X[,3])), type="l", lty=3)

#Plot the distances to centroid distribution
med.cent<-round(median(cent.dist), digit=3)
barplot(cent.dist, main=paste("median distance =", med.cent), width=0.5)

par(op)