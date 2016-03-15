# Testing significance

#############################
# Function for loading raw scripts from GitHub 
#############################
# Just loading the whole bunch of functions. I you want to look at them in detail, I've added the GitHub links for each individual function as I call them through this script
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


# Do the proper Multidimensional Scaling (corrected for negative eigenvalues)
pco <- cmdscale(trimmed.max.data$dist.matrix, k=nrow(trimmed.max.data$dist.matrix) - 2, add=T, eig=TRUE)
pco.data <- pco$points #generates k-2 eigenvectors (=n cladistic space dimensions)
pco.eigen<- pco$eig #generates k eigenvalues (including two null eigenvalues)

bins<-rev(seq(from=0, to=100, by=20))
pco_in_bins<-int.pco(pco.data, tree.data, bins, FAD_LAD=ages.data, diversity=TRUE)


#############################
# Testing significance
############################# 

#Creating three different disparities

disp_med.dist.cent.BS<-time.disparity(pco_in_bins$pco_intervals, method="centroid", relative=FALSE, verbose=TRUE, bootstraps=10, rarefaction=TRUE, centroid.type="median", central_tendency=mean, save.all=TRUE)
disp_obs.dist.cent<-time.disparity(pco_in_bins$pco_intervals, method="centroid", relative=FALSE, verbose=TRUE, bootstraps=0,  rarefaction=FALSE, centroid.type="full", central_tendency=median, save.all=TRUE)
disp_obs.dist.cent.BS<-time.disparity(pco_in_bins$pco_intervals, method="centroid", relative=FALSE, verbose=TRUE, bootstraps=10, rarefaction=TRUE, centroid.type="full", central_tendency=median, save.all=TRUE)

#Extract disparity from the rarefaction analysis
disp_med.dist.cent.BS_max<-extract.disp(disp_med.dist.cent.BS$quantiles, rarefaction="max")

op<-par(mfrow=c(1,3))
plot.disparity(extract.disp(disp_med.dist.cent.BS$quantiles, rarefaction="max"), ylim=c(0.8, 1.8))
plot.disparity(disp_obs.dist.cent$quantiles, ylim=c(0.8, 1.8))
plot.disparity(extract.disp(disp_obs.dist.cent.BS$quantiles, rarefaction="max"), ylim=c(0.8, 1.8))
par(op)


#ADONIS testing
disparity.test.time(pco_in_bins[[1]], method="euclidean", permutations=1000)

#If significant (or not) test differences between slices
#sequential (changes through time)
sequential_diff<-disparity.test(pco_in_bins[[1]], method="centroid", test="sequential", bootstraps=1000)

#reference (lag effect)
reference_diff<-disparity.test(pco_in_bins[[1]], method="centroid", test="reference", bootstraps=1000)

#Test before/after K-T
KT_test<-int.pco(pco.data, tree.data, intervals=c(100,75,0), FAD_LAD=ages.data, diversity=TRUE)
disparity.test.time(KT_test[[1]], method="euclidean", permutations=1000)
pair_test<-disparity.test(KT_test[[1]], method="centroid", test="pairwise", bootstraps=1000)




#############################
#Visualising the results
#############################

#Full cladisto-space
pco_data<-pco$points #For code later: can be extracted from intervals (unique) 
pco_data<-pco$points[,1:3] #3D version (easier to calculate)
intervals<-pco_in_bins$pco_intervals

#Calculating the overall centroid distances
centroid_overall<-apply(pco_data, 2, mean)

#Calculating the centroids per interval

centroids<-list()
cent.dist<-list()

for (int in 1:length(intervals)) {
    
    #Isolating one interval
    X<-pco_in_bins$pco_intervals[[int]][,1:3] #3D only!

    #Calculating it's centroid
    centroids[[int]]<-apply(X, 2, mean)
    
    #Calculating all the distances from this centroid
    Y<-NULL
    for (j in 1:nrow(X)){
        Y[j] <- dist(rbind(X[j,], centroids[[int]]), method="euclidean")
    }
    cent.dist[[int]]<-Y

}

#Calculating distances from time centroids to overall centroid
cent_to_overall<-NULL
cent_mat<-matrix(data=unlist(centroids), ncol=3, nrow=length(intervals), byrow=TRUE)
for (j in 1:nrow(cent_mat)){
    cent_to_overall[j] <- dist(rbind(cent_mat[j,], centroid_overall), method="euclidean")
}


#Visualizing the results in 3D
library(scatterplot3d)

for (int in 1:length(intervals)) {
    #Plot each observation
    s3d<-scatterplot3d(pco_data)
    #Add the centroid
    s3d$points3d(centroids[[int]][1],centroids[[int]][2],centroids[[int]][3], col="blue",pch=16)
    #Add the distances within centroids
    s3d$points3d(c(rbind(centroids[[int]][1],pco_in_bins$pco_intervals[[int]][,1])), c(rbind(centroids[[int]][2],pco_in_bins$pco_intervals[[int]][,2])), c(rbind(centroids[[int]][3],pco_in_bins$pco_intervals[[int]][,3])), type="l", lty=3, col="grey")
    #Pause
    Sys.sleep(1)
}

#Plot each observation
s3d<-scatterplot3d(pco_data)

#Add the centroid
for (int in 1:length(intervals)) {
    s3d$points3d(centroids[[int]][1],centroids[[int]][2],centroids[[int]][3], col="blue",pch=16)
}

#Add the overall centroid
s3d$points3d(centroid_overall[1],centroid_overall[2],centroid_overall[3], col="red",pch=16)
#Add the distances between centroids
s3d$points3d(c(rbind(centroid_overall[1],cent_mat[,1])), c(rbind(centroid_overall[2],cent_mat[,2])), c(rbind(centroid_overall[3],cent_mat[,3])), type="l", lty=1, col="red")

