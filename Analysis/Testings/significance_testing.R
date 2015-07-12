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


#Adding time category to each taxa
mat<-matrix(data=0, nrow=nrow(pco_data), ncol=length(intervals))
rownames(mat)<-rownames(pco_data)
#colnames(mat)<-names(intervals)
colnames(mat)<-c("A","B","C","D","E")
for (int in 1:length(intervals)) {
    mat[c(rownames(intervals[[int]])),int]<-1
}
mat<-as.data.frame(mat)

mat_full<-mat
for (column in 1:ncol(mat_full)) {
    mat_full[,column]<-1 ; mat_full[sample(1:nrow(mat_full), sample(1:nrow(mat_full), 1)), column]<-0
}
#Permanova
adonis(pco_data ~ A+B+C+D+E, data=mat, method='euclidean', permutations=1000)
adonis(pco_data ~ A:B:C:D:E, data=mat, method='euclidean', permutations=1000)

#Calculate the effect of different slices.
#OK

#Add a column for time
mat$time<-NA
for (row in 1:nrow(mat)) {
    mat$time[row]<-paste(colnames(mat)[which(mat[row,]==1)], collapse="")
}

boxplot(pco$points ~ mat$time)
#Permanova 
adonis(pco_data ~ time, data=mat, method='euclidean', permutations=1000)

#Number of overlapping cases = n*(n+1)/2

#Test with the distance matrix PROPER ONE
dist.matrix<-trimmed.max.data$dist.matrix
adonis(dist.matrix ~ time, data=mat, method='euclidean', permutations=1000)

#Test with the full pco ON THE PCO DATA
pco$points 
adonis(pco$points ~ time, data=mat, method='euclidean', permutations=1000)

#Test with time as a strata
adonis(pco$points ~ time, data=mat, strata=mat$time, method='euclidean', permutations=1000)


#GLM
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
print(d.AD <- data.frame(treatment, outcome))
glm.D93 <- glm(counts ~ outcome , family = poisson())
anova(glm.D93)
summary(glm.D93)

#treatments <- distances from centroids
#outcome <- time




#PROPER TESTING


eucl_list<-disp_obs.dist.cent$values
eucl_val<-unlist(eucl_list) ; names(eucl_val)<-NULL
time_counts<-unlist(lapply(eucl_list, lapply, length))
time<-as.factor(rep(names(eucl_list), time_counts))

#Univariate (on the euclidean distance)
bla<-aov(eucl_val~time)
summary(bla)
plot(TukeyHSD(bla))

#Multivariate (on the pco)
pco_mat<-rbind(intervals[[1]], intervals[[2]])
for (int in 3:length(intervals)) {
    pco_mat<-rbind(pco_mat, intervals[[int]])
}
nrow(pco_mat) == length(time)

adonis(pco_mat~time, method='euclidean', permutations=1000)


#add the

#colnames(mat)<-names(intervals)
colnames(mat)<-c("A","B","C","D","E")
for (int in 1:length(intervals)) {
    mat[c(rownames(intervals[[int]])),int]<-1
}
mat<-as.data.frame(mat)


#Do the t-tests


#Calculate manually
squares<-NULL
for (i in 1:nrow(pco_in_bins$pco_intervals[[1]])) {
    squares[i]<-(pco_in_bins$pco_intervals[[1]][i,1]-centroids[[1]][1])^2+(pco_in_bins$pco_intervals[[1]][i,2]-centroids[[1]][2])^2+(pco_in_bins$pco_intervals[[1]][i,3]-centroids[[1]][3])^2
}
SS_A<-sum(squares)/3

mean_SS_A
F.mod_A
R2_A
p_value


combinations<-function(k) {
    for (n in 1:k) {
        x<-combn(1:k, n)
        for(col in 1:ncol(x)) {
            cat(x[,col], sep="") ; cat(" ")
        }
    }
}









#DisparityTTest from Smith et al 

#This function uses the t-distribution to calculate significant differences in disparity among time bins. This script is based on the t-test script of Anderson and Friedman (2012).
#We added additional disparity metrics to their script.
#Output is a symmetric matrix of p-values. Note that the output is not corrected for multiple comparisons so the R function p.adjust() may be required.
#Depending on the number of comparisons you perform, you are likely to also need to do some type of p-value correction (e.g. Bonferroni). This can be achieved using the “p.adjust” function in R. 
#For example, say you were looking at four adjacent time bins (a, b, c, and d) calculated using the “DisparityTTest” function and were interested in the adjusted p-value between them (a vs. b, b vs. c, c vs. d). Say we obtained p-values of 0.035, 0.0012, and 0.081.
#By typing the following line of code  you can correct for multiple comparisons.
#p.adjust(c(0.035, 0.0012, 0.081), method="bonferroni")
#[1] 0.1050 0.0036 0.2430
#Each value here has been adjusted by the Bonferroni correction method. In essence, this method  multiplies the uncorrected  p-value by the number of comparisons (in this case, 3 comparisons were  made). 0.035 was corrected to 0.1050, 0.0012 was corrected to 0.0036, and so on
