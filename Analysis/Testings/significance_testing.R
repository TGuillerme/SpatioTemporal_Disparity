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
browseURL("https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Functions/disparity/R/MorphDistMatrix.verbose.R")

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
#Cailliez correction applied to negative eigenvalues: D' = -0.5*(D + 1.64469512246817 )^2, except diagonal elements

all(round(cmd$eig, digit=10) == round(pcoa$values$Corr_eig, digit=10))
as.matrix(round(cmd$points[,-13], digit=10)) == as.matrix(round(pcoa$vectors.cor, digit=10)) #weirdly the 8th and 9th dimension are reversed???

#########################
# VOLUME
#########################

#Dummy test
clad.mat<-matrix(data=c(0,0,0,1,0,0,1,1,0,1,1,1), nrow=3) ; rownames(clad.mat) <- letters[1:3]
dist.mat<-as.matrix(dist(clad.mat, diag=TRUE))
pco<-cmdscale(dist.mat, eig=TRUE)
pco.data<-pco$points
pco.eig<-pco$eig
plot(pco.data)

round(pco.eig[-3], digit=10) == round(apply(cov(pco.data),2,sum)*2, digit=10)

#volume.fast(pco.data, pco.eig)
#volume(pco.data)
disparity(pco.data, "volume")

#ok for this simple case

#Dummy test 2
clad.mat<-matrix(data=sample(0:1,1000, replace=TRUE), nrow=20)
dist.mat<-as.matrix(dist(clad.mat, diag=TRUE))
pco<-cmdscale(dist.mat, k=nrow(dist.mat)-1, eig=TRUE, add=F)
pco.data<-pco$points
pco.eig<-pco$eig
plot(pco.data)

round(pco.eig[-20], digit=10) == round(apply(cov(pco.data),2,sum)*19, digit=10)

#volume.fast(pco.data, pco.eig)
#volume(pco.data)
disparity(pco.data, "volume")

#volume increases with the number of characters (increased distances)

# Do the proper Multidimensional Scaling (corrected for negative eigenvalues)
pco <- cmdscale(trimmed.max.data$dist.matrix, k=nrow(trimmed.max.data$dist.matrix) - 2, add=T, eig=TRUE)
pco.data <- pco$points #generates k-2 eigenvectors (=n cladistic space dimensions)
pco.eigen<- pco$eig #generates k eigenvalues (including two null eigenvalues)

#Calculate the hyper-dimensional volume (i.e. size of the cladistic-space)

#volume.fast(pco.data, pco.eigen)
#volume(pco.data)
disparity(pco.data, "volume")

#Relation between eigen value and eigen vectors

set.seed(0)
X<-abs(matrix(rnorm(9),3,3))
diag(X)<-0
X[upper.tri(X)]<-X[lower.tri(X)]
first<-eigen(cov(X))
second<- apply(cov(first$vectors),2, sum)
#eigen(cov(X))$values != cmdscale(X, eig=TRUE)$eig #problem here?

eig.val<-cmdscale(X, eig=TRUE)$eig
eig.vec<-cmdscale(X, eig=TRUE)$points

round(eig.val[-3], digit=10) == round(apply(cov(eig.vec),2,sum)*2, digit=10)

#Works only for a distance matrix and when vectors are reported as in the cmdscale/pcoa algorithm?

# Let's first make a vector of intervals boundaries let's make it very simple: 3 bins with at least three taxa.
bins<-rev(seq(from=0, to=100, by=20))
# Note that the bins must be in reverse chronological order (time since the present)

# Now we can separate the pco.data in different bins
int.pco(pco.data, tree.data, bins)

# We can also make it more accurate by adding the FAD-LAD data
int.pco(pco.data, tree.data, bins, FAD_LAD=ages.data)

# We can also count the diversity
pco_in_bins<-int.pco(pco.data, tree.data, bins, FAD_LAD=ages.data, diversity=TRUE)
# The new object is a list of the binned pco data
pco_in_bins$pco_intervals
# And of the taxonomic counts per bins
pco_in_bins$diversity


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




disparity_full_ran_beck

#values for each slice
beck_65<-disparity_full_ran_beck$values$`65`
beck_60<-disparity_full_ran_beck$values$`60`
beck_55<-disparity_full_ran_beck$values$`55`
beck_50<-disparity_full_ran_beck$values$`50`
beck_45<-disparity_full_ran_beck$values$`45`
beck_40<-disparity_full_ran_beck$values$`40`
beck_35<-disparity_full_ran_beck$values$`35`
beck_30<-disparity_full_ran_beck$values$`30`

dis_pro_max_beck

#quantiles
dis_pro_max_beck<-extract.disp(disparity_full_ran_beck$quantiles, rarefaction="max")

adonis(beck_65~beck_60, permutations=1000, method="euclidean")
adonis(beck_65~beck_55, permutations=1000, method="euclidean")
adonis(beck_65~beck_50, permutations=1000, method="euclidean")
adonis(beck_65~beck_45, permutations=1000, method="euclidean")
adonis(beck_65~beck_40, permutations=1000, method="euclidean")
adonis(beck_65~beck_35, permutations=1000, method="euclidean")
adonis(beck_65~null_ran_centroid_ran[[2]]$values$`65`, permutations=1000, method="euclidean")


adonis(beck_65~beck_60+beck_55+beck_50,beck_45+beck_40+beck_35+beck_30, permutations=1000, method="euclidean")
#NPMANOVA of the PC axes  (e.g. Stayton 2005 and Ruta 2013)
 PC.man <- adonis(PC95axes~sp.fam$Family, data=sp.fam, permutations=999, method="euclidean")

adonis(dis_pro_max_beck[22,3]~dis_pro_max_beck[22,3], permutations=1000)

beck_test<-list(as.vector(beck_60),as.vector(beck_55),as.vector(beck_50),as.vector(beck_45),as.vector(beck_40),as.vector(beck_35),as.vector(beck_30))

bla<-lapply(beck_test, bhatt.coeff, y=as.vector(beck_65))
bhatt.coeff(as.vector(beck_65), as.vector(null_ran_centroid_ran[[2]]$values$`65`))

#Just do a Tukey HSD?