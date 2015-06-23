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