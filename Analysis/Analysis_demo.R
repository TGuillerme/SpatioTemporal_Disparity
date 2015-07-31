# This is a tidied up script going through each step of the disparity analysis on a subset of data for speeding up computational
# time.
# Please email me any questions or query (guillert@tcd.ie).

# This demo is divided in six parts:
# A - Preparing the data 
# B - Estimating the ancestral characters states
# C - Computing the distance matrix
# D - Measuring disparity through time
# E - Visualising disparity through time
# F - Measuring the effect of the K-Pg boundary
# The data necessary for all this analysis is located in the data folder from this repository. If you are running this script from
# the Analysis folder, and that you downloaded all the data you should have no problems running this script. Otherwise it is
# necessary to download both the tree from Beck & Lee 2014 and the the cladistic matrix. Both can be found at the following links
# @@@@@@@@@ ADD LINKS @@@@@@@@@
# The lines you might to modify the path will contain the following tag
    # ~~ PATH ~~

# Before the analysis, we need to download all the functions and packages used in this demo. Most functions are based on Graeme
# Lloyd's Claddis package:
browseURL("http://cran.r-project.org/web/packages/Claddis/")
# and on a draft version of a future package that I will develop along with Graeme Lloyd. For now this package is not properly tested
# and neither is it documented so I will describe the few used functions as we use them. This package, temporarily called 'disparity'
# can be downloaded from my gitHub repository. Because the developement of this package is not finished, there is no manuals
# associated to each function. However, you can find details on each functions in each function script here:
browseURL("https://github.com/TGuillerme/SpatioTemporal_Disparity/tree/master/Functions/disparity/R")

# We can use devtools to install it.
library(devtools)
install_github("TGuillerme/SpatioTemporal_Disparity/Functions/disparity")
library(disparity)

# To speed up the analysis, we are going to use only Beck & Lee (2014; Proceedings B) data.

#########################
#
#   A - Preparing the data
#
#########################

# The tree
Tree_data<-read.nexus("../Data/2014-Beck-ProcB-TEM.tre") # ~~ PATH ~~

# The cladistic matrix
Nexus_data<-ReadMorphNexus("../Data/2014-Beck-ProcB-matrix-morpho.nex") # ~~ PATH ~~
# Isolating the matrix
Nexus_matrix<-Nexus_data$matrix

# We then need to clean both the tree and the matrix using the clean.tree and clean.table function so the data matches exactly the tree.
# Cleaning the tree
Tree_data<-clean.tree(Tree_data, Nexus_matrix)
# Cleaning the matrix
Nexus_data$matrix<-clean.table(Nexus_matrix, Tree_data)
# Making sure they match
nrow(Nexus_data$matrix) == Ntip(Tree_data)

# A last step for the tree is to add the node labels and to make sure the tree has a root time information. We can do that using the
# lapply.root function
tree<-lapply.root(Tree_data)

# We can now visualise the tree using the geoscalePhylo function from the 'strap' package
geoscalePhylo(ladderize(tree), cex.age=0.6, cex.ts=0.7, cex.tip=0.5, units=c("Period","Epoch"), boxes="Epoch")
# We can also add the K-Pg boundary
abline(v=106.3, col="red")


#########################
#
#   B - Estimating the ancestral states
#
#########################

# We can estimate the ancestral states using the anc.state function. This function intakes:
# -A tree
# -A nexus data object
# -A specification on the method used: here we are going to use the ML-claddis method that uses the rerooting method to reconstruct ancestral
# states and Lloyd's method for dealing with missing data (see manuscript for more details).

# Warning, this task can be fairly long because it's estimating the states of 421 characters for 101 nodes. Skip this step if you are not patient
#anc_states<-anc.state(tree, Nexus_data, method='ML-claddis', verbose=TRUE)

# For speeding up this demonstration, we can directly load the previously calculated states.
anc_states<-get(load("../Data/Beck2014/Beck2014_ancestral_states.Rda")) # ~~ PATH ~~

# This object contains both the ancestral states for each nodes added to the cladistic matrix but also the marginal probability associated to
# each state.
tail(anc_states$state[,1:10])
tail(anc_states$prob[,1:10])

# Because the problems linked to estimating ancestral states, we can use the anc.unc function to set a threshold for the marginal likelihood
# under which estimations are replaced by missing data.
save_anc_states<-anc.unc(anc_states, 0.95, missing=NA)$state

# The cladistic matrix now contains missing data
tail(save_anc_states[,1:10])

# We can now create the estimated cladistic matrix (i.e. observed states for tips + estimated states for nodes).
estimated_cladistic_matrix<-Nexus_data
estimated_cladistic_matrix$matrix<-save_anc_states

#########################
#
#   C - Computing the distance matrix
#
#########################

# We can then use the full estimated cladistic matrix to compute the distance matrix between every tips and every nodes. To do so, we can use
# a slightly modified version from the 'MorphDistMatrix' function from 'Claddis' that as been modified to be verbose and is now called
# MorphDistMatrix.verbose"

# Warning this step also takes some times (calculating 203 pairwise distances) so skip this step if you are not patient.
#distance_matrix<-MorphDistMatrix.verbose(estimated_cladistic_matrix, verbose=TRUE)

# For speeding up this demonstration, we can directly load the previously calculated states.
distance_matrix<-get(load("../Data/Beck2014/Beck2014_distance-nodes95.Rda")) # ~~ PATH ~~

# The MorphDistMatrix produces 5 different matrices using 5 different distances. In this manuscript we are just going to use the gower distance
# matrix (i.e. the euclidean distance divided by the number of shared characters).
gower_distance_matrix<-distance_matrix$gower.dist.matrix


#########################
#
#   D - Measuring disparity through time
#
#########################

# Before measuring disparity we need to set up two time variables: the time of the slices (every 5 million years from 170 Mya to the present) and
# getting information on the fossil occurrence (FAD/LAD). The occurrence data has already been compiled in a .csv file.
slices=rev(seq(from=0, to=170, by=5))
FADLAD<-read.csv('../Data/Beck2014_FADLAD.csv', row.names=1) # ~~ PATH ~~

# The first step is to transform the distance matrix into a cladisto-space (i.e. ordination, e.g. pco) using the cmdscale function. Because our
# distance are not Euclidean we are using the Cailliez correction (option 'add=TRUE') and setting the dimensions to the number of tips + nodes -2
cladisto_space<-cmdscale(gower_distance_matrix, k=nrow(gower_distance_matrix) - 2, add=TRUE)$points

# We time-sample this cladisto-space using the 'slice.pco' function. This function needs:
# -The full cladisto-space
# -The tree
# -The slices
# -The slicing method (the evolutionary model)
# -The FAD/LAD data
# It can also take various options such as being verbose (verbose=TRUE) and counting the taxonomic diversity at each sub-sample (diversity=TRUE).
# Here we are going to use the gradual model (method="proximity") but other methods are available (see the function script header for info).
cladisto_space_sliced<-slice.pco(cladisto_space, tree, slices, method='proximity', FAD_LAD=FADLAD, verbose=TRUE, diversity=TRUE)

# We can then isolate the sliced cladisto_space and the taxonomic diversity
taxonomic_diversity<-cladisto_space_sliced$diversity
cladisto_space_sliced<-cladisto_space_sliced$pco_slices

# We can already visualise the taxonomic diversity
plot(taxonomic_diversity, type="l", xlab="sub-samples")


slices_nodes95_div<-pco_slices_nodes95_pro[[2]] ; pco_slices_nodes95_pro<-pco_slices_nodes95_pro[[1]]
#Disparity
disp_sli_nodes95_pro<-time.disparity(pco_slices_nodes95_pro, verbose=TRUE, rarefaction=TRUE, save.all=TRUE)
save(disp_sli_nodes95_pro, file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes95_pro.Rda', sep=''))
#Observed disparity
disp_sli_nodes95_pro_obs<-time.disparity(pco_slices_nodes95_pro, method='centroid', bootstraps=0, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
save(disp_sli_nodes95_pro_obs, file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes95_pro_obs.Rda', sep=''))
#Observed disparity (BS)
disp_sli_nodes95_pro_obs_BS<-time.disparity(pco_slices_nodes95_pro, method='centroid', bootstraps=1000, verbose=TRUE, rarefaction=TRUE, save.all=TRUE, centroid.type='full')
save(disp_sli_nodes95_pro_obs_BS,file=paste(data_path, chain_name, '/',chain_name,'-disp_sli_nodes95_pro_obs_BS.Rda', sep=''))
















#Loading the package
library(devtools)
install_github("TGuillerme/SpatioTemporal_Disparity/Functions/disparity")
library(disparity)

###################################
#
# MAIN FIGURE
#
###################################

######################################
# Isolating the data
######################################

#Setting the variables
#Constant
data_path<-"../Data/"
disaprity_pro<-"-disp_sli_nodes95_pro.Rda"
disparity_ran<-"-disp_sli_nodes95_ran.Rda"
diversity_ful<-"-slices_nodes95_div.Rda"
distance_gowr<-"_distance-nodes95.Rda"

#Data
chain_name<-c("Slater2013", "Beck2014")
file_matrix<-c("../Data/2013-Slater-MEE-matrix-morpho.nex", "../Data/2014-Beck-ProcB-matrix-morpho.nex")
file_tree<-c("../Data/2013-Slater-MEE-TEM.tre","../Data/2014-Beck-ProcB-TEM.tre")

#Loading the data
#Disparity + tree
slat_tmp1<-read.data(chain_name[1], data_path, file_matrix[1], file_tree[1], disaprity_pro)
slat_tmp2<-read.data(chain_name[1], data_path, file_matrix[1], file_tree[1], disparity_ran)
beck_tmp1<-read.data(chain_name[2], data_path, file_matrix[2], file_tree[2], disaprity_pro)
beck_tmp2<-read.data(chain_name[2], data_path, file_matrix[2], file_tree[2], disparity_ran)
#Diversity
slat_div<-load(paste(data_path, chain_name[1], "/", chain_name[1], diversity_ful, sep=""))
beck_div<-load(paste(data_path, chain_name[2], "/", chain_name[2], diversity_ful, sep=""))

#Extracting the data
#Tree
tree_slater<-slat_tmp1[[2]]
tree_beck  <-beck_tmp1[[2]]
#Disparity
disparity_full_pro_slater<-slat_tmp1[[3]]
disparity_full_ran_slater<-slat_tmp2[[3]]
disparity_full_pro_beck  <-beck_tmp1[[3]]
disparity_full_ran_beck  <-beck_tmp2[[3]]

#Diversity
diversity_full_slater<-get(slat_div)
diversity_full_beck  <-get(beck_div)

######################################
# Tree visualisation
######################################

#Plot the tree (slater)
geoscalePhylo(ladderize(tree_slater), cex.age=0.6, cex.ts=0.7, cex.tip=0.5, units=c("Period","Epoch"), boxes="Epoch")
abline(v=240, col="red")
dev.new()
geoscalePhylo(ladderize(tree_beck), cex.age=0.6, cex.ts=0.7, cex.tip=0.5, units=c("Period","Epoch"), boxes="Epoch")
abline(v=106, col="red")
dev.new()

######################################
# Disparity visualisation
######################################
dis_pro_max_beck<-extract.disp(disparity_full_ran_beck$quantiles, rarefaction="max")
dis_pro_max_slater<-extract.disp(disparity_full_ran_slater$quantiles, rarefaction="max")
dis_ran_max_beck<-extract.disp(disparity_full_pro_beck$quantiles, rarefaction="max")
dis_ran_max_slater<-extract.disp(disparity_full_pro_slater$quantiles, rarefaction="max")

#Pretty plot with the tree and the disparity
op<-par(mfrow=c(2,2), bty="n", mar=c(4,4,4,4))
plot.disparity(dis_ran_max_beck, diversity=dis_ran_max_beck$rarefaction, main="Eutherian (punctuated)", xlab="Time (Mya)", y2lab="", ylab="Median distance from centroid")
abline(v=22, col="red")
plot.disparity(dis_ran_max_slater, diversity=dis_ran_max_slater$rarefaction, main="Mammaliaformes (punctuated)", xlab="Time (Mya)", y2lab="Species richness", ylab="")
abline(v=22, col="red")
plot.disparity(dis_pro_max_beck, diversity=dis_pro_max_beck$rarefaction, main="Eutherian (gradual)", xlab="Time (Mya)", y2lab="", ylab="Median distance from centroid")
abline(v=22, col="red")
plot.disparity(dis_pro_max_slater, diversity=dis_pro_max_slater$rarefaction, main="Mammaliaformes (gradual)", xlab="Time (Mya)", y2lab="Species richness", ylab="")
abline(v=22, col="red")
par(op)

######################################
# Rarefaction
######################################
dis_pro_max_beck<-extract.disp(disparity_full_ran_beck$quantiles, rarefaction=8)
dis_pro_max_slater<-extract.disp(disparity_full_ran_slater$quantiles, rarefaction="min")
dis_ran_max_beck<-extract.disp(disparity_full_pro_beck$quantiles, rarefaction=8)
dis_ran_max_slater<-extract.disp(disparity_full_pro_slater$quantiles, rarefaction="min")

#Pretty plot with the tree and the disparity
op<-par(mfrow=c(2,2), bty="n", mar=c(4,4,4,4))
plot.disparity(dis_ran_max_beck, diversity=dis_ran_max_beck$rarefaction, main="Eutherian (punctuated)", xlab="Time (Mya)", y2lab="", ylab="Median distance from centroid")
abline(v=22, col="red")
plot.disparity(dis_ran_max_slater, diversity=dis_ran_max_slater$rarefaction, main="Mammaliaformes (punctuated)", xlab="Time (Mya)", y2lab="Species richness", ylab="")
abline(v=22, col="red")
plot.disparity(dis_pro_max_beck, diversity=dis_pro_max_beck$rarefaction, main="Eutherian (gradual)", xlab="Time (Mya)", y2lab="", ylab="Median distance from centroid")
abline(v=22, col="red")
plot.disparity(dis_pro_max_slater, diversity=dis_pro_max_slater$rarefaction, main="Mammaliaformes (gradual)", xlab="Time (Mya)", y2lab="Species richness", ylab="")
abline(v=22, col="red")
par(op)



















#Loading the package
library(devtools)
install_github("TGuillerme/SpatioTemporal_Disparity/Functions/disparity")
library(disparity)

###################################
#
# TESTING THE DIFFERENCES THROUGH TIME
#
###################################

######################################
# Loading the ordinated data (pco)
######################################

slices<-as.numeric(strsplit(c(noquote('170,165,160,155,150,145,140,135,130,125,120,115,110,105,100,95,90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,5,0')), split=',')[[1]])

#Setting the variables
#Constant
data_path<-"../Data/"
disaprity_pro<-"-disp_sli_nodes95_pro.Rda"
disparity_ran<-"-disp_sli_nodes95_ran.Rda"
diversity_ful<-"-slices_nodes95_div.Rda"
distance_gowr<-"_distance-nodes95.Rda"
FADLADbeck<-'../Data/Beck2014_FADLAD.csv'
FADLADslat<-'../Data/Slater2013_FADLAD.csv'

#FAD/LAD
FADLADbeck<-read.csv(FADLADbeck, row.names=1)
FADLADslat<-read.csv(FADLADslat, row.names=1)

#Data
chain_name<-c("Slater2013", "Beck2014")
file_matrix<-c("../Data/2013-Slater-MEE-matrix-morpho.nex", "../Data/2014-Beck-ProcB-matrix-morpho.nex")
file_tree<-c("../Data/2013-Slater-MEE-TEM.tre","../Data/2014-Beck-ProcB-TEM.tre")

#Loading the data
#Trees
slat_tmp1<-read.data(chain_name[1], data_path, file_matrix[1], file_tree[1], disaprity_pro)
beck_tmp1<-read.data(chain_name[2], data_path, file_matrix[2], file_tree[2], disaprity_pro)
#Distances
slat_dist<-load(paste(data_path, chain_name[1], "/", chain_name[1], distance_gowr, sep=""))
beck_dist<-load(paste(data_path, chain_name[2], "/", chain_name[2], distance_gowr, sep=""))

#Extracting the data
#Trees
#For slater, get the trimmed tree!
tree_slater<-get(load(paste(strsplit(file_tree[1], ".tre")[[1]], ".trimmed_nod", sep="")))
tree_beck  <-beck_tmp1[[2]]
#Distance
#For slater, get the trimmed tree!
distance_gower_slater<-get(load(paste(strsplit(paste(data_path, chain_name[1], "/", chain_name[1], distance_gowr, sep=""), ".Rda")[[1]], ".trimmed", sep="")))
distance_gower_beck  <-get(beck_dist)[[3]]

#Running the pco
pco_slater<-cmdscale(distance_gower_slater, k=nrow(distance_gower_slater) - 2, add=T)$points #-1 or -2? Check manuscript
pco_beck  <-cmdscale(distance_gower_beck  , k=nrow(distance_gower_beck  ) - 2, add=T)$points

#Slicing
pco_slice_slater_pro<-slice.pco(pco_slater, tree_slater, slices, method='proximity', FAD_LAD=FADLADslat, verbose=TRUE)
pco_slice_slater_ran<-slice.pco(pco_slater, tree_slater, slices, method='random'  , FAD_LAD=FADLADslat, verbose=TRUE)
pco_slice_beck_pro  <-slice.pco(pco_beck  , tree_beck  , slices, method='proximity', FAD_LAD=FADLADbeck, verbose=TRUE)
pco_slice_beck_ran  <-slice.pco(pco_beck  , tree_beck  , slices, method='random', FAD_LAD=FADLADbeck, verbose=TRUE)

######################################
# Significance testing
######################################

#PERMANOVA
permanova_pro_slater<-disparity.test.time(pco_slice_slater_pro, method="euclidean", permutations=1000)
permanova_pro_beck  <-disparity.test.time(pco_slice_beck_pro  , method="euclidean", permutations=1000)
permanova_ran_slater<-disparity.test.time(pco_slice_slater_ran, method="euclidean", permutations=10000)
permanova_ran_beck  <-disparity.test.time(pco_slice_beck_ran  , method="euclidean", permutations=1000)

#reference
reftest_pro_slater<-disparity.test(pco_slice_slater_pro[22:35], method="centroid", test="reference", bootstraps=1000)
reftest_pro_beck  <-disparity.test(pco_slice_beck_pro[22:35]  , method="centroid", test="reference", bootstraps=1000)
#reftest_ran_slater<-disparity.test(pco_slice_slater_ran[22:35], method="centroid", test="reference", bootstraps=1000)
reftest_ran_beck  <-disparity.test(pco_slice_beck_ran[22:35]  , method="centroid", test="reference", bootstraps=1000)

#reference (rarefied)
reftestRAR_pro_slater<-disparity.test(pco_slice_slater_pro[22:35], method="centroid", test="reference", bootstraps=1000, rarefaction=8)
reftestRAR_pro_beck  <-disparity.test(pco_slice_beck_pro[22:35]  , method="centroid", test="reference", bootstraps=1000, rarefaction=8)
#reftestRAR_ran_slater<-disparity.test(pco_slice_slater_ran[22:35], method="centroid", test="reference", bootstraps=1000, rarefaction=8)
reftestRAR_ran_beck  <-disparity.test(pco_slice_beck_ran[22:35]  , method="centroid", test="reference", bootstraps=1000, rarefaction=8)


#sequential
#seqtest_pro_slater<-disparity.test(pco_slice_slater_pro, method="centroid", test="sequential", bootstraps=1000)
#seqtest_pro_beck  <-disparity.test(pco_slice_beck_pro  , method="centroid", test="sequential", bootstraps=1000)


#KT test
#KTdata_pro_slater<-int.pco(pco_slater, tree_slater , c(170,65,0), FAD_LAD=FADLADbeck)
#KTtest_pro_slater<-disparity.test(KTdata_pro_slater, method="centroid", test="pairwise", bootstraps=1000)

#KTdata_pro_beck  <-int.pco(pco_beck  , tree_beck  , c(170,65,0), FAD_LAD=FADLADbeck)
#KTtest_pro_beck  <-disparity.test(KTdata_pro_beck, method="centroid", test="pairwise", bootstraps=1000)

######################################
# Making the xtable
######################################
library(xtable)

make.table.permanova<-function() {
    #Fixed rows
    #Empty matrix
    permanova_terms<-as.data.frame(matrix(" ", nrow=8, ncol=3))
    #Column names
    colnames(permanova_terms)<-c("Data", "model", "terms")#, c(colnames(permanova_pro_beck[[1]][1,])))
    #Data
    permanova_terms$Data<-c("Eutherian", rep(" ", 3), "Mammaliformes", rep(" ", 3))
    permanova_terms$model<-rep(c(c("gradual", rep(" ", 1)),c("punctuate", rep(" ", 1))),2)
    permanova_terms$terms<-rep(c("time", "Residuals"), 4)

    #Fixed rows
    #Empty matrix
    permanova_results<-as.data.frame(matrix(NA, nrow=8, ncol=6))
    #Column names
    colnames(permanova_results)<-c(colnames(permanova_pro_beck[[1]][1,]))

    #Variable rows
    permanova_results[1, ]<-as.vector(permanova_pro_beck[[1]][1,])
    permanova_results[2, ]<-as.vector(permanova_pro_beck[[1]][2,])
    permanova_results[3, ]<-as.vector(permanova_ran_beck[[1]][1,])
    permanova_results[4, ]<-as.vector(permanova_ran_beck[[1]][2,])
    permanova_results[5, ]<-as.vector(permanova_pro_slater[[1]][1,])
    permanova_results[6, ]<-as.vector(permanova_pro_slater[[1]][2,])
    permanova_results[7, ]<-as.vector(permanova_ran_slater[[1]][1,])
    permanova_results[8, ]<-as.vector(permanova_ran_slater[[1]][2,])

    #Rounding
    for (n in 1:6) {
        permanova_results[,n]<-round(permanova_results[,n], digit=n)
    }

    return(cbind(permanova_terms, permanova_results))
}

#permanova table
xtable(make.table.permanova(), digits=4)

#lag test table (to modify manually in LaTeX)
xtable(cbind(reftest_pro_beck[[1]][,-5], reftest_ran_beck[[1]][,-5]), digits=4)
#xtable(cbind(reftest_ran_beck[[1]][,-5], reftest_ran_slater[[1]][,-5]))
cat("add: & & gradual & & & & punctuated & & \\\\' in the header.")
cat("replace 'rrrrrrrrr' by 'rrrrr|rrrr'.")

xtable(cbind(reftestRAR_pro_beck[[1]][,-5], reftestRAR_ran_beck[[1]][,-5]), digits=4)
#xtable(cbind(reftest_ran_beck[[1]][,-5], reftest_ran_slater[[1]][,-5]))
cat("add: & & gradual & & & & punctuated & & \\\\' in the header.")
cat("replace 'rrrrrrrrr' by 'rrrrr|rrrr'.")


xtable(cbind(reftest_pro_slater[[1]][,-5], reftestRAR_pro_slater[[1]][,-5]), digits=4)
#xtable(cbind(reftest_ran_beck[[1]][,-5], reftest_ran_slater[[1]][,-5]))
cat("add: & & full data & & & & rarefied & & \\\\' in the header.")
cat("replace 'rrrrrrrrr' by 'rrrrr|rrrr'.")