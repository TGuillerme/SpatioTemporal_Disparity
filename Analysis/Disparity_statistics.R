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

#Distances
slat_dist<-load(paste(data_path, chain_name[1], "/", chain_name[1], diversity_ful, sep=""))
beck_dist<-load(paste(data_path, chain_name[2], "/", chain_name[2], diversity_ful, sep=""))





######################################
# Significance testing
######################################

#PERMANOVA
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


