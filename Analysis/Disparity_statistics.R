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
tree_slater<-slat_tmp1[[2]]
tree_beck  <-beck_tmp1[[2]]
#Distance
distance_gower_slater<-get(slat_dist)[[3]] #BUG WITH SLATER TREE!
distance_gower_beck  <-get(beck_dist)[[3]]

#Running the pco
pco_slater<-cmdscale(distance_gower_slater, k=nrow(distance_gower_slater) - 2, add=T)$points #-1 or -2? Check manuscript
pco_beck  <-cmdscale(distance_gower_beck  , k=nrow(distance_gower_beck  ) - 2, add=T)$points

#Slicing
pco_slice_slater_pro<-slice.pco(pco_slater, tree_slater, slices, method='proximity', FAD_LAD=FADLADslat, verbose=TRUE)
pco_slice_beck_pro  <-slice.pco(pco_beck  , tree_beck  , slices, method='proximity', FAD_LAD=FADLADbeck, verbose=TRUE)

######################################
# Significance testing
######################################

#PERMANOVA
permanova_pro_slater<-disparity.test.time(pco_slice_slater_pro, method="euclidean", permutations=1000)
permanova_pro_beck  <-disparity.test.time(pco_slice_beck_pro  , method="euclidean", permutations=1000)

#T-Test
#sequential
seqtest_pro_slater<-disparity.test(pco_slice_slater_pro, method="centroid", test="sequential", bootstraps=1000)
seqtest_pro_beck  <-disparity.test(pco_slice_beck_pro  , method="centroid", test="sequential", bootstraps=1000)

#reference
reftest_pro_slater<-disparity.test(pco_slice_slater_pro[22:35], method="centroid", test="reference", bootstraps=1000)
reftest_pro_beck  <-disparity.test(pco_slice_beck_pro[22:35]  , method="centroid", test="reference", bootstraps=1000)

#KT test
KTdata_pro_slater<-int.pco(pco_slater, tree_slater , c(170,65,0), FAD_LAD=FADLADbeck)
KTtest_pro_slater<-disparity.test(KTdata_pro_slater, method="centroid", test="pairwise", bootstraps=1000)

KTdata_pro_beck  <-int.pco(pco_beck  , tree_beck  , c(170,65,0), FAD_LAD=FADLADbeck)
KTtest_pro_beck  <-disparity.test(KTdata_pro_beck, method="centroid", test="pairwise", bootstraps=1000)
