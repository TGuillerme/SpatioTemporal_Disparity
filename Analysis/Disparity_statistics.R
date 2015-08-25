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
#permanova_pro_slater<-disparity.test.time(pco_slice_slater_pro, method="euclidean", permutations=1000)
#permanova_pro_beck  <-disparity.test.time(pco_slice_beck_pro  , method="euclidean", permutations=1000)
#permanova_ran_slater<-disparity.test.time(pco_slice_slater_ran, method="euclidean", permutations=1000)
#permanova_ran_beck  <-disparity.test.time(pco_slice_beck_ran  , method="euclidean", permutations=1000)

#reference
reftest_pro_slater<-disparity.test(pco_slice_slater_pro[21:28], method="centroid", test="reference", bootstraps=1000, correction="none")
reftest_pro_beck  <-disparity.test(pco_slice_beck_pro[21:28]  , method="centroid", test="reference", bootstraps=1000, correction="none")
reftest_ran_slater<-disparity.test(pco_slice_slater_ran[21:28], method="centroid", test="reference", bootstraps=1000, correction="none")
reftest_ran_beck  <-disparity.test(pco_slice_beck_ran[21:28]  , method="centroid", test="reference", bootstraps=1000, correction="none")

#reference (rarefied)
reftestRAR_pro_slater<-disparity.test(pco_slice_slater_pro[21:28], method="centroid", test="reference", bootstraps=1000, rarefaction=8, correction="none")
reftestRAR_pro_beck  <-disparity.test(pco_slice_beck_pro[21:28]  , method="centroid", test="reference", bootstraps=1000, rarefaction=8, correction="none")
reftestRAR_ran_slater<-disparity.test(pco_slice_slater_ran[21:28], method="centroid", test="reference", bootstraps=1000, rarefaction=8, correction="none")
reftestRAR_ran_beck  <-disparity.test(pco_slice_beck_ran[21:28]  , method="centroid", test="reference", bootstraps=1000, rarefaction=8, correction="none")


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

# make.table.permanova<-function() {
#     #Fixed rows
#     #Empty matrix
#     permanova_terms<-as.data.frame(matrix(" ", nrow=8, ncol=3))
#     #Column names
#     colnames(permanova_terms)<-c("Data", "model", "terms")#, c(colnames(permanova_pro_beck[[1]][1,])))
#     #Data
#     permanova_terms$Data<-c("Eutherian", rep(" ", 3), "Mammaliformes", rep(" ", 3))
#     permanova_terms$model<-rep(c(c("gradual", rep(" ", 1)),c("punctuate", rep(" ", 1))),2)
#     permanova_terms$terms<-rep(c("time", "Residuals"), 4)

#     #Fixed rows
#     #Empty matrix
#     permanova_results<-as.data.frame(matrix(NA, nrow=8, ncol=6))
#     #Column names
#     colnames(permanova_results)<-c(colnames(permanova_pro_beck[[1]][1,]))

#     #Variable rows
#     permanova_results[1, ]<-as.vector(permanova_pro_beck[[1]][1,])
#     permanova_results[2, ]<-as.vector(permanova_pro_beck[[1]][2,])
#     permanova_results[3, ]<-as.vector(permanova_ran_beck[[1]][1,])
#     permanova_results[4, ]<-as.vector(permanova_ran_beck[[1]][2,])
#     permanova_results[5, ]<-as.vector(permanova_pro_slater[[1]][1,])
#     permanova_results[6, ]<-as.vector(permanova_pro_slater[[1]][2,])
#     permanova_results[7, ]<-as.vector(permanova_ran_slater[[1]][1,])
#     permanova_results[8, ]<-as.vector(permanova_ran_slater[[1]][2,])

#     #Rounding
#     for (n in 1:6) {
#         permanova_results[,n]<-round(permanova_results[,n], digit=n)
#     }

#     return(cbind(permanova_terms, permanova_results))
# }

# #permanova table
# xtable(make.table.permanova(), digits=4)

#lag test table (to modify manually in LaTeX)
xtable(cbind(reftest_pro_beck[[1]][,-5], reftest_ran_beck[[1]][,-5]), digits=3)
#xtable(cbind(reftest_ran_beck[[1]][,-5], reftest_ran_slater[[1]][,-5]))
cat("add: & & gradual & & & & punctuated & & \\\\' in the header.")
cat("replace 'rrrrrrrrr' by 'rrrrr|rrrr'.")

#lag test table (to modify manually in LaTeX)
xtable(cbind(reftest_pro_slater[[1]][,-5], reftest_ran_slater[[1]][,-5]), digits=3)
#xtable(cbind(reftest_ran_slater[[1]][,-5], reftest_ran_slater[[1]][,-5]))
cat("add: & & gradual & & & & punctuated & & \\\\' in the header.")
cat("replace 'rrrrrrrrr' by 'rrrrr|rrrr'.")

xtable(cbind(reftestRAR_pro_beck[[1]][,-5], reftestRAR_ran_beck[[1]][,-5]), digits=3)
#xtable(cbind(reftest_ran_beck[[1]][,-5], reftest_ran_slater[[1]][,-5]))
cat("add: & & gradual & & & & punctuated & & \\\\' in the header.")
cat("replace 'rrrrrrrrr' by 'rrrrr|rrrr'.")

xtable(cbind(reftestRAR_pro_slater[[1]][,-5], reftestRAR_ran_slater[[1]][,-5]), digits=3)
#xtable(cbind(reftest_ran_slater[[1]][,-5], reftest_ran_slater[[1]][,-5]))
cat("add: & & gradual & & & & punctuated & & \\\\' in the header.")
cat("replace 'rrrrrrrrr' by 'rrrrr|rrrr'.")
