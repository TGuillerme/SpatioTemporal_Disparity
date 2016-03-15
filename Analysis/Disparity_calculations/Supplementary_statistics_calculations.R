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

#Setting the variables
#Constant
data_path<-"../Data/"
disparity_int<-"-disp_int_tips.Rda"
disparity_ino<-"-disp_int_nodes95.Rda"
disparity_ran<-"-disp_sli_nodes95_ran.Rda"
disparity_acc<-"-disp_sli_nodes95_acc.Rda"
disaprity_del<-"-disp_sli_nodes95_del.Rda"
disaprity_pro<-"-disp_sli_nodes95_pro.Rda"
distance_node<-"_distance-nodes95.Rda"
distance_inte<-"_distance-tips.Rda"

FADLADbeck<-'../Data/Beck2014_FADLAD.csv'
FADLADslat<-'../Data/Slater2013_FADLAD.csv'

#time
slices<-as.numeric(strsplit(c(noquote('170,165,160,155,150,145,140,135,130,125,120,115,110,105,100,95,90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,5,0')), split=',')[[1]])
intervals=as.numeric(strsplit(c(noquote('170.300,168.300,166.100,163.500,157.300,152.100,145.000,139.800,132.900,129.400,125.000,113.000,100.500,93.900,89.800,86.300,83.600,72.100,66.000,61.600,59.200,56.000,47.800,41.300,38.000,33.900,28.100,23.030,23.030,20.440,15.970,13.820,11.620,7.246,5.333,0.000')), split=',')[[1]])

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
distance_nodes_slater<-get(load(paste(strsplit(paste(data_path, chain_name[1], "/", chain_name[1], distance_node, sep=""), ".Rda")[[1]], ".trimmed", sep="")))
distance_nodes_beck  <-get(load(paste(data_path, chain_name[2], "/", chain_name[2], distance_node, sep="")))$gower.dist.matrix
distance_inter_slater<-get(load(paste(strsplit(paste(data_path, chain_name[1], "/", chain_name[1], distance_inte, sep=""), ".Rda")[[1]], ".trimmed", sep="")))
distance_inter_beck  <-get(load(paste(data_path, chain_name[2], "/", chain_name[2], distance_inte, sep="")))$gower.dist.matrix
distance_intno_slater<-get(load(paste(strsplit(paste(data_path, chain_name[1], "/", chain_name[1], distance_node, sep=""), ".Rda")[[1]], ".trimmed", sep="")))
distance_intno_beck  <-get(load(paste(data_path, chain_name[2], "/", chain_name[2], distance_node, sep="")))$gower.dist.matrix

#Extracting the data
#Trees
tree_beck_nod  <-read.data(chain_name[2], data_path, file_matrix[2], file_tree[2], disaprity_pro)[[2]]
tree_beck_tip  <-read.data(chain_name[2], data_path, file_matrix[2], file_tree[2], disparity_int)[[2]]
tree_slater_nod<-get(load(paste(strsplit(file_tree[1], ".tre")[[1]], ".trimmed_nod", sep="")))
tree_slater_tip<-get(load(paste(strsplit(file_tree[1], ".tre")[[1]], ".trimmed_tip", sep="")))

#Running the pco
pco_nodes_slater<-cmdscale(distance_nodes_slater, k=nrow(distance_nodes_slater) - 2, add=T)$points
pco_nodes_beck  <-cmdscale(distance_nodes_beck  , k=nrow(distance_nodes_beck  ) - 2, add=T)$points
pco_inter_slater<-cmdscale(distance_inter_slater, k=nrow(distance_inter_slater) - 2, add=T)$points
pco_inter_beck  <-cmdscale(distance_inter_beck  , k=nrow(distance_inter_beck  ) - 2, add=T)$points
pco_intno_slater<-cmdscale(distance_intno_slater, k=nrow(distance_intno_slater) - 2, add=T)$points
pco_intno_beck  <-cmdscale(distance_intno_beck  , k=nrow(distance_intno_beck  ) - 2, add=T)$points

#Intervals
pco_int_slater_int<-int.pco(pco_inter_slater, tree_slater_tip, intervals,  include.nodes=FALSE, FAD_LAD=FADLADslat)
pco_int_beck_int  <-int.pco(pco_inter_beck  , tree_beck_tip  , intervals,  include.nodes=FALSE, FAD_LAD=FADLADbeck)
pco_int_slater_ino<-int.pco(pco_intno_slater, tree_slater_nod, intervals,  include.nodes=TRUE , FAD_LAD=FADLADslat)
pco_int_beck_ino  <-int.pco(pco_intno_beck  , tree_beck_nod  , intervals,  include.nodes=TRUE , FAD_LAD=FADLADbeck)
#correcting the intervals (a minimum of three species per interval)
pco_int_slater_int<-cor.time.pco(pco_int_slater_int)
pco_int_slater_ino<-cor.time.pco(pco_int_slater_ino)
pco_int_beck_int  <-cor.time.pco(pco_int_beck_int  )
pco_int_beck_ino  <-cor.time.pco(pco_int_beck_ino  )

#Slicing
#rand
pco_slice_slater_ran<-slice.pco(pco_nodes_slater, tree_slater_nod, slices, method='random', FAD_LAD=FADLADslat, verbose=TRUE)
pco_slice_beck_ran  <-slice.pco(pco_nodes_beck  , tree_beck_nod  , slices, method='random', FAD_LAD=FADLADbeck, verbose=TRUE)
#acc
pco_slice_slater_acc<-slice.pco(pco_nodes_slater, tree_slater_nod, slices, method='acctran', FAD_LAD=FADLADslat, verbose=TRUE)
pco_slice_beck_acc  <-slice.pco(pco_nodes_beck  , tree_beck_nod  , slices, method='acctran', FAD_LAD=FADLADbeck, verbose=TRUE)
#del
pco_slice_slater_del<-slice.pco(pco_nodes_slater, tree_slater_nod, slices, method='deltran', FAD_LAD=FADLADslat, verbose=TRUE)
pco_slice_beck_del  <-slice.pco(pco_nodes_beck  , tree_beck_nod  , slices, method='deltran', FAD_LAD=FADLADbeck, verbose=TRUE)
#pro
pco_slice_slater_pro<-slice.pco(pco_nodes_slater, tree_slater_nod, slices, method='proximity', FAD_LAD=FADLADslat, verbose=TRUE)
pco_slice_beck_pro  <-slice.pco(pco_nodes_beck  , tree_beck_nod  , slices, method='proximity', FAD_LAD=FADLADbeck, verbose=TRUE)

######################################
# Significance testing
######################################

#PERMANOVA
#Intervals - Slater
permanova_int_slater  <-disparity.test.time(pco_int_slater_int, method="euclidean", permutations=1000)
permanova_ino_slater  <-disparity.test.time(pco_int_slater_ino, method="euclidean", permutations=1000)
#Slices - Slater
permanova_ran_slater  <-disparity.test.time(pco_slice_slater_ran, method="euclidean", permutations=1000)
permanova_acc_slater  <-disparity.test.time(pco_slice_slater_acc, method="euclidean", permutations=1000)
permanova_del_slater  <-disparity.test.time(pco_slice_slater_del, method="euclidean", permutations=1000)
permanova_pro_slater  <-disparity.test.time(pco_slice_slater_pro, method="euclidean", permutations=1000)
#Saving
permanovas_slater<-list("intervals"=permanova_int_slater, "int+nodes"=permanova_ino_slater, "punctuate"=permanova_ran_slater, "acctran"=permanova_acc_slater, "deltran"=permanova_del_slater, "gradual"=permanova_pro_slater)
save(permanovas_slater, file="permanova_slater.Rda")


#Intervals - Beck
permanova_int_beck  <-disparity.test.time(pco_int_beck_int  , method="euclidean", permutations=1000)
permanova_ino_beck  <-disparity.test.time(pco_int_beck_ino  , method="euclidean", permutations=1000)
#Slices - Beck
permanova_ran_beck  <-disparity.test.time(pco_slice_beck_ran  , method="euclidean", permutations=1000)
permanova_acc_beck  <-disparity.test.time(pco_slice_beck_acc  , method="euclidean", permutations=1000)
permanova_del_beck  <-disparity.test.time(pco_slice_beck_del  , method="euclidean", permutations=1000)
permanova_pro_beck  <-disparity.test.time(pco_slice_beck_pro  , method="euclidean", permutations=1000)
#Saving
permanovas_beck<-list("intervals"=permanova_int_beck, "int+nodes"=permanova_ino_beck, "punctuate"=permanova_ran_beck, "acctran"=permanova_acc_beck, "deltran"=permanova_del_beck, "gradual"=permanova_pro_beck)
save(permanovas_beck, file="permanova_beck.Rda")

##################################
#
#REFERENCE TESTING
#
##################################
#centroid
#################
#slater

reftest_centroid_int_slater<-disparity.test(pco_int_slater_int[17:19], method="centroid", test="reference", bootstraps=1000)
reftest_centroid_ino_slater<-disparity.test(pco_int_slater_ino[17:23], method="centroid", test="reference", bootstraps=1000)
reftest_centroid_ran_slater<-disparity.test(pco_slice_slater_ran[22:35], method="centroid", test="reference", bootstraps=1000)
reftest_centroid_acc_slater<-disparity.test(pco_slice_slater_acc[22:35], method="centroid", test="reference", bootstraps=1000)
reftest_centroid_del_slater<-disparity.test(pco_slice_slater_del[22:35], method="centroid", test="reference", bootstraps=1000)
reftest_centroid_pro_slater<-disparity.test(pco_slice_slater_pro[22:35], method="centroid", test="reference", bootstraps=1000)
#Saving
reftests_centroid_slater<-list("intervals"=reftest_centroid_int_slater, "int+nodes"=reftest_centroid_ino_slater, "punctuate"=reftest_centroid_ran_slater, "acctran"=reftest_centroid_acc_slater, "deltran"=reftest_centroid_del_slater, "gradual"=reftest_centroid_pro_slater)
save(reftests_centroid_slater, file="reftests_centroid_slater.Rda")

#beck
reftest_centroid_int_beck<-disparity.test(pco_int_beck_int[9:17], method="centroid", test="reference", bootstraps=1000)
reftest_centroid_ino_beck<-disparity.test(pco_int_beck_ino[13:21], method="centroid", test="reference", bootstraps=1000)
reftest_centroid_ran_beck<-disparity.test(pco_slice_beck_ran[22:35], method="centroid", test="reference", bootstraps=1000)
reftest_centroid_acc_beck<-disparity.test(pco_slice_beck_acc[22:35], method="centroid", test="reference", bootstraps=1000)
reftest_centroid_del_beck<-disparity.test(pco_slice_beck_del[22:35], method="centroid", test="reference", bootstraps=1000)
reftest_centroid_pro_beck<-disparity.test(pco_slice_beck_pro[22:35], method="centroid", test="reference", bootstraps=1000)
#Saving
reftests_centroid_beck<-list("intervals"=reftest_centroid_int_beck, "int+nodes"=reftest_centroid_ino_beck, "punctuate"=reftest_centroid_ran_beck, "acctran"=reftest_centroid_acc_beck, "deltran"=reftest_centroid_del_beck, "gradual"=reftest_centroid_pro_beck)
save(reftests_centroid_beck, file="reftests_centroid_beck.Rda")

#################
#sum.range
#################
#slater
reftest_sum.range_int_slater<-disparity.test(pco_int_slater_int[17:19], method="sum.range", test="reference", bootstraps=1000)
reftest_sum.range_ino_slater<-disparity.test(pco_int_slater_ino[17:23], method="sum.range", test="reference", bootstraps=1000)
reftest_sum.range_ran_slater<-disparity.test(pco_slice_slater_ran[22:35], method="sum.range", test="reference", bootstraps=1000)
reftest_sum.range_acc_slater<-disparity.test(pco_slice_slater_acc[22:35], method="sum.range", test="reference", bootstraps=1000)
reftest_sum.range_del_slater<-disparity.test(pco_slice_slater_del[22:35], method="sum.range", test="reference", bootstraps=1000)
reftest_sum.range_pro_slater<-disparity.test(pco_slice_slater_pro[22:35], method="sum.range", test="reference", bootstraps=1000)
#Saving
reftests_sum.range_slater<-list("intervals"=reftest_sum.range_int_slater, "int+nodes"=reftest_sum.range_ino_slater, "punctuate"=reftest_sum.range_ran_slater, "acctran"=reftest_sum.range_acc_slater, "deltran"=reftest_sum.range_del_slater, "gradual"=reftest_sum.range_pro_slater)
save(reftests_sum.range_slater, file="reftests_sum.range_slater.Rda")

#beck
reftest_sum.range_int_beck<-disparity.test(pco_int_beck_int[9:17], method="sum.range", test="reference", bootstraps=1000)
reftest_sum.range_ino_beck<-disparity.test(pco_int_beck_ino[13:21], method="sum.range", test="reference", bootstraps=1000)
reftest_sum.range_ran_beck<-disparity.test(pco_slice_beck_ran[22:35], method="sum.range", test="reference", bootstraps=1000)
reftest_sum.range_acc_beck<-disparity.test(pco_slice_beck_acc[22:35], method="sum.range", test="reference", bootstraps=1000)
reftest_sum.range_del_beck<-disparity.test(pco_slice_beck_del[22:35], method="sum.range", test="reference", bootstraps=1000)
reftest_sum.range_pro_beck<-disparity.test(pco_slice_beck_pro[22:35], method="sum.range", test="reference", bootstraps=1000)
#Saving
reftests_sum.range_beck<-list("intervals"=reftest_sum.range_int_beck, "int+nodes"=reftest_sum.range_ino_beck, "punctuate"=reftest_sum.range_ran_beck, "acctran"=reftest_sum.range_acc_beck, "deltran"=reftest_sum.range_del_beck, "gradual"=reftest_sum.range_pro_beck)
save(reftests_sum.range_beck, file="reftests_sum.range_beck.Rda")

#################
#product.range
#################
#slater
reftest_product.range_int_slater<-disparity.test(pco_int_slater_int[17:19], method="product.range", test="reference", bootstraps=1000)
reftest_product.range_ino_slater<-disparity.test(pco_int_slater_ino[17:23], method="product.range", test="reference", bootstraps=1000)
reftest_product.range_ran_slater<-disparity.test(pco_slice_slater_ran[22:35], method="product.range", test="reference", bootstraps=1000)
reftest_product.range_acc_slater<-disparity.test(pco_slice_slater_acc[22:35], method="product.range", test="reference", bootstraps=1000)
reftest_product.range_del_slater<-disparity.test(pco_slice_slater_del[22:35], method="product.range", test="reference", bootstraps=1000)
reftest_product.range_pro_slater<-disparity.test(pco_slice_slater_pro[22:35], method="product.range", test="reference", bootstraps=1000)
#Saving
reftests_product.range_slater<-list("intervals"=reftest_product.range_int_slater, "int+nodes"=reftest_product.range_ino_slater, "punctuate"=reftest_product.range_ran_slater, "acctran"=reftest_product.range_acc_slater, "deltran"=reftest_product.range_del_slater, "gradual"=reftest_product.range_pro_slater)
save(reftests_product.range_slater, file="reftests_product.range_slater.Rda")

#beck
reftest_product.range_int_beck<-disparity.test(pco_int_beck_int[9:17], method="product.range", test="reference", bootstraps=1000)
reftest_product.range_ino_beck<-disparity.test(pco_int_beck_ino[13:21], method="product.range", test="reference", bootstraps=1000)
reftest_product.range_ran_beck<-disparity.test(pco_slice_beck_ran[22:35], method="product.range", test="reference", bootstraps=1000)
reftest_product.range_acc_beck<-disparity.test(pco_slice_beck_acc[22:35], method="product.range", test="reference", bootstraps=1000)
reftest_product.range_del_beck<-disparity.test(pco_slice_beck_del[22:35], method="product.range", test="reference", bootstraps=1000)
reftest_product.range_pro_beck<-disparity.test(pco_slice_beck_pro[22:35], method="product.range", test="reference", bootstraps=1000)
#Saving
reftests_product.range_beck<-list("intervals"=reftest_product.range_int_beck, "int+nodes"=reftest_product.range_ino_beck, "punctuate"=reftest_product.range_ran_beck, "acctran"=reftest_product.range_acc_beck, "deltran"=reftest_product.range_del_beck, "gradual"=reftest_product.range_pro_beck)
save(reftests_product.range_beck, file="reftests_product.range_beck.Rda")

#################
#sum.variance
#################
#slater
reftest_sum.variance_int_slater<-disparity.test(pco_int_slater_int[17:19], method="sum.variance", test="reference", bootstraps=1000)
reftest_sum.variance_ino_slater<-disparity.test(pco_int_slater_ino[17:23], method="sum.variance", test="reference", bootstraps=1000)
reftest_sum.variance_ran_slater<-disparity.test(pco_slice_slater_ran[22:35], method="sum.variance", test="reference", bootstraps=1000)
reftest_sum.variance_acc_slater<-disparity.test(pco_slice_slater_acc[22:35], method="sum.variance", test="reference", bootstraps=1000)
reftest_sum.variance_del_slater<-disparity.test(pco_slice_slater_del[22:35], method="sum.variance", test="reference", bootstraps=1000)
reftest_sum.variance_pro_slater<-disparity.test(pco_slice_slater_pro[22:35], method="sum.variance", test="reference", bootstraps=1000)
#Saving
reftests_sum.variance_slater<-list("intervals"=reftest_sum.variance_int_slater, "int+nodes"=reftest_sum.variance_ino_slater, "punctuate"=reftest_sum.variance_ran_slater, "acctran"=reftest_sum.variance_acc_slater, "deltran"=reftest_sum.variance_del_slater, "gradual"=reftest_sum.variance_pro_slater)
save(reftests_sum.variance_slater, file="reftests_sum.variance_slater.Rda")

#beck
reftest_sum.variance_int_beck<-disparity.test(pco_int_beck_int[9:17], method="sum.variance", test="reference", bootstraps=1000)
reftest_sum.variance_ino_beck<-disparity.test(pco_int_beck_ino[13:21], method="sum.variance", test="reference", bootstraps=1000)
reftest_sum.variance_ran_beck<-disparity.test(pco_slice_beck_ran[22:35], method="sum.variance", test="reference", bootstraps=1000)
reftest_sum.variance_acc_beck<-disparity.test(pco_slice_beck_acc[22:35], method="sum.variance", test="reference", bootstraps=1000)
reftest_sum.variance_del_beck<-disparity.test(pco_slice_beck_del[22:35], method="sum.variance", test="reference", bootstraps=1000)
reftest_sum.variance_pro_beck<-disparity.test(pco_slice_beck_pro[22:35], method="sum.variance", test="reference", bootstraps=1000)
#Saving
reftests_sum.variance_beck<-list("intervals"=reftest_sum.variance_int_beck, "int+nodes"=reftest_sum.variance_ino_beck, "punctuate"=reftest_sum.variance_ran_beck, "acctran"=reftest_sum.variance_acc_beck, "deltran"=reftest_sum.variance_del_beck, "gradual"=reftest_sum.variance_pro_beck)
save(reftests_sum.variance_beck, file="reftests_sum.variance_beck.Rda")

#################
#product.variance
#################
#slater
reftest_product.variance_int_slater<-disparity.test(pco_int_slater_int[17:19], method="product.variance", test="reference", bootstraps=1000)
reftest_product.variance_ino_slater<-disparity.test(pco_int_slater_ino[17:23], method="product.variance", test="reference", bootstraps=1000)
reftest_product.variance_ran_slater<-disparity.test(pco_slice_slater_ran[22:35], method="product.variance", test="reference", bootstraps=1000)
reftest_product.variance_acc_slater<-disparity.test(pco_slice_slater_acc[22:35], method="product.variance", test="reference", bootstraps=1000)
reftest_product.variance_del_slater<-disparity.test(pco_slice_slater_del[22:35], method="product.variance", test="reference", bootstraps=1000)
reftest_product.variance_pro_slater<-disparity.test(pco_slice_slater_pro[22:35], method="product.variance", test="reference", bootstraps=1000)
#Saving
reftests_product.variance_slater<-list("intervals"=reftest_product.variance_int_slater, "int+nodes"=reftest_product.variance_ino_slater, "punctuate"=reftest_product.variance_ran_slater, "acctran"=reftest_product.variance_acc_slater, "deltran"=reftest_product.variance_del_slater, "gradual"=reftest_product.variance_pro_slater)
save(reftests_product.variance_slater, file="reftests_product.variance_slater.Rda")

#beck
reftest_product.variance_int_beck<-disparity.test(pco_int_beck_int[9:17], method="product.variance", test="reference", bootstraps=1000)
reftest_product.variance_ino_beck<-disparity.test(pco_int_beck_ino[13:21], method="product.variance", test="reference", bootstraps=1000)
reftest_product.variance_ran_beck<-disparity.test(pco_slice_beck_ran[22:35], method="product.variance", test="reference", bootstraps=1000)
reftest_product.variance_acc_beck<-disparity.test(pco_slice_beck_acc[22:35], method="product.variance", test="reference", bootstraps=1000)
reftest_product.variance_del_beck<-disparity.test(pco_slice_beck_del[22:35], method="product.variance", test="reference", bootstraps=1000)
reftest_product.variance_pro_beck<-disparity.test(pco_slice_beck_pro[22:35], method="product.variance", test="reference", bootstraps=1000)
#Saving
reftests_product.variance_beck<-list("intervals"=reftest_product.variance_int_beck, "int+nodes"=reftest_product.variance_ino_beck, "punctuate"=reftest_product.variance_ran_beck, "acctran"=reftest_product.variance_acc_beck, "deltran"=reftest_product.variance_del_beck, "gradual"=reftest_product.variance_pro_beck)
save(reftests_product.variance_beck, file="reftests_product.variance_beck.Rda")

##################################
#
#SEQUENTIAL TESTING
#
##################################
#centroid
#################
#slater
seqtest_centroid_int_slater<-disparity.test(pco_int_slater_int, method="centroid", test="sequential", bootstraps=1000)
seqtest_centroid_ino_slater<-disparity.test(pco_int_slater_ino, method="centroid", test="sequential", bootstraps=1000)
seqtest_centroid_ran_slater<-disparity.test(pco_slice_slater_ran, method="centroid", test="sequential", bootstraps=1000)
seqtest_centroid_acc_slater<-disparity.test(pco_slice_slater_acc, method="centroid", test="sequential", bootstraps=1000)
seqtest_centroid_del_slater<-disparity.test(pco_slice_slater_del, method="centroid", test="sequential", bootstraps=1000)
seqtest_centroid_pro_slater<-disparity.test(pco_slice_slater_pro, method="centroid", test="sequential", bootstraps=1000)
#Saving
seqtests_centroid_slater<-list("intervals"=seqtest_centroid_int_slater, "int+nodes"=seqtest_centroid_ino_slater, "punctuate"=seqtest_centroid_ran_slater, "acctran"=seqtest_centroid_acc_slater, "deltran"=seqtest_centroid_del_slater, "gradual"=seqtest_centroid_pro_slater)
save(seqtests_centroid_slater, file="seqtests_centroid_slater.Rda")

#beck
seqtest_centroid_int_beck<-disparity.test(pco_int_beck_int, method="centroid", test="sequential", bootstraps=1000)
seqtest_centroid_ino_beck<-disparity.test(pco_int_beck_ino, method="centroid", test="sequential", bootstraps=1000)
seqtest_centroid_ran_beck<-disparity.test(pco_slice_beck_ran, method="centroid", test="sequential", bootstraps=1000)
seqtest_centroid_acc_beck<-disparity.test(pco_slice_beck_acc, method="centroid", test="sequential", bootstraps=1000)
seqtest_centroid_del_beck<-disparity.test(pco_slice_beck_del, method="centroid", test="sequential", bootstraps=1000)
seqtest_centroid_pro_beck<-disparity.test(pco_slice_beck_pro, method="centroid", test="sequential", bootstraps=1000)
#Saving
seqtests_centroid_beck<-list("intervals"=seqtest_centroid_int_beck, "int+nodes"=seqtest_centroid_ino_beck, "punctuate"=seqtest_centroid_ran_beck, "acctran"=seqtest_centroid_acc_beck, "deltran"=seqtest_centroid_del_beck, "gradual"=seqtest_centroid_pro_beck)
save(seqtests_centroid_beck, file="seqtests_centroid_beck.Rda")

#################
#sum.range
#################
#slater
seqtest_sum.range_int_slater<-disparity.test(pco_int_slater_int, method="sum.range", test="sequential", bootstraps=1000)
seqtest_sum.range_ino_slater<-disparity.test(pco_int_slater_ino, method="sum.range", test="sequential", bootstraps=1000)
seqtest_sum.range_ran_slater<-disparity.test(pco_slice_slater_ran, method="sum.range", test="sequential", bootstraps=1000)
seqtest_sum.range_acc_slater<-disparity.test(pco_slice_slater_acc, method="sum.range", test="sequential", bootstraps=1000)
seqtest_sum.range_del_slater<-disparity.test(pco_slice_slater_del, method="sum.range", test="sequential", bootstraps=1000)
seqtest_sum.range_pro_slater<-disparity.test(pco_slice_slater_pro, method="sum.range", test="sequential", bootstraps=1000)
#Saving
seqtests_sum.range_slater<-list("intervals"=seqtest_sum.range_int_slater, "int+nodes"=seqtest_sum.range_ino_slater, "punctuate"=seqtest_sum.range_ran_slater, "acctran"=seqtest_sum.range_acc_slater, "deltran"=seqtest_sum.range_del_slater, "gradual"=seqtest_sum.range_pro_slater)
save(seqtests_sum.range_slater, file="seqtests_sum.range_slater.Rda")

#beck
seqtest_sum.range_int_beck<-disparity.test(pco_int_beck_int, method="sum.range", test="sequential", bootstraps=1000)
seqtest_sum.range_ino_beck<-disparity.test(pco_int_beck_ino, method="sum.range", test="sequential", bootstraps=1000)
seqtest_sum.range_ran_beck<-disparity.test(pco_slice_beck_ran, method="sum.range", test="sequential", bootstraps=1000)
seqtest_sum.range_acc_beck<-disparity.test(pco_slice_beck_acc, method="sum.range", test="sequential", bootstraps=1000)
seqtest_sum.range_del_beck<-disparity.test(pco_slice_beck_del, method="sum.range", test="sequential", bootstraps=1000)
seqtest_sum.range_pro_beck<-disparity.test(pco_slice_beck_pro, method="sum.range", test="sequential", bootstraps=1000)
#Saving
seqtests_sum.range_beck<-list("intervals"=seqtest_sum.range_int_beck, "int+nodes"=seqtest_sum.range_ino_beck, "punctuate"=seqtest_sum.range_ran_beck, "acctran"=seqtest_sum.range_acc_beck, "deltran"=seqtest_sum.range_del_beck, "gradual"=seqtest_sum.range_pro_beck)
save(seqtests_sum.range_beck, file="seqtests_sum.range_beck.Rda")

#################
#product.range
#################
#slater
seqtest_product.range_int_slater<-disparity.test(pco_int_slater_int, method="product.range", test="sequential", bootstraps=1000)
seqtest_product.range_ino_slater<-disparity.test(pco_int_slater_ino, method="product.range", test="sequential", bootstraps=1000)
seqtest_product.range_ran_slater<-disparity.test(pco_slice_slater_ran, method="product.range", test="sequential", bootstraps=1000)
seqtest_product.range_acc_slater<-disparity.test(pco_slice_slater_acc, method="product.range", test="sequential", bootstraps=1000)
seqtest_product.range_del_slater<-disparity.test(pco_slice_slater_del, method="product.range", test="sequential", bootstraps=1000)
seqtest_product.range_pro_slater<-disparity.test(pco_slice_slater_pro, method="product.range", test="sequential", bootstraps=1000)
#Saving
seqtests_product.range_slater<-list("intervals"=seqtest_product.range_int_slater, "int+nodes"=seqtest_product.range_ino_slater, "punctuate"=seqtest_product.range_ran_slater, "acctran"=seqtest_product.range_acc_slater, "deltran"=seqtest_product.range_del_slater, "gradual"=seqtest_product.range_pro_slater)
save(seqtests_product.range_slater, file="seqtests_product.range_slater.Rda")

#beck
seqtest_product.range_int_beck<-disparity.test(pco_int_beck_int, method="product.range", test="sequential", bootstraps=1000)
seqtest_product.range_ino_beck<-disparity.test(pco_int_beck_ino, method="product.range", test="sequential", bootstraps=1000)
seqtest_product.range_ran_beck<-disparity.test(pco_slice_beck_ran, method="product.range", test="sequential", bootstraps=1000)
seqtest_product.range_acc_beck<-disparity.test(pco_slice_beck_acc, method="product.range", test="sequential", bootstraps=1000)
seqtest_product.range_del_beck<-disparity.test(pco_slice_beck_del, method="product.range", test="sequential", bootstraps=1000)
seqtest_product.range_pro_beck<-disparity.test(pco_slice_beck_pro, method="product.range", test="sequential", bootstraps=1000)
#Saving
seqtests_product.range_beck<-list("intervals"=seqtest_product.range_int_beck, "int+nodes"=seqtest_product.range_ino_beck, "punctuate"=seqtest_product.range_ran_beck, "acctran"=seqtest_product.range_acc_beck, "deltran"=seqtest_product.range_del_beck, "gradual"=seqtest_product.range_pro_beck)
save(seqtests_product.range_beck, file="seqtests_product.range_beck.Rda")

#################
#sum.variance
#################
#slater
seqtest_sum.variance_int_slater<-disparity.test(pco_int_slater_int, method="sum.variance", test="sequential", bootstraps=1000)
seqtest_sum.variance_ino_slater<-disparity.test(pco_int_slater_ino, method="sum.variance", test="sequential", bootstraps=1000)
seqtest_sum.variance_ran_slater<-disparity.test(pco_slice_slater_ran, method="sum.variance", test="sequential", bootstraps=1000)
seqtest_sum.variance_acc_slater<-disparity.test(pco_slice_slater_acc, method="sum.variance", test="sequential", bootstraps=1000)
seqtest_sum.variance_del_slater<-disparity.test(pco_slice_slater_del, method="sum.variance", test="sequential", bootstraps=1000)
seqtest_sum.variance_pro_slater<-disparity.test(pco_slice_slater_pro, method="sum.variance", test="sequential", bootstraps=1000)
#Saving
seqtests_sum.variance_slater<-list("intervals"=seqtest_sum.variance_int_slater, "int+nodes"=seqtest_sum.variance_ino_slater, "punctuate"=seqtest_sum.variance_ran_slater, "acctran"=seqtest_sum.variance_acc_slater, "deltran"=seqtest_sum.variance_del_slater, "gradual"=seqtest_sum.variance_pro_slater)
save(seqtests_sum.variance_slater, file="seqtests_sum.variance_slater.Rda")

#beck
seqtest_sum.variance_int_beck<-disparity.test(pco_int_beck_int, method="sum.variance", test="sequential", bootstraps=1000)
seqtest_sum.variance_ino_beck<-disparity.test(pco_int_beck_ino, method="sum.variance", test="sequential", bootstraps=1000)
seqtest_sum.variance_ran_beck<-disparity.test(pco_slice_beck_ran, method="sum.variance", test="sequential", bootstraps=1000)
seqtest_sum.variance_acc_beck<-disparity.test(pco_slice_beck_acc, method="sum.variance", test="sequential", bootstraps=1000)
seqtest_sum.variance_del_beck<-disparity.test(pco_slice_beck_del, method="sum.variance", test="sequential", bootstraps=1000)
seqtest_sum.variance_pro_beck<-disparity.test(pco_slice_beck_pro, method="sum.variance", test="sequential", bootstraps=1000)
#Saving
seqtests_sum.variance_beck<-list("intervals"=seqtest_sum.variance_int_beck, "int+nodes"=seqtest_sum.variance_ino_beck, "punctuate"=seqtest_sum.variance_ran_beck, "acctran"=seqtest_sum.variance_acc_beck, "deltran"=seqtest_sum.variance_del_beck, "gradual"=seqtest_sum.variance_pro_beck)
save(seqtests_sum.variance_beck, file="seqtests_sum.variance_beck.Rda")

#################
#product.variance
#################
#slater
seqtest_product.variance_int_slater<-disparity.test(pco_int_slater_int, method="product.variance", test="sequential", bootstraps=1000)
seqtest_product.variance_ino_slater<-disparity.test(pco_int_slater_ino, method="product.variance", test="sequential", bootstraps=1000)
seqtest_product.variance_ran_slater<-disparity.test(pco_slice_slater_ran, method="product.variance", test="sequential", bootstraps=1000)
seqtest_product.variance_acc_slater<-disparity.test(pco_slice_slater_acc, method="product.variance", test="sequential", bootstraps=1000)
seqtest_product.variance_del_slater<-disparity.test(pco_slice_slater_del, method="product.variance", test="sequential", bootstraps=1000)
seqtest_product.variance_pro_slater<-disparity.test(pco_slice_slater_pro, method="product.variance", test="sequential", bootstraps=1000)
#Saving
seqtests_product.variance_slater<-list("intervals"=seqtest_product.variance_int_slater, "int+nodes"=seqtest_product.variance_ino_slater, "punctuate"=seqtest_product.variance_ran_slater, "acctran"=seqtest_product.variance_acc_slater, "deltran"=seqtest_product.variance_del_slater, "gradual"=seqtest_product.variance_pro_slater)
save(seqtests_product.variance_slater, file="seqtests_product.variance_slater.Rda")

#beck
seqtest_product.variance_int_beck<-disparity.test(pco_int_beck_int, method="product.variance", test="sequential", bootstraps=1000)
seqtest_product.variance_ino_beck<-disparity.test(pco_int_beck_ino, method="product.variance", test="sequential", bootstraps=1000)
seqtest_product.variance_ran_beck<-disparity.test(pco_slice_beck_ran, method="product.variance", test="sequential", bootstraps=1000)
seqtest_product.variance_acc_beck<-disparity.test(pco_slice_beck_acc, method="product.variance", test="sequential", bootstraps=1000)
seqtest_product.variance_del_beck<-disparity.test(pco_slice_beck_del, method="product.variance", test="sequential", bootstraps=1000)
seqtest_product.variance_pro_beck<-disparity.test(pco_slice_beck_pro, method="product.variance", test="sequential", bootstraps=1000)
#Saving
seqtests_product.variance_beck<-list("intervals"=seqtest_product.variance_int_beck, "int+nodes"=seqtest_product.variance_ino_beck, "punctuate"=seqtest_product.variance_ran_beck, "acctran"=seqtest_product.variance_acc_beck, "deltran"=seqtest_product.variance_del_beck, "gradual"=seqtest_product.variance_pro_beck)
save(seqtests_product.variance_beck, file="seqtests_product.variance_beck.Rda")
