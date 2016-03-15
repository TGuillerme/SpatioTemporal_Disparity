#Phylo corret disparity... Is that a thing?

#Loading the package
library(devtools)
install_github("TGuillerme/SpatioTemporal_Disparity/Functions/disparity")
library(disparity)

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
#Diversity
slat_div<-load(paste(data_path, chain_name[1], "/", chain_name[1], diversity_ful, sep=""))

#Extracting the data
#Tree
tree_slater<-slat_tmp1[[2]]
#Disparity
disparity_full_pro_slater<-slat_tmp1[[3]]
disparity_full_ran_slater<-slat_tmp2[[3]]

#Diversity
diversity_full_slater<-get(slat_div)


######################################
# Significance testing (classic)
######################################

#PERMANOVA
permanova_pro_slater<-disparity.test.time(pco_slice_slater_pro, method="euclidean", permutations=1000)
permanova_ran_slater<-disparity.test.time(pco_slice_slater_ran, method="euclidean", permutations=1000)

#reference
reftest_pro_slater<-disparity.test(pco_slice_slater_pro[22:35], method="centroid", test="reference", bootstraps=1000)
reftest_ran_slater<-disparity.test(pco_slice_slater_ran[22:35], method="centroid", test="reference", bootstraps=1000)

#reference (rarefied)
reftestRAR_pro_slater<-disparity.test(pco_slice_slater_pro[22:35], method="centroid", test="reference", bootstraps=1000, rarefaction=8)
reftestRAR_ran_slater<-disparity.test(pco_slice_slater_ran[22:35], method="centroid", test="reference", bootstraps=1000, rarefaction=8)

######################################
# Including phylogeny (non-independent)
######################################

# The idea here is that in the PERMANOVA (pco ~ time), both pco and time are not independent between each other and within each other.
# This is probably not super cool stats-wise...
# Correcting for time might be tricky since it kind of always correlates with everything by essence, is uniderectional, etc. Maybe doing equal (or random?) time slicing might be sufficient.
        # By the way, maybe random time sampling might be a overall better idea! - no to get uniform distribution needs to resample > 10 times subsamples. That's 350 slices. Might as well just do 170 slices directly.
# Correcting for pco might be "easier":
# Each pco sub-sample (i.e. slice) are dependant but we have a phylogeny for all of them. Maybe the idea could be to get a pco average in each sub-sample, corrected for phylogeny, and test for average_pco ~ time instead?