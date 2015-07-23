#Quick analysis to look at the effect of rarefaction

#Loading the package
library(devtools)
install_github("TGuillerme/SpatioTemporal_Disparity/Functions/disparity")
library(disparity)

###################################
#
# LOAD THE DATA
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

###################################
#
# CALCULATE THE DISTRIBUTIONS OVERLAP
#
###################################

######################################
# Function
######################################

#select the slices with highest taxonomic diversity
bhat.coeff.table<-function(disparity, diversity) {
    results<-as.data.frame(matrix(NA, nrow=length(diversity), ncol=7))
    colnames(results)<-c("max. taxa", "no overlap", "BC","%", "overlap", "BC","%")
    rownames(results)<-names(diversity)
    results[,1]<-as.vector(diversity)

    for (slice in 1:nrow(results)) {
        if(length(disparity$values[[slice]]) == 1) {
            results[slice, 2:7]<-NA
        } else {
            #select the max slice
            tot_slice<-disparity$values[[slice]][[length(disparity$values[[slice]])]]
            #select the slice that is significantly different
            min_slice<-1
            BC_min<-bhatt.coeff(tot_slice, disparity$values[[slice]][[min_slice]])
            if(BC_buffer > 0.05) {
                min_slice<-min_slice
                BC_min<-BC_min
            } else {
                while(BC_min < 0.05) {
                    min_slice<-min_slice+1
                    BC_min<-bhatt.coeff(tot_slice, disparity$values[[slice]][[min_slice]])
                }
            }
            #select the slice that is significantly similar
            max_slice<-1
            BC_max<-bhatt.coeff(tot_slice, disparity$values[[slice]][[max_slice]])
            while(BC_max < 0.95) {
                max_slice<-max_slice+1
                BC_max<-bhatt.coeff(tot_slice, disparity$values[[slice]][[max_slice]])
            }  

            #Saving the results
            results[slice, 2]<-min_slice+2
            results[slice, 3]<-BC_min
            results[slice, 4]<-(min_slice+2)/(length(disparity$values[[slice]])+2)
            results[slice, 5]<-max_slice+2
            results[slice, 6]<-BC_max
            results[slice, 7]<-(max_slice+2)/(length(disparity$values[[slice]])+2)
        }
    }
    return(results)
}

warning("Problem with calculating the minimum slice. Seems to be stuck at 3.")

######################################
# Beck
######################################

bhat.coeff.table(disparity_full_ran_beck, diversity_full_beck)
bhat.coeff.table(disparity_full_pro_beck, diversity_full_beck)

######################################
# Slater
######################################

bhat.coeff.table(disparity_full_ran_slater, diversity_full_slater)
bhat.coeff.table(disparity_full_pro_slater, diversity_full_slater)
