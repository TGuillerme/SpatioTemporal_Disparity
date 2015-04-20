#wd
setwd('~/PhD/Projects/SpatioTemporal_Disparity/Cluster')

#Load package
library(disparity)

#Set variables
chain_name<-"Beck2014"
save_path<-"../Data/"
methods<-c("int", "sli_ran", "sli_acc", "sli_del", "sli_pro")
null<-c("random", "sim.char")
disparity<-c("centroid", "product.range", "product.var", "sum.range", "sum.var")

#Looping through every combination
for(null_type in 1:length(null)) {
    #null_type
    for(methods_type in 1:length(methods)) {
        #methods_type
        for(disparity_type in 1:length(disparity)) {
            #disaprity_type

            #Extract all the simulations results
            #Setting the parameters
            chain_path<-paste(chain_name, "/", sep="")
            chain_pattern<-paste(chain_name, "-diparity_null-", methods[methods_type], "-", null[null_type], "_", disparity[disparity_type], sep="")
            list_match<-list.files(path=chain_path, pattern=chain_pattern)

            message(paste("Starting ", chain_pattern, ":", sep=""), appendLF=FALSE)
            #Extracting the values
            diversity_out<-list() ; disparity_out<-list() ; values_out<-list()
            for (file in 1:length(list_match)) {
                file_load<-load(paste(chain_path, "/",list_match[file], sep=""))
                diversity_out[[file]]<-get(file_load)[[1]]
                message('.', appendLF=FALSE)
                disparity_out[[file]]<-get(file_load)[[2]][[1]]
                message('.', appendLF=FALSE)
                values_out[[file]]<-get(file_load)[[2]][[2]]
                message('.', appendLF=FALSE)
            }

            #Combining the results
            disparity_null<-list("diversity"=diversity_out[[1]], "disaprity"=combine.disp(disparity_out), "values"=values_out)

            #Saving the results
            save(disparity_null, file=paste(save_path, chain_name, "/", chain_name, "-disparity_null-", methods[methods_type], "-", null[null_type], "-", disparity[disparity_type], ".Rda", sep=""))
            message("Done.\n", appendLF=FALSE)
        }
    }
}