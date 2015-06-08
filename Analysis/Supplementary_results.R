#Loading the package
library(devtools)
install_github("TGuillerme/SpatioTemporal_Disparity/Functions/disparity")
library(disparity)

###################################
#
# SUPPLEMENTARY FIGURES
#
###################################

#SLATER

#Selecting the slater data
chain_name<-"Slater2013"
data_path<-"../Data/"
file_matrix<-"../Data/2013-Slater-MEE-matrix-morpho.nex"
file_tree<-"../Data/2013-Slater-MEE-TEM.tre"
disparity_data<-"-disp_sli_nodes95_pro.Rda"

#SLICES

#Isolating the slater diversity
diversity_data<-"-diversity_slice.Rda"
div<-load(paste(data_path, chain_name, "/", chain_name, diversity_data, sep=""))
diversity_full_slater<-get(div)$diversity

#Extracting all the data
tmp<-read.data(chain_name, data_path, file_matrix, file_tree, disparity_data)
#Isolating the slater data
disparity_full_pro_slater<-tmp[[3]]
#Isolating the slater tree
tree_slater<-tmp[[2]]

#Extracting the random slicing data
disparity_data<-"-disp_sli_nodes95_ran.Rda"
tmp<-read.data(chain_name, data_path, file_matrix, file_tree, disparity_data)
#Isolating the slater data
disparity_full_ran_slater<-tmp[[3]]

#Extracting the random slicing data
disparity_data<-"-disp_sli_nodes95_acc.Rda"
tmp<-read.data(chain_name, data_path, file_matrix, file_tree, disparity_data)
#Isolating the slater data
disparity_full_acc_slater<-tmp[[3]]

#Extracting the random slicing data
disparity_data<-"-disp_sli_nodes95_del.Rda"
tmp<-read.data(chain_name, data_path, file_matrix, file_tree, disparity_data)
#Isolating the slater data
disparity_full_del_slater<-tmp[[3]]

#INTERVALS

#Extracting the intervals diversities
diversity_data<-"-div_int_tip.Rda"
div<-load(paste(data_path, chain_name, "/", chain_name, diversity_data, sep=""))
diversity_full_tip_slater<-get(div)

diversity_data<-"-div_int_nod.Rda"
div<-load(paste(data_path, chain_name, "/", chain_name, diversity_data, sep=""))
diversity_full_nod_slater<-get(div)

#Extracting the tip intervals
disparity_data<-"-disp_int_tips.Rda"
tmp<-read.data(chain_name, data_path, file_matrix, file_tree, disparity_data)
#Isolating the beck data
disparity_full_tip_slater<-tmp[[3]]

#Extracting the tip+nodes intervals
disparity_data<-"-disp_int_nodes.Rda"
tmp<-read.data(chain_name, data_path, file_matrix, file_tree, disparity_data)
#Isolating the beck data
disparity_full_node_slater<-tmp[[3]]

#BECK

#Selecting the Beck data
chain_name<-"Beck2014"
data_path<-"../Data/"
file_matrix<-"../Data/2014-Beck-ProcB-matrix-morpho.nex"
file_tree<-"../Data/2014-Beck-ProcB-TEM.tre"
disparity_data<-"-disp_sli_nodes95_pro.Rda"

#SLICES

#Isolating the beck diversity
diversity_data<-"-diversity_slice.Rda"
div<-load(paste(data_path, chain_name, "/", chain_name, diversity_data, sep=""))
diversity_full_beck<-get(div)$diversity

#Extracting all the data
tmp<-read.data(chain_name, data_path, file_matrix, file_tree, disparity_data)
#Isolating the beck data
disparity_full_pro_beck<-tmp[[3]]
#Isolating the beck tree
tree_beck<-tmp[[2]]

#Extracting the random slicing data
disparity_data<-"-disp_sli_nodes95_ran.Rda"
tmp<-read.data(chain_name, data_path, file_matrix, file_tree, disparity_data)
#Isolating the beck data
disparity_full_ran_beck<-tmp[[3]]

#Extracting the random slicing data
disparity_data<-"-disp_sli_nodes95_acc.Rda"
tmp<-read.data(chain_name, data_path, file_matrix, file_tree, disparity_data)
#Isolating the beck data
disparity_full_acc_beck<-tmp[[3]]

#Extracting the random slicing data
disparity_data<-"-disp_sli_nodes95_del.Rda"
tmp<-read.data(chain_name, data_path, file_matrix, file_tree, disparity_data)
#Isolating the beck data
disparity_full_del_beck<-tmp[[3]]

#INTERVALS

#Extracting the intervals diversities
diversity_data<-"-div_int_tip.Rda"
div<-load(paste(data_path, chain_name, "/", chain_name, diversity_data, sep=""))
diversity_full_tip_beck<-get(div)

diversity_data<-"-div_int_nod.Rda"
div<-load(paste(data_path, chain_name, "/", chain_name, diversity_data, sep=""))
diversity_full_nod_beck<-get(div)

#Extracting the tip intervals
disparity_data<-"-disp_int_tips.Rda"
tmp<-read.data(chain_name, data_path, file_matrix, file_tree, disparity_data)
#Isolating the beck data
disparity_full_tip_beck<-tmp[[3]]

#Extracting the tip+nodes intervals
disparity_data<-"-disp_int_nodes.Rda"
tmp<-read.data(chain_name, data_path, file_matrix, file_tree, disparity_data)
#Isolating the beck data
disparity_full_node_beck<-tmp[[3]]

######################################
#Rarefaction analysis on both data sets
######################################

#Slater
dis_ran_max_slater<-extract.disp(disparity_full_pro_slater$quantiles, rarefaction="max")
dis_ran_min_slater<-extract.disp(disparity_full_pro_slater$quantiles, rarefaction="min")
dis_ran_mod_slater<-extract.disp(disparity_full_pro_slater$quantiles, rarefaction=mode.val(diversity_full_slater))

#Plot
op<-par(mfrow=c(3,1), bty="n", mar=c(4,4,4,4))
plot.disparity(dis_ran_max_slater, diversity=log(dis_ran_max_slater$rarefaction), main="Maximum (all taxa)", xlab="Time (Mya)", y2lab="Diversity (log)")
abline(v=22, col="red")
plot.disparity(dis_ran_mod_slater, diversity=log(dis_ran_mod_slater$rarefaction), main=paste("Mode (", mode.val(diversity_full_slater), " taxa)", sep=""), xlab="Time (Mya)", y2lab="Diversity (log)")
abline(v=22, col="red")
plot.disparity(dis_ran_min_slater, diversity=log(dis_ran_min_slater$rarefaction), main=paste("Minimum (", unique(dis_ran_min_slater$rarefaction), " taxa)", sep=""), xlab="Time (Mya)", y2lab="Diversity (log)")
abline(v=22, col="red")
par(op)

#Beck
dis_ran_max_beck<-extract.disp(disparity_full_pro_beck$quantiles, rarefaction="max")
dis_ran_min_beck<-extract.disp(disparity_full_pro_beck$quantiles, rarefaction="min")
dis_ran_mod_beck<-extract.disp(disparity_full_pro_beck$quantiles, rarefaction=mode.val(diversity_full_beck))

#Plot
op<-par(mfrow=c(3,1), bty="n", mar=c(4,4,4,4))
plot.disparity(dis_ran_max_beck, diversity=log(dis_ran_max_beck$rarefaction), main="Maximum (all taxa)", xlab="Time (Mya)", y2lab="Diversity (log)")
abline(v=22, col="red")
plot.disparity(dis_ran_mod_beck, diversity=log(dis_ran_mod_beck$rarefaction), main=paste("Mode (", mode.val(diversity_full_beck), " taxa)", sep=""), xlab="Time (Mya)", y2lab="Diversity (log)")
abline(v=22, col="red")
plot.disparity(dis_ran_min_beck, diversity=log(dis_ran_min_beck$rarefaction), main=paste("Minimum (", unique(dis_ran_min_slater$rarefaction), " taxa)", sep=""), xlab="Time (Mya)", y2lab="Diversity (log)")
abline(v=22, col="red")
par(op)

######################################
#Evolutionary models variation + metrics variation
######################################

#-------------------
#SLATER
#-------------------

dis_tips_slater <-extract.disp(disparity_full_tip_slater$quantiles, rarefaction="max")
dis_nodes_slater<-extract.disp(disparity_full_node_slater$quantiles, rarefaction="max")
dis_ran_slater  <-extract.disp(disparity_full_ran_slater$quantiles, rarefaction="max")
dis_acc_slater  <-extract.disp(disparity_full_acc_slater$quantiles, rarefaction="max")
dis_del_slater  <-extract.disp(disparity_full_del_slater$quantiles, rarefaction="max")
dis_pro_slater  <-extract.disp(disparity_full_pro_slater$quantiles, rarefaction="max")


#The following is a "BIG" plot comparing all methods/metrics for Slater
quartz(width = 15.6, height = 11.2) #A5 landscape
#Windows dimensions
op<-par(mfrow=c(5, 6), bty="n", mar=c(4,4,4,4))# oma=c(bottom, left, top, right)
#Centroid
plot.disparity(dis_tips_slater, xlab="", ylab="Distance from centroid", measure="Cent.dist", main="Intervals (tips only)", diversity=log(dis_tips_slater$rarefaction), y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_nodes_slater, xlab="", ylab="", measure="Cent.dist", main="Intervals (tips and nodes)", diversity=log(dis_nodes_slater$rarefaction), y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_ran_slater, xlab="", ylab="", measure="Cent.dist", main="Slices (punctuated)", diversity=log(diversity_full_slater), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_slater, xlab="", ylab="", measure="Cent.dist", main="Slices (punctuated:acctran)", diversity=log(diversity_full_slater), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_slater, xlab="", ylab="", measure="Cent.dist", main="Slices (punctuated:deltran)", diversity=log(diversity_full_slater), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_slater, xlab="", ylab="", measure="Cent.dist", main="Slices (gradual)", diversity=log(diversity_full_slater), y2lab="Diversity (log)")
abline(v= 22, col="red")

#Sum of ranges
plot.disparity(dis_tips_slater, xlab="", ylab="Sum of ranges", measure="Sum.range", main="Intervals (tips only)", diversity=log(dis_tips_slater$rarefaction), y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_nodes_slater, xlab="", ylab="", measure="Sum.range", main="Intervals (tips and nodes)", diversity=log(dis_nodes_slater$rarefaction), y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_ran_slater, xlab="", ylab="", measure="Sum.range", main="Slices (punctuated)", diversity=log(diversity_full_slater), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_slater, xlab="", ylab="", measure="Sum.range", main="Slices (punctuated:acctran)", diversity=log(diversity_full_slater), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_slater, xlab="", ylab="", measure="Sum.range", main="Slices (punctuated:deltran)", diversity=log(diversity_full_slater), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_slater, xlab="", ylab="", measure="Sum.range", main="Slices (gradual)", diversity=log(diversity_full_slater), y2lab="Diversity (log)")
abline(v= 22, col="red")

#Sum of variance
plot.disparity(dis_tips_slater, xlab="", ylab="Sum of variances", measure="Sum.var", main="Intervals (tips only)", diversity=log(dis_tips_slater$rarefaction), y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_nodes_slater, xlab="", ylab="", measure="Sum.var", main="Intervals (tips and nodes)", diversity=log(dis_nodes_slater$rarefaction), y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_ran_slater, xlab="", ylab="", measure="Sum.var", main="Slices (punctuated)", diversity=log(diversity_full_slater), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_slater, xlab="", ylab="", measure="Sum.var", main="Slices (punctuated:acctran)", diversity=log(diversity_full_slater), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_slater, xlab="", ylab="", measure="Sum.var", main="Slices (punctuated:deltran)", diversity=log(diversity_full_slater), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_slater, xlab="", ylab="", measure="Sum.var", main="Slices (gradual)", diversity=log(diversity_full_slater), y2lab="Diversity (log)")
abline(v= 22, col="red")

#Product of ranges
plot.disparity(dis_tips_slater, xlab="", ylab="Product of ranges", measure="Prod.range", main="Intervals (tips only)", diversity=log(dis_tips_slater$rarefaction), y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_nodes_slater, xlab="", ylab="", measure="Prod.range", main="Intervals (tips and nodes)", diversity=log(dis_nodes_slater$rarefaction), y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_ran_slater, xlab="", ylab="", measure="Prod.range", main="Slices (punctuated)", diversity=log(diversity_full_slater), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_slater, xlab="", ylab="", measure="Prod.range", main="Slices (punctuated:acctran)", diversity=log(diversity_full_slater), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_slater, xlab="", ylab="", measure="Prod.range", main="Slices (punctuated:deltran)", diversity=log(diversity_full_slater), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_slater, xlab="", ylab="", measure="Prod.range", main="Slices (gradual)", diversity=log(diversity_full_slater), y2lab="Diversity (log)")
abline(v= 22, col="red")

#Product of variance
plot.disparity(dis_tips_slater, xlab="Time bins (Mya)", ylab="Product of variances", measure="Prod.var", main="Intervals (tips only)", diversity=log(dis_tips_slater$rarefaction), y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_nodes_slater, xlab="Time bins (Mya)", ylab="", measure="Prod.var", main="Intervals (tips and nodes)", diversity=log(dis_nodes_slater$rarefaction), y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_ran_slater, xlab="Time (Mya)", ylab="", measure="Prod.var", main="Slices (punctuated)", diversity=log(diversity_full_slater), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_slater, xlab="Time (Mya)", ylab="", measure="Prod.var", main="Slices (punctuated:acctran)", diversity=log(diversity_full_slater), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_slater, xlab="Time (Mya)", ylab="", measure="Prod.var", main="Slices (punctuated:deltran)", diversity=log(diversity_full_slater), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_slater, xlab="Time (Mya)", ylab="", measure="Prod.var", main="Slices (gradual)", diversity=log(diversity_full_slater), y2lab="Diversity (log)")
abline(v= 22, col="red")

par(op)

#-------------------
#BECK
#-------------------


dis_tips_beck <-extract.disp(disparity_full_tip_beck$quantiles, rarefaction="max")
dis_nodes_beck<-extract.disp(disparity_full_node_beck$quantiles, rarefaction="max")
dis_ran_beck  <-extract.disp(disparity_full_ran_beck$quantiles, rarefaction="max")
dis_acc_beck  <-extract.disp(disparity_full_acc_beck$quantiles, rarefaction="max")
dis_del_beck  <-extract.disp(disparity_full_del_beck$quantiles, rarefaction="max")
dis_pro_beck  <-extract.disp(disparity_full_pro_beck$quantiles, rarefaction="max")


#The following is a "BIG" plot comparing all methods/metrics for Slater
quartz(width = 15.6, height = 11.2) #A5 landscape
#Windows dimensions
op<-par(mfrow=c(5, 6), bty="n", mar=c(4,4,4,4))# oma=c(bottom, left, top, right)
#Centroid
plot.disparity(dis_tips_beck, xlab="", ylab="Distance from centroid", measure="Cent.dist", main="Intervals (tips only)", diversity=log(dis_tips_beck$rarefaction), y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_nodes_beck, xlab="", ylab="", measure="Cent.dist", main="Intervals (tips and nodes)", diversity=log(dis_nodes_beck$rarefaction), y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_ran_beck, xlab="", ylab="", measure="Cent.dist", main="Slices (punctuated)", diversity=log(diversity_full_beck), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_beck, xlab="", ylab="", measure="Cent.dist", main="Slices (punctuated:acctran)", diversity=log(diversity_full_beck), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_beck, xlab="", ylab="", measure="Cent.dist", main="Slices (punctuated:deltran)", diversity=log(diversity_full_beck), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_beck, xlab="", ylab="", measure="Cent.dist", main="Slices (gradual)", diversity=log(diversity_full_beck), y2lab="Diversity (log)")
abline(v= 22, col="red")

#Sum of ranges
plot.disparity(dis_tips_beck, xlab="", ylab="Sum of ranges", measure="Sum.range", main="Intervals (tips only)", diversity=log(dis_tips_beck$rarefaction), y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_nodes_beck, xlab="", ylab="", measure="Sum.range", main="Intervals (tips and nodes)", diversity=log(dis_nodes_beck$rarefaction), y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_ran_beck, xlab="", ylab="", measure="Sum.range", main="Slices (punctuated)", diversity=log(diversity_full_beck), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_beck, xlab="", ylab="", measure="Sum.range", main="Slices (punctuated:acctran)", diversity=log(diversity_full_beck), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_beck, xlab="", ylab="", measure="Sum.range", main="Slices (punctuated:deltran)", diversity=log(diversity_full_beck), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_beck, xlab="", ylab="", measure="Sum.range", main="Slices (gradual)", diversity=log(diversity_full_beck), y2lab="Diversity (log)")
abline(v= 22, col="red")

#Sum of variance
plot.disparity(dis_tips_beck, xlab="", ylab="Sum of variances", measure="Sum.var", main="Intervals (tips only)", diversity=log(dis_tips_beck$rarefaction), y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_nodes_beck, xlab="", ylab="", measure="Sum.var", main="Intervals (tips and nodes)", diversity=log(dis_nodes_beck$rarefaction), y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_ran_beck, xlab="", ylab="", measure="Sum.var", main="Slices (punctuated)", diversity=log(diversity_full_beck), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_beck, xlab="", ylab="", measure="Sum.var", main="Slices (punctuated:acctran)", diversity=log(diversity_full_beck), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_beck, xlab="", ylab="", measure="Sum.var", main="Slices (punctuated:deltran)", diversity=log(diversity_full_beck), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_beck, xlab="", ylab="", measure="Sum.var", main="Slices (gradual)", diversity=log(diversity_full_beck), y2lab="Diversity (log)")
abline(v= 22, col="red")

#Product of ranges
plot.disparity(dis_tips_beck, xlab="", ylab="Product of ranges", measure="Prod.range", main="Intervals (tips only)", diversity=log(dis_tips_beck$rarefaction), y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_nodes_beck, xlab="", ylab="", measure="Prod.range", main="Intervals (tips and nodes)", diversity=log(dis_nodes_beck$rarefaction), y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_ran_beck, xlab="", ylab="", measure="Prod.range", main="Slices (punctuated)", diversity=log(diversity_full_beck), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_beck, xlab="", ylab="", measure="Prod.range", main="Slices (punctuated:acctran)", diversity=log(diversity_full_beck), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_beck, xlab="", ylab="", measure="Prod.range", main="Slices (punctuated:deltran)", diversity=log(diversity_full_beck), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_beck, xlab="", ylab="", measure="Prod.range", main="Slices (gradual)", diversity=log(diversity_full_beck), y2lab="Diversity (log)")
abline(v= 22, col="red")

#Product of variance
plot.disparity(dis_tips_beck, xlab="Time bins (Mya)", ylab="Product of variances", measure="Prod.var", main="Intervals (tips only)", diversity=log(dis_tips_beck$rarefaction), y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_nodes_beck, xlab="Time bins (Mya)", ylab="", measure="Prod.var", main="Intervals (tips and nodes)", diversity=log(dis_nodes_beck$rarefaction), y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_ran_beck, xlab="Time (Mya)", ylab="", measure="Prod.var", main="Slices (punctuated)", diversity=log(diversity_full_beck), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_beck, xlab="Time (Mya)", ylab="", measure="Prod.var", main="Slices (punctuated:acctran)", diversity=log(diversity_full_beck), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_beck, xlab="Time (Mya)", ylab="", measure="Prod.var", main="Slices (punctuated:deltran)", diversity=log(diversity_full_beck), y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_beck, xlab="Time (Mya)", ylab="", measure="Prod.var", main="Slices (gradual)", diversity=log(diversity_full_beck), y2lab="Diversity (log)")
abline(v= 22, col="red")

par(op)


######################################
#Wills 1994 measures
######################################

#Slater
op<-par(mfrow=c(2,2), bty="n", mar=c(4,4,4,4))
plot.disparity(dis_ran_max_slater, diversity=log(dis_ran_max_slater$rarefaction), measure="Sum.range", main="Sum of ranges", xlab="Time (Mya)", y2lab="")
abline(v=22, col="red")
plot.disparity(dis_ran_max_slater, diversity=log(dis_ran_max_slater$rarefaction), measure="Sum.var", main="Sum of variance", xlab="Time (Mya)", y2lab="Diversity (log)")
abline(v=22, col="red")
plot.disparity(dis_ran_max_slater, diversity=log(dis_ran_max_slater$rarefaction), measure="Prod.range", main="Product of ranges", xlab="Time (Mya)", y2lab="")
abline(v=22, col="red")
plot.disparity(dis_ran_max_slater, diversity=log(dis_ran_max_slater$rarefaction), measure="Prod.var", main="Product of variance", xlab="Time (Mya)", y2lab="Diversity (log)")
abline(v=22, col="red")
par(op)

#beck
op<-par(mfrow=c(2,2), bty="n", mar=c(4,4,4,4))
plot.disparity(dis_ran_max_beck, diversity=log(dis_ran_max_beck$rarefaction), measure="Sum.range", main="Sum of ranges", xlab="Time (Mya)", y2lab="")
abline(v=22, col="red")
plot.disparity(dis_ran_max_beck, diversity=log(dis_ran_max_beck$rarefaction), measure="Sum.var", main="Sum of variance", xlab="Time (Mya)", y2lab="Diversity (log)")
abline(v=22, col="red")
plot.disparity(dis_ran_max_beck, diversity=log(dis_ran_max_beck$rarefaction), measure="Prod.range", main="Product of ranges", xlab="Time (Mya)", y2lab="")
abline(v=22, col="red")
plot.disparity(dis_ran_max_beck, diversity=log(dis_ran_max_beck$rarefaction), measure="Prod.var", main="Product of variance", xlab="Time (Mya)", y2lab="Diversity (log)")
abline(v=22, col="red")
par(op)


#####################################
#
