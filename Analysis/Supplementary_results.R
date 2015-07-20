#Loading the package
library(devtools)
install_github("TGuillerme/SpatioTemporal_Disparity/Functions/disparity")
library(disparity)

###################################
#
# SUPPLEMENTARY FIGURES
#
###################################


######################################
# Isolating the data
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
diversity_ful<-"-slices_nodes95_div.Rda"
distance_gowr<-"_distance-nodes95.Rda"
intervals=as.numeric(strsplit(c(noquote('170.300,168.300,166.100,163.500,157.300,152.100,145.000,139.800,132.900,129.400,125.000,113.000,100.500,93.900,89.800,86.300,83.600,72.100,66.000,61.600,59.200,56.000,47.800,41.300,38.000,33.900,28.100,23.030,23.030,20.440,15.970,13.820,11.620,7.246,5.333,0.000')), split=',')[[1]])
FADLADbeck<-'../Data/Beck2014_FADLAD.csv'
FADLADslat<-'../Data/Slater2013_FADLAD.csv'

#Data
chain_name<-c("Slater2013", "Beck2014")
file_matrix<-c("../Data/2013-Slater-MEE-matrix-morpho.nex", "../Data/2014-Beck-ProcB-matrix-morpho.nex")
file_tree<-c("../Data/2013-Slater-MEE-TEM.tre","../Data/2014-Beck-ProcB-TEM.tre")

#Loading the data
#Disparity + tree
#slater
slat_tmp1<-read.data(chain_name[1], data_path, file_matrix[1], file_tree[1], disparity_int)
slat_tmp2<-read.data(chain_name[1], data_path, file_matrix[1], file_tree[1], disparity_ino)
slat_tmp3<-read.data(chain_name[1], data_path, file_matrix[1], file_tree[1], disparity_ran)
slat_tmp4<-read.data(chain_name[1], data_path, file_matrix[1], file_tree[1], disparity_acc)
slat_tmp5<-read.data(chain_name[1], data_path, file_matrix[1], file_tree[1], disaprity_del)
slat_tmp6<-read.data(chain_name[1], data_path, file_matrix[1], file_tree[1], disaprity_pro)
#beck
beck_tmp1<-read.data(chain_name[2], data_path, file_matrix[2], file_tree[2], disparity_int)
beck_tmp2<-read.data(chain_name[2], data_path, file_matrix[2], file_tree[2], disparity_ino)
beck_tmp3<-read.data(chain_name[2], data_path, file_matrix[2], file_tree[2], disparity_ran)
beck_tmp4<-read.data(chain_name[2], data_path, file_matrix[2], file_tree[2], disparity_acc)
beck_tmp5<-read.data(chain_name[2], data_path, file_matrix[2], file_tree[2], disaprity_del)
beck_tmp6<-read.data(chain_name[2], data_path, file_matrix[2], file_tree[2], disaprity_pro)


#Extracting the data
#Tree
tree_slater<-slat_tmp1[[2]]
tree_beck  <-beck_tmp1[[2]]

#FAD/LAD
FADLADslat<-read.csv(FADLADslat, row.names=1)
FADLADbeck<-read.csv(FADLADbeck, row.names=1)

#Disparity
#Slater
disparity_full_int_slater<-slat_tmp1[[3]]
disparity_full_ino_slater<-slat_tmp2[[3]]
disparity_full_ran_slater<-slat_tmp3[[3]]
disparity_full_acc_slater<-slat_tmp4[[3]]
disparity_full_del_slater<-slat_tmp5[[3]]
disparity_full_pro_slater<-slat_tmp6[[3]]
#beck
disparity_full_int_beck  <-beck_tmp1[[3]]
disparity_full_ino_beck  <-beck_tmp2[[3]]
disparity_full_ran_beck  <-beck_tmp3[[3]]
disparity_full_acc_beck  <-beck_tmp4[[3]]
disparity_full_del_beck  <-beck_tmp5[[3]]
disparity_full_pro_beck  <-beck_tmp6[[3]]

#Diversity
#Slice nodes
slat_div_nod<-get(load(paste(data_path, chain_name[1], "/", chain_name[1], diversity_ful, sep="")))
beck_div_nod<-get(load(paste(data_path, chain_name[2], "/", chain_name[2], diversity_ful, sep="")))

#Intervals tips
load(paste(data_path, chain_name[1], '/', chain_name[1], '_distance-tips.Rda', sep=''))
trimmed_max_data_tips<-TrimMorphDistMatrix(dist_tips$gower.dist.matrix)
tree_tips<-drop.tip(tree_slater, trimmed_max_data_tips$removed.taxa) ; tree_tips$root.time<-max(tree.age(tree_tips)[,1])
pco_data_tips<-cmdscale(trimmed_max_data_tips$dist.matrix, k=nrow(trimmed_max_data_tips$dist.matrix) - 2, add=T)$points
slat_div_int<-int.pco(pco_data_tips, tree_tips, intervals, include.nodes=FALSE, FAD_LAD=FADLADslat, diversity=TRUE)[[2]]

load(paste(data_path, chain_name[2], '/', chain_name[2], '_distance-tips.Rda', sep=''))
trimmed_max_data_tips<-TrimMorphDistMatrix(dist_tips$gower.dist.matrix)
tree_tips<-drop.tip(tree_beck  , trimmed_max_data_tips$removed.taxa) ; tree_tips$root.time<-max(tree.age(tree_tips)[,1])
pco_data_tips<-cmdscale(trimmed_max_data_tips$dist.matrix, k=nrow(trimmed_max_data_tips$dist.matrix) - 2, add=T)$points
beck_div_int<-int.pco(pco_data_tips, tree_tips, intervals, include.nodes=FALSE, FAD_LAD=FADLADbeck, diversity=TRUE)[[2]]

#Intervals nodes
load(paste(data_path, chain_name[1], '/', chain_name[1], '_distance-nodes95.Rda', sep=''))
trimmed_max_data_tips<-TrimMorphDistMatrix(dist_tips$gower.dist.matrix)
tree_tips<-drop.tip(tree_slater, trimmed_max_data_tips$removed.taxa) ; tree_tips$root.time<-max(tree.age(tree_tips)[,1])
pco_data_tips<-cmdscale(trimmed_max_data_tips$dist.matrix, k=nrow(trimmed_max_data_tips$dist.matrix) - 2, add=T)$points
slat_div_ino<-int.pco(pco_data_tips, tree_tips, intervals, include.nodes=FALSE, FAD_LAD=FADLADslat, diversity=TRUE)[[2]]

load(paste(data_path, chain_name[2], '/', chain_name[2], '_distance-nodes95.Rda', sep=''))
trimmed_max_data_tips<-TrimMorphDistMatrix(dist_tips$gower.dist.matrix)
tree_tips<-drop.tip(tree_beck  , trimmed_max_data_tips$removed.taxa) ; tree_tips$root.time<-max(tree.age(tree_tips)[,1])
pco_data_tips<-cmdscale(trimmed_max_data_tips$dist.matrix, k=nrow(trimmed_max_data_tips$dist.matrix) - 2, add=T)$points
beck_div_ino<-int.pco(pco_data_tips, tree_tips, intervals, include.nodes=FALSE, FAD_LAD=FADLADbeck, diversity=TRUE)[[2]]


######################################
#Rarefaction analysis on beck
######################################

#Gradual
#Beck
dis_ran_max_beck  <-extract.disp(disparity_full_ran_beck  $quantiles, rarefaction="max")
dis_pro_max_beck  <-extract.disp(disparity_full_pro_beck  $quantiles, rarefaction="max")
dis_ran_min_beck  <-extract.disp(disparity_full_ran_beck  $quantiles, rarefaction="min")
dis_pro_min_beck  <-extract.disp(disparity_full_pro_beck  $quantiles, rarefaction="min")
dis_ran_mod_beck  <-extract.disp(disparity_full_ran_beck  $quantiles, rarefaction=8)
dis_pro_mod_beck  <-extract.disp(disparity_full_pro_beck  $quantiles, rarefaction=8)

#Plot
op<-par(mfrow=c(3,2), bty="n", mar=c(4,4,4,4))
plot.disparity(dis_ran_max_beck  , diversity=dis_ran_max_beck  $rarefaction, main="Punctuated (all taxa)", xlab="Time (Mya)", y2lab="", ylab="Distance from centroid")
abline(v=22, col="red")
plot.disparity(dis_pro_max_beck  , diversity=dis_pro_max_beck  $rarefaction, main="Gradual (all taxa)", xlab="Time (Mya)", y2lab="Species richness", ylab="")
abline(v=22, col="red")

plot.disparity(dis_ran_mod_beck  , diversity=dis_ran_mod_beck  $rarefaction, main="Punctuated (8 taxa)", xlab="Time (Mya)", y2lab="", ylab="Distance from centroid")
abline(v=22, col="red")
plot.disparity(dis_pro_mod_beck  , diversity=dis_pro_mod_beck  $rarefaction, main="Gradual (8 taxa)", xlab="Time (Mya)", y2lab="Species richness", ylab="")
abline(v=22, col="red")

plot.disparity(dis_ran_min_beck  , diversity=dis_ran_min_beck  $rarefaction, main="Punctuated (3 taxa)", xlab="Time (Mya)", y2lab="", ylab="Distance from centroid")
abline(v=22, col="red")
plot.disparity(dis_pro_min_beck  , diversity=dis_pro_min_beck  $rarefaction, main="Gradual (3 taxa)", xlab="Time (Mya)", y2lab="Species richness", ylab="")
abline(v=22, col="red")
par(op)


######################################
#Evolutionary models variation + metrics variation
######################################

#-------------------
#SLATER
#-------------------

dis_tips_slater <-extract.disp(disparity_full_int_slater$quantiles, rarefaction="max")
dis_nodes_slater<-extract.disp(disparity_full_ino_slater$quantiles, rarefaction="max")
dis_ran_slater  <-extract.disp(disparity_full_ran_slater$quantiles, rarefaction="max")
dis_acc_slater  <-extract.disp(disparity_full_acc_slater$quantiles, rarefaction="max")
dis_del_slater  <-extract.disp(disparity_full_del_slater$quantiles, rarefaction="max")
dis_pro_slater  <-extract.disp(disparity_full_pro_slater$quantiles, rarefaction="max")


#The following is a "BIG" plot comparing all methods/metrics for Slater
quartz(width = 15.6, height = 11.2) #A5 landscape
#Windows dimensions
op<-par(mfrow=c(5, 6), bty="n", mar=c(4,4,4,4))# oma=c(bottom, left, top, right)
#Centroid
plot.disparity(dis_tips_slater, xlab="", ylab="Distance from centroid", measure="Cent.dist", main="Intervals (tips only)", diversity=dis_tips_slater$rarefaction, y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_nodes_slater, xlab="", ylab="", measure="Cent.dist", main="Intervals (tips and nodes)", diversity=dis_nodes_slater$rarefaction, y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_ran_slater, xlab="", ylab="", measure="Cent.dist", main="Slices (punctuated)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_slater, xlab="", ylab="", measure="Cent.dist", main="Slices (punctuated:acctran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_slater, xlab="", ylab="", measure="Cent.dist", main="Slices (punctuated:deltran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_slater, xlab="", ylab="", measure="Cent.dist", main="Slices (gradual)", diversity=slat_div_nod, y2lab="Species richness")
abline(v= 22, col="red")

#Sum of ranges
plot.disparity(dis_tips_slater, xlab="", ylab="Sum of ranges", measure="Sum.range", main="Intervals (tips only)", diversity=dis_tips_slater$rarefaction, y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_nodes_slater, xlab="", ylab="", measure="Sum.range", main="Intervals (tips and nodes)", diversity=dis_nodes_slater$rarefaction, y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_ran_slater, xlab="", ylab="", measure="Sum.range", main="Slices (punctuated)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_slater, xlab="", ylab="", measure="Sum.range", main="Slices (punctuated:acctran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_slater, xlab="", ylab="", measure="Sum.range", main="Slices (punctuated:deltran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_slater, xlab="", ylab="", measure="Sum.range", main="Slices (gradual)", diversity=slat_div_nod, y2lab="Species richness")
abline(v= 22, col="red")

#Sum of variance
plot.disparity(dis_tips_slater, xlab="", ylab="Sum of variances", measure="Sum.var", main="Intervals (tips only)", diversity=dis_tips_slater$rarefaction, y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_nodes_slater, xlab="", ylab="", measure="Sum.var", main="Intervals (tips and nodes)", diversity=dis_nodes_slater$rarefaction, y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_ran_slater, xlab="", ylab="", measure="Sum.var", main="Slices (punctuated)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_slater, xlab="", ylab="", measure="Sum.var", main="Slices (punctuated:acctran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_slater, xlab="", ylab="", measure="Sum.var", main="Slices (punctuated:deltran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_slater, xlab="", ylab="", measure="Sum.var", main="Slices (gradual)", diversity=slat_div_nod, y2lab="Species richness")
abline(v= 22, col="red")

#Product of ranges
plot.disparity(dis_tips_slater, xlab="", ylab="Product of ranges", measure="Prod.range", main="Intervals (tips only)", diversity=dis_tips_slater$rarefaction, y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_nodes_slater, xlab="", ylab="", measure="Prod.range", main="Intervals (tips and nodes)", diversity=dis_nodes_slater$rarefaction, y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_ran_slater, xlab="", ylab="", measure="Prod.range", main="Slices (punctuated)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_slater, xlab="", ylab="", measure="Prod.range", main="Slices (punctuated:acctran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_slater, xlab="", ylab="", measure="Prod.range", main="Slices (punctuated:deltran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_slater, xlab="", ylab="", measure="Prod.range", main="Slices (gradual)", diversity=slat_div_nod, y2lab="Species richness")
abline(v= 22, col="red")

#Product of variance
plot.disparity(dis_tips_slater, xlab="", ylab="Product of variances", measure="Prod.var", main="Intervals (tips only)", diversity=dis_tips_slater$rarefaction, y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_nodes_slater, xlab="", ylab="", measure="Prod.var", main="Intervals (tips and nodes)", diversity=dis_nodes_slater$rarefaction, y2lab="", cex.xaxis=0.7)
abline(v= 17.5, col="red")
plot.disparity(dis_ran_slater, xlab="Time (Mya)", ylab="", measure="Prod.var", main="Slices (punctuated)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_slater, xlab="Time (Mya)", ylab="", measure="Prod.var", main="Slices (punctuated:acctran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_slater, xlab="Time (Mya)", ylab="", measure="Prod.var", main="Slices (punctuated:deltran)", diversity=slat_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_slater, xlab="Time (Mya)", ylab="", measure="Prod.var", main="Slices (gradual)", diversity=slat_div_nod, y2lab="Species richness")
abline(v= 22, col="red")

par(op)

#-------------------
#BECK
#-------------------


dis_tips_beck <-extract.disp(disparity_full_int_beck$quantiles, rarefaction="max")
dis_nodes_beck<-extract.disp(disparity_full_ino_beck$quantiles, rarefaction="max")
dis_ran_beck  <-extract.disp(disparity_full_ran_beck$quantiles, rarefaction="max")
dis_acc_beck  <-extract.disp(disparity_full_acc_beck$quantiles, rarefaction="max")
dis_del_beck  <-extract.disp(disparity_full_del_beck$quantiles, rarefaction="max")
dis_pro_beck  <-extract.disp(disparity_full_pro_beck$quantiles, rarefaction="max")


#The following is a "BIG" plot comparing all methods/metrics for Slater
quartz(width = 15.6, height = 11.2) #A5 landscape
#Windows dimensions
op<-par(mfrow=c(5, 6), bty="n", mar=c(4,4,4,4))# oma=c(bottom, left, top, right)
#Centroid
plot.disparity(dis_tips_beck, xlab="", ylab="Distance from centroid", measure="Cent.dist", main="Intervals (tips only)", diversity=dis_tips_beck$rarefaction, y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_nodes_beck, xlab="", ylab="", measure="Cent.dist", main="Intervals (tips and nodes)", diversity=dis_nodes_beck$rarefaction, y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_ran_beck, xlab="", ylab="", measure="Cent.dist", main="Slices (punctuated)", diversity=beck_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_beck, xlab="", ylab="", measure="Cent.dist", main="Slices (punctuated:acctran)", diversity=beck_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_beck, xlab="", ylab="", measure="Cent.dist", main="Slices (punctuated:deltran)", diversity=beck_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_beck, xlab="", ylab="", measure="Cent.dist", main="Slices (gradual)", diversity=beck_div_nod, y2lab="Species richness")
abline(v= 22, col="red")

#Sum of ranges
plot.disparity(dis_tips_beck, xlab="", ylab="Sum of ranges", measure="Sum.range", main="Intervals (tips only)", diversity=dis_tips_beck$rarefaction, y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_nodes_beck, xlab="", ylab="", measure="Sum.range", main="Intervals (tips and nodes)", diversity=dis_nodes_beck$rarefaction, y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_ran_beck, xlab="", ylab="", measure="Sum.range", main="Slices (punctuated)", diversity=beck_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_beck, xlab="", ylab="", measure="Sum.range", main="Slices (punctuated:acctran)", diversity=beck_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_beck, xlab="", ylab="", measure="Sum.range", main="Slices (punctuated:deltran)", diversity=beck_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_beck, xlab="", ylab="", measure="Sum.range", main="Slices (gradual)", diversity=beck_div_nod, y2lab="Species richness")
abline(v= 22, col="red")

#Sum of variance
plot.disparity(dis_tips_beck, xlab="", ylab="Sum of variances", measure="Sum.var", main="Intervals (tips only)", diversity=dis_tips_beck$rarefaction, y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_nodes_beck, xlab="", ylab="", measure="Sum.var", main="Intervals (tips and nodes)", diversity=dis_nodes_beck$rarefaction, y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_ran_beck, xlab="", ylab="", measure="Sum.var", main="Slices (punctuated)", diversity=beck_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_beck, xlab="", ylab="", measure="Sum.var", main="Slices (punctuated:acctran)", diversity=beck_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_beck, xlab="", ylab="", measure="Sum.var", main="Slices (punctuated:deltran)", diversity=beck_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_beck, xlab="", ylab="", measure="Sum.var", main="Slices (gradual)", diversity=beck_div_nod, y2lab="Species richness")
abline(v= 22, col="red")

#Product of ranges
plot.disparity(dis_tips_beck, xlab="", ylab="Product of ranges", measure="Prod.range", main="Intervals (tips only)", diversity=dis_tips_beck$rarefaction, y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_nodes_beck, xlab="", ylab="", measure="Prod.range", main="Intervals (tips and nodes)", diversity=dis_nodes_beck$rarefaction, y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_ran_beck, xlab="", ylab="", measure="Prod.range", main="Slices (punctuated)", diversity=beck_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_beck, xlab="", ylab="", measure="Prod.range", main="Slices (punctuated:acctran)", diversity=beck_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_beck, xlab="", ylab="", measure="Prod.range", main="Slices (punctuated:deltran)", diversity=beck_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_beck, xlab="", ylab="", measure="Prod.range", main="Slices (gradual)", diversity=beck_div_nod, y2lab="Species richness")
abline(v= 22, col="red")

#Product of variance
plot.disparity(dis_tips_beck, xlab="", ylab="Product of variances", measure="Prod.var", main="Intervals (tips only)", diversity=dis_tips_beck$rarefaction, y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_nodes_beck, xlab="", ylab="", measure="Prod.var", main="Intervals (tips and nodes)", diversity=dis_nodes_beck$rarefaction, y2lab="", cex.xaxis=0.8)
abline(v= 10.5, col="red")
plot.disparity(dis_ran_beck, xlab="Time (Mya)", ylab="", measure="Prod.var", main="Slices (punctuated)", diversity=beck_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_acc_beck, xlab="Time (Mya)", ylab="", measure="Prod.var", main="Slices (punctuated:acctran)", diversity=beck_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_del_beck, xlab="Time (Mya)", ylab="", measure="Prod.var", main="Slices (punctuated:deltran)", diversity=beck_div_nod, y2lab="")
abline(v= 22, col="red")
plot.disparity(dis_pro_beck, xlab="Time (Mya)", ylab="", measure="Prod.var", main="Slices (gradual)", diversity=beck_div_nod, y2lab="Species richness")
abline(v= 22, col="red")

par(op)


######################################
#Wills 1994 measures
######################################

#Slater
op<-par(mfrow=c(2,2), bty="n", mar=c(4,4,4,4))
plot.disparity(dis_ran_max_slater, diversity=dis_ran_max_slater$rarefaction, measure="Sum.range", main="Sum of ranges", xlab="Time (Mya)", y2lab="")
abline(v=22, col="red")
plot.disparity(dis_ran_max_slater, diversity=dis_ran_max_slater$rarefaction, measure="Sum.var", main="Sum of variance", xlab="Time (Mya)", y2lab="Species richness")
abline(v=22, col="red")
plot.disparity(dis_ran_max_slater, diversity=dis_ran_max_slater$rarefaction, measure="Prod.range", main="Product of ranges", xlab="Time (Mya)", y2lab="")
abline(v=22, col="red")
plot.disparity(dis_ran_max_slater, diversity=dis_ran_max_slater$rarefaction, measure="Prod.var", main="Product of variance", xlab="Time (Mya)", y2lab="Species richness")
abline(v=22, col="red")
par(op)

#beck
op<-par(mfrow=c(2,2), bty="n", mar=c(4,4,4,4))
plot.disparity(dis_ran_max_beck, diversity=dis_ran_max_beck$rarefaction, measure="Sum.range", main="Sum of ranges", xlab="Time (Mya)", y2lab="")
abline(v=22, col="red")
plot.disparity(dis_ran_max_beck, diversity=dis_ran_max_beck$rarefaction, measure="Sum.var", main="Sum of variance", xlab="Time (Mya)", y2lab="Species richness")
abline(v=22, col="red")
plot.disparity(dis_ran_max_beck, diversity=dis_ran_max_beck$rarefaction, measure="Prod.range", main="Product of ranges", xlab="Time (Mya)", y2lab="")
abline(v=22, col="red")
plot.disparity(dis_ran_max_beck, diversity=dis_ran_max_beck$rarefaction, measure="Prod.var", main="Product of variance", xlab="Time (Mya)", y2lab="Species richness")
abline(v=22, col="red")
par(op)


#####################################
#
