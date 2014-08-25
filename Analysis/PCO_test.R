#Header
setwd('~/PhD/Projects/SpatioTemporal_Disparity/Analysis')
library(ape)

#data
tree<-read.tree('BDtree.tre')
table<-read.dna('MorphoMat.phylip', as.character=TRUE, as.matrix=TRUE)

#Euclidian matrix
eucl.table<-dist(table, method = "euclidean")

#pcoa
pco<-pcoa(eucl.table)

#plot
pco$values
biplot(pco)
#or nicer
clade1<-extract.clade(tree, 88)$tip.label
clade2<-extract.clade(tree, 54)$tip.label
pco.scores<-as.data.frame(pco$vectors)
clade.col<-ncol(pco.scores)+1
pco.scores[, clade.col]<-"NA"
pco.scores[match(clade1, row.names(pco.scores)),clade.col]<-"clade1"
pco.scores[match(clade2, row.names(pco.scores)),clade.col]<-"clade2"
plot(pco$vectors[,1], pco$vectors[,2], xlab="PC1", ylab="PC2", col=c("red", "blue")[as.factor(pco.scores[,clade.col])])


#Really dummy data
species<-LETTERS[1:10]
t1<-c(rep(0,5), rep(1,5))
t2<-c(rep(0,6), rep(1,4))
t3<-c(rep(0,7), rep(1,3))
t4<-c(rep(0,8), rep(1,2))
t5<-c(rep(0,9), rep(1,1))
t6<-c(rep(1,5), rep(0,5))
t7<-c(rep(1,6), rep(0,4))
t8<-c(rep(1,7), rep(0,3))
t9<-c(rep(1,8), rep(0,2))
t10<-c(rep(1,9), rep(0,1))
table<-data.frame(row.names=species, c1=rep(NA,10), c2=rep(NA,10), c3=rep(NA,10), c4=rep(NA,10), c5=rep(NA,10), c6=rep(NA,10), c7=rep(NA,10), c8=rep(NA,10), c9=rep(NA,10), c10=rep(NA,10))
table[1,]<-t1
table[2,]<-t2
table[3,]<-t3
table[4,]<-t4
table[5,]<-t5
table[6,]<-t6
table[7,]<-t7
table[8,]<-t8
table[9,]<-t9
table[10,]<-t10

#Euclidian matrix
eucl.table<-dist(table, method = "euclidean")

#pcoa
pco<-pcoa(eucl.table)

#plot
pco$values
biplot(pco)

clade1<-LETTERS[1:5]
clade2<-LETTERS[6:10]
pco.scores<-as.data.frame(pco$vectors)
clade.col<-ncol(pco.scores)+1
pco.scores[, clade.col]<-"NA"
pco.scores[match(clade1, row.names(pco.scores)),clade.col]<-"clade1"
pco.scores[match(clade2, row.names(pco.scores)),clade.col]<-"clade2"
plot(pco$vectors[,1], pco$vectors[,2], xlab="PC1", ylab="PC2", col=c("red", "blue")[as.factor(pco.scores[,clade.col])])

#Not less dummy
species<-LETTERS[1:10]
t1<-c(rep(1,1), rep(0,5), rep(1,4))
t2<-c(rep(1,2), rep(0,5), rep(1,3))
t3<-c(rep(1,3), rep(0,5), rep(1,2))
t4<-c(rep(1,4), rep(0,5), rep(1,1))
t5<-c(rep(1,5), rep(0,5))
table<-data.frame(row.names=species, t1=t1, t2=t2, t3=t3, t4=t4, t5=t5)

#Euclidian matrix
eucl.table<-dist(table, method = "euclidean")

#pcoa
pco<-pcoa(eucl.table)

#plot
pco$values
biplot(pco)

clade1<-LETTERS[1:5]
clade2<-LETTERS[6:10]
pco.scores<-as.data.frame(pco$vectors)
clade.col<-ncol(pco.scores)+1
pco.scores[, clade.col]<-"NA"
pco.scores[match(clade1, row.names(pco.scores)),clade.col]<-"clade1"
pco.scores[match(clade2, row.names(pco.scores)),clade.col]<-"clade2"
plot(pco$vectors[,1], pco$vectors[,2], xlab="PC1", ylab="PC2", col=c("red", "blue")[as.factor(pco.scores[,clade.col])])
