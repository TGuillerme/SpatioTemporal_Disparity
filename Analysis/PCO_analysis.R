#Header
if(grep("TGuillerme", getwd())) {
    setwd('~/PhD/Projects/SpatioTemporal_Disparity/Analysis')
} else {
    warning("You might have to change the directory!")
}
source("functions.R")

#data='Slater'
data='Beck'

if (data == 'Slater') {
    #Data input (Slater)
    Slater.table<-read.table("../Data/2013-Slater-MEE-morpho.table", header=F, sep=" ", row.names=1) 
    Slater.tree<-read.nexus('../Data/2013-Slater-MEE-TEM.tre')

    #Remove species with only missing data before hand
    if (any(apply(as.matrix(Beck.table), 1, function(x) levels(as.factor((x)))) == "?")) {
        Slater.table<-Slater.table[-c(as.vector(which(apply(as.matrix(Slater.table), 1, function(x) levels(as.factor(x))) == "?"))),]
    }
    Slater.table<-Slater.table[-c(grep("Aegialodon", row.names(Slater.table)), grep("Murtoilestes", row.names(Slater.table))),] #Aegialodon and Murtoilestes are bugged

    #Cleaning the tree and the table
    tree<-clean.tree(Slater.tree, Slater.table)
    table<-clean.table(Slater.table, Slater.tree)
    #Making the tree binary
    tree<-bin.tree(tree)
    #adding node names
    tree$node.label<-paste("n",seq(1:Nnode(tree)), sep="")

    #test
    expect_equal(Ntip(tree), nrow(table))
    suppressWarnings({eucl.table<-dist(table, method = "euclidean")})
    pco<-pcoa(eucl.table)
    expect_is(pco, "pcoa")
}

if (data == 'Beck')
    #Data input (Beck)
    Beck.table<-read.table("../Data/2014-Beck-ProcB-morpho.table", header=F, sep=" ", row.names=1) 
    Beck.tree<-read.nexus('../Data/2014-Beck-ProcB-TEM.tre')

    #Remove species with only missing data before hand
    if (any(apply(as.matrix(Beck.table), 1, function(x) levels(as.factor((x)))) == "?")) {
        Beck.table<-Beck.table[-c(as.vector(which(apply(as.matrix(Beck.table), 1, function(x) levels(as.factor(x))) == "?"))),]
    }

    #Cleaning the tree and the table
    tree<-clean.tree(Beck.tree, Beck.table)
    table<-clean.table(Beck.table, Beck.tree)
    #Making the tree binary
    tree<-bin.tree(tree)
    #adding node names
    tree$node.label<-paste("n",seq(1:Nnode(tree)), sep="")   

    #test
    expect_equal(Ntip(tree), nrow(table))
    suppressWarnings({eucl.table<-dist(table, method = "euclidean")})
    pco<-pcoa(eucl.table)
    expect_is(pco, "pcoa")
}

#Ancestral state matrix
anc.matrix.save<-anc.state(tree, table, model='ML', verbose=TRUE)

#Test: replacing "?" by NAs
#anc.matrix<-anc.matrix.save
#for (character in 1:ncol(anc.matrix$state)) {
#    for (taxa in 1:nrow(anc.matrix$state)) {
#        if(as.character(anc.matrix$state[taxa, character]) == "?") {
#            anc.matrix$state[taxa, character] <- NA
#        }
#    }
#}


#Recalculating the matrix with a 0.95 probability lower limit
anc.matrix<-anc.unc(anc.matrix.save, 0.95)

#PCO
dist.matrix<-dist(anc.matrix$state, method="euclidian")
pco<-pcoa(dist.matrix, correction="lingoes") #or "cailliez", see ?pcoa
pco.scores<-pco$vectors

dat3b<-scale(dat3[-c(218:220),c(3:9)])
s3b<-svd(cor(dat3b))
sum<-s3b$d/sum(s3b$d)
proj3<-dat3b%*%s3b$u


#svd?
dm.svd<-svd(dist.matrix)
dm.svd.sca<-svd(scale(dist.matrix))
dm.svd.sco<-svd(cor(scale(dist.matrix)))
#Cumulative variance
cum.dm.svd<-cumsum(dm.svd$d/sum(dm.svd$d))
cum.dm.svd.sca<-cumsum(dm.svd.sca$d/sum(dm.svd.sca$d))
cum.dm.svd.sco<-cumsum(dm.svd.sco$d/sum(dm.svd.sco$d))

#Checking the axis variance
op<-par(mfrow=c(3,2))
barplot(dm.svd$d/sum(dm.svd$d), main="distance matrix")
hist(cum.dm.svd, main="distance matrix")
barplot(dm.svd.sca$d/sum(dm.svd.sca$d), main="scale(dist.mat)")
hist(cum.dm.svd.sca, main="scale(dist.mat)")
barplot(dm.svd.sco$d/sum(dm.svd.sco$d), main="cor(scale(dist.mat))")
hist(cum.dm.svd.sco, main="cor(scale(dist.mat))")
par(op)


#Clades Slater
#Australosphenids (monotremes and relatives) node 160
#pco.scores<-set.group(tree, pco.scores, type='clade', node=160, name='Australosphenids')
#Marsupialiomorphs (marsupials and relatives) node 122
#pco.scores<-set.group(tree, pco.scores, tax.col="taxonomy", type='clade', node=122, name='Marsupialiomorphs')
#Placentaliomorphs (placentals and relatives) node 106
#pco.scores<-set.group(tree, pco.scores, tax.col="taxonomy", type='clade', node=106, name='Placentaliomorphs')
#Stem mammaliforms node 91 to 101
#pco.scores<-set.group(tree, pco.scores, tax.col="taxonomy", type='grade', node=c(91, 101), name='Stem_mammaliforms')
#Stem theriforms node 102 105
#pco.scores<-set.group(tree, pco.scores, tax.col="taxonomy", type='grade', node=c(102, 105), name='Stem_theriforms')

#Clades Beck
#Placental
pco.scores<-set.group(tree, pco.scores, type='clade', node=153, name='Placental')
#Stem-placental
pco.scores<-set.group(tree, pco.scores, type='grade', node=c(103,153), name='Stem-placental')

#Full pco plot
plot.pco(pco.scores, "taxonomy", main="Entire \"morphospace\"", legend=TRUE)

#lim
xlim=c(min(pco.scores[,1]), max(pco.scores[,1]))
ylim=c(min(pco.scores[,2]), max(pco.scores[,2]))

#Time slices pco plots Beck
pco.slice(tree, pco.scores, c(0, 40, 50, 60, 70, 80, 90, 100, 110), 'ACCTRAN', tax.col="taxonomy", legend=FALSE, pars=c(3,3), xlim=xlim, ylim=ylim)
dev.new()
pco.slice(tree, pco.scores, c(0, 40, 50, 60, 70, 80, 90, 100, 110), 'DELTRAN', tax.col="taxonomy", legend=FALSE, pars=c(3,3), xlim=xlim, ylim=ylim)

#Time slices pco plots Slate
pco.slice(tree, pco.scores, c(0,25,50,75,100,125,150,175,200), 'ACCTRAN', tax.col="taxonomy", legend=FALSE, pars=c(3,3), xlim=xlim, ylim=ylim)
pco.slice(tree, pco.scores, c(0,25,50,75,100,125,150,175,200), 'DELTRAN', tax.col="taxonomy", legend=FALSE, pars=c(3,3), xlim=xlim, ylim=ylim)
