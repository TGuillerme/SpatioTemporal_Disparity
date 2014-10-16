#Data input (Slater)
#Slater.table<-read.table("https://raw.githubusercontent.com/TGuillerme/SpatioTemporal_Disparity/master/Data/2013-Slater-MEE-morpho.table", header=F, sep=" ", row.names=1) 
Slater.table<-read.table("../Data/2013-Slater-MEE-morpho.table", header=F, sep=" ", row.names=1) 
#Slater.tree<-read.nexus('https://raw.githubusercontent.com/TGuillerme/SpatioTemporal_Disparity/master/Data/2013-Slater-MEE-TEM.tre')
Slater.tree<-read.nexus('../Data/2013-Slater-MEE-TEM.tre')

#Remove species with only missing data before hand
Slater.table<-Slater.table[-c(as.vector(which(apply(as.matrix(Slater.table), 1, function(x) levels(as.factor(x))) == "?"))),]
#Slater.table<-Slater.table[-c(grep("Aegialodon", row.names(Slater.table)), grep("Murtoilestes", row.names(Slater.table)))] #Somehow Aegialodon and Murtoilestes are bugged and Adelobasileus seems to produce NAs in the distance matrix

#Cleaning the tree and the table
tree<-clean.tree(Slater.tree, Slater.table)
table<-clean.table(Slater.table, Slater.tree)
cat("Created the morphological characters table as:\ntable\n")

#Making the tree binary
tree<-bin.tree(tree)
#adding node names
tree$node.label<-paste("n",seq(1:Nnode(tree)), sep="")
cat("Created the tree as:\ntree\n")

#test
expect_equal(Ntip(tree), nrow(table))

#Characters list
Mandible<-seq(1:39)
Premolars<-seq(40:56)
Molar_morpho<-seq(57:126)
Molar_wear<-seq(127:141)
Dental_other<-seq(142:173)
Vertebrae<-seq(174:185)
Shoulder<-seq(186:208)
Forelimb<-seq(209:224)
Pelvic<-seq(225:237)
Hindlimb<-seq(238:287)
Postcranial_other<-seq(288:292)
Basicranium<-seq(293:365)
Middle_ear<-seq(366:387)
Cranial_other<-seq(388:436)
Cranial_vault<-seq(437:443)
Soft_tissue<-seq(444:445)
character.details<-c("Mandible","Premolars","Molar_morpho","Molar_wear","Dental_other","Vertebrae","Shoulder","Forelimb","Pelvic","Hindlimb","Postcranial_other","Basicranium","Middle_ear","Cranial_other","Cranial_vault","Soft_tissue")

#Groups
Dental<-c(Premolars, Molar_morpho, Molar_wear, Dental_other)
Cranial<-c(Mandible, Basicranium, Middle_ear, Cranial_other, Cranial_vault)
PostCranial<-c(Vertebrae, Shoulder, Forelimb, Pelvic, Hindlimb, Postcranial_other)
cat("Created the different morphological characters categories as:\nDental\nCranial\nPostCranial\nCharacters details list is available in the object:\ncharacter.details\n")

#taxonomy.list
taxonomy.list<-list("Australosphenids"=164, "Marsupialiomorphs"=125, "Placentaliomorphs"=108, "Stem_mammaliforms"=c(93, 103), "Stem_theriforms"=c(104,107)) #Slater
cat("Created the taxonomical list as:\ntaxonomy.list\n")

#Age slices
slices<-c(0,25,50,75,100,125,150,175,200)
cat("Created the age slices list as:\nslices\n")

#load("https://raw.githubusercontent.com/TGuillerme/SpatioTemporal_Disparity/master/Data/slater.ML.mat.rda")
load("../Data/slater.ML.mat.rda")
anc.matrix.save<-anc.matrix.slater.ML
cat("Created the ancestral matrix in \'ML\' as:\nanc.matrix.save\n")

