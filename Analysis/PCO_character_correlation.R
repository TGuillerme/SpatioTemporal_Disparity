##########################
#CHARACTER CORRELATION ANALYSIS
##########################

#Correlation matrix?
matrix<-table #Not on the estimated characters!
#NA removal from pco.std()
for (character in 1:ncol(matrix)) {
    for (taxa in 1:nrow(matrix)) {
        if(as.character(matrix[taxa, character]) == "?") {
            matrix[taxa, character] <- NA
        }
    }
}
#Transforming the matrix to numeric
matrix<-as.data.frame(matrix)
for (character in 1:ncol(matrix)) {
    matrix[,character]<-as.numeric(matrix[,character])
}

#Scaling
scaled.matrix<-scale(matrix, scale=TRUE, center=TRUE)
#Calculating the distance matrix
distance.matrix<-vegdist(scaled.matrix, method="euclidean", na.rm=TRUE)

#correlation in the numeric matrix
cor.mat<-cor(matrix[], use = "pairwise")
diag(cor.mat)<-NA

tabNAs<-length(which(is.na(matrix)))/(ncol(matrix)*nrow(matrix)) # 50% NAs
corNAs<-length(which(is.na(cor.mat)))/length(cor.mat) # 5% NAs
hist(cor.mat)

#Isolating completely correlated characters
#Positive
cor.pairs.pos<-list()
for(character in 1:ncol(cor.mat)) {
    if(length(which(cor.mat[,character]==1)) == 0){
        cor.pairs.pos[[character]]<-NA
    } else {
       cor.pairs.pos[[character]]<-which(cor.mat[,character]==1) 
    }
}
names(cor.pairs.pos)<-paste("V",seq(1:ncol(cor.mat)), sep="")
cor.pairs.pos<-cor.pairs.pos[-which(is.na(cor.pairs.pos))]
#Negative
cor.pairs.neg<-list()
for(character in 1:ncol(cor.mat)) {
    if(length(which(cor.mat[,character]==-1)) == 0){
        cor.pairs.neg[[character]]<-NA
    } else {
       cor.pairs.neg[[character]]<-which(cor.mat[,character]==-1) 
    }
}
names(cor.pairs.neg)<-paste("V",seq(1:ncol(cor.mat)), sep="")
cor.pairs.neg<-cor.pairs.neg[-which(is.na(cor.pairs.neg))]

#Visualise the correlated characters
hist(unlist(cor.pairs.pos), breaks=421, main="positively correlated characters")
hist(unlist(cor.pairs.neg), breaks=421, main="negatively correlated characters")

image(cor.mat) #hooooooo