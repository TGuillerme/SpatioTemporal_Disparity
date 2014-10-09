#TEST ANC STAT

#functions

#Testing pco.std
#example
tree<-rtree(10)
matrix<-data.frame(row.names=tree$tip.label, 'char1'=sample(c(0,1), Ntip(tree), replace=TRUE), 'char2'=sample(c(0,1), Ntip(tree), replace=TRUE))
anc.matrix<-anc.state(tree, matrix, method='ML', verbose=FALSE)

#test
#pco output
expect_is(pco.std(anc.matrix), 'pcoa') ; message('.', appendLF=FALSE)
expect_equal(names(pco.std(anc.matrix)), c("correction","note","values","vectors","trace")) ; message('.', appendLF=FALSE)
#errors (matrix)
expect_error(pco.std(c(1,2))) ; message('.', appendLF=FALSE)

#errors (vegdist)
expect_error(pco.std(anc.matrix, distance="far away")) ; message('.', appendLF=FALSE)
expect_error(pco.std(anc.matrix, na.rm="yes please")) ; message('.', appendLF=FALSE)
#options
expect_is(pco.std(anc.matrix, diag=TRUE), 'pcoa') ; message('.', appendLF=FALSE)

#errors (pcoa)
expect_error(pco.std(anc.matrix, correction="yes please")) ; message('.', appendLF=FALSE)
