#TEST ANC.UNC

#functions

#Testing anc.unc
#example
tree<-rtree(10)
matrix.state<-data.frame(row.names=tree$tip.label, 'char1'=sample(c(0,1), Ntip(tree), replace=TRUE), 'char2'=sample(c(0,1), Ntip(tree), replace=TRUE))
set.seed(10)
#prob matrix with only two cells == 1
matrix.prob<-data.frame(row.names=tree$tip.label, "char1"=rnorm(10), "char2"=rnorm(10))
matrix.prob=matrix.prob^2 ; matrix.prob[,1]<-matrix.prob[,1]/max(matrix.prob[,1]) ; matrix.prob[,2]<-matrix.prob[,2]/max(matrix.prob[,2])
anc.mat<-list("state"=matrix.state, "prob"=matrix.prob)
test<-anc.unc(anc.mat, threshold=1)

#test
test_that('Testing anc.state_ace()', {
    #anc.mat
    expect_is(anc.mat, 'list')
    expect_that(length(anc.mat), equals(2))
    expect_that(names(anc.mat), equals(c("state", "prob")))

    #object
    expect_is(test, 'list')
    expect_that(length(test), equals(2))
    expect_that(names(test), equals(c("state", "prob")))

    #only two cells in the right threshold
    expect_equal(which(is.numeric(suppressWarnings(as.numeric(test$state[,1])))), 1)
    expect_true(any(match("?", test$state[,1])))
    expect_equal(which(is.numeric(suppressWarnings(as.numeric(test$state[,2])))), 1)
    expect_true(any(match("?", test$state[,2])))

    #errors
    expect_error(anc.unc(anc.mat, threshold='YES'))
    expect_error(anc.unc(anc.mat, threshold=1.8))
    expect_error(anc.unc("hahaha"))
})
message('.', appendLF=FALSE)
