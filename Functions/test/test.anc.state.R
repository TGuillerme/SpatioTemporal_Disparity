#TEST ANC STAT

#functions

#Testing anc.state_ace
#example
tree<-rtree(10)
matrix<-data.frame(row.names=tree$tip.label, 'char1'=sample(c(0,1), Ntip(tree), replace=TRUE), 'char2'=sample(c(0,1), Ntip(tree), replace=TRUE))
test<-anc.state_ace(tree, matrix, method='ML', verbose=FALSE)

#test
test_that('Testing anc.state_ace()', {
    #class
    expect_is(test, 'list')
    expect_is(test[[1]], 'ace')
    expect_is(test[[2]], 'ace')
    #number of reconstructions
    expect_that(length(test), equals(ncol(matrix)))
    #verbose
    expect_message(bla<-anc.state_ace(tree, matrix, method='ML', verbose=TRUE), 'Estimating the ancestral states for 2 characters:')
    expect_message(bla<-anc.state_ace(tree, matrix, method='ML', verbose=TRUE), '.')
    expect_message(bla<-anc.state_ace(tree, matrix, method='ML', verbose=TRUE), 'Done.')
})
message('.', appendLF=FALSE)

#Testing anc.state_prob
#example
tree<-rtree(10)
matrix<-data.frame(row.names=tree$tip.label, 'char1'=sample(c(0,1), Ntip(tree), replace=TRUE), 'char2'=sample(c(0,1), Ntip(tree), replace=TRUE))
test<-anc.state_prob(tree, matrix, anc.state_ace(tree, matrix, method='ML', verbose=FALSE))

#test
test_that('Testing anc.state_prob()', {
    #class
    expect_is(test, 'matrix')
    #number of columns
    expect_equal(ncol(test), ncol(matrix))
    #number of rows
    expect_equal(nrow(test), (Ntip(tree)+Nnode(tree)))
    #probabilities
    expect_that(max(test), is_less_than(1.000001))
    expect_that(min(test), is_more_than(-0.000001))
})
message('.', appendLF=FALSE)


#Testing anc.state_state
tree<-rtree(10)
matrix<-data.frame(row.names=tree$tip.label, 'char1'=sample(c(0,1), Ntip(tree), replace=TRUE), 'char2'=sample(c(0,1), Ntip(tree), replace=TRUE))
test<-anc.state_state(tree, matrix, anc.state_ace(tree, matrix, method='ML', verbose=FALSE))

#test
test_that('Testing anc.state_state()', {
    #class
    expect_is(test, 'matrix')
    #number of columns
    expect_equal(ncol(test), ncol(matrix))
    #number of rows
    expect_equal(nrow(test), (Ntip(tree)+Nnode(tree)))
    #states
    expect_is(test[1,1], 'character')
})
message('.', appendLF=FALSE)


#Testing anc.states
tree<-rtree(10)
matrix<-data.frame(row.names=tree$tip.label, 'char1'=sample(c(0,1), Ntip(tree), replace=TRUE), 'char2'=sample(c(0,1), Ntip(tree), replace=TRUE))
test<-anc.state(tree, matrix, method='ML', verbose=FALSE)

#test
test_that('Testing anc.state()', {
    #class
    expect_is(test, 'list')
    #length
    expect_equal(length(test), 2)
    #elements
    expect_is(test[[1]], 'matrix')
    expect_is(test[[2]], 'matrix')
    #names
    expect_equal(names(test), c('state', 'prob'))
    #correct levels
    for(character in 1:ncol(matrix)) {
        expect_that(levels(as.factor(matrix[,character])) %in% as.vector(test$state[,character]), equals(rep(TRUE, length(levels(as.factor(matrix[,character]))))))
    }
    #verbose
    expect_message(bla<-anc.state(tree, matrix, method='ML', verbose=TRUE), 'Estimating the ancestral states for 2 characters:')
    expect_message(bla<-anc.state(tree, matrix, method='ML', verbose=TRUE), '.')
    expect_message(bla<-anc.state(tree, matrix, method='ML', verbose=TRUE), 'Done.')
})
message('.', appendLF=FALSE)
