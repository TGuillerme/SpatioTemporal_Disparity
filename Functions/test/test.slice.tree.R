#TEST SLICE TREE

#Testing slice.tree_parent.node
#example
tree<-read.tree(text="(((((A:1,B:1):2,C:3):1,D:1):1,E:5):1,F:3);")
tree$node.label<-as.character(seq(1:5))

#test
test_that("Testing slice.tree_parent.node()", {
    #class
    expect_error(slice.tree_parent.node(tree, '1'))
    expect_is(slice.tree_parent.node(tree, 'B'), 'character')
    #length
    expect_error(length(slice.tree_parent.node(tree, '1')))
    expect_equal(length(slice.tree_parent.node(tree, 'B')), 1)
    #null
    expect_error(slice.tree_parent.node(rtree(5), '1'))
    expect_error(slice.tree_parent.node(rtree(5), 'B'))
    #sister nodes/taxa
    expect_equal(slice.tree_parent.node(tree, 'A'), slice.tree_parent.node(tree, 'B'))
    expect_equal(slice.tree_parent.node(tree, 'E'), slice.tree_parent.node(tree, '3'))
})
message('.', appendLF=FALSE)

#Testing slice.tree_offspring.node
#example
tree<-read.tree(text="(((((A:1,B:1):2,C:3):1,D:1):1,E:5):1,F:3);")
tree$node.label<-as.character(seq(1:5))

#test
test_that("Testing slice.tree_offspring.node()", {
    #error
    expect_error(slice.tree_offspring.node(tree, '1'))
    expect_error(slice.tree_offspring.node(tree, 'E', 'E'))
    expect_error(slice.tree_offspring.node(tree, '5', 'E'))
    expect_error(slice.tree_offspring.node(tree, '1', '2'))
    #class
    expect_is(slice.tree_offspring.node(tree, '1', 'A'), 'character')
    #length
    expect_equal(length(slice.tree_offspring.node(tree, '1', 'A')), 1)
    #sister nodes/taxa
    expect_equal(slice.tree_offspring.node(tree, '1', 'A'), slice.tree_offspring.node(tree, '1', 'E'))
    expect_equal(slice.tree_offspring.node(tree, '4', 'A'), slice.tree_offspring.node(tree, '4', 'B'))
})
message('.', appendLF=FALSE)

#Testing slice.tree_DELTRAN
#example
tree<-read.tree(text="(((((A:1,B:1):2,C:3):1,D:1):1,E:5):1,F:3);")
tree$node.label<-as.character(seq(1:5))
slice_tree<-suppressMessages(tree_slice<-timeSliceTree(tree, 3, drop.extinct=TRUE, plot=FALSE))
test<-slice.tree_DELTRAN(tree, 'A', tree_slice)

#test
test_that("Testing slice.tree_DELTRAN()", {
    #class
    expect_is(test, 'character')
    #length
    expect_equal(length(test), 1)
    #result (node 3)
    expect_equal(test, '3')
})
message('.', appendLF=FALSE)

#Testing ACCTRAN
#example
tree<-read.tree(text="(((((A:1,B:1):2,C:3):1,D:1):1,E:5):1,F:3);")
tree$node.label<-as.character(seq(1:5))
slice_tree<-suppressMessages(tree_slice<-timeSliceTree(tree, 3, drop.extinct=TRUE, plot=FALSE))
test<-slice.tree_ACCTRAN(tree, 'A', tree_slice)

#test
test_that("Testing slice.tree_ACCTRAN()", {
    #class
    expect_is(test, 'character')
    #length
    expect_equal(length(test), 1)
    #result (node 4)
    expect_equal(test, '4')
})
message('.', appendLF=FALSE)

#Testing slice.tree
tree<-read.tree(text="(((((A:1,B:1):2,C:3):1,D:1):1,E:5):1,F:3);")
tree$node.label<-as.character(seq(1:5))

#test
test_that("Testing slice.tree()", {
    #class
    expect_is(slice.tree(tree, 0, 'ACCTRAN'), 'phylo')
    #number of tips at time 0 (4)
    expect_equal(Ntip(slice.tree(tree, 0, 'ACCTRAN')), 4)
    #tips at time 0 ABCE
    expect_equal(sort(slice.tree(tree, 0, 'ACCTRAN')$tip.label), LETTERS[c(1:3,5)])
    expect_equal(sort(slice.tree(tree, 0, 'DELTRAN')$tip.label), LETTERS[c(1:3,5)])
    #to old
    expect_error(slice.tree(tree, 7, 'ACCTRAN'))
    #to few taxa
    expect_error(slice.tree(tree, 5.5, 'DELTRAN'))
})
message('.', appendLF=FALSE)