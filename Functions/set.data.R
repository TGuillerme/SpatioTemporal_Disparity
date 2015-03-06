##########################
#Setting up the data for disparity analysis
##########################
#Wrapper for setting up the data for disparity analysis
#v0.1
##########################
#SYNTAX :
#<matrix_path> the path to the nexus matrix
#<tree_path> the path to the nexus tree
#<ace> whether to estimate ancestral states
#<dist> whether to calculate the distance matrix
#<save> whether to create a read.data file and save the results
#<verbose> whether to be verbose or not
##########################
#----
#guillert(at)tcd.ie 06/03/2015
##########################

set.data<-function(matrix_path, tree_path, ace=TRUE, dist=TRUE, save=TRUE, verbose=TRUE) {
    #SANITIZING
    #matrix_path
    check.class(matrix_path, 'character', ' must be the path to the matrix data.')
    check.length(matrix_path, 1, ' must be the path to the matrix data.', errorif=FALSE)

    #tree_path
    check.class(tree_path, 'character', ' must be the path to the tree data.')
    check.length(tree_path, 1, ' must be the path to the matrix data.', errorif=FALSE)

    #ace
    check.class(ace, 'logical', ' must be logical.')

    #dist
    check.class(dist, 'logical', ' must be logical.')

    #save
    check.class(save, 'logical', ' must be logical.')

    #verbose
    check.class(verbose, 'logical', ' must be logical.')

    

}