#SANITYZING FUNCTIONS


#Checking the class of an object and returning an error message if != class
check.class<-function(object, class, msg, errorif=FALSE) {
    #Set msg if missing
    if(missing(msg)) {
        msg<-paste(" must be ", class, ".", sep="")
    }
    #check if object is class.
    if(length(class) == 1) {
        if(errorif==FALSE) {
            if(class(object) != class) {
                stop(as.character(substitute(object)), msg , call.=FALSE)
            }
        } else {
            if(class(object) == class) {
                stop(as.character(substitute(object)), msg , call.=FALSE)
            }        
        }
    } else {
    #check if object is class in a cascade (class[1] else class[2] else class[3], etc..)
    #returns error only if object is not of any class
        for (i in 1:length(class)) {
            if(class(object) == class[i]) {
                class.test<-class[i]
            }        
        }
        if(exists(as.character(quote(class.test)))) {
            return(class.test)
        } else {
            stop(as.character(substitute(object)), msg , call.=FALSE)
        }
    }
}


#Checking the class of an object and returning an error message if != class
check.length<-function(object, length, msg, errorif=FALSE) {
    if(errorif==FALSE) {
        if(length(object) != length) {
            stop(as.character(substitute(object)), msg , call.=FALSE)
        }
    } else {
        if(length(object) == length) {
            stop(as.character(substitute(object)), msg , call.=FALSE)
        }        
    }
}


#Cleaning a tree so that the species match with the ones in a table
clean.tree<-function(tree, table, verbose=FALSE) {
    missing.species<-comparative.data(tree, data.frame("species"=row.names(table), "dummy"=rnorm(nrow(table)), "dumb"=rnorm(nrow(table))), "species")$dropped
    if(length(missing.species$tips) != 0) {
        tree.tmp<-drop.tip(tree, missing.species$tips)
        if (verbose==TRUE) {
            cat("Dropped tips:\n")
            cat(missing.species$tips, sep=", ")
        }
        tree<-tree.tmp
    }

    return(tree)
}

#Cleaning a table so that the species match with the ones in the tree
clean.table<-function(table, tree, verbose=FALSE) {
    missing.species<-comparative.data(tree, data.frame("species"=row.names(table), "dummy"=rnorm(nrow(table)), "dumb"=rnorm(nrow(table))), "species")$dropped
    if(length(missing.species$unmatched.rows) != 0) {
        table.tmp<-table[-c(match(missing.species$unmatched.rows, rownames(table))),]
        if (verbose==TRUE) {
            cat("Dropped rows:\n")
            cat(missing.species$unmatched.rows, sep=", ")
        }
        table<-table.tmp
    }
    return(table)
}

#Transforming a tree to binary with no 0 branch length.

bin.tree<-function(tree){
    if(!is.binary.tree(tree)) {
        tree<-multi2di(tree)
        warning('tree is now binary.' , call.=FALSE)
    }
    #Null branch length?
    if(any(tree$edge.length == 0)){
        tree$edge.length[which(tree$edge.length == 0)]<-min(tree$edge.length[-which(tree$edge.length == 0)])*0.01
        warning('New branches length generated are set to 1% of the minimum branch length.' , call.=FALSE)
    }
    return(tree)
}

#Replacing a value to be NA
replace.na<-function(x, y="?") {
    x[which(x == y)] <- NA
    return(x)
}

MorphDistMatrix.verbose <- function(morph.matrix, transform.proportional.distances="arcsine_sqrt", verbose=TRUE) {

  #verbose
  if(verbose == TRUE) {
    message("This is a verbose version of Claddis::MorphDistMatrix function v0.1.\n", appendLF=FALSE)
  }

  # Check format of transform.proportional.distances:
  if(transform.proportional.distances != "none" && transform.proportional.distances != "sqrt" && transform.proportional.distances != "arcsine_sqrt") {

    # Give error if something other than three possible settings is given:
    stop("ERROR: transform.proportional.distances must be one of \"none\", \"sqrt\", or \"arcsine_sqrt\".")

  }

  # Isolate ordering element of morphology matrix:
  ordering <- morph.matrix$ordering

  # Isolate max values element of morphology matrix:
  max.vals <- morph.matrix$max.vals
  
  # Isolate min values element of morphology matrix:
  min.vals <- morph.matrix$min.vals
  
  # Isolate weighting element of morphology matrix:
  weights <- morph.matrix$weights
  
  # Isolate character-taxon matrix element of morphology matrix:
  morph.matrix <- morph.matrix$matrix

  # Create empty vectors to store S and W value for Wills 2001 equations (1 and 2):
  differences <- maximum.differences <- vector(mode="numeric")
    
  # Distance matrices for storing:
  comp.char.matrix <- gower.dist.matrix <- max.dist.matrix <- dist.matrix <- matrix(0, nrow=length(rownames(morph.matrix)), ncol=length(rownames(morph.matrix)))
    
  # Fill comparable characters diagonal:
  comp.fill<-function(X) {
    return(length(X) - length(grep(TRUE, is.na(X))))
  }
  diag(comp.char.matrix)<-apply(morph.matrix, 1, comp.fill)

  #for(i in 1:length(morph.matrix[, 1])) {
  #  comp.char.matrix[i,i] <- length(morph.matrix[i, ]) - length(grep(TRUE, is.na(morph.matrix[i, ])))
  #}

  # Set up empty matrix for storing data to calculate the Generalised Euclidean Distance of Wills (2001):
  GED.data <- matrix(nrow=0, ncol=ncol(morph.matrix))

  #verbose
  if(verbose == TRUE) {
    message("Calculating the distances:", appendLF=FALSE)
  }

  # Go through matrix rows:
  for(i in 1:(length(morph.matrix[, 1]) - 1)) {
        
    # Go through matrix columns:
    for(j in (i + 1):length(morph.matrix[, 1])) {
            
      # Get just the comparable characters (those coded for both taxa):
      compchar <- intersect(grep(TRUE, !is.na(morph.matrix[rownames(morph.matrix)[i], ])), grep(TRUE, !is.na(morph.matrix[rownames(morph.matrix)[j], ])))
            
      # Get comparable characters for ith taxon:
      firstrow <- morph.matrix[rownames(morph.matrix)[i], compchar]

      # Get comparable characters for jth taxon:
      secondrow <- morph.matrix[rownames(morph.matrix)[j], compchar]
            
      # Deal with polymorphic characters (if present):
      if(length(grep("&", unique(c(firstrow, secondrow)))) > 0) {
                
        # Find ampersands (polymorphisms):
        ampersand.elements <- sort(c(grep("&", firstrow), grep("&", secondrow)))
                
        # Go through each polymorphic character:
        for(k in 1:length(ampersand.elements)) {
                    
          # Find out if two codings overlap once all polymorphism resolutions are considered:
          intersection.value <- intersect(strsplit(firstrow[ampersand.elements[k]], "&")[[1]], strsplit(secondrow[ampersand.elements[k]], "&")[[1]])
                    
          # Case if polymorphic and non-polymorphic values overlap:
          if(length(intersection.value) > 0) {
                        
            # Set ith value as zero (no difference):
            firstrow[ampersand.elements[k]] <- 0

            # Set jth value as zero (no difference)
            secondrow[ampersand.elements[k]] <- 0

          }
                    
          # Case if polymorphic and non-polymorphic values do not overlap:
          if(length(intersection.value) == 0) {
                        
            # Case if character is unordered (max difference is 1):
            if(ordering[compchar[ampersand.elements[k]]] == "unord") {
                            
              # Set ith value as zero:
              firstrow[ampersand.elements[k]] <- 0

              # Set jth value as 1 (making the ij difference equal to one):
              secondrow[ampersand.elements[k]] <- 1

            }
                        
            # Case if character is ordered (max difference is > 1):
            if(ordering[compchar[ampersand.elements[k]]] == "ord") {
                            
              # Get first row value(s):
              firstrowvals <- as.numeric(strsplit(firstrow[ampersand.elements[k]], "&")[[1]])
                            
              # Get second row value(s):
              secondrowvals <- as.numeric(strsplit(secondrow[ampersand.elements[k]], "&")[[1]])
                            
              # Make mini distance matrix:
              poly.dist.mat <- matrix(0, nrow=length(firstrowvals), ncol=length(secondrowvals))
                            
              # Go through each comparison:
              for(l in 1:length(firstrowvals)) {
                                
                # Record absolute difference:
                for(m in 1:length(secondrowvals)) poly.dist.mat[l, m] <- sqrt((firstrowvals[l] - secondrowvals[m]) ^ 2)

              }
                            
              # Set first value as zero:
              firstrow[ampersand.elements[k]] <- 0
                            
              # Set second value as minimum possible difference:
              secondrow[ampersand.elements[k]] <- min(poly.dist.mat)

            }

          }

        }

      }
      #End dealing with polymorphic characters (if present)

            
      # Get the absolute difference between the two rows:
      raw.diffs <- diffs <- abs(as.numeric(firstrow) - as.numeric(secondrow))
            
      # If there are differences greater than 1 for unordered characters then rescore as 1:
      if(length(grep(TRUE, diffs > 1)) > 0) diffs[grep(TRUE, diffs > 1)[grep(TRUE, ordering[compchar[grep(TRUE, diffs > 1)]] == "unord")]] <- 1

      # Find the incomparable characters:
      incompchar <- setdiff(1:ncol(morph.matrix), compchar)

      # Store data for GED with NAs for missing distances:
      GED.data <- rbind(GED.data, rbind(c(diffs, rep(NA, length(incompchar))), c(weights[compchar], weights[incompchar])))

      # Get weighted differences:
      diffs <- as.numeric(weights[compchar]) * diffs

      # Get raw Euclidean distance:
      raw.dist <- dist(rbind(diffs, rep(0, length(diffs))), method="euclidean")

      # Work out maximum difference (again, checked against ordering) using compchar characters only:
      raw.maxdiffs <- maxdiffs <- as.numeric(max.vals[compchar]) - as.numeric(min.vals[compchar])

      # Correct maximum possible differences for unordered characters:
      if(length(grep(TRUE, maxdiffs > 1)) > 0) maxdiffs[grep(TRUE, maxdiffs > 1)[grep(TRUE, ordering[compchar[grep(TRUE, maxdiffs > 1)]] == "unord")]] <- 1

      # Get vector of maximum differences (corrected for character weights):
      maxdiffs <- as.numeric(weights[compchar]) * maxdiffs

      # Store raw distance:
      dist.matrix[i, j] <- dist.matrix[j, i] <- raw.dist

      # Store Gower distance:
      gower.dist.matrix[i, j] <- gower.dist.matrix[j, i] <- sum(diffs) / sum(weights[compchar])

      # Store maximum-rescaled distance:
      max.dist.matrix[i, j] <- max.dist.matrix[j, i] <- sum(diffs) / sum(maxdiffs)

      # Store N comparable characters:
      comp.char.matrix[i, j] <- comp.char.matrix[j, i] <- length(compchar)

      # Add to maximum differences (S_ijk * W_ijk in equation 1 of Wills 2001):
      differences <- c(differences, diffs)

      # Add to maximum differences (S_ijk_max * W_ijk in equation 1 of Wills 2001):
      maximum.differences <- c(maximum.differences, maxdiffs)
      
      #verbose
      if(verbose == TRUE) {
        message(".", appendLF=FALSE)
      }

    }

  }
  #End of the loop for calculating the distance
  #verbose
  if(verbose == TRUE) {
    message("Done.\n", appendLF=FALSE)
  }

  # Calculated weighted mean univariate distance for calculating GED (equation 2 in Wills 2001):
  S_ijk_bar <- sum(differences) / sum(maximum.differences)

  # Replace missing distances with S_ijk_bar (i.e., results of equation 2 in Wills 2001 into equation 1 of Wills 2001):
  GED.data[is.na(GED.data)] <- S_ijk_bar

  # Isolate the distances:
  S_ijk <- GED.data[grep(TRUE, (1:nrow(GED.data) %% 2) == 1), ]

  # Isolate the weights:
  W_ijk <- GED.data[grep(TRUE, (1:nrow(GED.data) %% 2) == 0), ]

  # Calculate the GED (equation 1 of Wills 2001) for each pairwise comparison (ij):
  GED_ij <- sqrt(apply(W_ijk * (S_ijk ^ 2), 1, sum))

  # Create empty GED distance matrix:
  GED.dist.matrix <- matrix(0, nrow=nrow(morph.matrix), ncol=nrow(morph.matrix))

  # Set initial value for counter:
  counter <- 1

  #verbose
  if(verbose == TRUE) {
    message("Storing the distances:", appendLF=FALSE)
  }

  # Go through matrix rows:
  for(i in 1:(length(morph.matrix[, 1]) - 1)) {

    # Go through matrix columns:
    for(j in (i + 1):length(morph.matrix[, 1])) {

      # Store distance:
      GED.dist.matrix[i, j] <- GED.dist.matrix[j, i] <- GED_ij[counter]

      # Update counter:
      counter <- counter + 1

      #verbose
      if(verbose == TRUE) {
        message(".", appendLF=FALSE)
      }
    }

  }

  #verbose
  if(verbose == TRUE) {
    message("Done.\n", appendLF=FALSE)
  }
    
  # Set diagonals as zero:
  diag(gower.dist.matrix) <- diag(max.dist.matrix) <- 0

  # Add row and column names (taxa) to distance matrices:
  rownames(comp.char.matrix) <- colnames(comp.char.matrix) <- rownames(GED.dist.matrix) <- colnames(GED.dist.matrix) <- rownames(gower.dist.matrix) <- colnames(gower.dist.matrix) <- rownames(max.dist.matrix) <- colnames(max.dist.matrix) <- rownames(dist.matrix) <- colnames(dist.matrix) <- rownames(morph.matrix)

  # If transformation option is not "none":
  if(transform.proportional.distances != "none") {

    #verbose
    if(verbose == TRUE) {
      message("Transforming the proportional distances:...", appendLF=FALSE)
    }

    # If transformation option is "sqrt":
    if(transform.proportional.distances == "sqrt") {

      # Replace NaN with NA for Gower distances and take square root:
      gower.dist.matrix <- matrix(as.numeric(gsub(NaN, NA, sqrt(gower.dist.matrix))), nrow=nrow(gower.dist.matrix), dimnames=list(rownames(gower.dist.matrix), rownames(gower.dist.matrix)))

      # Replace NaN with NA for Max distances and take square root:
      max.dist.matrix <- matrix(as.numeric(gsub(NaN, NA, sqrt(max.dist.matrix))), nrow=nrow(max.dist.matrix), dimnames=list(rownames(max.dist.matrix), rownames(max.dist.matrix)))

    # If transformation option is "arcsine_sqrt":
    } else {

      # Establish correction factor to ensure Gower data is proportional:
      gower.correction <- max(c(max(sort(gower.dist.matrix)), 1))

      # Ensure all Gower values are on 0 to 1 scale then take arcsine of sqrt to get values that better approximate a normal distribution:
      gower.dist.matrix <- matrix(as.numeric(gsub(NaN, NA, asin(sqrt(gower.dist.matrix / gower.correction)))), nrow=nrow(gower.dist.matrix), dimnames=list(rownames(gower.dist.matrix), rownames(gower.dist.matrix)))

      # Take arcsine square root of all MOD dist values:
      max.dist.matrix <- matrix(as.numeric(gsub(NaN, NA, asin(sqrt(max.dist.matrix)))), nrow=nrow(max.dist.matrix), dimnames=list(rownames(max.dist.matrix), rownames(max.dist.matrix)))

    }
    
    #verbose
    if(verbose == TRUE) {
      message("Done.\n", appendLF=FALSE)
    }

  # If transformation option is "none":
  } else {

    # Replace NaN with NA for Gower distances:
    gower.dist.matrix <- matrix(as.numeric(gsub(NaN, NA, gower.dist.matrix)), nrow=nrow(gower.dist.matrix), dimnames=list(rownames(gower.dist.matrix), rownames(gower.dist.matrix)))

    # Replace NaN with NA for Max distances:
    max.dist.matrix <- matrix(as.numeric(gsub(NaN, NA, max.dist.matrix)), nrow=nrow(max.dist.matrix), dimnames=list(rownames(max.dist.matrix), rownames(max.dist.matrix)))

  }

  # Compile results as a list:
  result <- list(dist.matrix, GED.dist.matrix, gower.dist.matrix, max.dist.matrix, comp.char.matrix)
    
  # Add names to results list:
  names(result) <- c("raw.dist.matrix", "GED.dist.matrix", "gower.dist.matrix", "max.dist.matrix", "comp.char.matrix")
    
  # Output result:
  return(result)

}

#FUNCTIONS FOR DISPARITY

#Performs bootstrap and eventual rarefaction
Bootstrap.rarefaction<-function(data, bootstraps, rarefaction) {

    #Set rarefaction (or not)
    if(rarefaction == TRUE) {
        rarefaction_max<-seq(1:nrow(data))
    } else {
        rarefaction_max<-nrow(data)
    }
    #Rarefaction
    result<-NULL
    BSresult<-NULL
    for(rare in rarefaction_max){
        #Bootstraps
        for(BS in 1:bootstraps){ #bootstraps -> bootstraps
            #Bootstrap
            output<-as.matrix(data[sample(1:nrow(data),rare,TRUE),])
            result[BS] <- list(output)
        }
        #Rarefaction + BS results
        BSresult[rare]<-list(result)
    }

    #Remove first element if rarefaction
    if(rarefaction == TRUE) {
        BSresult[[1]]=NULL
    } else {
        #Removing the n-1 first elements
        BSresult<-BSresult[-c(1:(rarefaction_max-1))]
    }
    return(BSresult)
}

#Range Calculations
range.calc<-function(list_table) {
    #Empty matrix (for output)
    output<-matrix(nrow=length(list_table), ncol=ncol(list_table[[1]]))
    #Looping through columns and rows
    for(row in 1:length(list_table)) { #Rows are bootstraps
        for(column in 1:ncol(list_table[[1]])) { #Columns are axes
            output[row,column]<-max(list_table[[row]][,column])-min(list_table[[row]][,column])
        }
    }
    return(output)
}

#Variance calculation
variance.calc<-function(list_table) {
    #Empty matrix (for output)
    output<-matrix(nrow=length(list_table), ncol=ncol(list_table[[1]]))
    #Looping through columns and rows
    for(row in 1:length(list_table)) { #Rows are bootstraps
        for(column in 1:ncol(list_table[[1]])) { #Columns are axes
            output[row,column]<-var(list_table[[row]][,column])
        }
    }
    return(output)
}

#Set-up for the NthRoot function in order to scale your products correctly.
nth.root<-function(x, n){
    x^(1/n)
}

#Apply loop for calculating the product
prod.apply<-function(X) {
    output<-nth.root(apply(X,1,prod), n=ncol(X))
    return(output)
}

#Apply loop for calculating the sum
sum.apply<-function(X) {
    output<-apply(X,1,sum)
    return(output)
}

#No apply (does nothin)
no.apply<-function(X) {
    return(X)
}

#Apply loop for calculating the centroid
centroid.apply<-function(X) {
    #FUNCTION FROM SIVE, ADD DOI
    #Centroid (mean score of each PC axis)
    centroid<-apply(X, 2, mean)
    #Euclidean distances to the centroid
    cent.dist<-NULL
    for (j in 1:nrow(X)){
        cent.dist[j] <- dist(rbind(X[j,], centroid), method="euclidean")
    }
    return(cent.dist)
}

#Lapply loop for calculating the centroid
centroid.calc<-function(X) {
    Y<-lapply(X, centroid.apply)
    return(matrix(nrow=10, data=unlist(Y), byrow=TRUE))
}

#Converts one or more CI into a quantile probabilities
CI.converter<-function(CI) {
    sort(c(50-CI/2, 50+CI/2)/100)
}

#Calculate product for the central tendency and the CIs for variance or range
Disparity.measure.table<-function(type_function, distribution_variable, central_tendency, CI, save.all) {

    #Products/Sum of distribution_variable (correct by the NthRoot)
    Disparity_measure<-lapply(distribution_variable, type_function)
    #Confidence intervals for the products/sum of distribution_variable
    Disparity_measure_quantile<-lapply(Disparity_measure, quantile, probs=CI.converter(CI))
    #Calculate the central tendency
    Disparity_measure_central<-lapply(Disparity_measure, central_tendency)

    #Transform the results into a data.frame
    #Add just the first column (central tendency)
    Disparity_measure_table<-data.frame("Central"=unlist(Disparity_measure_central))
    #Add the CIs
    Disparity_measure_table<-cbind(Disparity_measure_table, matrix(ncol=length(CI)*2, data=unlist(Disparity_measure_quantile), byrow=TRUE))
    #Add the CIs names
    colnames(Disparity_measure_table)[-1]<-paste(CI.converter(CI)*100, "%", sep="")

    #return only the quantile table
    if(save.all == FALSE) {
        return(Disparity_measure_table)
    } else {
        output<-list("quantiles"=Disparity_measure_table, "values"=Disparity_measure[[1]])
        return(output)
    }

}

##########################
#Disparity functions
##########################
#Calculate the disparity as the distance from centroid
#This function is based on DisparityCalc() from Smith et al. 2014 - Evolution (http://dx.doi.org/10.1111/evo.12435) http://datadryad.org/resource/doi:10.5061/dryad.d380g 
#v0.2.2
##########################
#SYNTAX :
#<distance> the distance matrix
#<method> the method for calculating the disparity. Can be any of the following: "centroid", "sum.range", "product.range", "sum.variance", "product.variance"
#<CI> confidence intervals (default=c(50,95))
#<bootstraps> the number of boostrap replicates (default=1000)
#<central_tendency> any function for calculating the central tendency
#<verbose> whether to be verbose or not
#<rm.last.axis> whether to remove the last axis from the pco data. Can be a threshold value.
#<save.all> whether to save all the disparity values (TRUE) or just the quantiles (FALSE (default)).
##########################
#----
#guillert(at)tcd.ie 16/03/2015
##########################

disparity<-function(data, method=c("centroid", "sum.range", "product.range", "sum.variance", "product.variance"), CI=c(50, 95), bootstraps=1000, central_tendency=median, rarefaction=FALSE, verbose=FALSE, rm.last.axis=FALSE, save.all=FALSE) {

    #SANITIZING
    #distance
    check.class(data, "matrix", " must be a distance matrix.")

    #Test if applicable (> 2 rows)
    if(nrow(data) < 2) {
        stop("Disparity can not be calculated because less than two taxa are present in the data!")
    } 

    #method
    check.class(method, "character", " must be 'centroid', 'sum.range', 'product.range', 'sum.variance' or/and 'product.variance'.")
    methods_list<-c("centroid", "sum.range", "product.range", "sum.variance", "product.variance")
    if(all(is.na(match(method, methods_list)))) {
        stop("method must be 'centroid', 'sum.range', 'product.range', 'sum.variance' or/and 'product.variance'.")
    }

    #Bootstrap
    check.class(bootstraps, "numeric", " must be a single (entire) numerical value.")
    check.length(bootstraps, 1, " must be a single (entire) numerical value.")
    #Make sure the bootstrap is a whole number
    bootstraps<-round(abs(bootstraps))

    #CI
    #only relevant if bootstrap != 0)
    if(bootstraps != 0) {
        check.class(CI, "numeric", " must be any value between 1 and 100.")
        #remove warnings
        options(warn=-1)
        if(any(CI) < 1) {
            stop("CI must be any value between 1 and 100.")
        }
        if(any(CI) > 100) {
            stop("CI must be any value between 1 and 100.")
        }
        options(warn=0)
    }
    
    #Central tendency
    check.class(central_tendency, "function", " must be either a function (e.g. 'mean' or 'median'.")

    #rarefaction
    check.class(rarefaction, "logical", " must be logical.")

    #verbose
    check.class(verbose, "logical", " must be logical.")

    #rm.last.axis
    if(class(rm.last.axis) == "logical") {
        if(rm.last.axis == FALSE) {
            rm.axis<-FALSE
        } else {
            rm.axis<-TRUE
            last.axis<-0.95
        }
    } else {
        check.class(rm.last.axis, "numeric", " must be logical or a probability threshold value.")
        check.length(rm.last.axis, 1, " must be logical or a probability threshold value.", errorif=FALSE)
        if(rm.last.axis < 0) {
            stop("rm.last.axis must be logical or a probability threshold value.")
        } else {
            if(rm.last.axis > 1) {
                stop("rm.last.axis must be logical or a probability threshold value.")
            } else {
                rm.axis<-TRUE
                last.axis<-rm.last.axis
            }
        }
    }

    #verbose
    check.class(save.all, "logical", " must be logical.")

    #CALCULATING THE DISPARITY

    #Removing the last pco axis
    if(rm.axis==TRUE) {
        #calculate the cumulative variance per axis
        scree_data<-cumsum(apply(data, 2, var) / sum(apply(data, 2, var)))
        #extract the axis  below the threshold value
        axis_selected<-length(which(scree_data < last.axis))
        #remove the extra axis
        data<-data[,c(1:axis_selected)]
        #warning
        message(paste("The", length(scree_data)-axis_selected, "last axis have been removed from the pco data."))
    }

    #Bootstraping the matrix
    #verbose
    if(verbose==TRUE) {
        message("Bootstraping...", appendLF=FALSE)
    }
    BSresult<-Bootstrap.rarefaction(data, bootstraps, rarefaction)
    if(verbose==TRUE) {
        message("Done.", appendLF=TRUE)
    }

    #CENTROID
    #Distance form centroid
    if(any(method == 'centroid')) {
        #Calculate the distance from centroid for the rarefaction and the bootstrapped matrices
        if(verbose==TRUE) {
            message("Calculating distance from centroid...", appendLF=FALSE)
        }
        centroids<-lapply(BSresult, centroid.calc)
        #Distance to centroid
        Centroid_dist_table<-Disparity.measure.table(type_function=no.apply, centroids, central_tendency, CI, save.all)
        #Results type
        if(save.all == FALSE) {
            #Renaming the column
            colnames(Centroid_dist_table)[1]<-"Cent.dist"
        } else {
            #Isolating the data parts
            Centroid_values<-Centroid_dist_table[[2]]
            Centroid_dist_table<-Centroid_dist_table[[1]]
            colnames(Centroid_dist_table)[1]<-"Cent.dist"
        }
        if(verbose==TRUE) {
            message("Done.", appendLF=TRUE)
        }
    }

    #RANGES
    if(any(grep("range", method))) {
        #Calculate the range for the rarefaction and the bootstrapped matrices
        if(verbose==TRUE) {
            message("Calculating ranges...", appendLF=FALSE)
        }
        ranges<-lapply(BSresult, range.calc)

        #Sum of ranges
        if(any(method == 'sum.range')) {
            #Sum of range
            Sum_range_table<-Disparity.measure.table(type_function=sum.apply, ranges, central_tendency, CI, save.all)

            #Results type
            if(save.all == FALSE) {
                #Renaming the column
                colnames(Sum_range_table)[1]<-"Sum.range"
            } else {
                #Isolating the data parts
                Sum_range_values<-Sum_range_table[[2]]
                Sum_range_table<-Sum_range_table[[1]]
                colnames(Sum_range_table)[1]<-"Sum.range"
            }

        }

        #Product of ranges
        if(any(method == 'product.range')) {
            #Product of range
            Product_range_table<-Disparity.measure.table(type_function=prod.apply, ranges, central_tendency, CI, save.all)

            #Results type
            if(save.all == FALSE) {
                #Renaming the column
                colnames(Product_range_table)[1]<-"Prod.range"
            } else {
                #Isolating the data parts
                Prod_range_values<-Product_range_table[[2]]
                Product_range_table<-Product_range_table[[1]]
                colnames(Product_range_table)[1]<-"Prod.range"
            }

        }
        if(verbose==TRUE) {
            message("Done.", appendLF=TRUE)
        }
    }

    #VARIANCE
    if(any(grep("variance", method))) {
        #Calculate the variance for the rarefaction and the bootstrapped matrices
        if(verbose==TRUE) {
            message("Calculating variance...", appendLF=FALSE)
        }
        variances<-lapply(BSresult, variance.calc)

        #Sum of variance
        if(any(method == 'sum.variance')) {
            #Sum of variance
            Sum_variance_table<-Disparity.measure.table(type_function=sum.apply, variances, central_tendency, CI, save.all)

            #Results type
            if(save.all == FALSE) {
                #Renaming the column
                colnames(Sum_variance_table)[1]<-"Sum.var"
            } else {
                #Isolating the data parts
                Sum_variance_values<-Sum_variance_table[[2]]
                Sum_variance_table<-Sum_variance_table[[1]]
                colnames(Sum_variance_table)[1]<-"Sum.var"
            }
  
        }

        #Product of variance
        if(any(method == 'product.variance')) {
            #Product of variance
            Product_variance_table<-Disparity.measure.table(type_function=prod.apply, variances, central_tendency, CI, save.all)

            #Results type
            if(save.all == FALSE) {
                #Renaming the column
                colnames(Product_variance_table)[1]<-"Prod.var"
            } else {
                #Isolating the data parts
                Prod_variance_values<-Product_variance_table[[2]]
                Product_variance_table<-Product_variance_table[[1]]
                colnames(Product_variance_table)[1]<-"Prod.var"
            }
        
        }
        if(verbose==TRUE) {
            message("Done.", appendLF=TRUE)
        }
    }

    #OUTPUT
    #Empty output table
    if(rarefaction==FALSE) {
        output<-matrix(nrow=1, data=rep(NA, 1))
        colnames(output)[1]<-"rarefaction"
    } else {
        output<-matrix(nrow=(nrow(data)-1), data=seq(from=2, to=nrow(data)))
        colnames(output)[1]<-"rarefaction"
    }

    #Distance form centroid
    if(any(method == 'centroid')) {
        #Combine the results
        output<-cbind(output, Centroid_dist_table)
    }
    #Sum of ranges
    if(any(method == 'sum.range')) {
        #Combine the results
        output<-cbind(output, Sum_range_table)
    }
    #Product of ranges
    if(any(method == 'product.range')) {
        #Combine the results
        output<-cbind(output, Product_range_table)
    }
    #Sum of variance
    if(any(method == 'sum.variance')) {
        #Combine the results
        output<-cbind(output, Sum_variance_table)  
    }
    #Product of variance
    if(any(method == 'product.variance')) {
        #Combine the results
        output<-cbind(output, Product_variance_table)   
    }

    if(save.all == FALSE) {
        #Quantiles only
        return(output)
    } else {
        #Quantiles and values
        output<-list("table"=output)
        #Add the values of each metric
        #centroid
        if(any(method == 'centroid')) {
            output[[length(output)+1]]<-Centroid_values
            names(output)[length(output)]<-"centroid"
        }
        #Sum of sum.range
        if(any(method == 'sum.range')) {
            output[[length(output)+1]]<-Sum_range_values
            names(output)[length(output)]<-"sum.range"
        }
        #Product of ranges
        if(any(method == 'product.range')) {
            output[[length(output)+1]]<-Prod_range_values
            names(output)[length(output)]<-"product.range"
        }
        #Sum of variance
        if(any(method == 'sum.variance')) {
            output[[length(output)+1]]<-Sum_variance_values
            names(output)[length(output)]<-"sum.variance"
        }
        #Product of variance
        if(any(method == 'product.variance')) {
            output[[length(output)+1]]<-Prod_variance_values
            names(output)[length(output)]<-"product.variance"
        }

        return(output)
    }

#End
}

##########################
#int.pco
##########################
#Select the number of taxa per intervals
#v0.2.2
##########################
#SYNTAX :
#<pco_data> the pco data to split in intervals.
#<tree> a 'phylo' object. The tree must be dated.
#<intervals> a series of intervals breaks limits.
#<FAD_LAD> a data.frame containing the first and last apparition datums. If none is provided, or if taxa are missing, taxa are assumed to have the same FAD and LAD.
#<include.nodes> logical, whether to include nodes or not in the intervals. default = FALSE. If TRUE, the nodes must be the same name in the pco_data and in the tree.
#<diversity> logical, whether to count the number of taxa in each interval.
##########################
#Update: fixed FAD_LAD to me more plastic: if input FAD_LAD contains extra taxa, they are now being discarded from the analysis.
#----
#guillert(at)tcd.ie 19/03/2014
##########################

int.pco<-function(pco_data, tree, intervals, FAD_LAD, include.nodes=FALSE, diversity=FALSE) {

    #SANITIZING
    #pco
    check.class(pco_data, 'matrix', ' must be a pco scores matrix.')

    #tree
    check.class(tree, 'phylo', ' must be a phylo object.')
    #the tree must be dated
    if(length(tree$root.time)==0){
        stop("Tree must be a dated tree with $root.time.")
    }

    #intervals
    if(class(intervals) != 'numeric') {
        if(class(intervals) != 'integer') {
            stop("intervals must be numeric.")
        }
    }
    #length must be greater than one
    if(length(intervals) < 2) {
        stop("At least two breaks should be specified for the intervals.")
    }

    #FAD_LAD
    if(missing(FAD_LAD)) {
        #Create the FAD_LAD table
        FAD_LAD<-data.frame("FAD"=tree.age(tree)[1:Ntip(tree),1], "LAD"=tree.age(tree)[1:Ntip(tree),1], row.names=tree.age(tree)[1:Ntip(tree),2])
        message("No FAD_LAD table has been provided so every tip is assumed to interval single points in time.")
    } else {
        #Check if the FAD_LAD contains all taxa
        if(any(tree$tip.label %in% rownames(FAD_LAD) == FALSE)) {
            message("Some tips have FAD/LAD and are assumed to interval single points in time.")
            #If not generate the FAD_LAD for the missing taxa
            missing_FADLAD<-which(is.na(match(tree$tip.label, rownames(FAD_LAD))))
            add_FAD_LAD<-data.frame(tree.age(tree)[missing_FADLAD,1], tree.age(tree)[missing_FADLAD,1], row.names=tree.age(tree)[missing_FADLAD,2])
            colnames(add_FAD_LAD)<-colnames(FAD_LAD)
            FAD_LAD<-rbind(FAD_LAD, add_FAD_LAD)
        }
        #Remove FAD_LAD taxa not present in the tree
        if(nrow(FAD_LAD) != Ntip(tree)) {
            FAD_LAD<-FAD_LAD[-c(which(is.na(match(rownames(FAD_LAD), tree$tip.label)))),]
        }

    }

    #include.nodes
    check.class(include.nodes, 'logical', " must be logical.")
    #Check if nodes are present in the pco_data object and in the tree
    if(include.nodes == TRUE) {
        #Check if node labels are present
        if(length(tree$node.label) == 0) {
            stop("Provided tree has no nodes labels.")
        } else {
            for (node in 1:length(tree$node.label)) {
                if(length(grep(tree$node.label[node], rownames(pco_data))) == 0) {
                    stop(paste("node", tree$node.label[node], "not found."))
                }
            }
        }

        #Check if the pco_data contains more nodes
        if(nrow(pco_data) > (Ntip(tree)+Nnode(tree))) {
            message("Some rows in pco_data are not present in the tree!")
        }

    } else {
        #Check if the pco_data contains more nodes
        if(nrow(pco_data) > (Ntip(tree)+Nnode(tree))) {
            message("Some rows in pco_data are not present in the tree!")
        }
    }

    #diversity
    check.class(diversity, "logical", " must be logical.")


    #BINING THE PCO
    #ages of tips/nodes + FAD/LAD
    ages_tree_FAD<-tree.age(tree)
    ages_tree_LAD<-tree.age(tree)
    #Change the age if FAD or LAD are higher/lower than the age of the tip
    for(tip in 1:nrow(FAD_LAD)) {
        #Replace age of the tip if FAD is higher
        if(FAD_LAD[tip,1] > ages_tree_FAD$ages[which(ages_tree_FAD$edges == rownames(FAD_LAD)[tip])]) {
            ages_tree_FAD$ages[which(ages_tree_FAD$edges == rownames(FAD_LAD)[tip])]<-FAD_LAD[tip,1]
        }
        #Replace age of the tip if LAD is lower
        if(FAD_LAD[tip,2] < ages_tree_LAD$ages[which(ages_tree_LAD$edges == rownames(FAD_LAD)[tip])]) {
            ages_tree_LAD$ages[which(ages_tree_LAD$edges == rownames(FAD_LAD)[tip])]<-FAD_LAD[tip,2]
        }
    }

    #Empty list element per interval
    int_elements<-NULL
    int_elements<-list()

    #Attribute each taxa/node to it's interval
    for (interval in 1:(length(intervals)-1)) {
        #Select the elements of one interval
        int_elements[[interval]]<-ages_tree_FAD$edges[which(ages_tree_FAD$ages >= intervals[interval+1] & ages_tree_LAD$ages <= intervals[interval])]
    }
    
    #Remove the nodes (if necessary)
    if(include.nodes==FALSE) {
        for (interval in 1:length(int_elements)) {
        #Remove nomatch with tree$tip.label
            int_elements[[interval]]<-int_elements[[interval]][match(tree$tip.label, int_elements[[interval]])[-which(is.na(match(tree$tip.label, int_elements[[interval]])))]]
        }
    }

    #Making the pco interval list
    pco_intervals<-NULL
    pco_intervals<-list()

    for (interval in 1:length(int_elements)) {
        #Matching list
        matching<-match(as.character(int_elements[[interval]]),as.character(rownames(pco_data)))
        #If only one taxa is matching, make sure it's not a vector
        if(length(matching) == 1) {
            pco_intervals[[interval]]<-matrix(data=pco_data[matching,], nrow=1)
            rownames(pco_intervals[[interval]])<-rownames(pco_data)[matching]
        } else {
            pco_intervals[[interval]]<-pco_data[matching,]
        }
    }

    #Naming the intervals
    name_list<-NULL
    for(interval in 1:length(int_elements)) {
        name_list[interval]<-paste(intervals[interval], intervals[interval+1], sep="-")
    }
    
    #If interval is empty, send warning and delete the interval
    #list of empty intervals (empty)
    empty_intervals<-NULL
    for (interval in 1:length(pco_intervals)) {
        if(length(pco_intervals[[interval]]) == 0) {
            #Remove the interval
            empty_intervals[interval]<-interval
            #Select the empty interval
            empty_interval<-paste(intervals[interval], intervals[interval+1], sep="-")
            message("The following interval is empty: ", empty_interval, ".")
        }
    }

    #If any empty intervals
    if(!is.null(empty_intervals)) {
        #NA removal from empty_intervals vector (if any)
        if(any(is.na(empty_intervals))) {
            empty_intervals<-empty_intervals[-which(is.na(empty_intervals))]
        }
        #Removing the empty intervals
        pco_intervals<-pco_intervals[c(-empty_intervals)]
        #Removing the empty intervals names
        name_list<-name_list[-empty_intervals]
    }
    
    names(pco_intervals)<-name_list

    #Diversity
    if(diversity == TRUE) {
        #count the elements per intervals
        if(!is.null(empty_intervals)) {
            diversity_counts<-unlist(lapply(int_elements[-empty_intervals], length))
        } else {
            diversity_counts<-unlist(lapply(int_elements, length))
        }
        #add the interval names
        names(diversity_counts)<-name_list

        #Output
        output<-list("pco_intervals"=pco_intervals, "diversity"=diversity_counts)
        return(output)
    
    } else {
        return(pco_intervals)
    }
}

##########################
#Tree ages
##########################
#Extract the node and the tips ages of a chronogram
#v1.1
#Update: Syntax and typos
##########################
#SYNTAX :
#<tree> a 'phylo' object
#<scale> the scale of the tree (i.e. the age of the root)
#<type> either 'past' if the tree must be a dated tree (units = time to past; tips=0), or 'present' if the tree is an absolute dated tree (units = time to present; root=0)
##########################
#----
#guillert(at)tcd.ie 30/06/2014
#Modified from [R-sig-phylo] nodes and taxa depth II - 21/06/2011 - Paolo Piras - ppiras(at)uniroma3.it
##########################


tree.age<-function(tree, scale, type='past'){

#SANITYZING

    #tree
    check.class(tree, 'phylo', ' must be a phylo object.')

    #scale
    if(missing(scale)) {
      #Using the tree height as scale if scale is missing
        scale=max(dist.nodes(tree)[, Ntip(tree)+1])
    }
    check.class(scale, 'numeric', ' must be a numerical value.')
    check.length(scale, '1', ' must a a single value.')

    #type
    check.class(type, 'character', ' must be \'past\' or \'present\'.')
    if(type !='past') {
        if(type !='present') {
            stop('type must be \'past\' or \'present\'.')
        }
    }

#CALCULATE THE EDGES AGE

    if(scale == 0) {
        ages.table<-tree.age_table(tree)
    } else {
        ages.table<-tree.age_scale(tree.age_table(tree), scale)
    }

    #Type
    if(type == 'past'){
        tree.height<-max(ages.table$ages)
        ages.table$ages<-round(abs(ages.table$ages-tree.height), digit=3)
    } else {
        ages.table$ages<-round(ages.table$ages, digit=7)
    }

    #Output
    #ages.table<-round(ages.table[1,], digit=3)
    return(ages.table)

    #Example
    example=FALSE
    if(example == TRUE){
        #Examples
        ##Generate a birth-death tree with fossil and living species
        library(diversitree)
        tree<-tree.bd(c(1,0.3), max.taxa=20, include.extinct=TRUE)

        ##Calculate the edges age by setting the root a 65 Mya
        Treeage(tree, scale=65)

        ##Ploting the distribution of the node ages
        hist(ages[-c(1:length(tree$tip.label)),1], xlab="Time", main="Divergence frequency per time")

        ##Calculate when the fossil went extinct
        ages.fossil<-Treeage(tree, scale=65, type='past')
        tree.fossil<-tree
        tree.fossil$tip.label[grep("ex",ages.fossil[,2])]<-paste(tree.fossil$tip.label[grep("ex",ages.fossil[,2])],round(ages.fossil[grep("ex",ages.fossil[,2]),1], digit=2), "Mya", sep=" ")
        plot(tree.fossil)

        ##Ploting the node age from root
        ages.nodes<-Treeage(tree, scale=65, type='present') #change 'present' into 'past' to plot a classical chronogram
        tree.nodes<-tree
        tree.nodes$node.label<-round(ages.nodes[-c(1:length(tree.nodes$tip.label)),1], digit=2)
        plot(tree.nodes)
        nodelabels(tree.nodes$node.label, cex=0.6)
    }
}


#FUNCTIONS FOR tree.age


#Extract the ages table from a tree
#Calculating the tips and the edges age
tree.age_table<-function(tree){
    ages<-dist.nodes(tree)[length(tree$tip.label)+1,]
    tip.names<-tree$tip.label
    if(is.null(tree$node.label)) {
        nod.names<-c((length(tree$tip.label)+1):length(dist.nodes(tree)[,1]))
    } else {
        nod.names<-tree$node.label
    }
    edges<-c(tip.names, nod.names)
    ages.table<-data.frame(ages=ages,edges=edges)
    return(ages.table)
}


#Scaling the ages from tree.age_table if scale != 1
tree.age_scale<-function(ages.table, scale){
    ages.table$ages<-ages.table$ages/max(ages.table$ages)
    ages.table$ages<-ages.table$ages*scale
    return(ages.table)
}

##########################
#time.disparity
##########################
#Calculates the disparity for interval pco.data and output a interval.disparity table object
#v0.2.3
##########################
#SYNTAX :
#<time_pco> time intervals or slices from a pco
#<...> disparity arguments (see ?disparity for information)
##########################
#----
#guillert(at)tcd.ie 16/03/2014
##########################

time.disparity<-function(time_pco, method=c("centroid", "sum.range", "product.range", "sum.variance", "product.variance"), CI=c(50, 95), bootstraps=1000, central_tendency=median, rarefaction=FALSE, verbose=FALSE, rm.last.axis=FALSE, save.all=FALSE) {
    #SANITIZING
    #time_pco
    check.class(time_pco, "list", " must be a list of time sections of pco data.")
    if(length(names(time_pco))!=length(time_pco)) {
        stop("time_pco data must have time sections names.")
    }

    #rarefaction
    if(rarefaction == TRUE) {
        message("Rarefaction is calculated and slows down the disparity calculation.\nUse Rarefaction=FALSE to speed up the calculations.")
    }

    #Managing bins with only one data point
    while(any(unlist(lapply(time_pco, nrow)) < 2)) {
        wrong_intervals<-which(unlist(lapply(time_pco, nrow)) < 2)
        #Moving the first wrong interval to the next interval in time (unless the wrong interval is the last one)
        if(wrong_intervals[1] != length(time_pco)) {
            host_interval<-wrong_intervals[1]+1
            message("Intervals ", names(time_pco)[wrong_intervals[1]], " and ", names(time_pco)[host_interval], " are combined due to insufficient data.")
        } else {
            #Moving the wrong interval in the preceding one
            host_interval<-wrong_intervals[1]-1
            message("Intervals ", names(time_pco)[host_interval], " and ", names(time_pco)[wrong_intervals[1]], " are combined due to insufficient data.")
        }
        #Creating the new interval
        new_interval<-rbind(time_pco[[wrong_intervals[1]]], time_pco[[host_interval]])
        #Creating the new time_pco data
        new_time_pco<-time_pco ; names(new_time_pco)<-names(time_pco)
        #replacing the wrong interval
        new_time_pco[[host_interval]]<-new_interval
        #renaming the interval
        if(wrong_intervals[1] != length(time_pco)) {
            names(new_time_pco)[host_interval]<-paste(strsplit(names(new_time_pco)[wrong_intervals[1]], split="-")[[1]][1],strsplit(names(new_time_pco)[host_interval], split="-")[[1]][2],sep="-")
        } else {
            names(new_time_pco)[host_interval]<-paste(strsplit(names(new_time_pco)[host_interval], split="-")[[1]][1],strsplit(names(new_time_pco)[wrong_intervals[1]], split="-")[[1]][2],sep="-")
        }
        #removing empty interval
        new_time_pco[[wrong_intervals[1]]]<-NULL
        time_pco<-new_time_pco
    }

    #CALCULATING THE DISPARITY FOR EACH BIN
    disparity_interval<-lapply(time_pco, disparity, method=method, CI=CI, bootstraps=bootstraps, central_tendency=central_tendency, rarefaction=rarefaction, verbose=verbose, rm.last.axis=rm.last.axis, save.all=save.all)

    #Return the table only
    if(save.all == FALSE) {
        #Sorting the data as a table
        if(rarefaction == FALSE) {
            #create the table's first row
            disparity_intervals_table<-disparity_interval[[1]]
            #Loop through the other elements of the table
            for(interval in 2:length(disparity_interval)) {
                disparity_intervals_table<-rbind(disparity_intervals_table, disparity_interval[[interval]])
            }
            #Renaming the rarefaction column interval
            colnames(disparity_intervals_table)[1]<-"time"
            #Saving the interval names
            disparity_intervals_table[,1]<-names(time_pco)

        } else {

            #If rarefaction has been calculated, only get the last element of each rarefaction table

            #create an interval row
            interval_row<-matrix(nrow=(nrow(disparity_interval[[1]])), data=rep(names(time_pco)[[1]]))
            #add the disparity results (with rarefaction)
            interval_tab<-cbind(interval_row, disparity_interval[[1]])
            #binding the interval table
            disparity_intervals_table<-rbind(interval_tab)

            #Loop through the other intervals
            for(interval in 2:length(time_pco)) {
                #create an interval row
                interval_row<-matrix(nrow=(nrow(disparity_interval[[interval]])), data=rep(names(time_pco)[[interval]]))
                #add the disparity results (with rarefaction)
                interval_tab<-cbind(interval_row, disparity_interval[[interval]])
                #binding the interval table
                disparity_intervals_table<-rbind(disparity_intervals_table, interval_tab)
            }

            #Renaming the rarefaction column interval
            colnames(disparity_intervals_table)[1]<-"time"
        }

        return(disparity_intervals_table)
    
    } else {

        #Sorting the data as a table
        if(rarefaction == FALSE) {
            #Creating the quantile table
            #create the table's first row
            disparity_intervals_table<-disparity_interval[[1]][[1]]
            #Loop through the other elements of the table
            for(interval in 2:length(disparity_interval)) {
                disparity_intervals_table<-rbind(disparity_intervals_table, disparity_interval[[interval]][[1]])
            }
            #Renaming the rarefaction column interval
            colnames(disparity_intervals_table)[1]<-"time"
            #Saving the interval names
            disparity_intervals_table[,1]<-names(time_pco)

        } else {

            #If rarefaction has been calculated, only get the last element of each rarefaction table

            #create an interval row
            interval_row<-matrix(nrow=(nrow(disparity_interval[[1]][[1]])), data=rep(names(time_pco)[[1]]))
            #add the disparity results (with rarefaction)
            interval_tab<-cbind(interval_row, disparity_interval[[1]][[1]])
            #binding the interval table
            disparity_intervals_table<-rbind(interval_tab)

            #Loop through the other intervals
            for(interval in 2:length(time_pco)) {
                #create an interval row
                interval_row<-matrix(nrow=(nrow(disparity_interval[[interval]][[1]])), data=rep(names(time_pco)[[interval]]))
                #add the disparity results (with rarefaction)
                interval_tab<-cbind(interval_row, disparity_interval[[interval]][[1]])
                #binding the interval table
                disparity_intervals_table<-rbind(disparity_intervals_table, interval_tab)
            }

            #Renaming the rarefaction column interval
            colnames(disparity_intervals_table)[1]<-"time"
        }

        #saving the results per time section
        disparity_intervals_values<-list()
        #Loop through all the elements of the list to extract the values
        for(interval in 1:length(disparity_interval)) {
            disparity_intervals_values[[interval]]<-disparity_interval[[interval]][[2]]
        }
        #Renaming the list elements
        names(disparity_intervals_values)<-names(time_pco)

        #output
        output<-list("quantiles"=disparity_intervals_table, "values"=disparity_intervals_values)
        return(output)
    }
}

##########################
#Plotting disparity results
##########################
#Plots the disparity results
#v0.2
##########################
#SYNTAX :
#<disparity> disparity data
#<measure> the name of the column containing the disparity measurement. If set to 'default' the measure will be the first measure (second column) of the table.
#<rarefaction> whether to plot the rarefaction results or not
#<diversity> optional. Must be a vector of the same length as disparity_data.
#<add> optional. Whether to add the a previous called graph
##########################
#----
#guillert(at)tcd.ie 19/03/2015
##########################

plot.disparity<-function(disparity_data, measure="default", rarefaction=FALSE, xlab="default", ylab="default", col="default", las=2, diversity, add=FALSE, ylim, ...){
    #SANITIZING
    #Disparity
    check.class(disparity_data, 'data.frame', " must be a disparity data.frame")
    if(length(disparity_data) < 4) {
        stop("Disparity data.frame must have in the following order:\na 'rarefaction' column, a 'measurement' column and at least two 'Confidence Interval' columns.\nUse the disparity() function to generate the proper formatted data.frame.")
    }

    #Measure
    check.class(measure, 'character', " must be 'default' or one of the names of the columns in the disparity data.frame.")
    check.length(measure, 1, " must be 'default' or one of the names of the columns in the disparity data.frame.", errorif=FALSE)
    #Get the right column number
    if(measure == 'default') {
        measure_col<-2
    } else {
        measure_col<-grep(measure, colnames(disparity_data))
        if(length(measure_col) != 1) {
            stop("measure column not found in disparity_data.\nUse the disparity() function to generate the proper formatted data.frame.")
        }
    }

    #Get the Confidence intervals columns
    measure_col_tmp<-measure_col
    while(length(grep("%", colnames(disparity_data)[(measure_col_tmp+1)]))==1) {
        measure_col_tmp<-measure_col_tmp+1
    }
    #Get the column position for the different CI_values
    CI_length<-(measure_col_tmp-measure_col)
    CI_min<-measure_col+1
    CI_max<-measure_col+CI_length
    #Get the eventual intermediate CIs
    #number of CI values
    n_CI<-CI_length/2
    #extracting all the CI column number pairs in a table
    CI_pairs<-matrix(nrow=n_CI, ncol=2)
    for(n in 1:n_CI) {
        CI_pairs[n,1]<-CI_min+n-1
        CI_pairs[n,2]<-CI_max-n+1
    }
    #Setting the line types for the CIs
    lty_list<-c(44,33,22,21,12)

    #Rarefaction
    check.class(rarefaction, 'logical')
    #Check if rarefaction data is available
    options(warn=-1)
    if(rarefaction == TRUE & is.na(disparity_data[,1])) {
        stop("No rarefaction data available.\nUse the disparity() function to generate the proper formatted data.frame with the option 'rarefaction=TRUE'.")
    }
    #If no rarefaction, check if disparity data is > 1
    if(rarefaction == FALSE & nrow(disparity_data) < 2) {
        warning("Only one disparity point is available.")
    }
    options(warn=0)

    #xlab
    check.class(xlab, "character", " must be a character string.")
    check.length(xlab, 1, " must be a character string.", errorif=FALSE)
    #ylab
    check.class(ylab, "character", " must be a character string.")
    check.length(ylab, 1, " must be a character string.", errorif=FALSE)
    #col
    check.class(col, "character", " must be a character string.")

    #diversity
    if(missing(diversity)) {
        plot.diversity<-FALSE
    } else {
        plot.diversity<-TRUE
        #check.class(diversity, "integer", " must be a numeric vector of the same number of rows as disparity_data.")
        check.length(diversity, nrow(disparity_data), " must be a numeric vector of the same number of rows as disparity_data.")
    }

    #add
    check.class(add, "logical")

    #PLOTTING THE DISPARITY RESULTS
    if(add == FALSE) {
        if(rarefaction == TRUE) {
            #ylim
            if(missing(ylim)) {
                ylim=c(min(disparity_data[,CI_min]),max(disparity_data[,CI_max]))
            }
            #Plotting the rarefaction curve
            plot(disparity_data[,1], disparity_data[,measure_col], type='l', ylim=ylim , ...)
            #Add the CIs
            for (n in 1:(CI_length/2)) {
                #Add both lines
                lines(disparity_data[,1], disparity_data[,CI_pairs[n,1]], type='l', lty=lty_list[n+1])
                lines(disparity_data[,1], disparity_data[,CI_pairs[n,2]], type='l', lty=lty_list[n+1])
            }

        } else {
            #Plotting the disparity curve
            if(nrow(disparity_data) == 1) {
                #ylim
                if(missing(ylim)) {
                    ylim=c(min(disparity_data[,CI_min]),max(disparity_data[,CI_max]))
                }
                #If only one data point is available, do box plot style
                plot(1,1, xlab='', ylab='', ylim=ylim, type='n', xaxt='n')
                points(1,disparity_data[,measure_col], pch=19)
                #line types for this one
                lty_list2<-c(44,1,1,1,1,1)
                for (n in 1:(CI_length/2)) {
                    #Add CIs lines
                    lines(c(1,1), c(disparity_data[,CI_pairs[n,1]], disparity_data[,CI_pairs[n,2]]), lwd=1+(n-1)*3, lty=lty_list2[n])
                }
            } else {
                #Setting plot options (defaults)
                #xlab
                if(xlab=="default") {
                    xlab="bins"
                }
                #ylab
                if(ylab=="default") {
                    ylab=names(disparity_data)[measure_col]
                }
                #colors
                if(col=="default") {
                    #line color
                    line_color<-"black"
                    #polygon_colors
                    polygon_colors<-c("lightgrey", "grey")
                } else {
                    #line color
                    line_color<-col[1]
                    #polygon_colors
                    polygon_colors<-col
                }

                #ylim
                if(missing(ylim)){
                    ylim=c(min(disparity_data[,CI_min]),max(disparity_data[,CI_max]))
                }

                #Plotting the curve
                plot(seq(from=1, to=nrow(disparity_data)), disparity_data[,measure_col], type='l', 
                    ylim=ylim ,col="white", ylab=ylab, xlab=xlab, xaxt='n' , ...)
                    #ylim=ylim ,col="white", ylab=ylab, xlab=xlab, xaxt='n'); warning("debug")
                if(class(disparity_data[,1]) == "character") {
                    axis(side = 1, 1:nrow(disparity_data), disparity_data[,1], las=las)
                } else {
                    axis(side = 1, 1:nrow(disparity_data))
                }
                #Add the polygons
                for (n in 1:(CI_length/2)) {
                    polygon(c(seq(from=1, to=nrow(disparity_data)), seq(from=nrow(disparity_data), to=1)),
                        c(disparity_data[,CI_pairs[n,1]], rev(disparity_data[,CI_pairs[n,2]])),
                        col=polygon_colors[n], density=100-(100/(CI_length/2)/2.5*n))
                }
                #ylim
                if(missing(ylim)){
                    ylim=c(min(disparity_data[,CI_min]),max(disparity_data[,CI_max]))
                }
                #Add the central tendency line
                lines(seq(from=1, to=nrow(disparity_data)), disparity_data[,measure_col], type='l', ylim=ylim, col=line_color)

                #Add the diversity (optional)
                if(plot.diversity==TRUE) {
                    par(new=TRUE)
                    plot(diversity, type="l", lty=2, xaxt="n",yaxt="n",xlab="",ylab="")
                    axis(4)
                }
            }
        }
    } else {

        #Adding to the existing plot
        stop("add=TRUE option in development.")
    }
}