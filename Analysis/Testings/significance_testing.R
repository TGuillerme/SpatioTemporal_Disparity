# Testing significance

#############################
# Function for loading raw scripts from GitHub 
#############################
# Just loading the whole bunch of functions. I you want to look at them in detail, I've added the GitHub links for each individual function as I call them through this script
library(devtools)
install_github("TGuillerme/SpatioTemporal_Disparity/Functions/disparity")
library(disparity)

#############################
# Data loading
#############################
# This is just exactly your workshop script up until the distance matrix calculations (without the "install.packages")

# Load the morphological matrix
nexus.data <- ReadMorphNexus("http://www.graemetlloyd.com/nexus/Cullen_etal_2013a.nex")

# Load the tree
tree.data <- read.tree("http://www.graemetlloyd.com/firstmpt/Cullen_etal_2013a.tre")

# Fix polytomies.
tree.data <- multi2di(tree.data)

# Read the age data
ages.data <- read.table("http://www.graemetlloyd.com/teaching/RE2014/Cullenages.txt", row.names=1, sep="\t", header=T)

# Building the tree
tree.data <- timePaleoPhy(tree.data, ages.data, type="mbl", vartime=2)

# Safe taxonomic reduction
safe.data <- SafeTaxonomicReduction(nexus.data)

#############################
# Distance matrix
#############################

# Distance matrix
dist.data <- MorphDistMatrix(nexus.data)

# For small datasets this is instantly but I'm playing around with reasonably big ones (100 taxa * 400 characters) and some can take a day or so to compute.
# To make sure the script is not frozen, I've just added a verbose function.
# I clumsily renamed it MorphDistMatrix.verbose to avoid any confusion.
#browseURL("https://github.com/TGuillerme/SpatioTemporal_Disparity/blob/master/Functions/disparity/R/MorphDistMatrix.verbose.R")

# The function is exactly the same but can take an optional additional argument: verbose (TRUE (default) or FALSE)
# With verbose = FALSE, the function is exactly the same
dist.data.verbose <- MorphDistMatrix.verbose(nexus.data, verbose=FALSE)

#The verbose option just adds some dots and some text (probably worth modifying but that's the general idea)
dist.data.verbose <- MorphDistMatrix.verbose(nexus.data, verbose=TRUE)

# Just for the sake of testing the differences
for (different_distances in 1:length(dist.data)) {
    print( all( (dist.data.verbose[[different_distances]] == dist.data[[different_distances]]), na.rm=TRUE) )
}

# Let's just continue through your workshop to get to the disparity part

# Trim inapplicable data
trimmed.max.data <-TrimMorphDistMatrix(dist.data$max.dist.matrix)

# Remove the trimmed taxa from the tree for later
tree.data<-drop.tip(tree.data, trimmed.max.data$removed.taxa)


# Do the proper Multidimensional Scaling (corrected for negative eigenvalues)
pco <- cmdscale(trimmed.max.data$dist.matrix, k=nrow(trimmed.max.data$dist.matrix) - 2, add=T, eig=TRUE)
pco.data <- pco$points #generates k-2 eigenvectors (=n cladistic space dimensions)
pco.eigen<- pco$eig #generates k eigenvalues (including two null eigenvalues)

bins<-rev(seq(from=0, to=100, by=20))
pco_in_bins<-int.pco(pco.data, tree.data, bins, FAD_LAD=ages.data, diversity=TRUE)


#############################
# Testing significance
############################# 

#Creating three different disparities

disp_med.dist.cent.BS<-time.disparity(pco_in_bins$pco_intervals, method="centroid", relative=FALSE, verbose=TRUE, bootstraps=10, rarefaction=TRUE, centroid.type="median", central_tendency=mean, save.all=TRUE)
disp_obs.dist.cent<-time.disparity(pco_in_bins$pco_intervals, method="centroid", relative=FALSE, verbose=TRUE, bootstraps=0,  rarefaction=FALSE, centroid.type="full", central_tendency=median, save.all=TRUE)
disp_obs.dist.cent.BS<-time.disparity(pco_in_bins$pco_intervals, method="centroid", relative=FALSE, verbose=TRUE, bootstraps=10, rarefaction=TRUE, centroid.type="full", central_tendency=median, save.all=TRUE)

#Extract disparity from the rarefaction analysis
disp_med.dist.cent.BS_max<-extract.disp(disp_med.dist.cent.BS$quantiles, rarefaction="max")

op<-par(mfrow=c(1,3))
plot.disparity(extract.disp(disp_med.dist.cent.BS$quantiles, rarefaction="max"), ylim=c(0.8, 1.8))
plot.disparity(disp_obs.dist.cent$quantiles, ylim=c(0.8, 1.8))
plot.disparity(extract.disp(disp_obs.dist.cent.BS$quantiles, rarefaction="max"), ylim=c(0.8, 1.8))
par(op)



#Before and after KT difference
KT<-65

preKT<-unlist(c(disp_med.dist.cent.BS$values$`100-80`, disp_med.dist.cent.BS$values$`80-60`))
posKT<-unlist(c(disp_med.dist.cent.BS$values$`60-40`, disp_med.dist.cent.BS$values$`40-20`, disp_med.dist.cent.BS$values$`20-0`))


wilcox.test(preKT, posKT)

data(dune)
data(dune.env)
adonis(dune ~ Management*A1, data=dune.env, permutations=99)


### Example of use with strata, for nested (e.g., block) designs.

dat <- expand.grid(rep=gl(2,1), NO3=factor(c(0,10)),field=gl(3,1) )
dat
Agropyron <- with(dat, as.numeric(field) + as.numeric(NO3)+2) +rnorm(12)/2
Schizachyrium <- with(dat, as.numeric(field) - as.numeric(NO3)+2) +rnorm(12)/2
total <- Agropyron + Schizachyrium
dotplot(total ~ NO3, dat, jitter.x=TRUE, groups=field,
        type=c('p','a'), xlab="NO3", auto.key=list(columns=3, lines=TRUE) )

Y <- data.frame(Agropyron, Schizachyrium)
mod <- metaMDS(Y)
plot(mod)
### Hulls show treatment
with(dat, ordihull(mod, group=NO3, show="0"))
with(dat, ordihull(mod, group=NO3, show="10", col=3))
### Spider shows fields
with(dat, ordispider(mod, group=field, lty=3, col="red"))

### Correct hypothesis test (with strata)
adonis(Y ~ NO3, data=dat, strata=dat$field, perm=999)

### Incorrect (no strata)
adonis(Y ~ NO3, data=dat, perm=999)



disparity_full_ran_beck

#values for each slice
beck_65<-disparity_full_ran_beck$values$`65`
beck_60<-disparity_full_ran_beck$values$`60`
beck_55<-disparity_full_ran_beck$values$`55`
beck_50<-disparity_full_ran_beck$values$`50`
beck_45<-disparity_full_ran_beck$values$`45`
beck_40<-disparity_full_ran_beck$values$`40`
beck_35<-disparity_full_ran_beck$values$`35`
beck_30<-disparity_full_ran_beck$values$`30`

dis_pro_max_beck

#quantiles
dis_pro_max_beck<-extract.disp(disparity_full_ran_beck$quantiles, rarefaction="max")

adonis(beck_65~beck_60, permutations=1000, method="euclidean")
adonis(beck_65~beck_55, permutations=1000, method="euclidean")
adonis(beck_65~beck_50, permutations=1000, method="euclidean")
adonis(beck_65~beck_45, permutations=1000, method="euclidean")
adonis(beck_65~beck_40, permutations=1000, method="euclidean")
adonis(beck_65~beck_35, permutations=1000, method="euclidean")
adonis(beck_65~null_ran_centroid_ran[[2]]$values$`65`, permutations=1000, method="euclidean")


adonis(beck_65~beck_60+beck_55+beck_50,beck_45+beck_40+beck_35+beck_30, permutations=1000, method="euclidean")
#NPMANOVA of the PC axes  (e.g. Stayton 2005 and Ruta 2013)
 PC.man <- adonis(PC95axes~sp.fam$Family, data=sp.fam, permutations=999, method="euclidean")

adonis(dis_pro_max_beck[22,3]~dis_pro_max_beck[22,3], permutations=1000)

beck_test<-list(as.vector(beck_60),as.vector(beck_55),as.vector(beck_50),as.vector(beck_45),as.vector(beck_40),as.vector(beck_35),as.vector(beck_30))

bla<-lapply(beck_test, bhatt.coeff, y=as.vector(beck_65))
bhatt.coeff(as.vector(beck_65), as.vector(null_ran_centroid_ran[[2]]$values$`65`))

#Just do a Tukey HSD?



#DisparityTTest from Smith et al 

#This function uses the t-distribution to calculate significant differences in disparity among time bins. This script is based on the t-test script of Anderson and Friedman (2012).
#We added additional disparity metrics to their script.
#Output is a symmetric matrix of p-values. Note that the output is not corrected for multiple comparisons so the R function p.adjust() may be required.
#Depending on the number of comparisons you perform, you are likely to also need to do some type of p-value correction (e.g. Bonferroni). This can be achieved using the “p.adjust” function in R. 
#For example, say you were looking at four adjacent time bins (a, b, c, and d) calculated using the “DisparityTTest” function and were interested in the adjusted p-value between them (a vs. b, b vs. c, c vs. d). Say we obtained p-values of 0.035, 0.0012, and 0.081.
#By typing the following line of code  you can correct for multiple comparisons.
#p.adjust(c(0.035, 0.0012, 0.081), method="bonferroni")
#[1] 0.1050 0.0036 0.2430
#Each value here has been adjusted by the Bonferroni correction method. In essence, this method  multiplies the uncorrected  p-value by the number of comparisons (in this case, 3 comparisons were  made). 0.035 was corrected to 0.1050, 0.0012 was corrected to 0.0036, and so on

DisparityTTest<-function(data, replicates=1000, upper.bound=0.975, lower.bound=0.025, method=c("Sum of Ranges","Sum of Variances","Product of Ranges","Product of Variances")){


#data = Matrix: A PCO matrix with taxa in the first column, time bins in the second column, and the PCO data in the subsequent columns
#replicates = 1000 Numeric: Number of iterations to calculate a t-distribution
#upper/lower bound = CI
#method = disparity methods

bin1<-pco_in_bins[[1]][[1]]
bin2<-pco_in_bins[[1]][[2]]
bin3<-pco_in_bins[[1]][[5]]

#Creating the "name" vector (col1)
bin1[,1]<-row.names(bin1)
#Creating the "age" vector (col2)
bin1[,2]<-sample(c(1:2), nrow(bin1), replace=TRUE)

data<-bin1


##Sorting out bins
groupsVector<<-as.character(as.matrix(data[ ,2])) #Selecting the bins


#create vector for every "group" indicating which taxa belong to each group
groupBins <- list()
allLevelsList <- strsplit(groupsVector,split="\\|")
allLevels <- unlist(allLevelsList)
allLevels <- allLevels[allLevels!="\\|"]
groupLevels <- unique(allLevels)

for(currentLevel in groupLevels){
  tempBin <- NULL
  for(taxon in 1:length(allLevelsList)){
    if(any(allLevelsList[[taxon]]==currentLevel))
      tempBin <- c(tempBin, taxon)
  }
  groupBins[[currentLevel]] <- tempBin
}
groupBins<<-groupBins #Selecting which taxa is in which group (bin)


bins <- groupLevels
number.bins <- length (groupLevels)
standard.error <- array (dim = number.bins)
sample.size <- array (dim = number.bins)

#specify arrays to hold dates and disparity values
complete.variance <- array (dim = c(replicates, number.bins))
mean.variance <- array (dim = number.bins)
t.statistics <- array (dim = c(number.bins, number.bins))
degrees.freedom <- array (dim = c (number.bins, number.bins))
significance.values <- array (dim = c(number.bins, number.bins))
number.axes <- ncol (data) - 2 #calculates the number of axes to be analyzed
CI <- array (dim = c(number.bins, 2))

#Set up for the product disparity metrics to scale by the nth root.
NthRoot<-function(x=data, n=number.axes){
    x^(1/n)
}

Methods<-sort(method)
IList<-list()



####Sum of Ranges Statistics####

#TG: this is just a normal t-test comparing means...

if (any(Methods == "Sum of Ranges")){
    for (i in 1: replicates){
        for (k in 1: number.bins){
            bin.data <- data[groupBins[[k]],]

            dimensions <- dim (bin.data)
            height <- dimensions [1]
            #Bootstrapping the data
            bin.data <- bin.data [sample (1:height, replace = TRUE), ]
            sample.size [k] <- height
            variance.holder <- 0
            
            for (m in 1:number.axes){
                variance.holder <- variance.holder + ((max(as.numeric(bin.data[, (m+2)])))-(min(as.numeric(bin.data[, (m+2)]))))
            }

            complete.variance [i, k] <- variance.holder #Different sum of ranges
        }
    }

    for (j in 1:number.bins){
        #Getting the CIs?
        standard.error [j] <- sd (complete.variance [, j])
        mean.variance [j] <- mean (complete.variance [,j])
        CI [j, 1] <- quantile (complete.variance [, j], probs = lower.bound)
        CI [j, 2] <- quantile (complete.variance [, j], probs = upper.bound)
    }


    for (x in 1: number.bins){ # go down rows
        for (q in 1:number.bins){ #goes across columns
            denominator <- (mean.variance[x]-mean.variance[q])
            numerator.one <- ((sample.size[x]-1)*(sample.size[x])*((standard.error[x])^2) + (sample.size[q]-1)*(sample.size[q])*((standard.error[q])^2))/(sample.size[x]+sample.size[q]+2) 
            numerator.two <- (sample.size[x] + sample.size [q])/(sample.size[x] * sample.size [q]) 
            t.statistics [x, q] <- denominator/((numerator.one*numerator.two)^.5)
            degrees.freedom [x, q] <- sample.size[x] + sample.size [q] - 2
            significance.values [x, q] <- 1-pt(t.statistics [x, q], df = degrees.freedom [x, q])
            
            #make test two-tailed
            if (significance.values [x,q] > 0.5)
            {
                significance.values [x,q] <- 2*(1-significance.values[x,q])
            }
            else if (significance.values [x,q] < 0.5)
            {
                significance.values [x,q] <- 2*(significance.values[x,q])   
            }
            else if (significance.values [x,q] == 0.5)
            {
                significance.values [x,q] <- 1  

            }   
        }
    }
 }


set.seed(1)
X1<-rnorm(100, 0)
X2<-rnorm(100, 1)
pt<-t.test(X1,X2)
plot(density(X1), col="red", xlim=c(min(c(X1,X2)), max(c(X1,X2))))
lines(density(X2), col="blue")