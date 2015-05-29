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

#quantiles
dis_pro_max_beck<-extract.disp(disparity_full_ran_beck$quantiles, rarefaction="max")

adonis(beck_65~beck_60, permutations=1000, method="euclidean")
adonis(beck_65~beck_55, permutations=1000, method="euclidean")
adonis(beck_65~beck_50, permutations=1000, method="euclidean")
adonis(beck_65~beck_45, permutations=1000, method="euclidean")
adonis(beck_65~beck_40, permutations=1000, method="euclidean")
adonis(beck_65~beck_35, permutations=1000, method="euclidean")
adonis(beck_65~null_ran_centroid_ran[[2]]$values$`65`, permutations=1000, method="euclidean")

#NPMANOVA of the PC axes  (e.g. Stayton 2005 and Ruta 2013)
 PC.man <- adonis(PC95axes~sp.fam$Family, data=sp.fam, permutations=999, method="euclidean")


beck_test<-list(as.vector(beck_60),as.vector(beck_55),as.vector(beck_50),as.vector(beck_45),as.vector(beck_40),as.vector(beck_35),as.vector(beck_30))

bla<-lapply(beck_test, bhatt.coeff, y=as.vector(beck_65))
bhatt.coeff(as.vector(beck_65), as.vector(null_ran_centroid_ran[[2]]$values$`65`))

#Just do a Tukey HSD?