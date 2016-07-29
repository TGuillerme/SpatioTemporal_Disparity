##########################
#Reading landmark data (TPS style)
##########################
#Reading landmark coordinates output by TPS
#v.0.1
##########################
#SYNTAX :
#<data> the path to a csv file.
#
##########################
#----
#guillert(at)tcd.ie - 14/07/2016
##########################
#Requirements:
#-R 3


#This amounted to 47 total points collected: 7 ‘internal’ landmarks (figure 2a), 20 points around the talonid, and 20 points around the trigonid

#######
# Prepare the data
#######

# Reading in raw data
raw_coordinates <- read.csv("../Data/molar_LM_coords_GrossnickleNewham2016.csv", header = F, stringsAsFactors = F)

# Removing the taxa with no coordinates
raw_coordinates <- raw_coordinates[-which(raw_coordinates[,1] == 0),]
# Removing the associated ID and LMs
LMs <- which(raw_coordinates[,1] == "LM=47")
# Which of these have an NA in second position
removed_taxa <- which(is.na(raw_coordinates[LMs+1,2]))
raw_coordinates <- raw_coordinates[-c(LMs[removed_taxa], LMs[removed_taxa]+1),]

# Getting the species names
species_names_raws <- grep("ID=", raw_coordinates[,1])
species_names_list <- as.list(as.character(raw_coordinates[species_names_raws, 1]))
species_names <- unlist(lapply(species_names_list, strsplit, split = "ID="))
species_names <- species_names[-which(species_names == "")]

# Getting the number of specimens
number_specimen <- length(species_names)

# Getting landmarks function
get.landmarks <- function(raw_data, number_of_specimens, shift_size, landmarks) {
    list_out <- list()
    for(spec in 1:number_specimen) {
        #specimen shifter
        shifter <- (shift_size)*(spec-1)
        #storing the internal landmarks in the list
        list_out[[spec]] <- data.matrix(raw_data[(landmarks[1]+shifter):(landmarks[2]+shifter),])
    }
    return(list_out)
}

# Getting the internal landmarks
internal_landmarks <- get.landmarks(raw_coordinates, number_specimen, 49, c(2,8))
talonid_slides <- get.landmarks(raw_coordinates, number_specimen, 49, c(9,28))
trigonid_slides <- get.landmarks(raw_coordinates, number_specimen, 49, c(29,48))

# Adding the specimens names
names(internal_landmarks) <- names(talonid_slides) <- names(trigonid_slides) <- species_names

# Reorient all the landmarks and curves
# The first landmark (paraconid) must have the lowest x value
# The sixth landmark (hypoconulid) must have the highest x value


# for(which_outline in 1:length(internal_landmarks)) {

#     one_landmark_set <- internal_landmarks[[which_outline]]
#     #Set the detection range
#     detection_x <- range(one_landmark_set[,1])[2] - range(one_landmark_set[,1])[1]
#     detection_x <- detection_x - 0.1*detection_x #with 10% tolerance
#     detection_y <- range(one_landmark_set[,2])[2] - range(one_landmark_set[,2])[1]
#     detection_y <- detection_y - 0.1*detection_y #with 10% tolerance

#     if(one_landmark_set[6,1] < one_landmark_set[1,1] && abs(one_landmark_set[1,1]-one_landmark_set[6,1]) > detection_x) {
#         internal_landmarks[[which_outline]] <- rotate.coordinates(internal_landmarks[[which_outline]], mirror = "horizontal")
#         talonid_slides[[which_outline]] <- rotate.coordinates(talonid_slides[[which_outline]], mirror = "horizontal")
#         trigonid_slides[[which_outline]] <- rotate.coordinates(trigonid_slides[[which_outline]], mirror = "horizontal")
#     } else {
#         if(one_landmark_set[6,2] < one_landmark_set[1,2] && abs(one_landmark_set[6,2]-one_landmark_set[1,2]) > detection_y) {
#             internal_landmarks[[which_outline]] <- rotate.coordinates(internal_landmarks[[which_outline]], rotation = 270)
#             talonid_slides[[which_outline]] <- rotate.coordinates(talonid_slides[[which_outline]], rotation = 270)
#             trigonid_slides[[which_outline]] <- rotate.coordinates(trigonid_slides[[which_outline]], rotation = 270)
#         } else {
#             if(one_landmark_set[6,2] > one_landmark_set[1,2] && abs(one_landmark_set[6,2]-one_landmark_set[1,2]) > detection_y) {
#                 internal_landmarks[[which_outline]] <- rotate.coordinates(internal_landmarks[[which_outline]], rotation = 90)
#                 talonid_slides[[which_outline]] <- rotate.coordinates(talonid_slides[[which_outline]], rotation = 90)
#                 trigonid_slides[[which_outline]] <- rotate.coordinates(trigonid_slides[[which_outline]], rotation = 90)
#             }
#         }
#     }
# }

rotate.coordinates <- function(coordinates, rotation, mirror) {

    # warning("proportions not preserved!")
    rotated_coords <- coordinates

    if(!missing(rotation) && rotation != 0) {
        # rotation_matrix <- matrix(c(cos(rotation), -sin(rotation), sin(rotation), cos(rotation)), 2, 2)
        # rotated_coords <- coordinates %*% rotation_matrix

        rotated_coords[,1] = coordinates[,1] * cos(rotation) - coordinates[,2] * sin(rotation)
        rotated_coords[,2] = coordinates[,2] * cos(rotation) + coordinates[,1] * sin(rotation)
        # x' = x cos f - y sin f
        # y' = y cos f + x sin f

    }
    if(!missing(mirror)) {
        if(mirror == "vertical") {
            rotated_coords[,1] <- rotated_coords[,1]*-1
        }
        if(mirror == "horizontal") {
            rotated_coords[,2] <- rotated_coords[,2]*-1
        }
    }

    return(rotated_coords)
}


# Plotting the teeth!
plot.teeth <- function(which_outline, talonid_slides, trigonid_slides, internal_landmarks, legend = TRUE, axis = TRUE) {
    op <- par(bty = "n")
    # Plots the talonid
    if(axis) {
        plot(rbind(talonid_slides[[which_outline]], trigonid_slides[[which_outline]], internal_landmarks[[which_outline]]), col = "white", xlab = "", ylab = "")
    } else {
        plot(rbind(talonid_slides[[which_outline]], trigonid_slides[[which_outline]], internal_landmarks[[which_outline]]), col = "white", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
    }
    lines(talonid_slides[[which_outline]], col = "blue")
    # Plots the trigonid
    lines(trigonid_slides[[which_outline]], col = "red")
    # Plots the cusp (landmarks)
    points(internal_landmarks[[which_outline]], pch = 19)
    # Legend
    if(legend) {
        limits <- apply(rbind(talonid_slides[[which_outline]], trigonid_slides[[which_outline]], internal_landmarks[[which_outline]]), 2, range)
        legend(limits[1,1], limits[2,2], legend = c("Talonid", "Trigonid", "Cusps"), lty = 1, pch = 19, col = c("blue", "red", "black"), cex = 0.75, bty = "n")
    }
    par(op)
}

par(mfrow = c(5,7))
for(i in 1:35) {
    plot.teeth(i, talonid_slides, trigonid_slides, internal_landmarks, legend = FALSE, axis = FALSE)
}


# plot(rotate.coordinates(bla, mirror = "vertical"), type = "l")


# Combining the talonid and trigonid slides
combined_slides <- mapply(rbind, talonid_slides, trigonid_slides, SIMPLIFY = FALSE)
# Combining all the data
slides_n_landmarks <- mapply(rbind, internal_landmarks, combined_slides, SIMPLIFY = FALSE)

# Set data as an array
list.to.array <- function(list) {
    return(array(unlist(list), dim = c(nrow(list[[1]]), ncol(list[[1]]), length(list))))
}

# Setting the different arrays
internal_landmarks_array <- list.to.array(internal_landmarks)
combined_slides_array <- list.to.array(combined_slides)
slides_n_landmarks_array <- list.to.array(slides_n_landmarks)

# Curves data for procrustes analysis

# Procrustes analysis
library(geomorph)
proc_internal_landmarks <- gpagen(internal_landmarks_array)

ordination <- prcomp(proc_internal_landmarks$data)$x

#######
# 2D Geomorph
#######


#######
# 2D Fourier
#######


