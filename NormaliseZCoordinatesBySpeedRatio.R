# ------------------------------------------------
# Normalise Z coordinates by Speed Ratio
# ------------------------------------------------
# Code provided as is and can be used or modified freely. 
#
# The code is from the following paper:
#
# Demsar U and Long JA, 2019, 
# Potential Path Volume - a geometric estimator for space use in 3D
# 
# ------------------------------------------------
# Author: Urska Demsar
# University of St Andrews
# St Andrews, Scotland, UK
# http://udemsar.com
# ------------------------------------------------
# This file re-scales Z coordinates by ration of v_z_maz/v_xy_max
# ------------------------------------------------
# Supporting R files:
# DistanceTwo3Dpoints.R
# WriteCSVtableFromRectangular3DArrays.R
# ------------------------------------------------
# Input: data with a list of trajectories
# This must be a csv file including the following variables:
# - TrajID - ID of each trajectory
# - X - x coordinate of each trajectory point in some projected
# coordinate system (in metres)
# - Y - y coordinate of each trajectory point in some projected
# coordinate system (in metres)
# - Z - z coordinate of each trajectory point in some projected
# coordinate system (in metres)
# - T - time  
# ------------------------------------------------
# Output: a csv file where Z coordinates are re-scaled
# ------------------------------------------------


# ------------------------------------------------
# Parameters to be set/changed by the user

# Set working directory
setwd('.')

# Name of input and output files
# Input
trajectories.file <- 'One4Dsegment.csv'

# Output
output.file <- 'One4Dsegment_Rescaled.csv'

# ------------------------------------------------
# General setup & read data

# Read data
dfAll <- read.csv(trajectories.file,stringsAsFactors=FALSE)
head(dfAll)

# ------------------------------------------------
# Outer Loop across trajectories

# Find unique trajectories and count how many there are
listOfTrajIDs <- unique(dfAll$TrajID)
noOfTraj <- length(listOfTrajIDs)

# Initialise segment-based geometric parameters
dfAll$dz <- 0     # distance on segment in z direction  
dfAll$dxy <- 0    # distance on segment on xy plane
dfAll$dt <- 0     # time on segment: dt_i = T_i - T_(i-1) 
dfAll$vz <- 0     # z speed on segment: v_i = d_i / dt_i
dfAll$vxy <- 0    # xy speed on segment: v_i = d_i / dt_i

# Start outer loop: for each trajectory
print('Loop 1 - parameter calculation')

for (shk in 1:noOfTraj) {
  #for (shk in 1:1) {
  
  # How long is this trajectory?
  trajIndices <- which(dfAll$TrajID == listOfTrajIDs[shk])
  trajLength <- length(trajIndices)
  
  # ----------------------------------------
  # Inner Loop 1: for each segment calculate
  # physical parameters of movement along each trajectory
  
  # For all trajectory points except the first one, 
  # calculate physical parameters along the segment that ends in point i
  # i.e. segment P_i-1 to P_i
 
  for (i in 2:trajLength) {
    
    print(c(shk,i))

    # distance between the two 3D points that define the segment: d_i = d(P_i, P_(i-1))
    dfAll$dz[trajIndices[i]] <- abs(dfAll$Z[trajIndices[i-1]]-dfAll$Z[trajIndices[i]])
    dfAll$dxy[trajIndices[i]] <- sqrt((dfAll$X[trajIndices[i-1]]-dfAll$X[trajIndices[i]])^2+
                                  (dfAll$Y[trajIndices[i-1]]-dfAll$Y[trajIndices[i]])^2)
    
    # time on segment: dt_i = T_i - T_(i-1)
    dfAll$dt[trajIndices[i]] <- dfAll$T[trajIndices[i]]-dfAll$T[trajIndices[i-1]]

    # speeds on segment: v_i = d_i / dt_i
    dfAll$vz[trajIndices[i]] <- dfAll$dz[trajIndices[i]]/dfAll$dt[trajIndices[i]]
    dfAll$vxy[trajIndices[i]] <- dfAll$dxy[trajIndices[i]]/dfAll$dt[trajIndices[i]]
        
  } # Closing inner loop 1: for (i in 2:trajLength)
  
  # ----------------------------------------
  # Calculate max velocity on this trajectory (Long and Nelson 2012)
  
  # Find the largest velocities
  vz_m <- max(dfAll$vz[trajIndices])
  vxy_m <- max(dfAll$vxy[trajIndices])

  # Calculate the scaling factor: if vz = 0, then factor is 1
  # if vxy = 0, then factor is vz,
  # otherwise take the ratio
  if (vz_m == 0) { v_sf <-1
  } else if (vxy_m == 0) { v_sf <- v_z 
  } else { v_sf <- vz_m/vxy_m } 
  
  # Rescale Z coordinates
  dfAll$Z[trajIndices] <- dfAll$Z[trajIndices]/v_sf
  
} # Closing outer loop: for (shk in 1:noOfTraj)

# Export re-scaled coordinates

write.csv(dfAll[,1:5],output.file,row.names=FALSE)
