# ------------------------------------------------
# Potential Path Volume around a set of 3D trajectories
# ------------------------------------------------
# Code provided as is and can be used or modified freely. 
#
# The code is from the following paper:
#
# Demsar U and Long JA, 2018, 
# Potential Path Volume - a geometric estimator for space use in 3D
# 
# ------------------------------------------------
# Author: Urska Demsar
# University of St Andrews
# St Andrews, Scotland, UK
# http://udemsar.com
# ------------------------------------------------
# This is the main file that implements the PPV calculation
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
# The following also need to be manually specified: 
# - volume extent (min/max X, min/max Y, min/max T)
# - voxel size
# ------------------------------------------------
# Output: a csv file with PPV volume, i.e. with four columns:
# x, y, z, PPV value (1 if voxel is inside the volume, 0 otherwise)
# ------------------------------------------------
# We also provide two test files, one with a trajectory with
# one segment only (One4Dsegment.csv) and one with three 
# segments (Three4DSegments.csv).
# ------------------------------------------------


# ------------------------------------------------
# Parameters to be set/changed by the user

# Set working directory
setwd('.')

# Name of input and output files
# Input
trajectories.file <- 'One4Dsegment.csv'
# trajectories.file <- 'Three4Dsegments.csv' 

# Output
output.file <- 'One4Dsegment_PPV.csv'
# output.file <- 'Three4Dsegments_PPV.csv'

# Extent: coordinates of eight points that define the extent - for Voxler
# extent.file <- 'One4Dsegment_Extent.csv'
# extent.file <- 'Three4Dsegments_Extent.csv'
# extent.file <- '3DRandomWalk_40points_Extent.csv'
extent.file <- '4155724_Extent.csv'

# Output file with trajectories and all parameters
parameter.file <- 'One4Dsegment_Parameters.csv'
# parameter.file <- 'Three4Dsegments_Parameters.csv'

# Output file with no. of voxels, voxelSize and volume of PPV (voxelSize^3*noOfvoxels)
volume.file <- 'One4Dsegment_Values.csv'
#volume.file <- 'Three4Dsegments_Values.csv'

# ------------------------------------------------
# General setup & read data

# Read fuctions in additional files
# Function for saving volume as as csv coordinate file
source("WriteCSVtableFromRectangular3DArrays.R") 
# Functions for geometric and kernel calculations
source("DistanceTwo3Dpoints.R")

# Package for 3D plotting & colour schemes
library(plot3D)
#library(RColorBrewer)
#library(plotly)

# Read data
dfAll <- read.csv(trajectories.file,stringsAsFactors=FALSE)
head(dfAll)

# Set extent of the volume 
# Read extent from data 
minXcoord <- min(dfAll$X) # westernmost point 
maxXcoord <- max(dfAll$X) # easternmost point
minYcoord <- min(dfAll$Y) # southernmost point
maxYcoord <- max(dfAll$Y) # northernmost point
minZcoord <- min(dfAll$Z) # lowest point
maxZcoord <- max(dfAll$Z) # highest point

# ------------------------------------------------
# Outer Loop across trajectories

# Find unique trajectories and count how many there are
listOfTrajIDs <- unique(dfAll$TrajID)
noOfTraj <- length(listOfTrajIDs)

# Initialise segment-based geometric parameters
dfAll$d <- 0      # distance d_i = d(P_i, P_(i-1)), length of segment  
dfAll$dt <- 0     # time on segment: dt_i = T_i - T_(i-1) 
dfAll$v <- 0      # velocity on segment: v_i = d_i / dt_i
dfAll$a <- 0      # major axis of ellipsoid a_i
dfAll$b <- 0      # minor axes of ellipsoid b_i 
dfAll$xc <- 0     # central point on segment P_C, X coordinate 
dfAll$yc <- 0     # dentral point on segment P_C, Y coordinate 
dfAll$zc <- 0     # central point on segment P_C, Z coordinate 
dfAll$alpha <- 0  # first rotation angle, alpha_i 
dfAll$beta <- 0   # second rotation angle, beta_i 
dfAll$direction <- 0 # +/- for direction angle

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
    dfAll$d[trajIndices[i]] <- DistanceTwo3Dpoints(dfAll$X[trajIndices[i-1]],dfAll$Y[trajIndices[i-1]],dfAll$Z[trajIndices[i-1]],
                                                   dfAll$X[trajIndices[i]],dfAll$Y[trajIndices[i]],dfAll$Z[trajIndices[i]])
    
    # time on segment: dt_i = T_i - T_(i-1)
    dfAll$dt[trajIndices[i]] <- dfAll$T[trajIndices[i]]-dfAll$T[trajIndices[i-1]]

    # velocity on segment: v_i = d_i / dt_i
    dfAll$v[trajIndices[i]] <- dfAll$d[trajIndices[i]]/dfAll$dt[trajIndices[i]]
    
  } # Closing inner loop 1: for (i in 2:trajLength)
  
  # ----------------------------------------
  # Calculate max velocity on this trajectory (Long and Nelson 2012)
  
  # Find the largest velocity
  v_m <- max(dfAll$v[trajIndices])

  # Find next largest velocity to five decimal places of precision
  vv <- dfAll$v[trajIndices]
  #vv <- as.data.frame(dfAll$v[trajIndices])
  sort.v <- sort(vv, decreasing = TRUE) # sort in descending order
  where.next <- 1 
  while (round(sort.v[where.next],digits=5) == round(sort.v[1], digits=5)) {where.next <- where.next + 1}
  v_m1 <- sort.v[where.next]
  
  # Calculate v_max for PPV modeling
  v_max <- 2 * v_m - v_m1
  
  # ----------------------------------------
  # Inner Loop 2: for each segment 
  # Within this loop there are two steps: 
  # Step 1: calculate ellipsoid parameters and rotation angles
  # Step 2: Inner loop 3: for each voxel in currentPPV decide if it's inside or outside the ellipsoid
  
  # For all trajectory points except the first one, 
  # calculate ellipsoid parameters and rotation angles along the segment that ends in point i
  # i.e. segment P_i-1 to P_i

  for (i in 2:trajLength) {
    
    # Step 1: calculate ellipsoid parameters and rotation angles
    # Major and minor axes of ellipsoid, a_i, b_i
    dfAll$a[trajIndices[i]] <- (v_max * dfAll$dt[trajIndices[i]])/2 
    dfAll$b[trajIndices[i]] <- sqrt(dfAll$a[trajIndices[i]]^2-(dfAll$d[trajIndices[i]]^2/4))  
    
    # Coordinates of central point on segment P_C: xc_i, yc_i, zc_i
    dfAll$xc[trajIndices[i]] <- (dfAll$X[trajIndices[i]]+dfAll$X[trajIndices[i-1]])/2     
    dfAll$yc[trajIndices[i]] <- (dfAll$Y[trajIndices[i]]+dfAll$Y[trajIndices[i-1]])/2      
    dfAll$zc[trajIndices[i]] <- (dfAll$Z[trajIndices[i]]+dfAll$Z[trajIndices[i-1]])/2 
    
    # Rotation angles alpha and beta 
    dfAll$alpha[trajIndices[i]] <- atan((dfAll$Y[trajIndices[i]]-dfAll$Y[trajIndices[i-1]])/(dfAll$X[trajIndices[i]]-dfAll$X[trajIndices[i-1]]))   
    dfAll$beta[trajIndices[i]] <- asin((dfAll$Z[trajIndices[i]]-dfAll$Z[trajIndices[i-1]])/dfAll$d[trajIndices[i]])
    
    # Rotation direction for beta depends on the direction of movement (up or down) along the segment and dx
    dZ <- dfAll$Z[trajIndices[i]]-dfAll$Z[trajIndices[i-1]]
    dX <- dfAll$X[trajIndices[i]]-dfAll$X[trajIndices[i-1]]
    dfAll$direction[trajIndices[i]] <- dZ * dX
    
  } # Closing inner loop 2: for (i in 2:trajLength)
      
} # Closing outer loop: for (shk in 1:noOfTraj)

# ----------------------------------------
# Calculate the largest extent, depending on the largest b

# Add a buffer in distance to extent to fit ellipsiods around end points
# Buffer is the size of the largest b across all trajectories and segments

distBuffer <- max(dfAll$b)
#distBuffer <- 0.3 * max(distX,distY,distZ)
distX <- maxXcoord - minXcoord
distY <- maxYcoord - minYcoord
distZ <- maxZcoord - minZcoord
minXcoord <- minXcoord - distBuffer
maxXcoord <- maxXcoord + distBuffer
minYcoord <- minYcoord - distBuffer
maxYcoord <- maxYcoord + distBuffer
minZcoord <- minZcoord - distBuffer
maxZcoord <- maxZcoord + distBuffer

# Generate eight points defining extent for visualising in Voxler and export file
xE <- seq(minXcoord,maxXcoord,by=distX+2*distBuffer)
yE <- seq(minYcoord,maxYcoord,by=distY+2*distBuffer)
zE <- seq(minZcoord,maxZcoord,by=distZ+2*distBuffer)
MExtent <- mesh(xE,yE,zE)
xEc <- MExtent$x
yEc <- MExtent$y
zEc <- MExtent$z
WriteCSVtableFromThreeRectangular3DArrays(xEc,yEc,zEc,extent.file)

# Set automatic voxel size depending on max distance between points, so that there
# are 100 voxels on the side that is longest
#voxelSize <- max(distX,distY,distZ)/100

# For one/three segments, we can take a fixed voxel size of 0.2 for testing
voxelSize <- 0.2  

# ------------------------------------------------
# Build three 3D arrays of x, y, z coordinates and initialise the total density volume

startX <- floor(minXcoord/voxelSize)*voxelSize
startY <- floor(minYcoord/voxelSize)*voxelSize
startZ <- floor(minZcoord/voxelSize)*voxelSize
endX <- ceiling(maxXcoord/voxelSize)*voxelSize
endY <- ceiling(maxYcoord/voxelSize)*voxelSize
endZ <- ceiling(maxZcoord/voxelSize)*voxelSize

# No of voxels
xvoxels <- (endX-startX)/voxelSize+1
yvoxels <- (endY-startY)/voxelSize+1
zvoxels <- (endZ-startZ)/voxelSize+1

# Build the volume
x <- seq(startX,endX,by=voxelSize)
y <- seq(startY,endY,by=voxelSize)
z <- seq(startZ,endZ,by=voxelSize)
M <- mesh(x,y,z)
xcoord <- M$x
ycoord <- M$y
zcoord <- M$z

# Initialise the total Accessibility Volume of PPVs as zeros everywhere in a 3D array of the same 
# size as the xcoord/ycoord/vcoord volumes
totalPPV <- array(data=0,dim=dim(xcoord))



# ----------------------------------------
# Restart Loop 1: for each trajectory calculate PPV
print('Loop 2 - PPVs')

for (shk in 1:noOfTraj) {
      #for (shk in 1:1) {
      
  # How long is this trajectory?
  trajIndices <- which(dfAll$TrajID == listOfTrajIDs[shk])
  trajLength <- length(trajIndices)
      
  # Restart inner loop 2: for each segment 
  for (i in 2:trajLength) {
    
    # Step 2: Inner loop 3: for each voxel in currentPPV decide if it's inside or outside the ellipsoid
    # For all trajectory points except the first one, 
    # calculate ellipsoid parameters and rotation angles along the segment that ends in point i
    # i.e. segment P_i-1 to P_i
    for (ii in 1:xvoxels) { print(c(shk,i,ii))
      for (jj in 1:yvoxels){
        for (kk in 1:zvoxels){
          
          # Find coordinates of the voxel
          x_v <- xcoord[ii,jj,kk]
          y_v <- ycoord[ii,jj,kk]
          z_v <- zcoord[ii,jj,kk]
          
          # Transform the coordinates by 1) translation to PC, 2) rotation for alpha, 3)rotation for beta
          # 1) Translation to PC
          x_v_1 <- x_v - dfAll$xc[trajIndices[i]]
          y_v_1 <- y_v - dfAll$yc[trajIndices[i]]
          z_v_1 <- z_v - dfAll$zc[trajIndices[i]]
          
          # 2) Rotation for alpha, in x-y plane, counter-clockwise around z axis
          x_v_2 <- cos(dfAll$alpha[trajIndices[i]])*x_v_1 + sin(dfAll$alpha[trajIndices[i]])*y_v_1
          y_v_2 <- - sin(dfAll$alpha[trajIndices[i]])*x_v_1 + cos(dfAll$alpha[trajIndices[i]])*y_v_1
          z_v_2 <- z_v_1

          # 3)Rotation for beta, rotating in z direction, around y axis, x-z plane needs to stay fixed (y = fixed)
          y_new <- y_v_2
          
          if(dfAll$direction[trajIndices[i]]>=0) { # dz*dx is more than zero -> we turn up or down based on dz
          
            if (dZ>0) { # Moving up: z(i)>z(i-1) -> Rotation 1
              x_new <- cos(dfAll$beta[trajIndices[i]])*x_v_2 + sin(dfAll$beta[trajIndices[i]])*z_v_2
              z_new <- - sin(dfAll$beta[trajIndices[i]])*x_v_2 + cos(dfAll$beta[trajIndices[i]])*z_v_2
            } else if (dZ<0) { # Moving down: z(i)<z(i-1) -> Rotation 2
              x_new <- cos(dfAll$beta[trajIndices[i]])*x_v_2 - sin(dfAll$beta[trajIndices[i]])*z_v_2
              z_new <- sin(dfAll$beta[trajIndices[i]])*x_v_2 + cos(dfAll$beta[trajIndices[i]])*z_v_2
            } else { # Same level: z(i)=z(i-1) 
              x_new <- x_v_2
              z_new <- z_v_2
            }
            
          } else if (dfAll$direction[trajIndices[i]]<0) { # dz*dx is less than zero -> we turn down or up based on dz
            
            if (dZ>0) { # Moving up: z(i)>z(i-1) -> Rotation 2
              x_new <- cos(dfAll$beta[trajIndices[i]])*x_v_2 - sin(dfAll$beta[trajIndices[i]])*z_v_2
              z_new <- sin(dfAll$beta[trajIndices[i]])*x_v_2 + cos(dfAll$beta[trajIndices[i]])*z_v_2
            } else if (dZ<0) { # Moving down: z(i)<z(i-1) -> Rotation 1
              x_new <- cos(dfAll$beta[trajIndices[i]])*x_v_2 + sin(dfAll$beta[trajIndices[i]])*z_v_2
              z_new <- - sin(dfAll$beta[trajIndices[i]])*x_v_2 + cos(dfAll$beta[trajIndices[i]])*z_v_2
            } else { # Same level: z(i)=z(i-1) 
              x_new <- x_v_2
              z_new <- z_v_2
            }
          } # else if
          
          # Is the new point (x_new, y_new, z_new) inside the ellipsoid?
          inside <- (x_new/dfAll$a[trajIndices[i]])^2+(y_new/dfAll$b[trajIndices[i]])^2+(z_new/dfAll$b[trajIndices[i]])^2
          # It is, if inside <= 1, then we assign value 1 to the voxel in totalPPV: this will
          # produce a union of all the new voxels that were inside the PPV for this trajectory:
          # add them to the totalPPV volume:
          if (inside<=1) {totalPPV[ii,jj,kk] <- 1}
          #df
        } # for kk
      } # for jj
    } # for ii
    # Closing inner loop 3: for all voxels in the volume

  } # Closing inner loop 2: for (i in 2:trajLength)
  
} # Closing outer loop: for (shk in 1:noOfTraj)

# Export PPV into a csv file
WriteCSVtableFromFourRectangular3DArrays(xcoord,ycoord,zcoord,totalPPV,output.file)

# Export all parameters for each segment
write.csv(dfAll,parameter.file,row.names=FALSE)

# Calculate volume inside the PPV
nonZeroVoxels <- length(which(totalPPV>0))
volume <- nonZeroVoxels * voxelSize^3
noOfVoxels <- xvoxels*yvoxels*zvoxels
volumeSize <- cbind(xvoxels,yvoxels,zvoxels,voxelSize,noOfVoxels,nonZeroVoxels,volume)
write.csv(volumeSize,volume.file,row.names=FALSE)
