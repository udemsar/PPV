# ------------------------------------------------
# DistanceTwo3Dpoints.R
# ------------------------------------------------
# Code provided as is and can be used or modified freely. 
#
# We would however appreciate if you cite the following paper:
#
# Demsar U and Long J, Time-Geography in Four Dimensions: 
# Potential Path Volumes around 3D trajectories, GIScience 2016
# ------------------------------------------------
# Author: Urska Demsar
# University of St Andrews
# St Andrews, Scotland, UK
# http://udemsar.com
# Date: 11 Aug 2016
# ------------------------------------------------


# ------------------------------------------------
# Calculates distance between 2 points in 3D, 
# given as (x,y,z) and (x1,y1,z1).
# ------------------------------------------------

DistanceTwo3Dpoints <- function(x,y,z,x1,y1,z1) {
  
  d2 <- (x-x1)^2+(y-y1)^2+(z-z1)^2

  f <- sqrt(d2)
  
  return(f)  
  }