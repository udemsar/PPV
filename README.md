# PPV
Potential Path Volume - a geometric estimator for space use in 3D

This is an implementation of the algorithm for the Potential Path Volume (PPV), a 
geometric estimator of space use around a trajecotry where the location is measured
in three physical dimensions.

The code is free and open and provided as supplementary information to the following paper:

Demšar U and Long JA, 2018, Potential Path Volume (PPV) - a geometric estimator for space use in 3D
(Currently under review)

How to run the code:
The main file is PPVolumesFrom4DTrajectories.R, which runs a worked example on a trajectory 
with one segment, provided as One4Dsegment.csv. Another example is a trajecotry with three 
segments, given in the file Three4Dsegments.csv.

The repository also contains results from these two examples (all remaining .csv files) 
as well as a Voxler file PPV_Test_Segments.voxb (zipped), which visualises the results using the
3D volumetric visualisation environment Voxler (Golden Software).

A faster implementation of the code is integrated into the wildlifeTG R package, available at: 
https://github.com/jedalong/wildlifeTG
