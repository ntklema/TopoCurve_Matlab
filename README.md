This code was used to generate the results presented in 

Klema, N., Karlstrom, L., and Roering, J.: Discrete differential geometry of fluvial landscapes, https://doi.org/10.5194/egusphere-2025-4431, 2025

The code uses a formal surface theory approach to calculate geometry metrics on gridded DEMs that are not sensative to projection distortion.

Geometric calculations are done within the CurveObj object class, which also stores geometry metrics that are derived through various binning and analysis functions. 
The script ESurf_Figs.m shows the basic workflow used to generate manuscript figures, and will give basic intuition around application of the CurveObj class. The package 
relies heavily on the TopoToolbox package for basic processing of geospatial datasets.  Topotoolbox can be downloaded at https://topotoolbox.wordpress.com/download/ and 
can be added through the built in MATLAB package manager.

Please reach out to ntklema@fortlewis.edu with questions.
