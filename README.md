# Spatial_RedSnapper_OM_and_SM
Spatially explicit population simulator, and data simulator, for Red Snapper which was used in Fisch et al. (2021) - Fisheries Research and Fisch et al. (2023) - IJMS

Data_Processing.R is a file which processes some raw data into files then used in the Spatially explict OM (e.g., substrate and depth files)

SpatialModel_Function.R runs the Spatially Explicit Operating Model 
   Note that it is currently set up to run 50 unfished years, and either 40 or 80 fished years (formulation from Fisch et al., 2023)
   This could be changed, but the operator would need to alter some aspects of the function (e.g., effort for the time series)

Get_Data.R runs a sampling model for a Spatially Explicit Operating Model that has already been run. 

