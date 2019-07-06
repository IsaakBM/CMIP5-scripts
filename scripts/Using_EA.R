# files.nc <- list.files("netCDF/rcp85/", pattern = "*.nc", full.names = TRUE) #from 2006-2100
# names(files.nc) <- c("CanESM2", "CNRM-CM5", "CSIRO-Mk3", "GFDL-CM3", "GFDL-ESM2G",
#                      "GFDL-ESM2M", "GISS-E2-H", "GISS-E2-R", "MIROC5", "MPI-ESM-LR", "MRI-CGCM3")
# 
# for (i in 1:length(files.nc)) {
#   test <- rs.ncdf(nc = files.nc[i], model = as.character(names(files.nc[i])), rcp = "rcp85", outdir = "griFiles/rcp85/")
# }

# using previous scripts
source("scripts/ensemble_average.R")
sst_rcp45 <- model.mean("griFiles/rcp45/")
writeRaster(sst_rcp45[[1]], "griFiles/thetao_surface_MonthlyEnsembleMean_rcp45_r1i1p1_2006-2100.grd")
writeRaster(sst_rcp45[[2]], "griFiles/thetao_surface_MonthlyEnsembleSD_rcp45_r1i1p1_2006-2100.grd")



