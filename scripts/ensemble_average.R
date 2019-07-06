

library(raster)
library(rasterImage)
library(VoCC)
library(sf)

# after using this scritp, save your output and close R (it takes a lot of memory from R and it will slow your sesion)
# this function was thinking and wrote for a parallel structure (ask for the other code...)
model.mean <- function(path) {
  # Packages
    library(raster)
  
  # Establish path for files
    files_multi <- list.files(path = path, pattern = ".grd", full.names = TRUE)
  # Reads and saves every model into a list class object
    rs_multi <- list()
      for (k in 1:length(files_multi)) {
        names.yrs <- names(files_multi[[k]])
        rs <- stack(files_multi[k]) # first step read every file
          if (res(rs)[1] != 1) { # == 0.25 if models came in high-res, standardised them at 1x1 degree of latitude aggregating by the mean
            rs1 <- raster::aggregate(rs, 4)
          } else if (nlayers(rs) > 1140) { # models with more years? I only want 2006-2100
              rs1 <- raster::subset(rs, names.yrs)
          } else {rs1 <- rs}
        
        rs_multi[k] <- rs1 # add those rasters to a list element
      }    

  # Create rasterstacks of models by "number of years"
    names.yrs <- names(rs_multi[[1]])
    st <- stack()
    rs_multi02 <- list()
      for (i in 1:length(names.yrs)) {
        for (j in 1:length(rs_multi)) {
          if (j == 1) {
          dt <- subset(rs_multi[[j]], names.yrs[i]) # subset every rasterstacks element of the list
          st1 <- dt
          } else {
          dt <- subset(rs_multi[[j]], names.yrs[i]) # same again!
          st1 <- stack(st1, dt)
          }
        print(paste0(i, " of ", length(names.yrs)))
        }
      rs_multi02[i] <- st1
      }
  
  # Estimate mean and SD by year to get our model-ensemble means
    model.mean <- list()
    model.sd <- list()
      for (l in 1:length(rs_multi02)) {
        model.mean[l] <- stackApply(rs_multi02[[l]], indices = rep(1, nlayers(rs_multi02[[l]])), fun = mean)
        model.sd[l] <- stackApply(rs_multi02[[l]], indices = rep(1, nlayers(rs_multi02[[l]])), fun = sd)
      }
        
  # Function to stack a list of rasters     
    staking <- function(dat) {
      d <- stack(dat[[1]])
        for (m in 2:length(dat)) { 
          d <- addLayer(d, dat[[m]])
        }
      names(d)<- names.yrs
      crs(d) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
      return(d)
    }
        
  # Create a ensemble-model means raster and the SD
    final.mean <- staking(model.mean)
    final.sd <- staking(model.sd)

    final.ls <- list(final.mean, final.sd)
    return(final.ls)
    
}   
        
