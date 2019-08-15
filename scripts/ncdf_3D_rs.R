# This code was written by Isaac Brito-Morales (i.britomorales@uq.edu.au)
# Please do not distribute this code without permission.
# NO GUARANTEES THAT CODE IS CORRECT
# Caveat Emptor!

# AIM: Function that reads a 4 dimensional netCDF file (long, lat, depth, time) and returns a list of rasterstacks
# nc = netCDF dile directory
# this function needs to be runned with par_rs_ncdf script
  
  # without interpolation process and for monthly data
  ncdf_3D_rs <- function(nc, v = "thetao", x = "lon", y = "lat", depth = "lev") { # v = "thetao" for RCPs and "temperature" for EN4
                                                                                      # depth = "lev" for RCPs and "depth"for EN4
    library(akima)
    library(raster)
    library(ncdf4)
    library(ncdf4.helpers)
    library(PCICt)
    
    # ifelse(names(nc$variable) == "temperature", "thetao", "thetao") # some thing like that to uniform the criteria
    # Extract data from the netCDF file  
    nc <- nc_open(nc)
    dat <- ncvar_get(nc, v) # x, y, depth, year 
    dat[] <- dat-273 # from Kelvin to Celcius
    depth <- ncvar_get(nc, depth)
    # X <- ncvar_get(nc, x)
    X <- dim(dat)[1]
    # Y <- ncvar_get(nc, y)
    Y <- dim(dat)[2]
    tt <- nc.get.time.series(nc, v = "time", time.dim.name = "time") # from packages ncdf4.helpers&PCICt
      tt <- as.POSIXct(tt)
      tt <- as.Date(tt)
    nc_close(nc)
    # range(Y); range(X); range(dat, na.rm = TRUE)
    rs <- raster(nrow = Y, ncol = X) # Make a raster with the right dims to fill with lat&lon
    # Fix orientation of original data [and then create a raster with this fix orientation and paste deths and time...]
    drs <- data.frame(coordinates(rs))
    # drs$x <- drs$x + 180 # to get 0-360 # for GISS-E2-H & GISS-E2-H-CC models
    # drs$x <- ifelse(drs$x < 180, drs$x, drs$x-360) # for GISS-E2-H & GISS-E2-H-CC models
    
    # Create rasters stacks of depths for every month
    rs_list <- list() # empty list to allocate results
    for (i in 1:length(depth)) {
      for (j in 1:length(tt)) {
        if(j == 1) {
          dt1 <- rasterFromXYZ(cbind(drs, as.vector(dat[,, i, j]))) # create a raster from the original raster (already fixed it) with the 1 depth and 1 month
            dt1[]<- ifelse(dt1[]<= -2, NA, dt1[]) # for some models that have weird temperatures (0-300K)
            dt1[]<- ifelse(dt1[] >= 40, NA, dt1[])
          # extent(dt1) <- extent(-180, 180, -84, 89) # for GISS-E2-H & GISS-E2-H-CC models
          
          rs.final <- flip(dt1, 2)
        } else {
          dt1 <- rasterFromXYZ(cbind(drs, as.vector(dat[,, i, j]))) # create a raster from the original raster (already fixed it) with the 1 depth and 1 month
            dt1[]<- ifelse(dt1[]<= -2, NA, dt1[]) # for some models that have weird temperatures (0-300K)
            dt1[]<- ifelse(dt1[] >= 40, NA, dt1[])
          # extent(dt1) <- extent(-180, 180, -84, 89) # for GISS-E2-H & GISS-E2-H-CC models
          rs.final <- stack(rs.final, flip(dt1, 2)) 
        }
        print(paste0(i, " of ", length(tt)))
      }
      names(rs.final) <- paste0(as.character(tt))
      rs_list[i] <- rs.final # save it on a list
    }
    
    # Working to get ocean layers 
      # 1.Unlist the list new element ([model files] x [depth layers] and [yearly stack layers])
        rs_final <- unlist(rs_list)
        names(rs_final) <- as.character(round(depth)) # name every element according with depth
  
    return(rs_final)
    
  }
  