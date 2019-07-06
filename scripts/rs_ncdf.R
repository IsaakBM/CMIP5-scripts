# This code was written by Isaac Brito-Morales (i.britomorales@uq.edu.au)
# Please do not distribute this code without permission.
# NO GUARANTEES THAT CODE IS CORRECT
# Caveat Emptor!

# AIM: Function that reads a 3 dimensional netCDF file (long, lat, time) and returns a list of rasterstacks
# nc = netCDF dile directory
# this function needs to be runned with par_rs_ncdf script
  
  # without interpolation process and for monthly data
  rs.ncdf <- function(nc, model, rcp, outdir, v = "thetao", x = "lon", y = "lat", depth = "lev", rrun = "r1i1p1", yrs = "2006-2100") {
                                                                                      
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
    X <- dim(dat)[1]
    Y <- dim(dat)[2]
    tt <- nc.get.time.series(nc, v = "time", time.dim.name = "time") # from packages ncdf4.helpers&PCICt
      tt <- as.POSIXct(tt)
      tt <- as.Date(tt)
    nc_close(nc)
    rs <- raster(nrow = Y, ncol = X) # Make a raster with the right dims to fill with lat&lon
    # Fix orientation of original data [and then create a raster with this fix orientation and paste deths and time...]
    drs <- data.frame(coordinates(rs))
    # Create rasters stacks of depths for every month
    rs_list <- list() # empty list to allocate results
    st <- stack()
      for (i in 1:length(tt)) {
        dt1 <- rasterFromXYZ(cbind(drs, as.vector(dat[,, i])))
          dt1[]<- ifelse(dt1[] <= -2, NA, dt1[]) # for some models that have weird temperatures (-273)
          dt1[]<- ifelse(dt1[] >= 40, NA, dt1[])
        st <- addLayer(st, flip(dt1, 2))
        print(paste0(i, " of ", length(tt)))
      }
    
    # names(st) <- paste0(as.character(tt))
    names(st) <- seq(as.Date("2006/1/1"), as.Date("2100/12/1"), by = "mon") # to standardize days (data came in same year+month but different days)
    crs(st) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    name.rs <- paste(v, model, rcp, rrun, yrs, ".grd", sep = "_")
    writeRaster(st, paste(outdir, name.rs, sep = ""), overwrite = TRUE)
    
    return(st)
    
  }
  