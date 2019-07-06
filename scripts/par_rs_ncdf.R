# This code was written by Isaac Brito-Morales (i.britomorales@uq.edu.au)
# Please do not distribute this code without permission.
# NO GUARANTEES THAT CODE IS CORRECT
# Caveat Emptor!

# AIM: A parallelized function of the rs_ncdf05 script that returns a list of RasterStacks by depth
# nc = netCDF dile directory
# this function needs to be runned with par_rs_ncdf script


# Libraries
  library(doParallel)
  library(foreach)
  library(raster)

# 1. Define the structure that will be parallelized and loop all files
  rm(list = ls())
  # Establish directory
    source("scripts/rs_ncdf05.R")
    files.nc <- list.files("netCDF02/MRI-CGCM3_rcp45_b(1deg)/", pattern = "*.nc", full.names = TRUE) #from 2006-2100
    files_list <- list() # to allocate results
  # Begin the parallel structure
    UseCores <- detectCores() -1 # Define how many cores
    cl <- makeCluster(UseCores)  
    registerDoParallel(cl) # Register CoreCluster
  # Create big RasterStack list
    rs_list <- foreach(i = 1:length(files.nc)) %dopar% { # ~20 minutes for 19 files at 0.25deg (takes less at 1deg)
      library(raster)
      library(ncdf4)
      
      files_list[i] <- rs.ncdf03(files.nc[i])
      # here the other function?
    }
  
  stopCluster(cl) # alwasys stop the cluster

# 2. Manipulate RasterStack list 
  # Function to stack lists of rasters
    rs_stack <- function(data, depth_layer) { # A generic function that uses staking function (add argument time_period = 2006:2100)
      
      staking <- function(dat) { # function that allows to stack a list of RasterStacks
        d <- stack(dat[[1]])
        for (i in 2:length(dat)) {
          d <- addLayer(d, dat[[i]])
        }
        names(d)<- paste(rep(2006:2100, each = 1)) #1960:2017 EN4 #2005-2110 for HadGEM2-ES and CC # 1955-2005 for historical
        crs(d) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
        return(d)
      }
      
      slayer.list <- list()
      for (j in 1:length(data)) {
        slayer <- data[[j]][[depth_layer]] # it will stak every depth from the list
        slayer.list[[j]] <- slayer
      }
      rs.slayer <- staking(slayer.list) # this function uses the staking function previously defined
      return(rs.slayer)
    }

# 3. Loop to get a final list in which every element is a depth and inside there is a temporal yearly RasterStack
  n_depths <- round(as.numeric(names(rs_list[[1]]))) # depths
  layers <- seq(1, length(n_depths), 1) # sequence to loop
  final_list <- list() # allocate results
  for (k in 1:length(layers)) {
    final_list[k] <- rs_stack(data = rs_list, depth_layer = k)
  }
  names(final_list) <- as.character(n_depths)
  
# 4. Define ocean's layers/zones by depth 
  ocean_layers <- function(data) { # data is a big raw rasterstacks (final_list) change for data in case of memory (add argument time_period = 2006:2100)
    
    # Standard ocean's layers
      s_list <- list()
      m_list <- list()
      b_list <- list()
      bt_list <- list()
        for (i in 1:length(data)) {
          if (as.numeric(names(data[i])) >= 0 & as.numeric(names(data[i])) <= 6) { # surface (5 m)
            s_list <- unlist(list(s_list, data[i]))
            # s_list2 <- s_list[2:length(s_list)]
          } else if (as.numeric(names(data[i])) > 200 & as.numeric(names(data[i])) <= 1000) { # mesopelagic (up to 1000 m)
            m_list <- unlist(list(m_list, data[i]))
            # m_list2 <- m_list[2:length(m_list)]
          } else if (as.numeric(names(data[i])) > 1000 & as.numeric(names(data[i])) <= 2000) { # bathy EMUs (up to 2000 m)
            b_list <- unlist(list(b_list, data[i]))
            # b_list2 <- b_list[2:length(b_list)]
          } else if (as.numeric(names(data[i])) > 2000) { # bottom EMUs (higher than 2000 m)
            bt_list <- unlist(list(bt_list, data[i]))
            # bt_list2 <- bt_list[2:length(bt_list)]
          }
        }
    
    # Alternatives ocena's layers
      s2_list <- list()
      b2_list <- list()
      bt2_list <- list()
        for (j in 1:length(data)) {
          if (as.numeric(names(data[j])) >= 0 & as.numeric(names(data[j])) <= 200) { # surface (up to 200 m)
            s2_list <- unlist(list(s2_list, data[j]))
            # s2_list2 <- s2_list[1:length(s2_list)]
          } else if (as.numeric(names(data[j])) > 1000 & as.numeric(names(data[j])) <= 4000) { # bathy NOAA (up to 4000 m)
            b2_list <- unlist(list(b2_list, data[j]))
            # b2_list2 <- b2_list[1:length(b2_list)]
          } else if (as.numeric(names(data[j])) > 4000) { # bottom NOAA (higher than 4000 m)
            bt2_list <- unlist(list(bt2_list, data[j]))
            # bt2_list2 <- bt2_list[1:length(bt2_list)]
          }
        }
      
        bt3_list <- list()
          for (k in 1:length(data)) {
            if (as.numeric(names(data[k])) > 1000) { # bathy + bottom ( any layer > 1000 m)
              bt3_list <- unlist(list(bt3_list, data[k]))
              # s2_list2 <- s2_list[1:length(s2_list)]
            } 
          }
    
    # Function to establish depth_weights
      depthW <- function(top, slices, bottom) {
        if(length(slices) > 1) {
          l <- diff(c(top, slices, bottom)) # Differences between consecutive values
          ll <- rep(.5, length(l)) # A vector to divide differenes in half
            ll[1] <- ll[length(l)] <- 1 # Don't divide differences between top or bottom from adjacent values
          ll <- l*ll # Do the division for all non-top/bottom differences
          d <- numeric(length(ll)) # Catch results
          for(i in 1:(length(ll))) {d[i] <- ll[i] + ll[i+1]} # Add consecutive layer widths
            return(as.numeric(na.omit(d))) # Dump any NAs
          } else {
            return(1) # If there is only one slice in the depth category, no need to weight
          }
      }
      
      # Use the function to create weights by each depth layer category
        dw_s <- depthW(0, as.numeric(names(s_list)), 6)
        dw_m <- depthW(200, as.numeric(names(m_list)), 1000)
        dw_b <- depthW(1000, as.numeric(names(b_list)), 2000)
        dw_bt <- depthW(2000, as.numeric(names(bt_list)), 6000)
          dw_s2 <- depthW(0, as.numeric(names(s2_list)), 200)
          dw_b2 <- depthW(1000, as.numeric(names(b2_list)), 4000)
          dw_bt2 <- depthW(4000, as.numeric(names(bt2_list)), 6000)
            dw_bt3 <- depthW(1000, as.numeric(names(bt3_list)), 6000) # 6500 for MRI-CGCM3 model
        
    # Create rasterstacks of models by "number of years"
      fx_rs <- function(olayer, depth_weights) {
        # Define how many years
          names.yrs <- paste("X", as.character(seq(2006, 2100)), sep = "")
        # Loop through
          st <- stack()
          rs_multi02 <- list()
          for (i in 1:length(names.yrs)) {
            for (j in 1:length(olayer)) {
              if (j == 1) {
                dt <- subset(olayer[[j]], names.yrs[i]) # subset every rasterstacks element of the list
                st1 <- dt
              } else {
                dt <- subset(olayer[[j]], names.yrs[i]) # same again!
                st1 <- stack(st1, dt)
              }
              print(paste0(i, " of ", length(names.yrs)))
            }
            rs_multi02[i] <- st1
          }
      
      # Estimate mean of depths to get ocean's layers
        model.mean <- list()
        for (l in 1:length(rs_multi02)) {
          model.mean[l] <- stackApply(rs_multi02[[l]], indices = rep(1, nlayers(rs_multi02[[l]])), 
                                      fun = function(x, ...) raster::weighted.mean(x, depth_weights, na.rm = TRUE))
        }
      # To stack those raster
        staking <- function(dat) {
          d <- stack(dat[[1]])
          for (m in 2:length(dat)) { 
            d <- addLayer(d, dat[[m]])
          }
          names(d)<- paste(rep(2006:2100, each = 1))
          crs(d) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
          return(d)
        }
      
      # Final stack
        final.stack <- staking(model.mean)
        return(final.stack)
    }
    
    
    # DO THE ANALYSIS FOR ALL LAYERS [to much memory for high-res rasters]
      rs_stack <- list(fx_rs(olayer = s_list, depth_weights = dw_s), fx_rs(olayer = s2_list, dw_s2), # surface (5 m) & surface (up to 200 m)
                       fx_rs(olayer = m_list, depth_weights = dw_m), # mesopelagic (up to 1000 m)
                       fx_rs(olayer = b_list, depth_weights = dw_b), fx_rs(olayer = b2_list, depth_weights = dw_b2), # bathy EMUs (up to 2000 m) & # bathy NOAA (up to 4000 m)
                       fx_rs(olayer = bt_list, depth_weights = dw_bt), fx_rs(olayer = bt2_list, depth_weights = dw_bt2), 
                       fx_rs(olayer = bt3_list, depth_weights = dw_bt3)) # bottom EMUs (higher than 2000 m), # bottom NOAA (higher than 4000 m), # bathy + bottom ( any layer > 1000 m)
    
    return(rs_stack)
    
  }
  
  system.time(model <- ocean_layers(final_list)) # ~ 40 minutes(server) ~ 15 minutes(mac)
    # plot(model[[2]]$X2005, useRaster = FALSE) # check if all is OK? YES
    # summary(model[[2]]$X2005[])
    # summary(model[[3]]$X2005[])
    # summary(model[[5]]$X2005[])
    # summary(model[[7]]$X2005[])
    
  
  # Function to write climate models
    write.rs <- function(data, outdir, var, model, rcp, rrun = "r1i1p1", yrs = "2006-2100") {
      for (i in 1:length(data)) {
        single <- data[[i]]
        name.rs <- paste(var, i, model, rcp, rrun, yrs, ".grd", sep = "_")
        writeRaster(single, paste(outdir, name.rs, sep = ""), overwrite = TRUE)
      }
    }
  
    write.rs(data = model, outdir = "Data03/MRI-CGCM3_rcp45_b(1deg)/", var = "thetao", model = "MRI-CGCM3", rcp = "rcp45")
      