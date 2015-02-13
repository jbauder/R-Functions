# --- HOMERANGE FUNCTIONS ------------------------------------------------------
# General functions for creation, analysis, and simulation of home ranges
# ------------------------------------------------------------------------------

# CenterXYInCell Function ------------------------------------------------------

###  Centers x and y values into the center of a raster cell based on the xmin, 
###    ymin, and cellsize parameters
###  Usage: CenterXYInCell(x, y, xmin, ymin, cellsize)
###  Arguments: x = x value
###             y = y value
###             xmin = minimum value of x that establishes the grid arrangement
###             ymin = minimum value of y that establishes the grid arrangement
###             cellsize = cell size in units of the x and y values
###  Returns: vector of centered x and y (i.e., (x,y))   
###  Notes:
###  Blake Massey
###  2015.01.07

CenterXYInCell <- function(x, 
                           y, 
                           xmin, 
                           ymin, 
                           cellsize) { 
    x <- xmin + (floor(((x-xmin)/cellsize))*cellsize) + (cellsize/2)
    y <- ymin + (floor(((y-ymin)/cellsize))*cellsize) + (cellsize/2)
    output <- c(x, y)
    return(output)
}

# ConvertAngle Function --------------------------------------------------------

###  Converts a radian angle that is outside of the Unit Circle range (i.e., 
###    [0, 2pi]) to the equivalent value that within the Unit Circle range. 
###  Usage: ConvertAngle(x)
###  Arguments: x = radian value
###  Returns: vector of one value 
###  Notes: Used in several simulations so that the angles never go outside of 
###    the Unit Circle range when angles are being added or subtracted from each 
###    other over and over during a simulation. 
###  Blake Massey
###  2015.01.07

ConvertAngle <- function(x) {
  if (x > 2*pi) x <- x-(2*pi)
  if (x < 0) x <- x+(2*pi)
  return(x)
}

# CreateHomeRangeKernels Function ----------------------------------------------

###  Creates RasterLayers of homerange kernels based on home centroids
###  Usage: CreateHomeRangeKernels(df_all, df_home, base, max_r, 
###   home_inflection, home_scale, avoid_inflection, avoid_scale, output_dir, 
###   write_dist, write_homerange)
###  Arguments: df_all = dataframe of all homerange centroids 
###             df_home = dataframe of homerange centroids to calculate
###               homerange kernels (these must be a subset of the df_all 
###               dataframe), Default is to use 'df_all' dataframe
###             base = base Raster that sets the projection, extent, and 
###               dimensions of the study area 
###             max_r = maximum radius to calculate the homerange raster from 
###               each df_home centroid
###             home_inflection = inflection point of the Logistic function that
###               governs the home kernel
###             home_scale = scale parameter of the Logistic function that 
###               governs the scale parameter
###             avoid_inflection = inflection point of the Logistic function
###               that governs the conspecific avoidance kernel
###             avoid_scale = scale parameter of the Logistic function that 
###               governs the conspecific avoidance kernel  
###  Returns: A list containing homerange kernel Rasters for all the df_home 
###    centroids
###  Notes: 
###  Blake Massey and Javan Bauder
###  2015.02.06 

CreateHomeRangeKernels <- function(df_all,
                                   df_home = df_all,
                                   base,
                                   max_r,
                                   home_inflection, 
                                   home_scale,
                                   avoid_inflection, 
                                   avoid_scale,
                                   output_dir,
                                   write_dist=FALSE,
                                   write_homerange=FALSE) {
  cellsize <- res(base)[1]
  max_r_cells <- ceiling(max_r/cellsize)
  xmin <- xmin(base)
  ymin <- ymin(base)
  df_sp <- SpatialPointsDataFrame(df_all[,c("x","y")], df_all, 
    proj4string=crs(base))
  writeLines(noquote(paste("Calculating global distance")))
  if (write_distance == TRUE) {
    ifelse(exists("i"), i <- i , i <- 1)
    filename <- paste0(output_dir,"/global_dist_", sprintf("%03d", i), ".tif")
    global_dist <- distanceFromPoints(base, df_sp, filename=filename, 
        overwrite=TRUE)
    writeLines(noquote(paste("Writing:", filename)))
  } else {
    global_dist <- distanceFromPoints(base, df_sp)
  }
  homerange_ids <- df_home$id
  total <- length(homerange_ids)
  homerange_kernels <- as.list(setNames(rep(NA,length(homerange_ids)), 
    homerange_ids), homerange_ids)   
  for (j in 1:nrow(df_home)) {
    writeLines(noquote(paste("Calculating homerange", j, "of", total)))
    home <- df_home[j,]
    home_sp <- SpatialPointsDataFrame(home[,c("x","y")], home, 
      proj4string=crs(base))
    xy <- CenterXYInCell(home_sp@coords[,"x"], home_sp@coords[,"y"], 
      xmin, ymin, cellsize)
    cell_extent <- extent(xy[1]-(cellsize/2), xy[1]+(cellsize/2), xy[2]-
      (cellsize/2), xy[2]+(cellsize/2))
    cell <- setValues(raster(cell_extent, crs=projection(base), res=cellsize),j)
    home_ext <- extend(cell, c(max_r_cells, max_r_cells), value=NA)
    home_dist <- distanceFromPoints(home_ext, home[,c("x","y")])
    home_kern <- calc(home_dist, fun = function(x){(1/(exp((-(x - 
      home_inflection)) / home_scale) + 1))})
    global_dist_crop <- crop(global_dist, home_dist)
    cent_dist <- overlay(home_dist, global_dist_crop, fun=function (x,y)
      {ifelse(x != y, NA, x)})
    cent_bounds <- boundaries(cent_dist)
    cent_bounds <- subs(cent_bounds, data.frame(from=c(0,1), to=c(NA,1)))
    edge_dist_abs <- distance(cent_bounds)
    edge_dist <- overlay(cent_dist, edge_dist_abs, fun=function(x,y)
      {ifelse(!is.na(x), y*-1, y*1)})    
    avoid_kern <- calc(edge_dist, fun = function(x){(1/(exp((-(x - 
      avoid_inflection)) / avoid_scale) + 1))})
    homerange_kern <- overlay(avoid_kern, home_kern, fun=function(x,y){
      p <- y*x #    y+x-1
    #       ifelse(p<=0, 0, p)
      })
    if (write_distance == TRUE) {
      k <- home[,"id"]
      filename <- paste0(output_dir,"/homerange_", k, "_",sprintf("%03d",
        i), ".tif")  
      writeLines(noquote(paste0("Writing: ", filename)))
      writeRaster(homerange_kern, filename=filename, overwrite=TRUE)
    }
    homerange_kernels[[j]] <- homerange_kern
    names(homerange_kernels[[j]]) <- paste0("homerange_", home[,"id"])
  }
  return(homerange_kernels)
}

# CreateHomeRangeKernelRaster Function (DEPRECATED) ----------------------------

###  Create a homerange kernel raster based on homerange and conspecific centers
###  Usage: CreateHomeRangeKernelRaster(home, cons, base, max_dist, location, 
###    scale)
###  Arguments: home = two-column dataframe with x,y location 
###             cons = two-column dataframe with x,y locations
###             base = raster that is used to rasterize() the input data
###             max_dist = maximum distance that the homerange will extent
###             location = location of logisitic function that determines the
###               function of the kernel interactions, default is 10 
###             scale = scale of logistic function that determines the function
###               of the kernel interactions, default is -15
###  Returns: raster object  
###  Notes:
###  Blake Massey
###  2015.01.27

# CreateHomeRangeKernelRaster <- function(home, 
#                                         cons, 
#                                         base, 
#                                         max_dist, 
#                                         location = 10, 
#                                         scale = -15) {
#   suppressPackageStartupMessages(require(raster))
#   source('C:/Work/R/Functions/gis.R')
#   home_df <- data.frame(id=1, x=home[,1], y=home[,2], stringsAsFactors=FALSE)
#   cons_df <- data.frame(id=99, x=cons[,1], y=cons[,2], stringsAsFactors=FALSE)
#   home <- CreateRasterFromPointsAndBase(home_df, "id", "x", "y", base=base)
#   cons <- CreateRasterFromPointsAndBase(cons_df, "id", "x", "y", base=base)
#   home_dist <- distance(home)
#   cons_dist <- distance(cons)
#   home_max <- reclassify(home_dist, c(max_dist, Inf, NA))
#   cons_max <- reclassify(cons_dist, c(max_dist, Inf, NA))
#   home_one <- (home_max-cellStats(home_max, 'min')) /
#     (cellStats(home_max, 'max') - cellStats(home_max, 'min'))
#   cons_one <- (cons_max-cellStats(cons_max, 'min')) /
#     (cellStats(cons_max, 'max') - cellStats(cons_max, 'min'))
#   home_kernel <- calc(home_one, fun=function(x) exp(location + scale*x) /
#     (1 + exp(location + scale*x)))
#   cons_kernel <- calc(cons_one, fun=function(x) exp(location + scale*x) / 
#     (1 + exp(location + scale*x)))
#   home_kernel[is.na(home_kernel)] <- 0
#   cons_kernel[is.na(cons_kernel)] <- 0
#   
#   homerange_kernel <- overlay(home_kernel, cons_kernel, fun=function(x,y) 
#     {return(x-y)})
#   
#   homerange_kernel[homerange_kernel <= 0] <- 0
#   return(homerange_kernel)
# }

# CreateParetoKernel Function --------------------------------------------------

###  Creates a probability matrix based on a Pareto distribution.
###  Usage: CreateParetoKernel(scale, shape, max_r, cellsize)
###  Arguments: scale = scale parameter of Pareto distribution
###             shape = shape parameter of Pareto distribution
###             max_r = maximum radius of kernel in meters, default = 100
###             cellsize = cell size in meters, default = 1
###  Returns: RasterLayer  
###  Notes:
###  Javan Bauder and Blake Massey
###  2015.01.03

CreateParetoKernel <- function(scale, 
                               shape, 
                               max_r = 100, 
                               cellsize = 1) {
  require("VGAM")
  max_r_cells <- max_r/cellsize
  size = ceiling(max_r_cells) * 2 + 1
  center = ceiling(max_r_cells) + 1
  kernel <- new("matrix", 0, size, size)
  for (i in 1:size) for (j in 1:size) {
    r = sqrt((i - center)^2 + (j - center)^2) * cellsize  
    if (r <= max_r) 
      kernel[i, j] <- VGAM::dgpd(r,scale=scale,shape=shape,log=FALSE)
  }
  kernel[center, center] <- 1/scale  
  kernel <- kernel / sum(kernel)
  # This last part deletes the cells at the edge if they are all zero
  if (all(kernel[1, ] == 0, kernel[, 1] == 0, 
          kernel[nrow(kernel),] == 0, kernel[, ncol(kernel)] == 0)) {
    kernel <- kernel[2:(nrow(kernel) - 1), 2:(ncol(kernel) - 1)]
  }
  return(kernel)
}

# CreateRedistKernel Function -------------------------------------------------

###  Create a redistribution kernel matrix based on a wrapped Cauchy 
###    distribution for direction and a Pareto distribution for distance.
###  Usage: CreateRedistKernel(max_r, cellsize, mu, rho, shape, scale, 
###    ignore_cauchy, ignore_pareto)
###  Arguments: max_r = maximum radius of kernel in meters, default = 300
###             cellsize = cell size in meters, default = 30
###             mu = mu parameter of wrapped Cauchy distribution, 0 radians is 
###               due east because everything is based on the Unit Circle
###             rho = rho parameter of wrapped Cauchy distribution
###             shape = shape parameter of Pareto distribution
###             scale = scale parameter of Pareto distribution
###             ignore_cauchy = logical, removes cauchy kernel's contribution to 
###               output raster. Default is FALSE.
###             ignore_pareto = logical, removes pareto kernel's contribution to
###               output raster. Default is FALSE.
###  Returns: matrix
###  Notes:
###  Javan Bauder and Blake Massey
###  2015.01.03

CreateRedistKernel <- function(max_r = 300,
                               cellsize = 30,
                               mu,
                               rho,
                               shape,
                               scale, 
                               ignore_cauchy = FALSE,
                               ignore_pareto = FALSE) {
  suppressPackageStartupMessages(require(circular))
  suppressPackageStartupMessages(require(raster))
  suppressPackageStartupMessages(require(texmex))
  AngleToPoint <- function(origin_x, 
                         origin_y, 
                         target_x, 
                         target_y){
    dx <- c(target_x - origin_x)
    dy <- c(target_y - origin_y)
    abs_angle <- atan2(dy, dx)
    abs_angle <- ifelse(abs_angle < 0, (2*pi) + abs_angle, abs_angle)
  } 
  # Create the empty kernel objects
  max_r_cells <- max_r/cellsize
  size <- ceiling(max_r_cells) * 2 + 1
  center <- ceiling(max_r_cells) + 1
  wrpc_kernel <- new("matrix", 0, size, size)
  gpd_kernel <- new("matrix", 0, size, size)
  for (i in 1:size){
    for(j in 1:size){
      r = sqrt((i - center)^2 + (j - center)^2) * cellsize  
      b = AngleToPoint(center, center, j, i)
      if(r <= max_r){
        wrpc_kernel[i, j] <- round(suppressWarnings(dwrappedcauchy(b, mu=mu, 
          rho=rho)), 5)
        gpd_kernel[i, j] <- dgpd(r, sigma=scale, xi=shape, log=FALSE)
      }
    }
  } 
  wrpc_kernel <- apply(wrpc_kernel, 2, rev)
  gpd_kernel[center, center] <- 1/scale
  # This last part deletes the cells at the edge if they are all zero
  if (all(wrpc_kernel[1, ] == 0, wrpc_kernel[, 1] == 0, 
    wrpc_kernel[nrow(wrpc_kernel),] == 0, wrpc_kernel[, ncol(wrpc_kernel)] ==0))
    wrpc_kernel <- wrpc_kernel[2:(nrow(wrpc_kernel) - 1), 2:(ncol(wrpc_kernel) 
      - 1)]
  if (all(gpd_kernel[1, ] == 0, gpd_kernel[, 1] == 0, 
    gpd_kernel[nrow(gpd_kernel),] == 0, gpd_kernel[, ncol(gpd_kernel)] == 0))
    gpd_kernel <- gpd_kernel[2:(nrow(gpd_kernel) - 1), 2:(ncol(gpd_kernel) - 1)]
  # Multiply the two kernels together and re-normalize
  if (ignore_cauchy) wrpc_kernel <- 1
  if (ignore_pareto) gpd_kernel <- 1
  redist_kernel <- gpd_kernel*wrpc_kernel
  redist_kernel <- redist_kernel/sum(redist_kernel)
  return(redist_kernel)
}

# CreateRedistKernelVonMises Function ------------------------------------------

###  Create a redistribution kernel matrix based on a generalized von Mises
###    distribution for direction and a Pareto distribution for distance.
###  Usage: CreateRedistKernel2(max_r, cellsize, mu1, mu2, kappa1, kappa2,
###    shape, scale, ignore_mises, ignore_pareto)
###  Arguments: max_r = maximum radius of kernel in meters, default = 300
###             cellsize = cell size in meters, default = 30
###             mu1 = mu1 parameter of generalized von Mises distribution, 
###               0 radians is due east because everything is based on the 
###               Unit Circle
###             mu1 = mu2 parameter of generalized von Mises distribution
###             kappa1 = first concentration parameter of the generalized von 
###               Mises distribution
###             kappa2 = second concentration parameter of the generalized von 
###               Mises distribution
###             shape = shape parameter of Pareto distribution
###             scale = scale parameter of Pareto distribution
###             ignore_vonmises = logical, removes von Mises kernel's 
###               contribution to output raster. Default is FALSE.
###             ignore_pareto = logical, removes pareto kernel's contribution to
###               output raster. Default is FALSE.
###  Returns: matrix
###  Notes:
###  Javan Bauder and Blake Massey
###  2015.01.03

CreateRedistKernelVonMises <- function(max_r = 300,
                                       cellsize = 30,
                                       mu1,
                                       mu2,
                                       kappa1,
                                       kappa2,
                                       shape,
                                       scale, 
                                       ignore_vonmises = FALSE,
                                       ignore_pareto = FALSE) {
  suppressPackageStartupMessages(require(circular))
  suppressPackageStartupMessages(require(raster))
  suppressPackageStartupMessages(require(texmex))
  source('C:/Work/R/Functions/pars.R')
  AngleToPoint <- function(origin_x, 
                         origin_y, 
                         target_x, 
                         target_y){
    dx <- c(target_x - origin_x)
    dy <- c(target_y - origin_y)
    abs_angle <- atan2(dy, dx)
    abs_angle <- ifelse(abs_angle < 0, (2*pi) + abs_angle, abs_angle)
  } 
  # Create the empty kernel objects
  max_r_cells <- max_r/cellsize
  size <- ceiling(max_r_cells) * 2 + 1
  center <- ceiling(max_r_cells) + 1
  wrpc_kernel <- new("matrix", 0, size, size)
  gpd_kernel <- new("matrix", 0, size, size)
  for (i in 1:size){
    for(j in 1:size){
      r = sqrt((i - center)^2 + (j - center)^2) * cellsize  
      b = AngleToPoint(center, center, j, i)
      if(r <= max_r){
        wrpc_kernel[i, j] <- DensGenVonMises(b, mu1=mu1, mu2=mu2, kappa1=kappa1, 
          kappa2=kappa2)
        gpd_kernel[i, j] <- dgpd(r, sigma=scale, xi=shape, log=FALSE)
      }
    }
  } 
  wrpc_kernel <- apply(wrpc_kernel, 2, rev)
  gpd_kernel[center, center] <- 1/scale
  # This last part deletes the cells at the edge if they are all zero
  if (all(wrpc_kernel[1, ] == 0, wrpc_kernel[, 1] == 0, 
    wrpc_kernel[nrow(wrpc_kernel),] == 0, wrpc_kernel[, ncol(wrpc_kernel)] ==0))
    wrpc_kernel <- wrpc_kernel[2:(nrow(wrpc_kernel) - 1), 2:(ncol(wrpc_kernel) 
      - 1)]
  if (all(gpd_kernel[1, ] == 0, gpd_kernel[, 1] == 0, 
    gpd_kernel[nrow(gpd_kernel),] == 0, gpd_kernel[, ncol(gpd_kernel)] == 0))
    gpd_kernel <- gpd_kernel[2:(nrow(gpd_kernel) - 1), 2:(ncol(gpd_kernel) - 1)]
  # Multiply the two kernels together and re-normalize
  if (ignore_vonmises) wrpc_kernel <- 1
  if (ignore_pareto) gpd_kernel <- 1
  redist_kernel <- gpd_kernel*wrpc_kernel
  redist_kernel <- redist_kernel/sum(redist_kernel)
  return(redist_kernel)
}

# ExtractHomeConsDF Function ---------------------------------------------------

###  Extracts matrix of longitude and longitude values (as x,y) for home site or 
###    conspecifics within a distance buffer around the home site
###  Usage: ExtractHomeConsDF(df, home_id, id_var, buffer, proj, x, y, output)
###  Arguments: df = dataframe of home and conspecific locations
###             home_id = id of home site 
###             id_var = column name of id variable, default is "ExtJoinVal"
###             buffer = buffer distance around home 
###             proj = projection of 
###             x = column name of longitude
###             y = column name of latitude
###             output = "home" or "cons", 
###  Returns: matrix of latitude and longitude (i.e., (x,y)) values for either 
###    the home site or conspecific locations
###  Notes:
###  Blake Massey
###  2015.02.02

ExtractHomeConsDF <- function(df,
                              home_id,
                              id_var = "ExtJoinVal", 
                              buffer,
                              proj, 
                              x,
                              y,
                              output = "home") {
  suppressPackageStartupMessages(require(raster))
  suppressPackageStartupMessages(require(rgeos))
  suppressPackageStartupMessages(require(sp))
  df <- df
  df$id <- df[ ,id_var]
  df$x <- as.numeric(df[ ,x])
  df$y <- as.numeric(df[ ,y])
  home_xy <- subset(df, id == home_id, select = c("x", "y", "id"))
  if (output == "home") output_df <- home_xy
  if (output == "cons") { 
    
    coordinates(home_xy) <- ~x+y
    
    home_sp <- SpatialPoints(home_xy, crs(proj))
    
    
    home_sp_buff <- gBuffer(home_sp, width=buffer)
    cons_xy <- subset(df, id != home_id, select = c("x", "y", "id"))
    cons_xy <- coordinates(~x+y)
    cons_sp <- SpatialPoints(cons_xy, crs(proj))
    cons_sp_crop <- crop(cons_sp, extent(home_sp_buff))
    cons <- as.data.frame(cons_sp_crop)
    output_df <- cons
  }
  return(output_df)
}

# FitTurnAngleParsToDataLogistic Function --------------------------------------

###  Fits turn angle distribution data to turn angle function using logisitic
###  Usage: FitTurnAngleParsToDataLogisitic (df, by, nest_dist, con_dist, 
###    abs_angle, nest_angle, low_rho, low_location_nest, low_scale_nest,
###    low_location_con, low_scale_con )
###  Arguments: df = dataframe 
###             by = column to subset data 
###             var = variable to fit Wrapped Normal distribution
###             mu = starting value of mu, default is 1
###             rho = starting value of rho, default is 0
###  Returns: dataframe of "by" variable and parameters  
###  Notes:
###  Blake Massey
###  2014.12.05

FitTurnAngleParsToDataLogistic <- function(df,
                                           by = NULL,
                                           nest_dist = "nest_dist",
                                           con_dist =  "con_dist",
                                           abs_angle = "abs_angle",
                                           nest_angle = "nest_angle") {  
  suppressPackageStartupMessages(require(bbmle))
  df <- df
  ifelse(is.null(by), df$by <- "all", df$by <- df[,by])
  ifelse(is.null(by), by <- "by", by <- by)
  vars <- unique(df$by)
  pars <- data.frame(by=vars, rho=NA, location_nest=NA, scale_nest=NA, 
    location_con=NA, scale_con=NA, stringsAsFactors=FALSE)  
  for (i in vars){
    df_i <- subset(df, by == i)
    df_i <- df_i[!is.null(df_i), ]
    df_i$prev_abs_angle <- c(NA, df_i[1:nrow(df_i)-1, abs_angle])    
    df_i <- df_i[-1]
    df_i <- df_i[complete.cases(df_i[,c(nest_dist, con_dist, abs_angle, 
      nest_angle)]),]  
    start <- list(rho = 0, 
      location_nest = 0, 
      scale_nest = 0, 
      location_con = 0, 
      scale_con = 0)
    fit_i <- mle2(NLLTurnAngleLogistic, start=start,
      data=list(nest_dist=df_i[,nest_dist], con_dist=df_i[,con_dist],
      abs_angle=df_i[,abs_angle], nest_angle=df_i[,nest_angle], 
      prev_abs_angle=df_i[,"prev_abs_angle"]), method="Nelder-Mead") 
    pars[which(pars$by==i), "rho"] = as.numeric(coef(fit_i)["rho"])
    pars[which(pars$by==i), "location_nest"] = as.numeric(coef(fit_i)
      ["location_nest"])
    pars[which(pars$by==i), "scale_nest"] = as.numeric(coef(fit_i)
      ["scale_nest"])
    pars[which(pars$by==i), "location_con"] = as.numeric(coef(fit_i)
      ["location_con"])
    pars[which(pars$by==i), "scale_con"] = as.numeric(coef(fit_i)
      ["scale_con"])
  }
  names(pars)[names(pars) == 'by'] <- by
  return(pars)
}

# FitTurnAngleParsToDataTanH Function ---------------------------------------

###  Fits turn angle distribution data to turn angle function using TanH
###  Usage: FitTurnAngleParsToDataTanH(df, by, nest_dist, con_dist, abs_angle,
###    nest_angle, low_rho, low_b_nest, low_b_con, low_c_nest, low_c_con, 
###    upper_rho, upper_b_nest, upper_b_con, upper_c_nest, upper_c_con)
###  Arguments: df = dataframe 
###             by = column to subset data 
###             nest_dist = distance to nest
###             con_dist = distance to conspecific
###             abs_angle = absolute angle
###             nest_angle = angle to nest
###             low_rho = lower parameter limit for rho
###             low_b_nest = lower parameter limit for b_nest
###             low_b_con = lower parameter limit for b_con
###             low_c_nest = lower parameter limit for c_nest
###             low_c_con = lower parameter limit for c_con
###             upper_rho = upper parameter limit for rho
###             upper_b_nest = upper parameter limit for b_nest
###             upper_b_con = upper parameter limit for b_con
###             upper_c_nest = upper parameter limit for c_nest
###             upper_c_con = upper parameter limit for c_con
###  Returns: dataframe of "by" variable and parameters  
###  Notes:
###  Blake Massey
###  2014.12.05

FitTurnAngleParsToDataTanH <- function(df,
                                       by = NULL,
                                       nest_dist = "nest_dist",
                                       con_dist =  "con_dist",
                                       abs_angle = "abs_angle",
                                       nest_angle = "nest_angle",
                                       low_rho = 0,
                                       low_b_nest = .05,
                                       low_b_con = .5,
                                       low_c_nest = 0,
                                       low_c_con = -1,
                                       upper_rho = 1,
                                       upper_b_nest = .5,  
                                       upper_b_con = 5,
                                       upper_c_nest = 1,
                                       upper_c_con = 0) {  
  suppressPackageStartupMessages(require(bbmle))
  df <- df
  ifelse(is.null(by), df$by<-"all", df$by <- df[,by])
  ifelse(is.null(by), by <- "by", by <- by)
  vars <- unique(df$by)
  pars <- data.frame(by=vars, rho=NA, b_nest=NA, b_con=NA, c_nest=NA, c_con=NA, 
    stringsAsFactors=FALSE)  
  for (i in vars){
    df_i <- subset(df, by == i)
    df_i <- df_i[!is.null(df_i), ]
    df_i$prev_abs_angle <- c(NA, df_i[1:nrow(df_i)-1, abs_angle])    
    df_i <- df_i[-1]
    df_i <- df_i[complete.cases(df_i[,c(nest_dist, con_dist, abs_angle, 
      nest_angle)]),]  
    low_pars <- c(low_rho, low_b_nest, low_b_con, low_c_nest, low_c_con)
    upper_pars <- c(upper_rho, upper_b_nest, upper_b_con, upper_c_nest, 
      upper_c_con)
    start <- list(rho = (low_rho+upper_rho)/2, 
      b_nest = (low_b_nest+upper_b_nest)/2, 
      c_nest = (low_c_nest+upper_c_nest)/2, 
      b_con = (low_b_con+upper_b_con)/2, 
      c_con = (low_c_con+upper_c_con)/2)
    fit_i <- mle2(NLLTurnAngleTanH, start=start,
      data=list(nest_dist=df_i[,nest_dist], con_dist=df_i[,con_dist],
      abs_angle=df_i[,abs_angle], nest_angle=df_i[,nest_angle], 
      prev_abs_angle=df_i[,"prev_abs_angle"]), method="L-BFGS-B", 
      lower = low_pars, upper = upper_pars)
    pars[which(pars$by==i), "rho"] = as.numeric(coef(fit_i)["rho"])
    pars[which(pars$by==i), "b_nest"] = as.numeric(coef(fit_i)["b_nest"])
    pars[which(pars$by==i), "c_nest"] = as.numeric(coef(fit_i)["c_nest"])
    pars[which(pars$by==i), "b_con"] = as.numeric(coef(fit_i)["b_con"])
    pars[which(pars$by==i), "c_con"] = as.numeric(coef(fit_i)["c_con"])
  }
  names(pars)[names(pars) == 'by'] <- by
  return(pars)
}

# NLLTurnAngleLogistic Function ------------------------------------------------

###  Negative log-likelihood function for turn angle using logistic function
###  Usage: NLLTurnAngleLogistic(nest_dist, con_dist, abs_angle, nest_angle,
###    prev_abs_angle, rho, location_nest, scale_nest, location_con, scale_con) 
###  Arguments: nest_dist = distance to nest
###             con_dist = distance to conspecific nest
###             abs_angle = absolute angle
###             nest_angle = absolute angle to nest
###             prev_abs_angle = previous abs angle
###             rho = rho
###             location_nest = location parameter for nest
###             scale_con = scale parameter for conspecific
###             location_nest = location parameter for nest
###             scale_con = scale parameter for conspecific
###  Returns: NLL 
###  Notes: 
###  Blake Massey
###  2015.01.02

NLLTurnAngleLogistic <- function(nest_dist,
                                 con_dist,
                                 abs_angle,
                                 nest_angle, 
                                 prev_abs_angle,
                                 rho,
                                 location_nest,
                                 scale_nest,
                                 location_con,
                                 scale_con){  
  beta_nest = exp(location_nest + (scale_nest*nest_dist))/
    (1 + exp(location_nest + (scale_nest*nest_dist)))
  beta_con = exp(location_con + (scale_con*con_dist))/
    (1 + exp(location_con + (scale_con*df$con_dist)))
  beta = max(beta_nest, beta_con)
  exp_angle = ((1-beta)*(prev_abs_angle)) + (beta*nest_angle)
  -sum(log(dwrappedcauchy(abs_angle, mu=exp_angle, rho=rho)))
}

# NLLTurnAngleTanH Function ----------------------------------------------------

###  Negative log-Likelihood function for turn angle using TanH function
###  Usage: NLLTurnAngleTanH(data, nest_dist, con_dist, abs_angle, nest_angle,
###    prev_abs_angle, rho, b_nest, b_con, c_nest, c_con) 
###  Arguments: nest_dist = distance to nest
###             con_dist = distance to conspecific nest
###             abs_angle = absolute angle
###             nest_angle = absolute angle to nest
###             prev_abs_angle = previous abs angle
###             rho = rho
###             b_nest = b parameter for nest
###             b_con = b parameter for conspecific
###             c_nest = c parameter for nest
###             c_con = c parameter for conspecific
###  Returns: NLL of for turn angle
###  Notes: 
###  Blake Massey
###  2014.11.29

NLLTurnAngleTanH <- function(nest_dist,
                             con_dist,
                             abs_angle,
                             nest_angle, 
                             prev_abs_angle,
                             rho,
                             b_nest,
                             c_nest,
                             b_con,                         
                             c_con){  
  beta_nest = tanh(b_nest*nest_dist^c_nest)
  beta_con = tanh(b_con*con_dist^c_con)
  beta = max(beta_nest, beta_con)
  exp_angle = ((1-beta)*(prev_abs_angle)) + (beta*nest_angle)
  -sum(log(dwrappedcauchy(abs_angle, mu=exp_angle, rho=rho)))
}

# Plot3DRaster Function --------------------------------------------------------

###  Wrapper function for hist3D() and plotrgl() that plots an rgl plot
###    from a Raster layer
###  Usage: Plot3DRaster(raster, azimuth, colaltitude, col, border, x_lab, 
###    y_lab, z_lab, z_lim, main, legend_lab, rgl, rgl_window, spin, movie, 
###    movie_name, ...)
###  Arguments: raster = Raster layer to plot
###             azimuth = azimuth angle, default is 45
###             colalitude  = coaltitude angle, default is 30
###             col = color palette
###             border = color of the lines drawn around surface facets, default 
###               is "black".
###             x_lab = label for x-axis, default is "Longitude"
###             y_lab = label for y-axis, default is "Latitude"
###             z_lab = label for z-axis, default is ""
###             z_lim = range of values for z axis, default is range(z)
###             main = plot title, default is raster object name 
###             legend_lab = label for legend, default is z_lab
###             rgl = logical, create an interactive rgl object, default is TRUE
###             rgl_window = sets rgl window size, either "screen" or "image", 
###               optimized for viewing on screen or as an 1024x768 pixel image, 
###               respectively. Default is "screen".
###             spin = logical, spin the plot 360 degrees once, default is TRUE
###             movie = logical, create a .gif movie from the rgl plot. Default 
###               is FALSE.
###             movie_name = name of output gif movie. Default is saved in  
###               in working directroy as "RasterSpin.gif"
###             ... = additional arguments for the hist3D() function
###  Returns: Plot of Raster layer in RStudio, 3d plot in interactive rgl device
###    (optional), and a .gif movie of the plot rotating 360 degrees (optional)  
###  Notes: For additional arguments see ?persp3D. If a movie is made, a new 
###    rgl window will open set with the proper dimensions, record the movie, 
###    then automatically close. All NA values are converted to 0 because the 
###    hist3d() plots were not working with NA values included in the matrix.
###  Blake Massey
###  2015.01.18

Plot3DRaster <- function(raster,
                         azimuth = 45,
                         coaltitude = 30, 
                         col = NULL,
                         border="black",
                         x_lab = NULL,
                         y_lab = NULL,
                         z_lab = NULL,
                         z_lim = NULL,
                         main = NULL,
                         legend_lab = NULL,
                         rgl = TRUE, 
                         rgl_window = "screen",
                         spin = FALSE,
                         movie = FALSE,
                         movie_name = "RasterSpin",
                         ...) {
  suppressPackageStartupMessages(library(CircStats)) 
  suppressPackageStartupMessages(library(plot3D)) 
  suppressPackageStartupMessages(library(plot3Drgl)) 
  x <- xFromCol(raster, col=1:ncol(raster))
  y <- yFromRow(raster, row=1:nrow(raster))
  z <- t(as.matrix(raster))
  z[is.na(z)] <- 0  # otherwise the hist3d() plot does not work properly
  if (is.null(col)) col <- gg.col(length(unique(z)))
  if (is.null(x_lab)) x_lab <- "Longitude"
  if (is.null(y_lab)) y_lab <- "Latitude"
  if (is.null(z_lab)) z_lab <- ""
  if (is.null(z_lim)) z_lim <- range(z, na.rm = TRUE)
  if (is.null(main)) main <- deparse(substitute(raster))
  if (is.null(legend_lab)) legend_lab <- z_lab
  hist3D(x=x, y=y, z=z, shade=0, nticks=5, ticktype="detailed", col=col,
    bty="b2", expand=.25, phi=coaltitude, theta=azimuth, border=border, 
    facets=TRUE, axes=TRUE, image=FALSE, contour=FALSE, panel.first=NULL,
    ltheta=-135, lphi=0, space=0, add=FALSE, plot=TRUE, clab=legend_lab, 
    main=main, xlab="Longitude", ylab="Latitude", zlab="", zlim=z_lim, 
    colkey=list(side=4, line.clab=1,  length=.5, width=.5, adj.clab=0.1, 
      dist=-.03) , ...)
  ResetGraphics <- function(){
    rgl.clear(type = "bboxdeco")
    text_ids <- subset(rgl.ids(), type=="text", select="id") 
    for (i in 1:nrow(text_ids)){
      rgl.pop(id=text_ids[i,"id"])
    }   
    par(mar=c(2, 2, 2, 2)+.01, las=2)
    axis3d('x--', ntick=7)  # can be adjusted to add more or fewer tick marks
    axis3d('y+-', ntick=7)  # can be adjusted to add more or fewer tick marks
    axis3d('z--', ntick=4)  # can be adjusted to add more or fewer tick marks
    mtext3d(x_lab, edge='x--', line=2)
    mtext3d(y_lab, edge='y+-', line=2)
    mtext3d(z_lab, edge='z--', line=2.5)
    r1 <- rotationMatrix((coaltitude+270)*(pi/180), 1,0,0)  # 
    r2 <- rotationMatrix(-azimuth*pi/180, 0,0,1)  #
    r <- r1 %*% r2
    rgl.viewpoint(interactive=TRUE, userMatrix=r) # rotate
    observer3d(-0.075, -0.15, 3) 
    Sys.sleep(.5)
    bgplot3d({ 
      par(omd=c(0.75, 1.0000000, 0, 0.3), ps=18)
  #    par(cra=2)
      colkey(side = 4, clim = c(0, max(z)), add = FALSE, cex.clab = 1,  
        line.clab=.75, width = 2, length = 1.5, clab = legend_lab, 
        col=col, adj.clab = 0.05, cex.axis = .85)
      par(omd = c(0, 1, 0, .975), ps=35)
      title(main=main, font.main=2, cex.main=1)
    })
  }
  
  if (rgl == TRUE) { 
    plotrgl(new = TRUE) # new window
    if (rgl_window == "image") par3d(windowRect=c(150, 22, 1174, 790))  #      
    if (rgl_window == "screen") par3d(windowRect=c(0, 22, 1366, 728))  #
    ResetGraphics()
    if (spin == TRUE) { 
      play3d(spin3d(axis=c(0,0,1), rpm=6), duration=10)
    }
    if (movie == TRUE) { 
      if (rgl_window == "image") {     
        cat("Creating a movie file, this will take a few seconds", "\n") 
        movie3d(spin3d(axis=c(0,0,1), rpm=6), fps = 32, duration=10, 
        movie=movie_name, dir=getwd(), clean=TRUE) 
        cat(paste0("Created movie file: ", movie_name, ".gif"), "\n")
      }
      if (rgl_window == "screen") { 
        org <- as.numeric(rgl.cur())
        cat("Opening new rgl device with proper dimensions for a movie.", "\n")
        plotrgl(new = TRUE)
        par3d(windowRect=c(150, 22, 1174, 790))  # dimensions: 1028 X 768
        ResetGraphics()
        cat("Creating a movie file - this will take a few seconds.") 
        movie3d(spin3d(axis=c(0,0,1), rpm=6), fps = 32, duration=10,
          movie=movie_name, dir=getwd(), clean=TRUE, verbose=FALSE)
        cat(paste0("Created movie file: ", movie_name, ".gif"), "\n")
        cat("Returning to previous rgl device.")       
        rgl.close()
        rgl.set(which=org)
      }
    }
  }
}

# PlotAngleBiasTanHPDF Function ------------------------------------------------

###  Plot Angle Bias Probability Distribution Function using TanH
###  Usage: PlotAngleBiasPDF(b, c, x_max)
###  Arguments: b = b parameter
###             c = c parameter
###             x_max = maximum value on x scale, default is 500
###  Returns: plot of probability distribution
###  Notes:
###  Blake Massey
###  2014.12.15

PlotAngleBiasTanHPDF <- function(b = 1, 
                                 c = 1, 
                                 x_max = 500){
  suppressPackageStartupMessages(require(ggplot2))
  x <- seq(0, x_max, length=x_max)
  hx <- tanh(b*x^c)
  df <- as.data.frame(cbind(x, hx))
  main <- paste0("Angle Bias Distribution (b = ", signif(b,4), ", c = ", 
    signif(c,4), ")")
  g <- ggplot(df, aes(x, hx)) +
    geom_line(colour="dark green", size=1.5) + 
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black"))
  g + labs(x='Value', y='Probability Density', title=main)
}

# PlotSimHomeRange (DEPRECATED) Function ---------------------------------------

###  Plot locations, pathways, and bias values of SimulateHomeRangeLogistic(), 
###    and SimulateHomeRangeTanH() output dataframes
###  Usage: PlotSimHomeRange(df, label, x_lim, y_lim)
###  Arguments: df = dataframe with x, y, nest_x, and nest_y columns
###             x_lim = maximum value of x axis, default is NULL
###             y_lim = maximum value on y axis, default is NULL
###  Returns: plot of locations and pathways
###  Notes:
###  Blake Massey
###  2015.01.02

# PlotSimHomeRange <- function (df,
#                               x_lim = NULL,
#                               y_lim = NULL){
#   by_colors = c("nest"="green2", "con"="purple") 
#   g <- ggplot(df) + xlab("X") + ylab("Y") +  
#   geom_text(aes(x=x,y=y,label=signif(beta,3)), vjust=1.5, hjust=-.5, size=3) +
#   geom_point(data=df[1,], aes(x=nest_x, y=nest_y), size=8, color="black", 
#     fill="yellow",shape=22) + 
#   geom_point(data=df[1,], aes(x=con_x, y=con_y), size=8, color="black", 
#     fill="purple",shape=22) + 
#   geom_point(data=df[nrow(df),], aes(x=x, y=y),size=8, color="black", 
#     fill="red", shape=23) +     
#   geom_path(data=df, aes(x=x,y=y), arrow = arrow()) +
#   geom_point(data=df, aes(x=x,y=y, fill=factor(beta_name)), size=4, 
#     color="black", shape=21) +  
#   scale_fill_manual(values = by_colors) +
#   coord_fixed(ratio = 1) +
#   theme(legend.position="none") + theme(text=element_text(size=20, 
#     colour="black")) + theme(axis.text=element_text(colour="black")) + 
#   theme(axis.title.x = element_text(angle = 0, vjust = 0, hjust=0.5)) +
#   theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) 
#   if (!is.null(x_lim)) {
#     g <- g + scale_x_continuous(limits=c(0, x_lim))
#   }
#   if (!is.null(y_lim)) { 
#     g <- g + scale_y_continuous(limits=c(0, y_lim))
#   }
#   return(g)
# }

# PlotSimHomeRange Function ----------------------------------------------------

###  Plot locations, pathways, and bias values of SimulateHomeRange() 
###    output dataframes
###  Usage: PlotSimHomeRange(sim, home_df, homerange, col, label, x_range, 
###    y_range)
###  Arguments: sim = dataframe with x, y columns
###             home_df = 
###             homerange = 
###             col = color palette for base raster, default is jet2.col()
###             label = column name for point labels
###             x_range = x-axis range, default is NULL
###             y_range = y-axis range, default is NULL
###  Returns: plot of locations and pathways
###  Notes:
###  Blake Massey
###  2015.02.09

PlotSimHomeRange <- function (sim,
                              df_home,
                              homerange,
                              col = NULL,
                              label = NULL,
                              x_range = NULL,
                              y_range = NULL){
  require(grid)
  sim <- sim
  if (!is.null(label)) sim$label <- sim[,label]
  df_home <- df_home
  library(ggplot2)
  if (!is.null(homerange)) {
    homerange_df <- data.frame(rasterToPoints(homerange))
    colnames(homerange_df)[3] <- "Probability"
    if(is.null(col)) col <- jet2.col(length(unique(homerange)))
    g <- ggplot(homerange_df, aes(x=x, y=y))+ geom_raster(aes(fill=Probability)                                          , 
      interpolate=TRUE) + scale_fill_gradientn(colours=col)
  } else {
    g <- ggplot(data=sim) + theme(legend.position="none")
  } 
  g <- g + xlab("X") + ylab("Y") +   
  geom_point(data=sim[1,], aes(x=x, y=y),size=3, color="green", 
    fill="white", shape=18) +
  geom_point(data=sim[nrow(sim), ], aes(x=x, y=y),size=3, color="red", 
    fill="black", shape=18) +
  geom_path(data=sim, aes(x=x, y=y), color= "gray90", arrow = arrow(angle=20, 
    length=unit(.01, "npc"))) + 
  geom_point(data=sim, aes(x=x, y=y), size=2, color="gray20", shape=20) +   
  coord_fixed(ratio = 1) +
  theme(text=element_text(size=20, colour="black")) + 
  theme(axis.text=element_text(colour="black")) + 
  theme(axis.title.x = element_text(angle = 0, vjust = 0, hjust=0.5)) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=0.5))  
#  theme(panel.background = element_rect(fill = "white"))
  if (!is.null(label)) {   
    g <- g + geom_text(aes(x=x, y=y, label=signif(label,3)), vjust=1.5, 
      hjust=-.5, size=3)
  }
  if (!is.null(x_range)) {
    g <- g + scale_x_continuous(limits=c(x_range))
  }
  if (!is.null(y_range)) { 
    g <- g + scale_y_continuous(limits=c(y_range))
  }
  return(g)
}

# SimulateHomeRangeLogistic Function -------------------------------------------

###  Runs a simulation of home range behavior based on nest and conspecific 
###    locations and a logistic function for angle bias
###  Usage: SimulateHomeRangeLogisitic(n, nest_x, nest_y, con_x, con_y, mu, rho, 
###    pareto_scale, pareto_shape, location_nest, scale_nest, location_con, 
###    scale_con, nest_prob)
###  Arguments: n = number of locations
###             nest_x = longitude of nest
###             nest_y = latitude of nest
###             con_x = longitude of conspecific
###             con_y = latitude of conspecific
###             mu = wrapped cauchy mu parameter
###             rho = wrapped cauch rho parameter
###             pareto_scale = pareto scale parameter
###             pareto_shape = pareto shape parameter
###             location_nest = location parameter for nest 
###             scale_nest = scale c parameter for nest
###             location_con = location b parameter for con
###             scale_con = scale c parameter for con
###             nest_prob = proability to return to nest at any time step
###  Returns: dataframe with simulation data  
###  Notes: still a work in progress
###  Blake Massey
###  2015.01.02

SimulateHomeRangeLogistic <- function(n,
                                      nest_x, 
                                      nest_y, 
                                      con_x,
                                      con_y,
                                      mu = 0,  # Warpped Cauchy parameter 
                                      rho = .2,  # Wrapped Cauchy parameter
                                      pareto_scale = 10,
                                      pareto_shape = 5,
                                      location_nest = 0, #.21  
                                      scale_nest = 0, #.5
                                      location_con = 0, #.001
                                      scale_con = 0,
                                      nest_chance = 0) { 
  library(circular)
  library(grid)
  library(ggplot2)
  library(VGAM)
  source('C:/Work/R/Functions/gps.R')
  n <- n
  df <- data.frame(n = seq(1:n))
  df$x <- 0
  df$y <- 0
  df$abs_angle <- 0
  df$exp_angle <- 0 
  df$beta <- 0
  df$turn_angle <- 0
  df$diff_angle <- 0
  df$nest_angle <- 0
  df$con_angle <- 0
  df$step_length <- 0
  df$nest_dist <- 0
  df$con_dist <- 0
  df$nest_x <- nest_x
  df$nest_y <- nest_y
  df$con_x <- con_x
  df$con_y <- con_y
  ConvertAngle <- function(x) {
    if (x > 2*pi) x <- x-(2*pi)
    if (x < 0) x <- x+(2*pi)
    return(x)
  }
  for (i in 1:n){
    if (i == 1){
      df$x[i] <- df$nest_x[i]
      df$y[i] <- df$nest_y[i]
      df$abs_angle[i] <- sample(x=seq(from=0, to=(2*pi), by=(2*pi/360)), size=1)
      df$exp_angle[i] <- NA
      df$beta <- NA
      df$beta_name <- NA
      df$corr_turn_angle <- NA
      df$nest_turn_angle <- NA
      df$con_turn_angle <- NA
      df$turn_angle[i] <- NA
      df$diff_angle[i] <- NA
      df$nest_angle[i] <- NA 
      df$con_angle[i] <- NA
      df$nest_dist[i] <- as.integer(sqrt(sum((c(df$x[i],df$y[i]) - 
        c(df$nest_x[i], df$nest_y[i]))^2)))
      df$con_dist[i] <- as.integer(sqrt(sum((c(df$x[i],df$y[i]) - c(df$con_x[i], 
        df$con_y[i]))^2)))
      df$step_length[i] <- rpareto(1, pareto_scale, pareto_shape)
      df$x[i+1] <- df$x[i] + df$step_length[i]*(cos(df$abs_angle[i]))
      df$y[i+1] <- df$y[i] + df$step_length[i]*(sin(df$abs_angle[i]))
      df$diff_angle <- NA
    }
    if (i > 1 && i < n){
      nest <- rbinom(1, 1, nest_chance) 
      if (nest == TRUE) {
        df$nest_dist[i] <- as.integer(sqrt(sum((c(df$x[i],df$y[i]) - 
          c(df$nest_x[i], df$nest_y[i]))^2)))
        df$con_dist[i] <- as.integer(sqrt(sum((c(df$x[i],df$y[i]) - 
          c(df$con_x[i], df$con_y[i]))^2))) 
        df$nest_angle[i] <- CalculateAngleToPoint(df[i,], "x", "y", "nest_x",
          "nest_y")
        df$step_length[i] <- df$nest_dist[i]  # forces distance to nest
        beta_nest = exp(location_nest + scale_nest*df$nest_dist[i])/
          (1 + exp(location_nest + scale_nest*df$nest_dist[i]))
        beta_con = exp(location_con + scale_con*df$con_dist[i])/
          (1 + exp(location_con + scale_con*df$con_dist[i]))
        df$beta[i] = max(beta_nest, beta_con)      
        df$beta_name[i] <- names(which.max(c(nest=beta_nest, con=beta_con)))
        df$corr_turn_angle[i] = suppressWarnings(rwrappedcauchy(1, mu=mu,r=rho))
        df$corr_turn_angle[i] <- ConvertAngle(df$corr_turn_angle[i])
        df$nest_turn_angle[i] <-  df$nest_angle[i] - df$abs_angle[i-1] 
        df$nest_turn_angle[i] <- ConvertAngle(df$nest_turn_angle[i])
        df$turn_angle[i] <- df$nest_turn_angle[i]  # forces turn to nest
        df$turn_angle[i] <- ConvertAngle(df$turn_angle[i])
        df$abs_angle[i] <- df$abs_angle[i-1] + df$turn_angle[i]
        df$abs_angle[i] <- ConvertAngle(df$abs_angle[i])
        df$x[i+1] <- df$x[i] + df$step_length[i]*(cos(df$abs_angle[i]))
        df$y[i+1] <- df$y[i] + df$step_length[i]*(sin(df$abs_angle[i]))
        df$turn_angle[i] <- df$abs_angle[i-1] - df$abs_angle[i]
        df$turn_angle[i] <- ConvertAngle(df$turn_angle[i])
        df$diff_angle[i] <- df$abs_angle[i] - df$nest_angle[i]
        df$diff_angle[i] <- ConvertAngle(df$diff_angle[i])
      } else {
        df$nest_dist[i] <- as.integer(sqrt(sum((c(df$x[i],df$y[i]) - 
          c(df$nest_x[i], df$nest_y[i]))^2)))
        df$con_dist[i] <- as.integer(sqrt(sum((c(df$x[i],df$y[i]) - 
          c(df$con_x[i], df$con_y[i]))^2))) 
        df$nest_angle[i] <- CalculateAngleToPoint(df[i,], "x", "y", "nest_x",
          "nest_y")
        df$step_length[i] <- rpareto(1, pareto_scale, pareto_shape)  
        
        beta_nest = exp(location_nest + (scale_nest*df$nest_dist[i]))/
          (1 + exp(location_nest + (scale_nest*df$nest_dist[i])))
        beta_con = exp(location_con + (scale_con*df$con_dist[i]))/
          (1 + exp(location_con + (scale_con*df$con_dist[i])))
        df$beta[i] = max(beta_nest, beta_con)   
        df$corr_turn_angle[i] <- suppressWarnings(rwrappedcauchy(1,mu=mu,r=rho))
        df$corr_turn_angle[i] <- ConvertAngle(df$corr_turn_angle[i])
        df$beta_name[i] <- names(which.max(c(nest = beta_nest, con=beta_con)))
        df$nest_turn_angle[i] <-  df$nest_angle[i] - df$abs_angle[i-1] 
        df$nest_turn_angle[i] <- ConvertAngle(df$nest_turn_angle[i])
        df$turn_angle[i] = ((1-df$beta[i])*(df$corr_turn_angle[i])) + 
          (df$beta[i]*df$nest_turn_angle[i])
        df$turn_angle[i] <- ConvertAngle(df$turn_angle[i])
        df$abs_angle[i] <- df$abs_angle[i-1] + df$turn_angle[i]
        df$abs_angle[i] <- ConvertAngle(df$abs_angle[i])
        df$x[i+1] <- df$x[i] + df$step_length[i]*(cos(df$abs_angle[i]))
        df$y[i+1] <- df$y[i] + df$step_length[i]*(sin(df$abs_angle[i]))
        df$turn_angle[i] <- df$abs_angle[i-1] - df$abs_angle[i]
        df$turn_angle[i] <- ConvertAngle(df$turn_angle[i])
        df$diff_angle[i] <- df$abs_angle[i] - df$nest_angle[i]
        df$diff_angle[i] <- ConvertAngle(df$diff_angle[i])
      }
    }
    if (i == n){
      df$nest_dist[i] <- as.integer(sqrt(sum((c(df$x[i],df$y[i]) - 
        c(df$nest_x[i], df$nest_y[i]))^2)))
      df$con_dist[i] <- as.integer(sqrt(sum((c(df$x[i],df$y[i]) - c(df$con_x[i], 
        df$con_y[i]))^2))) 
      df$nest_angle[i] <- CalculateAngleToPoint(df[i,], "x", "y", "nest_x",
        "nest_y")
      df$con_angle[i] <- NA
      df$step_length[i] <- 0
      df$exp_angle[i] <- NA 
      df$abs_angle[i] <- NA
      df$turn_angle[i] <- NA
      df$diff_angle[i] <- NA
    }
  }
  return(df)
}

# SimulateHomeRange Function ---------------------------------------------------

###  Runs a simulation of home range behavior based on nest and conspecific 
###    locations and a logistic function for angle bias
###  Usage: SimulateHomeRange(n, home, cons, mu, rho, pareto_scale, 
###    pareto_shape, max_r, start_x, start_y, home_prob, save_redist_plots, 
###    save_arrow_plots, alpha)
###  Arguments: steps = number of locations
###             home = dataframe of home centroid (with "x", "y", and "id" 
###               columns)
###             homerange = homerange_kernel
###             mu = wrapped Cauchy mu parameter, default is 0
###             rho = wrapped Cauchy rho parameter
###             pareto_scale = Pareto scale parameter
###             pareto_shape = Pareto shape parameter
###             max_r = maximum radius of redistribution kernels. Default is 
###               NULL and automatically set to .995 quantile of the step length 
###               Pareto probability distribution 
###             start_xy = starting longitude and latitude location. Default is 
###               NULL and automatically set to the home location's x.
###             home_prob = proability to return to home at any time step. 
###               Default is 0.
###             save_redist_plots = save all of the redistribution kernel 
###               probability plots in the working directory. Default is FALSE.
###             save_arrow_plots = save all of the redistribution kernel 
###               probability plots with movement arrows in the working 
###               directory. Default is FALSE.
###             alpha = alpha (i.e., transparency) value for plotting the 
###               redistribution kernels. Default is .9 
###  Returns: Dataframe with simulation data
###  Notes: Still a work in progress. May modify code so all of the 
###    redistribution kernels are saved as a RasterStack. This would allow the 
###    legend scale to be set at fixed values and to potentially create 
###    ArcGIS rasters and KMLs after the simultion is run. 
###  Blake Massey
###  2015.01.08

SimulateHomeRange <- function(steps,
                              home, 
                              homerange, 
                              mu = 0, # Warpped Cauchy parameter 
                              rho,  # Wrapped Cauchy parameter
                              pareto_scale,
                              pareto_shape,
                              max_r = NULL,
                              start_xy = NULL,
                              home_prob = 0, 
                              x_range,
                              y_range,
                              save_redist_plots = FALSE,
                              save_arrow_plots = FALSE,
                              alpha = .9) {
  suppressPackageStartupMessages(require(circular))
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(grid))
  suppressPackageStartupMessages(require(rasterVis))
  suppressPackageStartupMessages(require(sampling))
  suppressPackageStartupMessages(require(VGAM))
  source('C:/Work/R/Functions/gen.R')
  source('C:/Work/R/Functions/gps.R')  
  source('C:/Work/R/Functions/home.R') 
  source('C:/Work/R/Functions/pars.R')  
  cellsize <- res(homerange)[1]
  home_xy <- as.numeric(home[,c("x", "y")])
  home_xy <- CenterXYInCell(home_xy[1], home_xy[2], xmin(base), ymin(base), 
   cellsize)
  n <- steps
  df <- data.frame(step = seq(1:n))
  df$x <- NA
  df$y <- NA
  df$abs_angle <- NA
  df$exp_angle <- NA
  if (is.null(max_r)) { 
    max_r <- cellsize*ceiling(qgpd(.995, pareto_scale, pareto_shape)/cellsize)
  }
  for (i in 1:n) {
    if (i == 1) {
      df$x[i] <- CenterXYInCell(start_xy[1], start_xy[2], xmin(homerange), 
        ymin(base), cellsize)[1] 
      df$y[i] <- CenterXYInCell(start_xy[1], start_xy[2], xmin(homerange), 
        ymin(base), cellsize)[2] 
      df$abs_angle[i] <- sample(x=seq(from=0, to=(2*pi), by=(2*pi/360)), size=1)     
      redist <- CreateRedistKernel(max_r=max_r, cellsize=cellsize, 
        mu=df$abs_angle[i], rho=rho, shape=pareto_shape, scale=pareto_scale)
      r <- (cellsize*((nrow(redist)-1)/2))+(cellsize/2)
      redist_raster <- raster(redist, xmn=-r, xmx=r, ymn=-r, ymx=r)
      redist_shift <- shift(redist_raster, x=df$x[i], y=df$y[i])
### ADDING IN OTHER PROBABILITY LAYERS
                
        homerange_crop <- crop(homerange, redist_shift)
        prob_raster <- overlay(redist_shift, homerange_crop, fun=function(x,y) 
          {return(x*y)}, recycle=FALSE) 
                
### END OF OTHER PROBABILITY LAYERS
      writeLines(noquote(paste0("Calculating prob raster for step ", i)))
      if (save_redist_plots == TRUE) {
        prob_raster_NA <- prob_raster
        prob_raster_NA[prob_raster_NA <= 0] <- NA
        g <- gplot(prob_raster_NA) + geom_raster(aes(fill=value), alpha=alpha) +
          coord_equal() +
          scale_x_continuous(limits=c(x_range[1], x_range[2])) +
          scale_y_continuous(limits=c(y_range[1], y_range[2])) +
          scale_fill_gradient(low="yellow", high="red", na.value = NA, 
            breaks=c(minValue(prob_raster_NA), maxValue(prob_raster_NA)), 
            labels=c("min", "max")) +
          labs(fill="p") +
          theme(plot.title=element_text(size=22)) +
          theme(text=element_text(size=20, colour="black")) +
          theme(axis.text=element_text(colour="black")) + 
          xlab("X") + ylab("Y") +
          labs(title =  paste0("Step ", i))
        SaveGGPlot(filename = paste0("Step", sprintf("%03d", i), ".jpeg"))
      }
      destination_cell <- suppressWarnings(strata(data=data.frame(cell=
        1:ncell(prob_raster)), stratanames=NULL, size=1, method="systematic", 
        pik=prob_raster@data@values))
      destination_xy <- as.vector(xyFromCell(prob_raster,destination_cell[1,1]))
      df$x[i+1] <- destination_xy[1]
      df$y[i+1] <- destination_xy[2]
      if (save_arrow_plots == TRUE) {
        df2 <- cbind.data.frame(x=df$x[i], y=df$y[i], xend=df$x[i+1], 
          yend=df$y[i+1])
        g2 <- g + geom_segment(data=df2, aes(x=x, y=y, xend=xend, yend=yend), 
          arrow = arrow(length = unit(0.01, "npc")))
        SaveGGPlot(filename = paste0("Step", sprintf("%03d", i), "a.jpeg"))
      }
    }
    if (i > 1 && i < n) {
      go_home <- rbinom(1, 1, home_prob) 
      if (go_home == TRUE) {
        df$x[i+1] <- home_xy[1]
        df$y[i+1] <- home_xy[2]
        df$abs_angle[i] <- CalculateAngleToPoint2(df$x[i], df$y[i], df$x[i+1],
          df$y[i+1])    
      } else {
        df$exp_angle[i] <- df$abs_angle[i-1]
        redist <- CreateRedistKernel(max_r=max_r, cellsize=cellsize, 
          mu=df$exp_angle[i], rho=rho, shape=pareto_shape, scale=pareto_scale)
        r <- (cellsize*((nrow(redist)-1)/2))+(cellsize/2)
        redist_raster <- raster(redist, xmn=-r, xmx=r, ymn=-r, ymx=r)
        redist_shift <- shift(redist_raster, x=df$x[i], y=df$y[i])
        
### ADDING IN OTHER PROBABILITY LAYERS
                
        homerange_crop <- crop(homerange, redist_shift, snap="out")
        prob_raster <- overlay(redist_shift, homerange_crop, fun=function(x,y) 
          {return(x*y)}, recycle=FALSE) 
                
### END OF OTHER PROBABILITY LAYERS

        writeLines(noquote(paste0("Calculating prob raster for step ", i)))
        if (save_redist_plots == TRUE) {
          prob_raster_NA <- prob_raster
          prob_raster_NA[prob_raster_NA <= 0] <- NA
          g <- gplot(prob_raster_NA) + geom_raster(aes(fill=value),alpha=alpha)+
            coord_equal() +
            scale_x_continuous(limits=c(x_range[1], x_range[2])) +
            scale_y_continuous(limits=c(y_range[1], y_range[2])) +
            scale_fill_gradient(low="yellow", high="red", na.value = NA, 
              breaks=c(minValue(prob_raster_NA), maxValue(prob_raster_NA)), 
              labels=c("min", "max")) +
            labs(fill = "p") +
            theme(plot.title=element_text(size=22)) +
            theme(text=element_text(size=20, colour="black")) +
            theme(axis.text=element_text(colour="black")) + 
            xlab("X") + ylab("Y") +
            labs(title =  paste0("Step ", i))
            SaveGGPlot(filename = paste0("Step", sprintf("%03d", i), ".jpeg"))
        }
        destination_cell <- suppressWarnings(strata(data=data.frame(cell=
          1:ncell(prob_raster)), stratanames=NULL, size=1, method="systematic", 
          pik=prob_raster@data@values))
        destination_xy <- xyFromCell(prob_raster, destination_cell[1,1])
        df$x[i+1] <- destination_xy[1]
        df$y[i+1] <- destination_xy[2]
        df$abs_angle[i] <- CalculateAngleToPoint2(df$x[i], df$y[i], df$x[i+1],
          df$y[i+1])
        if (save_arrow_plots == TRUE) {
          df2 <- cbind.data.frame(x=df$x[i], y=df$y[i], xend=df$x[i+1], 
            yend=df$y[i+1]) 
          g2 <- g + geom_segment(data=df2, aes(x=x, y=y, xend=xend, yend=yend),
            arrow = arrow(length = unit(0.01, "npc")))
          SaveGGPlot(filename = paste0("Step", sprintf("%03d", i), "a.jpeg"))
        }
      }
    }
    if (i == n) {
      df$exp_angle[i] <- NA 
      df$abs_angle[i] <- NA
    }
  } 
  return(df)
}

# SimulateHomeRangeTanH Function -----------------------------------------------

###  Runs a simulation of home range behavior based on nest and conspecific 
###    locations and a tanh function for angle bias
###  Usage: SimulateHomeRangeTanH(n, nest_x, nest_y, con_x, con_y, mu, rho, 
###    pareto_scale, pareto_shape, b_nest, c_nest, b_con, c_con, nest_prob)
###  Arguments: n = number of locations
###             nest_x = longitude of nest
###             nest_y = latitude of nest
###             con_x = longitude of conspecific
###             con_y = latitude of conspecific
###             mu = wrapped cauchy mu parameter
###             rho = wrapped cauch rho parameter
###             pareto_scale = pareto scale parameter
###             pareto_shape = pareto shape parameter
###             b_nest = tanh b parameter for nest 
###             c_nest = tanh c parameter for nest
###             b_con = tanh b parameter for con
###             c_con = tanh c parameter for con
###             nest_prob = proability to return to nest at any time step
###  Returns: dataframe with simulation data  
###  Notes: still a work in progress
###  Blake Massey
###  2015.01.02

SimulateHomeRangeTanH <- function(n,
                                  nest_x, 
                                  nest_y, 
                                  con_x,
                                  con_y,
                                  mu = 0,  # Cauchy 
                                  rho = .2,  # Cauchy
                                  pareto_scale = 10,
                                  pareto_shape = 5,
                                  b_nest = 0,
                                  c_nest = 0, 
                                  b_con = 0,
                                  c_con = 0,
                                  nest_prob = 0) { 
  library(circular)
  library(grid)
  library(ggplot2)
  library(VGAM)
  source('C:/Work/R/Functions/gps.R')
  n <- n
  df <- data.frame(n = seq(1:n))
  df$x <- 0
  df$y <- 0
  df$abs_angle <- 0
  df$exp_angle <- 0 
  df$beta <- 0
  df$turn_angle <- 0
  df$diff_angle <- 0
  df$nest_angle <- 0
  df$con_angle <- 0
  df$step_length <- 0
  df$nest_dist <- 0
  df$con_dist <- 0
  df$nest_x <- nest_x
  df$nest_y <- nest_y
  df$con_x <- con_x
  df$con_y <- con_y
  for (i in 1:n){
    if (i == 1){
      df$x[i] <- df$nest_x[i]
      df$y[i] <- df$nest_y[i]
      df$abs_angle[i] <- sample(x=seq(from=0, to=(2*pi), by=(2*pi/360)), size=1)
      df$exp_angle[i] <- NA
      df$beta <- NA
      df$beta_name <- NA
      df$corr_turn_angle <- NA
      df$nest_turn_angle <- NA
      df$con_turn_angle <- NA
      df$turn_angle[i] <- NA
      df$diff_angle[i] <- NA
      df$nest_angle[i] <- NA 
      df$con_angle[i] <- NA
      df$nest_dist[i] <- as.integer(sqrt(sum((c(df$x[i],df$y[i]) - 
        c(df$nest_x[i], df$nest_y[i]))^2)))
      df$con_dist[i] <- as.integer(sqrt(sum((c(df$x[i],df$y[i]) - c(df$con_x[i], 
        df$con_y[i]))^2)))
      df$step_length[i] <- rpareto(1, pareto_scale, pareto_shape)
      df$x[i+1] <- df$x[i] + df$step_length[i]*(cos(df$abs_angle[i]))
      df$y[i+1] <- df$y[i] + df$step_length[i]*(sin(df$abs_angle[i]))
      df$diff_angle <- NA
    }
    if (i > 1 && i < n){
      nest <- rbinom(1, 1, nest_prob) 
      if (nest == TRUE) {
        df$nest_dist[i] <- as.integer(sqrt(sum((c(df$x[i],df$y[i]) - 
          c(df$nest_x[i], df$nest_y[i]))^2)))
        df$con_dist[i] <- as.integer(sqrt(sum((c(df$x[i],df$y[i]) - 
          c(df$con_x[i], df$con_y[i]))^2))) 
        df$nest_angle[i] <- CalculateAngleToPoint(df[i,], "x", "y", "nest_x",
          "nest_y")
        df$step_length[i] <- df$nest_dist[i]  # forces distance to nest
        beta_nest <- tanh(b_nest*(df$nest_dist[i])^c_nest)
        beta_con <- tanh(b_con*(df$con_dist[i])^c_con)
        df$beta[i] <- max(beta_nest, beta_con)      
        df$corr_turn_angle[i] <- suppressWarnings(rwrappedcauchy(1,mu=mu,r=rho))
        df$corr_turn_angle[i] <- ConvertAngle(df$corr_turn_angle[i])
        df$nest_turn_angle[i] <- df$nest_angle[i] - df$abs_angle[i-1] 
        df$nest_turn_angle[i] <- ConvertAngle(df$nest_turn_angle[i])
        df$turn_angle[i] <- df$nest_turn_angle[i]  # forces turn to nest
        df$turn_angle[i] <- ConvertAngle(df$turn_angle[i])
        df$abs_angle[i] <- df$abs_angle[i-1] + df$turn_angle[i]
        df$abs_angle[i] <- ConvertAngle(df$abs_angle[i])
        df$x[i+1] <- df$x[i] + df$step_length[i]*(cos(df$abs_angle[i]))
        df$y[i+1] <- df$y[i] + df$step_length[i]*(sin(df$abs_angle[i]))
        df$turn_angle[i] <- df$abs_angle[i-1] - df$abs_angle[i]
        df$turn_angle[i] <- ConvertAngle(df$turn_angle[i])
        df$diff_angle[i] <- df$abs_angle[i] - df$nest_angle[i]
        df$diff_angle[i] <- ConvertAngle(df$diff_angle[i])
      } else {
        df$nest_dist[i] <- as.integer(sqrt(sum((c(df$x[i],df$y[i]) - 
          c(df$nest_x[i], df$nest_y[i]))^2)))
        df$con_dist[i] <- as.integer(sqrt(sum((c(df$x[i],df$y[i]) - 
          c(df$con_x[i], df$con_y[i]))^2))) 
        df$nest_angle[i] <- CalculateAngleToPoint(df[i,], "x", "y", "nest_x",
          "nest_y")
        df$step_length[i] <- rpareto(1, pareto_scale, pareto_shape)
        beta_nest = tanh(b_nest*(df$nest_dist[i])^c_nest)
        beta_con = tanh(b_con*(df$con_dist[i])^c_con)
        df$beta[i] = max(beta_nest, beta_con)   
        df$corr_turn_angle[i] <- suppressWarnings(rwrappedcauchy(1,mu=mu,r=rho))
        df$corr_turn_angle[i] <- ConvertAngle(df$corr_turn_angle[i])
        df$beta_name[i] <- names(which.max(c(nest = beta_nest, con=beta_con)))
        df$nest_turn_angle[i] <-  df$nest_angle[i] - df$abs_angle[i-1] 
        df$nest_turn_angle[i] <- ConvertAngle(df$nest_turn_angle[i])
        df$turn_angle[i] = ((1-df$beta[i])*(df$corr_turn_angle[i])) + 
          (df$beta[i]*df$nest_turn_angle[i])
        df$turn_angle[i] <- ConvertAngle(df$turn_angle[i])
        df$abs_angle[i] <- df$abs_angle[i-1] + df$turn_angle[i]
        df$abs_angle[i] <- ConvertAngle(df$abs_angle[i])
        df$x[i+1] <- df$x[i] + df$step_length[i]*(cos(df$abs_angle[i]))
        df$y[i+1] <- df$y[i] + df$step_length[i]*(sin(df$abs_angle[i]))
        df$turn_angle[i] <- df$abs_angle[i-1] - df$abs_angle[i]
        df$turn_angle[i] <- ConvertAngle(df$turn_angle[i])
        df$diff_angle[i] <- df$abs_angle[i] - df$nest_angle[i]
        df$diff_angle[i] <- ConvertAngle(df$diff_angle[i])
      }
    }
    if (i == n){
      df$nest_dist[i] <- as.integer(sqrt(sum((c(df$x[i],df$y[i]) - 
        c(df$nest_x[i], df$nest_y[i]))^2)))
      df$con_dist[i] <- as.integer(sqrt(sum((c(df$x[i],df$y[i]) - c(df$con_x[i], 
        df$con_y[i]))^2))) 
      df$nest_angle[i] <- CalculateAngleToPoint(df[i,], "x", "y", "nest_x",
        "nest_y")
      df$con_angle[i] <- NA
      df$step_length[i] <- 0
      df$exp_angle[i] <- NA 
      df$abs_angle[i] <- NA
      df$turn_angle[i] <- NA
      df$diff_angle[i] <- NA
    }
  }
  return(df)
}
