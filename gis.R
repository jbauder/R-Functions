# --- GIS FUNCTIONS ------------------------------------------------------------
# Functions for importing, manipulating, and displaying GIS data 
# ------------------------------------------------------------------------------

# AddLandscapeValues Function --------------------------------------------------

###  Adds the cell values from a list of landscape rasters to an input locations
###    dataset, usually a subset of BAEA.        
###  Usage: AddLandscapeValues(df, raster_stack, long, lat, clean_up)
###  Arguments: df = dataframe of location data
###             raster_stack = a RasterStack of landscape data
###             long = column name of longitude data, default = lat_utm
###             lat = column name of latitude data, default = long_utm
###             clean_up = logical, will run a series of steps to clean up the 
###               returned output (e.g. removing "_50mc" from column names and 
###               deleting unneeded landcover "definition" column).
###  Returns: The original locations df with additional landscape data columns.             
###  Notes: The locations dataset and rasters should have the same CRS, 
###    otherwise function stops and returns an error message. 
###  Blake Massey
###  2014.09.09

AddLandscapeValues <- function(df,
                               raster_stack,
                               long = "long_utm",
                               lat = "lat_utm",
                               clean_up = TRUE){
  suppressPackageStartupMessages(require(plyr))
  suppressPackageStartupMessages(require(raster))
  df <- df
  raster_stack <- raster_stack
  locs_xy <- cbind(df[,long], df[,lat])
  raster_crs <- sapply(raster_stack@layers, crs)  
  compare <- function(v) all(sapply(v[-1], FUN=function(z)compareCRS(z,v[[1]])))
  if (compare(raster_crs) == FALSE) {
     stop("Input rasters have different coordinate reference systems.")
  }
  col_names_list <- {}
  for (i in 1:nlayers(raster_stack)){
    raster_layer <- raster_stack[[i]]
    col_name <- raster_layer@data@names
    if (clean_up == TRUE){
      col_name <- gsub("_50mc", "", col_name)
    } 
    df[,col_name] <- extract(raster_layer, locs_xy)
    col_names_list <- append(col_names_list, col_name)
  }
  if (clean_up == TRUE){
    if ("elev" %in% colnames(df)){
      df$agl <- df$alt - df$elev
      col_names_list <- append(col_names_list, "agl")
    }
    if ("lc" %in% colnames(df)) {
      nlcd_classes <- read.csv(file.path("C:/ArcGIS/Data/Landcover",
        "NCLD_Landcover_Class_Definitions.csv"), header=TRUE, as.is=TRUE, 
        na.strings = "") 
   #   if(is.integer(nlcd_classes$lc)) nlcd_classes$lc <- as.numeric(nlcd_classes$lc)
      df <- plyr::join(x=df, y=nlcd_classes, by="lc")
      df$definition <- NULL
    }
  }
  cat("The following columns were added to the dataframe:",sep="\n")
  names_sorted<-rev(sort(col_names_list, decreasing=TRUE))
  results <- sapply(names_sorted, function(i) paste(" ", i))
  cat(results, sep="\n")
 return(df)
}

# CreateArcGISMaps Function ----------------------------------------------------

###  A wrapper function for creating an ArcGIS map from location data
###  Usage: CreateArcGISMaps(df)
###  Arguments: df = dataframe with locations
###  Returns: overwrites files in "C:/Work/Python" folders          
###  Notes: internal parameters of Python script specific to my BAEA GPS data
###    Identical to CreateArcGISandPDFMaps(), except it does not create PDFs. 
###  Blake Massey
###  2014.09.30

CreateArcGISMaps <- function(df = df){
  write.csv(df, file="C:/Work/Python/Data/CSV/BAEA.csv", row.names=FALSE)
  system('python C:/Work/Python/Scripts/BAEA/Create_ArcGIS_Maps.py')
}

# CreateArcGISandPDFMaps Function ----------------------------------------------

###  A wrapper function for creating an ArcGIS map and a series of PDFs
###  Usage: CreateArcGISandPDFMaps(df)
###  Arguments: df = dataframe with locations
###  Returns: overwrites files in "C:/Work/Python" folders          
###  Notes: internal parameters of Python script specific to my BAEA GPS data.
###    Identical to CreateArcGISMaps(), except it also creates PDFs. 
###  Blake Massey
###  2014.09.30

CreateArcGISandPDFMaps <- function(df = df){
  write.csv(df, file="C:/Work/Python/Data/CSV/BAEA.csv", row.names=FALSE)
  system('python C:/Work/Python/Scripts/BAEA/Create_ArcGISandPDF_Maps.py')
}

# CreateCategoricalLegend Function ---------------------------------------------

###  Creates a legend for categorical data
###  Usage: CreateCategoricalLegend (metadata, metadata_layer, color_alpha, 
###    pos_x, pos_y, main, main_cex, main_col, lab_col, signif_digits, 
###    left, plot_new)
###  Arguments: metadata = location of metadata .csv file. Metadata file must 
###               contain "id" and "icon_color" columns.
###             metadata_layer = not required, column name in metadata used 
###               for legend labels
###             color_alpha = two-digit rbga alpha value, "00" to "ff".
###             pos_x = relative position of left and right edge of color bar on 
###               first axis, [0,1]. Default is c(0.5, 0.55).
###             pos_y = relative position on lower and upper edge of colar bar 
###               on second axis, [0,1]. Default is c(0.05, 0.9). 
###             main = main title, written above the color bar. Default is NA.
###             main_cex = relative size of main title. Default is 1.
###             main_col = color of main title, default is "black". 
###             lab_col = color labels, default is "black". 
###             left = logical indicating whether to put the labels on the 
###               left (TRUE) or on the right (FALSE). Default is FALSE. 
###             plot_new = logical indicating whether to create a new plot, 
###               default is TRUE.
###  Returns: Plot of legend for categorical data 
###  Notes: Used in ExportKMLPolygon()
###  Blake Massey
###  2014.09.13

CreateCategoricalLegend <-function(metadata,
                                   metadata_layer = NULL,
                                   color_alpha = 99,
                                   pos_x = c(0.25, 0.3), 
                                   pos_y = c(0.25, 0.75), 
                                   main = NULL,
                                   main_cex = 1, 
                                   main_col = "black", 
                                   lab_col = "black",
                                   left = FALSE,
                                   plot_new = TRUE){
  suppressPackageStartupMessages(require(shape))
  suppressPackageStartupMessages(require(plotKML))
  if (plot_new == TRUE){
  plot.new()
  } 
  suppressMessages(par(new = TRUE))  # gave warning about calling par w/o a plot
  omar <- nmar <- par("mar")
  nmar[c(2, 4)] <- 0
  par(mar = nmar)
  emptyplot()
  pars <- par("usr")
  par(plt = c(0.01, .99, 0.01, 0.99)) # altered to extend graph plotting area 
  dim_x <- pars[2] - pars[1]
  x_min <- pars[1] + pos_x[1] * dim_x
  x_max <- pars[1] + pos_x[2] * dim_x
  dim_y <- pars[4] - pars[3]
  y_min <- pars[3] + pos_y[1] * dim_y
  y_max <- pars[3] + pos_y[2] * dim_y
  metadata<-read.csv(metadata, header=TRUE, as.is=TRUE, na.strings = "")
  if(!is.null(metadata_layer)){
  metadata$id<-metadata[,metadata_layer]
  }
  colors_kml <- metadata$icon_color
  colors <- colors_kml
  legend_interval_seq <- metadata$id  # this can be character or number 
  colors <- sub("#ff", color_alpha, colors_kml, ignore.case =FALSE)
  color_seq <- 1:length(legend_interval_seq) 
  z_min <- color_seq[1]
  z_max <- color_seq[length(color_seq)]
  y_seq <- seq(y_min, y_max, length.out = length(colors_kml) +1)
  rect(x_min, y_seq[-(length(colors_kml)+1)], x_max, y_seq[-1], col =colors_kml, 
    border = "black") # this draws the vertical color ramp
  dim_x <- (x_max - x_min)
  dim_y <- (y_max - y_min)
  if (left) {
    x_width <- -dim_x
    pos <- 2
    x_pos <- x_min + x_width * 0.1
  } else {
    x_width <- +dim_x
    pos <- 4
    x_pos <- x_max + x_width * 0.1
  }
  midpoints <- function(x){
    midpoints <-{}
    for (i in 1:length(x[-1])){
      lower <- x[i]
      upper <- x[i+1]
      midpoints <- append(midpoints,((lower+upper)/2))
    }
    return(midpoints)
  }
  y_pos <- midpoints(y_seq)    
  text(x_pos,  y_pos, legend_interval_seq, pos = 4, 
     col = lab_col)   
  if (!is.null(main)) {
    for (i in length(main):1) { 
      text(x = mean(c(x_min, x_max)), y = y_max + 0.075 * (length(main) - i+1), 
      labels = main[i], adj = c(0, 0.5), cex = main_cex, col = main_col)
    }
  }
  par(new = FALSE)
  par(mar = omar)
}

# CreateColorIntervalSequence Function -----------------------------------------

###  Creates an interval sequence used for assigning entities to a color ramp
###  Usage: CreateColorIntervalSequence(color_range, color_alpha, 
###    color_increment,log)
###  Arguments: color_range = two-valued vector, the minimum and maximum values.
###             color_levels = one value, number of color levels. Ignored if 
###               color_interval not equal to NA.
###             color_increment = one value, increment value of color 
###               breakpoints, one value.
###             log = logical of whether to log transform.
###  Returns: Vector of numbers
###  Notes: Parameter defaults are object names because this function is used 
###    internally in CreateColorPaletteLegend() & ExportKMLPolygon() to ensure 
###    the interval sequence selection is identical. The code uses pretty(), so 
###    the returned color levels and legend levels may not be identical to the 
###    input parameter value. 
###  Blake Massey
###  2014.08.12

CreateColorIntervalSequence<- function(color_range = color_range, 
                                       color_levels = color_levels, 
                                       color_increment = color_increment, 
                                       log = log) {
  if (log == TRUE && (min(color_range) <= 0) == TRUE ){  
    stop ("Unable to log-transform data: min(color_range) <= 0") 
  }
  if(!is.null(color_increment)){
    if (log){
      pretty_color_range <- 10^pretty(log10(color_range), n = (color_levels +1)) 
      color_interval_seq <- seq(log10(min(pretty_color_range)), 
        log10(max(pretty_color_range)), by = log10(color_increment)) 
    }
    if (!log){
      pretty_color_range <- pretty(color_range, n = (color_levels +1)) 
      color_interval_seq <- seq(min(color_range), max(color_range), 
        by = color_increment)
    }
  } else {
    if (log){
      color_interval_seq <- pretty(log10(color_range), n = (color_levels + 1)) 
    }
    if (!log){
      pretty_color_range <- pretty(color_range, n = (color_levels + 1))
      color_interval_seq <- seq(min(pretty_color_range),max(pretty_color_range), 
        length.out = color_levels + 1) 
    }
  }
  return(color_interval_seq)
}  

# CreateColorPaletteLegend Function --------------------------------------------

###  Creates a display of selected color palettes
###  Usage: CreateColorPaletteLegend (color_pal, color_range, color_alpha, 
###    color_levels, color_increment, legend_levels, log, pos_x, pos_y, main, 
###    main_cex, main_col, lab_col, signif_digits, left, plot_new)
###  Arguments: color_pal = color palette to be used, also allowed are two 
###               extremes or one value. Default is: c("yellow","red").
###             color_range = two-valued vector, the minimum and maximum values.
###             color_alpha = two-digit rbga alpha value, "00" to "ff", 
###               default is "ff". Note: ExportKMLPOlygon() uses default "cc"
###             color_levels = one value, number of color levels. Ignored if 
###               color_interval not equal to NA. Default is 5.
###             color_increment = one value, increment value of color 
###               breakpoints, one value. Default is NA. 
###             legend_levels = one value, increment in legend values, ignored  
###               if legend_values not equal to NA. Default is 5.
###             legend_values = vector of specific legend labels. Default is NA. 
###             log =  logical of whether to log transform. Default is FALSE. 
###             pos_x = relative position of left and right edge of color bar on 
###               first axis, [0,1]. Default is c(0.5, 0.55).
###             pos_y = relative position on lower and upper edge of colar bar 
###               on second axis, [0,1]. Default is c(0.05, 0.9). 
###             main = main title, written above the color bar. Default is NA.
###             main_cex = relative size of main title. Default is 1.
###             main_col = color of main title. Default is "black". 
###             lab_col = color labels, default is "black". 
###             signif_digits = integer, number of signifcant digits 
###             left = logical indicating whether to put the labels on the 
###               left (TRUE) or on the right (FALSE). Default is FALSE. 
###             plot_new = logical indicating whether to create a new plot. 
###               Default is TRUE.
###  Returns: Plot of legend with color palette  
###  Notes: Used in ExportKMLPolygon, make sure code used for color range, 
###    selction, log-tranformation, and color palette creation is identical or 
###    produces same results.  Color levels and legend levels use the pretty(), 
###    so the actual number of levels may not be identical to parameter values.
###  Blake Massey
###  2014.08.12

CreateColorPaletteLegend <- function (color_pal = c("yellow","red"), 
                                      color_range, 
                                      color_alpha = "ff",
                                      color_levels = 10,
                                      color_increment = NULL, 
                                      legend_levels = 10,
                                      legend_values = NULL,
                                      log = FALSE,
                                      pos_x = c(0.3, 0.35), 
                                      pos_y = c(0.0, 1), 
                                      main = NULL,
                                      main_cex = 1, 
                                      main_col = "black", 
                                      lab_col = "black",
                                      signif_digits = 2,
                                      left = FALSE,
                                      plot_new = TRUE){
  suppressPackageStartupMessages(require(shape))
  suppressPackageStartupMessages(require(plotKML))
  if (plot_new == TRUE){
  plot.new()
  } 
  suppressMessages(par(new = TRUE))  # gave warning about calling par w/o a plot
  omar <- nmar <- par("mar")
  nmar[c(2, 4)] <- 0
  par(mar = nmar)
  emptyplot()
  pars <- par("usr")
  par(plt = c(0.01, .99, 0.01, 0.99)) # altered to extend graph plotting area 
  dim_x <- pars[2] - pars[1]
  x_min <- pars[1] + pos_x[1] * dim_x
  x_max <- pars[1] + pos_x[2] * dim_x
  dim_y <- pars[4] - pars[3]
  y_min <- pars[3] + pos_y[1] * dim_y
  y_max <- pars[3] + pos_y[2] * dim_y
  color_interval_seq <- CreateColorIntervalSequence(color_range=color_range, 
    color_levels=color_levels, color_increment=color_increment, log=log)
  colors_rbg <- colorRampPalette(color_pal)(length(color_interval_seq)-1) 
  colors_rbg <- paste(colors_rbg, color_alpha, sep= "") 
  if (all(!is.null(legend_values))) {
    legend_interval_seq <- legend_values
  } else { # if only legend_levels is set
    if (log){
      legend_interval_seq <- pretty(log10(color_range), n = (legend_levels + 1)) 
    }
    if (!log){
      pretty_legend_range <- pretty(color_range, n = (legend_levels + 1))  
      legend_interval_seq <- seq(min(pretty_legend_range), 
        max(pretty_legend_range), length.out = (legend_levels + 1)) 
    }
  }
  z_min <- min(color_interval_seq)
  z_max <- max(color_interval_seq)
  y_seq <- seq(y_min, y_max, length.out = length(colors_rbg)+1 )
  rect(x_min, y_seq[-(length(colors_rbg)+1)], x_max, y_seq[-1], col=colors_rbg, 
    border=NA) 
  rect(x_min, y_min, x_max, y_max, border = lab_col)
    dim_x <- (x_max - x_min)
    dim_y <- (y_max - y_min)
    if (left) {
      x_width <- -dim_x
      pos <- 2
      x_pos <- x_min + x_width * 0.5
    } else {
      x_width <- +dim_x
      pos <- 4
      x_pos <- x_max + x_width * 0.5
    }
    y_pos <- y_min + (legend_interval_seq - z_min)/(z_max - z_min) * dim_y
    segments(x_min, y_pos, x_max, y_pos, col = lab_col)
    segments(x_pos + x_width * 0.2, y_pos, x_min, y_pos, col = lab_col)
  if(log==TRUE) { 
    legend_interval_labels <- sapply(legend_interval_seq, function(i)
      as.expression(bquote(10^ .(i))))
    } else {      
      legend_interval_labels <- format(legend_interval_seq, nsmall=0)
    }
  text(x_pos, y_pos, legend_interval_labels, pos = pos, col = lab_col)
  if (!is.null(main)) {
    for (i in length(main):1) { 
      text(x = mean(c(x_min, x_max)), y = y_max + 0.1 * (length(main) - i + 1), 
      labels = main[i], adj = c(0, 0.5), cex = main_cex, col = main_col)
    }
  }
  par(new = FALSE)
  par(mar = omar)
}

# CreateExtentBuffer Function --------------------------------------------------

###  CreateExtent Buffer function 
###  Usage: CreateExtentBuffer(df, buffer) 
###  Arguments: df = dataframe with columns: long_utm, lat_utm
###             buffer = minimum buffer extent, default = 500 
###  Returns: vector of extent (long_min, long_max, lat_min, lat_max)
###  Notes: finds lat and long extents, applies minimum buffer, and rounds to 
##    the same significance value as the buffer
###  Blake Massey
###  2014.07.04

CreateExtentBuffer <- function(df = df, 
                               buffer = 500) {
  suppressPackageStartupMessages(require(plyr))
  extent_matrix <- c(round_any(min(df$long_utm) - buffer, buffer, floor), 
    round_any(max(df$long_utm) + buffer, buffer, ceiling), 
    round_any(min(df$lat_utm) - buffer, buffer, floor),
    round_any(max(df$lat_utm) + buffer, buffer, ceiling))
  extent(extent_matrix)
}

# CreateKDEPoints Function -----------------------------------------------------

###  Creates a 'kde' object using a df of location points to calculate 
###    probability estimates for an input raster's cell centers        
###  Usage: CreateKDERaster(df, df_lat, df_long, blank_raster)
###  Arguments: df = dataframe of location data
###             df_long = df col name of longitude values, default is "long_utm"
###             df_lat = df col name of latitude values, default is "lat_utm"
###             buffer = value used in CreateExtentBuffer() around df locations,
###               default is 500 m 
###             blank_raster = file name or object name of raster file used to 
###               create output. This raster is cropped based on the df's extent
###               and kde estimates are based on each cell's center. Default 
###               file is: "C:/ArcGIS/Data/BlankRaster/maine_50mc.tif".
###  Returns: A 'kde' object with a kde estimate for the center point of each
###    blank raster cell that falls within the buffer around the df locations            
###  Notes: Ouput RasterLayer has the same CRS as the blank_raster input 
###  Blake Massey
###  2014.09.02

CreateKDEPoints <- function(df = df, 
                            df_long = "long_utm", 
                            df_lat = "lat_utm", 
                            buffer = 500, 
                            blank_raster = file.path("C:/ArcGIS/Data",
                              "BlankRaster/maine_50mc.tif")) {
  suppressPackageStartupMessages(require(ks))
  suppressPackageStartupMessages(require(raster))
  df_extent <- CreateExtentBuffer(df=df, buffer=500)
  if (is.raster(blank_raster) == TRUE){ 
    blank_raster <- blank_raster
  } else {
    blank_raster <- raster(blank_raster) 
  }
  blank_cropped <- crop(blank_raster, df_extent, snap="out")
  blank_points <- rasterToPoints(blank_cropped)
  df_xy <- cbind(df[,df_long], df[,df_lat])
  raster_xy <- cbind(blank_points[,1], blank_points[,2])
  hpi_df <- Hpi(x=df_xy, pilot="unconstr")
  kde_points <- kde(x=df_xy, H=hpi_df, eval.points=raster_xy, verbose=TRUE)
  return (kde_points)
}

# CreateKDEProbs Function ------------------------------------------------------

###  Calculates kernel probability levels based on 'kde' object and use them to 
###    recalssify the values in an kde 'RasterLayer' to the probability levels.
###  Usage: ExportKMLProbs(kde_df, kde_raster, probs)
###  Arguments: kde_object = 'kde' object, usually made by CreateKDEPoints()
###             kde_raster = 'RasterLayer' created using estimates from a 'kde' 
###               object, usually made by CreateKDERaster()
###             probs = vector [0,1] of probability values, default value is 
###               seq(0,1,by=.1)
###  Returns: A 'RasterLayer' with reclassified probabilities from using 
###    contourLevels() in the 'ks' package.           
###  Notes: To be used between CreateKDERaster() and ExportKMLProbs()
###  Blake Massey
###  2014.09.02

CreateKDEProbs <- function(kde_object = kde_object,
                           kde_raster = kde_raster, 
                           probs = seq(0,1,by=.1)) {
  suppressPackageStartupMessages(require(ks))
  suppressPackageStartupMessages(require(raster))
  kde_object <- kde_object
  kde_raster <- kde_raster
  probs <- probs
  cl <- (contourLevels(kde_object, prob=probs, approx=FALSE))
  cl2 <- c(cl[-length(cl)],1)  # top contour level will include all > probs
  min <- min(unlist(kde_raster@data@values))  # min raster value, not always 0
  m <- matrix(nrow=length(probs), ncol=3)
  if (1 %in% probs){
    m[,2]<-c(cl2[-1],1)
  } else {
    m[,2]<-c(cl[-1],1)
  }
  if (0 %in% probs){
    if (1 %in% probs){
      m[,1] <- c(min, cl2[-1])
    } else {
      m[,1] <- c(min,cl[-1])  
    }
  } else {
    if (1 %in% probs) {
      m[,1] <- cl2
    } else {
      m[,1] <- cl
    }
  }
  m[,3]<-c(probs)
  probs_raster <- reclassify(kde_raster, m)
  if (0 %in% probs == FALSE) {  
  m2 <- c(0, cl[1], NA)  # removes all raster values below first contour level 
  probs_raster <- reclassify(probs_raster, m2)
  }
  return(probs_raster)
}
  

# CreateKDERaster Function -----------------------------------------------------

###  Creates a 'RasterLayer' based on a 'kde' object of gridded point location 
###   probability estimates
###  Usage: CreateKDERaster(kde_points, cell_size, crs)
###  Arguments: kde_points = dataframe of kde estimates at gridded points
###             cell_size = cell size of output raster, default = 50 
###             crs = crs of output raster, default is UTM 19N, NAD83
###  Returns: A 'RasterLayer' object with kde estimates for each of the blank 
###    raster cells that fall within the buffered df locations.            
###  Notes: The kde_points file must have "eval.points" and "estimate", as 
###    created by CreateKDEPoints(). Returned 'RasterLayer' has the same CRS as 
###    blank_raster input. 
###  Blake Massey
###  2014.09.02  
  
CreateKDERaster <- function(kde_points = kde_points,
                            cell_size = 50,
                            crs = paste("+proj=utm +zone=19 +datum=NAD83",
                                        "+units=m +no_defs +ellps=GRS80", 
                                        "+towgs84=0,0,0")) {
  suppressPackageStartupMessages(require(raster))
  kde_points <- kde_points
  kde_xyz <- data.frame(kde_points["eval.points"], kde_points["estimate"])
  kde_raster <- rasterFromXYZ(kde_xyz, res=c(cell_size,cell_size), crs=crs, 
    digits=5)
  return(kde_raster)
}

# CreateProbIsoplethRaster Function --------------------------------------------

###  A wrapper function of used to create a RasterLayer of probabilities based 
###    on a df of location points      
###  Usage: CreateProbIsoplethRaster(df, buffer, cell_size, probs)
###  Arguments: df = dataframe of location data
###             buffer = value used in CreateExtentBuffer() around df locations,
###               default is 500 m
###             cell_size = cell size of output raster, default = 50
###             probs = vector [0,1] of probability values, default value is 
###               seq(0,1,by=.1)
###  Returns: A Raster object with a probability estimate for each
###    blank raster cell that falls within the buffer around the df locations            
###  Notes: Ouput RasterLayer has the same CRS as the blank_raster input 
###  Blake Massey
###  2014.09.30

CreateProbIsoplethRaster <- function(df,
                                     buffer = 500,
                                     cell_size = 50,
                                     probs = seq(0,1, by=.1)) {
  source('C:/Work/R/Functions/gis.R')
  kde_points <- CreateKDEPoints(df = df, df_long = "long_utm", 
    df_lat = "lat_utm", buffer = buffer, blank_raster = file.path("C:/ArcGIS",
    "Data/BlankRaster/maine_50mc.tif"))
  kde_raster <- CreateKDERaster(kde_points = kde_points, cell_size = cell_size,
    crs = paste("+proj=utm +zone=19 +datum=NAD83", "+units=m +no_defs",
    "+ellps=GRS80 +towgs84=0,0,0"))
  kde_probs <- CreateKDEProbs(kde_object = kde_points, kde_raster = kde_raster, 
    probs = probs) 
  return(kde_probs)
}

# CreateSpatialLines Function --------------------------------------------------

###  A wrapper function for creating a 'SpatialLines' object
###  Usage: CreateSpatialLines(df)
###  Arguments: df = dataframe with locations
###             long = df column name of longitude, default is "long_utm"
###             lat = df column name of latitude, default is "long_utm"
###  Returns: a 'SpatialLines'object          
###  Notes: 
###  Blake Massey
###  2014.10.13

CreateSpatialLines <- function (df = df,
                                long = "long_utm",
                                lat = "lat_utm"){
  suppressPackageStartupMessages(require(sp))
  df <- df
  coordinates(df) <- c(long, lat)
  coords <- coordinates(df)
  spatial_lines <- SpatialLines(list(Lines(list(Line(coords)), "1")))
  return(spatial_lines)
}

# ExportKMLPolygon Function ----------------------------------------------------

###  Create a Google Earth KML file from a SpatialPolygonsDataFrame
###  Usage: ExportKMLPolygon(object, object_layer, outfile, kml_name, 
###    categorical, metadata, metadata_layer, color_pal, color_alpha, 
###    color_range, color_min, color_max, color_levels, color_increment, 
###    legend_levels, legend_values, log, signif_digtits, outline, alt_mode, 
###    extrude, labelscale, create_kmz, ...) 
###  Arguments: object = a 'SpatialPolygonsDataFrame' object
###             object_layer = layer in object that is used for display
###             outfile = location of output KML file. Extensions (.kml or 
###               .kmz) will automatically determine file type.
###             kml_name = name of folder in KML file
###             categorical = logical, whether the data is categorical, default 
###               is false.
###             metadata = location of metadata .csv file for categorical data. 
###               Metadata file must contain "id" and "icon_color" columns.
###             metadata_layer = not required, column name in metadata used 
###               for legend labels instead of "id"
###             color_pal = color palette, can be a color ramp (i.e. c("white",
###               "red") or a specific palette ("SAGA_pal[[1]]")
###             color_alpha = display and legend alpha value, default "cc"
###             color_range = range of object values to create color palette 
###             color_min = min object values to create color palette
###             color_max = max object value to create color palette 
###             color_levels = number of breaks in color palette, ignored if 
###               color_interval not equal to NA 
###             color_increment = intervals for color palette breaks
###             legend_levels = number of breakpoints in legend, ignored if 
###               legend_interval not equal to NA 
###             legend_values = vector of values in legend 
###             log = logical of whether to log transform. Default is FALSE. 
###             signif_digits = number of signifcant digits for polygon labels
###             outline = 1 or 0, whether to draw an outline around each polygon
###             alt_mode = based on KML code: "absolute","clampedToGround",
###               "relativeToGround" (see KML documentation for description).
###               Default is "clampedToGround"
###             extrude = either 0 (default) or 1: 0 if for no line, 
###               1 extends a line from the point to the ground. 
###             labelscale = adjusts the size of the Google Earth location 
###               point labels. Default is 0, which hides the labels. To show 
###               labels, change to a value between 0.7-1.
###             create_kmz = will always create a KMZ, default is FALSE. 
###             ... = parameters passed to CreateCreateCategoricalLegend() or 
###               CreateColorPaletteLegend(), e.g., legend_levels, legend_values
###  Returns: KML or KMZ of polygons with an associated legend      
###  Notes: uses CreateCategoricalLegend() or CreateColorPaletteLegend() 
###    depending on if categorical = TRUE or FALSE 
###  Blake Massey
###  2014.09.15

ExportKMLPolygon <- function (object,
                              object_layer = NULL,
                              outfile = NULL,
                              kml_name = NULL,         
                              categorical = FALSE, 
                              metadata = NULL,
                              metadata_layer = NULL,
                              color_pal = c("yellow","red"),
                              color_alpha = "cc",
                              color_range = NULL,
                              color_min = NULL,
                              color_max = NULL,
                              color_levels = 5,
                              color_increment = NULL,                          
                              legend_levels = 5,
                              legend_values = NULL,
                              log = FALSE,
                              signif_digits = 3,
                              outline = 0,
                              alt_mode = "clampToGround",
                              extrude = 0,
                              labelscale = 0,
                              create_kmz = FALSE,
                              ...){   
  suppressPackageStartupMessages(require(maptools))
  suppressPackageStartupMessages(require(plotKML))
  suppressPackageStartupMessages(require(plyr))
  suppressPackageStartupMessages(require(tools))
  obj <- object
  if (class(obj) == "RasterLayer"){
    obj <- as(obj, 'SpatialPolygonsDataFrame')
  }
  ifelse(!is.null(object_layer), display_layer <- obj@data[[object_layer]], 
  display_layer <- obj@data[[1]])
  if(categorical == TRUE) {
    metadata_df<-read.csv(metadata, header=TRUE, as.is=TRUE, na.strings = "")
    if (!is.null(metadata_layer)){
      categorical_layer = metadata_df[,metadata_layer]
      metadata_display <- data.frame(display_layer = metadata_df$id,
        categorical_layer = categorical_layer)
    } else {
      metadata_display <- data.frame(display_layer = metadata_df$id)
    }
    display_df <- data.frame(display_layer)
    suppressMessages(display_layer <- join(display_df, metadata_display))
      #  join is based on match btwn. "display_layer" and "metadata_df$id"
  }
  obj <- spCbind(obj, display_layer)  # used to subset object
  if (is.null(kml_name) == TRUE) {
    if (!is.null(outfile)) {
      kml_name <- basename(outfile)
      kml_name <- sub(".kml", "", kml_name, ignore.case =TRUE)
      kml_name <- sub(".kmz", "", kml_name, ignore.case =TRUE)
    } else {
      kml_name <- deparse(substitute(object))
    }
  }
  if (is.null(outfile)) {
    if (!is.null(kml_name)) {
      outfile <- paste(getwd(), "/", kml_name, sep="")
    } else {
      outfile <- paste(getwd(), "/", deparse(substitute(object)), sep="")
    }
  } else {
  if (dirname(outfile) == "."){
    outfile <- paste(getwd(), "/", outfile, sep="") 
  }  
  }
  if (file_ext(outfile) == "kmz" |file_ext(outfile) == "KMZ") {
    create_kmz <- TRUE  # an outfile with a .kmz overrides create_kmz=FALSE
    outfile <- sub(".kmz", ".kml", outfile, ignore.case = TRUE)  # zip filename
  }
  if (file_ext(outfile) == "") {
    outfile <- paste(outfile, ".kml", sep="")  # if object has no extension
  }  
  if (categorical == FALSE) {  
  if (is.null(color_range[1])){  
    ifelse(!is.null(color_min), color_range[1] <- color_min, 
      color_range[1] <- min(obj$display_layer))
    ifelse(!is.null(color_max), color_range[2] <- color_max, 
      color_range[2] <- max(obj$display_layer)) 
  }
  color_interval_seq <- CreateColorIntervalSequence(color_range=color_range, 
    color_levels=color_levels, color_increment=color_increment, log=log)
    # function used so the results are identical to CreateColorPaletteLegend()
  colors_rbg <- colorRampPalette(color_pal)(length(color_interval_seq)-1) 
  colors_kml <- col2kml(colors_rbg)
  colors <- sub("#ff", color_alpha, colors_kml, ignore.case =FALSE)
  interval_num <- as.integer(seq(1, length(color_interval_seq)-1, 1))
  interval_colors <- data.frame(interval_num, colors)
  if(log){
  obj$display_layer <- log10(obj$display_layer) 
  }
  obj <- obj[obj$display_layer >= min(color_range) &
    obj$display_layer <= max(color_range),]  
  obj_display_layer <- obj$display_layer
  poly_style_num <- findInterval(obj_display_layer, color_interval_seq,
    rightmost.closed=TRUE)
  if (log == TRUE) { 
    obj_display_layer <- format(obj$display_layer, nsmall=0)
  } else {      
    obj_display_layer <- format(obj$display_layer, nsmall=0)
  }
  }
  prj.check <- check_projection(obj, control = TRUE)
  if (!prj.check) {
    obj <- reproject(obj)
  }
  pv <- length(obj@polygons) 
  pvn <- lapply(lapply(obj@polygons, slot, "Polygons"), length)
  coords <- rep(list(NULL), pv)
  hole <- rep(list(NULL), pv)  # currently not used.
  display_values <- rep(list(NULL), pv)
  categorical_values <- rep(list(NULL), pv)
  labpts <- rep(list(NULL), pv)
  for (i in 1:pv) {
    for (k in 1:pvn[[i]]) {
      poly_xyz <- slot(slot(obj@polygons[[i]], "Polygons")[[k]], "coords")
      pt_xyz <- slot(slot(obj@polygons[[i]], "Polygons")[[k]], "labpt")
      if (ncol(poly_xyz) == 2) {
        poly_xyz <- cbind(poly_xyz, rep(0, nrow(poly_xyz)))
      }
      hole[[i]][[k]] <- slot(slot(obj@polygons[[i]], "Polygons")[[k]], "hole")
      coords[[i]][[k]] <- paste(poly_xyz[, 1], ",", poly_xyz[, 2], ",", 
        poly_xyz[, 3], collapse = "\n ", sep = "")
      display_values[[i]][[k]] <- as.numeric(obj[["display_layer"]][[i]])
        # as numeric should remove any factor 
      if (categorical == TRUE) {
        categorical_values[[i]][[k]] <- 
          as.character(obj[["categorical_layer"]][[i]])      
      }
      labpts[[i]][[k]] <- paste(pt_xyz[1], ",", pt_xyz[2], ",", 0, 
      collapse = "\n ", sep = "")
    }
  }  
  poly_coords <- rep(list(""), pv)
  hole_coords <- rep(list(""), pv) 
  lab_pts <- rep(list(NA), pv)  
  poly_names <- rep(matrix(NA), pv)
  for (i in 1:pv) {
    for (k in 1:pvn[[i]]) {
     lab_pts[i] <- paste(unlist(labpts[i]), collapse = "\n ", sep = "")
     if (hole[[i]][[k]] == TRUE) {
       hole_coords[[i]][[k]] <- paste("\t\t\t\t\t<LinearRing>\n",
       "\t\t\t\t<altitudeMode>", alt_mode, "</altitudeMode>\n",
       "\t\t\t\t\t\t<coordinates>\n", coords[[i]][[k]], "\n",
       "\t\t\t\t\t\t</coordinates>\n",
       "\t\t\t\t\t</LinearRing>\n")
     }        
     if (hole[[i]][[k]] == FALSE) {
       poly_coords[[i]][[k]] <- paste("\t\t\t\t\t<LinearRing>\n",
       "\t\t\t\t\t<altitudeMode>", alt_mode, "</altitudeMode>\n",
       "\t\t\t\t\t\t<coordinates>\n", coords[[i]][[k]], "\n",
       "\t\t\t\t\t\t</coordinates>\n",
       "\t\t\t\t\t</LinearRing>\n")
     }
     if(is.numeric(display_values[[i]][[1]])){
       poly_names[[i]] <- signif(display_values[[i]][[1]], digits=signif_digits)    
     } else {
       poly_names[[i]] <- display_values[[i]][[1]]
     }
    }
  }
  hole_coords2<- lapply(hole_coords, paste, collapse="")    
  poly_coords2<- lapply(poly_coords, paste, collapse="")  
  if (categorical == TRUE) {
  poly_names <- categorical_values
  poly_style_num <- display_values
  interval_num <- metadata_df$id
  colors <- sub("#ff", color_alpha, col2kml(metadata_df$icon_color)) 
  }
  placemarks <- sprintf(paste(
  "\t<Placemark>\n",
  "\t\t<name>%s</name>\n",  
  "\t\t\t<styleUrl>#Poly%s</styleUrl>\n", 
  "\t\t\t<MultiGeometry>\n",
  "\t\t\t<Point>\n",  
  "\t\t\t\t\t\t<coordinates>\n%s\n",
  "\t\t\t\t\t\t</coordinates>\n",
  "\t\t\t</Point>\n",  
  "\t\t\t<Polygon>\n",
  "\t\t\t\t<extrude>", extrude, "</extrude>\n",
  "\t\t\t\t<outerBoundaryIs>\n%s",
  "\t\t\t\t</outerBoundaryIs>\n",
  "\t\t\t\t<innerBoundaryIs>\n%s",
  "\t\t\t\t</innerBoundaryIs>\n",
  "\t\t\t</Polygon>\n",
  "\t\t\t</MultiGeometry>\n",
  "\t</Placemark>\n", sep=""), poly_names, poly_style_num, lab_pts,
  poly_coords2, hole_coords2)  
  ## Icon Style Section ##
  hi_icon_label_scale <- 0.75   
  icon_scale <- 0.7
  ball_bg_color <- "ff333333"
  ball_text_color <- "ffffffff"
  stylemaps <- sprintf(paste(
  "\t<StyleMap id=\"Poly%s\">\n",
  "\t\t<Pair>\n",
  "\t\t\t<key>normal</key>\n",
  "\t\t\t\t<Style>\n",  
  "\t\t\t\t\t<LabelStyle>\n",      
  "\t\t\t\t\t<scale>",labelscale,"</scale>\n",  # to show label 
  "\t\t\t\t\t</LabelStyle>\n", 
  "\t\t\t\t\t<IconStyle>\n",
  "\t\t\t\t\t\t<scale>0</scale>\n", 
  "\t\t\t\t\t</IconStyle>\n", 
  "\t\t\t\t\t<BalloonStyle>\n",
  "\t\t\t\t\t<text>$[description]</text>\n", 
  "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>\n",
  "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>\n",
  "\t\t\t\t\t</BalloonStyle>\n",
  "\t\t\t\t\t<PolyStyle>\n",
  "\t\t\t\t\t\t<outline>",outline,"</outline>\n",
  "\t\t\t\t\t\t<color>%s</color>\n",
  "\t\t\t\t\t\t<colorMode>normal</colorMode>\n",
  "\t\t\t\t\t</PolyStyle>\n",
  "\t\t\t\t</Style>\n",
  "\t\t</Pair>\n",        
  "\t\t<Pair>\n",
  "\t\t\t<key>highlight</key>\n",
  "\t\t\t\t<Style>\n",
  "\t\t\t\t\t<LabelStyle>\n",     
  "\t\t\t\t\t<scale>",hi_icon_label_scale,"</scale>\n",  # to show label 
  "\t\t\t\t\t</LabelStyle>\n",
  "\t\t\t\t\t<IconStyle>\n",
  "\t\t\t\t\t\t<scale>0</scale>\n", 
  "\t\t\t\t\t</IconStyle>\n",
  "\t\t\t\t\t<BalloonStyle>\n",
  "\t\t\t\t\t<text>$[description]</text>\n", 
  "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>", "\n",
  "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>", "\n",
  "\t\t\t\t\t</BalloonStyle>\n",
  "\t\t\t\t\t<PolyStyle>\n",
  "\t\t\t\t\t\t<outline>1</outline>\n",
  "\t\t\t\t\t\t<color>%s</color>\n",
  "\t\t\t\t\t\t<colorMode>normal</colorMode>\n",
  "\t\t\t\t\t</PolyStyle>\n",
  "\t\t\t\t</Style>\n", 
  "\t\t</Pair>\n",
  "\t</StyleMap>\n", 
  sep = ""), interval_num, colors, colors)
  if (create_kmz == TRUE){
    base_file_name <- basename(outfile)
    org_outfile <- outfile
    temp_dir <- file.path(dirname(outfile), "TEMP")
    temp_files_dir <- file.path(temp_dir, "files")
    dir.create(file.path(temp_dir), showWarnings = FALSE)
    dir.create(file.path(temp_files_dir), showWarnings = FALSE)
    png_name <- sub(".kml", " - Legend.png", base_file_name, ignore.case =TRUE)
    href_png <- file.path ("files", png_name)    
    png_outfile <- file.path(temp_dir, href_png)
    kmz_outfile<- file.path(temp_dir, base_file_name)
    outfile <- file.path(temp_dir, base_file_name)
  } else {
    png_name <- sub(".kml", " - Legend.png", basename(outfile), 
      ignore.case = TRUE)
    png_outfile <- file.path(dirname(outfile), png_name)
    href_png <- png_outfile  # for noticeablility, legend in same folder as .kml
  }
  png(png_outfile, width = 350, height = 600, units = "px", bg = "transparent",
    res = 250, pointsize = 5.5, type = "cairo")  # sets png output parameters
  if(categorical == TRUE){
    CreateCategoricalLegend (metadata=metadata, metadata_layer=metadata_layer, 
      color_alpha = color_alpha, pos_x = c(0.05, 0.15), pos_y = c(0.00, 0.99), 
      main = kml_name, main_cex = 1.1, main_col = "white", lab_col = "white",
      ...)
  } else {  
    CreateColorPaletteLegend(color_pal = colors_rbg, color_alpha = color_alpha, 
      color_range = color_range, color_levels = color_levels, 
      color_increment = color_increment, legend_levels = legend_levels, 
      legend_values = legend_values, log=log, signif_digits = signif_digits, 
      pos_x = c(0.05, 0.15), pos_y = c(0.00, 0.99), main = kml_name, 
      main_cex = 1.1, main_col = "white", lab_col ="white", ...) 
  }
  dev.off()  # writes png_ouput file 
  screenoverlay <- paste("\t<ScreenOverlay>\n",
  "\t\t<name>Legend</name>\n",
  "\t\t<Icon>\n",
  "\t\t<href>",href_png,"</href>\n",
  "\t\t</Icon>\n",
  "\t\t<overlayXY x=\"0\" y=\"1\" xunits=\"fraction\" yunits=\"fraction\"/>\n",
  "\t\t<screenXY x=\"0\" y=\"1\" xunits=\"fraction\" yunits=\"fraction\"/>\n",
  "\t\t<rotationXY x=\"0.5\" y=\"0.5\" xunits=\"fraction\" yunits=\"fraction\"",
  "/>\n", "\t\t<size x=\"-1\" y=\"-1\" xunits=\"pixels\" yunits=\"pixels\"/>\n",
  "\t</ScreenOverlay>\n", sep="")
  if (file.exists(outfile)) file.remove(outfile)  # delete KML if already exists
  cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n",      
  "<kml xmlns=\"http://www.opengis.net/kml/2.2\"\n",      
  "xmlns:gx=\"http://www.google.com/kml/ext/2.2\"\n",
  "xmlns:atom=\"http://www.w3.org/2005/Atom\">\n\n",
  "<Document>\n", "\t<name>",kml_name,"</name>\n",file = outfile, 
  append = FALSE, sep = "") 
  cat(screenoverlay, file = outfile, append = TRUE, sep = "")  # all stylemaps  
  cat(stylemaps,file = outfile, append = TRUE, sep = "")  # all stylemaps    
  cat("\t<Folder>\n","\t<name>",kml_name,"</name>\n", "\t<open>0</open>\n", 
    file = outfile, append = TRUE, sep = "")    
  cat(placemarks,file = outfile, append = TRUE, sep = "")  # all placemarks   
  cat("\t</Folder>\n", file = outfile, append = TRUE, sep = "")
  cat("</Document>\n</kml>", file = outfile, append = TRUE)
  if (create_kmz == TRUE){
    zip_file <- sub(".kml", ".zip", org_outfile, ignore.case =TRUE)
    kmz_file <- sub(".kml", ".kmz", org_outfile, ignore.case =TRUE)
    command <- paste("7z a ","\"",zip_file,"\" ", "\"",outfile,"\" ", "\"", 
      temp_files_dir,"\"",  sep="")
    system(command)  # runs as if from the fcommand prompt 
    file.rename(zip_file, kmz_file)
    do.call(file.remove,list(list.files(temp_files_dir, full.names=TRUE)))
    unlink(temp_files_dir, recursive=TRUE)
    do.call(file.remove,list(list.files(temp_dir, full.names=TRUE)))    
    unlink(temp_dir, recursive=TRUE)
    writeLines(noquote(c("Writing: ", kmz_file)))
  } else {
  png_file <- file.path(dirname(outfile), png_outfile)
  writeLines(noquote(c("Writing: ", outfile, png_outfile)))
  }
}

# ExportKMLProbContour Function ------------------------------------------------

###  Create a Google Earth KML or KMZ file based on a 'RasterLayer' of kernel 
###    probability levels
###  Usage: ExportKMLProbContour(kde_raster, kml_name, dissolve, ...)
###  Arguments: probs_raster = RasterLayer of kde probability levels, usually 
###               made with CreateHRProbContourRaster()
###             kml_name = name of folder in KML file, default= "prob_countours"
###             dissolve = logical, TRUE = dissolve adjacent like polygons, 
###               FALSE = each cell remains intact; default is FALSE
###             ... = other arguments to pass to ExportKMLPolygon() 
###  Returns: Exports a KML or KMZ of the probability contour levels           
###  Notes: Use the dissolve functionality with caution - it is unable to 
###    properly handle polygon holes. This results from a known problem 
###    regarding holes in 'Polygon' objects. Therefore, at a given contour level
###    holes may be improperly filled in.
###  Blake Massey
###  2014.09.02
  
ExportKMLProbContour <- function(probs_raster = probs_raster,
                                 kml_name = "prob_contours",
                                 dissolve = FALSE, 
                                 ...){
  suppressPackageStartupMessages(require(raster))
  probs_raster <- probs_raster
  probs_polys <- rasterToPolygons(probs_raster, dissolve=dissolve) 
    # 0 !in probs -> huge 
  poly <- {}
  ID_list <- {}
  for (k in 1:length(probs_polys@polygons)){
    poly[k] <- lapply(probs_polys@polygons[k], slot, "Polygons")
    ID_list[k] <- probs_polys@data$layer[k] 
  }
  poly_ID_list <- {}
  for (k in 1:length(ID_list)){
    poly_ID_list <- append(poly_ID_list, rep(ID_list[k], length(poly[[k]])))
    }
  poly_list <- unlist(poly) 
  polys <- {}
  for (k in 1:length(poly_list)){
    polys <- c(polys, Polygons(poly_list[k],k))
  }
  # this may be the spot to add checkPolygonsHoles (but only if dissolve = TRUE)
  spolys = SpatialPolygons(unlist(polys), 1:length(unlist(polys)))
  spolysdf <- SpatialPolygonsDataFrame(spolys, data = data.frame(data = 
    poly_ID_list))
  projection(spolysdf) <- projection(probs_raster)
  ExportKMLPolygon(spolysdf, kml_name = kml_name, ...)
}

# ExportKMLRaster Function -----------------------------------------------------

###  Create a Google Earth KML file from a RasterLayer or a set of KMLs from a 
###    RasterStack, or RasterBrick. Reliant on ExportKMLPolygon(). 
###  Usage: ExportKMLRaster(object, object_layer, outfile, kml_name,categorical,
###    metadata, metadata_layer, color_pal, color_alpha, color_range, color_min,
###    color_max, color_levels, color_increment, outline, alt_mode, extrude, 
###    labelscale, create_kmz) 
###  Arguments: object = a 'RasterLayer' or 'SpatialPolygonsDataFrame' object
###             object_layer = layer in object that is used for display
###             outfile = location of output KML file. Extensions (.kml or 
###               .kmz) will automatically determine file type.
###             kml_name = name of folder in KML file             
###             categorical = logical, whether the data is categorical, default 
###               is false.
###             metadata = location of metadata .csv file for categorical data. 
###               Metadata file must contain "id" and "icon_color" columns.
###             metadata_layer = not required, column name in metadata used 
###               for legend labels instead of "id"
###             color_pal = color palette, can be a color ramp (i.e. c("white",
###               "red") or a specific palette ("SAGA_pal[[1]]"). Atuomatically 
###               set as Saga_pal[i] for RasterStack and RasterBrick. 
###             color_alpha = display and legend alpha value, default "cc"
###             color_range = range of object values to create color palette 
###             color_min = min object values to create color palette
###             color_max = max object value to create color palette 
###             color_levels = number of breaks in color palette, ignored if 
###               color_interval not equal to NA 
###             color_increment = intervals for color palette breaks
###             legend_levels = number of breakpoints in legend, ignored if 
###               legend_interval not equal to NA 
###             legend_values = vector of values in legend 
###             log = logical of whether to log transform. Default is FALSE. 
###             signif_digits = number of signifcant digits for polygon labels
###             outline = 1 or 0, whether to draw an outline around each polygon
###             alt_mode = based on KML code: "absolute","clampedToGround",
###               "relativeToGround" (see KML documentation for description).
###               Default is "clampedToGround"
###             extrude = either 0 (default) or 1: 0 for no line, 1 extends a l
###               ine from the point to the ground. 
###             labelscale = adjusts the size of the Google Earth location 
###               point labels. Default is 0, which hides the labels. To show 
###               labels, change to a value between 0.7-1.
###             create_kmz = will always create a KMZ, default is FALSE. 
###  Returns: KML or KMZ of polygons with an associated legend      
###  Notes: uses ExportKMLPolygon() and CreateColorPaletteLegend() 
###  Blake Massey
###  2014.09.15

ExportKMLRaster <- function (object = object,
                              object_layer = NULL,
                              outfile = NULL,
                              kml_name = NULL,
                              categorical = FALSE,
                              metadata = NULL,
                              metadata_layer = NULL,
                              color_pal = c("yellow","red"),
                              color_alpha = "cc",
                              color_range = NULL,
                              color_min = NULL,
                              color_max = NULL,
                              color_levels = 5,
                              color_increment = NULL,                          
                              legend_levels = 5,
                              legend_values = NULL,
                              log = FALSE,
                              signif_digits = 3,
                              outline = 0,
                              alt_mode = "clampToGround",
                              extrude = 0,
                              labelscale = 0,
                              create_kmz = FALSE) {   
  suppressPackageStartupMessages(require(maptools))
  suppressPackageStartupMessages(require(plotKML))
  suppressPackageStartupMessages(require(plyr))
  suppressPackageStartupMessages(require(tools))
  object <- object
  object_layer_org <- object_layer
  outfile_org <- outfile
  kml_name_org <- kml_name
  categorical_org <- categorical
  metadata_org <- metadata
  metadata_layer_org <- metadata_layer
  color_pal_org <- color_pal
  color_alpha_org <- color_alpha
  color_range_org <- color_range
  color_min_org <- color_min
  color_max_org <- color_max
  color_levels_org <- color_levels
  color_increment_org <- color_increment                          
  legend_levels_org <- legend_levels
  legend_values_org <- legend_values
  log_org <- log
  signif_digits_org <- signif_digits
  outline_org <- outline
  alt_mode_org <- alt_mode
  extrude_org <- extrude
  labelscale_org <- labelscale 
  create_kmz_org <- create_kmz
  specific_raster_list <- c("elev_50mc", "hydro_dir_50mc", "hydro_dist_50mc", 
    "lc_50mc", "maine_50mc")
  for (i in 1:nlayers(object)) {
    if (names(object[[i]]) %in% specific_raster_list) {
      if (names(object[[i]]) == "elev_50mc") {
        color_pal = R_pal[[3]]
        color_min = 1
        color_max = 500  # Mt Kathadin is around 5500
        color_levels = 50
        legend_levels = 10
      }
      if (names(object[[i]]) == "hydro_dir_50mc") {
        color_pal = R_pal[[10]]
        color_min = 0
        color_max = 360
        color_levels = 12
        legend_levels = 12
      }      
      if (names(object[[i]]) == "hydro_dist_50mc") {
        color_pal = R_pal[[9]]
        color_min = 0
        color_max = 2000
        color_levels = 20
        legend_levels = 10
      }
      if (names(object[[i]]) == "lc_50mc") {
        categorical = TRUE 
        metadata = "C:/Work/R/Data/Mapping/lc_50mc.csv"
        metadata_layer = "Land_Cover"
      }
      if (names(object[[i]]) == "maine_50mc") {
        color_pal = SAGA_pal[[22]]
        color_min = 0
        color_max = 10
        color_levels = 10
        legend_levels = 10
      }      
      ExportKMLPolygon(object = object[[i]], object_layer = object_layer, 
        outfile = outfile, kml_name = names(object[[i]]), categorical = 
        categorical, metadata = metadata, metadata_layer, color_pal = color_pal, 
        color_alpha = color_alpha, color_range = color_range, color_min = 
        color_min, color_max = color_max, color_levels = color_levels, 
        color_increment = color_increment, legend_levels = legend_levels, 
        legend_values = legend_values, log = log, signif_digits = signif_digits, 
        outline = outline, alt_mode = alt_mode, extrude = extrude, 
        labelscale = labelscale, create_kmz = create_kmz)
      # This section returns all parameter values back to their original values
      object_layer <- object_layer_org
      outfile <- outfile_org
      kml_name <- kml_name_org
      categorical <- categorical_org
      metadata <- metadata_org
      metadata_layer <- metadata_layer_org
      color_pal <- color_pal_org
      color_alpha <- color_alpha_org
      color_range <- color_range_org
      color_min <- color_min_org
      color_max <- color_max_org
      color_levels <- color_levels_org
      color_increment <- color_increment_org                          
      legend_levels <- legend_levels_org
      legend_values <- legend_values_org
      log <- log_org
      signif_digits <- signif_digits_org
      outline <- outline_org
      alt_mode <- alt_mode_org
      extrude <- extrude_org
      labelscale <- labelscale_org 
      create_kmz <- create_kmz_org
    } else {
      if (nlayers(object) > 1) { 
        ExportKMLPolygon(object = object[[i]], object_layer = object_layer, 
          outfile = outfile, kml_name = names(object[[i]]), categorical = 
          categorical, metadata = metadata, metadata_layer, color_pal = 
          SAGA_pal[[i]], color_alpha = color_alpha, color_range = color_range, 
          color_min=color_min, color_max = color_max, color_levels = 
          color_levels, color_increment = color_increment, legend_levels = 
          legend_levels, legend_values = legend_values, log = log, 
          signif_digits = signif_digits,outline = outline, alt_mode = alt_mode, 
          extrude = extrude, labelscale = labelscale, create_kmz = create_kmz)
      } else {
        ExportKMLPolygon(object = object, object_layer = object_layer, outfile = 
          outfile, kml_name = kml_name, categorical = categorical, metadata =
          metadata, metadata_layer, color_pal = color_pal, color_alpha = 
          color_alpha, color_range = color_range, color_min=color_min, 
          color_max = color_max, color_levels = color_levels, color_increment = 
          color_increment, legend_levels = legend_levels, legend_values = 
          legend_values, log = log, signif_digits = signif_digits, outline = 
          outline, alt_mode = alt_mode, extrude = extrude, labelscale =
          labelscale, create_kmz = create_kmz)
      }
    }
  }
}

# ExportKMLRasterOverlay Function ----------------------------------------------

###  Export KML Raster function 
###  Usage: ExportKMLRasterOverlay(x, outfile, outfolder) 
###  Arguments: x = Raster* object
###             color_pal = color palette, can be a color ramp (i.e. c("white",
###               "red") or a specific palette ("SAGA_pal[[1]]")
###             outfile = name of KML. default is to use name of raster 
###             outfolder = folder location. default = "C:/ArcGIS/Data/R_Output"
###  Returns: KML of a Raster
###  Notes: 
###  Blake Massey
###  2014.07.04

ExportKMLRasterOverlay<-function(x = x, 
                                 color_pal = c("white","red"), 
                                 outfile = NULL, 
                                 outfolder="C:/ArcGIS/Data/R_Output") {
  suppressPackageStartupMessages(require(raster))
  nm <-deparse(substitute(x))
  if(is.null(outfile)){
  outfile <- paste(outfolder,"/",nm, ".kml", sep ="")
  } else {
  nm <- outfile
  outfile <- paste(outfolder,"/",nm, ".kml", sep ="")  
  }
  x <- projectRaster(x, crs="+proj=longlat +datum=WGS84", method='ngb')
  unique_x <- length(unique(round(getValues(x))))
  cols <- colorRampPalette(color_pal)(unique_x)
  KML(x, file=outfile, col=cols, colNA=NA, maxpixels=500000, blur=10, zip='', 
    overwrite=TRUE)
}

# # ExportKMLTelemetry (moved) Function --------------------------------------------------
# 
# ###  Create a Google Earth KML file (points and multitrack) from lat/long 
# ###    coordinates
# ###  Usage: ExportKMLTelemetry(df, metadata, outfile, kml_folder, id, datetime, 
# ###    lat, long, speed, alt, alt_mode, behavior, behavior_color, path, extrude, 
# ###    labelscale, dateformat, timeformat, datetimeformat) 
# ###  Arguments: df = input dataframe; must have ID, lat, long, and datetime 
# ###             metadata = location of metadata .csv file. Metadata file must 
# ###               contain columns for animal ID, hexadecimal colors, and
# ###               additional icon data.
# ###             outfile = location of output KML file
# ###             kml_folder = name for parent KML folder
# ###             id = column name of unique identifier
# ###             datetime = column name of datetime in POSIXct format or as a 
# ###               character in the format (%Y/%m/%d %H:%M)
# ###             lat = column name of latitude coordinates (WGS84, dec. degree)
# ###             long = column name of longitude coordinates (WGS84, dec. degree)
# ###             speed = input dataframe column name for speed. Optional
# ###             alt = input dataframe column name for altitude(m). Optional.
# ###             alt_mode = based on KML code: "absolute","clampedToGround",
# ###               "relativeToGround" (see KML documentation for description).
# ###               Default is "clampedToGround"
# ###             behavior = input dataframe column name for behavior. Optional
# ###             behavior_color = colors are determined by behavior, shape by sex
# ###             path = to create paths. Default is TRUE.
# ###             extrude = either 0 (default) or 1: 0 if for no line, 
# ###               1 extends a line from the point to the ground. 
# ###             labelscale = adjusts the size of the Google Earth location 
# ###               point labels. Default is 0, which hides the labels. To show 
# ###               labels, change to a value between 0.7-1.
# ###             dateformat = changes the format of the date in the Google Earth 
# ###               location pop-up windows. Default is "%Y/%m/%d".
# ###             timeformat = changes the format of the time in the Google Earth  
# ###               locations pop-up windows. Default is "%I:%M %p".        
# ###             datetimeformat = changes the datetime format for the label of 
# ###               highlighted points. Default is "%Y/%m/%d %I:%M %p"
# ###  Returns: KML of points and multitracks          
# ###  Notes: 
# ###  Blake Massey
# ###  2014.06.02
# 
# ExportKMLTelemetry <- function (df = df,
#                                 metadata = file.path("C:/Work/R/Data",
#                                   "BAEA/BAEA_gps_data.csv"),
#                                 outfile ="C:/Users/Blake/Desktop/BAEA Data.kml",
#                                 kml_folder = "BAEA Data",
#                                 id = "id",                       
#                                 datetime = "datetime",
#                                 lat = "lat",
#                                 long = "long",
#                                 speed = NULL,
#                                 alt = NULL,
#                                 alt_mode = "clampToGround",
#                                 behavior = NULL,
#                                 behavior_color = FALSE,
#                                 path = TRUE,
#                                 extrude = 0,
#                                 labelscale = 0, 
#                                 dateformat = "%Y-%m-%d", 
#                                 timeformat = "%I:%M %p",
#                                 datetimeformat = "%Y-%m-%d %I:%M %p") {
#   suppressPackageStartupMessages(require(plotKML))
#   df$id <- df[ ,id]  # write id to the "id" column
#   ifelse(is.numeric(df$id), id_type<-"serial", id_type<-"deploy_location")
#   df$lat <- df[ ,lat]  # writes lat to the "lat"column
#   df$long <- df[ ,long]  # writes longitude to the "long" column
#   if (alt_mode == "absolute" | alt_mode == "relativeToGround") {
#   df$alt <- df[ ,alt]  # writes altitude to the "alt" column
#   alt1 <- '\t\t\t\t\tAltitude: '  # writes first part of "Altitude" description
#   alt2 <- '\n'  # writes second part of "Altitude" description
#   } else {
#   df$alt <- ""  # makes the altitude column a vector of blank values  
#   alt1 <- NULL  # prevents "Altitude" description from being written
#   alt2 <- NULL  # prevents "Altitude" description from being written    
#   }
#   if (!is.null(speed)) {
#   df$speed <- df[,speed]  # writes speed to the "speed" column
#   spd1 <- '\t\t\t\t\tSpeed: '  # writes first part of the "Speed" description
#   spd2 <- '\n'  # writes second part of the "Speed" description
#   } else {
#   df$speed <- ""  # makes the speed column a vector of blank values  
#   spd1 <- NULL  # prevents "Speed" description from being written
#   spd2 <- NULL  # prevents "Speed" description from being written    
#   }
#   if (!is.null(behavior)) {
#   df$behavior <- df[,behavior]  # writes behavior to the "behavior" column
#   beh1 <- '\t\t\t\t\tBehavior: '  # first part of the "Behavior" description
#   beh2 <- '\n'  # second part of the "Behavior" description
#   } else {
#   df$behavior <- ""  # makes the behavior column a vector of blank values  
#   beh1 <- NULL  # prevents "Behavior" description from being written
#   beh2 <- NULL  # prevents "Behavior" description from being written    
#   }
#   df$datetime <- df[,datetime]
#   df$datetime <- as.character(df$datetime)  # needed for KML parsing
#   df$datetimebegin <- df$datetime
#   end_times <- function(data) {   # returns "datetimeend" 
#     data$datetime2 <- data$datetime[c(2:length(data$datetime), 
#     length(data$datetime))]
#   } 
#   id <- as.character(df$id)  # vector of ids, as.character removes factor levels
#   df_split <- split(df, id)  # divides data by ids 
#   df_split <- lapply(df_split, end_times)  # 'end_times' on each section of data
#   datetimeend <- unsplit(df_split, id)  # returns array of returned values
#   df <- cbind(df, datetimeend)  # adds datetimeend column to original baea data
#   #if (behavior_color == TRUE) {
#   PlacemarkPoint <- function(plname, X, Y, Z, SP, BH, ID, SD, ST, ED, ET, DA,
#     TM) {
#   ifelse(behavior_color == TRUE, icon_style <- paste(sex,"-",BH,sep=""), 
#     icon_style <- unique.id)
#   cat("\t<Placemark>\n",
#   "\t\t<name>", plname, "</name>\n",
#   "\t\t<TimeSpan>\n",
#   "\t\t\t<begin>",SD ,"T" ,ST ,"</begin> " , "\n",
#   "\t\t\t<end>", ED, "T", ET, "</end> ", "\n",
#   "\t\t</TimeSpan>\n",
#   "\t\t\t<Snippet></Snippet>", "\n",
#   "\t\t\t\t<description>\n",
#   "\t\t\t\t\tID: ", ID, "\n",
#   "\t\t\t\t\tDate: ", DA, "\n",
#   "\t\t\t\t\tTime: ", TM, "\n",
#   "\t\t\t\t\tLongitude: ", X, "\n",
#   "\t\t\t\t\tLatitude: ", Y, "\n",
#   alt1, Z, alt2,  # written when alt is "absolute" or "relativeToGround"
#   spd1, SP, spd2,  # written when speed exists
#   beh1, BH, beh2, # written when behavior exists
#   "\t\t\t\t</description>", "\n",  
#   "\t\t\t<styleUrl>#Point_",icon_style,"</styleUrl>\n", 
#   "\t\t\t<Point>\n",
#   "\t\t\t\t<extrude>", extrude, "</extrude>\n",
#   "\t\t\t\t<altitudeMode>", alt_mode, "</altitudeMode>\n",
#   "\t\t\t\t\t<coordinates>", X, ",", Y, ",", Z, "</coordinates>\n",
#   "\t\t\t</Point>\n",
#   "\t</Placemark>\n",
#   file = outfile, append = TRUE, sep = "")  
#   }    
#   if (file.exists(outfile)) file.remove(outfile)  # delete KML if already exists
#   writeLines(noquote(c("Writing: ", outfile)))
#   ## Title Section ##
#   cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n",      
#   "<kml xmlns=\"http://www.opengis.net/kml/2.2\"\n",      
#   "xmlns:gx=\"http://www.google.com/kml/ext/2.2\"\n",
#   "xmlns:atom=\"http://www.w3.org/2005/Atom\">\n\n",
#   "<Document>\n", "\t<name>",kml_folder,"</name>\n",file = outfile, 
#   append = FALSE, sep = "") 
#   ## Icon Style Section ##
#   metadata<-read.csv(metadata, header=TRUE, as.is=TRUE, na.strings = "")
#   if (!("id" %in% colnames(metadata))){  
#     metadata$id <- metadata[,id_type]   # id_type = serial or deploy_location
#   if (id_type == "serial"){ # specific to BAEA project
#     metadata <- subset(metadata, is.na(deploy_location))  # not deployed 
#   } else {  # specific to BAEA project
#     metadata <- subset(metadata, !is.na(deploy_location))  # only deployed
#   }
#   row.names(metadata) <- NULL
#   }  
#   hi_icon_label_scale <- 0.75
#   icon_scale <- 0.7
#   ball_bg_color <- "ff333333"
#   ball_text_color <- "ffffffff"    
#   mt_icon_href <-
#     "http://earth.google.com/images/kml-icons/track-directional/track-0.png"
#   if (behavior_color == TRUE) {
#   behavior_colors_org <- CreateColorsByBehavior(output=TRUE)
#   behavior_colors <- col2kml(behavior_colors_org)
#   names(behavior_colors) <- names(behavior_colors_org)  # converts to KML color
#   behavior_colors <- gsub("#ff", "", behavior_colors)  # changes Alpha level
#   sex <- rep(c("male", "female"), each=length(behavior_colors))  
#   behavior <- rep(names(behavior_colors), 2)
#   icon_color <- rep(unname(behavior_colors),2)
#   href<-rep(c("http://maps.google.com/mapfiles/kml/shapes/placemark_square.png",
#     "http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png"), 
#     each=length(behavior_colors))
#   behaviors <- data.frame(sex, behavior, icon_color,href,stringsAsFactors=FALSE)
#   for (i in 1:nrow(behaviors)){
#   row <- behaviors[i, ]
#   sex <- row[, "sex"]
#   behavior <- row[, "behavior"]
#   icon_color <- row[ ,"icon_color"]
#   icon_href <- row[, "href"]
#   icon_m <- "http://maps.google.com/mapfiles/kml/shapes/placemark_square.png"
#   icon_f <- "http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png"
#   icon_na <- "http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png"
#   icon_na_color <- "50969696"
#   icon_na_color_hi <- "75969696"
#   icon_na_scale <- 0.3
#   cat("\t<StyleMap id=\"Point_",sex,"-",behavior,"\">\n",
#   "\t\t<Pair>\n",
#   "\t\t\t<key>normal</key>\n",
#   "\t\t\t\t<Style>\n",  
#   "\t\t\t\t\t<LabelStyle>\n",      
#   "\t\t\t\t\t<scale>",labelscale,"</scale>\n",  # to show label 
#   "\t\t\t\t\t</LabelStyle>\n",
#   "\t\t\t\t\t<IconStyle>\n",
#   "\t\t\t\t\t\t<color>FF",icon_color,"</color>\n",
#   "\t\t\t\t\t\t<scale>",icon_scale,"</scale>\n", 
#   "\t\t\t\t\t<Icon>\n",
#   "\t\t\t\t\t<href>",icon_href,"</href>\n",
#   "\t\t\t\t\t</Icon>\n",
#   "\t\t\t\t\t</IconStyle>\n",      
#   "\t\t\t\t\t<BalloonStyle>\n",
#   "\t\t\t\t\t<text>$[description]</text>\n", 
#   "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>\n",
#   "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>\n",
#   "\t\t\t\t\t</BalloonStyle>\n",
#   "\t\t\t\t</Style>\n",
#   "\t\t</Pair>\n",
#   "\t\t<Pair>\n",
#   "\t\t\t<key>highlight</key>\n",
#   "\t\t\t\t<Style>\n",
#   "\t\t\t\t\t<LabelStyle>\n",     
#   "\t\t\t\t\t<scale>",hi_icon_label_scale,"</scale>\n",  # to show label 
#   "\t\t\t\t\t</LabelStyle>\n",
#   "\t\t\t\t\t<IconStyle>\n",
#   "\t\t\t\t\t\t<color>FF",icon_color,"</color>\n", 
#   "\t\t\t\t\t\t<scale>0.8</scale>\n", 
#   "\t\t\t\t\t<Icon>\n",
#   "\t\t\t\t\t<href>",icon_href,"</href>\n",
#   "\t\t\t\t\t</Icon>\n",
#   "\t\t\t\t\t</IconStyle>\n",      
#   "\t\t\t\t\t<BalloonStyle>\n",
#   "\t\t\t\t\t<text>$[description]</text>\n", 
#   "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>", "\n",
#   "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>", "\n",
#   "\t\t\t\t\t</BalloonStyle>\n",
#   "\t\t\t\t</Style>\n", 
#   "\t\t</Pair>\n",
#   "\t</StyleMap>\n",
#   file = outfile, append = TRUE, sep = "")
#   }         
#   cat("\t<StyleMap id=\"Point_m-NA\">\n",
#   "\t\t<Pair>\n",
#   "\t\t\t<key>normal</key>\n",
#   "\t\t\t\t<Style>\n",  
#   "\t\t\t\t\t<LabelStyle>\n",      
#   "\t\t\t\t\t<scale>",labelscale,"</scale>\n",  # to show label 
#   "\t\t\t\t\t</LabelStyle>\n",
#   "\t\t\t\t\t<IconStyle>\n",
#   "\t\t\t\t\t\t<color>",icon_na_color,"</color>\n",
#   "\t\t\t\t\t\t<scale>",icon_na_scale,"</scale>\n", 
#   "\t\t\t\t\t<Icon>\n",
#   "\t\t\t\t\t<href>",icon_na,"</href>\n",
#   "\t\t\t\t\t</Icon>\n",
#   "\t\t\t\t\t</IconStyle>\n",      
#   "\t\t\t\t\t<BalloonStyle>\n",
#   "\t\t\t\t\t<text>$[description]</text>\n", 
#   "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>\n",
#   "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>\n",
#   "\t\t\t\t\t</BalloonStyle>\n",
#   "\t\t\t\t</Style>\n",
#   "\t\t</Pair>\n",
#   "\t\t<Pair>\n",
#   "\t\t\t<key>highlight</key>\n",
#   "\t\t\t\t<Style>\n",
#   "\t\t\t\t\t<LabelStyle>\n",     
#   "\t\t\t\t\t<scale>",hi_icon_label_scale,"</scale>\n",  # to show label 
#   "\t\t\t\t\t</LabelStyle>\n",
#   "\t\t\t\t\t<IconStyle>\n",
#   "\t\t\t\t\t\t<color>",icon_na_color_hi,"</color>\n", 
#   "\t\t\t\t\t\t<scale>",icon_na_scale,"</scale>\n", 
#   "\t\t\t\t\t<Icon>\n",
#   "\t\t\t\t\t<href>",icon_na,"</href>\n",
#   "\t\t\t\t\t</Icon>\n",
#   "\t\t\t\t\t</IconStyle>\n",      
#   "\t\t\t\t\t<BalloonStyle>\n",
#   "\t\t\t\t\t<text>$[description]</text>\n", 
#   "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>", "\n",
#   "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>", "\n",
#   "\t\t\t\t\t</BalloonStyle>\n",
#   "\t\t\t\t</Style>\n", 
#   "\t\t</Pair>\n",
#   "\t</StyleMap>\n",  
#   "\t<StyleMap id=\"Point_f-NA\">\n",
#   "\t\t<Pair>\n",
#   "\t\t\t<key>normal</key>\n",
#   "\t\t\t\t<Style>\n",  
#   "\t\t\t\t\t<LabelStyle>\n",      
#   "\t\t\t\t\t<scale>",labelscale,"</scale>\n",  # to show label 
#   "\t\t\t\t\t</LabelStyle>\n",
#   "\t\t\t\t\t<IconStyle>\n",
#   "\t\t\t\t\t\t<color>",icon_na_color,"</color>\n",
#   "\t\t\t\t\t\t<scale>",icon_na_scale,"</scale>\n", 
#   "\t\t\t\t\t<Icon>\n",
#   "\t\t\t\t\t<href>",icon_na,"</href>\n",
#   "\t\t\t\t\t</Icon>\n",
#   "\t\t\t\t\t</IconStyle>\n",      
#   "\t\t\t\t\t<BalloonStyle>\n",
#   "\t\t\t\t\t<text>$[description]</text>\n", 
#   "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>\n",
#   "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>\n",
#   "\t\t\t\t\t</BalloonStyle>\n",
#   "\t\t\t\t</Style>\n",
#   "\t\t</Pair>\n",
#   "\t\t<Pair>\n",
#   "\t\t\t<key>highlight</key>\n",
#   "\t\t\t\t<Style>\n",
#   "\t\t\t\t\t<LabelStyle>\n",     
#   "\t\t\t\t\t<scale>",hi_icon_label_scale,"</scale>\n",  # to show label 
#   "\t\t\t\t\t</LabelStyle>\n",
#   "\t\t\t\t\t<IconStyle>\n",
#   "\t\t\t\t\t\t<color>",icon_na_color_hi,"</color>\n", 
#   "\t\t\t\t\t\t<scale>",icon_na_scale,"</scale>\n", 
#   "\t\t\t\t\t<Icon>\n",
#   "\t\t\t\t\t<href>",icon_na,"</href>\n",
#   "\t\t\t\t\t</Icon>\n",
#   "\t\t\t\t\t</IconStyle>\n",      
#   "\t\t\t\t\t<BalloonStyle>\n",
#   "\t\t\t\t\t<text>$[description]</text>\n", 
#   "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>", "\n",
#   "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>", "\n",
#   "\t\t\t\t\t</BalloonStyle>\n",
#   "\t\t\t\t</Style>\n", 
#   "\t\t</Pair>\n",
#   "\t</StyleMap>\n",
#   "\t<StyleMap id=\"Track_male\">\n",
#   "\t\t<Pair>\n",
#   "\t\t\t<key>normal</key>\n",
#   "\t\t\t\t<Style>\n", 
#   "\t\t\t\t\t<LabelStyle>\n",
#   "\t\t\t\t\t<scale>0</scale>\n",  # to show label 
#   "\t\t\t\t\t</LabelStyle>\n",
#   "\t\t\t<IconStyle>\n",
#   "\t\t\t\t<scale>1</scale>\n",
#   "\t\t\t\t<Icon>\n",
#   "\t\t\t\t\t<href>",mt_icon_href,"</href>\n",
#   "\t\t\t\t</Icon>\n",
#   "\t\t\t</IconStyle>\n",
#   "\t\t\t\t\t<BalloonStyle>\n",
#   "\t\t\t\t\t<text>Male - Path</text>\n",       
#   "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>\n",
#   "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>\n",
#   "\t\t\t\t\t</BalloonStyle>\n",
#   "\t\t\t<LineStyle>\n",
#   "\t\t\t\t<color>99F01E14</color>\n",
#   "\t\t\t\t<width>1</width>\n",
#   "\t\t\t</LineStyle>\n",
#   "\t\t\t\t</Style>\n",
#   "\t\t</Pair>\n",
#   "\t\t<Pair>\n",
#   "\t\t\t<key>highlight</key>\n",
#   "\t\t\t\t<Style>\n", 
#   "\t\t\t\t\t<LabelStyle>\n",      
#   "\t\t\t\t\t<scale>0</scale>\n",  # to show label 
#   "\t\t\t\t\t</LabelStyle>\n",
#   "\t\t\t<IconStyle>\n",
#   "\t\t\t\t<scale>1</scale>\n",
#   "\t\t\t\t<Icon>\n",
#   "\t\t\t\t\t<href>",mt_icon_href,"</href>\n",
#   "\t\t\t\t</Icon>\n",
#   "\t\t\t</IconStyle>\n",
#   "\t\t\t\t\t<BalloonStyle>\n",
#   "\t\t\t\t\t<text>Male - Path</text>\n",       
#   "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>\n",
#   "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>\n",
#   "\t\t\t\t\t</BalloonStyle>\n",
#   "\t\t\t<LineStyle>\n",
#   "\t\t\t\t<color>99F01E14</color>\n",  # male path color 
#   "\t\t\t\t<width>1</width>\n",
#   "\t\t\t</LineStyle>\n",
#   "\t\t\t\t</Style>\n",
#   "\t\t</Pair>\n",
#   "\t</StyleMap>\n",      
#   ## Style Map for Track_f ##  
#   "\t<StyleMap id=\"Track_female\">\n",
#   "\t\t<Pair>\n",
#   "\t\t\t<key>normal</key>\n",
#   "\t\t\t\t<Style>\n", 
#   "\t\t\t\t\t<LabelStyle>\n",      
#   "\t\t\t\t\t<scale>0</scale>\n",  # to show label 
#   "\t\t\t\t\t</LabelStyle>\n",
#   "\t\t\t<IconStyle>\n",
#   "\t\t\t\t<scale>1</scale>\n",
#   "\t\t\t\t<Icon>\n",
#   "\t\t\t\t\t<href>",mt_icon_href,"</href>\n",
#   "\t\t\t\t</Icon>\n",
#   "\t\t\t</IconStyle>\n",
#   "\t\t\t\t\t<BalloonStyle>\n",
#   "\t\t\t\t\t<text>Female - Path</text>\n",       
#   "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>\n",
#   "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>\n",
#   "\t\t\t\t\t</BalloonStyle>\n",
#   "\t\t\t<LineStyle>\n",
#   "\t\t\t\t<color>99A078F0</color>\n",  # female path color
#   "\t\t\t\t<width>1</width>\n",
#   "\t\t\t</LineStyle>\n",
#   "\t\t\t\t</Style>\n",
#   "\t\t</Pair>\n",
#   "\t\t<Pair>\n",
#   "\t\t\t<key>highlight</key>\n",
#   "\t\t\t\t<Style>\n", 
#   "\t\t\t\t\t<LabelStyle>\n",      
#   "\t\t\t\t\t<scale>0</scale>\n",  # to show label 
#   "\t\t\t\t\t</LabelStyle>\n",
#   "\t\t\t<IconStyle>\n",
#   "\t\t\t\t<scale>1</scale>\n",
#   "\t\t\t\t<Icon>\n",
#   "\t\t\t\t\t<href>",mt_icon_href,"</href>\n",
#   "\t\t\t\t</Icon>\n",
#   "\t\t\t</IconStyle>\n",
#   "\t\t\t\t\t<BalloonStyle>\n",
#   "\t\t\t\t\t<text>Female - Path</text>\n",       
#   "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>\n",
#   "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>\n",
#   "\t\t\t\t\t</BalloonStyle>\n",
#   "\t\t\t<LineStyle>\n",
#   "\t\t\t\t<color>99A078F0</color>\n",  # female path color
#   "\t\t\t\t<width>1</width>\n",
#   "\t\t\t</LineStyle>\n",
#   "\t\t\t\t</Style>\n",
#   "\t\t</Pair>\n",
#   "\t</StyleMap>\n",   
#   file = outfile, append = TRUE, sep = "") 
#   } else {
#   df_ids <- unique(df$id)
#   metadata <- subset(metadata, id %in% df_ids) 
#   icon_href <- "http://maps.google.com/mapfiles/kml/shapes/placemark_square.png"
#   for (i in 1:nrow(metadata)) {
#   row <- metadata[i, ]
#   id <- row[, "id"]
#   icon_color <- col2kml(row[ ,"icon_color"])
#   icon_color <- gsub("#ff", "", icon_color)
#   icon_scale <- 0.7
#   cat("\t<StyleMap id=\"Point_",id,"\">\n",
#   "\t\t<Pair>\n",
#   "\t\t\t<key>normal</key>\n",
#   "\t\t\t\t<Style>\n",  
#   "\t\t\t\t\t<LabelStyle>\n",      
#   "\t\t\t\t\t<scale>",labelscale,"</scale>\n",  # to show label 
#   "\t\t\t\t\t</LabelStyle>\n",
#   "\t\t\t\t\t<IconStyle>\n",
#   "\t\t\t\t\t\t<color>FF",icon_color,"</color>\n",
#   "\t\t\t\t\t\t<scale>",icon_scale,"</scale>\n", 
#   "\t\t\t\t\t<Icon>\n",
#   "\t\t\t\t\t<href>",icon_href,"</href>\n",
#   "\t\t\t\t\t</Icon>\n",
#   "\t\t\t\t\t</IconStyle>\n",      
#   "\t\t\t\t\t<BalloonStyle>\n",
#   "\t\t\t\t\t<text>$[description]</text>\n", 
#   "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>\n",
#   "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>\n",
#   "\t\t\t\t\t</BalloonStyle>\n",
#   "\t\t\t\t</Style>\n",
#   "\t\t</Pair>\n",        
#   "\t\t<Pair>\n",
#   "\t\t\t<key>highlight</key>\n",
#   "\t\t\t\t<Style>\n",
#   "\t\t\t\t\t<LabelStyle>\n",     
#   "\t\t\t\t\t<scale>",hi_icon_label_scale,"</scale>\n",  # to show label 
#   "\t\t\t\t\t</LabelStyle>\n",
#   "\t\t\t\t\t<IconStyle>\n",
#   "\t\t\t\t\t\t<color>FF",icon_color,"</color>\n", 
#   "\t\t\t\t\t\t<scale>0.8</scale>\n", 
#   "\t\t\t\t\t<Icon>\n",
#   "\t\t\t\t\t<href>",icon_href,"</href>\n",
#   "\t\t\t\t\t</Icon>\n",
#   "\t\t\t\t\t</IconStyle>\n",      
#   "\t\t\t\t\t<BalloonStyle>\n",
#   "\t\t\t\t\t<text>$[description]</text>\n", 
#   "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>", "\n",
#   "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>", "\n",
#   "\t\t\t\t\t</BalloonStyle>\n",
#   "\t\t\t\t</Style>\n", 
#   "\t\t</Pair>\n",
#   "\t</StyleMap>\n",
#   ## Style Map for Track ##
#   "\t<StyleMap id=\"Track_",id,"\">\n",
#   "\t\t<Pair>\n",
#   "\t\t\t<key>normal</key>\n",
#   "\t\t\t\t<Style>\n", 
#   "\t\t\t\t\t<LabelStyle>\n",      
#   "\t\t\t\t\t<scale>0</scale>\n",  # to show label 
#   "\t\t\t\t\t</LabelStyle>\n",
#   "\t\t\t<IconStyle>\n",
#   "\t\t\t\t<scale>1</scale>\n",
#   "\t\t\t\t<Icon>\n",
#   "\t\t\t\t\t<href>",mt_icon_href,"</href>\n",
#   "\t\t\t\t</Icon>\n",
#   "\t\t\t</IconStyle>\n",
#   "\t\t\t\t\t<BalloonStyle>\n",
#   "\t\t\t\t\t<text>",id,"- Path</text>\n",       
#   "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>\n",
#   "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>\n",
#   "\t\t\t\t\t</BalloonStyle>\n",
#   "\t\t\t<LineStyle>\n",
#   "\t\t\t\t<color>dd",icon_color,"</color>\n",
#   "\t\t\t\t<width>1</width>\n",
#   "\t\t\t</LineStyle>\n",
#   "\t\t\t\t</Style>\n",
#   "\t\t</Pair>\n",    
#   "\t\t<Pair>\n",
#   "\t\t\t<key>highlight</key>\n",
#   "\t\t\t\t<Style>\n", 
#   "\t\t\t\t\t<LabelStyle>\n",      
#   "\t\t\t\t\t<scale>0</scale>\n",  # to show label 
#   "\t\t\t\t\t</LabelStyle>\n",
#   "\t\t\t<IconStyle>\n",
#   "\t\t\t\t<scale>1</scale>\n",
#   "\t\t\t\t<Icon>\n",
#   "\t\t\t\t\t<href>",mt_icon_href,"</href>\n",
#   "\t\t\t\t</Icon>\n",
#   "\t\t\t</IconStyle>\n",
#   "\t\t\t\t\t<BalloonStyle>\n",
#   "\t\t\t\t\t<text>",id," - Path</text>\n",       
#   "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>\n",
#   "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>\n",
#   "\t\t\t\t\t</BalloonStyle>\n",
#   "\t\t\t<LineStyle>\n",
#   "\t\t\t\t<color>ee",icon_color,"</color>\n",
#   "\t\t\t\t<width>1</width>\n",
#   "\t\t\t</LineStyle>\n",
#   "\t\t\t\t</Style>\n",
#   "\t\t</Pair>\n",
#   "\t</StyleMap>\n",              
#   file = outfile, append = TRUE, sep = "")
#   }  
#   }  # end of icon loops 
#   ids <- as.character(unique(df$id))  # as.character removes factor levels
#   for (i in ids){
#   sv = df$id %in% i
#   unique.id <- as.character(unique(df$id[sv]))    
#   sex <- as.character(unique(df$sex[sv]))
#   cat("<Folder>\n","<name>",unique.id,"</name>\n","<open>0</open>\n",
#   file = outfile, append = TRUE, sep = "")
#   cat("\t<Folder>\n","\t<name>",unique.id," - Locations</name>\n",
#     "\t<open>0</open>\n", file = outfile, append = TRUE, sep = "")
#   locs <- subset(df, id==unique.id)
#   for (i in 1:nrow(locs)){
#   loc <- locs[i,]
#   Xs <- loc[, "long"]
#   Ys <- loc[, "lat"]
#   Zs <- loc[, "alt"]
#   SPs <- loc[, "speed"]
#   BHs <- loc[, "behavior"]
#   IDs <- unique.id 
#   SDs <- substring(loc$datetime, 1,10) #start date
#   STs <- substring(loc$datetime, 12,16) #start time 
#   EDs <- substring(loc$datetimeend, 1,10) #end date
#   ETs <- substring(loc$datetimeend, 12,19) #end time 
#   DAs <- strftime(loc[, "datetimebegin"], dateformat)
#   TMs <- strftime(loc[, "datetimebegin"], timeformat)      
#   plnms <- strftime(loc[, "datetimebegin"], datetimeformat)
#   PlacemarkPoint(plnms, Xs, Ys, Zs, SPs, BHs, IDs, SDs, STs, EDs, ETs, DAs, TMs)     
#   }      
#   cat("\t</Folder>\n", file = outfile, append = TRUE, sep = "")
#   locs$Ts <- "T"
#   locs$Zs <- "Z"
#   locs$datetimedate <-substring(locs$datetime, 1,10) #start date
#   locs$datetimetime <- substring(locs$datetime, 12,16) #start time
#   whens <- locs[, c("datetimedate","Ts","datetimetime", "Zs")]
#   sgmts <- locs[, c("long","lat","alt")]
#   unique.id <- unique(locs$id)
#   bloc2 <- NULL
#   if (path == TRUE) {
#   if (behavior_color == TRUE) {
#   bloc2 <- c(bloc2, paste(
#   "\t<Placemark>\n",
#   "\t\t<name>",unique.id," - Path</name>\n",
#   "\t\t<styleUrl>#Track_",sex,"</styleUrl>\n",
#   "\t\t<gx:balloonVisibility>0</gx:balloonVisibility>\n",
#   "\t\t<gx:Track>\n",
#   "\t\t<altitudeMode>",alt_mode,"</altitudeMode>\n",
#   paste(paste("\t\t\t\t\t<when>", apply(whens, 1, paste, collapse = ""), 
#     sep = ""), "</when>", collapse = "\n"),"\n", 
#   paste(paste("\t\t\t\t\t<gx:coord>", apply(sgmts, 1, paste, collapse = " "), 
#     sep = ""),"</gx:coord>", collapse = "\n"),"\n", 
#   "\t\t</gx:Track>\n",
#   "\t</Placemark>\n",
#   "\t</Folder>\n",                  
#   sep = ""))  
#   } else {
#   bloc2 <- c(bloc2, paste(
#   "\t<Placemark>\n",
#   "\t\t<name>",unique.id," - Path</name>\n",
#   "\t\t<styleUrl>#Track_",unique.id,"</styleUrl>\n",
#   "\t\t<gx:balloonVisibility>0</gx:balloonVisibility>\n",
#   "\t\t<gx:Track>\n",
#   "\t\t<altitudeMode>",alt_mode,"</altitudeMode>\n",
#   paste(paste("\t\t\t\t\t<when>", apply(whens, 1, paste, collapse = ""), 
#     sep = ""), "</when>", collapse = "\n"),"\n", 
#   paste(paste("\t\t\t\t\t<gx:coord>", apply(sgmts, 1, paste, collapse = " "), 
#     sep = ""),"</gx:coord>", collapse = "\n"),"\n", 
#   "\t\t</gx:Track>\n",
#   "\t</Placemark>\n",
#   "\t</Folder>\n",                  
#   sep = ""))  
#   }
#   cat(bloc2, file = outfile, append = TRUE)                                      
#   }
#   }
#   cat("</Document>\n</kml>", file = outfile, append = TRUE)
# }
 

# ExportKMLWindProjects Function -----------------------------------------------

###  Exports KML of wind project locations
###  Usage: ExportKMLWindProjects(df, labelscale, outfile, metadata)
###  Arguments: df = the location of the .csv file
###             labelscale = adjusts the size of the Google Earth location 
###               point labels. Default is 0, which hides the labels. To show 
###               labels, change to a value between 0.7-1.
###             outfile = location of output KML file
###             metadata = location of metadata .csv file. Metadata file 
###               must contain columns for project ID, hexadecimal colors, 
###               and additional icon data.
###  Returns: creates KML     
###  Notes: defaults are specific to my file directories and locations
###  Blake Massey
###  2014.05.05

ExportKMLWindProjects <- function(df = file.path("C:/Work/R/Data", 
                                    "Turbines/Maine Wind Projects.csv"), 
                                  labelscale = 0, 
                                  outfile = file.path("C:/Users/Blake/Desktop",
                                    "Maine Wind Projects.kml"),
                                  metadata = file.path("C:/Work/R/Data/",
                                    "Turbines/Project_records.csv")) {
  plmark.pnt.fun <- function(project_name, developer, status, model, 
                             num_turbine, total_mw, hub_height, rotor_diam, 
                             rotor_top, rotor_btm, X, Y) {
    cat("\t<Placemark>\n",
        "\t\t<name>",project_name,"</name>\n",
        "<visibility>1</visibility>\n",
        "\t\t\t<Snippet></Snippet>", "\n",
        "\t\t\t\t<description>\n",
        "\t\t\t\t\tDeveloper: ",developer, "\n",
        "\t\t\t\t\tStatus: ",status, "\n",
        "\t\t\t\t\tTurbine Models: ",model, "\n",
        "\t\t\t\t\tNumber of Turbines: ",num_turbine, "\n",
        "\t\t\t\t\tTotal MW: ",total_mw, "\n",
        "\t\t\t\t\tHub Height: ",hub_height, " meters\n",
        "\t\t\t\t\tRotor Diameter: ",rotor_diam," meters\n",        
        "\t\t\t\t\tRotor Sweep Top: ",rotor_top, " meters\n",          
        "\t\t\t\t\tRotor Sweep Bottom: ",rotor_btm, " meters\n",  
        "\t\t\t\t\tLongitude: ",X, "\n",
        "\t\t\t\t\tLatitude: ",Y, "\n", 
        "\t\t\t\t</description>\n",
        "\t\t\t<styleUrl>#Style_",status,"</styleUrl>\n",
        "\t\t\t<Point>\n",
        "\t\t\t\t<extrude>0</extrude>\n",
        "\t\t\t\t<altitudeMode>clampedtoGround</altitudeMode>\n",
        "\t\t\t\t\t<coordinates>",X,",",Y,"0 </coordinates>\n",
        "\t\t\t</Point>\n",        
        "\t</Placemark>\n",
        file = outfile, append = TRUE, sep = "")  
  }
  if (file.exists(outfile)) file.remove(outfile)  # overwrite KML file
  writeLines(noquote(c("Writing: ",outfile)))  # print "Writing ... .kml"
  ## Title Section ##
  cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n",      
      "<kml xmlns=\"http://www.opengis.net/kml/2.2\"\n", 
      "xmlns:kml=\"http://www.opengis.net/kml/2.2\"\n",
      "xmlns:gx=\"http://www.google.com/kml/ext/2.2\"\n",
      "xmlns:atom=\"http://www.w3.org/2005/Atom\">\n\n",
      "<Document>\n",
      "\t<name>Maine Wind Projects</name>\n",
      file = outfile, append = FALSE, sep = "") 
  ## Icon Style Section ##
  metadata<-read.csv(metadata, header=TRUE, as.is=TRUE)
  for (i in 1:nrow(metadata)){
    row <- i
    sn<-metadata[row,]
    status <- sn[, "status"]
    icon_label_scale <- sn[, "icon_label_scale"]
    icon_color <-sn[, "icon_color"]
    icon_scale <-sn[, "icon_scale"]
    icon_href <-as.character(sn[, "icon_href"])
    icon_ball_bg_color <- sn[, "icon_ball_bg_color"]
    icon_ball_text_color <- sn[, "icon_ball_text_color"]    
    hi_icon_label_scale <- sn[, "hi_icon_label_scale"]
    hi_icon_color <- sn[, "hi_icon_color"]
    hi_icon_scale <- sn[, "hi_icon_scale"]
    hi_icon_href <- sn[, "hi_icon_href"]
    hi_icon_ball_bg_color <- sn[, "hi_icon_ball_bg_color"]
    hi_icon_ball_text_color <- sn[, "hi_icon_ball_text_color"] 
    cat("\t<StyleMap id=\"Style_",status,"\">\n",
        "\t\t<Pair>\n",
        "\t\t\t<key>normal</key>\n",
        "\t\t\t\t<Style>\n",  
        "\t\t\t\t\t<LabelStyle>\n",      
        "\t\t\t\t\t<scale>", icon_label_scale, "</scale>\n",  # to show label 
        "\t\t\t\t\t</LabelStyle>\n",
        "\t\t\t\t\t<IconStyle>\n",
        "\t\t\t\t\t\t<color>", icon_color, "</color>\n",
        "\t\t\t\t\t\t<scale>", icon_scale, "</scale>\n", 
        "\t\t\t\t\t<Icon>\n",
        "\t\t\t\t\t<href>", icon_href, "</href>\n",
        "\t\t\t\t\t</Icon>\n",
        "\t\t\t\t\t</IconStyle>\n",      
        "\t\t\t\t\t<BalloonStyle>\n",
        "\t\t\t\t\t<text>$[description]</text>\n", 
        "\t\t\t\t\t\t<bgColor>", icon_ball_bg_color, "</bgColor>\n", 
        "\t\t\t\t\t\t<textColor>", icon_ball_text_color, "</textColor>\n", 
        "\t\t\t\t\t</BalloonStyle>\n",
        "\t\t\t\t</Style>\n",
        "\t\t</Pair>\n",
        "\t\t<Pair>\n",
        "\t\t\t<key>highlight</key>\n",
        "\t\t\t\t<Style>\n",
        "\t\t\t\t\t<LabelStyle>\n",     
        "\t\t\t\t\t<scale>", hi_icon_label_scale, "</scale>\n",  # to show label 
        "\t\t\t\t\t</LabelStyle>\n",
        "\t\t\t\t\t<IconStyle>\n",
        "\t\t\t\t\t\t<color>", hi_icon_color, "</color>\n", 
        "\t\t\t\t\t\t<scale>", hi_icon_scale, "</scale>\n", 
        "\t\t\t\t\t<Icon>\n",
        "\t\t\t\t\t<href>", hi_icon_href, "</href>\n",
        "\t\t\t\t\t</Icon>\n",
        "\t\t\t\t\t</IconStyle>\n",      
        "\t\t\t\t\t<BalloonStyle>\n",
        "\t\t\t\t\t<text>$[description]</text>\n", 
        "\t\t\t\t\t\t<bgColor>", hi_icon_ball_bg_color, "</bgColor>", "\n", 
        "\t\t\t\t\t\t<textColor>", hi_icon_ball_text_color, "</textColor>","\n",
        "\t\t\t\t\t</BalloonStyle>\n",
        "\t\t\t\t</Style>\n", 
        "\t\t</Pair>\n",
        "\t</StyleMap>\n",
        file = outfile, append = TRUE, sep = "")
  }
  ## For Loop Section ##  
  df <- read.csv(df, header=TRUE, as.is=TRUE)
  df <- subset(df, show == 1)  # only include projects with 1 in "show" column
  status <- unique(df$status)
  for (i in status){
    sv = df$status %in% i
    status <- as.character(unique(df$status[sv]))
    cat("<Folder>\n",
        "<name>", status, "</name>\n",
        "<open>0</open>\n",
        file = outfile, append = TRUE, sep = "")
    projects <- subset(df, status==i)
    for (i in 1:nrow(projects)){
      row <- i
      project <- projects[row,]
      project_names <- project[, "project"]
      developers <- project[, "developer"]
      statuses <- project[, "status"]
      models <- project [, "turb_model"]
      num_turbines <- project [, "num_turbine"]       
      total_mws <- project [, "total_mw"] 
      hub_heights <- project[, "hub_ht_m"]
      rotor_diams <- project [, "rotor_diam"]
      rotor_tops <- project [, "rotor_top"]
      rotor_btms <- project [, "rotor_btm"]
      Xs <- project [, "long"]
      Ys <- project [, "lat"]
      plmark.pnt.fun(project_names, developers, statuses, models, num_turbines, 
                     total_mws, hub_heights, rotor_diams, rotor_tops, 
                     rotor_btms, Xs, Ys)      
    }                                         
    cat("</Folder>\n", file = outfile, append = TRUE, sep = "")
  }                                         
  cat("</Document>\n</kml>", file = outfile, append = TRUE)
}

# ExportKMLWindTurbines Function -----------------------------------------------

###  Exports KML of wind turbine locations
###  Usage: ExportKMLWindTurbines(df, labelscale, outfile, metadata)
###  Arguments: df = the location of the .csv file
###             labelscale = adjusts the size of the Google Earth location 
###               point labels. Default is 0, which hides the labels. To show 
###               labels, change to a value between 0.7-1.
###             outfile = location of output KML file
###             metadata = location of metadata .csv file. Metadata file 
###               must contain columns for project ID, hexadecimal colors, 
###               and additional icon data.
###  Returns: creates KML         
###  Notes: defaults are specific to my file directories and locations
###  Blake Massey
###  2014.05.05

ExportKMLWindTurbines <- function (df = file.path("C:/Work/R/Data",
                                     "Turbines/Maine Wind Turbines.csv"), 
                                   labelscale = 0, 
                                   outfile = file.path("C:/Users/Blake/Desktop",
                                     "Maine Wind Turbines.kml"),
                                   metadata = file.path("C:/Work/R/Data",
                                     "Turbines/Turbine_records.csv")) {
  plmark.pnt.fun <- function(plname, X, Y, pad, model, hub, rotor, sweep_top, 
                             sweep_btm, XP, XN, YP, YN) {
    cat("\t<Placemark>\n",
        "\t\t<name>Turbine ",plname,"</name>\n",
        "<visibility>1</visibility>\n",
        "\t\t\t<Snippet></Snippet>", "\n",
        "\t\t\t\t<description>\n",
        "\t\t\t\t\tTurbine ",plname,"\n",
        "\t\t\t\t\tLongitude: ",X, "\n",
        "\t\t\t\t\tLatitude: ",Y, "\n",        
        "\t\t\t\t\tPad Elevation: ",pad, " meters\n", 
        "\t\t\t\t\tModel: ",model, "\n",
        "\t\t\t\t\tHub Height: ",hub, " meters\n",
        "\t\t\t\t\tRotor Diameter: ",rotor," meters\n",        
        "\t\t\t\t\tRotor Sweep Top: ",sweep_top, " meters\n",          
        "\t\t\t\t\tRotor Sweep Bottom: ",sweep_btm, " meters\n",  
        "\t\t\t\t</description>\n",
        "\t\t\t<styleUrl>#Style_",model,"</styleUrl>\n",
        "\t\t\t<MultiGeometry>\n",
        "\t\t\t<Model>\n",
        "\t\t\t\t<range>500</range>\n",
        "\t\t\t\t<altitudeMode>clampedToGround</altitudeMode>\n",
        "\t\t\t\t<Location>\n",
        "\t\t\t\t\t<longitude>",X,"</longitude>\n",
        "\t\t\t\t\t<latitude>",Y,"</latitude>\n",
        "\t\t\t\t\t<altitude>0</altitude>\n",
        "\t\t\t\t</Location>\n",
        "\t\t\t\t<Orientation>\n",
        "\t\t\t\t\t<heading>120</heading>\n",
        "\t\t\t\t\t<tilt>0</tilt>\n",
        "\t\t\t\t\t<roll>0</roll>\n",
        "\t\t\t\t</Orientation>\n",
        "\t\t\t\t<Scale>\n",
        "\t\t\t\t\t<x>1</x>\n",
        "\t\t\t\t\t<y>1</y>\n",
        "\t\t\t\t\t<z>1</z>\n",
        "\t\t\t\t</Scale>\n",
        "\t\t\t\t<Link>\n",
        "\t\t\t\t\t<href>C:/ArcGIS/Data/Wind/Turbine Models/Collada Models/",
          model,".dae</href>\n",
        "\t\t\t\t</Link>\n",
        "\t\t\t</Model>\n",
        "\t\t\t<Polygon>\n",
        "\t\t\t\t<tesselate>1</tesselate>\n",
        "\t\t\t\t\t<outerBoundaryIs>\n",
        "\t\t\t\t\t\t<LinearRing>\n",
        "\t\t\t\t\t\t\t<coordinates>\n",
        "\t\t\t\t\t\t\t\t",XN,",",YP,",0 ",XN,",",YN,",0 ", XP,",",YN,",0 ",XP,
          ",",YP,",0 ",XN,",",YP,",0 \n",
        "\t\t\t\t\t\t\t</coordinates>\n",        
        "\t\t\t\t\t\t</LinearRing>\n",        
        "\t\t\t\t\t</outerBoundaryIs>\n",
        "\t\t\t</Polygon>\n",
        "\t\t\t</MultiGeometry>\n",        
        "\t</Placemark>\n",
        file = outfile, append = TRUE, sep = "")  
  }
  if (file.exists(outfile)) file.remove(outfile)  # overwrites KML file
  writeLines(noquote(c("Writing: ",outfile)))  # print "Writing ... .kml"
  ## Title Section ##
  cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n",      
      "<kml xmlns=\"http://www.opengis.net/kml/2.2\"\n", 
      "xmlns:kml=\"http://www.opengis.net/kml/2.2\"\n",
      "xmlns:gx=\"http://www.google.com/kml/ext/2.2\"\n",
      "xmlns:atom=\"http://www.w3.org/2005/Atom\">\n\n",
      "<Document>\n",
      "\t<name>Maine Wind Turbines</name>\n",
      file = outfile, append = FALSE, sep = "") 
  ## Icon Style Section ##
  metadata <- read.csv(metadata, header=TRUE, as.is=TRUE)
  for (i in 1:nrow(metadata)){
    row <- i
    sn <- metadata[row,]
    id <- sn[, "id"]
    model <- sn[, "model"]
    icon_label_scale <- sn[, "icon_label_scale"]
    icon_ball_bg_color <- sn[, "icon_ball_bg_color"]
    icon_ball_text_color <- sn[, "icon_ball_text_color"]    
    hi_icon_label_scale <- sn[, "hi_icon_label_scale"]
    hi_icon_ball_bg_color <- sn[, "hi_icon_ball_bg_color"]
    hi_icon_ball_text_color <- sn[, "hi_icon_ball_text_color"] 
    cat("\t<StyleMap id=\"Style_", model, "\">\n",
        "\t\t<Pair>\n",
        "\t\t\t<key>normal</key>\n",
        "\t\t\t\t<Style>\n",  
        "\t\t\t\t\t<LabelStyle>\n",      
        "\t\t\t\t\t<scale>", icon_label_scale, "</scale>\n",  # to show label
        "\t\t\t\t\t</LabelStyle>\n",
        "\t\t\t\t\t<IconStyle>\n",
        "\t\t\t\t\t\t<color>ff00aaff</color>\n",
        "\t\t\t\t\t\t<scale>1</scale>\n",
        "\t\t\t\t\t\t<Icon>\n",
        "\t\t\t\t\t\t<href>",
          "http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n",
        "\t\t\t\t\t\t</Icon>\n",
        "\t\t\t\t\t</IconStyle>\n",
        "\t\t\t\t\t<BalloonStyle>\n",
        "\t\t\t\t\t<text>$[description]</text>\n",
        "\t\t\t\t\t\t<bgColor>", icon_ball_bg_color, "</bgColor>\n",
        "\t\t\t\t\t\t<textColor>", icon_ball_text_color, "</textColor>\n",
        "\t\t\t\t\t</BalloonStyle>\n",
        "\t\t\t\t\t<LineStyle>\n",
        "\t\t\t\t\t\t<color>01ffffff</color>\n",
        "\t\t\t\t\t</LineStyle>\n",        
        "\t\t\t\t\t<PolyStyle>\n",
        "\t\t\t\t\t\t<color>01ffffff</color>\n",
        "\t\t\t\t\t</PolyStyle>\n",         
        "\t\t\t\t</Style>\n",
        "\t\t</Pair>\n",
        "\t\t<Pair>\n",
        "\t\t\t<key>highlight</key>\n",
        "\t\t\t\t<Style>\n",
        "\t\t\t\t\t<LabelStyle>\n",     
        "\t\t\t\t\t<scale>", hi_icon_label_scale, "</scale>\n", # to show label
        "\t\t\t\t\t</LabelStyle>\n",
        "\t\t\t\t\t<IconStyle>\n",
        "\t\t\t\t\t\t<color>ff00aaff</color>\n",
        "\t\t\t\t\t\t<scale>1</scale>\n",
        "\t\t\t\t\t\t<Icon>\n",
        "\t\t\t\t\t\t<href>",
          "http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n",
        "\t\t\t\t\t\t</Icon>\n",
        "\t\t\t\t\t</IconStyle>\n",
        "\t\t\t\t\t<BalloonStyle>\n",
        "\t\t\t\t\t<text>$[description]</text>\n", 
        "\t\t\t\t\t\t<bgColor>", hi_icon_ball_bg_color, "</bgColor>", "\n", 
        "\t\t\t\t\t\t<textColor>", hi_icon_ball_text_color,"</textColor>", "\n",
        "\t\t\t\t\t</BalloonStyle>\n",
        "\t\t\t\t\t<LineStyle>\n",
        "\t\t\t\t\t\t<color>01ffffff</color>\n",
        "\t\t\t\t\t</LineStyle>\n",        
        "\t\t\t\t\t<PolyStyle>\n",
        "\t\t\t\t\t\t<color>01ffffff</color>\n",
        "\t\t\t\t\t</PolyStyle>\n",
        "\t\t\t\t</Style>\n", 
        "\t\t</Pair>\n",
        "\t</StyleMap>\n",
        file = outfile, append = TRUE, sep = "") 
  }  
  ## For Loop Section ##  
  df <- read.csv(df, header=TRUE, as.is=TRUE)
  developer <- unique(df$developer)
  for (i in developer){
    sv = df$developer %in% i
    developer <- as.character(unique(df$developer[sv]))
    cat("<Folder>\n",
        "<name>", developer, "</name>\n",
        "<open>0</open>\n",
        file = outfile, append = TRUE, sep = "")
    projects <- subset(df, developer==i)
    project <- unique(projects$project)
    for (i in project){
      sv = projects$project %in% i
      project_name <- as.character(unique(projects$project[sv]))
      cat("<Folder>\n",
          "<name>", project_name, " Project</name>\n",
          "<open>0</open>\n",
          file = outfile, append = TRUE, sep = "")
      turbines <- subset(projects, project==i)
      for (i in 1:nrow(turbines)){
        row <- i
        turbine <- turbines[row,]
        plnms <- turbine[, "turbine_id"]
        Xs <- turbine[, "long"]
        XPs <- Xs + .001
        XNs <- Xs - .001
        Ys <- turbine[, "lat"]
        YPs <- Ys +.001
        YNs <- Ys -.001
        pads <- turbine[, "pad_elev_m"]
        models <- turbine[, "model"]
        hubs <- turbine[, "hub_ht_m"]
        rotors <- turbine [, "rotor_diam"]
        sweep_tops <- turbine [, "sweep_top"]
        sweep_btms <- turbine [, "sweep_btm"]
        plmark.pnt.fun(plnms, Xs, Ys, pads, models, hubs, rotors, sweep_tops, 
          sweep_btms, XPs, XNs, YPs, YNs)      
      }                                         
      cat("</Folder>\n", file = outfile, append = TRUE, sep = "")
    }                                         
    cat("</Folder>\n", file = outfile, append = TRUE, sep = "")
  }
  cat("</Document>\n</kml>", file = outfile, append = TRUE)
}

# ExportShapefileDF Function ###################################################

###  Exports shapefile of a dataframe's location data
###  Usage: ExportShapefileDF(df, layer, folder, overwrite_layer)
###  Arguments: df = dataframe with locations
###             layer = name of shapefile layer
###             folder = location for shapefile files
###             overwrite_layer = default FALSE, use with caution 
###  Returns: creates shapefile
###  Notes: requires "lat" and "long" columns
###  Blake Massey
###  2014.05.05

ExportShapefileDF<-function(df = df, 
                            layer = "baea", 
                            folder = "C:/Work/R/Data/Output",
                            overwrite_layer = FALSE){
  suppressPackageStartupMessages(require(rgdal))
  suppressPackageStartupMessages(require(maptools))
  xy<-(cbind(df$long, df$lat))
  classes <- as.character(sapply(df, class))
  colClasses <- which(classes == "c(\"POSIXct\", \"POSIXt\")")
  df[, colClasses] <- sapply(df[, colClasses], as.character)
  df_sp <- SpatialPointsDataFrame(xy, df, coords.nrs = numeric(0),
    proj4string = CRS("+proj=longlat +datum=WGS84"), match.ID = TRUE, 
    bbox = NULL)  # fails if there are spaces in CRS 
  writeLines(noquote(paste("Writing: ", folder, "/", layer, ".shp", sep="")))
  writeOGR(df_sp, folder, layer, driver = "ESRI Shapefile",
    overwrite_layer = overwrite_layer)
}

# ImportLandscapeRasterStack Function ------------------------------------------

###  Imports raster layers and combines them into a RasterStack       
###  Usage: ImportLandscapeRasterStack()
###  Arguments: none
###  Returns: The original locations df with additional landscape data columns.             
###  Notes: This function is specific to my directory and datafiles
###  Blake Massey
###  2014.09.09

ImportLandscapeRasterStack <- function(){
  suppressPackageStartupMessages(require(raster))
  hydro_50mc <- raster("C:/ArcGIS/Data/Hydrology/NHD_rasters/hydro_50mc.tif") 
  hydro_dir_50mc <- raster(file.path("C:/ArcGIS/Data/Hydrology/NHD_rasters",
    "hydro_dir_50mc.tif"))
  hydro_dist_50mc <- raster(file.path("C:/ArcGIS/Data/Hydrology/NHD_rasters",
    "hydro_dist_50mc.tif"))
  elev_50mc <- raster("C:/ArcGIS/Data/Elevation/Elevation/elev_50mc.tif")
  lc_50mc <- raster("C:/ArcGIS/Data/Landcover/Landcover/lc_50mc.tif")
  maine_50mc <- raster("C:/ArcGIS/Data/BlankRaster/maine_50mc.tif")
  raster_stack <- stack(elev_50mc, hydro_50mc, hydro_dir_50mc, hydro_dist_50mc,
    lc_50mc, maine_50mc)
  cat("The following layers were used to create the RasterStack:",sep="\n")
  names <- rev(sort(names(raster_stack), decreasing=TRUE))
  layers <- sapply(names, function(i) paste(" ", i))
  cat(layers, sep="\n")
  return(raster_stack)
}

# PrintRasterNames Function ----------------------------------------------------

###  Prints the position and name of the rasters in a RasterStack or RasterBrick       
###  Usage: PrinttRasterNames(raster)
###  Arguments: raster = RasterStack or RasterBrick
###  Returns: Prints a list of positions and names             
###  Notes: This function
###  Blake Massey
###  2014.09.10

PrintRasterNames <- function(raster){
  raster <- raster
  for (i in 1:nlayers(raster)){
  results <- paste(i,") ", names(raster[[i]]), sep="")
  cat(results, sep="\n")
  }
}
