# CreateColorsByMetadata Function ----------------------------------------------

###  Creates and/or displays dataframe of IDs and their associated colors 
###  Usage: CreateColorsByMetadata(file = ,output, display)
###  Arguments: file = CSV file with columns for "id" and "icon_color"
###             id = column name for unique identifier, default is "id"
###  Returns: df with ID names and hexidecimal colors          
###  Notes: Used in several other functions
###  Blake Massey
###  2014.10.21

CreateColorsByMetadata <- function(file,
                                   metadata_id = "id"){
  source('C:/Work/R/Functions/gen.R')
  source('C:/Work/R/Functions/gis.R')
  metadata <- read.csv(file, header=TRUE, as.is=TRUE, na.strings = "")  
  metadata$id <- metadata[ ,metadata_id] 
  id_colors <- metadata$icon_color   
  names(id_colors) <- metadata$id 
  return(id_colors)
}

# CreateColorsByVar Function ---------------------------------------------------

###  Creates dataframe of a variable and the associated colors 
###  Usage: CreateColorsByVar(df, by, b_pal, r_pal_num, pal, ouput, display)
###  Arguments: df = dataframe that contains the variable column 
###             by = column name of variable used to create colors. Creates a   
###               palette based on the unique number of factors in the variable 
###               and assigned colors based on the 'pal', 'b_pal', or 'r_pal' 
###               parameters.
###             pal = name of color palette funtions (e.g., rainbow, 
###               heat.colors, terrain.colors, topo.colors, or cm.colors)
###               used to create colors. This parameter has priority over the 
###               other color palette parameters. Default is NULL
###             r_pal = Specifc number of 'R_pal' color palette from the PlotKML
###               Package (e.g., 1 = R_pal[[1]]). This parameter has priority 
###               over the 'b_pal' parameter for setting the colors. Default is 
###               NULL.
###             b_pal = color palette name from 'RColorBrewer' package, default  
###               is "Set1". Automatically adjusts number of colors to match
###               the unique number of factors in the 'by' parameter. 
###  Returns: df with "by" names and a unique hexidecimal colors for each 
###  Notes: Used in  other functions. Color palette can either be a typical 
###    palette (e.g. rainbow, heat.heat.colors, terrain,colors, ) using the
###    'pal' parameter or a RColorBrewer palette using the 'r_pal' parameter, or
###    a RColorBrewer palette using the "b_pal" parameter (e.g. "Accent"). The 
###    number of colors is always automatically adjusted to match the unique 
###    number of factors in the "by" parameter column of the df.
###  Blake Massey
###  2014.11.26


CreateColorsByVar <- function (df,
                               by = NULL,
                               pal = NULL,
                               r_pal = NULL,
                               b_pal = "Set1",
                               output = TRUE,
                               plot = FALSE) {
  suppressPackageStartupMessages(require(graphics))
  suppressPackageStartupMessages(require(plotKML))
  suppressPackageStartupMessages(require(RColorBrewer))
  source('C:/Work/R/Functions/gen.R')
  if(is.null(by)){ 
    var_colors <- c("all" = "#377EB8")  # blue from RColorBrewer's "Set1"
  } else {
  vars <- unique(df[,by])
  vars_n <- as.numeric(length(unique(df[,by])))
  PalFunction <- function(x, fun) {
    fun(x)
  }
  if (!is.null(pal)) var_colors <- PalFunction(vars_n, pal)
  if (is.null(pal) && (!is.null(r_pal))) { 
    var_colors <- colorRampPalette(R_pal[[r_pal]])(vars_n)
  }
  if (is.null(pal) && (is.null(r_pal))) {
    ifelse(vars_n <= 8, var_colors <- brewer.pal(vars_n,b_pal), var_colors <- 
      colorRampPalette(brewer.pal(vars_n,b_pal))(vars_n))
  }
  for (i in 1:length(vars)) {
    names(var_colors)[i] <- vars[i] 
  }
  }
  return(var_colors)
}

# ExportKMLTelemetry Function --------------------------------------------------

###  Create a Google Earth KML file (points and multitrack) from lat/long 
###    coordinates
###  Usage: ExportKMLTelemetry(df, id, datetime, lat, long, speed, alt,alt_mode, 
###    alt_mode, agl, behavior, point_color, point_metadata, point_pal, 
###    point_r_pal, point_b_pal, path_color, path_metadata, path_pal 
###    path_r_pal, path_b_pal, extrude, labelscale, dateformat, 
###    timeformat, datetimeformat) 
###  Arguments: df = input dataframe, must have id, lat, long, and datetime 
###             id = column name of unique identifier, data is split into unique
###               paths and separate folders based on this parameter
###             datetime = column name of datetime in POSIXct format or as a 
###               character in the format (%Y/%m/%d %H:%M)
###             lat = column name of latitude coordinates (WGS84, dec. degree)
###             long = column name of longitude coordinates (WGS84, dec. degree)
###             speed = input dataframe column name for speed. Optional
###             alt = input dataframe column name for altitude(m). Optional.
###             alt_mode = based on KML code: "absolute","clampedToGround",
###               "relativeToGround" (see KML documentation for description).
###               Default is "clampedToGround".
###             agl =  input dataframe column name for "altitude above ground 
###               level", optional
###             behavior = input dataframe column name for behavior. Optional
###             point_color = column name that determines the color for each 
###               point, may be same as 'id' parameter, but may also be sex, 
###               behavior, season, etc. Default is 'id' parameter
###             point_metadata = location of metadata .csv file. Metadata file 
###               must have a column that matches name of 'point_color'
###               parameter and "icon_color" column with hexadecimal colors.
###             point_pal = name of color palette funtions (e.g., rainbow, 
###               heat.colors, terrain.colors, topo.colors, cm.colors
###               used to create colors. This parameter has priority over the 
###               other point color palette parameters. Default is NULL.
###             point_r_pal = Specifc number of 'R_pal' color palette from the 
###               'PlotKML' Package (e.g., 1 = R_pal[[1]]). This parameter has 
###               priority over the 'b_pal' parameter for setting the colors. 
###               Default is NULL.
###             point_b_pal = color palette name from RColorBrewer package, 
###               default is "Set1". Automatically adjusts number of colors to 
###               match the unique number of factors in the 'point_color' 
###               column of the input dataframe.
###             extrude = logical, either FALSE (default) for no line, or TRUE 
###               which extends a line from the point to the ground.
###             path = logical, to create Track paths. Default is TRUE.
###             path_color = similar to 'point_color' parameter, but the value
###               must have the same factor level structure as the id file, 
###               because each path is constructed for each id factor.
###               Default will use 'id' parameter.
###             path_metadata = location of metadata .csv file. Metadata file 
###               must have a column that matches name of 'path_color'
###               parameter and an "icon_color" column with hexadecimal colors.
###             path_pal = name of color palette funtions (e.g., rainbow, 
###               heat.colors, terrain.colors, topo.colors, cm.colors
###               used to create colors. This parameter has priority over the 
###               other point color palette parameters. Default is NULL.
###             path_r_pal = Specifc number of 'R_pal' color palette from the 
###               'PlotKML' Package (e.g., 1 = R_pal[[1]]). This parameter has 
###               priority over the 'b_pal' parameter for setting the colors. 
###               Default is NULL.
###             path_b_pal = color palette name from RColorBrewer package, 
###               default is "Set1". Automatically adjusts number of colors to 
###               match the unique number of factors in the 'point_color' 
###               column of the input dataframe.
###             kml_folder = name for folder in the KML file, default is name of 
###               'df' parameter
###             outfile = filepath of output KML file, default is working 
###               directory and name of 'df' parameter
###             labelscale = adjusts the size of the Google Earth location 
###               point labels. Default is 0, which hides the labels. To show 
###               labels, change to a value between 0.7-1.
###             dateformat = changes the format of the date in the Google Earth 
###               location pop-up windows. Default is "%Y/%m/%d".
###             timeformat = changes the format of the time in the Google Earth  
###               locations pop-up windows. Default is "%I:%M %p".        
###             datetimeformat = changes the datetime format for the label of 
###               highlighted points. Default is "%Y/%m/%d %I:%M %p"
###  Returns: KML of points and multitracks          
###  Notes: 
###  Blake Massey
###  2014.06.02

ExportKMLTelemetry <- function (df,
                                id = "id",                                  
                                datetime = "datetime",
                                lat = "lat",
                                long = "long", 
                                alt = NULL,
                                alt_mode = "clampToGround",
                                speed = NULL,
                                agl = NULL,
                                behavior = NULL,
                                point_color = NULL,
                                point_metadata = NULL,
                                point_pal = NULL,
                                point_r_pal = NULL,
                                point_b_pal = "Set1", 
                                extrude = FALSE,
                                path = TRUE,
                                path_color = NULL,
                                path_metadata = NULL,
                                path_pal = NULL,
                                path_r_pal = NULL,
                                path_b_pal = NULL,
                                arrow = TRUE,
                                icon_by_sex = FALSE,
                                labelscale = 0, 
                                dateformat = "%Y-%m-%d", 
                                timeformat = "%I:%M %p",
                                datetimeformat = "%Y-%m-%d %I:%M %p",
                                outfile = NULL,
                                kml_folder = NULL) {
  suppressPackageStartupMessages(require(plotKML))
  suppressPackageStartupMessages(require(tools))
  if (is.null(kml_folder) == TRUE) {
    if (!is.null(outfile)) {
      kml_folder <- basename(outfile)
      kml_folder <- sub(".kml", "", kml_folder, ignore.case =TRUE)
      kml_folder <- sub(".kmz", "", kml_folder, ignore.case =TRUE)
    } else {
      kml_folder <- deparse(substitute(df))
    }
  }
  if (is.null(outfile)) {
    if (!is.null(kml_folder)) {
      outfile <- paste(getwd(), "/", kml_folder, sep="")
      } else {
        outfile <- paste(getwd(), "/", deparse(substitute(df)), sep="")
      }
  } else {
    if (dirname(outfile) == "."){
      outfile <- paste(getwd(), "/", outfile, sep="") 
    }
  }
  if (file_ext(outfile) == "") {
    outfile <- paste(outfile, ".kml", sep="")  # if object has no extension
  } 
  df <- df 
  df$id <- df[ ,id]  
  df$lat <- df[ ,lat]  
  df$long <- df[ ,long] 
  if (!is.null(alt)){
    df$desc_alt <- df[ ,alt]  # writes altitude to the "alt" column
    alt1 <- '\t\t\t\t\tAltitude: '  # first part of "Altitude" description
    alt2 <- '\n'  # second part of "Altitude" description
  } else {
    df$desc_alt <- ""  # makes the altitude column a vector of blank values  
    alt1 <- NULL  # prevents "Altitude" description from being written
    alt2 <- NULL  # prevents "Altitude" description from being written    
  }
  if (!is.null(agl)) {
    df$desc_agl <- df[ ,agl]  
    agl1 <- '\t\t\t\t\tAGL: '  # first part of the "Altitdue" description
    agl2 <- '\n'  # second part of "Altitude" description
  } else {
    df$desc_agl <- ""  # makes the altitude column a vector of blank values  
    agl1 <- NULL  # prevents "Altitude Above Ground Level" from being written
    agl2 <- NULL  # prevents "Altitude Above Ground Level" from being written
  }
  if (!is.null(speed)) {
    df$desc_speed <- df[,speed]  # writes speed to the "speed" column
    spd1 <- '\t\t\t\t\tSpeed: '  # writes first part of the "Speed" description
    spd2 <- '\n'  # writes second part of the "Speed" description
  } else {
    df$desc_speed <- ""  # makes the speed column a vector of blank values  
    spd1 <- NULL  # prevents "Speed" description from being written
    spd2 <- NULL  # prevents "Speed" description from being written    
  }
  if (!is.null(behavior)) {
    df$desc_behavior <- df[,behavior]  # writes behavior to "behavior" column
    beh1 <- '\t\t\t\t\tBehavior: '  # first part of the "Behavior" description
    beh2 <- '\n'  # second part of the "Behavior" description
  } else {
    df$desc_behavior <- ""  # makes the behavior column a vector of blank values  
    beh1 <- NULL  # prevents "Behavior" description from being written
    beh2 <- NULL  # prevents "Behavior" description from being written    
  }
  df$datetime <- df[,datetime]
  df$datetime <- as.character(df$datetime)  # needed for KML parsing
  df$datetimebegin <- df$datetime
  EndTimes <- function(data) { # locates last time in data 
    data$datetime2 <- data$datetime[c(2:length(data$datetime), 
    length(data$datetime))]
  } 
  ids <- as.character(df$id)  # as.character removes factor levels
  df_split <- split(df, ids)  # divides data by ids 
  df_split <- lapply(df_split, EndTimes)
  datetimeend <- unsplit(df_split, ids)  # returns array of returned values
  df <- cbind(df, datetimeend)  # adds datetimeend column to original baea data
  ifelse(extrude == TRUE, extrude <- 1, extrude <- FALSE)
  PlacemarkPoint <- function(PN, X,  Y, Z, ZD, 
                             AG, SP, BH, SX, PS, 
                             ID, SD, ST, ED, ET, 
                             DA, TM) {
    if (icon_by_sex == TRUE) PS <- paste0(PS,"-",SX)
    cat("\t<Placemark>\n",
      "\t\t<name>",PN, "</name>\n",
      "\t\t<TimeSpan>\n",
      "\t\t\t<begin>",SD ,"T" ,ST ,"</begin> " , "\n",
      "\t\t\t<end>", ED, "T", ET, "</end> ", "\n",
      "\t\t</TimeSpan>\n",
      "\t\t\t<Snippet></Snippet>", "\n",
      "\t\t\t\t<description>\n",
      "\t\t\t\t\tID: ", ID, "\n",
      "\t\t\t\t\tDate: ", DA, "\n",
      "\t\t\t\t\tTime: ", TM, "\n",
      "\t\t\t\t\tLongitude: ", X, "\n",
      "\t\t\t\t\tLatitude: ", Y, "\n",
      alt1, ZD, alt2,  # written when !is.null(agl)
      agl1, AG, agl2, # written when !is.null(agl)
      spd1, SP, spd2,  # written when !is.null(speed)
      beh1, BH, beh2, # written when !is.null(behavior)
      "\t\t\t\t</description>", "\n",  
      "\t\t\t<styleUrl>#Point_",PS,"</styleUrl>\n", 
      "\t\t\t<Point>\n",
      "\t\t\t\t<extrude>", extrude, "</extrude>\n",
      "\t\t\t\t<altitudeMode>", alt_mode, "</altitudeMode>\n",
      "\t\t\t\t\t<coordinates>", X, ",", Y, ",", Z, "</coordinates>\n",
      "\t\t\t</Point>\n",
      "\t</Placemark>\n",
      file = outfile, append = TRUE, sep = "")  
    }  
  if (file.exists(outfile)) file.remove(outfile)  # delete KML if already exists
  writeLines(noquote(c("Writing: ", outfile)))
  ## Title Section ##
  cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n",      
  "<kml xmlns=\"http://www.opengis.net/kml/2.2\"\n",      
  "xmlns:gx=\"http://www.google.com/kml/ext/2.2\"\n",
  "xmlns:atom=\"http://www.w3.org/2005/Atom\">\n\n",
  "<Document>\n", "\t<name>",kml_folder,"</name>\n",file = outfile, 
  append = FALSE, sep = "") 
  ## Icon Style Section ##  
  if (is.null(point_color)) point_color <- id
  df$point_color <- df[ ,point_color]
  if (!is.null(point_metadata)) {
    point_colors <- CreateColorsByMetadata(file=point_metadata, 
      metadata_id=point_color)
    point_colors <- subset(point_colors, names(point_colors) %in% 
      unique(df$point_color))
  } else {
    suppressWarnings(point_colors <- CreateColorsByVar(by=point_color, df=df, 
      pal=point_pal, r_pal=point_r_pal, b_pal=point_b_pal))
  }  
  if (icon_by_sex == TRUE) {  
    point_colors_names <- c(sapply(names(point_colors), paste0, "-female"), 
      sapply(names(point_colors), paste0,"-male"))
    point_colors <- rep(point_colors, 2)
    names(point_colors) <- point_colors_names  
  }
  point_colors <- sapply(point_colors, col2kml)
  point_colors <- sapply(point_colors, substring, 4, 9)
  df$point_color <-df[,point_color]
  icon_scale <- 0.7
  hi_icon_label_scale <- 0.75
  ball_bg_color <- "ff333333"
  ball_text_color <- "ffffffff"    
  mt_icon_href <- file.path("http://earth.google.com/images/kml-icons",
    "track-directional/track-0.png")
  icon_href <- "http://maps.google.com/mapfiles/kml/shapes/placemark_square.png"
  
  for (i in 1:length(point_colors)) {  
    if (icon_by_sex == TRUE){  
      if (grepl("male", names(point_colors)[i]) == TRUE)  icon_href <- 
        "http://maps.google.com/mapfiles/kml/shapes/placemark_square.png"  
      if (grepl("female", names(point_colors)[i]) == TRUE) icon_href <- 
        "http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png"  
    }
    cat("\t<StyleMap id=\"Point_",names(point_colors)[i],"\">\n",
      "\t\t<Pair>\n",
      "\t\t\t<key>normal</key>\n",
      "\t\t\t\t<Style>\n",  
      "\t\t\t\t\t<LabelStyle>\n",      
      "\t\t\t\t\t<scale>",labelscale,"</scale>\n",  # to show label 
      "\t\t\t\t\t</LabelStyle>\n",
      "\t\t\t\t\t<IconStyle>\n",
      "\t\t\t\t\t\t<color>FF",point_colors[i],"</color>\n",
      "\t\t\t\t\t\t<scale>",icon_scale,"</scale>\n", 
      "\t\t\t\t\t<Icon>\n",
      "\t\t\t\t\t<href>",icon_href,"</href>\n",
      "\t\t\t\t\t</Icon>\n",
      "\t\t\t\t\t</IconStyle>\n",      
      "\t\t\t\t\t<BalloonStyle>\n",
      "\t\t\t\t\t<text>$[description]</text>\n", 
      "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>\n",
      "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>\n",
      "\t\t\t\t\t</BalloonStyle>\n",
      "\t\t\t\t</Style>\n",
      "\t\t</Pair>\n",
      "\t\t<Pair>\n",
      "\t\t\t<key>highlight</key>\n",
      "\t\t\t\t<Style>\n",
      "\t\t\t\t\t<LabelStyle>\n",     
      "\t\t\t\t\t<scale>",hi_icon_label_scale,"</scale>\n",  # to show label 
      "\t\t\t\t\t</LabelStyle>\n",
      "\t\t\t\t\t<IconStyle>\n",
      "\t\t\t\t\t\t<color>FF",point_colors[i],"</color>\n", 
      "\t\t\t\t\t\t<scale>0.8</scale>\n", 
      "\t\t\t\t\t<Icon>\n",
      "\t\t\t\t\t<href>",icon_href,"</href>\n",
      "\t\t\t\t\t</Icon>\n",
      "\t\t\t\t\t</IconStyle>\n",      
      "\t\t\t\t\t<BalloonStyle>\n",
      "\t\t\t\t\t<text>$[description]</text>\n", 
      "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>", "\n",
      "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>", "\n",
      "\t\t\t\t\t</BalloonStyle>\n",
      "\t\t\t\t</Style>\n", 
      "\t\t</Pair>\n",
      "\t</StyleMap>\n",
      file = outfile, append = TRUE, sep = "")
  } 
  if (path ==TRUE) {
    if (is.null(path_color)) path_color <- id
    if (is.null(path_b_pal)) path_b_pal <- point_b_pal
    if (!is.null(path_metadata)) {
      path_colors <- CreateColorsByMetadata(file=path_metadata, id=path_color)
      path_colors <- subset(path_colors, names(path_colors) %in% 
        unique(df[,path_color]))
    } else {
      suppressWarnings(path_colors <- CreateColorsByVar(by=path_color, df=df, 
        pal=path_pal, r_pal=path_r_pal, b_pal=path_b_pal))
    }  
    path_colors <- sapply(path_colors, col2kml)
    path_colors <- sapply(path_colors, substring, 4, 9)
    if (icon_by_sex == TRUE) {  
      path_colors_names <- c(sapply(names(path_colors), paste0, "-female"), 
      sapply(names(path_colors), paste0,"-male"))
      path_colors <- rep(path_colors, 2)
      names(path_colors) <- path_colors_names  
    }
    ifelse(arrow == TRUE, arrow <- 1, arrow <- 0)
  ## Style Map for Track ##
    for (i in 1:length(path_colors)) {    
      cat("\t<StyleMap id=\"Track_",names(path_colors)[i],"\">\n",
      "\t\t<Pair>\n",
      "\t\t\t<key>normal</key>\n",
      "\t\t\t\t<Style>\n", 
      "\t\t\t\t\t<LabelStyle>\n",      
      "\t\t\t\t\t<scale>0</scale>\n",  # to show label 
      "\t\t\t\t\t</LabelStyle>\n",
      "\t\t\t<IconStyle>\n",
      "\t\t\t\t<scale>",arrow,"</scale>\n",
      "\t\t\t\t<Icon>\n",
      "\t\t\t\t\t<href>",mt_icon_href,"</href>\n",
      "\t\t\t\t</Icon>\n",
      "\t\t\t</IconStyle>\n",
      "\t\t\t\t\t<BalloonStyle>\n",
      "\t\t\t\t\t<text>",names(path_colors)[i]," - Path</text>\n",       
      "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>\n",
      "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>\n",
      "\t\t\t\t\t</BalloonStyle>\n",
      "\t\t\t<LineStyle>\n",
      "\t\t\t\t<color>dd",path_colors[i],"</color>\n",
      "\t\t\t\t<width>1</width>\n",
      "\t\t\t</LineStyle>\n",
      "\t\t\t\t</Style>\n",
      "\t\t</Pair>\n",    
      "\t\t<Pair>\n",
      "\t\t\t<key>highlight</key>\n",
      "\t\t\t\t<Style>\n", 
      "\t\t\t\t\t<LabelStyle>\n",
      "\t\t\t\t\t<scale>0</scale>\n",  # to show label, change value to >= 0.7
      "\t\t\t\t\t</LabelStyle>\n",
      "\t\t\t<IconStyle>\n",
      "\t\t\t\t<scale>1</scale>\n",
      "\t\t\t\t<Icon>\n",
      "\t\t\t\t\t<href>",mt_icon_href,"</href>\n",
      "\t\t\t\t</Icon>\n",
      "\t\t\t</IconStyle>\n",
      "\t\t\t\t\t<BalloonStyle>\n",
      "\t\t\t\t\t<text>",names(path_colors)[i]," - Path</text>\n",
      "\t\t\t\t\t\t<bgColor>",ball_bg_color,"</bgColor>\n",
      "\t\t\t\t\t\t<textColor>",ball_text_color,"</textColor>\n",
      "\t\t\t\t\t</BalloonStyle>\n",
      "\t\t\t<LineStyle>\n",
      "\t\t\t\t<color>ee",path_colors[i],"</color>\n",
      "\t\t\t\t<width>1</width>\n",
      "\t\t\t</LineStyle>\n",
      "\t\t\t\t</Style>\n",
      "\t\t</Pair>\n",
      "\t</StyleMap>\n",
      file = outfile, append = TRUE, sep = "")
    }  # end of Track icon loop
  } # end of (path == TRUE) 

  ids <- as.character(unique(df$id))  # as.character removes factor levels
  for (i in ids) {
    sv = df$id %in% i
    unique_id <- as.character(unique(df$id[sv]))    
    cat("<Folder>\n","<name>",unique_id,"</name>\n","<open>0</open>\n",
      file = outfile, append = TRUE, sep = "")
    cat("\t<Folder>\n","\t<name>",unique_id," - Locations</name>\n",
      "\t<open>0</open>\n", file = outfile, append = TRUE, sep = "")
    locs <- subset(df, id == unique_id)
    for (i in 1:nrow(locs)){
      loc <- locs[i,]
      PNs <- strftime(loc[, "datetimebegin"], datetimeformat)
      Xs <- loc[, "long"]
      Ys <- loc[, "lat"]
      Zs <- loc[, "alt"]
      ZDs <- loc[, "desc_alt"]      
      SPs <- loc[, "desc_speed"]
      AGs <- loc[, "desc_agl"]
      BHs <- loc[, "desc_behavior"]
      SXs <- loc[, "sex"]
      PSs <- loc[, "point_color"] 
      IDs <- unique_id 
      SDs <- substring(loc$datetime, 1,10) #start date
      STs <- substring(loc$datetime, 12,16) #start time 
      EDs <- substring(loc$datetimeend, 1,10) #end date
      ETs <- substring(loc$datetimeend, 12,19) #end time 
      DAs <- strftime(loc[, "datetimebegin"], dateformat)
      TMs <- strftime(loc[, "datetimebegin"], timeformat)      
      PlacemarkPoint(PNs, Xs, Ys, Zs, ZDs,
                     AGs,  SPs, BHs, SXs, PSs, 
                     IDs,  SDs, STs, EDs,  ETs, 
                     DAs, TMs)     
    }      
    cat("\t</Folder>\n", file = outfile, append = TRUE, sep = "")
    locs$Ts <- "T"
    locs$Zs <- "Z"
    locs$datetimedate <-substring(locs$datetime, 1,10) #start date
    locs$datetimetime <- substring(locs$datetime, 12,16) #start time
    whens <- locs[, c("datetimedate","Ts","datetimetime", "Zs")]
    sgmts <- locs[, c("long","lat","alt")]
    unique_id <- unique(locs$id)
    ifelse(icon_by_sex == TRUE, path_id <- paste0(unique_id, "-", 
      unique(locs$sex)), path_id <- unique_id)
    bloc2 <- NULL
    if (path == TRUE) {
      bloc2 <- c(bloc2, paste(
        "\t<Placemark>\n",
        "\t\t<name>",unique_id," - Path</name>\n",
        "\t\t<styleUrl>#Track_",path_id,"</styleUrl>\n",
        "\t\t<gx:balloonVisibility>0</gx:balloonVisibility>\n",
        "\t\t<gx:Track>\n",
        "\t\t<altitudeMode>",alt_mode,"</altitudeMode>\n",
      paste(paste("\t\t\t\t\t<when>", apply(whens, 1, paste, collapse=""), 
        sep=""), "</when>", collapse="\n"), "\n", 
      paste(paste("\t\t\t\t\t<gx:coord>", apply(sgmts, 1, paste, collapse=" "), 
        sep = ""),"</gx:coord>", collapse = "\n"),"\n", 
        "\t\t</gx:Track>\n",
        "\t</Placemark>\n",                
        sep = ""))  
    }
    cat(bloc2, "\t</Folder>\n", file = outfile, append = TRUE)                                      
  }
  cat("</Document>\n</kml>", file = outfile, append = TRUE)
}