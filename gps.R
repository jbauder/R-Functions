# --- GPS FUNCTIONS ------------------------------------------------------------
# General functions for GPS telemetry data. Not specific to BAEA data.
# ------------------------------------------------------------------------------

# AddFirstLastDistance Function ------------------------------------------------

###  Adds columns with distance to first and last location of each day
###  Usage: AddFirstLastDistance(df, id, datetime, tz)
###  Arguments: df = dataframe with locations
###             by = column name of unique identifier, default is "id"
###             datetime = column name of datetime in POSIXct format
###               or as a character in the format (%Y/%m/%d %H:%M)
###             tz = return to this timezome. Default "Etc/GMT+5".
###  Returns: dataframe with "datetime", "first", "last"         
###  Notes: Caution with use of tz, may cause problems
###  Blake Massey
###  2014.11.11

AddFirstLastDistance <- function(df = df,
                                 by = "id",
                                 datetime = "datetime",
                                 tz = "Etc/GMT+5"){
  suppressPackageStartupMessages(require(plyr))  # for adply function
  df$by <- df[,by]
  df$datetime <- df[,datetime]
  df$date <- as.Date(df[,datetime], tz=tz)
  df$datetime <- as.character(df$datetime)  # convert to character
  df$date <- as.character(df$date)  # convert to character
  first <- ddply(df, .(date, by), function(x) x[1, ])  # first records
  last <- ddply(df, .(date, by), function(x) x[(nrow(x)), ])  # last records
  first$first <- "First"  # populates "first" column
  last$last <- "Last"  # populates "last" column
  first <- subset(first, select = c(by, first, date, datetime, long_utm, 
    lat_utm))  # subsets first records
  colnames(first) <- c("by", "first", "date", "datetime","long_utm_first",
    "lat_utm_first")  # need to id first locations
  last <- subset(last, select = c(by, last, date, datetime, long_utm,
    lat_utm))  # subsets last records
  colnames(last) <- c("by", "last", "date", "datetime","long_utm_last",
    "lat_utm_last")  # need to id last locations
  first_date <- subset(first, select=c(by, date, long_utm_first,
    lat_utm_first))  # all first locations
  last_date <- subset(last, select=c(by, date, long_utm_last, lat_utm_last))
    # all last locations
  first_last_date <- merge(first_date, last_date, by = c("by","date"),
    all.x = TRUE)  # first and last record locations
  first_last <- merge(first, last, all=TRUE)  # merges first and last records 
  first_last_datetime <- subset(first_last, select = c(by, first, last,
    datetime))  # only first and last records and their times
  df <- merge(df, first_last_datetime, by=c("by", "datetime"), all.x=TRUE)
  df <- merge(df, first_last_date, by=c("by", "date"), all.x=TRUE)
      #adds the lat and long for the first and last records to the dataframe
  df<-adply(df, 1, transform, dist_first = as.integer(sqrt(sum((c(long_utm, 
    lat_utm) - c(long_utm_first, lat_utm_first))^2))))  # Pythagorean theorem
  df<-adply(df, 1, transform, dist_last = as.integer(sqrt(sum((c(long_utm,
    lat_utm) - c(long_utm_last, lat_utm_last))^2)))) # Pythagorean theorem
  df$date <- as.Date(df$date, "%Y-%m-%d") #convert date into date format 
  df$datetime <- strptime(df$datetime, "%Y-%m-%d %H:%M:%S")  # back to POSIXct
  df$datetime <- as.POSIXct(df$datetime, tz=tz, usetz=FALSE)  # returns to tz
  drops <- c("long_utm_first", "lat_utm_first", "long_utm_last", 
             "lat_utm_last", "by")  # vector of columns to drop
  df<-df[ ,!(names(df) %in% drops)] 
  return(df)
}

# AddSegmentLengthTime (deprecated) Function -----------------------------------

###  Caculates distances and times between successive points for each individual 
###  Usage: AddSegmentLengthTime(df, id, datetime, long, lat)
###  Arguments: df = dataframe
###             id = "id"
###             datetime = "datetime" 
###             long = "long"
###             lat = "lat"
###  Returns: dataframe with "seg_length" and "time_diff" columns         
###  Notes: projects to UTM Zone 19 for Maine (from rgdal Package)
###  Blake Massey
###  2014.05.05

AddSegmentLengthTime <- function(df = df, 
                                 id = "id", 
                                 datetime = "datetime", 
                                 long = "long", 
                                 lat = "lat") {
  suppressPackageStartupMessages(require(circular))  # not sure this is needed
  suppressPackageStartupMessages(require(move))  # needed to create 'move' obj
  if(!("lat_utm" %in% colnames(df))) { 
    xy <- (cbind(df[ ,long], df[ ,lat]))  # for projection function in next step
    xy <- project(xy, "+proj=utm +zone=19 ellps=WGS84")  
    colnames(xy) <- c("long_utm", "lat_utm")
    xy <- round(xy)  # rounds lat and long to zero decimal places
    df <- cbind(df, xy)  # combines lat long with originial data  
  }
  SegmentLength <- function(df){ #caculates distances between successive points
    move_object <- move(x = df[ ,"long_utm"], y = df[, "lat_utm"], 
      time = df[,datetime], proj = CRS("+proj=utm +zone=19 ellps=WGS84"),
      animal = df[,id], sensor = "GPS") #convert baea to a move object
    seg_length <- as.integer(c(seglength(move_object), NA))  # NA at end  
  }
  SegmentTime <- function(df){  # calculates time between successive points
    t <- df[ ,datetime]
    diff <- difftime(head(t, -1), tail(t, -1))
    time_diff<-as.integer(c(-1*as.numeric(diff, units = "mins"), NA))
  }
  sn<-unique(df[,id])
  for (s in sn){
    sv = df[,id] %in% s
    individual <- df[df[ ,id] == s,]
    df$seg_length[sv] <- SegmentLength(individual)
    df$time_diff[sv] <- SegmentTime(individual)
  }     
  return(df) #returns the dataframe
}

# AddStepLengthAndAnglesFunction -----------------------------------------------

###  Calculates step length, absolute angle (N=0), and turn angle between 
###    successive point locations 
###  Usage: AddStepLengthAndAngles(df, by, datetime, long, lat)
###  Arguments: df = dataframe with locations and datetime data
###             by = column name to use to split, analyze, and merge data, 
###               default is "id" 
###             long = longitude, must be in UTM, default is "long_utm"
###             lat = latitude, must be in UTM, default is "lat_utm"
###             radian = logical, whether to return the angles in radian, 
###               default is FALSE
###  Returns: dataframe with "step_length", "abs_angle", and "turn_angle"
###    columns         
###  Notes: Coordinates must be in identically-scaled units (e.g. UTM meters).  
###   If data is in lat, long If needed, project coordinates to UTM with 'rgdal' 
###   package before running this function
###  Blake Massey
###  2014.11.11


AddStepLengthAndAngles <- function(df,
                                   by = "id", 
                                   long = "long_utm",
                                   lat = "lat_utm"){
  suppressPackageStartupMessages(require(plyr))
  df <- df
  ifelse(is.null(by), df$by <- "all", df$by <- df[,by])
  StepLengthAndAngles <- function(df=df, lat=lat, long=long){
    xy <- data.frame(x = df[,long], y = df[,lat])
    xy1 <- xy[-1, ]
    xy2 <- xy[-nrow(xy), ]
    step_length <- c(sqrt((xy1$x - xy2$x)^2 + (xy1$y - xy2$y)^2), NA)
    dx <- c(xy1$x - xy2$x, NA)
    dy <- c(xy1$y - xy2$y, NA)
    abs_angle <- ifelse(step_length < 1e-07, NA, (atan2(dy, dx)/(pi/180)))
      # if N=0, then in the previous line: atan2(dx,dy)
    abs_angle <- ifelse(abs_angle < 0, 360 + abs_angle, abs_angle)
    turn_angle <- c(abs_angle, 0) - c(0, abs_angle)
    turn_angle <- ifelse(turn_angle < 0, 360 + turn_angle, turn_angle)
    turn_angle <- turn_angle[-length(turn_angle)]
    turn_angle[1] <- NA 
    turn_angle_rad <- turn_angle*(pi/180)
#    turn_angle_180 <- ifelse(turn_angle>180, 360-turn_angle, turn_angle) 
#    turn_angle_180_rad <- turn_angle_180*(pi/180)
    out <- cbind.data.frame(dx=dx, dy=dy, step_length=step_length, 
      abs_angle=abs_angle, turn_angle=turn_angle, 
      turn_angle_rad=turn_angle_rad) #,
#      turn_angle_180=turn_angle_180, turn_angle_180_rad=
#      turn_angle_180_rad)
  }
  uniques <- unique(df[,"by"])
  out <- data.frame()
  for (j in uniques){
    sv = df[,by] %in% j
    data <- subset (df, by==j)
    df2<- cbind(data, StepLengthAndAngles(df=data, lat=lat, long=long)) 
    out <- rbind(out, df2)
  } 
  out$by<-NULL
  return(out)
}

# AddStepTime Function ---------------------------------------------------------

###  Calculates time between successive point locations
###  Usage: AddStepTime(df, by, datetime)
###  Arguments: df = dataframe with locations and datetime data
###             by = column name to use to split, analyze, and merge data, 
###               default is "id" 
###             datetime = datetime column name, default is "datetime"
###  Returns: dataframe with "step_time" column         
###  Notes: 
###  Blake Massey
###  2014.11.11

AddStepTime<- function(df, 
                       by = "id", 
                       datetime = "datetime"){
  df <- df
  df$by <- df[,by]
  StepTime <- function(df, datetime=datetime){
    t <- df[,datetime]
    diff <- difftime(head(t, -1), tail(t, -1))
    step_time <- as.integer(c(-1*as.numeric(diff, units = "mins"), NA))
  }
  
  uniques <- unique(df[,by])
  out <- data.frame()
  for (j in uniques){
    sv = df[,by] %in% j
    data <- subset (df, by==j)
    df2 <- cbind(data, step_time = StepTime(df=data, datetime=datetime)) 
    out <- rbind(out, df2)
  }
  out$by<-NULL
  return(out)
}

# AddSolarTimes Function -------------------------------------------------------

###  Add sunrise, sunset, and solarnoon to a dataframe of location data
###  Usage: AddSolarTimes(df, id, tz)
###  Arguments: df = dataframe
###             by = column name to use to split, analyze, and merge data, 
###               default is "id". 
###             tz = timezone (default = "Etc/GMT+5") 
###  Returns: df with sunrise, sunset, and solarnoon times          
###  Notes: "lat" and "long" must be WGS84, sunrise and solarnoon based on first
###    location, sunset based on last location.
###  Blake Massey
###  2014.05.05


AddSolarTimes <- function(df = df, 
                          by = "id", 
                          tz = "Etc/GMT+5"){  
  suppressPackageStartupMessages(require(lubridate)) # needed for tz function
  suppressPackageStartupMessages(require(maptools))
  suppressPackageStartupMessages(require(plyr))
  df <- df
  df$by <- df[,by]
  if( ! ("date" %in% colnames(df))) {
    df$date <- as.Date(df$datetime,tz = tz)
  }
  AddTimes<-function (df=df){
    first <- ddply(df, .(date), function(x) x[1, ])  # first records
    sunrise_coords <- cbind(first$long, first$lat)
    sunrise_datetime <- first$datetime
    sunrise <- sunriset(sunrise_coords, sunrise_datetime, proj4string = 
      CRS("+proj=longlat +datum=WGS84"), direction = "sunrise", 
      POSIXct.out= TRUE)
    solarnoon <- solarnoon(sunrise_coords, sunrise_datetime, proj4string = 
      CRS("+proj=longlat +datum=WGS84"), POSIXct.out = TRUE)
    sunrise$date <- as.Date(sunrise$time, tz = tz)
    sunrise$sunrise <- sunrise$time
    sunrise <- subset(sunrise, select = c(date, sunrise))
    solarnoon$date <- as.Date(solarnoon$time, tz=tz)
    solarnoon$solarnoon <- solarnoon$time
    solarnoon <- subset(solarnoon, select = c(date, solarnoon))
    last <- ddply(df, .(date), function(x) x[(nrow(x)), ])  # last records
    sunset_coords <- cbind(last$long, last$lat)
    sunset_datetime <- last$datetime
    sunset <- sunriset(sunset_coords, sunset_datetime, proj4string = 
      CRS("+proj=longlat +datum=WGS84"), direction = "sunset", 
      POSIXct.out= TRUE)
    sunset$date <- as.Date(sunset$time, tz=tz)
    sunset$sunset <- sunset$time
    sunset <- subset(sunset, select = c(date, sunset))
    df <- merge(df, sunrise, by="date", all.x = TRUE)
    df <- merge(df, solarnoon, by="date", all.x = TRUE)
    df <- merge(df, sunset, by="date", all.x = TRUE)
    df$hr_before_sunrise <- df$sunrise - 3600  # subtracts an hour from sunset
    tz(df$sunrise) <- tz  # sets timezone for sunrise times
    tz(df$hr_before_sunrise) <- tz  # sets timezone for sunrise times
    df$hr_after_sunset <- df$sunset + 3600  # adds an hour to sunset
    tz(df$sunset) <- tz  # sets timezone for sunrise times
    tz(df$hr_after_sunset) <- tz  # sets timezone for sunrise times
  return(df)
  }
  df2 <- ddply(df, .(by), AddTimes)
  df2$by <- NULL
  return(df2)
}

# CompileDownloads Function ----------------------------------------------------

###  Compiles the downloaded CTT files
###  Usage: CompileDownloads(units, compile, tz)
###  Arguments: units = "deployed" or "reserve", default = "deployed"
###             compile = "all" or "recent", default = "recent"
###             tz = timezone, default is "Etc/GMT+5"
###  Returns: df of compiled files          
###  Notes: internal parameters are set specifically for Maine data
###  Blake Massey
###  2014.05.12

CompileDownloads <- function(units = "deployed", 
                             compile = "all", 
                             tz = "Etc/GMT+5") {
  suppressPackageStartupMessages(require(rgdal))
  suppressPackageStartupMessages(require(sp))  
  infile <- paste("C:/Work/R/Data/Telemetry/", units, "/",
    compile, sep = "")
  filenames <- list.files(path=infile, full.names=TRUE)
  output <- paste("into ", units, "_", compile, sep="")
  writeLines(noquote(paste(c("Compiling files:", filenames, output))))
  suppressWarnings(df <- do.call("rbind", lapply(filenames, read.csv,colClasses= 
    c("character", "character", "character", "character", "character", 
    "character", "numeric", "integer", "integer", "numeric", 
    "numeric", "numeric", "character", "character","character", 
    "integer", "integer", "character", "character"), header = TRUE, 
    na.strings = "")))
  df <-  subset(df, select=serial:alt)
  date <- Sys.Date()
  if (compile == "recent" || compile == "all") {
    outfile <- paste("C:/Work/R/Data/Telemetry/", units, "/Archive/", compile, 
      "/",date,".csv", sep ="")
    if ( !file.exists(outfile)) {  
      writeLines(noquote(c("Writing: ", outfile)))
      write.csv(df, file=outfile, row.names=FALSE)
    }
  }
  df$serial <- as.integer(substring(df$serial, 16,20))  
    # removes first 15 digits in serial, then convert to integer
  colnames(df)[2] <- "date"  # "GPS_date_DDMMYYYY" to "date"
  colnames(df)[3] <- "time"  # "GPS_utc_HH:MM:SS" to "time"
  colnames(df)[4] <- "datetime"  # "GPS_YYYY..." to "datetime"
  df$datetimeUTC <- as.POSIXct(df$datetime, tz="GMT", usetz=FALSE)  
    # set time to UTC and convert to POSIXct format
  df$datetime <- format(df$datetimeUTC, tz=tz, usetz=FALSE)  # convert to EST
  df$datetime <- as.POSIXct(df$datetime, tz=tz, usetz=FALSE)  # convert to EST
  df$year <- as.numeric(strftime(df$datetime, format="%Y", usetz=FALSE))  # year
  df$date <- as.POSIXct(df$datetime, tz=tz, "%Y-%m-%d", usetz=FALSE)  # date
  df$date <- as.Date(df$datetime, tz=tz,"%Y-%m-%d", usetz=FALSE)  # only date
  if("lon" %in% colnames(df)){
     colnames(df)[which(names(df) == "lon")] <- "long"
  }
  df$long <- as.numeric(strtrim(df$long,9))  # removes N
  if (all(df$long >= 0)) df$long <- (df$long)*-1
  df$lat <- as.numeric(strtrim(df$lat,9))  # removes W
  gps_data <- read.csv("C:/Work/R/Data/BAEA/BAEA_gps_data.csv", 
    header=TRUE, as.is=TRUE, na.strings = "")
  gps_data <-  subset(gps_data, select=serial:notes)  # keeps id:notes col
  gps_data$id <- NA
  date_cols <- c("on_hand","deployed", "end_data",  "failed",  "removed",
      "recovered")
  for (i in date_cols) { 
    gps_data[,i] <- as.character(gps_data[,i]) 
    gps_data[,i] <- as.Date(gps_data[,i], "%Y%m%d")
  }
  gps_blank <- gps_data[0,]
  gps_blank[1:nrow(df),] <- NA
  gps_blank$serial <- NULL  # removes the redudant serial column 
  df <- cbind(df, gps_blank)
  if (units == "deployed") {
    deployed <- gps_data[which(!is.na(gps_data$deploy_location)),]  
    for (i in 1:nrow(deployed)) {
      record <- deployed[i,]  
    if (is.na(record$end_data)) {
      end_date <- Sys.Date() + 1
    } else {
      end_date <- record$end_data
    }
    sv <- df$serial == record$serial & df$date > record$deployed & 
      df$date < end_date
    record$id <- record$deploy_location
    record$serial <- NULL # prevents a redundant "serial.1" column
    df[sv, (ncol(df)-length(record)+1):ncol(df)] <- record[1,]
    }
    df <- subset(df, (!is.na(deploy_location)))
    }
  if (units == "reserve") {
    reserve <- gps_data
    for (i in 1:nrow(reserve)) {
      record <- reserve[i,]     
      if (is.na(record$end_data)) {
        end_date <- Sys.Date() + 1
      } else {
        end_date <- record$end_data
      }
      if (is.na(record$deployed)) {
        deploy_date <- Sys.Date() + 1
      } else {
        deploy_date <- record$deployed
      }
      sv <- df$serial == record$serial & df$date > record$on_hand & 
        df$date < deploy_date & df$date < end_date
      record$id <- record[1,"serial"]
      record$serial <- NULL  # prevents a redundant "serial.1" column
      record$deploy_location <- NA
      record$bird_ID <- NA
      record$sex <- NA
      record$deploy_seq <- NA
      df[sv, (ncol(df)-length(record)+1):ncol(df)] <- record[1,]
    }
    df <- subset(df, !is.na(id)) # still working at this point
  }
  xy <- cbind(df$long,df$lat) # 2 col for next step
  xy <- project(xy, "+proj=utm +zone=19 ellps=WGS84")  # projects to UTM Zone 19
  colnames(xy) <- c("long_utm", "lat_utm")  # name columns
  xy <- round(xy)  # rounds lat and long to zero decimal places
  df <- cbind(df, xy)  # combines lat long with data
  df$long_utm <- as.integer(df$long_utm)
  df$lat_utm <- as.integer(df$lat_utm)
  drops <- c("time", "utc", "cog", "data_voltage", "capacity", 
    "pow_voltage", "pow_timestamp", "barfcn", "dbm", "netnameasc", 
    "serv_timestamp", "id_long",
    "datetimeUTC")  # vector of columns to drop
  df<-df[ ,!(names(df) %in% drops)]  
  df<-unique(df)  # removes duplicate rows, found for 72179 2013-12-29 06:28:39 
  df<-df[,c("id",setdiff(names(df),"id"))]  # puts "id" column first 
  row.names(df)<-NULL
  return(df)
}

# CreateMove Function ----------------------------------------------------------

###  A wrapper function for creating a 'move' object
###  Usage: CreateMove(df)
###  Arguments: df = dataframe with locations
###  Returns: a 'move'object          
###  Notes: internal parameters set specifically for my BAEA GPS data
###  Blake Massey
###  2014.05.05

CreateMove <- function(df){
  suppressPackageStartupMessages(require(plyr))
  suppressPackageStartupMessages(require(circular))
  suppressPackageStartupMessages(require(move))
  df <- arrange(df, id, datetime)
  df$order <- NULL
  move <- move(x = df$long, y = df$lat, time = df$datetime, 
    proj = CRS("+proj=longlat"), animal = df$id, sensor = "GPS")
  return(move)
}  

# DownloadCTT Function ---------------------------------------------------------

###  Downloads data files from CTT's website
###  Usage: DownloadCTT(units, download)
###  Arguments: units = "deployed", "reserve", or "none"
###             download = "all" or "recent"
###  Returns: df of downloaded files
###  Notes: This is ENTIRELY DEPENDENT on my Python scripts/locations
###  Blake Massey
###  2014.05.12

DownloadCTT <- function(units="", 
                        download="recent") {
  if (units == "deployed" && download == "all") {  
    system('python C:/Work/Python/Scripts/BAEA_CTT/Import_Deployed_All.py')
    writeLines(noquote("Downloading all data for deployed units from CTT"))
  } 
  if (units == "deployed" && download == "recent") {
    system('python C:/Work/Python/Scripts/BAEA_CTT/Import_Deployed_Recent.py')
    writeLines(noquote("Downloading recent data for deployed units from CTT"))
  }
  if (units == "reserve" && download == "all") {
    system('python C:/Work/Python/Scripts/BEA_CTT/Import_Reserve_All.py')
    writeLines(noquote("Downloading all data for reserve units from CTT"))
  }
  if (units == "reserve" && download == "recent") {
    system('python C:/Work/Python/Scripts/BAEA_CTT/Import_Reserve_Recent.py')
    writeLines(noquote("Downloading recent data for reserve units from CTT"))
  }
  if (units == ""  | units == "none") {
    writeLines(noquote("Not downloading data from CTT"))
  }
}

# ImportUnits Function ---------------------------------------------------------

###  Imports previous records and merges with existing
###  Usage: ImportUnits(units, import)
###  Arguments: units = "deployed" or "reserve"
###             existing = exisiting file to merge with import
###             import = select to import previous records, default is TRUE
###  Returns: merged units file          
###  Notes: defaults are specific to my file directories and locations
###  Blake Massey
###  2014.05.12

ImportUnits <- function(units = "deployed", 
                        existing = NULL, 
                        import = TRUE) {
  if (import == TRUE) {    
    units_import <- read.csv(file = paste("C:/Work/R/Data/Telemetry/", units,
      "/", units, ".csv", sep =""), header=TRUE, stringsAsFactors=FALSE)
    units_import$datetime <- as.POSIXct(units_import$datetime,
      tz="Etc/GMT+5", usetz=FALSE) #convert to POSIXct in EST
  if (!is.null(existing)) {
    units_merge <- rbind(existing, units_import)
    units_full <- unique(units_merge)
    units_full <- units_full[with(units_full,order(id,datetime)),]
    row.names(units_full) <- NULL
    date <- Sys.Date()
    outfile <- paste("C:/Work/R/Data/Telemetry/", units, "/Archive/", units, 
      "/", units, "_", date, ".csv", sep ="")
  if (!file.exists(outfile)) {
    writeLines(noquote(paste("Merging existing and import")))
    write.csv(units_full, file=outfile, row.names=FALSE)
    writeLines(noquote(c("Writing: ", outfile, sep = ""))) 
    units_full_rewrite<-paste("C:/Work/R/Data/Telemetry/", units,
      "/", units, ".csv", sep ="")
    write.csv(units_full, file=units_full_rewrite, row.names=FALSE)  
      # rewrites import file
  } 
  } else {
    units_full <- units_import  
    date_cols<- c("date","on_hand","deployed", "end_data",  "failed",  "removed",
                "recovered")  # corresponds to list in CompileDownloads
   for (i in date_cols) { 
    units_full[,i] <- as.character(units_full[,i]) 
    units_full[,i] <- as.Date(units_full[,i], "%Y-%m-%d")
  }   
  }
  }
  if (import == FALSE) {
    writeLines(noquote(paste("Previous records NOT imported")))
  if (!is.null(existing)) {
    writeLines(noquote(paste("Coverted existing to output", sep="")))
    units_full <- existing
  }
  }
  return(units_full)
  return(units_full)
}

# Plot3DInteractive Function ---------------------------------------------------

###  Makes a 3D plot with 'rgl' package
###  Usage: Plot3DInteractive(df, x, y, z)
###  Arguments: df = dataframe with x, y, z data, Default is 'baea' 
###             x  = column name of x (longitude) locations
###             y  = column name of y (latitude) locations
###             z  = column name of z (elevation) locations
###  Returns: 3D plot from 'rgl'package         
###  Notes: can be any z, y, z, but labels are set for lat, long, elev data
###  Blake Massey
###  2014.05.05

Plot3DInteractive<-function(df = baea, 
                            x = "x", 
                            y = "y", 
                            z = "z"){
  suppressPackageStartupMessages(require(rgl))
  df$x <- df[,x]
  df$y <- df[,y]
  df$z <- df[,z]
  xyz <- cbind(df$x, df$y, df$z)  # df of x, y, z data
  colnames(xyz) <- c("x", "y", "z")
  open3d()  # open an RGL window
  bg3d("white")  # make the background white
  plot3d(xyz, aspect=c(1, 1, .5), col = "darkblue",
         xlab = "longitude", ylab = "latitude", zlab = "elevation(m)")
  lines3d(xyz, color="red")
}

# PlotLocationCount Function ---------------------------------------------------

###  Plot a linegraph of daily locations 
###  Usage: PlotLocationsCount(df, id, individual, color_factor, start,
###    end, breaks, addsolartimes, wrap)
###  Arguments: df = dataframe with id, datetime
###             id = column name of unique identifier
###             individual = individual/s (from id column) to keep,
###               format should be c(id, id), default is all          
###             color_factor = factor to determine location point color
###             start = filter start date, default is 1970-01-01   
###             end = filter end date, default is current date 
###             breaks = breaks on x-axis (i.e., "7 days")
###             tz = uses this timezome. Default "Etc/GMT+5". 
###             wrap = wraps ribbon of panels into 2d. Default TRUE. 
###  Returns: daily locations over a range of dates           
###  Notes: Using color_factor on "behavior"
###  Blake Massey
###  2014.05.20

PlotLocationCount <- function(df = df, 
                              id = "id", 
                              individual = "",
                              color_factor = NULL, 
                              start = NULL, 
                              end = NULL, 
                              breaks = "14 days",
                              tz = "Etc/GMT+5", 
                              wrap = FALSE) {
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(lubridate)) 
  suppressPackageStartupMessages(require(plyr)) 
  suppressPackageStartupMessages(require(scales))
  if(!is.null(color_factor)) {
    cf_name <- color_factor
    df$color_factor <- df[,color_factor]
  }
  if ( ! (individual == "" || individual == "all")){  
    df = df[which(df[,id] == individual), ]  # to extract individuals  
  }
  if(!is.null(start)) {
    start = as.Date(start)
  } else {
    start = min(df$date)
  }
  if(!is.null(end)) {
    end = as.Date(end)
  } else {
    end = max(df$date)
  }
  limits_x = c(start, end)
  sumstats <- SummarizeDailyLocations(df=df)
  p <- ggplot(data=sumstats, aes(x=date, y=total_loc)) +  
  labs(title = "Daily Locations", x="Date", y="Locations") + 
  scale_x_date(breaks=date_breaks(breaks), labels=date_format("%m/%d"), 
  limits=limits_x) 
  if(!is.null(color_factor)) {
    p <-  p + geom_line(aes(group=factor(id), colour = factor(id))) +
    geom_point(aes(group=factor(id), colour = factor(id))) +
    labs(shape=cf_name, colour=cf_name) 
  } else {
    p <- p + geom_line() + geom_point()
  }
  p <- p + theme_bw() + theme(plot.title = element_text(size = 22)) +
  theme(text = element_text(size = 18, colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black")) +    
  theme(text = element_text(size = 20, colour = "black")) +
  theme(axis.text = element_text(colour = "black"))  
  if (wrap == TRUE) {
  p + facet_wrap(~ id, scales= "free_x")
  } else {
  p
  }
}

# PlotLocationSunriseSunset Function -------------------------------------------

###  Plot diagram of locations in relation to sunrise, sunset, and solarnoon 
###  Usage: PlotLocationSunriseSunset(df, id, individual, color_factor, start,
###    end, breaks, addsolartimes, wrap)
###  Arguments: df = dataframe with id, lat, long, datetime
###             by = column name to use to split, analyze, and merge data, 
###               default is "id" 
###             individual = individual/s (from id column) to keep,
###               format should be c(id, id), default is all          
###             color_factor = factor to determine location point color
###             pal = color palette for CreateVarsColors(), default is NULL
###             b_pal = color palette from RColorBrewer for CreateVarsColors(),
###               default is "Set1"
###             start = filter start date, default is 1970-01-01   
###             end = filter end date, default is current date 
###             breaks = breaks on x-axis (i.e., "7 days")
###             tz = uses this timezome. Default "Etc/GMT+5".
###             addsolartimes = runs this function on data first. Default FALSE. 
###             wrap = wraps ribbon of panels into 2d. Default TRUE. 
###  Returns: plot of locations over a range of dates           
###  Notes: Using color_factor on "behavior"
###  Blake Massey
###  2014.10.06

PlotLocationSunriseSunset <- function(df = df, 
                                      by = "id", 
                                      individual = "",
                                      color_factor = NULL,
                                      pal = NULL,
                                      b_pal = "Set1",
                                      start = NULL, 
                                      end = NULL, 
                                      breaks = "1 days",
                                      tz = "Etc/GMT+5", 
                                      addsolartimes = FALSE, 
                                      wrap = TRUE) {

  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(scales))
  suppressPackageStartupMessages(require(lubridate)) 
  source('C:/Work/R/Functions/baea.R')
  source('C:/Work/R/Functions/gen.R')
  df=df
  df$by <- df[,by]
  if(!is.null(color_factor)) {
    cf_name <- color_factor
    df$color_factor<-df[,color_factor]
    by_colors <- CreateColorsByAny(by=color_factor, df=df, output=TRUE, pal=pal, 
    b_pal=b_pal)
  } else {
    ifelse(is.null(by),  df$by <- "all", df$by <- df[,by])     
    by_colors <- CreateColorsByAny(by=by, df=df, output=TRUE, pal=pal, 
    b_pal=b_pal)
  }
  if ( ! (individual == "" || individual == "all" || is.null(individual))){  
    df = df[which(df[,"id"] == individual), ]  # to extract individuals  
  }
  if (addsolartimes == TRUE) {
    df <- AddSolarTimes(df = df, by = by,  tz = tz)
    df$by <- df[,by]  # AddSolorTimes removes df$by  
  }
  df$loc_time <- format(df$datetime, format = "%H:%M:%S")
  df$loc_time <- as.POSIXct(df$loc_time, tz=tz, format = "%H:%M:%S")
  df$sunrise <- format(df$sunrise, format = "%H:%M:%S")
  df$sunrise <- as.POSIXct(df$sunrise, tz=tz, format = "%H:%M:%S")
  df$solarnoon <- format(df$solarnoon, format = "%H:%M:%S")
  df$solarnoon <- as.POSIXct(df$solarnoon, tz=tz, format ="%H:%M:%S")
  df$sunset <- format(df$sunset, format = "%H:%M:%S")
  df$sunset <- as.POSIXct(df$sunset, tz=tz, format = "%H:%M:%S")
  if(!is.null(start)) {
    start = as.POSIXct(start)
  } else {
    start = min(df$datetime)
  }
  if(!is.null(end)) {
    end = as.POSIXct(end)
  } else {
    end = max(df$datetime)
  }
  limits_x = c(start, end)
  p <- ggplot(data = df) +
    geom_line(aes(datetime, solarnoon), colour = "red", size = 2) +
    geom_line(aes(datetime, sunset), colour = "orange" , size = 2) +
    geom_line(aes(datetime, sunrise), colour = "orange", size = 2) +
    labs(title = "Locations in Relation to Sunrise and Sunset", 
    x="Date", y="Time") +
    scale_y_datetime(breaks=date_breaks("1 hour"), labels=date_format("%H")) +
    scale_x_datetime(breaks=date_breaks(breaks), labels=date_format("%m/%d"),
    limits=limits_x) 
  if(!is.null(color_factor)) {
    p <-  p + geom_point(aes(datetime, loc_time, colour = factor(color_factor), 
     shape = factor(color_factor)), size = 3) +
      labs(shape=cf_name, colour=cf_name) +  scale_fill_brewer(palette="Set1")
  } else {
    p <- p + geom_point(aes(datetime, loc_time))
  }
  p <- p + theme_bw() + theme(plot.title = element_text(size = 22)) +
    theme(text = element_text(size = 18, colour = "black")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black")) +    
    theme(text = element_text(size = 20, colour = "black")) +
    theme(axis.text = element_text(colour = "black"))  +
    theme(panel.background = element_rect(fill = "gray90")) +
    theme(panel.grid.major = element_line(color = "gray80")) +
    theme(panel.grid.minor = element_line(color = "gray80")) +
    scale_colour_manual(values = by_colors)
  if (wrap == TRUE) {
  p + facet_wrap(~ by)
  } else {
  p
  }
}

# PlotMove Function ------------------------------------------------------------

###  Wrapper for plotting a 'move' object
###  Usage: PlotMove(move)
###  Arguments: move = move object
###  Returns: plot of move object         
###  Notes: only work on move objects
###  Blake Massey
###  2014.05.05

PlotMove<-function(move = move){
  suppressPackageStartupMessages(require(move)) 
  units <- nlevels(idData(move))
  plot(move, type="o", col=seq(1,units,by=1), lwd=2, pch=20,
    xlab="Longitude", ylab="Latitude")
}

# PlotOvernightDistances Function ----------------------------------------------

###  Scatterplot of overnight distances. Non-sequential days are removed.
###  Usage: PlotOvernightDistances(df, id, ylim, point_size)
###  Arguments: df = dataframe of locations with "date", "last", "step_length"
###             id = column name of unique identifier
###             ylim = y-axis limit 
###             point_size = point size for time datapoints 
###  Returns: scatterplot          
###  Notes: 
###  Blake Massey
###  2014.05.05

PlotOvernightDistances<-function(df = baea, 
                                 individual = "",  
                                 ylim = NULL, 
                                 point_size = 3){
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(lubridate))
  suppressPackageStartupMessages(require(scales))
  df <- RemoveNonsequentialDays(df = df) 
  if (individual == ""|individual == "all"){
    df = df[which(df$last == "Last"),] 
  } else {
    df = df[which(df$last == "Last" & df$id == individual),]  
  }
  if (length(unique(df$id)) == 1) {
    title = paste("Overnight Distance by Date: Unit ", unique(df$id), sep = "")
    if(ylim == "" || is.null(ylim)){
      g<-ggplot(df, aes(x = date, y = step_length)) 
      g + geom_point(size = point_size, colour="dark blue") +
        theme(text = element_text(size = 20, colour="black")) +
        theme(axis.text.y = element_text(colour="black")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1,
          colour= "black")) + 
        scale_x_date(labels = date_format("%m/%d"), breaks="7 days") +
        labs(title=title, x="Date", y="Distance (m)", text="black") +
        guides(shape = guide_legend("ID #"))
    } else {
      title = paste("Plot of Overnight Distance by Date (y axis limited to ",
        ylim, "m)", sep = "")
      g <- ggplot(df, aes(x = date, y = step_length)) 
      g + ylim(0, ylim) +
        geom_point(size = point_size, colour = "darl blue") +
        theme(text = element_text(size = 20, colour = "black")) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1,
          colour= "black")) +  
        scale_x_date(labels = date_format("%m/%d"), breaks = "7 days") 
      labs(title = title, x = "Date", y = "Distance (m)", text="black") +
        guides(shape = guide_legend("ID #"))
    }
  } else {
    if(ylim == "" || is.null(ylim)){
      g <- ggplot(df, aes(x = date, y = step_length, colour = factor(id))) 
      g + geom_point(size = point_size) +
        theme(text = element_text(size = 20, colour = "black")) +
        theme(axis.text.y = element_text(colour = "black")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, 
          colour = "black")) + 
        scale_x_date(labels = date_format("%m/%d"), breaks = "7 days") +
        labs(title = "Overnight Distances by Date", x = "Date", 
             y = "Distance (m)", text="black") +
        guides(colour = guide_legend("ID #"))
    } else {
      title = paste("Plot of Overnight Distances by Date (y axis limited to ",
        ylim, "m)", sep = "")
      g <- ggplot(df, aes(x = date, y = step_length, colour = factor(id))) 
      g + ylim(0, ylim) +
        geom_point(size = point_size) +
        theme(text = element_text(size = 20, colour="black")) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1,
          colour= "black")) + 
        scale_x_date(labels = date_format("%m/%d"), breaks = "7 days") 
      labs(title="Overnight Distance by Date", x = "Date", y = "Distance (m)", 
        text = "black") +
        guides(colour = guide_legend("ID #")) 
    }
  }
}

# PlotOvernightDistancesECDF Function ------------------------------------------

###  Empirical cumulative distribution; non-sequential days removed
###  Usage: PlotOvernightDistancesECDF(df, id, ylim, point_size)
###  Arguments: df = dataframe of locations with "date", "last", "step_length"
###             id = column name of unique identifier
###             xlim = x-axis limit, default is NULL 
###             line_size = point size for time datapoints
###  Returns: scatterplot          
###  Notes: requires RemoveNonsequentialDays function
###  Blake Massey
###  2014.05.05

PlotOvernightDistancesECDF <- function(df = baea, 
                                       id = NULL, 
                                       xlim = NULL,
                                       line_size = 2){
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(lubridate))
  suppressPackageStartupMessages(require(scales))
  df <- RemoveNonsequentialDays(df = df)
  if (id == "" || is.null(id)){  
    df = df[which(df$last == "Last"), ] 
  } else {
    df = df[which(df$last == "Last" & df$id == id),]  
  }  
  if (is.null(xlim)){
    g <- ggplot(df, aes(step_length, colour = factor(id))) 
    g + stat_ecdf(size = line_size) +
      theme(text = element_text(size = 20, colour = "black")) +
      theme(axis.text = element_text(colour = "black")) +
      labs(title = "ECDF of Overnight Distances", x = "Distance (m)",
        y = "Proportion", text = "black") +
      guides(colour = guide_legend("ID #"))   
  } else {
    title = paste("Plot of Overnight Distance by Date (x axis limited to ", 
                xlim, "m)", sep = "")
    g <- ggplot(df, aes(step_length, colour = factor(id))) 
    g + stat_ecdf() +
      xlim(0, xlim) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"), 
        axis.text.y = element_text(colour = "black")) +
      theme(text = element_text(size = 20, colour = "black")) +
      theme(axis.text = element_text(colour = "black")) +
      labs(title=title, x="Distance (m)", y="Proportion",
           text="black") +
      guides(colour = guide_legend("ID #"))
  }
}

# RemoveNonsequentialDays Function --------------------------------------------- 

###  Removes days that do not have GPS data for the next day
###  Usage: RemoveNonsequentialDays(df)
###  Arguments:  df = database of location, default is "baea"
###  Returns: df
###  Notes: used in PlotOvernightDistancesECD()
###  Blake Massey
###  2014.05.05

RemoveNonsequentialDays<-function (df = baea){
  suppressPackageStartupMessages(require(plyr))
  suppressPackageStartupMessages(require(maptools))
  sumstats <- ddply(df, .(id, date), summarize, date = as.Date(unique(date)),
  total_loc = length(deploy_seq))
  nextday <- function(data = data){
    out <- sapply(2:nrow(data), function(i){data$date[i] - data$date[i-1]})
    next_day <- c(out, NA) 
    return(next_day)
  }
  list <- by(sumstats, sumstats$id, function(x) nextday(x))  # need a list
  sumstats$nextdayGPS <- unlist(list)  # list into array, combines with sumstats  
  df <- merge(df, sumstats, by = c("id","date"), all.x = TRUE)
  df <- subset(df, nextdayGPS == 1)
  df
}

# SummarizeDailyLocations Function ---------------------------------------------

###  Gets the daily number of locations for a location dataset
###  Usage: SummarizeDailyLocations(df)
###  Arguments: df = dataframe of locations, needs "id" and "date" columns
###  Returns: dataframe of total daily location, including 0 for days w/o data          
###  Notes: used in PlotDailyLocationsCount() 
###  Blake Massey
###  2014.05.21

SummarizeDailyLocations <-function (df = df){
  sumstats <- ddply(df, .(id, date), summarize, date = as.Date(unique(date)), 
    total_loc = length(id))
  sumstats_filled <- sumstats[0,]
  unique_id<- unique(df$id)
  for (i in unique_id) {
  df<-subset(sumstats, id == i)
  date_min <- min(df$date)
  date_max <- max(df$date)
  all_dates <- seq(date_min, date_max, by="day")
  all_dates_frame <- data.frame(list(date=all_dates))
  merged_frame <- merge(all_dates_frame, df, all=T)
  merged_frame$id <- i
  merged_frame$total_loc[which(is.na(merged_frame$total_loc))] <- 0
  sumstats_filled<-rbind(sumstats_filled,merged_frame)
  }
  return(sumstats_filled)
}

# SummarizeLocations Function --------------------------------------------------

###  Calculates Summarize stats for altitude, speed, moving and total locations
###  Usage: SummarizeLocations(df)
###  Arguments: df = dataframe with locations
###             pdf = whether to write PDF, default is FALSE
###  Returns: table and PDF (optional)         
###  Notes: default outputs are specific to my file directories and locations
###  Blake Massey
###  2014.05.05

SummarizeLocations <- function(df = df, 
                               pdf = FALSE) {
  suppressPackageStartupMessages(require(gridExtra))
  suppressPackageStartupMessages(require(plyr))
  sumstats <- ddply(df, .(id, date), summarize, date = as.Date(unique(date)), 
    max_alt = round(max(alt), 0), min_alt = round(min(alt), 0),
    max_speed = round(max(speed), 0), moving_loc = sum(speed>2),
    total_loc = length(id))
  sumstats_filled <- sumstats[0,]
  unique_id<- unique(df$id)
  for (i in unique_id) {
  df<-subset(sumstats, id == i)
  date_min <- min(df$date)
  date_max <- max(df$date)
  all_dates <- seq(date_min, date_max, by="day")
  all_dates_frame <- data.frame(list(date=all_dates))
  merged_frame <- merge(all_dates_frame, df, all=T)
  merged_frame$id <- i
  merged_frame$total_loc[which(is.na(merged_frame$total_loc))] <- 0
  sumstats_filled<-rbind(sumstats_filled,merged_frame)
  } 
  sumstats <-sumstats_filled
  write.csv(sumstats, "C:/Work/R/Data/Output/Summary Stats.csv", 
    row.names = FALSE) 
  if (pdf == TRUE){
  maxrow = 30 
  npages = ceiling(nrow(sumstats)/maxrow); #multi-page pdf 
  pdf("C:/Users/Blake/Desktop/Summary Stats.pdf", height=11, width=8.5) 
  for (i in 1:npages){idx = seq(1+((i-1)*maxrow), i*maxrow); grid.newpage();
    grid.table(sumstats[idx, ])}
  dev.off();
  } 
  return(sumstats)
}

# SummarizeOvernightDistances Function -----------------------------------------

###  Summarizes overnight distances
###  Usage: SummarizeOvernightDistances(df, individual, cuts)
###  Arguments: df = dataframe
###             individual = individual (from "id" column) to keep, optional 
###             cuts = proportional quantile cuts, default = 0.05
###  Returns: dataframe of quantiles       
###  Notes: requires RemoveNonsequentialDays function
###  Blake Massey
###  2014.05.05

SummarizeOvernightDistances <- function(df = df, 
                                        individual = NULL, 
                                        cuts = (0.05)) {
  df <- RemoveNonsequentialDays(df = df)
  if (individual == "" || is.null(individual)){
    df = df[which(df$last == "Last"),] 
  } else {
    df = df[which(df$last == "Last" && df$id == individual),]  
  }  
  quantiles <- quantile(df$step_length, seq(0, 1, cuts))
  return(quantiles)
}
