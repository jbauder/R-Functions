# Functions more specific to Bald Eagle behavior and data analysis 
# ------------------------------------------------------------------------------

# AddCruiseBehavior Function ---------------------------------------------------

###  Finds cruise behavior that meet given threshold parameters
###  Usage: AddCruiseBehavior(df, alt_abv_grd, min_speed)
###  Arguments: df = dataframe 
###             alt_abv_grd = distance (in meters) above ground
###             min_speed = minimum speed of bird
###  Returns: dataframe with behavior column that has "cruise"
###  Notes: should run AddNestBehavior, AddRoostBehavior, AddCruiseBehavior 
###    prior to this function
###  Blake Massey
###  2014.09.30

AddCruiseBehavior <- function(df = df, 
                              min_agl = 200,
                              min_speed = 5) {
  df <- df
  df$behavior <- ifelse(df$agl >= min_agl & df$speed >= min_speed & 
    is.na(df$behavior), "cruise", df$behavior)
  return(df)
}

# AddFlightData Function -------------------------------------------------------

###  Finds location data that meet given threshold parameters for flights, which
###    include the start point, mid-flight, and final location of a flight
###  Usage: AddFlightsData(df, alt_abv_grd, min_speed)
###  Arguments: df = dataframe 
###             min_speed = minimum speed of bird
###             min_step_length = distance (in meters) between locations
###             max_time_step = maximum time between locations
###  Returns: dataframe with columns for flight_index, flight_step, and
###    flight_lenght
###  Notes: should run AddLandscapeValues() prior to this function 
###  Blake Massey
###  2014.11.27

AddFlightData <- function(df, 
                          min_speed = 5,
                          min_step_length = 50,
                          max_step_time = 20){
  df_org <- df
  df <- df
  df$flight <- ifelse(df$speed > min_speed, TRUE, FALSE)
  df$movement <- ifelse(df$step_length > min_step_length, TRUE, FALSE)
  df$within_window <- ifelse(df$step_time < max_step_time, TRUE , FALSE)
  moved <- ifelse(df$step_length > min_step_length, TRUE, FALSE)
  df$moved <- c(NA, moved[-length(moved)]) 
  df$mid_flight <- ifelse(df$moved == TRUE & df$movement == TRUE & 
    df$flight == TRUE, TRUE, FALSE)
  df$index <- seq.int(nrow(df))
  before <- which(df$mid_flight == TRUE) - 1
  mid <- which(df$mid_flight == TRUE)
  after <- which(df$mid_flight == TRUE) + 1
  complete <- sort.int(unique(c(before, mid, after)))
  flights <- df[complete,c("id", "datetime", "index", "mid_flight")]
  index1 <- c(which(!diff(flights$index)==1), nrow(flights))
  index2 <- c(1, index1 + 1)  # index of first records of sequential groups 
  index3 <- diff(index2)  # length of each seq group
  index4 <- rep(seq(1,length(index3), by=1), times=index3)
  
  flights$flight_index <- index4
  flights$flight_step <- sequence(rle(flights$flight_index)$lengths)
  flights$flight_index <- factor(index4)  
  suppressPackageStartupMessages(require(plyr))
  flights <- ddply(flights, .(flight_index, id), transform, 
    flight_length=length(flight_step))
  df2 <- suppressMessages(plyr::join(df_org, flights))
  df2$index <- NULL
  return(df2)
}

# AddNestBehavior Function -----------------------------------------------------

###  Finds nest attendance behavior that meet given threshold parameters
###  Usage: AddNestBehavior(df, distance_to_nest)
###  Arguments: df = dataframe 
###             distance_to_nest = distance (in meters) from location to nest to
###               assign a location as "nest" behavior
###  Returns: dataframe with behavior column that has "nest"   
###  Notes: need to run AddNestData() prior to running this
###  Blake Massey
###  2014.06.11

AddNestBehavior <- function(df = df, 
                            distance_to_nest = 50) {
  df<-df
  df$behavior <- NA
  df$behavior <- ifelse(df$nest_dist <= distance_to_nest, "nest", df$behavior)
  return(df)
}

# AddNestConDist Function ------------------------------------------------------

###  Adds nest and conspecific distances to dataframe
###  Usage: AddNestConDist(df, nests_con_dist)
###  Arguments: df = dataframe of locations, must have coordinate columns
###             nest_con_dist = list of nest and conspecific distance rasters
###  Returns: dataframe with "nest_dist" and "con_dist" columns
###  Notes: The nest_con_dist has to have structure of list[[year]][[nest_id]]
###  Blake Massey
###  2014.12.11

AddNestConDist<- function(df, 
                          nest_con_dist = nest_con_dist) {
  df <- df
  nest_con_dist <- nest_con_dist
  AddDistances <- function(df){
    df <- df
    coordinates(df) <- cbind(df$long_utm, df$lat_utm)
    proj4string(df) <- CRS("+proj=utm +zone=19 +datum=NAD83")
    nest_id <- unique(df$nest_id)
    year <- unique(df$year)
    dists <- nest_con_dist[[as.character(year)]][[nest_id]]
    df$nest_dist <- round(as.vector(unlist(extract(dists, df, layer=1, 
      nl=1))))
    df$con_dist <-round(as.vector(unlist(extract(dists, df, layer=2, 
      nl=1))))
    output <- as.data.frame(df)
    return(output)
  }
  output <- ddply(df, .(year, id), AddDistances)
  output$x <- NULL
  output$y <- NULL
  return(output)
}

# AddNestData Function ---------------------------------------------------------

###  Adds nest locations to each record based on "date" and "id"
###  Usage: AddNestData(df, nests)
###  Arguments: df = dataframe 
###             nests = .csv file of nests associated with BAEA location data
###  Returns: dataframe with merged nest data from "nests" file   
###  Notes: 
###  Blake Massey
###  2014.12.10

AddNestData <- function(df = df, 
                        nests = "C:/Work/R/Data/BAEA/BAEA_nests.csv") {
  source('C:/Work/R/Functions/baea.R')
  suppressPackageStartupMessages(require(rgdal))
  suppressPackageStartupMessages(require(plyr))
  df <- df
  nests <- read.csv(nests, header=TRUE, stringsAsFactors=FALSE)
  date_cols<- c("use_start_date", "use_end_date", "clutch_initiation", 
    "breeding_end_date")
  for (i in date_cols) { 
    nests[,i] <- as.character(nests[,i]) 
    nests[,i] <- as.Date(nests[,i], "%Y%m%d")
  }
  nests_blank <- nests[0,]
  nests_blank[1:nrow(df),] <- NA
  df<-cbind(df, nests_blank)
  for (i in 1:nrow(nests)) {
    nest <- nests[i,] 
  if (is.na(nest$use_end_date)) {
    use_end_date <- Sys.Date() + 1
  } else {
    use_end_date <- nest$use_end_date
  }
  sv <- df$id == nest$eagle_id & df$date >= nest$use_start_date & 
    df$date <= use_end_date
  df[sv, (ncol(df)-length(nest)+1):ncol(df)] <- nest[1,]
  }
  drops <- c("eagle_id", "nest_long", "nest_lat", "active_2013", "active_2014", 
    "active_2015", "active_2016", "active_2017", "use_start_date", 
    "use_end_date", "clutch_initiation", "breeding_end_date", "clutch_size",
    "num_fledge", "nest_notes")  # vector of columns to drop
  df <- df[ ,!(names(df) %in% drops)]
  df$nest_angle <- CalculateAngleToPoint(df, "long_utm", "lat_utm", 
    "nest_long_utm", "nest_lat_utm")
#   xy <- cbind(df$nest_long,df$nest_lat) # 2 col for next step
#   xy <- project(xy, "+proj=utm +zone=19 ellps=WGS84")  # projects to UTM N19
#   colnames(xy) <- c("nest_long_utm", "nest_lat_utm")  # name columns
#   xy <- round(xy)  # rounds lat and long to zero decimal places
#   df <- cbind(df, xy)  # combines lat long with data
#   df<-adply(df, 1, transform, dist_nest = as.integer(sqrt(sum((c(long_utm,
#     lat_utm) - c(nest_long_utm, nest_lat_utm))^2)))) # Pythagorean theorem  
  return(df)
}

# AddPerchBehavior Function ----------------------------------------------------

###  Finds perch behavior that meet given threshold parameters
###  Usage: AddPerchBehavior(df, max_speed)
###  Arguments: df = dataframe
###             max_speed = max_speed to assign 
###  Returns: dataframe with behavior column that has "perch"   
###  Notes:
###  Blake Massey
###  2014.12.17

AddPerchBehavior <- function(df = df, 
                             max_speed = 5) {
  df<-df
  df$behavior <- NA
  df$behavior <- ifelse(df$speed < max_speed, "perch", df$behavior)
  return(df)
}

# AddRoostBehavior Function ----------------------------------------------------

###  Finds roost arrivals and departures that meet given threshold parameters
###  Usage: AddRoostBehavior(df, overnight_distance_threshold, 
###    at_roost_distance_threshold, daily_location_threshold)
###  Arguments: df = dataframe 
###             default_tz = used in IfElseTimedateNA/Compare functions
###             tz = timezone, default is "Etc/GMT+5"
###             overnight_distance_threshold = max overnight distance
###             at_roost_distance_threshold = max distance away from roost
###             depart_timediff_max = max diff between start of day and depart
###             arrive_timediff_max = max diff between end of day and arrival
###  Returns: dataframe with roost column that has arrive, depart, and roost   
###  Notes: automatically makes sure that data exists for the following day,
###    automatically checks that there are at least 7 locations in the 
###    first/last two hours of the day for departure/arrive
###  Blake Massey
###  2014.06.11

AddRoostBehavior <- function(df = baea, 
                             default_tz = "America/New_York",
                             tz = "Etc/GMT+5",
                             overnight_distance_threshold = 100,
                             at_roost_distance_threshold = 50,
                             depart_timediff_max = 1000, 
                             arrive_timediff_max = 1000){
  suppressPackageStartupMessages(require(lubridate))  # needed for tz function
  suppressPackageStartupMessages(require(plyr))
  suppressPackageStartupMessages(require(maptools))
  suppressPackageStartupMessages(require(lubridate))
  df$time_after_start <- as.integer(difftime(df$datetime, 
    df$hr_before_sunrise, tz=tz, units = ("mins"))) 
  df$time_before_end <- as.integer(difftime(df$hr_after_sunset, df$datetime, 
    tz=tz, units = ("mins")))
  df$two_hr_after_sunrise <- df$hr_before_sunrise + hours(2)
  df$two_hr_before_sunset <- df$hr_after_sunset - hours(2)
  df$sunrise_window_loc <- df$datetime <= df$two_hr_after_sunrise
  df$sunset_window_loc <- df$datetime >= df$two_hr_before_sunset
  sumstats<-ddply(df, .(id, date), summarize, 
              date = as.Date(unique(date)), total_loc = length(deploy_seq),
              am_loc = sum(sunrise_window_loc, na.rm=TRUE),
              pm_loc = sum(sunset_window_loc, na.rm=TRUE))
  nextday <- function(data = data){
    out <- sapply(2:nrow(data),function(i){data$date[i] - data$date[i-1]})
    next_day <- c(out, NA) 
    next_day
  }
  nextamloc <- function(data = data){
    out <- sapply(1:nrow(data),function(i) { data$am_loc[i+1] })
    next_am_loc <- c(out) 
    next_am_loc
  }
  list <- by(sumstats, sumstats$id, function(x) nextday(x))  # makes list
  sumstats$nextdayGPS <- unlist(list)  
  list <- by(sumstats, sumstats$id, function(x) nextamloc(x))  # makes list
  sumstats$next_am_loc <- unlist(list) 
  # At this point, sumstats has: "id", "date", "total_locs", "nextdayGPS", 
  # "am_loc", "pm_loc", and "next_am_loc"  
  # This has all the data needed to cull by nextDayGPS, and am/pm locations  
  df<-merge(df, sumstats, by = c("id", "date"), all.x = TRUE)
  last_threshold<- function(data, threshold=overnight_distance_threshold){ 
    subset (data, last == "Last" & step_length <= threshold & nextdayGPS == 1 &
    pm_loc >= 7 & next_am_loc >= 7, select=c("id", "date"))
  } 
  last_roost_confirmed <- last_threshold(data = df, 
    threshold=overnight_distance_threshold) 
  # This culls by overnight segment length, if next day has GPS locations, if
  # there are at least 7 locations within an hour of either side of sunset 
  # and if there are at least 7 locations within an hour of either side of 
  # sunrise on the following morning.
  row.names(last_roost_confirmed) <- NULL  # housekeeping
  first_roost_confirmed <- adply(last_roost_confirmed, 1, transform, 
    date = date+1)  # the mornings after "last_roost_confirmed" dates 
  roost_arrival_filtered <- join (df, last_roost_confirmed, type="inner") 
    # to "confirm" arrival based on overnight distance
  roost_departure_filtered <- join (df, first_roost_confirmed, 
    type = "inner")  # to "confirm" departure based on overnight distance  
  threshold_dist <- function(x, threshold=at_roost_distance_threshold) {
    x>threshold
  }  
  # This function sets the distance a location can still be considered at roost
  # based on its distance to the last and first locations
  depart <- ddply(roost_departure_filtered, .(date, id), function(x)   
    x[(Position(threshold_dist, x$dist_first, right = FALSE, nomatch = NULL) 
       - 1), c("id","date", "datetime")])
  arrive <- ddply(roost_arrival_filtered, .(date, id), function(x)   
    x[(Position(threshold_dist, x$dist_last, right = TRUE, nomatch = NULL) + 1), 
      c("id", "date", "datetime")])
  arrive$arr_dist_threshold_datetime <- arrive$datetime
  depart$dep_dist_threshold_datetime <- depart$datetime
  arrive$datetime <- NULL
  depart$datetime <- NULL
  threshold_time_depart <- function(x, threshold=depart_timediff_max) {
    x>threshold
  } 
  threshold_time_arrive <- function(x, threshold=arrive_timediff_max) {
    x>threshold
  } 
  depart_max <- ddply(roost_departure_filtered, .(date, id), 
    function(x) x[nrow(x), c("id","date", "datetime")])
  arrive_max <- ddply(roost_arrival_filtered, .(date, id), 
    function(x) x[1, c("id","date", "datetime")])
  depart_max$max_datetime <- depart_max$datetime
  arrive_max$max_datetime <- arrive_max$datetime
  depart_max$datetime <- NULL
  arrive_max$datetime <- NULL
  # max_datetime is the first/last record in the day
  depart_threshold <- ddply(roost_departure_filtered, .(date, id), function(x)   
    x[(Position(threshold_time_depart, x$time_after_start, right = FALSE, 
    nomatch = NULL) - 1), c("id","date", "datetime")])
  arrive_threshold <- ddply(roost_arrival_filtered, .(date, id), function(x)   
    x[(Position(threshold_time_arrive, x$time_before_end, right = TRUE, 
    nomatch = NULL) + 1), c("id","date", "datetime")])
  depart_threshold$threshold_datetime <- depart_threshold$datetime
  arrive_threshold$threshold_datetime <- arrive_threshold$datetime
  depart_threshold$datetime <- NULL
  arrive_threshold$datetime <- NULL
  # threshold_dateime is the first record within the threshold
  depart_max <- merge(depart_max, depart_threshold, by = c("date", "id"), 
    all.x= TRUE)
  arrive_max <- merge(arrive_max, arrive_threshold, by = c("date", "id"), 
    all.x= TRUE)
  rm(arrive_threshold, depart_threshold)
  depart_max <- IfElseTimedateNA(df=depart_max, col1="threshold_datetime",
    col2="max_datetime", result="max_threshold", default_tz=default_tz, tz=tz)
  arrive_max <- IfElseTimedateNA(df=arrive_max, col1="threshold_datetime",
    col2="max_datetime", result="max_threshold", default_tz=default_tz, tz=tz)
  depart <- merge(depart, depart_max, by = c("date", "id"), all.x= TRUE)
  arrive <- merge(arrive, arrive_max, by = c("date", "id"), all.x= TRUE)
  depart <- IfElseTimedateCompare(df=depart, col1="dep_dist_threshold_datetime",
    sign="<", col2="max_threshold", result="dep_datetime", 
    default_tz=default_tz, tz=tz)
  arrive <- IfElseTimedateCompare(df=arrive, col1="max_threshold", sign=">",
    col2="arr_dist_threshold_datetime", result="arr_datetime", 
    default_tz=default_tz, tz=tz)  
  roost_times <- merge(depart, arrive, by = c("date", "id"), all = TRUE) 
  roost_times <- subset(roost_times, select=c("date", "id", "dep_datetime", 
    "dep_dist_threshold_datetime", "arr_datetime", 
    "arr_dist_threshold_datetime"))
  df <- merge(df, roost_times, by = c("date", "id"), all= TRUE)
  df <- df[with(df, order(id, date, datetime)), ]                    
  df$loaf <- NA
  df$loaf <- ifelse(df$datetime <= df$dep_dist_threshold_datetime, "loaf", NA)
  df$loaf <- ifelse(is.na(df$loaf) & df$datetime >= 
    df$arr_dist_threshold_datetime, "loaf", df$loaf)
  df$depart <- ifelse(df$datetime == df$dep_datetime, "depart", NA)  
  df$arrive <- ifelse(df$datetime == df$arr_datetime, "arrive", NA)
  df$depart <- ifelse(is.na(df$depart) & !is.na(df$dep_datetime) &
                         df$datetime < df$dep_datetime, "roost", df$depart)
  df$arrive <- ifelse(is.na(df$arrive) & !is.na(df$arr_datetime) &
                         df$datetime > df$arr_datetime, "roost", df$arrive)
  df$roost <- ifelse(is.na(df$depart),df$arrive,df$depart)
  df$roost_loaf <- ifelse(!is.na(df$roost), df$roost, 
    df$loaf)
  if (!("behavior" %in% colnames(df))) df$behavior<-NA
  df$behavior <- ifelse(!is.na(df$roost_loaf) & is.na(df$behavior), 
    df$roost_loaf, df$behavior)
  drops <- c("arrive", "depart", "total_loc", "nextdayGPS", "dep_datetime",
    "two_hr_after_sunrise", "two_hr_before_sunset", "sunrise_window_loc",
    "arr_datetime", "sunset_window_loc", "am_loc", "pm_loc", "next_am_loc",
    "roost", "loaf", "roost_loaf" ,"dep_dist_threshold_datetime", 
    "arr_dist_threshold_datetime")
  df<-df[ ,!(names(df) %in% drops)] 
  row.names(df) <- NULL
  return(df)
}

# AddTimeStepProportion Function -----------------------------------------------

###  Adds time_steps, day_min, time_after_start, and time_proportion data
###  Usage: AddTimeStepProportion(df, time_step="15 min", tz = "Etc/GMT+5")
###  Arguments: df = dataframe 
###             time_step = time step length, based on lubriate times 
###             tz = timezone for dataset 
###  Returns: dataframe with time_steps, day_min, time_after_start, and 
###    time_proportion data
###  Notes: need to run AddSolarTimes() prior to running this function
###  Blake Massey
###  2014.06.11

AddTimeStepProportion <- function(df = df, 
                                  by = "id",
                                  time_step = "15 min",
                                  tz = "Etc/GMT+5") {
  suppressPackageStartupMessages(require(plyr))
  df <- df
  df$by <- df[,by]
  DailyTimeStepCount<-function (df=df){  
    days <- subset(df, select=c("id", "date", "hr_before_sunrise", 
      "hr_after_sunset"))
    days <- ddply(days, .(date), function(x) x[1, ])  # first records
    days$time_steps <- NA
    days$day_min <- NA
    for (i in 1:nrow(days)) {
      day <- days[i,]
      days[i,"time_steps"] <- length(seq(day[,"hr_before_sunrise"], 
        day[,"hr_after_sunset"], time_step))
      days[i,"day_min"] <- as.integer(difftime(day[,"hr_after_sunset"], 
        day[,"hr_before_sunrise"],tz=tz, units = ("mins")))
    }
  return(days)
  }
  df2 <- ddply(df, .(by), DailyTimeStepCount) 
  df2 <- merge(df,df2, all.x=TRUE)
  df2$time_after_start <- as.integer(difftime(df2$datetime, 
    df2$hr_before_sunrise, tz=tz, units = ("mins"))) 
  df2$time_proportion <- df2$time_after_start/df2$day_min
  df2$by <- NULL
  return(df2)
}

# ConvertNestIdToNum Function --------------------------------------------------

###  Converts alphanumeric "nest_id" column to a numeric "nest_id_num" column 
###  Usage:  ConvertNestIdToNum(df)
###  Arguments: df = input dataframe with "nest_id" column  
###  Returns: A dataframe with a numberic "nest_id_num" column       
###  Notes: If a "nest_id_num" column doesn't exist, the function automatically 
###    creates one  
###  Blake Massey
###  2014.12.08 
   
ConvertNestIdToNum <- function(df){ 
  suppressPackageStartupMessages(require(stringr))
  df <- df
  if(!"nest_id_num" %in% colnames(df)) {
    df$nest_id_num <- NA
  }
  nest_id <- df$nest_id  
  territory_number <- sapply(strsplit(nest_id, "[A-Z]"), "[", 1)  
  nest_letter <- str_extract(nest_id, "[A-Z]")
  nest_number <- sapply(strsplit(nest_id, "[A-Z]"), "[", 2)
  nest_number[is.na(nest_number)] <- ""
  letters <- LETTERS[1:26]  # had only gone to "k" in 2012
  numbers <- formatC(1:26, width = 2, format = "d", flag = "0") 
  for (i in 1:length(letters)){
    nest_letter <- gsub(letters[i], numbers[i], nest_letter)
  }
  df$nest_id_num <- as.numeric(sprintf("%s%s%s", territory_number, nest_letter,
    nest_number))
  return(df)
}

# ConvertNestNumToId Function --------------------------------------------------

###  Converts a numeric "nest_id_num" column to an alphanumeric "nest_id" column 
###  Usage:  ConvertNestIdToNum(df)
###  Arguments: df = input dataframe with "nest_id_num" column  
###  Returns: A dataframe with a numberic "nest_id" column       
###  Notes: If a "nest_id" column doesn't exist, the function automatically 
###    creates one  
###  Blake Massey
###  2014.12.08 

ConvertNestNumToId <- function(df){ 
  df <- df
  if(!"nest_id" %in% colnames(df)) {
    df$nest_id <- NA
  }
  nest_id_num <- df$nest_id_num
  letters <- LETTERS[1:26]  # had only gone to "k" in 2012
  numbers <- formatC(1:26, width = 2, format = "d", flag = "0") 
  for (i in 1:length(nest_id_num)){
    if (nchar(nest_id_num[i]) == 3){
      territory_number<- sprintf("%03d", as.numeric(substr(nest_id_num[i],1,1)))
      nest_letter <- substr(nest_id_num[i], 2, 3)
      for (j in 1:length(letters)){
        nest_letter <- gsub(numbers[j], letters[j], nest_letter)
      }
      df$nest_id[i] <- paste(territory_number, nest_letter, sep="") 
    }
    if (nchar(nest_id_num[i]) == 4){
      territory_number<- sprintf("%03d", as.numeric(substr(nest_id_num[i],1,2)))
      nest_letter <- substr(nest_id_num[i], 3, 4)
      for (j in 1:length(letters)){
        nest_letter <- gsub(numbers[j], letters[j], nest_letter)
      }
      df$nest_id[i] <- paste(territory_number, nest_letter, sep="") 
    }
    if (nchar(nest_id_num[i]) == 5){
      territory_number <- substr(nest_id_num[i], 1, 3)
      nest_letter <- substr(nest_id_num[i], 4, 5)
      for (j in 1:length(letters)){
        nest_letter <- gsub(numbers[j], letters[j], nest_letter)
      }
      df$nest_id[i] <- paste(territory_number, nest_letter, sep="") 
    }
    if (nchar(nest_id_num[i]) == 7){
      territory_number <- substr(nest_id_num[i], 1, 3)
      nest_letter <- substr(nest_id_num[i], 4, 5)
      nest_number <- substr(nest_id_num[i], 6, 7)
      for (j in 1:length(letters)){
        nest_letter <- gsub(numbers[j], letters[j], nest_letter)
      }
      df$nest_id[i] <- paste(territory_number, nest_letter, nest_number, 
        sep="")
    }
    if (nchar(nest_id_num[i]) == 8){
      territory_number <- substr(nest_id_num[i], 1, 4)
      nest_letter <- substr(nest_id_num[i], 5, 6)
      nest_number <- substr(nest_id_num[i], 7, 8)
      for (j in 1:length(letters)){
        nest_letter <- gsub(numbers[j], letters[j], nest_letter)
      }
      df$nest_id[i] <- paste(territory_number, nest_letter, nest_number, 
        sep="") 
    }
  }
  return(df)
}

# CreateColorsByAny Function ---------------------------------------------------

###  Creates and/or displays dataframe of "by" variable and associated colors 
###  Usage: CreateColorsByAny(by, df, pal, output, plot, ...)
###  Arguments: by = variable to determine colors. Specific outcomes for: 
###               "behavior", "id", or "sex"
###             df = dataframe with "by" variable - only required if "by" is not
###               "behavior", "id", or "sex"
###             pal = color palette for CreateVarsColors(), default is NULL
###             b_pal = color palette from RColorBrewer for CreateVarsColors(),
###               default is "Set1"
###             output = logical, whether or not to return the dataframe, 
###               default is TRUE 
###             plot = logical, whether or not to display names and colors,
###               default is FALSE
###  Returns: df with "by" variable and hexidecimal colors 
###  Notes: Used in several other functions. Color palettes determined 
###    automatically for "behavior", "id", and "sex". Others set by pal or 
###    b_pal.
###  Blake Massey
###  2014.11.04

CreateColorsByAny <- function (by,
                              df,
                              pal = NULL,
                              r_pal = NULL,
                              b_pal = "Set1", 
                              output = TRUE,
                              plot = FALSE,
                              ...) {
  source('C:/Work/R/Functions/kml.R')
  if (!is.null(by)){
  if (by == "behavior" || by == "id" || by == "sex") {
    if (by == "behavior") by_colors <- CreateColorsByMetadata(file=
      "C:/Work/R/Data/BAEA/BAEA_metadata.csv", metadata_id="behavior")
    if (by == "id") by_colors <- CreateColorsByMetadata(file=
      "C:/Work/R/Data/BAEA/BAEA_gps_data.csv", metadata_id="deploy_location")
    if (by == "sex") by_colors <- CreateColorsByMetadata(file=
      "C:/Work/R/Data/BAEA/BAEA_metadata.csv", metadata_id="sex")
  } else {
    by_colors <- CreateColorsByVar(df=df, by=by, pal=pal, r_pal = r_pal, b_pal = 
      b_pal)
  }
  } else {
    by_colors <- CreateColorsByVar(df=df, by=by, pal=pal, r_pal = r_pal, b_pal = 
      b_pal)
  }
  if (plot == TRUE) PlotColorPie(by_colors)
  if (output == TRUE) return(by_colors)  
}

# ExportKMLTelemetryBAEA Function ----------------------------------------------

###  Create a Google Earth KML file (points and multitrack) from lat/long 
###    coordinates, defaults are specifically set to my BAEA data
###  Usage: ExportKMLTelemetryBAEA(df, id, datetime, lat, long, speed, alt,
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
###  2014.12.10

ExportKMLTelemetryBAEA <- function (df,
                                    id = "id",                                  
                                    datetime = "datetime",
                                    lat = "lat",
                                    long = "long", 
                                    alt = "alt",
                                    alt_mode = "clampToGround",
                                    speed = "speed",
                                    agl = NULL,
                                    behavior = NULL,
                                    point_color = "deploy_location",
                                    point_metadata = file.path("C:/Work/R",
                                      "Data/BAEA/BAEA_gps_data.csv"),
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
                                    file = "BAEA Data.kml",
                                    kml_folder = "C:/Users/Blake/Desktop") {
  source('C:/Work/R/Functions/kml.R')
  if (point_color == "behavior" || point_color == "sex") {
      point_metadata = file.path("C:/Work/R/Data/BAEA/BAEA_metadata.csv")
  }
  if (!is.null(path_color) &&  path_color == "sex") {
      path_metadata = file.path("C:/Work/R/Data/BAEA/BAEA_metadata.csv")
  }
  ExportKMLTelemetry(df=df, id=id, datetime=datetime, lat=lat, long=long, 
    alt=alt, alt_mode=alt_mode, speed=speed, agl=agl, behavior=behavior, 
    point_color=point_color, point_metadata=point_metadata, point_pal=point_pal, 
    point_r_pal=point_r_pal, point_b_pal=point_b_pal, extrude=extrude, 
    path=path, path_color=path_color, path_metadata=path_metadata, 
    path_pal=path_pal, path_r_pal= path_r_pal, path_b_pal=path_b_pal, 
    arrow=arrow, icon_by_sex=icon_by_sex, labelscale=labelscale, 
    dateformat=dateformat, timeformat=timeformat, datetimeformat=datetimeformat,
    file=file, kml_folder=kml_folder)
}

# ExtractFlightData Function ---------------------------------------------------

###  Extracts the speed of flights for each bird  
###  Usage: ExtractFlightData(df)
###  Arguments: df = dataframe with locations
###             min_speed = minimum speed to be classified as flight
###  Returns: a dataframe of locations that are           
###  Notes: 
###  Blake Massey
###  2014.10.20

ExtractFlightSpeed <- function(df,
                               min_speed = 2) {
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(grid))
  suppressPackageStartupMessages(require(plyr))
  source("C://Work/R/Functions/gen.R")
  source("C://Work/R/Functions/baea.R")
  df <- df
  # get data that only fits the "flight" criteria
  flight <- df[which(df$speed>min_speed),]
  return(flight)
}

# ExtractMovements Function ----------------------------------------------------

###  Extracts detected movements for each bird  
###  Usage: ExtractMovements(df, min_step_length, max_step_time)
###  Arguments: df = dataframe with locations
###             min_step_length = threshold distance for movement classification
###             max_step_time = threshold difference in time for movement 
###               classification 
###             by = column name of unique identifier to split data, default is
###               "id"
###  Returns: a dataframe of selected movements
###  Notes: 
###  Blake Massey
###  2014.10.21

ExtractMovements <- function(df,
                             by = "id",
                             min_step_length = 50,
                             max_step_time = 20){
  df <- df 
  movements <- ddply(df, by, function(df){
    movements <- df[which(df$step_length >= min_step_length & 
    df$step_time <= max_step_time), ]}) 
  return(movements)
}

# ExtractRoostData Function ----------------------------------------------------

###  Returns a dataframe of sex, arr_diff_min, and dep_diff_min
###  Usage: ExtractRoostData(df)
###  Arguments: df = dataframe with "sex", "datetime", "behavior", 
###                    "hr_after_sunset", and "hr_before_sunrise
###  Returns: dataframe
###  Notes: used in other functions
###  Blake Massey
###  2014.06.02

ExtractRoostData <- function(df = df){
  arrive <- subset(df, behavior=="arrive")
  arrive$arr_diff_min <- as.integer(difftime(arrive$hr_after_sunset, 
    arrive$datetime))
  arrive <- subset(arrive, select=c("sex", "arr_diff_min"))
  arrive$dep_diff_min <- NA
  arrive$roost <- "arrive"
  depart <- subset(df, behavior=="depart")
  depart$dep_diff_min <- as.integer(difftime(depart$datetime, 
    depart$hr_before_sunrise))
  depart <- subset(depart, select=c("sex", "dep_diff_min"))
  depart$arr_diff_min <- NA
  depart$roost <- "depart"
  output <- rbind(arrive, depart)
  output$diff_min <- NA
  output$diff_min <- ifelse(!is.na(output$arr_diff_min), output$arr_diff_min, 
    output$dep_diff_min)
  return(output)
}

# FilterLocations Function -----------------------------------------------------

###  Filters location data by individual and dates  
###  Usage: FilterLocations(df, id, individual, start, end)
###  Arguments: df = dataframe with locations
###             id = column name of unique identifier
###             individual = individual/s (from id column) to keep,
###               format should be c(id, id), default is all
###             start = filter start date, default is 1970-01-01   
###             end = filter end date, default is current date
###             behavior = filter to specific behaviors
###  Returns: dataframe with subsetted rows         
###  Notes: defaults are specific to my file directories and locations
###  Blake Massey
###  2014.05.05

FilterLocations <- function(df = df, 
                            id = "id", 
                            individual = NULL, 
                            start = NULL, 
                            end = NULL, 
                            behavior = NULL){
  suppressPackageStartupMessages(require(rgdal))
  suppressPackageStartupMessages(require(sp))
  if (is.null(start) || start == ""){
    start <- "2013-01-01"
  }
  if (is.null(end) || end == ""){
    today <- Sys.Date()
    end <- format(today, format="%Y-%m-%d")
  }
  starts <- paste("Start date: ", start, sep="")
  ends <- paste("End date: ", end, sep="")
  writeLines(noquote(starts))
  writeLines(noquote(ends))
  end <- as.POSIXct(end)
  end <- trunc(end, "days") + 60*60*24
  df <- df[df$datetime >= as.POSIXct(start) & df$datetime <= end,]
  row.names(df) <- NULL
  if(all(is.null(individual)) || individual == ""){
    writeLines(noquote ("All individuals included")) 
  } else {
    writeLines(noquote(paste("Filtered to individual(s):", individual, sep="")))
    df <- df[df[,id] %in% individual,]
    row.names(df) <- NULL
  }
  if(is.null(behavior)){
    writeLines(noquote ("All behaviors included")) 
  } else {
    writeLines(noquote(paste("Filtered to behavior(s):", behavior, sep="")))
    df <- df[df[,"behavior"] %in% behavior,]
    row.names(df) <- NULL
  }
  return(df)
}

# IfElseTimedateCompare Function -----------------------------------------------

### Examine two time columns, returns the greater value
### Usage: IfElseTimedateCompare(df, col1, col2, result, default_tz, tz)
### Arguments: df = dataframe
###            col1 = col1 is compared to col2 
###            col2 = col2 is compared to col1 
###            result = result of selection between col1 and col2
###            default_tz = timezone that timedate reverts to, based on OS time
###            tz = timezone of original data (may be different from default_tz)
### Output: original df with a new "result" column 
### Notes: ensure that result is in the correct timezone. Used in
###   RoostArrivalDeparture().
### Blake Massey
### 2014.05.27

IfElseTimedateCompare <- function(df = df,
                                  col1 = "col1", 
                                  sign = ">", 
                                  col2 = "col2", 
                                  result = "result", 
                                  default_tz = default_tz, tz=tz) {  
  safe.ifelse <- function(cond, yes, no) {
    structure(ifelse(cond, yes, no), class = class(yes))
  }
  df$col1<-df[,col1]
  df$col2<-df[,col2]
  df$result <- NA
  if (sign == ">"){
    df$result <- safe.ifelse(df$col1 >= df$col2, df$col1, df$col2)
  }
  if (sign == "<"){
    df$result <- safe.ifelse(df$col1 <= df$col2, df$col1, df$col2)
  }
  df$result <- force_tz(df$result, tzone = default_tz) 
  df$result <- with_tz(df$result, tzone=tz)  
  df$col1 <- NULL
  df$col2 <- NULL
  df_length <- length(df)
  colnames(df)[df_length] <- result
  return(df)
}

# IfElseTimedateNA Function ----------------------------------------------------

### Examine two time columns, if one is NA, the other time is returned.
### Usage: IfElseTimedateNA(df, col1, col2, result, default_tz, tz)
### Arguments: df = dataframe
###            col1 = column to check if NA. If not NA, then col1 is returned
###            col2 = if col1 is NA, then col2 is returned
###            result = result of selection between col1 and col2
###            default_tz = timezone that timedate reverts to, based on OS time
###            tz = timezone of original data (may be different from default_tz)
### Output: original df with a new "result" column
### Notes: Ensure that result is in the correct timezone. Used in 
###   RoostArrivalDeparture().
### Blake Massey
### 2014.05.27

IfElseTimedateNA <- function(df = df, 
                             col1 = "col1", 
                             col2 = "col2", 
                             result = "result", 
                             default_tz = default_tz, 
                             tz = tz) {  
  safe.ifelse <- function(cond, yes, no) {
    structure(ifelse(cond, yes, no), class = class(yes))
  }
  df$col1<-df[,col1]
  df$col2<-df[,col2]
  df$result <- NA
  df$result <- safe.ifelse(is.na(df$col1), df$col2, df$col1)
  df$result <- force_tz(df$result, tzone = default_tz)
  df$result <- with_tz(df$result, tzone=tz)  
  df$col1 <- NULL
  df$col2 <- NULL
  df_length <- length(df)
  colnames(df)[df_length] <- result
  return(df)
}

# ImportBAEA Function ----------------------------------------------------------

###  Imports baea.csv and merges with existing
###  Usage: ImportBAEA(units, import, tz)
###  Arguments: existing = exisiting file to merge with baea. Can be NULL.
###             import = select to import baea. Default is TRUE.
###             tz = timezone. Default is "Etc/GMT+5".
###  Returns: merged baea file          
###  Notes: directory defaults are specific to my computer
###  Blake Massey
###  2014.05.20

ImportBAEA <- function(existing = deployed, 
                       import = TRUE, 
                       tz = "Etc/GMT+5") {
  suppressPackageStartupMessages(require(lubridate))
  if (import == TRUE) {
  baea <- read.csv(file="C:/Work/R/Data/BAEA/BAEA.csv", header = TRUE, 
    stringsAsFactors=FALSE)
  date_cols <- c("date","on_hand","deployed", "end_data",  "failed",  "removed",
    "recovered")
  for (i in date_cols) {  
    baea[,i] <- as.Date(baea[,i], "%Y-%m-%d", tz=tz)  # otherwise they are chars
  }
  datetime_cols <- c("datetime","sunset", "sunrise",  "solarnoon",  
    "hr_before_sunrise", "hr_after_sunset")
  for (i in datetime_cols) {  
    baea[,i] <- as.POSIXct(baea[,i], tz=tz, usetz=FALSE)
  }
  if (!is.null(existing)) {  
    existing <- subset(existing, date > (as.Date(max(baea$date)) - days(3)))
    baea <- subset(baea, date <= (as.Date(max(baea$date)) - days(3))) 
      max(baea$date)
      min(existing$date)
  # 3 days are removed from baea to ensure that AddSegmentTimeLength, etc. was 
  # done on a full dataset. The baea and existing datasets should not overlap.
  if(!("sunrise" %in% colnames(existing))) {
    existing <- AddSolarTimes(existing)
  }
  existing <- AddStepLengthAndAngles(existing)
  existing <- AddStepTime(existing)
  existing <- AddFirstLastDistance(existing)
  baea_full <- rbind(baea, existing)
  baea_full <- unique(baea_full)
  baea_full <- baea_full[with(baea_full,order(id,datetime)),]
  row.names(baea_full) <- NULL
  date <- Sys.Date()
  outfile <- paste("C:/Work/R/Data/BAEA/Archive/BAEA_", date, ".csv", sep ="")
  if (!file.exists(outfile)) {
  writeLines(noquote(paste("Merging existing and import")))
  write.csv(baea_full, file=outfile, row.names=FALSE)
  writeLines(noquote(c("Writing: ", outfile, sep = "")))
  }
  write.csv(baea_full, file="C:/Work/R/Data/BAEA/BAEA.csv", row.names=FALSE)  
    # rewrites import file
  writeLines(noquote("Writing: \"C:/Work/R/Data/BAEA/BAEA.csv\""))
  baea <- baea_full
  } 
  }
  if (import == FALSE) {
  writeLines(noquote(paste("BAEA.csv was NOT imported")))
  if (!is.null(existing)) {
  if(!("sunrise" %in% colnames(existing))) {
    existing <- AddSolarTimes(existing)
  }
  existing <- AddStepLengthAndAngles(existing, by = "id")
  existing < AddStepTime(existing, by = "id")
  existing <- AddFirstLastDistance(existing)       
  writeLines(noquote(paste("Coverted existing to baea", sep="")))
  baea <- existing
  } else {
  writeLines(noquote("Nothing imported or converted"))
  }
  }
  return(baea)
}

# PlotBehaviorProportionBar Function -------------------------------------------

###  Plots daily behavior distributions
###  Usage: PlotBehaviorProportionBar(df, breaks)
###  Arguments: df = dataframe with "sex", "behavior", and "time_proportion"
###             breaks = number of breaks in daily period
###  Returns: facetted plot of behavior proportion over daily period          
###  Notes: 
###  Blake Massey
###  2014.10.04

PlotBehaviorProportionBar<-function(df = df, 
                                    breaks = 20){
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(grid))
  suppressPackageStartupMessages(require(scales))
  suppressPackageStartupMessages(require(reshape))
  source('C:/Work/R/Functions/gps.R')
  behavior_colors <- CreateColorsByMetadata(file=
      "C:/Work/R/Data/BAEA/BAEA_metadata.csv", metadata_id="behavior")
  df$behavior <- factor(df$behavior)
  CutProportion <- function(data,breaks=breaks) {  
    b <- seq(0, 1, length=2*breaks+1)
    brk <- b[0:breaks*2+1]
    k <- cut(data, breaks=brk)
  }
  CutProportionMid <- function(data,breaks=breaks) {  
    b <- seq(0, 1, length=2*breaks+1)
    brk <- b[0:breaks*2+1]
    mid <- b[1:breaks*2]
    k <- cut(data, breaks=brk)
    mid[k]
  }
  df$bins <- CutProportion(df$time_proportion, breaks)
  df$bins_mid <- factor(CutProportionMid(df$time_proportion,breaks))
  melted <- melt(ddply(df, .(sex, bins_mid), 
    function(x){prop.table(table(x$behavior))}))
  names(melted)[names(melted) == 'variable'] <- 'behavior'
  melted$bins_mid <- as.numeric(as.character(melted$bins_mid))
  ggplot(melted, aes(x = bins_mid, y=value, ymax=1, fill= behavior)) +  
    facet_grid(~ sex) + geom_bar(stat="identity") +
    scale_fill_manual(values = behavior_colors) +
    scale_x_continuous(breaks=seq(0,1,.1)) + 
    theme(panel.margin = unit(1, "lines")) +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) +
    labs(x="Daily Period", 
    y="Behavior Proportion", 
    title="Daily Behavior Distributions")
}

# PlotBehaviorProportionLine Function ------------------------------------------

###  Plots daily behavior distributions
###  Usage: PlotBehaviorProportionLine(df, breaks)
###  Arguments: df = dataframe with "sex", "behavior", and "time_proportion"
###             breaks = number of breaks in daily period
###  Returns: facetted plot of behavior proportion over daily period          
###  Notes: 
###  Blake Massey
###  2014.10.03

PlotBehaviorProportionLine<-function(df = df, 
                                     breaks = 20){
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(grid))
  suppressPackageStartupMessages(require(scales))
  suppressPackageStartupMessages(require(reshape))
  df$behavior<-factor(df$behavior)
  source('C:/Work/R/Functions/gps.R')
  behavior_colors <- CreateColorsByMetadata(file=
      "C:/Work/R/Data/BAEA/BAEA_metadata.csv", metadata_id="behavior")
  CutProportion <- function(data,breaks=breaks) {  
    b <- seq(0, 1, length=2*breaks+1)
    brk <- b[0:breaks*2+1]
    k <- cut(data, breaks=brk)
  }
  CutProportionMid <- function(data,breaks=breaks) {  
    b <- seq(0, 1, length=2*breaks+1)
    brk <- b[0:breaks*2+1]
    mid <- b[1:breaks*2]
    k <- cut(data, breaks=brk)
    mid[k]
  }
  df$bins <- CutProportion(df$time_proportion, breaks)
  df$bins_mid <- factor(CutProportionMid(df$time_proportion,breaks))
  melted <- melt(ddply(df, .(sex, bins_mid), 
    function(x){prop.table(table(x$behavior))}))
  names(melted)[names(melted) == 'variable'] <- 'behavior'
  melted$bins_mid <- as.numeric(as.character(melted$bins_mid))
  ggplot(melted, aes(x = bins_mid, y=value, ymax=1, group= behavior, color=
    behavior)) +  facet_grid(~ sex) + theme(panel.margin=unit(1, "lines")) +
    scale_color_manual(values=behavior_colors)+
    geom_line(stat="identity", size=1.5) + 
    scale_x_continuous(breaks=seq(0,1,.1)) +
    theme(plot.title=element_text(size=22)) + 
    theme(text=element_text(size=20, colour="black")) + 
    theme(axis.text=element_text(colour="black")) + labs(x="Daily Period", 
    y="Behavior Proportion", title="Daily Behavior Distributions") 
}

# PlotStepLengths Function -----------------------------------------------------

###  Plots a histogram of the step-lengths for each bird  
###  Usage: PlotStepLengths(df)
###  Arguments: df = dataframe with flights
###             id = column name of unique identifier, default = "id"
###             xlim = x value limit, default is NULL
###             bin_width = bin size, default is: x-value range/30 
###  Returns: a plot of step-length distances          
###  Notes: 
###  Blake Massey
###  2014.10.21

PlotStepLengths <- function(df,
                            id = "id",
                            xlim = NULL,
                            bin_width = NULL){
  source('C:/Work/R/Functions/gen.R')
  source('C:/Work/R/Functions/baea.R')
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(grid))
  suppressPackageStartupMessages(require(plyr))
  df<-df
  sum_move <- SummarizeSE(df, "step_length", "id", na_rm=TRUE)
  id_colors <- CreateColorsByID(output=TRUE)  
  grid <- seq(min(df$step_length, na.rm=TRUE), max(df$step_length, na.rm=TRUE),
    length = 100)
  probs = c(0.5, 0.75, 0.95)
  sumstats <- ddply(df, id, function(df){
    quantile(df$step_length, probs = probs, na.rm=TRUE)
  })
  q_probs<-paste("q",probs, sep="")
  colnames(sumstats)[2:(length(probs)+1)] <- q_probs 
  if (is.null(xlim)) xlim <- max(df$step_length, na.rm=TRUE)
  g <- ggplot(df, aes(x=step_length, fill=id)) + facet_wrap( ~ id)  +
    scale_fill_manual(values=id_colors) +
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) + 
    xlab("Step length") + ylab("Density") + scale_x_continuous(limits=c(0,xlim))
  if (is.null(bin_width)) bin_width = xlim/30
  g <- g + geom_bar(aes(y = ..density.., fill=id), colour="black", 
    binwidth=bin_width)
  build<-ggplot_build(g)
  sumstats$xmax <- max(build$panel$ranges[[1]]$x.range)
  sumstats$ymax <- max(build$panel$ranges[[1]]$y.range)
  g + geom_vline(data=sumstats, aes(xintercept=q0.5), linetype="longdash", 
    size=1, colour="black")  +  
  geom_vline(data=sumstats, aes(xintercept=q0.75), linetype="dashed", 
    size=1, colour="grey20") +
  geom_vline(data=sumstats, aes(xintercept=q0.95), linetype="dashed", 
    size=1, colour="grey30") +
  geom_text(data=sumstats, aes(x=q0.5 + (xmax*0.02), y=ymax*.75, 
    label=paste("Median:", "\n", signif(q0.5,3), sep="")),
    color="black", hjust=0) +  
  geom_text(data=sumstats, aes(x=q0.75+(xmax*0.02), y=ymax*.5, 
    label=paste("75%:", "\n", signif(q0.75,3), sep="")),
    color = "grey20", hjust=0) +  
  geom_text(data=sumstats, aes(x=q0.95+(xmax*0.02), y=ymax*.25, 
    label=paste("95%:", "\n", signif(q0.95,3), sep="")),
    color = "grey30", hjust=0)
}

# PlotRoostECDF Function -------------------------------------------------------

###  Plots Roost empirical distribution funtion and fitted Weibull cumulative 
###    distribution function 
###  Usage: PlotRoostECDF(df, pars)
###  Arguments: df = dataframe of location data
###             pars = simulation parameters with male/female, 
###               arrive/depart, shape/scale 
###  Returns: facetted plot with roost data and fitted weibull distributions          
###  Notes: empirical distribution function extends to 15 min past the end of 
###    last time.
###  Blake Massey
###  2014.06.03

PlotRoostECDF <- function(df = baea, 
                          pars = baea_pars) {
  suppressPackageStartupMessages(require(ggplot2))
  df <- ExtractRoostData(df)
  df2 <- ExtractRoostPars(pars)
  df[df=="m"]<-"male"
  df2[df2=="m"]<-"male"  
  df[df=="f"]<-"female"
  df2[df2=="f"]<-"female"
  df2$diff_min<-.85*(max(df$diff_min)+15)
  df2$lab<-paste("shape: ",round(df2$shape,2))
  df2$lab2<-paste("scale: ", round(df2$scale,2))
  vec <- with(df, seq(0, max(diff_min)+15, length=max(diff_min)+16))
  weibull <- ddply(df2, .(sex,roost), function(df) {
    data.frame(diff_min=vec,
    density = pweibull(vec, shape=df$shape, scale=df$scale))
  })
  df2$density<-.25*(max(weibull$density))
  df2$density2<-.15*(max(weibull$density))
  ggplot(df,aes(x=diff_min)) +  
  stat_ecdf(aes(colour=sex), size=1) +   
  geom_line(aes(x=diff_min, y=density), data=weibull, 
            colour="orange", size=1) +
  geom_text(aes(x=diff_min, y=density, label=lab),
    data=df2) +
  geom_text(aes(x=diff_min, y=density2, label=lab2),
    data=df2) +
  facet_grid(roost ~ sex) +
  theme(plot.title=element_text(size=22)) +
  theme(text=element_text(size=20, colour="black")) +
  theme(axis.text=element_text(colour="black")) +
  labs(x='Time Difference from Start/End of Sim Day', 
   y='Cumulative Probability Density', 
   title="Roost Data with Fitted Weibull Distributions", 
   colour= "Sex")
}

# PlotRoostHistogram Function --------------------------------------------------

###  Plots roost times and fitted Weibull probability distributions
###  Usage: PlotRoostHistogram(df, pars)
###  Arguments: df = dataframe of location data
###             pars = simulation parameters with male/female, 
###               arrive/depart, shape/scale 
###  Returns: facetted plot with roost data and fitted weibull distributions          
###  Notes: weibull probability function extends 15 min past last time.
###  Blake Massey
###  2014.06.03

PlotRoostHistogram <- function(df, 
                               pars, 
                               binwidth = 15) {
  df <- ExtractRoostData(df)
  df2 <- ExtractRoostPars(pars)
  df[df=="m"] <- "male"
  df2[df2=="m"] <- "male"  
  df[df=="f"] <- "female"
  df2[df2=="f"] <- "female"
  df2$diff_min <- .85*max(df$diff_min)+15
  df2$lab <- paste("shape: ", round(df2$shape,2))
  df2$lab2 <- paste("scale: ", round(df2$scale,2))
  grid <- with(df, seq(0, max(diff_min)+15, length=max(diff_min)+16))
  weibull <- ddply(df2, .(sex,roost), function(df) {
    data.frame(diff_min=grid,
    density = dweibull(grid, shape=df$shape, scale=df$scale))
  })
  df2$density <- max(weibull$density)
  df2$density2 <- .85*df2$density
  ggplot(df) + 
  geom_histogram(aes(x=diff_min, y=..density..), binwidth=binwidth) +  
  geom_line(aes(x=diff_min, y=density), data=weibull, 
    colour = "orange", size=1) +   
  geom_text(aes(x=diff_min, y=density, label=lab),
    data=df2) +
  geom_text(aes(x=diff_min, y=density2, label=lab2),
    data=df2) +
  facet_grid(roost ~ sex) +
  theme(plot.title=element_text(size=22)) +
  theme(text=element_text(size=20, colour="black")) +
  theme(axis.text=element_text(colour="black")) +
  labs(x='Time Difference from Start/End of Sim Day', 
   y='Probability Density', 
   title="Histogram of Roost Data with Fitted Weibull Distributions")
}
