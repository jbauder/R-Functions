# --- DEPRECATED FUNCTIONS -----------------------------------------------------
# Deprecated functions. They are commented out because I do not want them to 
# load
# ------------------------------------------------------------------------------

# CreateSimlist Function -------------------------------------------------

###  Creates a simlist for a simulation
###  Usage: CreateSimlist(time_step, days, runs, tz, nests)
###  Arguments: time_step = time step length, based on lubriate times, default 
###               15 min 
###             days = number of days to run simulation, default 10
###             runs = number of runs of simulation (currently not used)
###             tz = timezone for dataset, default is "Etc/GMT+5"
###             nests = dataframe of nest locations used to create sim_list, 
###               default is "C:/Work/R/Data/Simulation/Sim_nests.csv"
###  Returns: a list object with time_steps, runs, and individual nest data
###  Notes: day 1 is clutch initiation date
###  Blake Massey
###  2014.10.11

CreateSimlist <- function(time_step = "15 min", 
                          days = 10, 
                          runs = 1, 
                          tz = "Etc/GMT+5",
                          nests = "C:/Work/R/Data/Simulation/Sim_nests.csv") { 
  suppressPackageStartupMessages(require(maptools)) 
  suppressPackageStartupMessages(require(lubridate))
  Make_dates <- function(df=df, 
                         days=days) {
    long <- df$nest_long  # sets longitude of nest
    lat <- df$nest_lat  # sets latitude of nest
    day_seq <- seq(from = as.POSIXct(df$clutch_initiation, "%m/%d/%Y", tz =tz),
                length.out = days, by = "days")  # creates sequence of dates
    nest_age <- seq(1, days, 1)  # day 1 is clutch inititation date
    coord <- matrix(c(x,y), nrow = 1)  # coordinates for calculating solartimes
    sunrise <- sunriset(coord, day_seq, direction="sunrise", POSIXct.out = TRUE)
    sunrise$sunrise <- as.POSIXct(sunrise$time, tz=tz) 
    sunset <- sunriset(coord, day_seq, direction="sunset",POSIXct.out = TRUE)
    sunset$sunset <- as.POSIXct(sunset$time, tz=tz) 
    solarnoon <- solarnoon(coord, day_seq, POSIXct.out=TRUE)
    solarnoon$solarnoon <- as.POSIXct(solarnoon$time, tz=tz) 
    dates <- data.frame(as.Date(day_seq),nest_age, sunrise$sunrise, 
      solarnoon$solarnoon, sunset$sunset)    
    colnames(dates) <- c("date", "nest_age", "sunrise", "solarnoon", "sunset")
    dates$start_time <- dates$sunrise - 3600  # to subtract an hour from sunrise
    dates$end_time <- dates$sunset + 3600  # to add an hour to sunset
    for (i in 1:nrow(dates)) {
      data <- dates[i,]  # selects each day
      steps <- seq(data[,"start_time"], data[,"end_time"], time_step) 
      dates[i,"steps"] <- length(steps)
      dates[i,"simday_min"] <- as.integer(difftime(data[,"end_time"],
        data[,"start_time"], tz=tz, units=("mins")))  
    }
    return(dates)  
  }
  Make_locations <- function(df=dates) {
    datetime <- numeric(0)  # need a blank vector
    days <- nrow(df)  # calculates number of rows in solartimes dataframe
    simtimes <- as.data.frame(NULL)
    for (i in 1:days) {
      data <- df[i,]  # selects each day
      time <- seq(data[,"start_time"], data[,"end_time"], time_step)
      datetime <- with_tz(time, tzone=tz)  # keeps the time as EST, not EDT
      date <- as.Date(datetime, tz=tz)
      time_after_start <- as.integer(difftime(time, data[,"start_time"],
        tz=tz, units = ("mins"))) 
      time_before_end <- as.integer(difftime(data[,"end_time"], time, tz=tz,
        units = ("mins")))
      simday_min <- as.integer(difftime(data[,"end_time"], data[,"start_time"],
        tz=tz, units = ("mins"))) 
      time_proportion <- time_after_start/simday_min
      timestep_after_start <- seq(0, length(time_after_start)-1, 1)
      timestep_before_end <- rev(timestep_after_start)
      timestep_proportion <- timestep_after_start/length(time_after_start)
      daily_simtimes <- data.frame(date, datetime, time_proportion,
        time_after_start, time_before_end, timestep_after_start,
        timestep_before_end, timestep_proportion)
      simtimes <- rbind(simtimes, daily_simtimes)
    }
    simtimes$datetime <- with_tz(simtimes$datetime, tzone = tz)
    simlist <- simtimes
    namevector <- c("behavior","forage_suc", "lat", "long", "alt")  # adds cols
    simlist[,namevector] <- NA  # adds columns and makes all values NA
    simlist <- simlist[with(simlist, order(datetime)),]
    row.names(simlist) <- NULL
    return (simlist)  
  }
  Make_individual <- function(df=df, 
                              locations=locations){
    id <- df[["nest_id"]]
    name <- df[["eagle_id"]]  
    id_f <- paste(id, "-1", sep="")
    name_f <- paste(name,"_f", sep="")
    sex_f <- "f"
    female <- list(id_f, name_f, sex_f, locations)    
    id_m <- paste(id, "-2", sep="")
    name_m <- paste(name,"_m", sep="")
    sex_m <- "m"
    male <- list(id_m, name_m,sex_m, locations)
    individuals <- list(female, male)
    return(individuals)
  }  
  sim_nests <- nests
  simlist <- list()  # creates a blank list object
  simlist$time_step <- time_step  # sets time step length
  simlist$runs <- runs  # sets number of simulation runs
  simlist$nest <- NULL #makes a blank dataframe 
  for (i in 1:nrow(sim_nests)){ 
    df <- sim_nests[i, ] #runs through each row of simnests
    core <- c(nest_id=df$nest_id, nest_site=df$nest_site, 
      long=as.numeric(df$nest_long), lat=as.numeric(df$nest_lat), 
      alt=as.numeric(df$nest_alt), 
      clutch_size=df$clutch_size, clutch_initiation=df$clutch_initiation)
    dates <- Make_dates(df=df,days=days)  # runs Make_date function
    locations <- Make_locations(df=dates)
    individual <- Make_individual(df=df, locations=locations)
    simlist$nest[[i]] <- list(core, dates, individual)
    names(simlist$nest[[i]]) <- c("attributes","dates","individuals")
    names(simlist[[3]][[i]][[3]]) <- c("female", "male")
    names(simlist[[3]][[i]][[3]][[1]]) <- c("bird_id","bird_name", "sex",
      "locations")
    names(simlist[[3]][[i]][[3]][[2]]) <- c("bird_id","bird_name", "sex",
      "locations")
  }
  names<-sim_nests$name
  names(simlist[[3]])<-names
  return (simlist)
}
