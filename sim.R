# --- SIMULATION FUNCTIONS -----------------------------------------------------
# General functions for simulation construction and analysis of simulation
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
    x <- df$nest_long  # sets longitude of nest
    y <- df$nest_lat  # sets latitude of nest
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

# ExtractSimlistDF Function ----------------------------------------------------

###  Compiles data from a simlist into a single dataframe
###  Usage: ExtractSimlistDF(ls)
###  Arguments: ls = list of nests in simualtion, default is "simlist_nests"
###  Returns: a dataframe of all the nest data
###  Notes: can be done before or after a simulation run
###  Blake Massey
###  2014.10.11

ExtractSimlistDF <- function(ls = simlist_nests) {
  df <- data.frame()
  for (i in 1:length(ls)) {
    id <- names(ls[i])
    df1 <- cbind(id, ls[[i]])    
    df <- rbind(df, df1)
  }
  return(df)
}

# PredictArrival Function ------------------------------------------------------

###  Predict roost arrival based on time until end of simulation day
###  Usage: PredictArrival(time, pars)
###  Arguments: time = time since being of simulation 
###             pars = vector containing Weibull shape and scale
###  Returns: 1 (arrive) or 0 (not arrive)   
###  Notes: need to run WeibullNLL prediction to get parameters
###  Blake Massey
###  2014.10.11

PredictArrival <- function(time_m = time_m, 
                          pars = pars) {
  p <- 1 - pweibull(q=time_m, shape=pars["shape"], scale=pars["scale"])
  z <- rbinom(1, size=1, prob=p)   
  return(z)
}

# PredictDayBehavior Function ------------------------------------------------------

###  Select behavior between roost departure and roost arrival
###  Usage: PredictDayBehavior(time_p, time_m, sim_pars)
###  Arguments: time_p = proportion of total time in simulated day 
###             time_m = time until end of simulation 
###             sim_pars = simulation parameters (not used yet)
###  Returns: a behavior ("forage", "loaf", "cruise")
###  Notes: very simplified version - no sim_pars used
###  Blake Massey
###  2014.10.11

PredictDayBehavior<-function(time_p = time_p, 
                             time_m = time_m, 
                             pars = pars){
  cruise <- dbeta(time_p, pars[["cruise_pars"]][["shape1"]], 
                          pars[["cruise_pars"]][["shape2"]])
  forage <- dbeta(time_p, pars[["forage_pars"]][["shape1"]], 
                          pars[["forage_pars"]][["shape2"]])
  loaf <- dbeta(time_p, pars[["loaf_pars"]][["shape1"]], 
                        pars[["loaf_pars"]][["shape2"]])
  nest <- dbeta(time_p, pars[["nest_pars"]][["shape1"]], 
                        pars[["nest_pars"]][["shape2"]])
  roost <- dbeta(time_p, pars[["roost_pars"]][["shape1"]], 
                         pars[["roost_pars"]][["shape2"]])
  territorial <- dbeta(time_p, pars[["territorial_pars"]][["shape1"]], 
                               pars[["territorial_pars"]][["shape2"]]
  )
  vec <- setNames(c(cruise, forage, loaf, nest, roost, territorial), 
                  c("cruise", "forage", "loaf", "nest", "roost", "territorial"))
  b <- names(which.max(vec))
  return(b)
}

# PredictDeparture Function ----------------------------------------------------
 
###  Predicts roost departure based on time since beginning of simulation day
###  Usage: PredictDeparture(time, pars)
###  Arguments: time = time since being of simulation 
###             pars = vector containing Weibull shape and scale
###  Returns: 1 (depart) or 0 (not depart)
###  Notes: need to run WeibullNLL prediction to get parameters 
###  Blake Massey
###  2013.10.11

PredictDeparture <- function(time, 
                          pars) {
    p <- pweibull(q=time, shape=pars["shape"], scale=pars["scale"])
    z <- rbinom(1, size=1, prob=p)   
    return(z)
}

# PredictDepartBehavior Function -----------------------------------------------

###  Select behavior after eagle departs roost/nest
###  Usage: PredictDepartBehavior(time, sim_pars)
###  Arguments: time = proportion of total time in simulated day, default is 
###               "time_proportion" 
###             sim_pars = simulation parameters (not used yet)
###  Returns: a behavior ("forage", "loaf", "cruise")
###  Notes: very simplified version - no sim_pars used
###  Blake Massey
###  2014.10.11

PredictDepartBehavior <- function(time_p = time_proportion) {  
  cruise <- dbeta(time_p, pars[["cruise_pars"]][["shape1"]], 
                          pars[["cruise_pars"]][["shape2"]])
  forage <- dbeta(time_p, pars[["forage_pars"]][["shape1"]], 
                          pars[["forage_pars"]][["shape2"]])
  loaf <- dbeta(time_p, pars[["loaf_pars"]][["shape1"]], 
                        pars[["loaf_pars"]][["shape2"]])
  nest <- dbeta(time_p, pars[["nest_pars"]][["shape1"]], 
                        pars[["nest_pars"]][["shape2"]])
  roost <- dbeta(time_p, pars[["roost_pars"]][["shape1"]], 
                         pars[["roost_pars"]][["shape2"]])
  territorial <- dbeta(time_p, pars[["territorial_pars"]][["shape1"]], 
                               pars[["territorial_pars"]][["shape2"]])
  vec <- setNames(c(cruise, forage, loaf, nest, roost, territorial), 
                  c("cruise", "forage", "loaf", "nest", "roost", "territorial"))
  b <- names(which.max(vec))
  return(b)
  return(b)
}

# RunSimulation Function -------------------------------------------------------

###  Runs a simulation on input simlist
###  Usage: RunSimulation(simlist, pars)
###  Arguments: simlist = simlist, created with CreateSimlist() 
###             pars = simulation parameters, default is 'pars'
###  Returns: input simlist with simulation data  
###  Notes: still a work in progress
###  Blake Massey
###  2014.05.11

RunSimulation <- function(simlist = simlist, 
                          pars = pars) {
  suppressPackageStartupMessages(require(plotKML))
  sim_pars <- pars
  for (i in 1:length(simlist[[3]])) {  # pulls out every nest
    nest <- simlist[[3]][[i]]  # pulls out an individual nest
    female <- nest[[3]][['female']][['locations']]
    male <- nest[[3]][['male']][['locations']]
    dates <- nest[[2]]  # pulls out nests dates
    date <- unique(dates[, "date"]) # pulls out the unique dates for the nest
  for (d in date) { 
    f_date <- female[female$date == d, ]  # creates a df for each date 
    m_date <- male[male$date == d, ]  # creates a df for each date
    nest_age <- dates[dates$date == d, "nest_age"]  # nest age for each date
    time <- unique(f_date$datetime)  # pulls of the unique time for each date
    t<-length(time)  # number of times in the day
    f_Depart <- 0  # to determine when the bird departs for day 
    f_Arrive <- 0  # to determine when the bird arrives for day
    m_Depart <- 0  # to determine when the bird departs for day
    m_Arrive <- 0  # to determine when the bird arrives for day
  for (t in 1:length(time)) {  # for each timestep of day
    datetime <- f_date[t, "datetime"]  # pulls out datetime
  if (nest_age == 1 && t == 1) {  # intitial data and location  
    f_date[1, "forage_suc"] <- as.integer(0) 
    m_date[1, "forage_suc"] <- as.integer(0)
    f_date[1, "behavior"] <- "nest"
    m_date[1, "behavior"] <- "roost"
    f_date[1, "lat"] <- as.numeric(nest$attributes['lat'])
    m_date[1, "lat"] <- as.numeric(nest$attributes['lat'])
    f_date[1, "long"] <- as.numeric(nest$attributes['long'])
    m_date[1, "long"] <- as.numeric(nest$attributes['long'])
    f_date[1, "alt"] <- as.integer(nest$attributes['alt']) + 15
    m_date[1, "alt"] <- as.integer(nest$attributes['alt']) + 15
  } else {
  ## Bring in Previous Location
  if (t == 1) {  # for first location of day
    f_prev_loc <- tail(female[female$date == d-1, ], 1)
    m_prev_loc <- tail(male[male$date == d-1, ], 1)
  }
  if(t != 1) {  # for all other locations, brings in previous location data 
    f_prev_loc <- f_date[t-1, ]
    m_prev_loc <- m_date[t-1, ]
  }  
  # Bring in Current Location
  f_current_loc <- f_date[t, ]
  m_current_loc <- m_date[t, ]
  # Select Behavior for female         
  # First location of day
  if (f_current_loc$timestep_after_start == 0) { 
    f_behavior <- f_prev_loc$behavior
  } else {
  # After first location
  if (f_Depart == 0) { 
    dep <- PredictDeparture(time=f_current_loc$time_after_start,
    pars=sim_pars[['female']][['depart_pars']])           
  if (dep == 0) {
    f_behavior <- f_prev_loc$behavior
  }
  if (dep == 1) {  
    f_behavior <- PredictDepartBehavior(time_p = f_current_loc$time_proportion)
    f_Depart <- 1
  }
  }
  if (f_Arrive == 1) {
    f_behavior <- "nest" 
  }
  if (f_Depart == 1 && f_Arrive == 0) {
    arr <- PredictArrival(time_m=f_current_loc$time_before_end,
      pars=sim_pars[['female']][['arrive_pars']])
  if (arr==0) {
    f_behavior <- PredictDayBehavior(time_p=f_current_loc$time_proportion,  
                              pars=sim_pars[["female"]])
  }
  if (arr==1) {
    f_behavior <- "nest"  # will need to force this to be nest or roost
    f_Arrive <- 1
  }
  }
  if (f_current_loc$timestep_before_end == 0 && f_behavior != "nest") {
    f_behavior <- "nest"  # will need to force this to be nest or roost
  }
  }
  # Select behavior for male
  # First location of day
  if (m_current_loc$timestep_after_start == 0) { 
    m_behavior <- m_prev_loc$behavior
  } else {
  # After first location
  if (m_Depart == 0) { 
    dep <- PredictDeparture(time = m_current_loc$time_after_start,
      pars = sim_pars[['male']][['depart_pars']])           
  if (dep == 0) {
    m_behavior <- m_prev_loc$behavior
  }
  if (dep == 1) {  
    m_behavior <- PredictDepartBehavior(time_p = m_current_loc$time_proportion)
    m_Depart<-1
  }
  }
  if (m_Arrive == 1){
    m_behavior<-"roost" 
  }
  if (m_Depart == 1 && m_Arrive == 0) {
    arr <- PredictArrival(time_m=m_current_loc$time_before_end, 
      pars=sim_pars[['male']][['arrive_pars']])
  if (arr == 0) {
    m_behavior <- PredictDayBehavior(time_p=m_current_loc$time_proportion, 
      pars=pars[["male"]])
  }
  if (arr == 1) {
    m_behavior <- "roost"  # needed to force this to be nest or roost
    m_Arrive <- 1
  }
  }
  if (m_current_loc$timestep_before_end == 0 && m_behavior != "roost") {
    m_behavior <- "roost(forced)"
  }
  }
  # Select locations
  if (f_behavior == "nest") {
    f_lat <- as.numeric(nest$attributes['lat'])
    f_long <- as.numeric(nest$attributes['long'])
    f_alt <- as.integer(nest$attributes['alt']) + 15
  }
  if (m_behavior == "nest") {
    m_lat <- as.numeric(nest$attributes['lat'])
    m_long <- as.numeric(nest$attributes['long'])
    m_alt<-as.integer(nest$attributes['alt']) + 15
  }
  if (f_behavior == "roost") {
    f_lat <- f_prev_loc$lat
    f_long <- f_prev_loc$long
    f_alt <- f_prev_loc$alt
  }
  if(m_behavior == "roost") {
    m_lat <- m_prev_loc$lat
    m_long <- m_prev_loc$long
    m_alt <- as.integer(m_prev_loc$alt + rnorm(1, 5, 5))
  }
  if(f_behavior == "cruise") {
    f_lat <- f_prev_loc$lat + rnorm(1, 0, .01)
    f_long <- f_prev_loc$long + rnorm(1, 0, .01)
    f_alt <- as.integer(f_prev_loc$alt + rnorm(1, 5, 5))
  }
  if (m_behavior == "cruise") {
    m_lat <- m_prev_loc$lat + rnorm(1, 0, .01)
    m_long <- m_prev_loc$long + rnorm(1, 0, .01)
    m_alt <- as.integer(m_prev_loc$alt + rnorm(1, 5, 5))
  }
  if (f_behavior == "forage") {
    f_lat <- f_prev_loc$lat+rnorm(1, 0, .01)
    f_long <- f_prev_loc$long+rnorm(1, 0, .01)
    f_alt <- as.integer(f_prev_loc$alt+rnorm(1, 5, 5))
  }
  if (m_behavior == "forage") {
    m_lat <- m_prev_loc$lat+rnorm(1, 0, .01)
    m_long <- m_prev_loc$long+rnorm(1, 0, .01)
    m_alt <- as.integer(m_prev_loc$alt + rnorm(1, 5, 5))
  }

  if(f_behavior == "territorial") {
    f_lat <- f_prev_loc$lat + rnorm(1, 0, .01)
    f_long <- f_prev_loc$long + rnorm(1, 0, .01)
    f_alt <- as.integer(f_prev_loc$alt + rnorm(1, 5, 5))
  }
  if (m_behavior == "territorial") {
    m_lat <- m_prev_loc$lat + rnorm(1 ,0 ,.01)
    m_long <- m_prev_loc$long + rnorm(1, 0, .01)
    m_alt <- as.integer(m_prev_loc$alt + rnorm(1, 5, 5))
  }    
  if (f_behavior == "loaf") {
    f_lat <- f_prev_loc$lat
    f_long <- f_prev_loc$long
    f_alt <- as.integer(f_prev_loc$alt)
  }
  if (m_behavior == "loaf") {
    m_lat <- m_prev_loc$lat
    m_long <- m_prev_loc$long
    m_alt <- as.integer(m_prev_loc$alt)
  }
  # Write data to object
  f_date[f_date$datetime == datetime, "forage_suc"] <- f_prev_loc$forage_suc + 
  as.integer(difftime(f_current_loc$datetime, f_prev_loc$datetime, 
                      units="min")) 
  m_date[m_date$datetime == datetime, "forage_suc"] <- m_prev_loc$forage_suc + 
  as.integer(difftime(m_current_loc$datetime, m_prev_loc$datetime, 
                      units="min")) 
  f_date[f_date$datetime == datetime, "behavior"] <- f_behavior
  m_date[m_date$datetime == datetime, "behavior"] <- m_behavior
  f_date[f_date$datetime == datetime, "lat"] <- f_lat
  m_date[m_date$datetime == datetime, "lat"] <- m_lat
  f_date[f_date$datetime == datetime, "long"] <- f_long
  m_date[m_date$datetime == datetime, "long"] <- m_long
  f_date[f_date$datetime == datetime, "alt"] <- f_alt
  m_date[m_date$datetime == datetime, "alt"] <- m_alt 
  # Delete behaviors
  rm(f_behavior)
  rm(m_behavior)
  }
  # After completing f_date/m_date, write to female/male
  sv <- female$date == d
  female[sv, ] <- f_date
  sv <- male$date == d
  male[sv, ] <- m_date
  }
  }
  simlist[[3]][[i]][[3]][[1]][[4]] <- female
  simlist[[3]][[i]][[3]][[2]][[4]] <- male
  }
  return(simlist)
}

# SummarizeSimlistNests Function -----------------------------------------------

###  Creates a list of the nests
###  Usage: SummarizeSimlistNests(simlist)
###  Arguments: ls = simlist
###  Returns: a lists of all the nest data
###  Notes: 
###  Blake Massey
###  2014.05.11

SummarizeSimlistNests<- function(simlist = simlist){
  all_nests<-simlist()
  for (i in 1:length(ls[[3]])){
    nest_f_name<-ls[[3]][[i]][[3]][[1]][["bird_name"]] #nest name
    nest_f_df<-ls[[3]][[i]][[3]][[1]][[4]]
    nest_f_df$sex<-"f"    
    nest_f_df<-nest_f_df[,c("sex",setdiff(names(nest_f_df),"sex"))] 
    nest_m_name<-ls[[3]][[i]][[3]][[2]][["bird_name"]] #nest name
    nest_m_df<-ls[[3]][[i]][[3]][[2]][[4]]
    nest_m_df$sex<-"m"
    nest_m_df<-nest_m_df[,c("sex",setdiff(names(nest_m_df),"sex"))] 
    nests<-list(nest_f_df, nest_m_df)
    names(nests)<-c(nest_f_name, nest_m_name) 
    all_nests<-c(all_nests, nests)
  }
  return(all_nests)
}

# SummarizeSimlist Function ----------------------------------------------------

###  Creates summary statistics for a simlist
###  Usage: SummarizeSimlist()
###  Arguments: simlist = simlist
###  Returns: a dataframe with the summary statistics
###  Notes: 
###  Blake Massey
###  2014.05.11

SummarizeSimlist <- function(simlist = simlist) {
  sumstats <- data.frame(NULL)
  for (i in 1:length(simlist[[3]])){
    sumstats[i,"nest"] <- simlist[[3]][[i]][[1]][[2]]  # nest name  
    sumstats[i, "male_locs"] <- length(simlist[[3]][[i]][[3]][[1]][[4]][[1]])
      # number of male locations
    sumstats[i, "female_locs"] <- length(simlist[[3]][[i]][[3]][[1]][[4]][[2]])
      # number of female locations
    sumstats[i, "start_date"] <- as.character(min(simlist[[3]][[i]][[2]][[1]]))
      # start date for sim
    sumstats[i, "end_date"] <- as.character(max(simlist[[3]][[i]][[2]][[1]])) 
      # end date for sim
    sumstats[i, "date_diff"] <- max(simlist[[3]][[i]][[2]][[1]]) - 
      min(simlist[[3]][[i]][[2]][[1]])
  }
  return(sumstats)
}
