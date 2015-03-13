# --- SIMULATION FUNCTIONS -----------------------------------------------------
# General functions for simulation construction and analysis of simulation
# Added a line
# ------------------------------------------------------------------------------

# CreateSumIntervals Function -------------------------------------------------

###  Creates summary intervals for the simulation input
###  Usage: CreateSumIntervals(sim_start, sim_period, sim_end, sum_period,
###    sum_interval, sum_interval_custom)
###  Arguments: sim_start = a POSIXct object for simulation's start, required
###             sim_period = a Period object (from 'lubridate' package), e.g.
###               period(10, "year"), period(10, "week"), period(10, "week")
###               that sets the length of time the simulation runs
###             sim_end = a POSIXct object for simulation's end. Setting this 
###               will override the 'sim_period' parameter. Default = NULL.
###             sum_period = a Period object that sets the summary interval  
###             sum_interval = vector or dateframe (with date in first column) 
###               of the format %B%d, %b%d, %d%B, or %B%d (see ?strptime for 
###               more details on formatting). Example: c("April01", "May15", 
###               "Sep1", "Oct15"). Overrides 'sum_period' parameter.
###             sum_interval_custom = a vector of POSIXct objects that set the 
###               start of the summary periods. If sim_start precedes the first
###               value, the first value will go until from start_sim to 
###               sum_interval_custom. Last interval will go from the last 
###               date in sum_interval_custom to end of sim_period or end_sim. 
###               Overrides 'sum_period' or 'sum_interval' parameters.
###  Returns: list of intervals  
###  Notes: Either 'sim_period' or 'sim_end' must be set, and either 
###    'sum_period', 'sum_interval', or 'sum_interval_custom' must be set. For
###    'sum_interval' the dates can be set as numeric for both day and month, 
###    but there is a possibility that the day and month will be inverted 
###    because dates with both values <12 can be confused (e.g. 02-10 can be 
###    either February 10th or October 2nd. Date format: "10Feb" or "02Oct" is 
###    preferred.  
###  Blake Massey
###  2015.03.10

CreateSumIntervals <- function(sim_start,
                               sim_period = NULL,
                               sim_end = NULL,
                               sum_period = NULL,
                               sum_interval = NULL,
                               sum_interval_custom = NULL) {
  suppressPackageStartupMessages(require(lubridate))
  if (is.null(sim_period) && is.null(sim_end)) {
    stop("must provide either 'sim_period' or 'sim_end' parameter value")
  }
  if (!is.null(sim_period))  sim_end <- sim_start + sim_period   
  if (!is.null(sum_period) && !is.null(sum_interval)){
    print("'sum_period' used instead of 'sum_interval'")
  }
  if (!is.null(sum_period)) {
    sum_period_unit <- ExtractUnitFromPeriod(sum_period) 
    first_int <- as.interval(sim_start, ceiling_date(sim_start+period(1, 
      "second"), unit=sum_period_unit))
    sum_intervals <- sim_start
    end_int <- first_int@start + first_int@.Data    
    sum_intervals <- with_tz(append(sum_intervals, end_int), tz(sim_start))      
    while (end_int < sim_end) {  
      end_int <- end_int + sum_period
      sum_intervals <- with_tz(append(sum_intervals, end_int), tz(sim_start))   
    }  
    if (tail(sum_intervals, 1) > sim_end) {
      sum_intervals[length(sum_intervals)] <- sim_end
    }
    sum_interval <- NULL
  } 
  if (!is.null(sum_interval)){
      if (is.data.frame(sum_interval)) sum_interval <- sum_interval[, 1]
      if (max(names(guess_formats(sum_interval, c("md", "dm")))) == "md") {
        md_format <- max(guess_formats(sum_interval, c("md", "dm")))
        sum_interval <- as.character(as.Date(sum_interval, format = md_format),
          "%d%b")  # date must be in day month order
      }                      
      dates <- dmy(paste0(sum_interval, 2000)) # creates POSIXct    
      sum_interval  <- sum_interval[order(dates)]  # orders sum_interval dates 
      first_sum_int <- FindFirstSumInterval(sim_start, sum_interval)
      end_int <- first_sum_int@start+first_sum_int@.Data 
      sum_intervals <- with_tz(append(sim_start, end_int), tz(sim_start))
      while (end_int < sim_end) {  
        end_int_jul <- yday(end_int)
        sum_int_jul <- yday(dmy(paste0(sum_interval, year(end_int))))
        sum_int_position <- match(end_int_jul, sum_int_jul)
        if (sum_int_position == length(sum_interval)) {
          next_int <- new_interval(end_int, dmy(paste0(sum_interval[1], 
            (year(end_int)))) + period(1, "year"))
        } else {
          next_int <- new_interval(end_int, dmy(paste0(sum_interval
            [sum_int_position+1], year(end_int))))
        }
        end_int <- next_int@start + next_int@.Data
        sum_intervals <- with_tz(append(sum_intervals, end_int), tz(sim_start))   
      }      
      if (tail(sum_intervals, 1) > sim_end) {
        sum_intervals[length(sum_intervals)] <- sim_end
      }
      sum_interval_custom <- NULL
  }
  if(!is.null(sum_interval_custom)){
    if (is.data.frame(sum_interval)) sum_interval <- sum_interval_custom[,1]
    sum_interval <- sum_interval_custom # needed if sum_interval_custom != df
    sum_intervals <- as.POSIXct(sum_interval, tz=tz(sim_start))
    if (tail(sum_intervals, 1) > sim_end && !is.null(sum_period)) {
      sum_intervals[length(sum_intervals)] <- sim_end
      warning("'sum_inverval_custom' exceeds 'sim_period' time length")
    }
    if (tail(sum_intervals, 1) > sim_end && is.null(sum_period)) {
      sum_intervals[length(sum_intervals)] <- sim_end
      warning("'sum_inverval_custom' exceeds 'sim_end' time length")
    }
    if (tail(sum_intervals, 1) < sim_end) {
      sum_intervals[length(sum_intervals)+1] <- sim_end
    }
  }
  sum_intervals_list <- list()
  for (i in 1:(length(sum_intervals)-1)) {
    interval <- as.interval(sum_intervals[i], sum_intervals[i+1])
    sum_intervals_list[[length(sum_intervals_list)+1]] <- interval
  }
  return(sum_intervals_list)
}

# ExtractFromPeriod Function ---------------------------------------------------

###  Helper function for CreateSumInterval(). Finds unit of Period object. 
###  Usage: ExtractUnitFromPeriod(period)
###  Arguments: sim_start = a POSIXct object for simulation's start, required
###  Returns: period object's unit as a chararcter 
###  Notes:
###  Blake Massey
###  2015.03.10

ExtractUnitFromPeriod <- function(period) {
  period <- period
  if (period@year > 0) unit <- "year"
  if (period@month > 0) unit <- "month"
  if (period@day == 7 ) unit <- "week"  
  if (period@day > 0 && period@day < 7) unit <- "day"
  if (period@day > 7) unit <- "day"
  if (period@hour > 0) unit <- "hour"  
  if (period@minute > 0) unit <- "minute"
  if (period@.Data > 0) unit <- "second"
  return(unit)
}

# FindFirstSumInterval Function ------------------------------------------------

###  Helper function for CreateSumInterval. Finds first interval based on  
###    sim_start datetime and sum_interval dates
###  Usage: FindFirstSumInterval(sim_start, sum_interval)
###  Arguments: sim_start = a POSIXct object for simulation's start, required
###             sum_interval = vector or dateframe (with date in first column) 
###               of the format %B%d, %b%d, %d%B, or %B%d (see ?strptime for 
###               more details on formatting). Example: c("April01", "May15", 
###               "Sep1", "Oct15"). Required.
###  Returns: an Interval object
###  Notes:
###  Blake Massey
###  2015.03.10

FindFirstSumInterval <- function(sim_start, 
                                 sum_interval){
  require(lubridate)
  start_year <- as.Date(0, origin=as.Date(floor_date(sim_start, "year")))
  if (length(sum_interval) == 1) {
    first_int <- as.interval(period(1, "year"), dmy(paste0(sum_interval[1], 
      year(start_year))))
    if (sim_start %within% first_int) { 
      first_sum_int <- new_interval(sim_start, dmy(paste0(sum_interval[1], 
        year(start_year + period(1, "year")))))
    } else {
      first_sum_int <- new_interval(sim_start, dmy(paste0(sum_interval[1], 
        year(start_year))))    
    }
  }
  if (length(sum_interval) == 2) {
    first_int <- as.interval(dmy(paste0("0101", year(start_year))),
      dmy(paste0(sum_interval[1], year(start_year))))
    if ((sim_start + period(1, "second")) %within% first_int) { 
      first_sum_int <- new_interval(sim_start, dmy(paste0(sum_interval[1], 
        year(start_year)))) 
    } 
    if (exists("first_sum_int") == FALSE) {
      second_int <- as.interval(dmy(paste0(sum_interval[1], year(start_year))),
        dmy(paste0(sum_interval[2], year(start_year))))
      if ((sim_start + period(1, "second")) %within% second_int) { 
        first_sum_int <- new_interval(sim_start, dmy(paste0(sum_interval[2], 
          year(start_year))))  
      } else { 
        first_sum_int <- new_interval(sim_start, dmy(paste0(sum_interval[1], 
          year(start_year+period(1, "year")))))   
      }
    }
  }
  if (length(sum_interval) > 2) { 
    first_int <- as.interval(dmy(paste0("0101", year(start_year))),
        dmy(paste0(sum_interval[1], year(start_year))))
    if ((sim_start + period(1, "second")) %within% first_int) { 
      first_sum_int <- new_interval(sim_start, dmy(paste0(sum_interval[1], 
        year(start_year))))
    }
    for (i in 2:(length(sum_interval)-1)) { 
      mid_int <- as.interval(dmy(paste0(sum_interval[i], year(start_year))),
        dmy(paste0(sum_interval[i+1], year(start_year))))
      if ((sim_start + period(1, "second")) %within% mid_int) { 
        first_sum_int <- new_interval(sim_start, dmy(paste0(sum_interval[i+1], 
          year(start_year))))  
      }
    }
    for (i in length(sum_interval)) { 
      last_int <- as.interval(dmy(paste0(sum_interval[i], year(start_year))),
        dmy(paste0("0101", (year(start_year)+1) )))
      if ((sim_start + period(1, "second")) %within% last_int) { 
        first_sum_int <- new_interval(sim_start, dmy(paste0(sum_interval[1], 
          (year(start_year)+1))))   
      }
    }
  }
  return(first_sum_int)
}

# CreateStepIntervals Function --------------------------------------------------

### Create step intervals within a summary interval
### Usage:     CreateStepIntervals(current_sum_interval = current_sum_interval,step_period = step_period)
### Arguments: current_sum_interval = an Interval object from lubridate that represents the
###             current summary interval
###            step_period = a Period object specifying the length of the step period. The
###             step_period should be of equal or greater length than the time_step_period. Note that if
###             current_sum_interval is not evenly divisible by step_period the last step_interval is removed
###             so that all step_intervals returned are of equal duration.
### Returns: a list of intervals
### Notes:
### Javan Bauder
### 2015.03.13

CreateStepIntervals <- function(current_sum_interval = current_sum_interval,
                              step_period = step_period) {
  suppressPackageStartupMessages(require(lubridate))
  step_intervals <- list()
  interval_counter <- 1
  current_start <- int_start(current_sum_interval)
  current_end <- (current_start + step_period)
  stop_point<-int_end(current_sum_interval)
  while(current_start < (stop_point)){
    current_end <- (current_start+step_period)
    step_intervals[[interval_counter]] <- new_interval(current_start,current_end)
    interval_counter <- interval_counter + 1
    current_start <- current_start + step_period
  }
  if(int_end(step_intervals[[length(step_intervals)]]) > stop_point){
    step_intervals[[length(step_intervals)]] <- NULL
  }
  return(step_intervals)
}



# CreateTimeSteps Function --------------------------------------------------

### Create the time steps within a step interval
### Usage:     CreateTimeSteps(current_step_interval = current_step_interval,time_step_period = time_step_period)
### Arguments: current_step_interval = an Interval object from lubridate that represents the
###             current step interval
###            time_step_period = a Period object specifying the length of the time step period. Note that if
###             the length  of the current_step_interval is equal to the time_step_period the function will return
###             a POSIXct object with the same date as the end of the current_step_interval.
### Returns: a list of intervals
### Notes:
### Javan Bauder
### 2015.03.13

CreateTimeSteps <- function(current_step_interval=current_step_interval,
                              time_step_period=time_step_period){
  suppressPackageStartupMessages(require(lubridate))
  options(lubridate.verbose=FALSE)
  time_steps <- list()
  step_start <- int_start(current_step_interval) + time_step_period
  time_steps[[1]] <- step_start
  stop_point <- int_end(current_step_interval)
  seconds_to_add <- period_to_seconds(time_step_period)
  t.step <- 2
  while(time_steps[[t.step-1]] != (stop_point)){
    time_steps[[t.step]] <- (step_start + seconds_to_add)
    t.step <- t.step + 1
    seconds_to_add <- seconds_to_add + period_to_seconds(time_step_period)
  }
  return(time_steps)
}

# CreateTimeStepsv2 Function --------------------------------------------------

### Create the time steps within a step interval
### Usage:     CreateTimeStepsv2(current_step_interval = current_step_interval,time_step_period = time_step_period)
### Arguments: current_step_interval = an Interval object from lubridate that represents the
###             current step interval
###            time_step_period = a Period object specifying the length of the time step period. Note that if
###             the length  of the current_step_interval is equal to the time_step_period the function will return
###             a list of length 1. In CreateTimeStepsv2, when the length  of the current_step_interval is 
###             equal to the time_step_period, the date of the POSIXct object of the start date of the 
###             current_step_interval.
### Returns: a list of intervals
### Notes:
### Javan Bauder
### 2015.03.13

CreateTimeStepsv2 <- function(current_step_interval=current_step_interval,
                       time_step_period=time_step_period){
  suppressPackageStartupMessages(require(lubridate))
  options(lubridate.verbose=FALSE)
  if(int_length(current_step_interval) == period_to_seconds(time_step_period)){
    time_steps <- int_start(current_step_interval)
  } else {
    time_steps <- list()
    step_start <- int_start(current_step_interval) + time_step_period
    time_steps[[1]] <- step_start
    stop_point <- int_end(current_step_interval)
    seconds_to_add <- period_to_seconds(time_step_period)
    t.step <- 2
    while(time_steps[[t.step-1]] != (stop_point)){
      time_steps[[t.step]] <- (step_start + seconds_to_add)
      t.step <- t.step + 1
      seconds_to_add <- seconds_to_add + period_to_seconds(time_step_period)
    }  
  }
  return(time_steps)
}
