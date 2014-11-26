# --- PARAMETERS FUNCTIONS -----------------------------------------------------
#  Functions for telemetry data and simulation parameter creation, extract, and 
#  plotting.
# ------------------------------------------------------------------------------

# CreateListIndex Function -----------------------------------------------------

###  Creates a "indexed" list by changing all the values into sequential numbers   
###  Usage: CreateListIndex(list)
###  Arguments: list = input list  
###  Returns: list with sequential numbers for all of the values    
###  Notes: Works with up to 3 nested levels (lists, sublists, subssublists)
###  Blake Massey
###  2014.10.14

CreateListIndex <- function(list){
  index_count <- 0
  for (i in 1:length(list)) {
    sublist  <- list[[i]]
    if (length(sublist) > 1) { 
      for (j in 1:length(sublist)) { 
        subsublist <- sublist[[j]]
        if (length(subsublist) > 1) {  
          for (k in 1:length(subsublist)) {
            index_count <- index_count + 1
            list[[i]][[j]][k] <- index_count
          }
        } else {
          index_count <- index_count + 1
          list[[i]][[j]] <- index_count    
        }
     }
     } else {
       index_count <- index_count + 1
       list[[i]] <- index_count    
    }
  }
  return(list)
}

# CreateNewPars Function ----------------------------------------------------

###  Creates a new list of behavior parameters with either NA or index values
###  Usage: CreateNewPars(index, random)
###  Arguments: index = logical, indicates if the parameter values should be 
###               numbers indicating their relative position. Default is TRUE.
###             random = logical, indicate if the parameter values should be
###               filled in using runif(1,1,2), random=TRUE overrides index=TRUE 
###               Default is FALSE.
###  Returns: a newly generated list of behavior parameters with NA or index 
###    numbers for each of the values        
###  Notes: Used to create a new behavior parameters list from scratch. 
###    Uses ExtractParsIndex().
###  Blake Massey
###  2014.10.14

CreateNewPars <- function(index=TRUE,
                          random=FALSE){
  pars <- list(arrive=list(shape=NA, scale=NA), 
               depart=list(shape=NA, scale=NA),
               cruise=list(shape1=NA, shape2=NA),
               forage=list(shape1=NA, shape2=NA),
               flight=list(shape=NA),
               loaf=list(shape1=NA, shape2=NA),
               nest=list(shape1=NA, shape2=NA),
               step=list(location=NA, scale=NA, shape=NA),
               territorial=list(shape1=NA, shape2=NA)
               )
  female<-list(female=pars)
  male<-list(male=pars)
  pars<-c(male,female)
  if (random == TRUE){
    pars <- rapply(pars, f=function(x) ifelse(is.na(x),runif(1,1,2),x), 
      how="replace" )
  }
  if (index == TRUE && random == FALSE){
    pars <- CreateListIndex(pars)
  }
  return(pars)
}

# ExtractParsBySex Function ---------------------------------------------------

###  Returns a dataframe of pars values 
###  Usage: ExtractParsBySex(pars)
###  Arguments: pars =  parameters
###             sex = sex to return: "male", "female", or NULL(both), default is 
###               NULL.
###  Returns: dataframe of behavior parameters      
###  Notes: 
###  Blake Massey
###  2014.10.14

ExtractParsBySex <- function(list,
                             sex = NULL){
  list <- list
  list2 <- sapply(list, FUN = function(X) unlist(X))
  df <- as.data.frame(list2)
  if (!is.null(sex) && sex == "male") df$female <- NULL
  if (!is.null(sex) && sex == "female") df$male <- NULL
  print(df)
  return(df)
}

# ExtractParsMatrix Function ---------------------------------------------------

###  Extracts a dataframe of parameter values with columns for "sex" and 
###    "behavior"  
###  Usage: ExtractParsMatrix(pars)
###  Arguments: pars = list of parameters
###  Returns: dataframe of behavior parameters  
###  Notes: Useful for comparing parameter values among behaviors
###  Blake Massey
###  2014.10.14

ExtractParsMatrix <- function(list){
  suppressPackageStartupMessages(require(plyr))
  df <- data.frame()
  for (i in 1:length(list)){  
    sublist  <- list[[i]]
    df_sub <- rbind.fill(lapply(sublist, as.data.frame))
    behavior <- names(sublist)
    sex <- (rep(names(list)[i], length(behavior)))
    df <- rbind(df,cbind(sex, behavior, df_sub))  
  }
  print(df)
  return(df)
}

# ExtractRoostPars Function ----------------------------------------------------

###  Returns df of roost parameters for male/female, arrive/depart, shape/scale 
###  Usage: ExtractRoostPars(pars)
###  Arguments: pars = simulation parameters with male/female, 
###               arrive/depart, shape/scale 
###  Returns: df of roost parameters for male/female, arrive/depart, shape/scale         
###  Notes: used in PlotRoostECDF
###  Blake Massey
###  2014.06.03

ExtractRoostPars <- function(pars = sim_pars){
  sex <- rep(c("f", "m"),each=2) 
  roost <- rep(c("depart", "arrive"),times=2) 
  shape <- rep(NA,length=4) 
  scale <- rep(NA,length=4) 
  output <- data.frame(sex, roost, shape, scale, stringsAsFactors=FALSE)
  for (i in 1:nrow(output)){
  row <- output[i,]
  sex <- row[, "sex"]
  p_sex <-ifelse(sex=="f", "female", "male")
  roost <- row [,"roost"]
  par <- paste(roost, "_pars", sep="")
  output[i, "shape"] <- 
    pars[[p_sex]][[par]][['shape']]      
  output[i, "scale"] <- 
    pars[[p_sex]][[par]][['scale']]    
  }  
  return(output)
}

# FitArrival Function --------------------------------------------------------

###  Fit Weibull parameters to roost arrival data
###  Usage: FitArrival(data, shape)
###  Arguments: data = dataframe with a "datetime", "hr_after_sunset" and 
###               "behavior" column
###             shape = starting value of shape, default is 10.
###  Returns: list(arr_scale, arr_shape)
###  Notes:
###  Blake Massey
###  2014.05.11

FitArrival <- function(data = data, 
                       shape = 10) {
  suppressPackageStartupMessages(require(bbmle))
  data <- subset(data, behavior=="arrive")
  data$arr_diff_min <- difftime(data$hr_after_sunset, data$datetime)
  data$arr_diff_min <- as.integer(as.numeric(data$arr_diff, units = "mins"))
  arr_mean <- mean(data$arr_diff_min) #calculate mean of arrival differences
  arrive_pars <- mle2(NLLWeibull, start=list(shape=shape, mu=arr_mean), 
    data=list(data=data$arr_diff_min)) #calculate weibull parameters
  arr_scale <- coef(arrive_pars)[2] / gamma(1 + (1 / coef(arrive_pars)[1]))
  names(arr_scale) <- "scale"
  arr_shape <- coef(arrive_pars)[1]
  names(arr_shape) <- "shape"
  z <- c(arr_scale, arr_shape)
  return(z)
}

# FitBehavior Function --------------------------------------------------------

###  Fit Beta parameters to behavior data
###  Usage: FitBehavior(df, shape1, shape2)
###  Arguments: df = dataframe with a "datetime", "hr_after_sunset" and 
###               "behavior" column
###             shape1 = starting value of shape1, default is 1.
###             shape2 = starting value of shape2, default is 1.
###  Returns: list(shape1, shape2)
###  Notes:
###  Blake Massey
###  2014.10.06

FitBehavior <- function(data, 
                        shape1 = 1,
                        shape2 = 1) {
  suppressPackageStartupMessages(require(bbmle))
  suppressPackageStartupMessages(require(plyr))
  suppressPackageStartupMessages(require(dplyr))
  suppressPackageStartupMessages(require(reshape2))
  suppressPackageStartupMessages(require(scales))
  if ("package:raster" %in% search()){
    detach("package:raster", unload=TRUE)  # masks select() in 'dplyr'
  }
  source('C:/Work/R/Functions/gps.R')
  load("C:/Work/R/Data/Simulation/blank_pars.RData")
  list_pars <- blank_pars
  breaks = 10
  df <- data
  df$sex <- gsub("m", "male", df$sex)
  df$sex <- gsub("f", "female", df$sex)
  df$behavior <- factor(df$behavior)
  source('C:/Work/R/Functions/gps.R')
  behavior_colors <- CreateColorsByBehavior(output=TRUE)
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
  cruise <- melted %>% filter(behavior == "cruise", sex=="male") %>% 
    select(bins_mid, value)
  cruise_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2), 
    data=list(data=cruise$value))   
  list_pars[["male"]][["cruise_pars"]][["shape1"]] <- coef(cruise_mle)[1]
  list_pars[["male"]][["cruise_pars"]][["shape2"]] <- coef(cruise_mle)[2]
  
  forage <- melted %>%
    filter(behavior == "forage", sex=="male") %>% select(bins_mid, value) 
  forage_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2), 
    data=list(data=forage$value)) 
  list_pars[["male"]][["forage_pars"]][["shape1"]] <- coef(forage_mle)[1]
  list_pars[["male"]][["forage_pars"]][["shape2"]] <- coef(forage_mle)[2]
  
  loaf <- melted %>%
    filter(behavior == "loaf", sex=="male") %>% select(bins_mid, value) 
  loaf_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2), 
    data=list(data=loaf$value))  
  list_pars[["male"]][["loaf_pars"]][["shape1"]] <- coef(loaf_mle)[1]
  list_pars[["male"]][["loaf_pars"]][["shape2"]] <- coef(loaf_mle)[2] 
  
  nest <- melted %>%
    filter(behavior == "nest", sex=="male") %>% select(bins_mid, value) 
  nest_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2), 
    data=list(data=nest$value))  
  list_pars[["male"]][["nest_pars"]][["shape1"]] <- coef(nest_mle)[1]
  list_pars[["male"]][["nest_pars"]][["shape2"]] <- coef(nest_mle)[2]
    
  roost <- melted %>%
    filter(behavior == "roost", sex=="male") %>% select(bins_mid, value) 
  roost_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2), 
    data=list(data=roost$value), optimizer = "nlminb")  
  list_pars[["male"]][["roost_pars"]][["shape1"]] <- coef(roost_mle)[1]
  list_pars[["male"]][["roost_pars"]][["shape2"]] <- coef(roost_mle)[2] 
  
  cruise <- melted %>%
    filter(behavior == "cruise", sex=="female") %>% select(bins_mid, value) 
  cruise_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2), 
    data=list(data=cruise$value))  
  list_pars[["female"]][["cruise_pars"]][["shape1"]] <- coef(cruise_mle)[1]
  list_pars[["female"]][["cruise_pars"]][["shape2"]] <- coef(cruise_mle)[2]
  
  forage <- melted %>%
    filter(behavior == "forage", sex=="female") %>% select(bins_mid, value) 
  forage_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2), 
    data=list(data=forage$value))
  list_pars[["female"]][["forage_pars"]][["shape1"]] <- coef(forage_mle)[1]
  list_pars[["female"]][["forage_pars"]][["shape2"]] <- coef(forage_mle)[2]
  
  loaf <- melted %>%
    filter(behavior == "loaf", sex=="female") %>% select(bins_mid, value) 
  loaf_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2), 
    data=list(data=loaf$value)) 
  list_pars[["female"]][["loaf_pars"]][["shape1"]] <- coef(loaf_mle)[1]
  list_pars[["female"]][["loaf_pars"]][["shape2"]] <- coef(loaf_mle)[2] 
  
  nest <- melted %>%
    filter(behavior == "nest", sex=="female") %>%
    select(bins_mid, value) 
  nest_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2), 
    data=list(data=nest$value)) #calculate weibull parameters 
  list_pars[["female"]][["nest_pars"]][["shape1"]] <- coef(nest_mle)[1]
  list_pars[["female"]][["nest_pars"]][["shape2"]] <- coef(nest_mle)[2]
    
  roost <- melted %>%
    filter(behavior == "roost", sex=="female") %>%
    select(bins_mid, value) 
  roost_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2), 
    data=list(data=roost$value), optimizer = "nlminb") 
  list_pars[["female"]][["roost_pars"]][["shape1"]] <- coef(roost_mle)[1]
  list_pars[["female"]][["roost_pars"]][["shape2"]] <- coef(roost_mle)[2] 
  
  return(list_pars) 
}

# FitDeparture Function --------------------------------------------------------

###  Fit Weibull parameters to roost departure data
###  Usage: FitDeparture(data, shape)
###  Arguments: df = dataframe with a "datetime", "hr_before_sunrise", and 
###               "behavior" column
###             shape = starting value of shape. Default is 10.
###  Returns: list(dep_scale, dep_shape)
###  Notes:
###  Blake Massey
###  2013.05.11

FitDeparture <- function(data=data, 
                         shape=10) {
  suppressPackageStartupMessages(require(bbmle))
  data <- subset(data, behavior=="depart")
  data$dep_diff_min <- difftime(data$datetime, data$hr_before_sunrise)
  data$dep_diff_min <- as.integer(as.numeric(data$dep_diff, units = "mins"))
  data <- subset(data, select = c("id", "date", "datetime", "behavior",
    "dep_diff_min"))
  dep_mean <- mean(data$dep_diff_min) # calculate mean of departure differences
  depart_pars <- mle2(NLLWeibull, start=list(shape=shape, mu=dep_mean), 
    data=list(data=data$dep_diff_min))  # calculate weibull parameters
  dep_scale <- coef(depart_pars)[2] / gamma(1 + (1 / coef(depart_pars)[1]))  
    # store departure weibull scale
  names(dep_scale)<-"scale"
  dep_shape <- coef(depart_pars)[1] #store arrive weibull shape
  z <- c(dep_scale, dep_shape)
  return(z)
}

# FitParetoParsToData Function -------------------------------------------------

###  Fits a generalized Pareto distribution to data
###  Usage: FitParetoParsToData(df, by, var, location, scale, shape)
###  Arguments: df = dataframe 
###             by = column to subset data 
###             var = variable to fit Pareto distribution
###             location = starting value of location, default is 
###               min(data$)-1
###             scale = starting value of scale, default is 1
###             shape = starting value of shape, default is 0
###  Returns: dataframe of "by" variable and parameters  
###  Notes:
###  Blake Massey
###  2014.11.06

FitParetoParsToData <- function(df, 
                                var,
                                by = NULL,
                                location = NULL,
                                scale = 1,
                                shape = 0) {
  suppressPackageStartupMessages(require(bbmle))
  suppressPackageStartupMessages(require(VGAM))
  df <- df
  ifelse(is.null(by), df$by <- "all", df$by <- df[,by])
  vars <- unique(df$by)
  pars <- data.frame(by=vars, location=NA, scale=NA, shape=NA, 
    stringsAsFactors=FALSE)
  if (is.null(location)) location <- min(df[,var],  na.rm = TRUE)-1  
  for (i in vars){
    data  <- subset(df, by == i, select = var, rm.na=TRUE)
    if (is.null(location)) location <- min(data,  na.rm = TRUE)-1 
    pars_i <- mle2(NLLPareto, start=list(location=location, scale=scale,
    shape=shape), data=list(data=data[,var]),  method="Nelder-Mead", 
    skip.hessian=TRUE) #calculates Pareto parameters
    pars[which(pars$by==i),"location"]=as.numeric(coef(pars_i)[1])
    pars[which(pars$by==i),"scale"]=as.numeric(coef(pars_i)[2])
    pars[which(pars$by==i),"shape"]=as.numeric(coef(pars_i)[3])
  }
  names(pars)[names(pars) == 'by'] <- by
  return(pars)
}

# FitParetoDensityToArray Function ---------------------------------------------

###  Fit Pareto Probability Density Function to an array of x values
###  Usage: FitParetoDensityToArray(df, var, pars, xlim)
###  Arguments: df = dataframe with data
###             var = variable fitted with Pareto distribution
###             pars = parameters of Pareto distribution, created by 
###               FitParetoParsToData()
###             xlim = sets xlim max, default is: max(df[,var])
###  Returns: density distributions of Pareto along an array of x values
###  Notes: x-axis is set to have 500 values, starting with: min(df[, var])+1
###  Blake Massey
###  2014.11.06

FitParetoDensityToArray <- function(df,
                                    var,
                                    pars, 
                                    xlim = NULL) {
  suppressPackageStartupMessages(require(stats))
  suppressPackageStartupMessages(require(VGAM))
  start <- min(df[, var])+1
  if (is.null(xlim)){
    end <- max(df[, var])
    x <- seq(start, end, length=500)
  } else {
    x <- seq(start, xlim, length=500)
  }
  dens <- data.frame(by=rep(pars[,1], each=length(x)), x=rep(x, 
    time=length(pars[,1])), y=NA, stringsAsFactors=FALSE)
  bys <- unique(pars[,1])
  for (i in bys){
  dens[which(dens$by == i), "y"] <- dgpd(x, location=pars[which(pars[,1] == i),
    "location"], scale=pars[which(pars[,1]==i),"scale"], shape=
    pars[which(pars[,1] == i),"shape"])
  }
  names(dens)[names(dens) == 'by'] <- names(pars)[1]
  return(dens)
}

# FitWrappedCauchyParsToData Function ------------------------------------------

###  Fits a Wrapped Cauchy distribution to data
###  Usage: FitWrappedCauchyParsToData(df, by, var, mu, rho)
###  Arguments: df = dataframe 
###             by = column to subset data
###             var = variable to fit Wrapped Cauchy distribution
###             mu = starting value of mu, default is 1
###             rho = starting value of rho, default is 0
###  Returns: dataframe of "by" variable and parameters  
###  Notes:
###  Blake Massey
###  2014.11.12

FitWrappedCauchyParsToData <- function(df, 
                                       var,
                                       by = NULL,
                                       mu = 1,
                                       rho = 0) {
  suppressPackageStartupMessages(require(CircStats))
  df <- df
  ifelse(is.null(by), df$by<-"all", df$by <- df[,by])
  vars <- unique(df$by)
  pars <- data.frame(by=vars, mu=NA, rho=NA, stringsAsFactors=FALSE)  
  for (i in vars){
    data  <- subset(df, by == i, select = var, rm.na=TRUE)
    data <- data[!is.na(data)]
    pars_i <- wrpcauchy.ml(data, mu=0, rho=.5, acc=1e-015)
    pars[which(pars$by==i), "mu"] = as.numeric(pars_i[1])
    pars[which(pars$by==i), "rho"] = as.numeric(pars_i[2])
  }
  names(pars)[names(pars) == 'by'] <- by
  return(pars)
}

# FitWrappedCauchyDensityToArray Function --------------------------------------

###  Fit Wrapped Cauchy Probability Density Function to an array of x values
###  Usage: FitWrappedCauchyDensityToArray(df, var, pars)
###  Arguments: var = variable fitted with Wrapped Cauchy distribution
###             pars = parameters of Wrapped Cauchy distribution, created by 
###               FitWrappedCauchyToData()
###  Returns: density distributions of Wrapped Cauchy along an array of x values
###  Notes: x-axis is set to have 359 values, starting with 1
###  Blake Massey
###  2014.11.12

FitWrappedCauchyDensityToArray <- function(var,
                                           pars) {
  suppressPackageStartupMessages(require(stats))
  suppressPackageStartupMessages(require(CircStats))
  x <- seq((2*pi/360), 2*pi, by=(2*pi/360))
  dens <- data.frame(by=rep(pars[,1], each=length(x)),
                   x=rep(x, time=length(pars[,1])),
                   y=NA, stringsAsFactors=FALSE)
  bys <- unique(pars[,1])
  for (i in bys){
    dens[which(dens$by==i),"y"] <- dwrpcauchy(x, mu=pars[which(pars[,1] == i),
    "mu"], rho=pars[which(pars[,1] == i),"rho"])
  }
  names(dens)[names(dens) == 'by'] <- names(pars)[1]
  return(dens)
}

# FitWrappedNormalParsToData Function ------------------------------------------

###  Fits a Wrapped normal distribution to data
###  Usage: FitWrappedNormalParsToData(df, by, var, mu, rho)
###  Arguments: df = dataframe 
###             by = column to subset data 
###             var = variable to fit Wrapped Normal distribution
###             mu = starting value of mu, default is 1
###             rho = starting value of rho, default is 0
###  Returns: dataframe of "by" variable and parameters  
###  Notes:
###  Blake Massey
###  2014.11.12

FitWrappedNormalParsToData <- function(df, 
                                       var,
                                       by = NULL,
                                       mu = 1,
                                       rho = 0) {
  suppressPackageStartupMessages(require(circular))
  df <- df
  ifelse(is.null(by), df$by<-"all", df$by <- df[,by])
  vars <- unique(df$by)
  pars <- data.frame(by=vars, mu=NA, rho=NA, stringsAsFactors=FALSE)  
  for (i in vars){
    data  <- subset(df, by == i, select = var, rm.na=TRUE)
    data <- data[!is.null(data)]
    suppressWarnings(pars_i <- mle.wrappednormal(data, mu=NULL, rho=NULL))
    pars[which(pars$by==i), "mu"] = as.numeric(pars_i$mu)
    pars[which(pars$by==i), "rho"] = as.numeric(pars_i$rho)
  }
  names(pars)[names(pars) == 'by'] <- by
  return(pars)
}

# FitWrappedNormalDensityToArray Function --------------------------------------

###  Fit Wrapped Normal Probability Density Function to an array of x values
###  Usage: FitWrappedNormalDensityToArray(df, var, pars)
###  Arguments: var = variable fitted with Wrapped Normal distribution
###             pars = parameters of Wrapped Normal distribution, created by 
###               FitWrappedNormalToData()
###  Returns: density distributions of Wrapped Normal along an array of x values
###  Notes: x-axis is set to have 359 values, starting with 1
###  Blake Massey
###  2014.11.12

FitWrappedNormalDensityToArray <- function(var,
                                           pars) {
  suppressPackageStartupMessages(require(stats))
  suppressPackageStartupMessages(require(circular))
  x <- seq((2*pi/360), 2*pi, by=(2*pi/360))
  dens <- data.frame(by=rep(pars[,1], each=length(x)),
                   x=rep(x, time=length(pars[,1])),
                   y=NA, stringsAsFactors=FALSE)
  bys <- unique(pars[,1])
  for (i in bys){
  dens[which(dens$by == i),"y"] <- 
    dwrappednormal(x, mu=pars[which(pars[,1] == i),"mu"], 
      rho=pars[which(pars[,1]==i),"rho"])
  }
  names(dens)[names(dens) == 'by'] <- names(pars)[1]
  return(dens)
}

# FitRoost Function ------------------------------------------------------------
###  Fits roost arrival and departure Weibull functions 
###  Usage: FitRoost(df, pars)
###  Arguments: df = dataframe of location data
###             pars = input parameters to have roost pars appended onto 
###  Returns: list of parameters          
###  Notes: uses ArriveFit() and DepartFit()
###  Blake Massey
###  2014.06.03

FitRoost <- function(df, 
                     pars=NULL){
  df_m <- subset(df, sex=="m")
  df_f <- subset(df, sex=="f")
  ifelse(!is.null(pars), pars <- pars, pars <- list())
  pars[['male']][['arrive_pars']] <- FitArrival(df_m)
  pars[['female']][['arrive_pars']] <- FitArrival(df_f)
  pars[['male']][['depart_pars']] <- FitDeparture(df_m)
  pars[['female']][['depart_pars']] <- FitDeparture(df_f)
  return(pars)
}

# FitStepLength Function -------------------------------------------------------

###  Fits a generalized Pareto distribution to step-length data
###  Usage: FitStepLength(data, location, scale, shape)
###  Arguments: data = dataframe with a "step_length" column
###             location = starting value of location, default is 
###               min(data$step_length)-1
###             scale = starting value of scale, default is 1
###             shape = starting value of shape, default is 0
###  Returns: list(step_location, step_scale, step_shape)
###  Notes:
###  Blake Massey
###  2014.10.03

FitStepLength <- function(df = df, 
                          location = NULL,
                          scale = 1,
                          shape = 0,
                          pars = NULL) {
  suppressPackageStartupMessages(require(bbmle))
  suppressPackageStartupMessages(require(VGAM))
  if (is.null(pars)){ 
    source('C:/Work/R/Functions/pars.R')
    pars <- CreateNewPars()
  }
  df <- df
  pars <- pars
  df$sex <- gsub("m", "male", df$sex)
  df$sex <- gsub("f", "female", df$sex)
  df_m<- subset(df, sex=="male")
  df_f<- subset(df, sex=="female") 
  if (is.null(location)) location <- min(df$step_length,  na.rm = TRUE)-1  
  step_pars_f <- mle2(NLLPareto, start=list(location=location, scale=scale,
    shape=shape), data=list(data=df_f$step_length),  method="Nelder-Mead", 
    skip.hessian=TRUE) #calculates Pareto parameters
  step_pars_m <- mle2(NLLPareto, start=list(location=location, scale=scale,
    shape=shape), data=list(data=df_m$step_length),  method="Nelder-Mead", 
    skip.hessian=TRUE) #calculates Pareto parameters
  pars$female$step = list(location=as.numeric(coef(step_pars_f)[1]), 
                           scale=as.numeric(coef(step_pars_f)[2]),   
                           shape = as.numeric(coef(step_pars_f)[3]))
  pars$male$step = list(location=as.numeric(coef(step_pars_m)[1]), 
                        scale=as.numeric(coef(step_pars_m)[2]),   
                        shape = as.numeric(coef(step_pars_m)[3]))
  return(pars)
}

# NLLBeta Function -------------------------------------------------------------

###  Beta Negative Log-Likelihood function
###  Usage: NLLBeta(data, shape1, shape2) 
###  Arguments: data = data
###             shape1 = shape1
###             shape2 = shape2 
###  Returns: NLL of Beta
###  Notes: used in BehaviorFit()
###  Blake Massey
###  2014.10.06

NLLBeta <- function(data, 
                    shape1, 
                    shape2) {  
  -sum(dbeta(data, shape=shape1, shape2=shape2, log=TRUE))
}

# NLLPareto Function -----------------------------------------------------------

###  Pareto Negative Log-Likelihood function
###  Usage: NLLPareto(data, shape, mu) 
###  Arguments: data = data
###             location = location ("left edge" of the probability density)
###             scale = scale ()
###             shape = shape
###  Returns: NLL of Pareto distribution
###  Notes: uses the generalized Pareto distribution from the VGAM package 
###  Blake Massey
###  2014.11.03

NLLPareto <- function(data, 
                      location, 
                      scale,
                      shape) {  
  suppressPackageStartupMessages(require(VGAM))
  -sum(dgpd(data, location=location, scale=scale, shape=shape, log=TRUE))
}

# NLLWeibull Function ----------------------------------------------------------

###  Weibull Negative Log-Likelihood function
###  Usage: NLLWeibull(data, shape, mu) 
###  Arguments: data = data
###             shape = shape
###             mu = mu 
###  Returns: NLL of Weibull distribution
###  Notes: used in DepartFit() and ArriveFit()
###  Blake Massey
###  2014.05.11

NLLWeibull <- function(data, 
                       shape, 
                       mu) {  
  scale <- mu/gamma(1+(1/shape))
  -sum(dweibull(data, shape=shape, scale=scale, log=TRUE))
}

# PlotBetaCDF Function ---------------------------------------------------------

###  Plot Beta Cumulative Probability Distribution Function
###  Usage: PlotBetaCDF(shape1, shape2) 
###  Arguments: shape1 = shape1, default is 5
###             shape2 = shape2, default is 5
###             colour = line color, default is "darkgreen"
###  Returns: plot of BetaCDF 
###  Notes: x-axis is set to go to 500
###  Blake Massey
###  2014.05.11

PlotBetaCDF <- function(shape1 = 5, 
                        shape2 = 5, 
                        color = "darkgreen") {
  suppressPackageStartupMessages(require(ggplot2))
  x <- seq(0, 1, length=500)
  y <- pbeta(x, shape=shape1, shape2=shape2)
  df <- as.data.frame(cbind(x, y))
  main <- paste("Beta Distribution (shape1 = ", signif(shape1,3), ", shape2 = ", 
    signif(shape2,3), ")",  sep="")
  g <- ggplot(df, aes(x, y)) +
    geom_line(color=color, size=1.5) + 
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black"))
  g+labs(x='Value', y='Cumulative Probability Density', title=main)
}

# PlotBetaPDF Function ---------------------------------------------------------

###  Plot Beta Probability Density Function
###  Usage: PlotBetaPDF(shape1, shape2)
###  Arguments: shape1 = shape1, default is 5
###             shape2 = shape2, default is 5
###             colour = line color, default is "darkgreen"
###  Returns: plot of BetaPDF
###  Notes: x-axis is set to have 500 values between 0 and 1 
###  Blake Massey
###  2014.06.13

PlotBetaPDF <- function(shape1 = 5, 
                        shape2 = 5,
                        color = "darkgreen") {
  suppressPackageStartupMessages(require(ggplot2))
  x <- seq(0, 1, length=500)
  y <- dbeta(x, shape=shape1, shape2=shape2)
  df <- as.data.frame(cbind(x, y))
  main <- paste("Beta Distribution (shape1 = ", signif(shape1,3), ", shape2 = ", 
    signif(shape2,3), ")",  sep="")
  g <- ggplot(df,aes(x,y)) +
    geom_line(color=color, size=1.5) + 
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black"))
  g+labs(x = 'Value', y= 'Probability Density', title=main)
}

# PlotCauchyPDF Function -------------------------------------------------------

###  Plot Cauchy Probability Density Function
###  Usage: PlotCauchyPDF(location, scale)
###  Arguments: location = location, default is 5
###             scale = scale, default is 1
###             xlim = sets axis limits of plot, default is: c(1st, 
###               99th quantile)
###             colour = line color, default is "darkgreen"
###  Returns: plot of CauchyPDF
###  Notes: x-axis is set to have 500 values 
###  Blake Massey
###  2014.10.13

PlotCauchyPDF <- function(location = 5, 
                          scale = 5, 
                          xlim = NULL,
                          color = "darkgreen") {
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(stats))
  if (is.null(xlim)){
    x <- seq(qcauchy(.01, location, scale), qcauchy(.99, location, scale), 
      length=500)
  } else {
    x <- seq(xlim[1], xlim[2], length=500)
  }
  y <- dcauchy(x, location=location, scale=scale)
  df <- as.data.frame(cbind(x, y))
  main <- paste("Cauchy Distribution (location = ", signif(location,3), 
    ", scale = ", signif(scale,3), ")",  sep="")
  g <- ggplot(df,aes(x,y)) +
    geom_line(color=color, size=1.5) + 
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black"))
  g + labs(x = 'Value', y= 'Probability Density', title=main) +
  geom_vline(xintercept = qcauchy(.95, location, scale), size=1, colour="red", 
    linetype = "longdash") +
  geom_vline(xintercept = qcauchy(.05, location, scale), size=1, colour="red", 
    linetype = "longdash")
}

# PlotDataAndPareto Function ---------------------------------------------------

###  Plots a histogram of the data with an option to fit a Pareto distribution  
###  Usage: PlotDataAndPareto(df, var, by, pars, xlim, bin_width, fit_pareto,
###           fit_color, hold_axes)
###  Arguments: df = dataframe of data
###             var = variable to fit Pareto distribution
###             by = optional, column name used to subset data, default = NULL 
###             pars = optional, a set of parameters used in place of pars 
###               generated within the function. Default is NULL.
###             xlim = x value limit, default is NULL
###             bin_width = bin size, default is: x-value range/30 
###             fit_pareto = logical, whether or not to fit and show Pareto 
###               distribution. Default is TRUE.
###             fit_color = color used for Pareto fit line, quantiles, and 
###               parameter value text. Default is "orangered"
###             hold_axes = logical, hold axes even if parameter distribution
###               goes off of the screen. Default is TRUE
###  Returns: a plot of the data with            
###  Notes: 
###  Blake Massey
###  2014.11.08

PlotDataAndPareto <- function(df, 
                              var,
                              by = NULL, 
                              pars = NULL, 
                              xlim = NULL, 
                              bin_width = NULL,             
                              fit_pareto = TRUE, 
                              fit_color = "orangered",
                              x_lab = NULL, 
                              hold_axes = TRUE) {  
  suppressPackageStartupMessages(require(ggplot2))
  source('C:/Work/R/Functions/gen.R')
  source('C:/Work/R/Functions/pars.R')  
  if (is.null(xlim)) xlim <- max(df[, var], na.rm=TRUE) 
  if (is.null(bin_width)) bin_width = xlim/30
  if (is.null(x_lab)) x_lab <- var
  ifelse(is.null(by), keep <- var, keep <- c(var,by))
  df <- subset(df, select = keep)
  ifelse(is.null(by),  df$by <- "all", df$by <- df[,by])
  df$var <- df[,var] 
  by_colors <- CreateColorsByAny(by=by, df=df) 
  ## hist plot
  g <- ggplot(df, aes(x=var)) +
    scale_fill_manual(values=by_colors) +
    theme(legend.position="none") +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) +
    xlab(x_lab) + ylab("Density") + 
    scale_x_continuous(limits=c(0,xlim)) +
  geom_bar(aes(y = ..density.., fill=by), color="black", binwidth=bin_width)
  if (!is.null(by)) g <- g + facet_wrap( ~ by)
  ## For "hold_axes"
  if (hold_axes == TRUE) {
  build <- ggplot_build(g)
  g <- g + coord_cartesian(ylim= c(min(build$panel$ranges[[1]]$y.range), 
                                   max(build$panel$ranges[[1]]$y.range)))  
  }
  ## Create and plot 'pars' (if fit_pareto = TRUE)
  if(fit_pareto == TRUE){
  if(is.null(pars)) pars <- FitParetoParsToData(df, var="var", by="by")
  if(is.null(by)){
    pars$by <- "all"
  } else {
    pars[,length(pars)+1] <- pars$by
    colnames(pars)[length(pars)] <- by
  } 
  dens <- FitParetoDensityToArray(df=df, var=var, pars=pars, xlim=xlim)
  if(is.null(by)){
    dens$by <- "all"
  } else {
    dens[,length(dens)+1] <- dens$by
    colnames(dens)[length(dens)] <- by
  }
  g <- g + geom_line(data=dens, aes(x=x, y=y), color=fit_color, size=1)
  }
  
  ## Set axes for the rest of the plots
  build <- ggplot_build(g)
  xmax <- max(build$panel$ranges[[1]]$x.range)
  ymax <- max(build$panel$ranges[[1]]$y.range) 
  
  ## Create and plot 'quantiles'  
  probs = c(0.5, 0.75, 0.95)
  quantiles <- ddply(df, "by", function(df) quantile(df$va,probs=probs,
    na.rm=TRUE)) 
  colnames(quantiles)[2:(length(probs)+1)] <- paste("q",probs, sep="")  
  quantiles$xmax <- xmax  # for geom_text 
  quantiles$ymax <- ymax  # for geom_text  
  g <- g + geom_text(data = quantiles, 
      aes(x = c(q0.5+(xmax*0.02), q0.75+(xmax*.02), q0.95+(xmax*.02)),
        y = c(ymax*.8, ymax*.6, ymax*.4),
        label = c(paste("Median:","\n",as.integer(q0.5)), 
        paste("75th:","\n",as.integer(q0.75)),
        paste("95th:","\n",as.integer(q0.95)))), 
      color = c("black", "gray20", "gray30"), hjust=0, vjust=1) +
    geom_vline(data = quantiles, aes(xintercept = c(q0.5, q0.75, q0.95)),
      color = c("black", "gray20", "gray30"),
      linetype = c("longdash", "dashed", "dotted"))  
    
  ## Plot 'pars' quantile lines and parameter value text
  if(fit_pareto == TRUE) {
  pars$xmax <- xmax  # for geom_text 
  pars$ymax <- ymax  # for geom_text 
  g <- g + geom_text(data = pars, aes(x = xmax*.95,y = ymax*.95,
    label = paste("Location: ", signif(location,3),"\n", 
    "Scale: ", signif(scale,3),"\n",
    "Shape: ", signif(shape,3)), sep=""), 
    color = fit_color, hjust=1, vjust=1) +
  geom_vline(data = pars, aes(xintercept = 
    c(qgpd(.5, location=location, scale=scale, shape=shape), 
      qgpd(.75, location=location, scale=scale, shape=shape),
      qgpd(.95, location=location, scale=scale, shape=shape))), 
    linetype=c("longdash", "dashed", "dotted"), colour=fit_color, size=1) +
  geom_text(data = pars,
    aes(x = c((qgpd(.5, location=location,scale=scale,shape=shape)+(xmax*0.02)),
      (qgpd(.75, location=location, scale=scale, shape=shape)+(xmax*0.02)),
      (qgpd(.95, location=location, scale=scale, shape=shape))+(xmax*0.02)), 
       y = c(ymax*.9, ymax*.7, ymax*.5),
    label = c(paste("Median:","\n",as.integer(qgpd(.5, location=location, 
      scale=scale, shape=shape))), 
    paste("75th:","\n",as.integer(qgpd(.75, location=location, 
      scale=scale, shape=shape))),
    paste("95th:","\n",as.integer(qgpd(.95, location=location, scale=scale, 
      shape=shape))))), 
  color= fit_color,hjust=0, vjust=1)
  }
  g
}

# PlotDataAndWrappedCauchy Function --------------------------------------------

###  Plots a histogram of the data with an option to fit a wrapped Cauchy 
###    distribution  
###  Usage: PlotDataAndWrappedCauchy (df, var, by, pars, bin_width, 
###    fit_wrap_cauchy, fit_color)
###  Arguments: df = dataframe of data
###             var = variable to fit wrapped Cauchy distribution
###             by = optional, column name used to subset data, default = NULL 
###             pars = optional, a set of parameters used in place of pars 
###               generated within the function. Default is NULL.
###             bin_width = bin size, default is: 15 degrees or 2*pi/24 radians 
###             fit_wrapped_cauchy = logical, whether or not to fit and show  
###               wrapped Cauchy distribution. Default is TRUE.
###             fit_color = color used for Pareto fit line, quantiles, and 
###               parameter value text. Default is "black"
###             x_lab = name for x-axis, default is 'var'
###  Returns: a plot of the data with a fitted wrapped Cauchy distribution            
###  Notes: Automatically adjusts plot for degrees or radians input, but all 
###    parameter estimates are based on radians 
###  Blake Massey
###  2014.11.14
                                                  
PlotDataAndWrappedCauchy<- function(df, 
                                    var,
                                    by = NULL, 
                                    pars = NULL, 
                                    bin_width = NULL,             
                                    fit_wrapped_cauchy = TRUE, 
                                    fit_color = "black",
                                    x_lab = NULL) {  
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(grid))
  source('C:/Work/R/Functions/baea.R')
  source('C:/Work/R/Functions/gen.R')
  source('C:/Work/R/Functions/pars.R')  
  ifelse(max(df[,var], na.rm=TRUE) <= 2*pi, radian <- TRUE, radian <- FALSE)
  if (is.null(x_lab)) x_lab <- var
  if (is.null(bin_width) & radian == TRUE) bin_width = (2*pi)/24
  if (radian == FALSE) {
    ifelse (is.null(bin_width), bin_width <- 15*(pi/180), bin_width <- 
      bin_width*(pi/180))
  }
  ifelse(is.null(by), keep <- var, keep <- c(var,by))
  df <- subset(df, select = keep)
  ifelse(is.null(by),  df$by <- "all", df$by <- df[,by])    
  ifelse(radian == TRUE, df$var <- df[,var], df$var <- df[,var]*(pi/180)) 
  by_colors <- CreateColorsByAny(by=by, df=df) 
  breaks <- seq(0, (2*pi), by=((2*pi)/12))
  breaks <- breaks[-length(breaks)]
  minor_breaks <- seq(0, 2*pi, by=bin_width)
  limits <- c(0, 2*pi)
  if (radian == TRUE){ 
    labels <- round(breaks, 2) 
  } else {
    labels <- round(breaks*(180/pi), 2)
  }
  g <- ggplot(df, aes(x=var)) +
    scale_fill_manual(values=by_colors) +
    theme(legend.position="none") +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) + 
    theme(axis.text.x = element_text(colour="grey20",size=12))+
    theme(axis.ticks = element_blank()) +
    xlab(x_lab) + ylab("Density") 
  g <- g + geom_bar(aes(y = ..density.., fill=by), color="black", 
    binwidth=bin_width) + coord_polar(start=(.5*pi)) +
    scale_x_continuous(limits=limits, breaks=breaks, minor_breaks=minor_breaks, 
      labels = labels) + scale_y_continuous(labels = NULL)
  if (!is.null(by)) g <- g + facet_wrap( ~ by)
  g
  if (fit_wrapped_cauchy == TRUE) {
    if (is.null(pars)) {
      pars <- FitWrappedCauchyParsToData(df=df, var="var", by="by")
      if(is.null(by)){
        pars$by <- "all"
      } else {
        pars[,length(pars)+1] <- pars$by
        colnames(pars)[length(pars)] <- by
      }
    }   
    dens <- FitWrappedCauchyDensityToArray(var=var, pars=pars)
    if(is.null(by)){
      dens$by <- "all"
    } else {
      dens[,length(dens)+1] <- dens$by
      colnames(dens)[length(dens)] <- by
    }
    g <- g + geom_line(data=dens, aes(x=x, y=y), color=fit_color, size=1)
    build <- ggplot_build(g)    
    pars$xmax <- max(build$panel$ranges[[1]]$theta.range)  # for geom_text 
    pars$ymax <- max(build$panel$ranges[[1]]$r.range)  # for geom_text 
    g <- g + geom_text(data = pars, aes(y = ymax, x = xmax*.875,
      label = paste("mu:", signif(mu,3),"\n","rho:", signif(rho,3)), sep=""), 
      size=4.5, color = fit_color, hjust=0, vjust=-1.5) 
  }
  g
}

# PlotDataAndWrappedNormal Function --------------------------------------------

###  Plots a histogram of the data with an option to fit a wrapped normal 
###    distribution  
###  Usage: PlotDataAndWrappedNormal (df, var, by, pars, bin_width, 
###    fit_wrap_cauchy, fit_color)
###  Arguments: df = dataframe of data
###             var = variable to fit wrapped Cauchy distribution
###             by = optional, column name used to subset data, default = NULL 
###             pars = optional, a set of parameters used in place of pars 
###               generated within the function. Default is NULL.
###             bin_width = bin size, default is: 15 degrees or 2*pi/24 radians 
###             fit_wrapped_normal = logical, whether or not to fit and show  
###               wrapped normal distribution. Default is TRUE.
###             fit_color = color used for Pareto fit line, quantiles, and 
###               parameter value text. Default is "black"
###             x_lab = name for x-axis, default is 'var'
###  Returns: a plot of the data with a fitted wrapped normal distribution            
###  Notes: Automatically adjusts plot for degrees or radians input, but all 
###    parameter estimates are based on radians 
###  Blake Massey
###  2014.11.14
                                                  
PlotDataAndWrappedNormal <- function(df, 
                                     var,
                                     by = NULL, 
                                     pars = NULL, 
                                     bin_width = NULL,             
                                     fit_wrapped_normal = TRUE, 
                                     fit_color = "black",
                                     x_lab = NULL) {  
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(grid))
  source('C:/Work/R/Functions/baea.R')
  source('C:/Work/R/Functions/gen.R')
  source('C:/Work/R/Functions/pars.R')  
  ifelse(max(df[,var], na.rm=TRUE) <= 2*pi, radian <- TRUE, radian <- FALSE)
  if (is.null(x_lab)) x_lab <- var
  if (is.null(bin_width) & radian == TRUE) bin_width = (2*pi)/24
  if (radian == FALSE) {
    ifelse (is.null(bin_width), bin_width <- 15*(pi/180), bin_width <- 
      bin_width*(pi/180))
  }
  ifelse(is.null(by), keep <- var, keep <- c(var,by))
  df <- subset(df, select = keep)
  ifelse(is.null(by),  df$by <- "all", df$by <- df[,by])    
  ifelse(radian == TRUE, df$var <- df[,var], df$var <- df[,var]*(pi/180)) 
  by_colors <- CreateColorsByAny(by=by, df=df) 
  breaks <- seq(0, (2*pi), by=((2*pi)/12))
  breaks <- breaks[-length(breaks)]
  minor_breaks <- seq(0, 2*pi, by=bin_width)
  limits <- c(0, 2*pi)
  if (radian == TRUE){ 
    labels <- round(breaks, 2) 
  } else {
    labels <- round(breaks*(180/pi), 2)
  }
  g <- ggplot(df, aes(x=var)) +
    scale_fill_manual(values=by_colors) +
    theme(legend.position="none") +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) + 
    theme(axis.text.x = element_text(colour="grey20",size=12))+
    theme(axis.ticks = element_blank()) +
    xlab(x_lab) + ylab("Density") 
  g <- g + geom_bar(aes(y = ..density.., fill=by), color="black", 
    binwidth=bin_width) + coord_polar(start=(.5*pi)) +
    scale_x_continuous(limits=limits, breaks=breaks, minor_breaks=minor_breaks, 
      labels = labels) + scale_y_continuous(labels = NULL)
  if (!is.null(by)) g <- g + facet_wrap( ~ by)
  if (fit_wrapped_normal == TRUE) {
    if (is.null(pars)) {
      pars <- FitWrappedNormalParsToData(df=df, var="var", by="by")
      if(is.null(by)){
        pars$by <- "all"
      } else {
        pars[,length(pars)+1] <- pars$by
        colnames(pars)[length(pars)] <- by
      }
    }   
    dens <- suppressWarnings(FitWrappedNormalDensityToArray(var=var, pars=pars))
    if(is.null(by)){
      dens$by <- "all"
    } else {
      dens[,length(dens)+1] <- dens$by
      colnames(dens)[length(dens)] <- by
    }
    g <- g + geom_line(data=dens, aes(x=x, y=y), color=fit_color, size=1)
    build <- ggplot_build(g)    
    pars$xmax <- max(build$panel$ranges[[1]]$theta.range)  # for geom_text 
    pars$ymax <- max(build$panel$ranges[[1]]$r.range)  # for geom_text 
    g <- g + geom_text(data = pars, aes(y = ymax, x = xmax*.875,
      label = paste("mu:", signif(mu,3),"\n","rho:", signif(rho,3)), sep=""), 
      size=4.5, color = fit_color, hjust=0, vjust=-1.5) 
  }
  g
}

# PlotParetoPDF Function -------------------------------------------------------

###  Plot Pareto Probability Density Function
###  Usage: PlotParetoPDF(location, scale, shape, xlim, color)
###  Arguments: location = location, default is 1
###             scale = scale, default is 1
###             shape = shape, default is 1
###             xlim = sets plot x lim max, default is 99th quantile
###             colour = line color, default is "darkgreen"
###  Returns: plot of Parteo PDF
###  Notes: x-axis is set to have 500 values. 5th and 95th percentiles are shown
###    with dashed red lines. 
###  Blake Massey
###  2014.10.26

PlotParetoPDF <- function(location = 1, 
                          scale = 1,
                          shape = 1,
                          xlim = NULL,
                          color = "darkgreen") {
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(stats))
  suppressPackageStartupMessages(require(VGAM))
  if (is.null(xlim)){
    x <- seq(location, qgpd(.99, location, scale, shape), length=501)
  } else {
    x <- seq(location, xlim, length=501)
  }
  x <- x[2:501]
  y <- dgpd(x, location=location, scale=scale, shape=shape)
  df <- as.data.frame(cbind(x, y))
  main <- paste("Pareto Distribution", "\n", "(location = ", signif(location,3), 
    ", scale = ", signif(scale,3),
    ", shape = ", signif(shape,3), ")",  sep="")
  g <- ggplot(df,aes(x,y)) +
    geom_line(color=color, size=1.5) + 
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) +
    scale_x_continuous (limits= c(0, max(x)))
  g + labs(x = 'Value', y= 'Probability Density', title=main) +
  geom_vline(xintercept = qgpd(.75, location=location, scale=scale, shape=shape), 
    size=1, colour="grey20", linetype = "longdash") +
  geom_vline(xintercept = qgpd(.95, location=location, scale=scale, shape=shape), 
    size=1, colour="grey30", linetype = "longdash")
}

# PlotParsBeta Function --------------------------------------------------------

###  Plot Beta probability distributions
###  Usage: PlotParsBeta(pars, length)
###  Arguments: pars = list, parameter values
###             length = number of x values used to plot lines, default is 100   
###  Returns: plot of BetaPDF with the parameter values of parameter list
###  Notes: automatically adds 1 to the length, since 0 is included
###  Blake Massey
###  2014.06.13

PlotParsBeta <- function(pars, 
                         length=100){
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(grid))
  suppressPackageStartupMessages(require(reshape2))
  source('c:/work/r/functions/baea.R')
  behavior_colors <- CreateColorsByBehavior(output=TRUE)
  male_pars <- pars[["male"]]
  female_pars <- pars[["female"]]
  length = length + 1 
  x <- seq(0, 1, length=length)
  behavior <- rep(NA, length)
  m <- rep("male",length)
  f <- rep("female", length)
  male <- data.frame(sex=m, x)
  female <- data.frame(sex=f, x)
  male$cruise <- dbeta(x, male_pars[["cruise_pars"]][["shape1"]], 
                          male_pars[["cruise_pars"]][["shape2"]])
  male$forage <- dbeta(x, male_pars[["forage_pars"]][["shape1"]], 
                          male_pars[["forage_pars"]][["shape2"]])
  male$nest <- dbeta(x, male_pars[["nest_pars"]][["shape1"]], 
                        male_pars[["nest_pars"]][["shape2"]])  
  male$loaf <- dbeta(x, male_pars[["loaf_pars"]][["shape1"]], 
                        male_pars[["loaf_pars"]][["shape2"]])  
  male$roost <- dbeta(x, male_pars[["roost_pars"]][["shape1"]], 
                         male_pars[["roost_pars"]][["shape2"]]) 
  male$territorial <- dbeta(x,male_pars[["territorial_pars"]][["shape1"]], 
                              male_pars[["territorial_pars"]][["shape2"]]) 
  female$cruise <- dbeta(x, female_pars[["cruise_pars"]][["shape1"]], 
                            female_pars[["cruise_pars"]][["shape2"]])
  female$forage <- dbeta(x, female_pars[["forage_pars"]][["shape1"]], 
                            female_pars[["forage_pars"]][["shape2"]])
  female$nest <- dbeta(x, female_pars[["nest_pars"]][["shape1"]], 
                          female_pars[["nest_pars"]][["shape2"]])  
  female$loaf <- dbeta(x, female_pars[["loaf_pars"]][["shape1"]], 
                          female_pars[["loaf_pars"]][["shape2"]])   
  female$roost <- dbeta(x, female_pars[["roost_pars"]][["shape1"]], 
                           female_pars[["roost_pars"]][["shape2"]]) 
  female$territorial <- dbeta(x,female_pars[["territorial_pars"]][["shape1"]], 
                                female_pars[["territorial_pars"]][["shape2"]]) 
  df <- as.data.frame(rbind(male, female))
  melted <- melt(df, c("sex","x"), 3:length(df), variable.name="behavior")
  ggplot(melted, aes(x = x, y=value, group= behavior, color=behavior)) +  
    facet_grid(~ sex) + theme(panel.margin=unit(1, "lines")) +
    scale_y_continuous(limits=c(0, 5)) +
    scale_color_manual(values=behavior_colors) +
    geom_line(stat="identity", size=1.5) + 
    scale_x_continuous(breaks=seq(0,1,.1)) +
    theme(plot.title=element_text(size=22)) + 
    theme(text=element_text(size=20, colour="black")) + 
    theme(axis.text=element_text(colour="black")) + labs(x="Daily Period", 
    y="Probability Density", title="Daily Behavior Probability Distributions") 
}

# PlotWeibullCDF Function ------------------------------------------------------

###  Plot Weibull Cumulative Probability Distribution Function
###  Usage: PlotWeibullCDF(data, shape, scale)
###  Arguments: shape = weibull shape, default is 1
###             scale = weibull scale, default is 1
###             max_x = maximum value on x scale, default is 250
###  Returns: plot of probability distribution
###  Notes: 
###  Blake Massey
###  2014.05.11

PlotWeibullCDF <- function(shape = 1, 
                           scale = 1, 
                           max_x = 250){
  suppressPackageStartupMessages(require(ggplot2))
  x <- seq(0, max_x, length=250)
  hx <- pweibull(x, shape=shape, scale=scale)
  df <- as.data.frame(cbind(x, hx))
  main <- paste("Weibull Distribution (shape = ", shape, ", scale = ", 
    scale, ")",  sep="")
  g <- ggplot(df, aes(x, hx)) +
    geom_line(colour="dark green", size=1.5) + 
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black"))
  g + labs(x='Value', y='Cumulative Probability Density', title=main)
}

# PlotWeibullPDF Function ------------------------------------------------------

###  Plot Weibull Probability Distribution Function
###  Usage: PlotWeibullPDF(data, shape, scale)
###  Arguments: shape = weibull shape, default is 1
###             scale = weibull scale, default is 1
###             max_x = maximum value on x scale, default is 250
###  Returns: plot of probability distribution
###  Notes:
###  Blake Massey
###  2014.05.11

PlotWeibullPDF <- function(shape = 1, 
                           scale = 1, 
                           max_x = 250){
  suppressPackageStartupMessages(require(ggplot2))
  x <- seq(0, max_X, length=250)
  hx <- dweibull(x, shape=shape, scale=scale)
  df <- as.data.frame(cbind(x, hx))
  main <- paste("Weibull Distribution (shape = ", shape, ", scale = ",
    scale, ")", sep="")
  g <- ggplot(df, aes(x, hx)) +
    geom_line(colour="dark green", size=1.5) + 
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black"))
  g + labs(x='Value', y='Probability Density', title=main)
}
