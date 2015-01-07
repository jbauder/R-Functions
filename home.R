# --- HOMERANGE FUNCTIONS ------------------------------------------------------
# General functions for creation and analysis of home ranges
# ------------------------------------------------------------------------------

  ConvertAngle <- function(x) {
    if (x > 2*pi) x <- x-(2*pi)
    if (x < 0) x <- x+(2*pi)
    return(x)
  }
  

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
                             cellsize){ 
    x <- xmin + ((floor(x/cellsize))*cellsize) + (cellsize/2)
    y <- ymin + ((floor(x/cellsize))*cellsize) + (cellsize/2)
    output <- c(x, y)
    return(output)
  }

# CreateRedistKernel Function -------------------------------------------------

###  Create a redistribution kernel raster based on a wrapped Cauchy 
###    distribution for direction and a Pareto distribution for distance.
###  Usage: CreateRedistKernel(max_r, cellsize, mu, rho, shape, scale)
###  Arguments: max_r = maximum radius of kernel in meters, default = 300
###             cellsize = cell size in meters, default = 30
###             mu = mu parameter of wrapped Cauchy distribution, 0 radians is 
###               due east because everything is based on the Unit Circle
###             rho = rho parameter of wrapped Cauchy distribution
###             shape = shape parameter of Pareto distribution
###             scale = scale parameter of Pareto distribution
###  Returns: raster object  
###  Notes:
###  Javan Bauder and Blake Massey
###  2015.01.03

CreateRedistKernel <- function(max_r = 300,
                               cellsize = 30,
                               mu,
                               rho,
                               shape,
                               scale) {
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
          rho=rho)), 3)
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
  redist_kernel <- gpd_kernel*wrpc_kernel
  redist_kernel <- redist_kernel/sum(redist_kernel)
  return(redist_kernel)
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

# PlotSimHomeRange Function ----------------------------------------------------

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

PlotSimHomeRange <- function (df,
                              x_lim = NULL,
                              y_lim = NULL){
  by_colors = c("nest"="green2", "con"="purple") 
  g <- ggplot(df) + xlab("X") + ylab("Y") +  
  geom_text(aes(x=x,y=y,label=signif(beta,3)), vjust=1.5, hjust=-.5, size=3) +
  geom_point(data=df[1,], aes(x=nest_x, y=nest_y), size=8, color="black", 
    fill="yellow",shape=22) + 
  geom_point(data=df[1,], aes(x=con_x, y=con_y), size=8, color="black", 
    fill="purple",shape=22) + 
  geom_point(data=df[nrow(df),], aes(x=x, y=y),size=8, color="black", 
    fill="red", shape=23) +     
  geom_path(data=df, aes(x=x,y=y), arrow = arrow()) +
  geom_point(data=df, aes(x=x,y=y, fill=factor(beta_name)), size=4, 
    color="black", shape=21) +  
  scale_fill_manual(values = by_colors) +
  coord_fixed(ratio = 1) +
  theme(legend.position="none") + theme(text=element_text(size=20, 
    colour="black")) + theme(axis.text=element_text(colour="black")) + 
  theme(axis.title.x = element_text(angle = 0, vjust = 0, hjust=0.5)) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) 
  if (!is.null(x_lim)) {
    g <- g + scale_x_continuous(limits=c(0, x_lim))
  }
  if (!is.null(y_lim)) { 
    g <- g + scale_y_continuous(limits=c(0, y_lim))
  }
  return(g)
}

# PlotSimHomeRangeRaster Function ----------------------------------------------

###  Plot locations, pathways, and bias values of SimulateHomeRangeRaster() 
###    output dataframes
###  Usage: PlotSimHomeRangeRaster(df, label, x_lim, y_lim)
###  Arguments: df = dataframe with x, y, nest_x, and nest_y columns
###             label = column name for point labels
###             x_lim = maximum value of x axis, default is NULL
###             y_lim = maximum value on y axis, default is NULL
###  Returns: plot of locations and pathways
###  Notes:
###  Blake Massey
###  2015.01.06

PlotSimHomeRangeRaster <- function (df,
                                    label = "home_dist",
                                    x_lim = NULL,
                                    y_lim = NULL){
  require(grid)
  df <- df
  df$label <- df[,label]
  by_colors = c("home"="green2") 
  g <- ggplot(df) + xlab("X") + ylab("Y") +  
  geom_text(aes(x=x, y=y, label=signif(label,3)), vjust=1.5, hjust=-.5, 
            size=3) +
  geom_point(data=df[1,], aes(x=home_x, y=home_y), size=8, color="black", 
    fill="yellow",shape=22) + 
  geom_point(data=df[1,], aes(x=x, y=y),size=8, color="black", 
    fill="blue", shape=23) + 
  geom_point(data=df[nrow(df),], aes(x=x, y=y),size=8, color="black", 
    fill="red", shape=23) +     
  geom_path(data=df, aes(x=x,y=y), arrow = arrow()) +
  geom_point(data=df, aes(x=x,y=y  #, fill=factor(df[,label])
                          ), 
    size=4, color="black", shape=21) +  
  scale_fill_manual(values = by_colors) +
  coord_fixed(ratio = 1) +
  theme(legend.position="none") + theme(text=element_text(size=20, 
    colour="black")) + theme(axis.text=element_text(colour="black")) + 
  theme(axis.title.x = element_text(angle = 0, vjust = 0, hjust=0.5)) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=0.5)) 
  if (!is.null(x_lim)) {
    g <- g + scale_x_continuous(limits=c(0, x_lim))
  }
  if (!is.null(y_lim)) { 
    g <- g + scale_y_continuous(limits=c(0, y_lim))
  }
  return(g)
}
PlotSimHomeRangeRaster(df, "home_dist")

# PlotTwoLogisticCDF Function --------------------------------------------------

###  Plot Logistic Cumulative Probability Distribution Function
###  Usage: PlotTwoLogisticCDF(location_1, scale_1, location_2, scale_2, 
###    x_label,x_max)
###  Arguments: location_1 = location parameter for distribution 1
###             scale_1 = scale parameter for distribution 1
###             location_2 = location parameter for distribution 2
###             scale_2 = scale parameter for distribution 2
###             x_label = label for x-axis, default if "X Value"
###             x_max = maximum value on x scale, default is 1000
###  Returns: plot of probability distribution
###  Notes:
###  Blake Massey
###  2014.12.16

PlotTwoLogisticCDF <- function(location_1, 
                               scale_1,
                               location_2,
                               scale_2,
                               x_label = "X Value",
                               x_max = 1000){
  suppressPackageStartupMessages(require(ggplot2))
  x <- seq(0, x_max, length=1000)
  y <- exp(location_1 + scale_1*x)/(1 + exp(location_1 + scale_1*x))
  y2 <- exp(location_2 + scale_2*x)/(1 + exp(location_2 + scale_2*x))
  df <- as.data.frame(cbind(x, y, y2))
  df$x_max <- x_max
  g <- ggplot(df) +
    geom_line(aes(x, y), colour="blue3", size=1.5) +     
    geom_line(aes(x, y2), colour="red3", size=1.5) + 
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) +
    labs(x=x_label, y='Probability Density') +
    annotate("text", x=x_max*.05, y=.5, label=paste0("Location: ",
      location_1,"\n", "Scale: ", scale_1), size=6, color="blue3", hjust=0, 
      vjust=1) +   
    annotate("text", x=x_max*.95, y=.5, label=paste0("Location: ",
      location_2,"\n", "Scale: ", scale_2), size=6, color="red3", hjust=1, 
      vjust=1)  
  g
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
                                 nest_prob = 0){ 
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
                                      nest_chance = 0){ 
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
