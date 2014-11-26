# --- GENERAL FUNCTIONS --------------------------------------------------------
# General descriptive, tranformation, and plotting functions
# ------------------------------------------------------------------------------

# AddNormDens Function ---------------------------------------------------------

###  Returns a vector of normal distribution density values based on the input 
###    values and input data's estimated parameters.
###  Usage: AddNormDens(x)
###  Arguments: x = input vector
###  Returns: a vector
###  Notes: Designed for use in plyr transform functions
###  Blake Massey
###  2014.10.20

AddNormDens <- function(x) {
   m <- mean(x)
   s <- sd(x)
   dnorm(x, m, s)
  }

# ExtractRandomSample Function -------------------------------------------------

###  Takes a random sample of rows from a dataframe
###  Usage: ExtractRandomSample(df, n)
###  Arguments: df = dataframe
###             n = sample size
###  Returns: dataframe with sample         
###  Notes: returns the dataframe with the original order sequence
###  Blake Massey
###  2014.05.16

ExtractRandomSample <- function(df = df,
                                n = n) { 
  df$rowID <- seq_len(nrow(df))
  df<- (df[sample(nrow(df), n),])
  df<-df[with(df, order(df$rowID)), ]
  df$rowID <- NULL 
  row.names(df)<-NULL
  return(df)
}

# PlotColorPalette Function -------------------------------------------------

###  Plots selected color palettes with indexes and names
###  Usage: PlotColorPalette(pal, select)
###  Arguments: pal = the color palette, not in quotes
###             select = layer(s) in color palette to be displayed
###  Returns: Plot of color palette(s)
###  Notes: For displaying packaged color palettes (e.g. R_pal, SAGA_pal) 
###  Blake Massey
###  2014.08.12

PlotColorPalette <- function (pal, 
                              sel = 1:length(pal)){
  suppressPackageStartupMessages(require(plotKML))
  for (i in 1:length(pal)){
    names(pal)[i] <- paste("[", i, "] ", names(pal)[i], sep = "")  # index num
  }
  if (length(sel) > 10){
    sel <- sel[1:10]
  }
  par(mfrow = c(length(sel), 1), mar = c(1.5, 0.8, 1.5, 0.5))
  for (j in sel){
    plot(y = rep(1, length(pal[[j]])), x = 1:length(pal[[j]]), 
    axes = FALSE, xlab = "", ylab = "", pch = 15, cex = 5, col = pal[[j]])
    mtext(names(pal)[j], cex = 1, line = .1, side = 3)
  }
}

# PlotColorPie Function -------------------------------------------------

###  Plots selected colors in a pie chart
###  Usage: PlotColorPie(pal, name)
###  Arguments: pal = the color palette
###             name = logical, show name assigned to each color, default is
###               TRUE.
###             radius = the pie is drawn centered in a square box whose sides 
###               range from -1 to 1. If the character strings labeling the 
###               slices are long it may be necessary to use a smaller radius.
###  Returns: Pie chart of colors
###  Notes: For displaying a color palette 
###  Blake Massey
###  2014.11.21

PlotColorPie <- function(pal, names = TRUE, radius = 100){
  plot.new()
  par(mfrow=c(1,1))  
  if (names == TRUE){
    pie(rep(1,length(pal)), labels = names(pal), col = pal, radius = radius) 
  } else { 
    pie(rep(1,length(pal)), col = pal, radius = radius) 
  }
}

# PlotHistogramByVar -----------------------------------------------------------

###  Plots a histogram of a variable by id  
###  Usage: PlotHistogramByVar(df)
###  Arguments: df = dataframe
###             var = column name with data to create histogram
###             id = column name of unique identifier, default = "id"
###  Returns: a histogram plot         
###  Notes: 
###  Blake Massey
###  2014.10.21    
  
PlotHistogramByVar <- function(df,
                               var,
                               id = "id") {
  df <- df
  sum_df <- SummarizeSE(df, var, id)
  id_colors <- CreateIDColors(output=TRUE)  
  grid <- seq(min(df[,var], na.rm = TRUE), max(df[,var], na.rm = TRUE), 
      length = 100)
  normaldens <- ddply(df, id, function(df) {data.frame(var = grid, 
    density = dnorm(grid, mean(df[,var]), sd(df[,var])))})  
  g <- ggplot(df, aes(x = var, fill=id)) + facet_wrap( ~ id)  +
    scale_fill_manual(values=id_colors) +
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) + 
    xlab("Speed") + ylab("Density") 
  g + geom_bar(aes(y = ..density.., fill=id), colour="black", binwidth = 2) +
#    geom_text(aes(label=paste("mean: ", signif(speed,3),
#      "\n","sd: ",signif(sd,3), sep="")), data=sum_df) +
    geom_line(aes(y = density), colour="black", size=1, data = normaldens)
}

# RemoveExcept Function --------------------------------------------------------

###  Removes all the objects in the environment except one
###  Usage: RemoveExcept(object)
###  Arguments: object = objects to keep
###  Returns: environment with only the kept objects         
###  Notes: 
###  Blake Massey
###  2014.05.19

RemoveExcept <- function(object = object){
  if (length(setdiff(ls(pos = .GlobalEnv), object)) > 0) {
    rm(list=setdiff(ls(pos = .GlobalEnv), object), pos = .GlobalEnv)
}
}

# ReplaceFilesText  Function -----------------------------------------------

###  Replaces a string within files
###  Usage: ReplaceFilesText(object)
###  Arguments: files = list of filenames 
###             text = string to search for
###             replace = replacement string  
###  Returns: environment with only the one object         
###  Notes: 
###  Blake Massey
###  2014.05.19

ReplaceFilesText <- function(files, 
                             text, 
                             replace) {
  for(i in files){
    x <- readLines(i)
    y <- gsub(text, replace, x)
    cat(y, file=i, sep="\n")
  }
}

# SavePlot Function ------------------------------------------------------------

###  A wrapper function for ggsave()
###  Usage: SavePlot(filename, path)
###  Arguments: filename = file name/filename of plot
###             path  = path to save plot to (if you just want to set path and 
###               not filename). If not set, reverts to working directory
###             width = width, default is 10
###             height = height, default is 7.5
###             units	= units for width and height when either one is explicitly
###               specified (in, cm, or mm), default is "in"
###             dpi = dpi to use for raster graphics, default = 300
###  Returns: Saves a jpeg file of the last displayed plot
###  Notes: Default output format is set for PowerPoint presentations
###  Blake Massey
###  2014.10.11

SavePlot <- function (filename, 
                      path = getwd(),
                      width = 10, 
                      height = 7.5,
                      units = "in",
                      dpi=300){
  suppressPackageStartupMessages(require(ggplot2))
  if(!exists("path")){ 
    path <- getwd()
  }
  ggsave(filename = filename, path = path, width=width, height=height, 
         units=units, dpi=dpi)
}  

# SummarizeSE Function ---------------------------------------------------------

###  Summarizes data with count, mean, standard deviation, standard error of the
###    mean, and confidence interval
###  Usage: SummarizeSE (df, var, groups, na_rm, conf_int)
###  Arguments: data = a data frame.
###             var = name of column that contains the variable to be summarized
###             groups = a vector containing names of columns that contain 
###               grouping variables
###             na_rm = logical that indicates whether to ignore NA's
###             conf_interval = the percent range of the confidence interval, 
###               default is 95%)
###  Returns: a dataframe
###  Notes: Based on example function from "Cookbook for R" website.
###  Blake Massey
###  2014.10.20


SummarizeSE <- function(data, 
                        var, 
                        by = NULL, 
                        na_rm = FALSE, 
                        conf_interval = .95) {
  suppressPackageStartupMessages(require(plyr))
  LengthNaRm <- function (x, na.rm=FALSE) {
  if (na_rm == TRUE){ 
    sum(!is.na(x))
  } else {       
    length(x)
  }
  }
  ifelse(is.null(by), data$by <- "all", data$by <- data[,by]) 
  data2 <- ddply(data, by, .fun = function(x, var) {c(n = LengthNaRm(x[[var]], 
    na.rm=na_rm), mean = mean(x[[var]], na.rm=na_rm), variance = var(x[[var]], 
    na.rm=na_rm), sd = sd(x[[var]], na.rm=na_rm))}, var)
  data2 <- rename(data2, c("mean" = var))
  data2$se <- data2$sd / sqrt(data2$n)  # Calculate standard error of the mean
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ci <- qt(conf_interval/2 + .5, data2$n-1)
  data2$ci <- data2$se * ci
  return(data2)
}
