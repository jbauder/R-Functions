# --- ALL FUNCTIONS ------------------------------------------------------------
# Sources all of the R functions needed for the BAEA Project 
# ------------------------------------------------------------------------------

filenames <- c('C:/Work/R/Functions/baea.R', 
               'C:/Work/R/Functions/gen.R',
               'C:/Work/R/Functions/gis.R', 
               'C:/Work/R/Functions/gps.R',
               'C:/Work/R/Functions/kml.R',
               'C:/Work/R/Functions/pars.R', 
               'C:/Work/R/Functions/sim.R')

sapply(filenames, source)

### Function to replace text string in multiple functions
#ReplaceFilesText(filenames, "multiTrack", "Track")
