# Setup ----

## Packages ----
# TODO: Remove packages that are not actually used like units, tealeaves, optimr
# (either copy functions from plantecophys or note them in the script!)

all_pkgs <- c(
  ## Coding and Wrangling
  "here", # For organizing working directory
  "tidyverse", # Cleaner coding
  # "magrittr", # Cleaner coding
  "lubridate", # Better date handling
  "renv",
  
  ## Visuals
  "scico",
  "visdat",    # Visualizing
  # "ggridges",  # Visualizing data
  "patchwork", # Visualizing data
  "ggpmisc",   # Adding linear regression info to ggplots
  # "ggrepel",   # Avoid overlapping label
  "rgdal",
  "raster",
  
  # "ggmap",
  # "spatialEco",
  
  ## Statistics
  "lme4", # Fitting non-linear model to data
  # "zoo",  # To calculate rolling averages
  
  ## Misc
  "plantecophys"
  # "rsofun", 
  # "ingestr", 
  # "rpmodel"
  )

for (i in all_pkgs){
  
  if (i %in% c("rsofun", "ingestr")) next
  
  if(! i %in% installed.packages()){
    message("[>] Installing ", i, "-------------------------------------------")
    install.packages(i, dependencies = TRUE)
  }
  library(i, character.only = T)
}

## Scripts

load_rfiles_in_dir <- function(dir_) {
  
  for (s in list.files(dir_)) {
    
    fl <- paste0(dir_, s)
    # print(fl)
    source(fl)
  }
}

load_rfiles_in_dir(here::here("R/functions/"))
