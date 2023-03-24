# Euler Code to Extract Forcing ----
## Disclaimer ----
#' To run this file, one needs to have access to worldclim, watch-wfde5, and etopo datasets.
#' Additionally, the {ingestr} packages is required to be installed and working.
#' To run this file, cd to location of this file on the cluster and use bash code below:
# bsub -u schneipa -J get_siteinfo -R "rusage[mem=25000]" "Rscript --vanilla rscript_ingest_siteinfo.R"

## Packages
library(tidyverse)
library(ingestr)

#__________________________________________________________________________________________________#
## Settings ----
## Set prefix to generalize code [!] <<<
prefix     <- "k19"
debug      <- F
start_at_i <- 14

#__________________________________________________________________________________________________#
## Make setup ----
## Load siteinfo data
siteinfo_all <- read_rds(paste0("~/data-local/", prefix, "_sitename_siteinfo_sitedata.rds")) %>%
    dplyr::select(sitename, siteinfo) %>% 
    unnest(siteinfo)

## For test run
if (debug) siteinfo_all <- siteinfo_all %>% slice(1:2)

#__________________________________________________________________________________________________#
## Start loop ----
for (i in 1:nrow(siteinfo_all)){
    
    ## Fix to start loop from where it was broken
    if (i < start_at_i) {
        next
    }

    ## Verbose
    siteinfo <- siteinfo_all[i, ]
    cat("\n > Site:", siteinfo$sitename, " | Nr.", i, "/", nrow(siteinfo_all), "\n")
    

    ## Debugging due to limited data availability:
    climate_source <- "wfde5"
    timescale     <- "h"
    
    yrs <- c(siteinfo[["year_start"]], siteinfo[["year_end"]])
    
    if (climate_source == "wfde5" & !all(yrs >= 2017)) {
        ## WFDE5 data is only available for ppfd at the moment
        getvars <- c("temp", "ppfd", "vpd", "wind")
    } else {
        getvars <- c("temp", "ppfd", "vpd", "prec", "wind")
    }
        
    ## Get elevation data
    df_etopo <- ingest(
        siteinfo,
        source = "etopo1",
        dir = "~/data/etopo/"
    ) %>% unnest()
    
    siteinfo <- left_join(siteinfo, df_etopo)
    
    ## Get meteorological data
    # Get meteo data
    df_meteo <- ingest(
        siteinfo = siteinfo,
        source    = climate_source,
        getvars   = getvars,
        dir       = "~/data/wfde5/",
        settings  = list(correct_bias = "worldclim", dir_bias = "~/data/worldclim"),
        timescale = timescale) %>%
        nest(data = names(.)[-1]) %>%
        unnest(data) %>%
        unnest(data)
    
    ## Get fAPAR
    ## Get uniformal 1.0 fapar dataframe
    df_fapar <- ingest(
        siteinfo  = siteinfo,
        source    = "fapar_unity",
        timescale = timescale) %>%
        unnest(data)
    
    ## Using remotely-sensed data: START ____
    # # Get settings
    # settings_modis <- get_settings_modis(
    #     bundle            = "modis_fpar",
    #     data_path         = "~/data/modis_subsets/",
    #     method_interpol   = "loess",
    #     network           = c("fluxnet","icos"),
    #     keep              = TRUE,
    #     overwrite_raw     = FALSE,
    #     overwrite_interpol= TRUE,
    #     n_focal           = 0
    # )
    # 
    # # Get data
    # df_modis_fpar <- ingest(
    #     fluxnet_sites,
    #     source = "modis",
    #     settings = settings_modis,
    #     parallel = FALSE,
    #     ncores = 1
    # )
    # 
    # # Rename data
    # df_modis_fpar <- df_modis_fpar %>%
    #     mutate(
    #         data = purrr::map(data, ~rename(., fapar = modisvar_filled))
    #     )
    # 
    ## Using remotely-sensed data: END ____
    
    # Merging
    df_forcing <-
        df_meteo %>%
        left_join(df_fapar) %>%
        nest(forcing = !all_of(c("sitename")))
    
    ## Get home temperature data
    df_tc_home <- ingest(
        siteinfo,
        source    = "worldclim",
        settings  = list(varnam = c("tmax")),
        dir       = "~/data/worldclim") %>%
        
        # Select mean maximum temperature of warmest month
        unnest(cols = c(data)) %>%
        rowwise() %>% 
        mutate(tc_home = max(c_across(tmax_01:tmax_12))) %>% 
        dplyr::select(!(starts_with("tmax")))
    
    ## Join dfs
    siteinfo <- 
        siteinfo %>% 
        left_join(df_tc_home) %>% 
        nest(siteinfo = !all_of(c("sitename"))) %>% 
        left_join(df_forcing)
    
    ## Save df
    saveRDS(siteinfo, paste0("~/data-local/tmp/", prefix, "_df_forc_", i, "-", nrow(siteinfo_all),".rds"))
}

#__________________________________________________________________________________________________#
## Save final df ----

## Run loop that binds all files together
df_out <- tibble()

for (i in 1:nrow(siteinfo_all)){
    ## Get filename
    fln <- paste0("~/data-local/tmp/", prefix, "_df_forc_", i, "-", nrow(siteinfo_all),".rds")
    
    ## Check if file exists
    if (!file.exists(fln)) stop("File does not exist: ", fln)
    
    ## Get file
    df_tmp <- readRDS(fln)
    
    ## Attach file
    df_out <- bind_rows(df_out, df_tmp)
    
    ## Save df_out
    saveRDS(df_out, paste0("~/data-local/", prefix, "_df_forc_from_cluster.rds"))
    
    ## Delete file
    file.remove(fln)
}
