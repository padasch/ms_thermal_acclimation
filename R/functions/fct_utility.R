# Recycled functions for analysis workflow ----

# . ----
# To wrangle data ----

## Dampen Climate Data ----
dampen_vec <- function( vec, tau ){
  
  if (length(vec)<365){
    
    rlang::abort("dampen_vec(): Aborting. Length of argument 'vec' must be at least 365.")
    
  } else {
    
    ## Add one year's data to the head of 'vec' to avoid "boundary effects"
    vec <- c(vec[1:365], vec)
    
    ## apply dampening
    vec_damped <- rep(NA, length(vec))
    vec_damped[1] <- vec[1] # this is the unwanted "boundary effect"
    for (idx in 2:length(vec)){
      dvar <- (1.0/tau) * (vec[idx] - vec_damped[idx - 1])
      vec_damped[idx] <- vec_damped[idx - 1] + dvar
    }
    
    ## remove the first year again, that was added artificially before
    vec_damped <- vec_damped[366:length(vec_damped)]
    
  }
  
  return( vec_damped )
}

## .................................................................................................
## Get growth conditions ----
change_time_res <- function(df_in,                    # input df_drivers()
                            vars    = "all",          # vector of variables to transform, or "all" numerics
                            res_in  = "h",            # input resolution of df
                            res_out = "d",            # output resolution of df
                            time_subset = 1:24,       # subset of hours to consider
                            method,            # Method for calculation
                            verbose = F
){ 
  
  ## Description:
  ## This function down-scales hourly input to a daily level, depending on chosen method.
  
  ## Debugging:
  # df_in <- df$forcing[[1]]
  # vars  <- df_in %>% select_if(~is.numeric(.x)) %>% names()
  
  ## Improvements:
  # - Could add sanity checks so that nrow(df_input)/24 = nrow(df_out)
  # - Could improve the max-function to return the mean of the highest four values of each variable
  
  ## Arguments Checks:
  if (!(res_in %in% c("hh", "h"))) {rlang::abort("change_time_res(): Only hourly or half_hourly input possible for now.")}
  if (res_out != "d") {rlang::abort("change_time_res(): Only daily output possible for now.")}
  if (vars == "all")  {vars <- df_in %>% select_if(~is.numeric(.x)) %>% names()}
  
  ## Timestep Checks: Assuming that second timestep in df_in is indicative for its resolution
  if(str_detect(df_in$date[2], "00:30:00")){midday_n <- 6} # half-hourly (half past midnight)
  if(str_detect(df_in$date[2], "01:00:00")){midday_n <- 3} # hourly      (1 a.m.)
  
  ## Method A: Daily Mean Values ----
  if (method == "A"){
    if(verbose) message("change_time_res(): Returning daily mean values without nighttime values for ppfd.")
    
    ## Separate ppfd from other variables to avoid dilution of mean by night values
    df_ppfd <-
      df_in %>% 
      dplyr::select(date, ppfd) %>% 
      dplyr::mutate(ppfd = ifelse(ppfd == 0, NA, ppfd), # Exclude night-time values from mean
                    date = lubridate::date(date)) %>%   # Make sure date is in date-format
      pivot_longer(cols = !all_of("date"), names_to = "nam", values_to = "val") %>% 
      group_by(date, nam) %>%
      summarise(val_mean = mean(val, na.rm = T), .groups = "keep") %>% 
      dplyr::mutate(val_mean = ifelse(is.nan(val_mean), 0, val_mean)) # Change all-night days back to 0 ppfd
    
    ## Get daily mean of remaining variables
    df_rem <-
      df_in %>% 
      dplyr::select(!all_of("ppfd")) %>% 
      mutate(date = lubridate::date(date)) %>% # Make sure date is in date-format
      pivot_longer(cols = !all_of("date"), names_to = "nam", values_to = "val") %>% 
      group_by(date, nam) %>%
      summarise(val_mean = mean(val, na.rm = T), .groups = "keep")
    
    ## Bind df together, ungroup, and widen to return final df
    df_out <-
      bind_rows(df_ppfd, df_rem) %>% 
      ungroup() %>% 
      pivot_wider(names_from = "nam", values_from = "val_mean") 
    
    return(df_out)
  }
  
  ## Method A_night: Daily Mean including nighttime values ----
  if (method == "A_night"){
    if(verbose) message("change_time_res(): Returning daily mean values without excluding night-time values for ppfd.")
    
    df_out <-
      df_in %>% 
      mutate(date = lubridate::date(date)) %>% # Make sure date is in date-format
      pivot_longer(cols = !all_of("date"), names_to = "nam", values_to = "val") %>% 
      group_by(date, nam) %>%
      summarise(val_mean = mean(val, na.rm = T), .groups = "keep") %>% 
      ungroup() %>% 
      pivot_wider(names_from = "nam", values_from = "val_mean") 
    
    return(df_out)
  }
  
  ## Method B: Midday values ----
  if (method == "B"){
    if(verbose) message("change_time_res(): Returning midday values for all variables.")
    
    df_out <-
      df_in %>%
      mutate(date = lubridate::date(date)) %>% 
      nest(data = !any_of("date")) %>%
      dplyr::select(date, data) %>% 
      mutate(data = purrr::map(data, ~slice_max(., ppfd, n = midday_n))) %>% 
      unnest(data) %>% 
      pivot_longer(cols = !any_of("date"), names_to = "nam", values_to = "val") %>% 
      group_by(date, nam) %>% 
      summarise(val_mean = mean(val, na.rm = T), .groups = "keep") %>% 
      ungroup() %>% 
      pivot_wider(names_from = "nam", values_from = "val_mean") 
    
    return(df_out)
  }
  
  ## Method C: Midday light and daily mean others ----
  if (method == "C"){
    if(verbose) message("change_time_res(): Returning midday for light, daily mean for rest.")
    
    ## Get midday light values first
    df_ppfd <-
      df_in %>% 
      dplyr::select(date, ppfd) %>%  
      mutate(date = lubridate::date(date)) %>% 
      nest(data = !any_of("date")) %>%
      mutate(data = purrr::map(data, ~slice_max(., ppfd, n = midday_n))) %>% 
      unnest(data) %>% 
      pivot_longer(cols = !any_of("date"), names_to = "nam", values_to = "val") %>% 
      group_by(date, nam) %>% 
      summarise(val_mean = mean(val, na.rm = T), .groups = "keep")
    
    ## Get daily mean values for remaining variables
    df_rem <-
      df_in %>% 
      mutate(date = lubridate::date(date)) %>% 
      pivot_longer(cols = !any_of("date"), names_to = "nam", values_to = "val") %>% 
      dplyr::filter(nam != "ppfd") %>% 
      group_by(date, nam) %>% 
      summarise(val_mean = mean(val, na.rm = T), .groups = "keep")
    
    df_out <-
      bind_rows(df_ppfd, df_rem) %>% 
      ungroup() %>% 
      pivot_wider(names_from = "nam", values_from = "val_mean") 
    
    return(df_out)
  }
  
  ## Method D: Individual daily max values ----
  if (method == "D"){
    if(verbose) message("change_time_res(): Returning daily max values.")
    
    df_out <-
      df_in %>% 
      mutate(date = lubridate::date(date)) %>% # Make sure date is in date-format
      pivot_longer(cols = !all_of("date"), names_to = "nam", values_to = "val") %>% 
      group_by(date, nam) %>%
      summarise(val_mean = max(val, na.rm = T), .groups = "keep") %>% 
      ungroup() %>% 
      pivot_wider(names_from = "nam", values_from = "val_mean") 
    
    return(df_out)
  }
}

# . ----
# To run model----
get_settings <- function(){
  settings <- list(
    
    ## Global options
    verbose = T,
    load_forcing = T,
    save_plots = T,
    dir_prefix = "base_plots",
    
    ## Daily forcing conditions
    daily_conditions  = "B",              # Options: A, B, C, D
    tau               = 15,               # Dampening timescale in days
    
    ## Acclimated P-Model setup
    # Methods
    method_optim      = "analytical",     # Options: analytical or numerical
    method_jmaxlim    = "smith37",        # Options: smith37 or farquhar89
    method_ftemp      = "kumarathunge19", # Options: kattge07 or kumarathunge2019
    method_eb         = "off",            # Options: off, plantecophys
    do_ftemp_kphio    = T,                # T or F to use thermally acclimated kphio
    do_ftemp_theta    = F,                # T or F to use Rogers' acclimation for theta
    num_jmax_via_ratio= T,                # Using ratio of Jmax:Vcmax to predict numerical Jmax
    
    # Parameters
    kphio_calib       = 0.04479,       # rsofun v3.3 GPP calibrated: 0.06946276, T_opt calibrated: 0.04479
    beta              = 146,
    theta             = 0.85,
    
    # Soil Stress
    do_soilmstress    = F,
    soilm             = 1,
    soilmstress       = 1,
    meanalpha         = 1,
    apar_soilm_calib  = 0.33349283,       # Options: numeric, calibrated
    bpar_soilm_calib  = 1.45602286,       # Options: numeric, calibrated
    
    # PTF
    c4                = FALSE,
    
    ## Instantaneous P-Model setup
    method_rd25       = "atkin15",        # Options: atkin15 or kumarathunge19
    method_rd         = "heskel16",       # Options: q10, arrhenius or heskel2016
    method_vcmax25    = "rpmodel_accl",   # Options: prescribed or rpmodel_accl
    method_ci         = "prentice14",     # Options: prescribed numeric or prentice14
    q10               = 2,                # Options: numeric
    
    ## Simulation P-Model setup
    rpmodel_exp = "metainfo",             # Options: "ambient", metainfo" or numeric in mol/m2/s
    vpd_inst_method = "climate_scaled",    # Combination of ("meas" or "climate") and ("const", "scaled", "interp")
    
    ## Evaluation P-Model 
    method_chi = "X21",
    which_leaf_width = 0.1      # "empirical" or numeric in [meter]
  )
  
  return(settings)
}

k19_create_df_forcing <- function(settings) {
  
  ## Check if down-scaled and dampened forcing already exists
  dir_tmp <- here("data/tmp")
  if (!dir.exists(dir_tmp)) dir.create(dir_tmp, recursive = T, showWarnings = F)
  
  fln <- paste0(here("data", "tmp"), "/k19_df_forc_", settings$daily_conditions, settings$tau, ".rds")
  
  ## Redo forcing if this function has been updated recently
  last_update <- ymd("2022-08-09")
  days_since  <- as.integer(paste(round(difftime(today(), "2022-08-09", units = "d")))) 
  
  if (file.exists(fln) && settings$load_forcing && days_since > 0) {
    ctime <- as.Date(file.info(fln)$ctime)
    ndays <- length(seq(from = as.Date(ctime), to = as.Date(now()), by = 'day'))
    rlang::inform(paste0("* Load pre-saved forcing from ", ndays, " days ago (", ctime, ")..."))
    df_forc <- readRDS(fln)
  } else {
    ## If no prepared forcing available, load forcing from scratch
    df_forc <- read_rds(here("data", "final", "k19_sitename_siteinfo_forcing_hh.rds")) %>%
      ## Reduce hourly to daily timescale
      mutate(forcing_d = purrr::map(forcing_hh, ~ change_time_res(., method = settings$daily_conditions))) %>%
      ## Attach growth conditions per day
      mutate(forcing_growth = purrr::map(forcing_d, ~ mutate(., across(!date, ~ dampen_vec(., tau = settings$tau)))))
    
    saveRDS(df_forc, fln)
  }
  
  return(df_forc)
}

k19_create_input_for_acc <- function(settings, df_forc = NA) {
  rlang::inform("* Wrangle Forcing Data...")
  
  ## Load forcing data
  if (identical(df_forc, NA)) {
    df_forc <- k19_create_df_forcing(settings)
  }
  
  ## Get climate vars (but not date)
  climate_vars <- names(df_forc$forcing_growth[[1]])
  climate_vars <- climate_vars[!climate_vars %in% c("date")]
  
  ## Load observational data
  df_tmp <- read_rds(here("data", "final", "k19_sitename_siteinfo_sitedata.rds"))
  
  ## Change from site level to site-date level
  df_site_date <-
    
    ## Unnest df with observational data
    df_tmp %>%
    unnest(sitedata) %>%
    rename(date = agg_date) %>%
    
    ## Drop original siteinfo because additional information from cluster was added
    dplyr::select(-siteinfo) %>% 
    
    ## Attach daily-forcing (keep siteinfo)
    left_join(
      df_forc %>% 
        dplyr::select(sitename, siteinfo, forcing_d) %>%
        unnest(forcing_d)) %>% 
    nest(forcing_d = all_of(climate_vars)) %>% 
    
    ## Attach growth-forcing
    left_join(
      df_forc %>% dplyr::select(sitename, forcing_growth) %>%
        unnest(forcing_growth)) %>% 
    nest(forcing_growth = all_of(climate_vars))
  
  ## Report if all sites have a forcing:
  df_site_date <-
    df_site_date %>%
    rowwise() %>%
    dplyr::mutate(qc = !is.null(forcing_d)) %>% 
    ungroup()
  
  qc <- df_site_date %>% dplyr::filter(qc == F)
  
  if (nrow(qc) != 0) {
    warning("[!] Sites without forcing were removed: \n", paste(qc$sitename, qc$date, " "))
  }
  
  df_site_date <-
    df_site_date %>%
    ungroup() %>% 
    dplyr::filter(qc == T) %>%
    dplyr::select(-qc)
  
  ## Correcting naming of leaf_size
  if (is.numeric(settings$which_leaf_width)) {
    df_site_date <- 
      df_site_date %>%
      mutate(forcing_d = purrr::map(forcing_d, ~mutate(., leaf_width = settings$which_leaf_width )),
             forcing_growth = purrr::map(forcing_growth, ~mutate(., leaf_width = settings$which_leaf_width )))
  } else {
    df_site_date <- 
      df_site_date %>%
      mutate(forcing_d = purrr::map(forcing_d, ~rename(., leaf_width = leaf_size )),
             forcing_growth = purrr::map(forcing_growth, ~rename(., leaf_width = leaf_size )))
    
  }
  
  return(df_site_date)
}

k19_run_acc <- function(df_site_date, settings) {
  
  #_____________________________________________________________________________
  rlang::inform("* Running acclimated P-Model...")
  
  
  if ("sensitivity_place_holder" %in% names(df_site_date)) {
    settings$kphio <- NULL
    settings$kphio_calib <- NULL
    settings$beta <- NULL
    settings$c_cost <- NULL
    message("!!! ATTENTION: Using kphio from df instead from settings, for sens. analys. Needs to be revised for normal runs! ")
  }
  
  df_tmp <- 
    df_site_date %>% 
    mutate(settings = list(as_tibble(settings))) %>% 
    dplyr::select(any_of(c("sitename", "date", "settings", "siteinfo", "forcing_growth"))) %>% 
    unnest(any_of(c("settings", "siteinfo", "forcing_growth")))
  
  if (("temp" %in% names(df_tmp))) {
    df_tmp <- df_tmp %>% rename(tc_air = temp) # This check is due to unproper naming in rpmodel routine
  }
  
  df_tmp <-
    df_tmp %>%
    nest(inputs = !all_of(c("sitename", "date")))
  
  df_acc <- 
    df_tmp %>%
    rowwise() %>%
    mutate(rpm_acc = rpmodel(inputs) %>%
             list()
    ) %>%
    ungroup() %>%
    dplyr::select(sitename, date, rpm_acc) %>%
    right_join(df_site_date)
  
  return(df_acc)
}

k19_get_df_inst <- function(df_in, settings, return_early = FALSE) {
  
  ## Check Method
  met <- c("climate_const", "climate_scaled", "meas_const", "meas_scaled", "meas_interp")
  if (!(settings$vpd_inst_method %in% met)) {
    stop("Selected method for vpd_inst_method cannot be used. You used: ", settings$vps_inst_method)
  }
  
  ## Define Variables
  tc_norm <- 25 # Changing this does not affect scaling of VPD
  step_size_tcair   <- 0.5
  
  ## Create df_inst
  df_inst <-  tibble(tc_leaf  = seq(step_size_tcair, 40, step_size_tcair), 
                     tc_air   = NA)
  df_inst$sitename <- unique(df_in$sitename)
  df_inst$date <- unique(df_in$date)
  
  # get patm (not reported in measurements, taking condition of the day)
  df_inst$patm <- mean(df_in$patm, na.rm = T)
  
  ## get ppfd (check if ppfd in meas_cond is available, else take metainfo)
  if (is.character(settings$rpmodel_exp)) {
    ppfd_all_na <- all(is.na(df_in$PARi))
    if (ppfd_all_na) {
      df_inst$ppfd <- mean(df_in$ppfd_measurement)/1e6
    } else {
      df_inst$ppfd <- mean(df_in$PARi, na.rm = T)/1e6
    }
  } else {
    df_inst$ppfd <- settings$rpmodel_exp ## If ppfd_measurement is prescribed in settings take that value instead
  }
  
  ## get co2
  ## From climate data (will be overwritten below if needed)
  df_inst$co2 <- mean(df_in$co2, na.rm = T) 
  
  # ## From meas data (outcommented to isolate effect of vpd)
  # if (str_detect(settings$vpd_inst_method, "meas")) {
  #   co2s_all_na <- all(is.na(df_in$CO2S))
  #   co2r_all_na <- all(is.na(df_in$CO2R))
  #   
  #   ## Use leaf-level surface data
  #   if (!(co2s_all_na)) {
  #     df_inst$co2 <- mean(df_in$CO2S, na.rm = T)
  #     ## Use reported reference co2
  #   } else if (co2r_all_na) {
  #     df_inst$co2 <- mean(df_in$co2)
  #   }
  # } 
  
  #_______________________________________________________________________________
  ## get vpd (first calculate all variations and select vpd to return in the end)
  df_tmp <- df_inst
  df_tmp$climate_tcair <- unique(df_in$temp)
  df_tmp$climate_vpd   <- unique(df_in$vpd)
  df_tmp$climate_patm  <- unique(df_in$patm)
  
  ## Calculate vpd_inst variations
  df_tmp$climate_vpd_scaled <-
    VPDleafToAir(
      df_tmp$climate_vpd / 1000,
      df_tmp$climate_tcair,
      df_tmp$tc_leaf,
      df_tmp$climate_patm / 1000
    ) * 1000
  
  ## Check if VpdL available (if not, skip leaf-level and take climate data directly)
  df_in_without_vpdna <- df_in %>% drop_na(VpdL)
  vpd_all_na <- nrow(df_in_without_vpdna)
  if (vpd_all_na == 0) {
    settings$vpd_inst_method <- str_replace(settings$vpd_inst_method, "meas", "climate")
    settings$vpd_inst_method <- str_replace(settings$vpd_inst_method, "interp", "scaled")
  } else {
    ## Leaf-level VPD data
    df_tmp$leaf_vpd_norm       <- mean(VPDairToLeaf(df_in_without_vpdna$VpdL, df_in_without_vpdna$Tleaf, tc_norm, Pa = df_in_without_vpdna$patm/1000), na.rm = T) 
    df_tmp$leaf_vpd_scaled_fun <- VPDairToLeaf(df_tmp$leaf_vpd_norm, tc_norm, df_tmp$tc_leaf, df_tmp$climate_patm/1000) * 1000
    df_tmp$leaf_vpd_meas       <- mean(df_in_without_vpdna$VpdL*1000, na.rm = T)
    
    ## Linear interpolation by binning Tleaf and take average vpd of each bin
    bin_size <- 2.5
    df_in_without_vpdna <- 
      df_in_without_vpdna %>% 
      mutate(tc_bin_1 = round(Tleaf / bin_size) * bin_size)
    
    df_2 <-
      aggregate(df_in_without_vpdna, 
                by = list(cut(df_in_without_vpdna$Tleaf, seq(0, 40, bin_size))), 
                mean) %>%
      mutate(
        tc_bin = round(Tleaf / bin_size) * bin_size,
        vpd_bin = VpdL
      )
    
    ## If only one bin available, skip and take climate data directly
    if (length(unique(df_2$tc_bin)) == 1) {
      settings$vpd_inst_method <- str_replace(settings$vpd_inst_method, "meas", "climate")
    } else {
      # Interpolate between each bin and use lowest value for bins below available Tleaf
      lin_int <- approx(df_2$tc_bin, df_2$vpd_bin, xout=df_tmp$tc_leaf, yleft = min(df_2$vpd_bin))
      df_tmp$leaf_vpd_scaled_lin <- lin_int$y * 1000
      
      ## Overwrite leaf_vpd_scaled_lin where NA
      min_ <- df_tmp %>%  slice_min(leaf_vpd_scaled_lin) %>% pull(tc_leaf)
      max_ <- df_tmp %>%  slice_max(leaf_vpd_scaled_lin) %>% pull(tc_leaf)
      
      df_tmp <-
        df_tmp %>% 
        mutate(
          # to use lower lm: leaf_vpd_scaled_lin = ifelse(tc_leaf < min_, y_low, leaf_vpd_scaled_lin),
          # to use lower bin_vpd (done above in approx()): leaf_vpd_scaled_lin = ifelse(tc_leaf < min_, min(df_tmp$leaf_vpd_scaled_lin, na.rm = T), leaf_vpd_scaled_lin),
          leaf_vpd_scaled_lin = ifelse(tc_leaf > max_, leaf_vpd_scaled_fun, leaf_vpd_scaled_lin))
    }
  }
  
  # Keep vpd of interest:
  if (str_detect(settings$vpd_inst_method, "climate")) {
    if (str_detect(settings$vpd_inst_method, "scaled")) {
      df_inst$vpd <- 
        ifelse(
          # df_tmp$climate_vpd_scaled > df_tmp$climate_vpd, df_tmp$climate_vpd_scaled,  df_tmp$climate_vpd # If scaled larger than climate, take scaled, else take climate
          df_tmp$climate_vpd_scaled > 200, df_tmp$climate_vpd_scaled,  200 # If scaled larger than 100, take scaled, else take 100
        )                             
    } else {
      df_inst$vpd <- df_tmp$climate_vpd
    }
    
  } else if (str_detect(settings$vpd_inst_method, "meas")) {
    if (str_detect(settings$vpd_inst_method, "scaled")) {
      df_inst$vpd <- df_tmp$leaf_vpd_scaled_fun
    } else if (str_detect(settings$vpd_inst_method, "interp")) {
      df_inst$vpd <- df_tmp$leaf_vpd_scaled_lin
    } else if (str_detect(settings$vpd_inst_method, "const")) {
      df_inst$vpd <- df_tmp$leaf_vpd_meas
    }
  }
  
  #_______________________________________________________________________________
  # Short-cut to return df for plotting outside
  if (return_early) {
    return_list <- 
      list(
        df_tmp = df_tmp,
        df_inst = df_inst,
        df_in = df_in
      )
    
    return(return_list)
  }
  
  #_______________________________________________________________________________
  ## Plot it
  p_vpd <- ggplot()
  
  if (str_detect(settings$vpd_inst_method, "meas")) {
    p_vpd <- 
      p_vpd +
      # geom_line(data = df_tmp,  aes(x = tc_leaf, y = leaf_vpd_meas,       color = "const", linetype = "meas")) +
      # geom_line(data = df_tmp,  aes(x = tc_leaf, y = leaf_vpd_scaled_fun, color = "scaled", linetype = "meas")) +
      geom_line(data = df_tmp,  aes(x = tc_leaf, y = leaf_vpd_scaled_lin, color = "interpolated", linetype = "meas"))
  }
  
  p_vpd <-
    p_vpd +
    geom_line(data = df_tmp,  aes(x = tc_leaf, y = climate_vpd,         color = "const", linetype = "climate")) +
    geom_line(data = df_tmp,  aes(x = tc_leaf, y = climate_vpd_scaled,  color = "scaled", linetype = "climate")) +
    # geom_line(data = df_inst, aes(x = tc_leaf, y = vpd, color = "final_choice", linetype = "final_choice")) +
    # geom_point(data = df_inst %>% dplyr::filter(row_number() %% 2 == 0), aes(x = tc_leaf, y = vpd, fill = "final_choice"), shape = 21) +
    geom_point(data = df_in, aes(x = Tleaf,   y = VpdL*1000, fill = "original_data"), shape = 21) +
    scale_color_manual(values = c("red",    "blue", "orange", "orange", "darkgreen"),
                       breaks = c("scaled", "const", "interpolated", "original_data", "final_choice")) +
    scale_fill_manual(values = c("red",    "blue", "orange", "orange", "darkgreen"),
                      breaks = c("scaled", "const", "interpolated", "original_data", "final_choice")) +
    scale_linetype_manual(values = c("dashed", "solid", "dotted"),
                          breaks = c("climate", "meas", "final_choice")) +
    geom_hline(yintercept = 0 ,color = "grey") +
    annotate(geom = "text", x = 0, y = 6000, hjust = 0, vjust = 1, 
             label = paste("ID:", unique(df_in$sitename), 
                           "\nAgg. Date:", unique(df_in$date))) +
    # labs(subtitle = paste0("Site-date: ", unique(df_in$sitename)," ", unique(df_in$date), "\n",
    #                        "Method chosen: ", settings$vpd_inst_method)) +
    xlim(0, 40) +
    ylim(0, 6000) +
    theme_classic() +
    guides(color = "none", fill = "none", linetype = "none") +
    xlab(bquote(T[leaf] ~ "[°C]")) +
    ylab(bquote(VPD[leaf] ~ "[Pa]")) 
  
  p_co2 <-
    ggplot() +
    # geom_line(data = df_inst, aes(tc_leaf, co2,  color = "final_choice")) +
    geom_line(data = df_in, aes(Tleaf, Ci, color = "meas_cond_ci")) +
    geom_line(data = df_in, aes(Tleaf, CO2S, color = "meas_cond_co2_leaf_surf")) +
    geom_line(data = df_in, aes(Tleaf, CO2R, color = "meas_cond_co2_air")) +
    geom_line(data = df_in, aes(Tleaf, co2, color = "daily_cond")) +
    geom_line(data = df_in, aes(Tleaf, co2, color = "daily_cond")) +
    scale_color_manual(values = c("darkgreen",     "blue",          "skyblue",               "purple",             "red"),
                       breaks = c("final_choice", "meas_cond_ci", "meas_cond_co2_leaf_surf", "meas_cond_co2_air" ,"daily_cond")) +
    ylim(150, 450) +
    xlim(0, 40) +    
    ylab("co2 [ppm]") +
    xlab("tc_leaf [degC]") +
    theme(legend.position = c(0.75, 0.5))
  
  p_ppfd <-
    ggplot() +
    # geom_line(data = df_inst, aes(tc_leaf, ppfd*1e6,  color = "final_choice")) +
    geom_line(data = df_in, aes(Tleaf, PARi, color = "meas_cond")) +
    geom_line(data = df_in, aes(Tleaf, ppfd*1e6, color = "daily_cond")) +
    scale_color_manual(values = c("darkgreen", "blue", "red"),
                       breaks = c("final_choice", "meas_cond", "daily_cond")) +
    ylim(500, 2500) +
    xlim(0, 40) +
    ylab("ppfd [umol/m2/s]") +
    xlab("tc_leaf [degC]") +
    theme(legend.position = c(0.75, 0.5))
  
  #_______________________________________________________________________________
  ## Define output
  df_inst <- 
    df_inst  %>% 
    dplyr::filter(vpd > 0)%>% # remove vpd below zero because impossible
    nest(forcing_inst = !all_of(c("sitename", "date"))) # nest for nice output
  
  out <- list(df_inst = df_inst,
              p_vpd   = p_vpd,
              p_co2   = p_co2,
              p_ppfd  = p_ppfd)
  
  return(out)
  
  #_______________________________________________________________________________
  ## Old Code for remembering:
  # 2. Make linear interpolation between bins
  # lin_int <- approx(df_2$tc_bin, df_2$vpd_bin, xout=df_tmp$tc_leaf,
  #                   yleft = min(df_2$vpd_bin),
  #                   yright = max(df_2$vpd_bin))
  # 3. Attach linear interpolation to df_tmp
  # df_tmp$leaf_vpd_scaled_lin <- lin_int$y * 1000
  # 4. Extra: Take linear regression of lowest and highest three bins and use them to extrapolate
  # Recode df_tmp to *not* hold left and right constant values
  # bins <- sort(unique(df_in$tc_bin_1))
  # bins_low <- bins[c(1,2)]
  # bins_upp <- bins[c(length(bins)-2,length(bins)-1,length(bins))]
  # 
  # df_in_low <- df_in %>% dplyr::filter(tc_bin_1 %in% bins_low)
  # df_in_up  <- df_in %>% dplyr::filter(tc_bin_1 %in% bins_upp)
  # 
  # lm             <- summary(lm(VpdL ~ Tleaf, data = df_in_up))
  # q              <- lm$coefficients[[1]]
  # m              <- lm$coefficients[[2]]
  # df_tmp$y_upp  <- (m * df_tmp$tc_leaf + q) * 1000
  # 
  # 
  # lm             <- summary(lm(VpdL ~ Tleaf, data = df_in_low))
  # q              <- lm$coefficients[[1]]
  # m              <- lm$coefficients[[2]]
  # df_tmp$y_low  <- (m * df_tmp$tc_leaf + q) * 1000
}

k19_run_inst <- function(df_acc, settings, sensana = FALSE) {
  
  if (sensana){ # Skip to df_inst directly if sensana
    df_tmp_1 <- df_acc %>% mutate(settings = list(as_tibble(settings)))
    
  } else  {  
    ## This function runs a thermal response curve under light saturation to estimate T_opt
    rlang::inform("* Running instantaneous P-Model...")
    
    ## Check if necessary variables are present
    if (!all(c("sitename", "date", "forcing_d", "forcing_growth", "siteinfo", "rpm_acc") %in% names(df_acc))){
      stop("One of the following variables is missing in the input: \n > sitename, date, forcing, forcing_growth, siteinfo, rpm_acc")
    } 
    
    ## Attach settings to tibble
    if (!("settings" %in% names(df_acc))) {
      df_tmp <-
        df_acc %>% 
        mutate(settings = list(as_tibble(settings)))
    }
    
    
    #_____________________________________________________________________________
    ## Create forcing_inst dataframe for simulations
    
    df_inst_all <- tibble()
    p_vpd_all   <- list()
    p_co2_all   <- list()
    p_ppfd_all  <- list()
    
    for (i in 1:nrow(df_tmp)) {
      # message(i)
      df_i <- 
        df_tmp %>%
        slice(i) %>%
        dplyr::select(sitename, siteinfo, date, meas_cond, forcing_d) %>% 
        unnest(c(meas_cond, forcing_d, siteinfo))
      
      out <- k19_get_df_inst(df_i, settings)
      
      df_inst_all <- rbind(df_inst_all, out$df_inst)
      p_vpd_all[[i]]  <- out$p_vpd
      p_ppfd_all[[i]] <- out$p_ppfd
      p_co2_all[[i]] <- out$p_co2
      
      ## Remove only x-axis
      if (i %in% c(1, 13, 25, 37, 49)) {
        p_vpd_all[[i]] <- p_vpd_all[[i]] + 
          xlab(NULL) +
          theme(axis.text.x = element_blank())
        
        ## Remove only y-axis
      } else if (i %in% c(62, 63, 64, 65, 66, 67, 56, 57, 58, 59, 60)) {
        p_vpd_all[[i]] <- p_vpd_all[[i]] +
          ylab(NULL) +
          theme(axis.text.y = element_blank())
        
        ## Remove both axis  
      } else if (i != 61) {
        p_vpd_all[[i]] <- p_vpd_all[[i]] + 
          xlab(NULL) + 
          ylab(NULL) +
          theme(axis.text.y = element_blank(),
                axis.text.x = element_blank())
      }
    }
    
    #_____________________________________________________________________________
    ## Make big forcing plot if requested
    
    if (T & settings$save_plots) { # Changed to F to speed things up
      p_vpd <- p_vpd_all[[1]]
      p_ppfd <- p_ppfd_all[[1]]
      p_co2 <- p_co2_all[[1]]
      
      for (p in 2:length(p_vpd_all)) {
        p_vpd  <- p_vpd  + p_vpd_all[[p]]  + theme(legend.position = "none")
        p_ppfd <- p_ppfd + p_ppfd_all[[p]] + theme(legend.position = "none")
        p_co2  <- p_co2  + p_co2_all[[p]]  + theme(legend.position = "none")
      }
      
      p_vpd  <- p_vpd  + plot_layout(guides = "collect", ncol = 12)
      p_ppfd <- p_ppfd + plot_layout(guides = "collect", ncol = 12)
      p_co2  <- p_co2  + plot_layout(guides = "collect", ncol = 12)
      
      message("[>] Saving plots to: ", settings$dir_now)
      ggsave(paste0(settings$dir_now, "/", "inst_vpd_",  settings$vpd_inst_method, ".pdf"), p_vpd,  height = 15, width = 32)
      # ggsave(paste0(settings$dir_now, "/", "inst_ppfd_", settings$vpd_inst_method, ".pdf"), p_ppfd, height = 18, width = 45)
      # ggsave(paste0(settings$dir_now, "/", "inst_co2_",  settings$vpd_inst_method, ".pdf"), p_co2,  height = 18, width = 45)
    }
    
    df_tmp_1 <- left_join(df_tmp, df_inst_all)
    
  }
  
  #_____________________________________________________________________________
  ## Extract the relevant variables needed for forcing the instant pmodel
  df_tmp_2 <-
    df_tmp_1 %>%
    mutate(tc_home           = purrr::map_dbl(siteinfo,       ~pull(., tc_home)),
           ppfd_measurement  = purrr::map_dbl(siteinfo,       ~pull(., ppfd_measurement)),
           tc_growth_air  = purrr::map_dbl(forcing_growth, ~pull(., temp)),
           tc_growth_leaf = purrr::map_dbl(rpm_acc,        ~pull(., tc_leaf)),
           kphio          = purrr::map_dbl(rpm_acc,        ~pull(., kphio  )),
           vcmax25        = purrr::map_dbl(rpm_acc,        ~pull(., vcmax25)),
           rd25           = purrr::map_dbl(rpm_acc,        ~pull(., rd25)),
           jmax25         = purrr::map_dbl(rpm_acc,        ~pull(., jmax25)),
           xi             = purrr::map_dbl(rpm_acc,        ~pull(., xi))) %>%
    dplyr::select(-any_of(c("meas_cond", "fit_opt", "rpm_acc", "forcing_growth",
                            "forcing_d", "data_org", "fit_opt", "siteinfo"))) %>% 
    unnest(c("settings", "forcing_inst"))  %>% 
    mutate(id = row_number(), # Have to add a row_id for 1-by-1 feeding of inputs
           method_eb = "off") # To ensure EB is *not* called here
  
  #_____________________________________________________________________________
  ## Run instant model
  df_rpm_inst <-
    df_tmp_2 %>% 
    nest(inputs = !any_of(c("sitename", "date", "id"))) %>%
    rowwise() %>%
    mutate(rpm_inst = rpmodel_inst(inputs) %>%
             as_tibble() %>%
             list()
    ) %>%
    ungroup()
  
  if (!sensana) {
    # Fast extraction for none edge conditions
    
    df_traits_inst <-
      df_rpm_inst %>%
      mutate(
        tc_growth_air  = purrr::map_dbl(inputs, ~pull(., tc_growth_air )),
        tc_growth_leaf = purrr::map_dbl(inputs, ~pull(., tc_growth_leaf)),
        # Round to nearest 0.5 to extract fitting A_growth
        tc_growth_air = round(tc_growth_air/0.5)*0.5,
        tc_growth_leaf = round(tc_growth_leaf/0.5)*0.5,
      ) %>%
      # dplyr::select(-inputs, -id) %>%
      unnest(rpm_inst) %>%
      nest(rpm_inst = !all_of(c("sitename", "date"))) %>%
      mutate(
        tc_opt        = purrr::map_dbl(rpm_inst, ~ slice_max(., anet)   %>% pull(tc_leaf)),
        tc_opt_agross = purrr::map_dbl(rpm_inst, ~ slice_max(., agross) %>% pull(tc_leaf)),
        tc_opt_ac     = purrr::map_dbl(rpm_inst, ~ slice_max(., ac)     %>% pull(tc_leaf)),
        tc_opt_aj     = purrr::map_dbl(rpm_inst, ~ slice_max(., aj)     %>% pull(tc_leaf)),
        
        min_a         = purrr::map_chr(rpm_inst, ~ slice_max(., anet)   %>% pull(min_a)),
        min_a         = as.factor(min_a),
        
        agross_opt    = purrr::map_dbl(rpm_inst, ~ slice_max(., agross) %>% pull(agross) * 1e6),
        anet_opt      = purrr::map_dbl(rpm_inst, ~ slice_max(., anet) %>% pull(anet) * 1e6),
        
        anet_growth   = 0, # purrr::map_dbl(rpm_inst, ~ dplyr::filter(., tc_growth_air == tc_leaf) %>% pull(anet) * 1e6),
        min_agrowth   = "ac", # purrr::map_chr(rpm_inst, ~ dplyr::filter(., tc_growth_air == tc_leaf) %>% pull(min_a)),
        min_agrowth   = "ac", # as.factor(min_agrowth),
        
        tspan_l       = purrr::map_dbl(rpm_inst,
                                       ~ dplyr::filter(., anet > (0.9 * max(anet))) %>%
                                         slice_min(tc_leaf) %>%
                                         pull(tc_leaf)),
        tspan_h       = purrr::map_dbl(rpm_inst,
                                       ~ dplyr::filter(., anet > (0.9 * max(anet))) %>%
                                         slice_max(tc_leaf) %>%
                                         pull(tc_leaf)),
        tspan          = tspan_h - tspan_l
      ) %>%
      nest(rpm_sim = !all_of(c("sitename", "date", "rpm_inst"))) %>%
      right_join(df_acc)
    
  } else {
    # Slower extraction to capture errors in sensitivity analysis
    df_traits_inst <-
      df_rpm_inst %>% 
      mutate(
        tc_growth_air  = purrr::map_dbl(inputs, ~pull(., tc_growth_air )),
        tc_growth_leaf = purrr::map_dbl(inputs, ~pull(., tc_growth_leaf)),
        # Round to nearest 0.5 to extract fitting A_growth
        tc_growth_air = round(tc_growth_air/0.5)*0.5,
        tc_growth_leaf = round(tc_growth_leaf/0.5)*0.5,
      ) %>%
      # dplyr::select(-inputs, -id) %>% 
      unnest(rpm_inst) %>% 
      nest(rpm_inst = !all_of(c("sitename", "date"))) %>% 
      mutate(
        anet_qc   = purrr::map_dbl(rpm_inst, ~ slice_max(., anet) %>% nrow()),
        agross_qc = purrr::map_dbl(rpm_inst, ~ slice_max(., agross) %>% nrow()),
        ac_qc     = purrr::map_dbl(rpm_inst, ~ slice_max(., ac) %>% nrow()),
        aj_qc     = purrr::map_dbl(rpm_inst, ~ slice_max(., aj) %>% nrow()),
        tspan_l_qc  = purrr::map_dbl(rpm_inst, 
                                     ~ dplyr::filter(., anet > (0.9 * max(anet))) %>%
                                       slice_min(tc_leaf) %>% 
                                       nrow(.)),
        tspan_h_qc  = purrr::map_dbl(rpm_inst, 
                                     ~ dplyr::filter(., anet > (0.9 * max(anet))) %>% 
                                       slice_max(tc_leaf) %>% 
                                       nrow(.))
      )
    
    for (i in 1:nrow(df_traits_inst)) {
      
      df_traits_inst$tc_opt[i]        <- NA
      df_traits_inst$tc_opt_agross[i] <- NA
      df_traits_inst$tc_opt_ac[i]     <- NA
      df_traits_inst$tc_opt_aj[i]     <- NA
      df_traits_inst$min_a[i]         <- NA
      df_traits_inst$agross_opt[i]    <- NA
      df_traits_inst$anet_opt[i]      <- NA
      
      df_traits_inst$tspan_h[i] <- NA
      df_traits_inst$tspan_l[i] <- NA
      df_traits_inst$tspan[i]   <- NA
      
      # Short debug for agrowth variables
      df_traits_inst$anet_growth[i] <- 0
      df_traits_inst$min_agrowth[i] <- "ac"
      df_traits_inst$min_agrowth[i] <- "ac"
      
      # Check for anet peaked related traits
      if (df_traits_inst$anet_qc[i] != 0 |
          df_traits_inst$agross_qc[i] != 0 |
          df_traits_inst$ac_qc[i] != 0 |
          df_traits_inst$aj_qc[i] != 0) {
        
        df_traits_inst$tc_opt[i]        <- df_traits_inst$rpm_inst[[i]] %>% slice_max(anet)   %>% pull(tc_leaf)
        df_traits_inst$tc_opt_agross[i] <- df_traits_inst$rpm_inst[[i]] %>% slice_max(agross)   %>% pull(tc_leaf)
        df_traits_inst$tc_opt_ac[i]     <- df_traits_inst$rpm_inst[[i]] %>% slice_max(ac)   %>% pull(tc_leaf)
        df_traits_inst$tc_opt_aj[i]     <- df_traits_inst$rpm_inst[[i]] %>% slice_max(aj)   %>% pull(tc_leaf)
        df_traits_inst$min_a[i]         <- df_traits_inst$rpm_inst[[i]] %>% slice_max(anet)   %>% pull(min_a)
        df_traits_inst$anet_opt[i]      <- df_traits_inst$rpm_inst[[i]] %>% slice_max(anet) %>% pull(anet) * 1e6
        df_traits_inst$agross_opt[i]    <- df_traits_inst$rpm_inst[[i]] %>% slice_max(agross) %>% pull(agross) * 1e6
        
      }
      
      # Check for tspan related traits
      if (df_traits_inst$tspan_h_qc[i] != 0 |
          df_traits_inst$tspan_l_qc[i] != 0) {
        
        df_traits_inst$tspan_h[i] <- 
          df_traits_inst$rpm_inst[[i]] %>% 
          dplyr::filter(anet > (0.9 * max(anet))) %>%
          slice_max(tc_leaf) %>% 
          pull(tc_leaf)
        
        df_traits_inst$tspan_l[i] <-
          df_traits_inst$rpm_inst[[i]] %>% 
          dplyr::filter(anet > (0.9 * max(anet))) %>%
          slice_min(tc_leaf) %>% 
          pull(tc_leaf)
        
        df_traits_inst$tspan[i]   <- df_traits_inst$tspan_h[i] - df_traits_inst$tspan_l[i]
      }
    }
    
    df_traits_inst <- 
      df_traits_inst %>% 
      mutate(min_a = as.factor(min_a)) %>% 
      dplyr::select(-contains("qc")) %>% 
      nest(rpm_sim = !all_of(c("sitename", "date", "rpm_inst"))) %>%
      right_join(df_acc)
  }
  
  return(df_traits_inst)
}

k19_create_final_df <- function(df_inst, settings) {
  
  rlang::inform("* Extract evluation data...")
  df_tmp <- df_inst %>% mutate(fit_thermal_response = list(tibble()))
  
  ## Attach fitted temperature response curve 
  for (i in 1:nrow(df_inst)) {
    df_tmp$fit_thermal_response[[i]] <- 
      tibble(
        tleaf = seq(0, 50, 0.1),
        anet = df_tmp$fit_opt[[i]]$aopt - ( df_tmp$fit_opt[[i]]$b * ( tleaf - df_tmp$fit_opt[[i]]$aopt) ^ 2),
        tc_growth_air = round(df_tmp$forcing_growth[[i]]$temp/0.1)*0.1)
  }
  
  ## Extract observed T_opt, A_opt, T_growth
  df_eval <-
    df_tmp %>%
    # slice(1:3) %>%
    dplyr::select(sitename, date, fit_opt, fit_thermal_response) %>%
    mutate(
      tc_opt       = purrr::map_dbl(fit_opt, ~ pull(., topt)),
      tc_opt_se    = purrr::map_dbl(fit_opt, ~ pull(., topt.se)),
      anet_opt     = purrr::map_dbl(fit_opt, ~ pull(., aopt)),
      anet_opt_se  = purrr::map_dbl(fit_opt, ~ pull(., aopt.se)),
      agrowth      = 0, #purrr::map_dbl(fit_thermal_response, ~ dplyr::filter(., tleaf == tc_growth_air) %>% pull(anet)),
      agrowth      = 0, #ifelse(agrowth > 0, agrowth, NA),
      
      tspan_l       = purrr::map_dbl(fit_thermal_response, 
                                     ~ dplyr::filter(., anet > (0.9 * max(anet))) %>% 
                                       slice_min(tleaf) %>% 
                                       pull(tleaf)),
      tspan_h       = purrr::map_dbl(fit_thermal_response, 
                                     ~ dplyr::filter(., anet > (0.9 * max(anet))) %>% 
                                       slice_max(tleaf) %>% 
                                       pull(tleaf)),
      tspan          = tspan_h - tspan_l
      
    ) %>%
    ## Nest evaluation data into eva_data
    nest(eva_data = any_of(c("tc_opt", "tc_opt_se", "anet_opt", "anet_opt_se", "agrowth", "tspan_l", "tspan_h", "tspan"))) %>%
    dplyr::select(sitename, date, eva_data) %>%
    right_join(df_inst)
  
  df_plot <-
    df_eval %>%
    unnest(rpm_sim,  names_sep = "_") %>%
    unnest(eva_data, names_sep = "_") %>%
    unnest(siteinfo) %>%
    unnest(forcing_growth) %>%
    
    ## Calculate biases
    mutate(bias_tc_opt =   (rpm_sim_tc_opt   - eva_data_tc_opt)    / (eva_data_tc_opt)   * 100)
  
  # saveRDS(df_plot, "~/projects/mscthesis/data/kumarathunge2019_acclimation/processed/df_plot.rds")
  
  return(df_plot)
}

k19_from_forc_to_plot <- function(settings, 
                                  debug_slice = F, 
                                  returndf = "") {
  
  ## Define var needed later
  skip_to_plotting <- F
  savep <- settings$save_plots
  
  ## Change settings if debug_slice is activated
  if (debug_slice) settings$save_plots <- F
  
  ## First check if this model setup has been run before today:
  dir_path  <- here("output", "base_plots","model_runs")     
  dir_setup <- paste0(
    "setup",
    "__FORC__", settings$daily_conditions, settings$tau,
    "__PARAM__", settings$kphio_calib, "_ftempvcmaxjmax-", settings$method_ftemp,
    "__ACC__", settings$method_jmaxli, "_ftempkphio-", settings$do_ftemp_kphio,
    "__INST__", settings$vpd_inst_method, "_ppfd-", settings$rpmodel_exp)
  
  # Make directory to save model    
  dir_model <- paste0(dir_path, "/", dir_setup)
  dir.create(dir_model, recursive = T, showWarnings = F)
  
  # Make directory to save basic plots
  dir_now <- here("output", "base_plots", today(), dir_setup)
  settings$dir_now <- dir_now
  
  dir.create(settings$dir_now, recursive = T, showWarnings = F)
  
  if (!debug_slice & settings$load_forcing & returndf %in% c("", "df_plot")) {
    fln <- paste0(dir_model, "/df_plot.rds")
    if (file.exists(fln)) {
      # If folder exists, load saved model
      
      ctime <- as.Date(file.info(fln)$ctime)
      ndays <- length(seq(from = as.Date(ctime), to = as.Date(now()), by = 'day')) - 1
      warning(paste0("* Load model from ", ndays, " days ago (", ctime, ")... \n", "  Setup: ", dir_setup))
      
      df_plot  <- readRDS(paste0(dir_model, "/df_plot.rds"))
      settings <- readRDS(paste0(dir_model, "/settings.rds")) 
      settings$dir_now <- dir_now
      
      skip_to_plotting <- T
      settings$save_plots <- savep 
    }
    # Create folder if not yet existing
    dir.create(dir_model)
  }
  
  # If setup does not exist,
  # Or if setup is not from today,
  # Or if debug_slice is requested, run model and save it (only if not debug_sliced)
  
  ## Make directory for this model run and save settings:
  if (settings$save_plots & !debug_slice & !skip_to_plotting) {
    
    sink(paste0(settings$dir_now, "/", "model_setup.txt"))
    cat(
      "\n Forcing: ", settings$daily_conditions, settings$tau,
      "\n\n",
      "Acc: ", 
      "| ", settings$method_optim,
      "| ", settings$method_jmaxlim,
      "| f_kphio", settings$do_ftemp_kphio,
      "| ", settings$method_ftemp,
      "| eb is", settings$method_eb,
      "| ftemp_theta", settings$do_ftemp_theta,
      "| num_jmax_via_ratio", settings$num_jmax_via_ratio, 
      "\n\n",
      "Params: ",
      "| kphio_calib is", settings$kphio_calib,
      "\n\n",
      "Inst: ",
      settings$vpd_inst_method,
      "\n\n",
      "Missing: More info on inst simulation and other settings!")
    sink()
  }
  
  if (returndf == "update_settings") {
    return(settings)
  }
  
  if (!skip_to_plotting) {
    
    df_site_date <- k19_create_input_for_acc(settings)
    if (debug_slice) df_site_date <- df_site_date %>% slice(round(seq(1, nrow(df_site_date), length.out = 20)))
    if (returndf == "df_site_date") return(list(df = df_site_date, set = settings))
    
    df_acc       <- k19_run_acc(df_site_date, settings)
    if (returndf == "df_acc") return(list(df = df_acc, set = settings))
    
    df_inst      <- k19_run_inst(df_acc, settings)
    if (returndf == "df_inst") return(list(df = df_inst, set = settings))
    
    df_plot      <- k19_create_final_df(df_inst, settings)
    
    ## Save df_plot for later use
    if (!debug_slice) {
      saveRDS(df_plot, paste0(dir_model, "/df_plot.rds"))
      saveRDS(settings, paste0(dir_model, "/settings.rds"))
    }
  }
  
  if (returndf == "df_plot") return(list(df = df_plot, set = settings))
  
  p_all_out <- plot_all_final_plots_from_df_plot(df_plot, settings)
  
  return(list(df = df_plot, set = settings, p = p_all_out))
  
}

# Additional functions
get_df_ref <- function(settings = NA){
  
  if(is.na(settings)[[1]]) {
    settings <- get_settings()
    message("Standard settings from get_settings()")
  }
  
  df_ref <- tibble(
    ## Environment:
    tc        = 25,
    tc_home   = tc,
    tc_growth = tc,
    temp      = tc, # needed for old formulation of rpmodel where temp is used instead of tc_air/tc_leaf, etc.
    vpd       = 1000,
    co2       = 400,
    ppfd      = 1500e-6,
    patm      = 101325,
    fapar     = 1,
    kphio     = settings$rpmodel_accl$kphio_calib,
    ci        = 275,
    wind      = 2,
    
    ## Photosynthetic variables
    kmm       = calc_kmm(25, 101325),
    gammastar = calc_gammastar(25, 101325),
    ns_star   = 1,
    
    ## Parameters
    nsteps      = 50,
    vcmax_start = 50e-6,
    jmax_start  = 100e-6,
    gs_start    = 5e-6,
    
    ## Added for sensitivity analysis
    xi          = 48,
    kphio       = settings$kphio_calib,
    kphio_calib = settings$kphio_calib,
    vcmax25     = vcmax_start,
    jmax25      = jmax_start,
    gs          = gs_start,
    ppfd_measurement = 1500e-6,
    tc_leaf     = tc,
    tc_growth_air = tc,
    tc_growth_leaf = tc,
    beta = 146) 
  
  return(df_ref)
}

# . ----
# To calculate statistics ----
standard_error <- function(x) {
  return(sd(x, na.rm = T) / sqrt(length(x)))
}

# . ----
# To make plots ----

## Anet-Tleaf ----
plot_topt_extraction <- function(df_in,
                                 add_simulation = F,
                                 dir = "automatic",
                                 make_one_plot = F,
                                 called_from_wrangling = F) {
  
  # This function plots the observed A_net-T curves,
  # based on the extracted parabola fit and its metrics.
  # Additionally, it plots the P-Model simulated A_net-T curves.
  
  # Create directory if inexistent
  if (dir == "automatic") {
    # Save it
    directory <- here("output", "wrangling_experimental_data",  "figures")
  } else {
    directory <- here(dir)
  }
  if (!dir.exists(directory)) dir.create(directory, recursive = T, showWarnings = F)
  
  # Make empty list for big plot
  p_list <- list()
  
  # Arrange df_in by climate zone - sitename - date
  if ("cl" %in% names(df_in)) df_in <- df_in %>% arrange(cl, sitename, date)
  
  for (i in 1:nrow(df_in)) {
    
    # Get a slice of data
    df_i <- df_in %>% slice(i) %>% unnest(fit_opt)
    
    # Check if t_opt was extracted from the slice
    # if (identical(df_i$topt, NA)) {
    #   warning("> T_opt extraction failed for: ", df_i$sitename, " ", df_i$agg_date)
    #   next
    # }
    
    # Rename if date column has "wrong" name
    if ("date" %in% names(df_i)) {
      df_i <- df_i %>% rename(agg_date = date)
    }
    
    # Unnest data
    df_i <- df_i %>% unnest(data_org)
    
    # Get info on site and extraction
    site       <- as.character(df_i$sitename[1])
    date       <- as.character(df_i$sample_date[1])
    agg_date   <- as.character(df_i$agg_date[1])
    fit_method <- as.character(df_i$fit_method[1])
    r2         <- round(df_i$r2[1], 2)
    pval       <- round(df_i$pvalue[1], 2)
    tc_opt     <- round(df_i$topt[1], 2)
    tc_opt_se  <- round(df_i$topt.se[1], 2)
    n_points   <- nrow(df_i)
    
    df_i$x <- as.double(df_i[["Tleaf"]])
    df_i$y <- as.double(df_i[["Photo"]])
    
    # Get extracted parabola
    df_2 <- tibble(
      tleaf = seq(0, 50, 0.1),
      anet = df_i$aopt[1] - ( df_i$b[1] * ( tleaf - df_i$topt[1]) ^ 2))
    
    # Define coef or easy plotting of T_opt:
    coef_1 <- df_i[1, ]
    ymax   <- 
      ifelse(
        is.na(coef_1$aopt),
        max(df_i$y),
        max(df_i$y, coef_1$aopt)
      ) + 3
    if (T) {
      ymax <- 30
      # warning("Fixing max ylim at", ymax)
    }
    
    
    if (is.na(df_i$topt[1])) {
      # Plot if extraction failed ______________________________________________
      p <-
        df_i %>% 
        ggplot() +
        geom_point(aes(x = x, y = y, color = as.factor(sample_date))) +
        xlim(0, 50) +
        ylim(0, ymax) +
        ylab(bquote(A[net] ~ "[µmol" ~ CO[2] ~ m^-2 ~ s^-1 ~ "]")) +
        xlab(bquote(T[leaf] ~ "[°C]")) +
        scale_color_brewer(palette = "Dark2") +
        labs(subtitle = paste0("Extraction failed for \n Site: ", df_i$sitename, 
                               ", Date: ", df_i$agg_date)) + 
        guides(linetype = "none") +
        theme_classic() +
        theme(legend.position = "right",
              legend.box = "vertical",
              plot.caption = element_text(hjust = 0))
    } else {
      
      # Plot extracted parabola __________________________________________________
      add_text <- paste0(
        "Site: ", df_i$sitename, 
        "\n Date: ", df_i$agg_date,
        "\n Model: ", fit_method,
        "\n p-value < ", pval,
        "\n r^2 = ", r2,
        "\n topt = ", tc_opt, " +- ", tc_opt_se, 
        "\n n = ", n_points)
      
      
      p <-
        ggplot() +
        geom_point(data = df_i, aes(x = x, y = y, color = as.factor(sample_date))) +
        geom_line(data = df_2, aes(x = tleaf, y = anet, linetype = "Fitted")) +
        geom_errorbarh(data = coef_1, aes(y = aopt, xmin = topt - topt.se, xmax = topt + topt.se), color = "red") +
        geom_errorbar(data = coef_1, aes(x = topt, ymin = aopt - aopt.se, ymax = aopt + aopt.se), color = "red") +
        geom_point(data = coef_1, aes(x = topt, y = aopt),
                   shape = 21, 
                   size = 3,
                   fill = "red") + # shape = 18, size = 4
        xlim(0, 50) +
        ylim(0, ymax) +
        ylab(bquote(A[net] ~ "[µmol" ~ CO[2] ~ m^-2 ~ s^-1 ~ "]")) +
        xlab(bquote(T[leaf] ~ "[°C]")) +
        # annotate("text", x = 25, y = ymax, label = add_text) +
        labs(
          color = "Sampling Date: ",
          # title = "Example for date-wise aggregation of measurements",
          subtitle = add_text
        ) +
        scale_color_brewer(palette = "Dark2") +
        scale_linetype_manual(
          "", 
          breaks = c("Fitted", "Simulated"),
          values = c("solid", "dashed")) +
        guides(linetype = "none") +
        theme_classic() +
        theme(legend.position = "right",
              legend.box = "vertical",
              plot.caption = element_text(hjust = 0))
      
      
      # Add simulated response on top of parabola ________________________________
      if (add_simulation) {
        limrate  <- paste0(unique(df_i$rpm_sim_min_a))
        sim_topt <- unique(df_i$rpm_sim_tc_opt)
        sim_aopt <- unique(df_i$rpm_sim_anet_opt)
        
        
        # for (d in 1:length(unique(df_i$sample_date))) {
        #   if (d == 1) {
        #     txt_sa <- paste0(unique(df_i$sample_date[d]))
        #   } else {
        #     txt_sa <- paste0(txt_sa, ", ", unique(df_i$sample_date[d]))
        #   }
        # }
        
        txt_ann <- 
          paste0("ID: ", df_i$sitename[1], " (", df_i$cl[1], ")",
                "\nAgg. Date: ", df_i$agg_date[1],
                # ", #Sampl. Days: ", length(unique(df_i$sample_date)), # txt_sa,
                "\nMod: T_opt = ", sim_topt, ", Lim: ", limrate,
                "\nObs: T_opt = ", tc_opt, " +- ", tc_opt_se, "\nFit: p < ", pval, ", R^2 = ", r2)
                 # "\n   Fit: ", fit_method, ", p < ", pval, ", R^2 = ", r2)
        
        p <-
          p +
          geom_line(
            data = df_i %>% unnest(rpm_inst),
            aes(x = tc_leaf, y = anet * 10^6, linetype = "Simulated"),
          ) +
          geom_point(data = df_i %>% slice(1),
                     aes(x = rpm_sim_tc_opt,
                         y = rpm_sim_anet_opt),
                     fill = "orange",
                     shape = 24,
                     size = 3) +
          annotate(geom = "text", x = 1, y = 30, vjust = 1, hjust = 0, label = txt_ann) +
          guides(linetype = "legend") +
          labs(linetype = "Line: ",
               subtitle = NULL)
      }
      
    }
    
    # Save plot for big plot
    if (make_one_plot) {
      
      plot_ncol <- 12
      n_cols <- nrow(df_in) - plot_ncol
      n_rows <- ceiling(nrow(df_in) /  plot_ncol)
      
      p_list[[i]] <- p + guides(color = "none", linetype = "none")
      
      
      ## Axis labelling for final plot
      ## Keep both axes for bottom left plot
      if (i == (n_rows -1) * plot_ncol + 1) {
        # Do nothing
        
      } else {
        
        ## Remove x axis, if not on bottom row 
        if (i %in% seq(1, n_cols)) {
          p_list[[i]] <- p_list[[i]] + 
            xlab(NULL) +
            theme(axis.text.x = element_blank())
        }
        
        ## Remove y axis if not first in row
        if (i %% plot_ncol != 1) {
          p_list[[i]] <- p_list[[i]] +
            ylab(NULL) +
            theme(axis.text.y = element_blank())
        }
      }
    }
  }
  
  if (make_one_plot) {
    
    # Make one big plot
    p_all <- p_list[[1]]
    for (pi in 2:length(p_list)) {
      p_all <- p_all + p_list[[pi]]
    }
    
    p_all <- p_all + plot_layout(plot_ncol)
    
    # Save Plot
    if (called_from_wrangling) {
      ggsave(paste0(directory, "/all_sites_in_one.pdf"), p_all, width = 32, height = 15, units = "in")
      message("Big plot was saved under: ", directory)
    }
    
    return(list(p_big_plot = p_all))
  }
}

plot_shape_thermal_response <- function(df_plot){
  
  # Original Code Below, working as it is
  ## Reduce dataframe to relevant variables
  df_1 <- 
    df_plot %>% 
    dplyr::select(sitename, date, rpm_inst, fit_opt, data_org, cl, rpm_sim_tc_opt) %>% 
    mutate(fit_inst = purrr::map(fit_opt, ~tibble(rep())))
  
  ## Loop through each simulated response and extract fit
  df_out <- tibble()
  
  for (i in 1:nrow(df_1)) {
    df_i <- df_1 %>% slice(i) %>% unnest(fit_opt)
    
    df_2 <- tibble(
      tc_leaf    = seq(0.5, 40, 0.5),
      anet_fit   = df_i$aopt[1] - ( df_i$b[1] * ( tc_leaf - df_i$topt[1]) ^ 2),
      anet_fit_s = anet_fit/df_i$aopt[1],
      tdiff_fit  = tc_leaf - df_i$topt[1])
    
    df_i <- 
      df_i %>% 
      mutate(fit_inst = list(df_2),
             data_org = purrr::map(data_org,
                                   ~mutate(.,
                                           tdiff_org   = Tleaf - df_i$topt[1],
                                           anet_org_s  = Photo/df_i$aopt[1] * 100)))
    
    
    df_out <- rbind(df_out, df_i)
  }
  
  ## Create dataframes for plotting
  # Get alphabetical order of cl for nice facet wrap
  cl_order <- sort(unique(df_plot$cl))
  
  df_unscaled <- 
    df_out %>% 
    dplyr::select(sitename, date, fit_inst) %>% 
    unnest(fit_inst) %>% 
    right_join(df_1 %>% unnest(rpm_inst)) %>% 
    rename(anet_rpm = anet) %>% 
    mutate(anet_rpm = anet_rpm * 1e6,
           id = paste0(sitename, date)) %>% 
    pivot_longer(cols = c('anet_rpm', 'anet_fit'), names_to = 'source', values_to = 'anet') %>% 
    dplyr::select(sitename, date, cl, tc_leaf, anet, source, id, rpm_sim_tc_opt) %>% 
    mutate(source = as.factor(source),
           cl = factor(cl, levels = cl_order),
           id = as.factor(id),
           sitename = as.factor(sitename))
  
  df_scaled <- 
    df_unscaled %>% 
    mutate(anet = anet/max(anet) * 100)
  
  ## Create function for plot basis
  make_pbase <- function(df_in){
    p <-
      df_in %>% 
      ggplot() +
      theme_linedraw() +
      theme(legend.position = c(0.9, 0.1),
            legend.justification = c(1, 0)) +
      scale_color_brewer(
        palette = "Dark2",
        labels = c("Observation", "Model")) + 
      scale_linetype_manual(
        values = c("solid", "dashed"),
        labels = c("Observation", "Model")) +
      labs(
        x = bquote(T[leaf] ~ "[°C]"),
        color = "Legend",
        linetype = "Legend"
      ) +
      facet_wrap(~cl) +
      geom_line(aes(
        x = tc_leaf, 
        y = anet, 
        group = interaction(id, source), 
        linetype = source,
        color = source
      )) 
    
    return(p)
  }
  
  ## Get unscaled plot
  p_unscaled <- 
    make_pbase(df_unscaled) +
    ylim(0, 30) +
    ylab(bquote(A[net] ~ "[µmol" ~ m ^-2 ~ s ^-1 ~ "]"))
  
  ## Get scaled plot
  p_scaled_old <- 
    make_pbase(df_scaled) +
    ylim(0, 100) +
    ylab(bquote(A[net] ~ "[%]"))
  
  #_____________________________________________________________________________
  # New approach for scaled plot
  
  df_scaled_new <- 
    df_out %>% 
    dplyr::select(sitename, date, fit_inst, data_org) %>% 
    unnest(fit_inst) %>% 
    right_join(df_1 %>% 
                 unnest(rpm_inst) %>% 
                 dplyr::select(-data_org)) %>% 
    rename(anet_rpm = anet) %>% 
    mutate(anet_rpm = anet_rpm * 1e6,
           id = paste0(sitename, date)) %>% 
    group_by(id) %>% 
    mutate(
      # Observed
      anet_fit_s = anet_fit_s / max(anet_fit_s) * 100,
      
      # Modelled
      anet_rpm_s = anet_rpm/max(anet_rpm) * 100,
      tdiff_rpm  = tc_leaf - rpm_sim_tc_opt,
      
      # Info
      cl = factor(cl, levels = cl_order),
      id = as.factor(id),
      sitename = as.factor(sitename)) %>% 
    ungroup() 
  
  df_org_s <-
    df_scaled_new %>% 
    dplyr::select(id, data_org, cl, sitename) %>% 
    distinct() %>% 
    unnest(data_org)
  
  #_____________________________________________________________________________
  # Facet by climate zone
  # Abbreviations:
  # org = original measurement data
  # fit = fitted to measurement data
  # rpm = modelled data
  
  p_scaled_cl <-
    df_scaled_new %>% 
    ggplot() +
    
    ## Adding measurement data as points  
    geom_point(data = df_org_s, 
               aes(tdiff_org, 
                   anet_org_s,
                   fill = "Measurement", 
                   color = "Measurement"
               ),
               alpha = 0.2,
               # size = 0.5,
               shape = 16) +
    
    ## Uncheck for densitiy plots on axes
    # geom_rug(data = df_org_s,
    #          aes(tdiff_org, anet_org_s),
    #          # fill = "grey",
    #          alpha = 0.2,
    #          sides = "b",
    #          # size = 0.5,
    #          color = "grey50",
    #          length = unit(0.02, "npc"),
    #          position = "jitter"
  #          ) +
  
  
  ## Uncheck for adding fitted parabolas
  # geom_line(data = df_scaled_new, 
  #           aes(tdiff_fit, anet_fit_s, group = id, linetype = "Measurement (loess)", color = "Measurement (loess)"),
  #           alpha = 0.5) +
  
  ## Add modelled response curves
  geom_line(data = df_scaled_new,
            aes(tdiff_rpm, 
                anet_rpm_s, 
                group = id, 
                linetype = "Model", 
                color = "Model",
                fill = "Model"),
            alpha = 0.5,
            size = 1) +
    
    ## Add smoother over measurement data
    geom_smooth(data = df_org_s, 
                aes(tdiff_org, 
                    anet_org_s, 
                    color = "Measurement (loess)", 
                    fill = "Measurement (loess)",
                    linetype = "Measurement (loess)"),
                # fill = "grey", 
                alpha = 0.5,
                size = 1,
                method = 'loess',
                level = 0.99,
                fullrange = T) +
    
    ## Aesthetics
    scale_linetype_manual(
      name = "",
      breaks = c("Model", "Measurement (loess)"),
      values = c("solid", "solid"),
      guide = "none") +
    scale_color_manual(
      name = "",
      breaks = c("Model", "Measurement (loess)", "Measurement"),
      values = c(RColorBrewer::brewer.pal(3, "Dark2")[1:2], "grey50")) +
    scale_fill_manual(
      name = "",
      breaks = c("Model", "Measurement (loess)", "Measurement"),
      values = c(RColorBrewer::brewer.pal(3, "Dark2")[1:2], "grey50"),
      guide = guide_legend(
        override.aes = list(shape = c(NA, NA, 19),
                            size = c(1, 1, 2.5),
                            linetype = c(1,1,0),
                            fill = c(NA, RColorBrewer::brewer.pal(3, "Dark2")[2], "white"),
                            color = c(RColorBrewer::brewer.pal(3, "Dark2")[1:2], "black")))) +
    facet_wrap(~cl, nrow = 2) +
    ylim(0, 140) +
    xlim(-20, 20) +
    labs(
      title = expression(bold("Observed and modelled temperature response curves")),
      y = bquote("Scaled " ~ A[net] ~ "[%]"),
      x = bquote(T[leaf] ~ "-" ~ T[opt] ~ "[°C]")) +
    theme_linedraw(base_size = 14) +
    theme(
      # legend.position = c(0.9, 0.05),
      # legend.justification = c(1, 0),
      legend.position = "bottom",
      # strip.background =element_rect(fill="white"),
      # strip.text = element_text(colour = 'black', face = "bold"),
      panel.grid = element_blank(),
      plot.subtitle = element_text(hjust = 0, face = "bold"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0, size = 14)) 
  
  
  #_____________________________________________________________________________
  # Facet by climate zone (but coarser cl level)
  p_scaled_cl_2 <-
    df_scaled_new %>% 
    mutate(cl = as.factor(str_extract(cl, "[A-Z]{1,2}"))) %>% 
    mutate(across(cl, factor, levels=c('ET', 'D', 'C', 'A'))) %>% 
    ggplot() +
    geom_point(data = df_org_s %>% 
                 mutate(cl = as.factor(str_extract(cl, "[A-Z]{1,2}"))) %>% 
                 mutate(across(cl, factor, levels=c('ET', 'D', 'C', 'A'))), 
               aes(tdiff_org, 
                   anet_org_s,
                   fill = "Measurement", 
                   color = "Measurement"
               ),
               alpha = 0.2,
               # size = 0.5,
               shape = 16) +
    # geom_rug(data = df_org_s,
    #          aes(tdiff_org, anet_org_s),
    #          # fill = "grey",
    #          alpha = 0.2,
    #          sides = "b",
    #          # size = 0.5,
    #          color = "grey50",
    #          length = unit(0.02, "npc"),
    #          position = "jitter"
    #          ) +
    # geom_line(data = df_scaled_new, 
  #           aes(tdiff_fit, anet_fit_s, group = id, linetype = "Measurement (loess)", color = "Measurement (loess)"),
  #           alpha = 0.5) +
  geom_line(data = df_scaled_new %>% mutate(cl = as.factor(str_extract(cl, "[A-Z]{1,2}"))),
            aes(tdiff_rpm, anet_rpm_s, group = id, 
                linetype = "Model", 
                color = "Model",
                fill = "Model"),
            alpha = 0.5,
            size = 1) +
    geom_smooth(data = df_org_s %>% mutate(cl = as.factor(str_extract(cl, "[A-Z]{1,2}"))), 
                aes(tdiff_org, anet_org_s, 
                    color = "Measurement (loess)", 
                    fill = "Measurement (loess)",
                    linetype = "Measurement (loess)"),
                # fill = "grey", 
                alpha = 0.5,
                size = 1,
                method = 'loess',
                level = 0.99,
                fullrange = T) +
    scale_linetype_manual(
      name = "",
      breaks = c("Model", "Measurement (loess)"),
      values = c("solid", "solid"),
      guide = "none") +
    scale_color_manual(
      name = "",
      breaks = c("Model", "Measurement (loess)", "Measurement"),
      values = c(RColorBrewer::brewer.pal(3, "Dark2")[1:2], "grey50")) +
    scale_fill_manual(
      name = "",
      breaks = c("Model", "Measurement (loess)", "Measurement"),
      values = c(RColorBrewer::brewer.pal(3, "Dark2")[1:2], "grey50"),
      guide = guide_legend(
        override.aes = list(shape = c(NA, NA, 19),
                            size = c(1, 1, 2.5),
                            linetype = c(1,1,0),
                            fill = c(NA, RColorBrewer::brewer.pal(3, "Dark2")[2], "white"),
                            color = c(RColorBrewer::brewer.pal(3, "Dark2")[1:2], "black")))) +
    facet_wrap(~cl, 
               nrow = 1,
               labeller = as_labeller(c(
                 'A'  = "Tropical",
                 'C'  = "Temperate",
                 'D'  = "Boreal",
                 'ET' = "Arctic"
               ))) +
    ylim(0, 140) +
    xlim(-20, 20) +
    labs(
      title = expression(bold("Observed and modelled temperature response curves")),
      y = bquote("Scaled " ~ A[net] ~ "[%]"),
      x = bquote(T[leaf] ~ "-" ~ T[opt] ~ "[°C]")) +
    theme_linedraw(base_size = 14) +
    theme(
      # legend.position = c(0.9, 0.05),
      # legend.justification = c(1, 0),
      legend.position = "bottom",
      # strip.background =element_rect(fill="white"),
      # strip.text = element_text(colour = 'black', face = "bold"),
      panel.grid = element_blank(),
      plot.subtitle = element_text(hjust = 0, face = "bold"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0, size = 14))
  
  #_______________________________________________________________________________
  # Facet by climate (coarse) with parabolas
  p_scaled_cl_2_parab <-
    df_scaled_new %>% 
    mutate(cl = as.factor(str_extract(cl, "[A-Z]{1,2}"))) %>% 
    ggplot() +
    
    ## Adding measurement data as points  
    geom_point(data = df_org_s %>% mutate(cl = as.factor(str_extract(cl, "[A-Z]{1,2}"))), 
               aes(tdiff_org, 
                   anet_org_s,
                   fill = "Measurement", 
                   color = "Measurement"
               ),
               alpha = 0.2,
               # size = 0.5,
               shape = 16) +
    
    ## Uncheck for densitiy plots on axes
    # geom_rug(data = df_org_s,
    #          aes(tdiff_org, anet_org_s),
    #          # fill = "grey",
    #          alpha = 0.2,
    #          sides = "b",
    #          # size = 0.5,
    #          color = "grey50",
    #          length = unit(0.02, "npc"),
    #          position = "jitter"
  #          ) +
  
  ## Add smoother over measurement data
  # geom_smooth(data = df_org_s, 
  #             aes(tdiff_org, anet_org_s, 
  #                 color = "Measurement (loess)", 
  #                 fill = "Measurement (loess)",
  #                 linetype = "Measurement (loess)"),
  #             # fill = "grey", 
  #             alpha = 0.5,
  #             size = 1,
  #             method = 'loess',
  #             level = 0.99,
  #             fullrange = T) +
  
  ## Add fitted parabolas
  geom_line(data = df_scaled_new %>% mutate(cl = as.factor(str_extract(cl, "[A-Z]{1,2}"))),
            aes(tdiff_fit, 
                anet_fit_s, 
                group = id, 
                linetype = "Measurement (fit)", 
                color = "Measurement (fit)",
                fill = "Measurement (fit)"),
            size = 1, 
            alpha = 0.5) +
    
    ## Add modelled response curves
    geom_line(data = df_scaled_new %>% mutate(cl = as.factor(str_extract(cl, "[A-Z]{1,2}"))),
              aes(tdiff_rpm, 
                  anet_rpm_s, 
                  group = id, 
                  linetype = "Model", 
                  color = "Model",
                  fill = "Model"),
              alpha = 0.5,
              size = 1) +
    
    ## Aesthetics
    scale_linetype_manual(
      name = "",
      breaks = c("Model", "Measurement (fit)"),
      values = c("solid", "solid"),
      guide = "none") +
    scale_color_manual(
      name = "",
      breaks = c("Model", "Measurement (fit)", "Measurement"),
      values = c(RColorBrewer::brewer.pal(3, "Dark2")[1:2], "grey50")) +
    scale_fill_manual(
      name = "",
      breaks = c("Model", "Measurement (fit)", "Measurement"),
      values = c(RColorBrewer::brewer.pal(3, "Dark2")[1:2], "grey50"),
      guide = guide_legend(
        override.aes = list(shape = c(NA, NA, 19),
                            size = c(1, 1, 2.5),
                            linetype = c(1,1,0),
                            fill = c(NA, RColorBrewer::brewer.pal(3, "Dark2")[2], "white"),
                            color = c(RColorBrewer::brewer.pal(3, "Dark2")[1:2], "black")))) +
    facet_wrap(~cl, 
               nrow = 1,
               labeller = as_labeller(c(
                 'A'  = "Tropical",
                 'C'  = "Temperate",
                 'D'  = "Boreal",
                 'ET' = "Arctic"
               ))) +
    ylim(0, 140) +
    xlim(-20, 20) +
    labs(
      title = expression(bold("Observed and modelled temperature response curves")),
      y = bquote("Scaled " ~ A[net] ~ "[%]"),
      x = bquote(T[leaf] ~ "-" ~ T[opt] ~ "[°C]")) +
    theme_linedraw(base_size = 14) +
    theme(
      # legend.position = c(0.9, 0.05),
      # legend.justification = c(1, 0),
      legend.position = "bottom",
      # strip.background =element_rect(fill="white"),
      # strip.text = element_text(colour = 'black', face = "bold"),
      panel.grid = element_blank(),
      plot.subtitle = element_text(hjust = 0, face = "bold"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0, size = 14))
  
  #_____________________________________________________________________________
  # Facet by site
  p_scaled_site <- p_scaled_cl + facet_wrap( ~ sitename, nrow = 2)
  
  #_____________________________________________________________________________
  ## Return
  return(
    list(
      p_unscaled   = p_unscaled,
      p_scaled_cl   = p_scaled_cl,
      p_scaled_cl_2 = p_scaled_cl_2,
      p_scaled_cl_2_parab = p_scaled_cl_2_parab,
      p_scaled_site = p_scaled_site, 
      p_scaled_old = p_scaled_old
    )
  )
}


## Standard x ~ y plot ----
plot_obs_vs_pred <- function(df_in, 
                             x, 
                             y,
                             x_se = NULL, 
                             y_se = NULL, 
                             fill = NULL,
                             shape = NULL,
                             yxlim = NULL,
                             xlim  = NULL,
                             ylim   = NULL,
                             plot11 = T,
                             add_lm_metrics = T) {
  
  ## Define local dataframe
  df_temp <- df_in
  
  ## Add generic variables to df
  df_temp$y <- df_temp[[y]]
  df_temp$x <- df_temp[[x]]
  
  ## Drop NA if needed
  if (nrow(df_temp) != nrow(df_temp %>% drop_na(x,y))) {
    warning("[!] NA values in either x or y variable! NA entries are removed")
    df_temp <- df_temp %>% drop_na(x,y)
  }
  
  ## Attach standard error variable either based on character naming
  if (isTRUE(y_se)) {
    warning("[!] y_se needs to be defined from now on!")
  }
  if (isTRUE(x_se)) {
    warning("[!] x_se needs to be defined from now on!")
  }
  
  if (!(is.null(y_se))) {
    df_temp$y_se <- df_temp[[paste0(y_se)]]
  } else {
    df_temp$y_se <- 0
    y_se <- "add_se"
  }
  
  if (!(is.null(x_se))) {
    df_temp$x_se <- df_temp[[paste0(x_se)]]
  } else {
    df_temp$x_se <- 0
    x_se <- "add_se"
  }
  
  # Check if fill is required, if not use basic fill
  if (!(is.null(fill))) {
    ## Re-leveling factor levels to ensure correct coloring
    if (str_detect(fill, "cl")) {
      df_temp$cl <- factor(df_temp$cl, levels = c(df_temp$cl %>% str_sort() %>% unique()))
    }
    if (!(fill %in% names(df_temp))) {
      # Check if given fill variable is in df_in
      rlang::warn("> Warning: Fill variable not found in data frame!")
      fill <- NULL
    } else {
      df_temp <- df_temp %>% arrange(eval(parse(text = fill)))
    }
  } else {
    # fill <- "black"
  }
  
  # Check if shape is required, if not, use basic shape
  if (!(is.null(shape))) {
    if (!(shape %in% names(df_temp))) {
      # Check if given shape variable is in df_in
      rlang::warn("> Warning: shape variable not found in data frame!")
      shape <- NULL
    }
  } else {
    # shape <- 21 
  }
  
  ## Overwrite shape metric if observed values are entered (no difference betwwen Ac and Aj limitation)
  if (str_detect(y, "eva_data") & x == "temp" & !(is.null(shape))) {
    df_temp[[shape]] <- "ac"
  }
  
  ## Set linetype for linear model (will be overwritten if lm is significant)
  lm_ltyp <-  "dashed"
  p_thresh <- 0.01 # global p value threshold
  ci_thresh <- 0.99 # global CI threshold
  
  # Linear regression metrics ----
  if (T) { #(plot11) {
    n <- nrow(df_temp)
    fit <- lm(y ~ x, data = df_temp)
    sry <- summary(fit)
    r2 <- sry$adj.r.squared %>% round(2)
    # r2      <- round(cor(df_temp$x, df_temp$y, method = "pearson"),2)
    
    # Catch error if there is no linear regression fitable
    sry_failed <- TRUE
    not_sign   <- TRUE
    
    if (length(sry$coefficients) > 4 & !is.nan(sry$coefficients[1, 4])) {
      
      sry_failed <- FALSE
      
      ## Get relevant metrics
      p_val <- round(sry$coefficients[[8]], 5)
      rmse <- sqrt((c(crossprod(sry$residuals)) / length(sry$residuals))) %>% round(2)
      slope <- sry$coefficients[2, 1] %>% round(2)
      b0 <- sry$coefficients[1, 1] %>% round(2)
      b1 <- sry$coefficients[2, 1] %>% round(2)
      b0_sign <- sry$coefficients[1, 4] <= p_thresh # T = Intercept is signficantly different from 0
      b1_sign <- !(between(1, confint(fit, level = ci_thresh)[2, 1] %>% round(2), confint(fit, level = ci_thresh)[2, 2] %>% round(2))) # T = Slope is significantly different from 1
      bias <- mean(df_temp$y - df_temp$x) %>% round(2)
      
      ## Get confidence intervall
      ci <- round(confint(fit, level = ci_thresh), 4)
      
      b0_ci_lo <- ci[1 ,1]
      b0_ci_up <- ci[1, 2]
      
      b1_ci_lo <- ci[2, 1]
      b1_ci_up <- ci[2, 2]
      
      
      ## Make tibble for linear regression plus confidence interval ribbon
      step_accuracy <- 0.01
      
      # data_lm <- 
      #   tibble(x = seq(0, 40, step_accuracy), 
      #          y = (b1 * x + b0),
      #          ylo = (b1_ci_lo * x + b0_ci_up),
      #          yup = (b1_ci_up * x + b0_ci_lo)
      #   ) %>% 
      #   rowwise() %>% 
      #   mutate(ymax = max(ylo, yup),
      #          ymin = min(ylo, yup)) %>% 
      #   ungroup()
      
      ## Dataframe with linear regression
      data_lm <- tibble(x = seq(0, 40, step_accuracy),
                        y = b1 * x + b0)
      
      data_lm$y <- round(data_lm$y/step_accuracy) * step_accuracy
      
      ## Add lower confidence interval line
      data_lm$y_ci_lo <- 
        predict(fit, 
                newdata = data.frame(
                  x = seq(0, 40, step_accuracy)
                ), 
                interval = "confidence", 
                level = ci_thresh)[, 2] %>% 
        as.vector()
      
      data_lm$y_ci_lo <- round(data_lm$y_ci_lo/step_accuracy) * step_accuracy
      
      ## Add upper confidence interval line
      data_lm$y_ci_up <- 
        predict(fit, 
                newdata = data.frame(
                  x = seq(0, 40, step_accuracy)
                ), 
                interval = "confidence", 
                level = ci_thresh)[, 3] %>% 
        as.vector()
      
      data_lm$y_ci_up <- round(data_lm$y_ci_up/step_accuracy) * step_accuracy
      
      ## Calculate intersection point with 1-1 line (x = y) and CI for it
      ip_xy <- 
        data_lm %>% 
        dplyr::filter(x == y) %>% 
        pull(x) %>% 
        mean()
      
      ip_xmin <-
        data_lm %>% 
        dplyr::filter(x == y_ci_lo) %>% 
        pull(x) %>% 
        min()
      
      ip_xmax <- 
        data_lm %>% 
        dplyr::filter(x == y_ci_up) %>% 
        pull(x) %>% 
        max()
      
      ci_lo_at_ip <-
        data_lm %>% 
        dplyr::filter(x == y) %>% 
        pull(y_ci_lo) %>% 
        min()
      
      ci_up_at_ip <- 
        data_lm %>% 
        dplyr::filter(x == y) %>% 
        pull(y_ci_up) %>% 
        max()
      
      ## Overwrite linetype if not significant
      if (p_val <= p_thresh) {
        lm_ltyp  <-  "solid"
        not_sign <- FALSE
      }
      
    } else {
      p_val <- NA
      rmse <- NA
      slope <- NA
      b0_sign <- NA
      b1_sign <- NA
      bias <-NA
      ip_xy <- NA
      ip_xmin <- NA
      ip_xmax <- NA
      ci_lo_at_ip <- NA
      ci_up_at_ip <- NA
      b0 <- NA
      b1 <- NA
      b0_ci_up <- NA
      b0_ci_lo <- NA
      b1_ci_up <- NA
      b1_ci_lo <- NA
    }
    
    df_metrics <- tibble(
      # Intercept
      b0_ci_up = b0_ci_up,
      b0 = b0,
      b0_ci_lo = b0_ci_lo,
      b0_sign_not_0 = b0_sign,
      
      # Slope
      b1_ci_lo = b1_ci_lo,
      slope = slope,
      b1 = b1,
      b1_ci_up = b1_ci_up,
      b1_sign_not_one = b1_sign,
      
      # Intersection
      ip_xmin = ip_xmin,
      ip_xy = ip_xy,
      ip_xmax = ip_xmax,
      ci_lo_at_ip = ci_lo_at_ip,
      ci_up_at_ip = ci_up_at_ip,
      
      # Metrics
      pval = p_val,
      r2 = r2,
      rmse = rmse,
      bias = bias,
      n = n
    )
  } else {
    df_metrics <- tibble()
  }
  
  ## Add information to df_metrics
  df_metrics$max_x <- round(max(df_temp$x), 2)
  df_metrics$max_y <- round(max(df_temp$y), 2)
  df_metrics$min_x <- round(min(df_temp$x), 2)
  df_metrics$min_y <- round(min(df_temp$y), 2)
  
  df_metrics$mean_x <- round(mean(df_temp$x), 2)
  df_metrics$sd_x <- round(sd(df_temp$x), 2)
  
  df_metrics$mean_y <- round(mean(df_temp$y), 2)
  df_metrics$sd_y <- round(sd(df_temp$y), 2)
  
  #_____________________________________________________________________________
  # Start Plot ----
  plot <- 
    df_temp %>%
    drop_na(x, y) %>%
    ggplot()
  
  ## Add 1:1 aesthetics ----
  if (plot11) {
    
    lim <- max(max(df_temp$y), max(df_temp$x))
    
    plot <- plot +
      ylim(0, ifelse(is.null(yxlim), lim, yxlim)) +
      xlim(0, ifelse(is.null(yxlim), lim, yxlim)) +
      geom_abline(linetype = "dotted")
    
  } else {
    plot <- plot +
      ylim(0, max(ylim, df_temp$y)) +
      xlim(0, max(xlim, df_temp$x))
  }
  
  ## Add uncertainty range of intersection point if required ----
  if (str_detect(y, "tc_opt") & x == "temp") {
    # plot <-
    #   plot +
    #   geom_vline(xintercept = ip_xy,
    #              color = "grey60",
    #              linetype = "longdash") +
    #   annotate('ribbon',
    #            y = c(-Inf, Inf),
    #            xmin = ip_xmin,
    #            xmax = ip_xmax,
    #            alpha = 0.1, fill = 'grey10')

    plot <-
      plot +
      geom_vline(xintercept = ip_xy,
                 color = "grey60",
                 linetype = "dashed") +
      annotate('ribbon',
               y = c(-Inf, Inf),
               xmin = ip_xmin,
               xmax = ip_xmax,
               alpha = 0.1,
               fill = 'grey10')
  }
  
  ## Add error bars before points ----
  if (!(is.null(y_se))) {
    if (is.null(fill)) {
      plot <- plot + geom_errorbar(aes(x = x, y = y, ymin = y - y_se, ymax = y + y_se), 
                                   width = 0,
                                   alpha = 0.8)
    } else {
      plot <- plot + geom_errorbar(aes(x = x, y = y, ymin = y - y_se, ymax = y + y_se, color = eval(parse(text = fill))), 
                                   width = 0,
                                   alpha = 0.8)
    }
  }
  
  if (!(is.null(x_se))) {
    if (is.null(fill)) {
      plot <- plot + geom_errorbarh(aes(y = y, x = x, xmin = x - x_se, xmax = x + x_se), 
                                    width = 0,
                                    alpha = 0.8)
    } else {
      plot <- plot + geom_errorbarh(aes(y = y, x = x, xmin = x - x_se, xmax = x + x_se, color = eval(parse(text = fill))), 
                                    width = 0,
                                    alpha = 0.8)
    }
  }
  
  ## Add points with right coloring and shape ----
  if (is.null(fill) & is.null(shape)) {
    plot <- plot +
      geom_point(aes(y = y, x = x), color = "black", size = 2.5, pch = 21)
    
  } else {
    if (str_detect(fill, "cl") & !is.null(shape)) {
      plot <-
        plot +
        geom_point(
          aes(
            y = y, x = x,
            fill = eval(parse(text = fill)),
            shape = eval(parse(text = shape))),
          color = "black",
          size = 2
        )
      
      ## Add geom_smooth only for significant relationships
      # if (lm_ltyp == "solid") plot <- plot + gg_smooth
      
    } else if (str_detect(fill, "cl")) {
      plot <-
        plot +
        geom_point(
          aes(
            y = y, x = x,
            fill = eval(parse(text = fill))),
          shape = 21,
          color = "black",
          size = 2
        )
      
      ## Add geom_smooth only for significant relationships
      # if (lm_ltyp == "solid") plot <- plot + gg_smooth
      
    } else {
      plot <-
        plot +
        geom_point(aes(y = y, x = x, fill = eval(parse(text = fill))), color = "black", size = 2.5, pch = 21)
      
    }
  }
  
  ## Add linear model with confidence interval for significant relationships ----
  if (!not_sign) {
    plot <-
      plot +
      geom_ribbon(
        data = data_lm,
        aes(x = x, 
            ymin = y_ci_lo,
            ymax = y_ci_up),
        alpha = 0.25
      ) +       
      geom_line(
        data = data_lm,
        aes(x = x, 
            y = y)
      )
  }
  
  
  ## Add text and color to plot ----
  
  ## Add labels
  plot <- plot +
    ylab(paste(y)) +
    xlab(paste(x)) +
    theme_classic() +
    theme(legend.position = "bottom")
  
  ## Add model metrics if needed
  caption11 <- ""
  
  if (plot11) {
    caption11 <- paste0("Is b0 ~ 0?: ", ifelse(b0_sign, "No", "Yes"), " | Is b1 ~ 1?: ", ifelse(b1_sign, "No", "Yes"), " | #N: ", n)
    
    if (add_lm_metrics & !is.na(p_val <= p_thresh | p_val)) {
      txt_b0   <- deparse(bquote("Intercept: " * .(round(b0, 2)) * " ["  * .(round(b0_ci_lo, 2))  * ", "  * .(round(b0_ci_up, 2)) * "]"))
      txt_b1   <- deparse(bquote("Slope: " * .(round(b1, 2)) * " ["  * .(round(b1_ci_lo, 2))  * ", "  * .(round(b1_ci_up, 2)) * "]"))
      # txt_pv   <- ifelse(p_val > 0.05, "p > 0.05 ", "p ≤ 0.05")
      txt_r2 <- deparse(bquote(R^2 * ": " * .(r2)))
      txt_rmse <- deparse(bquote("Bias: " * .(bias)))
      txt_bias <- deparse(bquote("RMSE: " * .(rmse)))
      
      ## Distribute text equally in top left third
      n_texts <- 5
      rel_pos <- yxlim/3/n_texts # relative difference between lines
      
      pos_0 <- yxlim - 0 * rel_pos
      pos_1 <- yxlim - 1 * rel_pos
      pos_2 <- yxlim - 2 * rel_pos
      pos_3 <- yxlim - 3 * rel_pos
      pos_4 <- yxlim - 4 * rel_pos
      
      plot <-
        plot +
        # annotate(geom = "text", x = 0, y = yxlim*0.8,  size = 3.5, vjust = 1, hjust = 0, label = txt_pv) +
        annotate(geom = "text", x = 0, y = pos_0,  size = 2.75, vjust = 1, hjust = 0, label = txt_b0, parse = T) +
        annotate(geom = "text", x = 0, y = pos_1,  size = 2.75, vjust = 1, hjust = 0, label = txt_b1, parse = T) +
        annotate(geom = "text", x = 0, y = pos_2,  size = 2.75, vjust = 1, hjust = 0, label = txt_r2, parse = T) +
        annotate(geom = "text", x = 0, y = pos_3,  size = 2.75, vjust = 1, hjust = 0, label = txt_rmse, parse = T) +
        annotate(geom = "text", x = 0, y = pos_4,  size = 2.75, vjust = 1, hjust = 0, label = txt_bias, parse = T)
    }
    
    plot <-
      plot +
      labs(
        title = paste0("Pred-Obs of ", str_remove(x, "rpm_acc_")),
        # subtitle = bquote(R^2 ~ " = " ~ .(r2) ~ " | bias = " ~ .(bias) ~ " | RMSE = " ~ .(rmse) ~ " | slope = " ~ .(slope)),
        # subtitle = bquote(R^2 ~ " = " ~ .(r2) ~ " | bias = " ~ .(bias) ~ " | RMSE = " ~ .(rmse)), # Without slope because it is in ggpmisc annotation
        caption = caption11
      )
  }
  
  ## Add climate zone coloring if needed
  if (!is.null(fill)) {
    if (fill == "cl") {
      
      ## - ORIGINAL KG ZONE COLORS
      ## Get subset of values and colors
      # all <- tibble(
      #   climate_colors = c("#960000", "#FF0000", "#FF6E6E", "#FFCCCC", "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64", "#007800", "#005000", "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00", "#550055", "#820082", "#C800C8", "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF", "#64FFFF", "#F5FFFF"),
      #   cl = c("Af", "Am", "As", "Aw", "BSh", "BSk", "BWh", "BWk", "Cfa", "Cfb", "Cfc", "Csa", "Csb", "Csc", "Cwa", "Cwb", "Cwc", "Dfa", "Dfb", "Dfc", "Dfd", "Dsa", "Dsb", "Dsc", "Dsd", "Dwa", "Dwb", "Dwc", "Dwd", "EF", "ET", "Ocean")
      # ) %>%
      #   arrange(cl)
      # sub_1 <- df_temp$cl %>%
      #   str_sort() %>%
      #   unique()
      # sub_2 <- all[which(all$cl %in% sub_1), ]
      # 
      # plot <-
      #   plot +
      #   scale_fill_manual(
      #     values = c(sub_2$climate_colors),
      #     name = "Köppen-Geiger Climate Zone",
      #     guide = guide_legend(direction = "horizontal",
      #                          title.position = "top",
      #                          ncol = 12,
      #                          override.aes = list(shape = 21),
      #                          reverse = TRUE)
      #   ) +
      #   scale_color_manual(
      #     values = c(sub_2$climate_colors),
      #     name = "Köppen-Geiger Climate Zone",
      #     guide = guide_legend(direction = "horizontal",
      #                          title.position = "top",
      #                          ncol = 12,
      #                          override.aes = list(shape = 21),
      #                          reverse = TRUE)
      #   )
      ## - ORIGINAL KG ZONE COLORS
      
      ## - COLORBLIND FRIENDLY COLORS
      # Brewer: RdYlBu has too little contrast.
      # Viridis: turbo does not work well
      plot <-
        plot +
       scale_color_viridis_d(
         option = "magma",
         direction = -1,
         # begin = 0.25, # For making better dark mode plots
         # end  = 0.85,  # Form making better light mode plots
         name = "Köppen-Geiger Climate Zone",
         guide = guide_legend(direction = "horizontal",
                              title.position = "top",
                              ncol = 12,
                              override.aes = list(shape = 21),
                              reverse = TRUE)
        ) +
        scale_fill_viridis_d(
          option = "magma",
          direction = -1,
          # begin = 0.25, # For making better dark mode plots
          # end   = 0.85,  # Form making better light mode plots
          name = "Köppen-Geiger Climate Zone",
          guide = guide_legend(direction = "horizontal",
                               title.position = "top",
                               ncol = 12,
                               override.aes = list(shape = 21),
                               reverse = TRUE)
        )
      
      # # Scico (suitable for dark plots):
      # plot <-
      #   plot +
      #   scico::scale_color_scico_d(
      #     palette = "roma",
      #     # direction = -1,
      #     # begin = 0.25, # For making better dark mode plots
      #     # end  = 0.85,  # Form making better light mode plots
      #     name = "Köppen-Geiger Climate Zone",
      #     guide = guide_legend(direction = "horizontal",
      #                          title.position = "top",
      #                          ncol = 12,
      #                          override.aes = list(shape = 21),
      #                          reverse = TRUE)
      #   ) +
      #   scico::scale_fill_scico_d(
      #     palette = "roma",
      #     # direction = -1,
      #     # begin = 0.25, # For making better dark mode plots
      #     # end   = 0.85,  # Form making better light mode plots
      #     name = "Köppen-Geiger Climate Zone",
      #     guide = guide_legend(direction = "horizontal",
      #                          title.position = "top",
      #                          ncol = 12,
      #                          override.aes = list(shape = 21),
      #                          reverse = TRUE)
      #   )
      
      ## - COLORBLIND FRIENDLY COLORS
      
    } else {
      plot <-
        plot +
        scale_fill_viridis_d(guide = guide_legend(direction = "horizontal", title.position = "top", ncol = 12, override.aes = list(shape = 21))) +
        scale_color_viridis_d(guide = guide_legend(direction = "horizontal", title.position = "top", ncol = 12, override.aes = list(shape = 21)))
    }
  }
  
  ## Add shape legend
  if (!(is.null(shape))) {
    
    df_metrics$n_aj <- NA
    df_metrics$n_ac <- NA
    
    if (!is.null(shape)) {
      if (str_detect(shape, "rpm_sim_min_a")) {
        
        if (all(c("ac", "aj") %in% df_temp[["rpm_sim_min_a"]])) {
          n_ac <- table(df_temp[[shape]])[["ac"]]
          if (n_ac == nrow(df_temp)) {
            n_aj <- 0
          } else {
            n_aj <- table(df_temp[[shape]])[["aj"]]
          }
        } else if ("ac" %in% df_temp[["rpm_sim_min_a"]]) {
          n_ac <- nrow(df_temp)
          n_aj <- 0
        } else if ("ac" %in% df_temp[["rpm_sim_min_a"]]) {
          n_ac <- 0
          n_aj <- nrow(df_temp)
        } else {
          stop("Something went wrong with adding shape to legend.")
        }
        
        df_metrics$n_aj <- n_aj
        df_metrics$n_ac <- n_ac
        
        plot <-
          plot +
          geom_point(
            data = tibble(x = -10, y = -10, s = as.factor(ifelse(n_ac == 0, "ac", "aj"))),
            aes(x, y, shape = s)) +
          scale_shape_manual(
            name = "Limit.",
            breaks = c("ac", "aj"),
            values = c(21, 23),
            labels = c("Ac", "Aj"),
            guide = guide_legend(
              direction = "horizontal", 
              title.position = "top",
              override.aes = list(
                shape = c(21, 23),
                size = c(2, 2)
              ))) +
          labs(caption = paste0(caption11, " | #Ac: ", n_ac, " | #Aj: ", n_aj))
      }
    }
  }
  
  out <- list(
    plot = plot,
    df_metrics = df_metrics,
    df_in = df_in
  )
  
  return(out)
}

## Collecting plots ----
plot_all_final_plots_from_df_plot <- function(df_plot, settings){

  #_____________________________________________________________________________
  ## Get plots
  ## Modobs and var ~ T_growth
  p_all    <- plot_all_for_one_model(df_plot)
  row_topt <- p_all$modobs_topt + p_all$p_topt_mod_tcair + p_all$p_topt_obs_tcair + plot_layout(nrow = 1) & theme(legend.position = "none") 
  row_aopt <- p_all$modobs_aopt + p_all$p_aopt_mod_tcair + p_all$p_aopt_obs_tcair + plot_layout(nrow = 1) & theme(legend.position = "bottom")
  
  row_agro <- p_all$modobs_agrowth + p_all$p_agrowth_mod_tcair + p_all$p_agrowth_obs_tcair + plot_layout(nrow = 1) & theme(legend.position = "bottom")
  row_tspan <- p_all$modobs_tspan  + p_all$p_tspan_mod_tcair    + p_all$p_tspan_obs_tcair    + plot_layout(nrow = 1) & theme(legend.position = "bottom")
  
  p_1 <- row_topt / row_aopt + patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")
  p_2 <- row_agro / row_tspan + patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  # THERMAL RESPONSES
  p_tr <- plot_shape_thermal_response(df_plot)
  
  ## Big plot with all thermal response curves
  p_big_plot <- plot_topt_extraction(df_plot, add_simulation = T, dir = settings$dir_now, make_one_plot = T)
  
  
  #_____________________________________________________________________________
  # Save plots only when requested
  if (F & settings$save_plots) { # Set to F to save only the big plots within the pmodel routine
  
    message(paste0("Saving to... ", settings$dir_now))
    
    ## Base Plots
    ggsave(paste0(settings$dir_now, "/all_final_topt_aopt.pdf"), p_1, height = 15, width = 25)
    ggsave(paste0(settings$dir_now, "/all_final_agro_tspan.pdf"), p_2, height = 15, width = 25)
    
    ## Big plot with all thermal response curves
    ggsave(paste0(settings$dir_now, "/all_response_curves.pdf"), p_big_plot$p_big_plot, height = 15, width = 32)
    
    ## Thermal Response Curves
    # By climate zone
    ggsave( paste0(settings$dir_now, "/thermal_response_scaled_cl.pdf"), p_tr$p_scaled_cl, height = 6, width = 12)
    # By climate zone (coarser)
    ggsave( paste0(settings$dir_now, "/thermal_response_scaled_cl_2.pdf"), p_tr$p_scaled_cl_2, height = 4, width = 16)
    # By climate zone (coarser) with fitted parabolas
    ggsave( paste0(settings$dir_now, "/thermal_response_scaled_cl_2_parabolas.pdf"), p_tr$p_scaled_cl_2_parab, height = 4, width = 16)
    # By site
    ggsave(paste0(settings$dir_now, "/thermal_response_scaled_site.pdf"),p_tr$p_scaled_site,height = 6,  width = 12)
    
    ## Traits
    ggsave(paste0(settings$dir_now, "/traits_vs_tgrowth.pdf"), p_all$p_traits,  height = 5, width = 15)
    
  }
  
  out <- append(append(p_all, p_tr), p_big_plot)
  
  return(out)
}

plot_all_for_one_model <- function(df_plot, title_1 = "model_1") {
  
  #___________________________________________________________________________
  # ├ Modobs ----
  ## T_opt
  tmp <-
    plot_obs_vs_pred(
      df_in = df_plot,
      x = "rpm_sim_tc_opt",
      y = "eva_data_tc_opt",
      x_se = NULL,
      y_se = "eva_data_tc_opt_se",
      fill = "cl",
      shape = NULL,
      yxlim = 40)
  
  p_modobs_topt <- 
    tmp$plot  +
    xlab(bquote("Modelled "   ~ T[opt] ~ " [°C]")) +
    ylab(bquote("Observed "  ~ T[opt] ~ " [°C]")) +
    ggtitle(bquote("Modobs " ~ T[opt]))
  
  metrics_modobs_topt <- tmp$df_metrics
  
  ## A_opt
  tmp <- 
    plot_obs_vs_pred(
      df_plot,
      x = "rpm_sim_anet_opt", 
      y = "eva_data_anet_opt",
      x_se = NULL, 
      y_se = "eva_data_anet_opt_se",
      fill = "cl", 
      shape = NULL,
      yxlim = 30)
  
  p_modobs_aopt <- 
    tmp$plot  +
    xlab(bquote("Modelled "   ~ A[opt] ~ "[µmol" ~ m^-2 ~ s^-1 ~ "]")) +
    ylab(bquote("Observed "  ~ A[opt] ~ "[µmol" ~ m^-2 ~ s^-1 ~ "]")) +
    ggtitle(bquote("Modobs " ~ A[opt]))
  
  metrics_modobs_aopt <- tmp$df_metrics
  
  ## A_growth
  tmp <- 
    plot_obs_vs_pred(
      df_in = df_plot %>% 
        ## Adding SE of 0 only for matching legends with other modobs plots
        ## This does not change any metrics
        mutate(eva_data_agrowth_se = 0),
      x = "rpm_sim_anet_growth", 
      y = "eva_data_agrowth",
      x_se = NULL, 
      y_se = "eva_data_agrowth_se",
      fill = "cl", 
      shape = NULL,
      yxlim = 30)
  
  p_modobs_agrowth <- 
    tmp$plot +
    xlab(bquote("Modelled "   ~ A[growth] ~ "[µmol m¯² s¯¹ ]")) +
    ylab(bquote("Observed "  ~ A[growth] ~ "[µmol m¯² s¯¹ ]")) +
    ggtitle(bquote("Modobs " ~ A[growth]))
  
  metrics_modobs_agrowth <- tmp$df_metrics
  
  ## T_span
  tmp <- 
    plot_obs_vs_pred(
      df_in = df_plot %>% 
        ## Adding SE of 0 only for matching legends with other modobs plots
        ## This does not change any metrics
        mutate(eva_data_tspan_se = 0), 
      x = "rpm_sim_tspan", 
      y = "eva_data_tspan",
      x_se = NULL, 
      y_se = "eva_data_tspan_se",
      fill = "cl", 
      shape = NULL,
      yxlim = 30)
  
  p_modobs_tspan <- 
    tmp$plot +
    xlab(bquote("Modelled "   ~ T[span] ~ " [°C]")) +
    ylab(bquote("Observed "  ~ T[span] ~ " [°C]")) +
    ggtitle(bquote("Modobs " ~ T[span]))
  
  metrics_modobs_tspan <- tmp$df_metrics
  
  
  #_____________________________________________________________________________
  # ├ Obs Thermal Acclimation ----
  
  ## T_opt
  tmp <-
    plot_obs_vs_pred(
      df_in = df_plot, 
      x = "temp", 
      x_se = NULL,
      y ="eva_data_tc_opt", 
      y_se = "eva_data_tc_opt_se", 
      fill = "cl", 
      shape = NULL,
      yxlim = 40,
      add_lm_metrics = F)
  
  p_topt_obs_tcair <-
    tmp$plot +
    xlab(bquote(T[growth] ~ " [°C]")) +
    ylab(bquote("Observed " ~ T[opt] ~ " [°C]")) +
    ggtitle(bquote("Observed Thermal Acclimation Potential of" ~ T[opt]), NULL)
  
  metrics_obstemp_topt <- tmp$df_metrics
  
  ## A_opt
  tmp <-
    plot_obs_vs_pred(
      df_in = df_plot, 
      x = "temp", 
      x_se = NULL,
      y ="eva_data_anet_opt", 
      y_se = "eva_data_anet_opt_se", 
      fill = "cl", 
      shape = NULL,
      plot11 = F,
      xlim = 40,
      ylim = 30,
      add_lm_metrics = F)
  
  p_aopt_obs_tcair <- 
    tmp$plot +
    xlab(bquote(T[growth] ~ " [°C]")) +
    ylab(bquote("Observed " ~ A[opt] ~ "[µmol" ~ m^-2 ~ s^-1 ~ "]")) +
    ggtitle(bquote("Observed Thermal Acclimation Potential of" ~ A[opt]), NULL)
  
  metrics_obstemp_aopt <- tmp$df_metrics
  
  ## A_growth
  tmp <-
    plot_obs_vs_pred(
      df_in = df_plot, 
      x = "temp", 
      x_se = NULL,
      y = "eva_data_agrowth", 
      y_se = NULL, 
      fill = "cl", 
      shape = NULL,
      plot11 = F,
      xlim = 40,
      ylim = 30,
      add_lm_metrics = F)
  
  p_agrowth_obs_tcair <- 
    tmp$plot +
    xlab(bquote(T[growth] ~ " [°C]")) +
    ylab(bquote("Observed " ~ A[growth] ~ "[µmol m¯² s¯¹ ]")) +
    ggtitle(bquote("Observed Thermal Acclimation Potential of" ~ A[growth]), NULL)
  
  metrics_obstemp_agrowth <- tmp$df_metrics
  
  ## T_span
  tmp <-
    plot_obs_vs_pred(
      df_in = df_plot, 
      x = "temp", 
      x_se = NULL,
      y = "eva_data_tspan", 
      y_se = NULL, 
      fill = "cl", 
      shape = NULL,
      plot11 = F,
      xlim = 40,
      ylim = 30,
      add_lm_metrics = F)
  
  p_tspan_obs_tcair <- 
    tmp$plot +
    xlab(bquote(T[growth] ~ " [°C]")) +
    ylab(bquote("Observed " ~ T[span] ~ "[°C]")) +
    ggtitle(bquote("Observed Thermal Acclimation Potential of" ~ T[span]), NULL)
  
  metrics_obstemp_tspan <- tmp$df_metrics
  
  #_____________________________________________________________________________
  # ├ Mod Thermal Acclimation ----
  ## T_opt
  tmp <-
    plot_obs_vs_pred(
      df_in = df_plot, 
      x = "temp", 
      x_se = NULL,
      y ="rpm_sim_tc_opt", 
      y_se = NULL, 
      fill = "cl", 
      shape = NULL,
      yxlim = 40,
      add_lm_metrics = F)
  
  p_topt_mod_tcair <- 
    tmp$plot +
    xlab(bquote(T[growth] ~ " [°C]")) +
    ylab(bquote("Modelled " ~ T[opt] ~ " [°C]")) +
    ggtitle(bquote("Modelled Thermal Acclimation Potential of" ~ T[opt]), NULL)
  
  metrics_modtemp_topt <- tmp$df_metrics
  
  ## A_opt
  tmp <-
    plot_obs_vs_pred(
      df_in = df_plot, 
      x = "temp", 
      x_se = NULL,
      y ="rpm_sim_anet_opt",  
      y_se = NULL, 
      fill = "cl", 
      shape = NULL,
      plot11 = F, 
      xlim = 40, 
      ylim = 30,
      add_lm_metrics = F)
  
  p_aopt_mod_tcair <- 
    tmp$plot +
    xlab(bquote(T[growth] ~ " [°C]")) +
    ylab(bquote("Modelled " ~ A[opt] ~ "[µmol" ~ m^-2 ~ s^-1 ~ "]")) +
    ggtitle(bquote("Modelled Thermal Acclimation Potential of" ~ A[opt]), NULL)
  
  metrics_modtemp_aopt <- tmp$df_metrics
  
  ## A_growth
  tmp <-
    plot_obs_vs_pred(
      df_in = df_plot, 
      x = "temp", 
      x_se = NULL,
      y = "rpm_sim_anet_growth", 
      y_se = NULL, 
      fill = "cl", 
      shape = NULL,
      plot11 = F,
      xlim = 40,
      ylim = 30,
      add_lm_metrics = F)
  
  p_agrowth_mod_tcair <- 
    tmp$plot +
    xlab(bquote(T[growth] ~ " [°C]")) +
    ylab(bquote("Modelled " ~ A[growth] ~ "[µmol m¯² s¯¹ ]")) +
    ggtitle(bquote("Modelled Thermal Acclimation Potential of" ~ A[growth]), NULL)
  
  metrics_modtemp_agrowth <- tmp$df_metrics
  
  ## T_span
  tmp <-
    plot_obs_vs_pred(
      df_in = df_plot, 
      x = "temp", 
      x_se = NULL,
      y ="rpm_sim_tspan",  
      y_se = NULL, 
      fill = "cl", 
      shape = NULL,
      plot11 = F, 
      xlim = 40, 
      ylim = 30,
      add_lm_metrics = F)
  
  p_tspan_mod_tcair <- 
    tmp$plot +
    xlab(bquote(T[growth] ~ " [°C]")) +
    ylab(bquote("Modelled " ~ T[span] ~ "[°C]")) +
    ggtitle(bquote("Modelled Thermal Acclimation Potential of" ~ T[span]), NULL)
  
  metrics_modtemp_tspan <- tmp$df_metrics
  
  #___________________________________________________________________________
  # ├ T_opt - T_growth ----
  
  #_______________________________________________________________________________
  # Wrangle input df
  df_red_1 <- 
    df_plot %>%
    mutate(pft = purrr::map_chr(data_org, ~pull(., PFT) %>% unique())) %>% 
    unnest(rpm_acc) %>%
    dplyr::select(
      sitename, date, tc_leaf, temp, pft,
      rpm_sim_tc_opt, rpm_sim_anet_opt,
      eva_data_tc_opt, eva_data_tc_opt_se,
      eva_data_anet_opt, eva_data_anet_opt_se,
      eb_convergence, cl, gs_c, vcmax25, jmax25, chi, xi, vcmax, jmax, rd25) %>%
    mutate(model = title_1, eb_convergence = ifelse(is.na(eb_convergence), T, eb_convergence)) %>%
    filter(temp > 0, eb_convergence != F) %>%
    mutate(gs_c = gs_c * 1e6,
           rd25 = rd25 * 1e6,
           vcmax25 = vcmax25*1e6,
           vcmax = vcmax*1e6,
           jmax25  = jmax25 *1e6,
           jmax  = jmax *1e6)
  
  # Get dynamic max and min for gs_c coloring
  gsc_max <- max(df_red_1$gs_c)
  gsc_min <- min(df_red_1$gs_c)
  
  #___________________________________________________________________________
  # LEAF AIR DECOUPLING PLOT
  p_leaf_air1 <-
    df_red_1 %>% 
    ggplot() +
    geom_point(aes(x = temp, y = tc_leaf, fill = gs_c), shape = 21, size = 2.5) +
    geom_smooth(aes(x = temp, y = tc_leaf), color = "red", size = 0.5) +
    scale_fill_continuous(low = "#122339", high = "#e0efff", limits=c(0, gsc_max)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
    xlim(0, 40) +
    ylim(0, 40) +
    xlab(expression("Observed Growth " ~ T[air])) +
    ylab(expression("Modelled Growth " ~ T[leaf])) +
    labs(subtitle = "tc_leaf ~ tc_air",
         fill = bquote(g[s,c] ~ "[µmol C" ~ m^-2 ~ s ^-1 ~ Pa^-1 ~ "]: ")) +
    theme_classic() +
    theme(legend.position = "bottom") 
  
  #___________________________________________________________________________
  # OBSERVED ~ OBSERVED GROWTH TC_AIR
  # T_OPT
  # color <- "deepskyblue3"
  # p_topt_obs_tcair <-
  #   df_red_1 %>% 
  #   ggplot() +
  #   aes(x = temp, y = eva_data_tc_opt) +
  #   geom_errorbar(aes(x = temp, ymin = eva_data_tc_opt - eva_data_tc_opt_se, ymax = eva_data_tc_opt + eva_data_tc_opt_se), color = color) +
  #   geom_point(aes(x = temp, y = eva_data_tc_opt), fill = color, shape = 21, size = 2.5) +
  #   geom_smooth(aes(x = temp, y = eva_data_tc_opt), color = color, method = "lm", fullrange = T) +
  #   ylab(expression("Observed " ~ T[opt])) + 
  #   xlab(expression("Observed Growth " ~ T[air])) + 
  #   xlim(0, 40) + 
  #   ylim(0, 40) +
  #   geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  #   ggtitle(expression("Observed " ~ T[opt] ~ "~" ~ "Observed Growth " ~ T[air])) +
  #   ggpmisc::stat_poly_eq(data = df_red_1,
  #                         formula = y ~ x,
  #                         coef.digits = 2,
  #                         method = "lm",
  #                         aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), 
  #                         parse = TRUE)
  # 
  # # A_OPT
  # color <- "deepskyblue3"
  # p_aopt_obs_tcair <-
  #   df_red_1 %>% 
  #   ggplot() +
  #   aes(x = temp, y = eva_data_anet_opt) +
  #   geom_errorbar(aes(x = temp, ymin = eva_data_anet_opt - eva_data_anet_opt_se, ymax = eva_data_anet_opt + eva_data_anet_opt_se), color = color) +
  #   geom_point(aes(x = temp, y = eva_data_anet_opt), fill = color, shape = 21, size = 2.5) +
  #   geom_smooth(aes(x = temp, y = eva_data_anet_opt), color = color, method = "lm", fullrange = T) +
  #   ylab(bquote("Observed " ~ A[opt] ~ "[µmol" ~ m^-2 ~ s ^-1 ~ "]")) + 
  #   xlab(expression("Observed Growth " ~ T[air] ~ "[°C]")) + 
  #   xlim(0, 40) + 
  #   ylim(0, 40) +
  #   geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  #   ggtitle(expression("Observed " ~ A[opt] ~ "~" ~ "Observed Growth " ~ T[air])) +
  #   ggpmisc::stat_poly_eq(data = df_red_1,
  #                         formula = y ~ x,
  #                         coef.digits = 2,
  #                         method = "lm",
  #                         aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), 
  #                         parse = TRUE)
  
  #___________________________________________________________________________
  # MODELLED ~ OBSERVED GROWTH TC_AIR
  # T_OPT
  # color <- "deepskyblue3"
  # p_topt_mod_tcair <-
  #   df_red_1 %>% 
  #   ggplot() +
  #   aes(x = temp, y = rpm_sim_tc_opt) +
  #   geom_point(aes(x = temp, y = rpm_sim_tc_opt), fill = color, shape = 21, size = 2.5) +
  #   geom_smooth(aes(x = temp, y = rpm_sim_tc_opt), color = color, method = "lm", fullrange = T) +
  #   ylab(expression("Modelled " ~ T[opt])) + 
  #   xlab(expression("Observed Growth " ~ T[air])) + 
  #   xlim(0, 40) + 
  #   ylim(0, 40) +
  #   geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  #   ggtitle(expression("Predicted " ~ T[opt] ~ "~" ~ "Growth " ~ T[air])) +
  #   ggpmisc::stat_poly_eq(data = df_red_1,
  #                         formula = y ~ x,
  #                         coef.digits = 2,
  #                         method = "lm",
  #                         aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), 
  #                         parse = TRUE)
  # 
  # # A_OPT
  # color <- "deepskyblue3"
  # p_aopt_mod_tcair <-
  #   df_red_1 %>% 
  #   ggplot() +
  #   aes(x = temp, y = rpm_sim_anet_opt) +
  #   geom_point(aes(x = temp, y = rpm_sim_anet_opt), fill = color, shape = 21, size = 2.5) +
  #   geom_smooth(aes(x = temp, y = rpm_sim_anet_opt), color = color, method = "lm", fullrange = T) +
  #   ylab(bquote("Modelled" ~ A[opt] ~ "[µmol" ~ m^-2 ~ s ^-1 ~ "]")) + 
  #   xlab(expression("Observed Growth " ~ T[air])) + 
  #   xlim(0, 40) + 
  #   ylim(0, 40) +
  #   geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  #   ggtitle(expression("Predicted " ~ T[opt] ~ "~" ~ "Growth " ~ T[air])) +
  #   ggpmisc::stat_poly_eq(data = df_red_1,
  #                         formula = y ~ x,
  #                         coef.digits = 2,
  #                         method = "lm",
  #                         aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), 
  #                         parse = TRUE)
  
  # __________________________________________________________________________
  # OBSERVED TOP_T  ~ PREDICTED TC_LEAF GROWTH
  color <- "darkgreen"
  p_topt_obs_tcleaf <-
    df_red_1 %>% 
    ggplot() +
    aes(x = tc_leaf, y = eva_data_tc_opt) +
    geom_errorbar(aes(x = tc_leaf, ymin = eva_data_tc_opt - eva_data_tc_opt_se, ymax = eva_data_tc_opt + eva_data_tc_opt_se), color = color) +
    geom_point(aes(x = tc_leaf, y = eva_data_tc_opt), fill = color, shape = 21, size = 2.5) +
    geom_smooth(aes(x = tc_leaf, y = eva_data_tc_opt), color = color, method = "lm", fullrange = T) +
    ylab(expression("Observed " ~ T[opt])) + 
    xlab(expression("Modelled Growth " ~ T[leaf])) + 
    xlim(0, 40) + 
    ylim(0, 40) +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
    ggtitle(expression("Modelled " ~ T[opt] ~ "~" ~ "Modelled Growth " ~ T[leaf])) +
    ggpmisc::stat_poly_eq(data = df_red_1,
                          formula = y ~ x,
                          coef.digits = 2,
                          method = "lm",
                          aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), 
                          parse = TRUE) +
    theme_classic()
  
  
  color <- "skyblue"
  p_topt_obs_tcair_monocol <-
    df_red_1 %>% 
    ggplot() +
    aes(x = temp, y = eva_data_tc_opt) +
    geom_errorbar(aes(x = temp, ymin = eva_data_tc_opt - eva_data_tc_opt_se, ymax = eva_data_tc_opt + eva_data_tc_opt_se), color = color) +
    geom_point(aes(x = temp, y = eva_data_tc_opt), fill = color, shape = 21, size = 2.5) +
    geom_smooth(aes(x = temp, y = eva_data_tc_opt), color = color, method = "lm", fullrange = T) +
    ylab(expression("Observed " ~ T[opt])) + 
    xlab(expression("Observed Growth " ~ T[air])) + 
    xlim(0, 40) + 
    ylim(0, 40) +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
    ggtitle(expression("Observed " ~ T[opt] ~ "~" ~ "Observed Growth " ~ T[air])) +
    ggpmisc::stat_poly_eq(data = df_red_1,
                          formula = y ~ x,
                          coef.digits = 2,
                          method = "lm",
                          aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), 
                          parse = TRUE) +
    theme_classic()
  
  # __________________________________________________________________________
  # ├ Traits ~ T_growth ----
  list_traits <- list()
  
  for (trait in c("chi", "vcmax25", "jmax25", "rd25", "vcmax", "jmax", "gs_c", "xi")) {
    # for (trait in c("vcmax25", "jmax25")) {
    xlim_val <- max(df_red_1[["temp"]]) + 5
    df_red_1$y <- df_red_1[[trait]]
    
    if (trait == "vcmax25") {
      ylab_txt <- bquote(V[cmax]^25 ~ "[µmol" ~ m^-2 ~ s ^-1 ~ "]")
      ylim_val <- 150
    } else if (trait == "rd25") {
      ylab_txt <- bquote(R[d]^25 ~ "[µmol" ~ m^-2 ~ s ^-1 ~ "]")
      ylim_val <- 150*0.015
    } else if (trait == "vcmax") {
      ylab_txt <- bquote(V[cmax]^acc ~ "[µmol" ~ m^-2 ~ s ^-1 ~ "]")
      ylim_val <- 150
    } else if (trait == "jmax25") {
      ylab_txt <-  bquote(J[max]^25 ~ "[µmol" ~ m^-2 ~ s ^-1 ~ "]")
      ylim_val <- 300
    } else if (trait == "jmax") {
      ylab_txt <-  bquote(J[max]^acc ~ "[µmol" ~ m^-2 ~ s ^-1 ~ "]")
      ylim_val <- 300
    } else if (trait == "chi") {
      ylab_txt     <- bquote(c[i] ~ "/" ~ c[a] ~ "[-]")
      ylim_val     <- 1
    } else if (trait == "gs_c") {
      ylab_txt    <- bquote(g[s,c] ~ "[µmol C" ~ m^-2 ~ s ^-1 ~ Pa^-1 ~ "]")
      ylim_val <- max(df_red_1$y) + round(max(df_red_1$y) * 0.1)
    } else if (trait == "xi") {
      ylab_txt    <- bquote("xi [" ~ Pa^-0.5 ~ "]")
      # ylim_val <- max(df_red_1$y) + round(max(df_red_1$y) * 0.1)
      ylim_val <- 150
    }
    
    list_traits[[trait]] <-
      df_red_1 %>% 
      ggplot() +
      aes(x = temp, y = y) +
      ## Color by PFT
      # geom_point(aes(x = temp, y = y, fill = pft), size = 2, shape = 21) +
      # scale_fill_viridis_d("PFT", option = "viridis", direction = -1) +
      
      ## Color by Climate Zone
      geom_point(aes(x = temp, y = y, fill = cl), color = "black", size = 2, shape = 21) +
      scico::scale_fill_scico_d(palette = "roma") +
      ylab("Trait Value") + 
      xlab(expression(T[growth] ~ "[°C]")) + 
      labs(subtitle=ylab_txt) +
      xlim(0, xlim_val) +
      ylim(0, ylim_val)
    
    if (trait %in% c("chi", "vcmax25", "jmax25", "rd25")) {
      list_traits[[trait]] <- 
        list_traits[[trait]] +
        geom_smooth(aes(x = temp, y = y), method = "lm", fullrange = T, color = "red") +
        ggpmisc::stat_poly_eq(data = df_red_1,
                              formula = y ~ x,
                              coef.digits = 2,
                              method = "lm",
                              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), 
                              parse = TRUE)  
    }
  }
  
  # 1x3
  p_traits <-
    (list_traits[["vcmax25"]]) +
    (list_traits[["rd25"]]) +
    (list_traits[["jmax25"]] + ylab(NULL)) +
    (list_traits[["xi"]] + ylab(NULL)) +
    plot_layout(nrow = 1, guides = "collect") &
    theme(legend.position = "bottom") &
    theme_classic() &
    plot_annotation("Traits along growth temperature gradient")
  
  
  
  # 2x3
  # p_traits <-
  #   (list_traits[["vcmax25"]] + xlab(NULL)) +
  #   (list_traits[["jmax25"]] + xlab(NULL) + ylab(NULL)) +
  #   (list_traits[["chi"]] + xlab(NULL) + ylab(NULL)) +
  #   list_traits[["vcmax"]] +
  #   (list_traits[["jmax"]] + ylab(NULL)) +
  #   (list_traits[["xi"]] + ylab(NULL)) +
  #   # list_traits[["gs_c"]] + ylab(NULL) +
  #   plot_layout(nrow = 2) &
  #   theme_classic() &
  #   plot_annotation("Traits along growth temperature gradient")
  # 
  
  # __________________________________________________________________________
  # FINAL PLOTS
  
  p_all1 <- p_modobs_topt + p_modobs_aopt + p_leaf_air1 + p_topt_obs_tcleaf + plot_annotation(title = title_1)
  
  out <- list(
    ## Summary Plots
    all = p_all1,
    p_traits = p_traits,
    p_traits_list = list_traits,
    
    ## Air-Leaf Plots
    leaf_air = p_leaf_air1,
    topt_obs_tcair_monocol = p_topt_obs_tcair_monocol,
    
    ## T_opt Plots
    modobs_topt = p_modobs_topt,
    p_topt_mod_tcair   = p_topt_mod_tcair,
    p_topt_obs_tcair   = p_topt_obs_tcair,
    p_topt_obs_tcleaf  = p_topt_obs_tcleaf,
    
    ## A_OPT PLOTS
    modobs_aopt = p_modobs_aopt,
    p_aopt_obs_tcair = p_aopt_obs_tcair,
    p_aopt_mod_tcair = p_aopt_mod_tcair,
    
    ## A_growth PLOTS
    modobs_agrowth      = p_modobs_agrowth,
    p_agrowth_obs_tcair = p_agrowth_obs_tcair,
    p_agrowth_mod_tcair = p_agrowth_mod_tcair,
    
    ## T_span
    modobs_tspan       = p_modobs_tspan,
    p_tspan_obs_tcair    = p_tspan_obs_tcair,
    p_tspan_mod_tcair    = p_tspan_mod_tcair,
    
    ## Metrics of modobs
    metrics_modobs_topt = metrics_modobs_topt,
    metrics_modobs_aopt = metrics_modobs_aopt,
    metrics_modobs_agrowth = metrics_modobs_agrowth,
    metrics_modobs_tspan = metrics_modobs_tspan,
  
    ## Observed thermal acclimation capacity
    metrics_obstemp_topt = metrics_obstemp_topt,
    metrics_obstemp_aopt = metrics_obstemp_aopt,
    metrics_obstemp_agrowth = metrics_obstemp_agrowth,
    metrics_obstemp_tspan = metrics_obstemp_tspan,
    
    ## Modelled thermal acclimation capacity
    metrics_modtemp_topt = metrics_modtemp_topt,
    metrics_modtemp_aopt = metrics_modtemp_aopt,
    metrics_modtemp_agrowth = metrics_modtemp_agrowth,
    metrics_modtemp_tspan = metrics_modtemp_tspan)
  
  return(out)
}


## Analysis-specific ----
### Isolated acclimation processes ----
plot_iso_acc_modobs <- function(trait,
                                p_noacc_scaled,
                                p_fullacc,
                                p_k19,
                                p_xi,
                                p_vj_phi,
                                p_k19_vj_phi,
                                p_xi_k19,
                                p_xi_vj_phi,
                                dir_figs,
                                modobs_or_modtcair) {
  
  if (modobs_or_modtcair == "modobs") tmp <- paste0("modobs_", trait)
  if (modobs_or_modtcair == "modtcair") tmp <- paste0("p_", trait, "_mod_tcair")
  obs_tcair <- paste0("p_", trait, "_obs_tcair")
  
  # ## Full vs no acclimation
  # p_fullno <-
  #   p_noacc_scaled[[tmp]] + labs(title = NULL, caption = NULL, subtitle = "No Acclimation") +
  #   p_fullacc[[tmp]] + labs(title = NULL, caption = NULL, subtitle = "Full Acclimation") +
  #   plot_layout(guides = "collect", nrow = 1) &
  #   theme(legend.position = "bottom")
  # 
  # ## 1/3 acclimated processes
  # p_single <-
  #   p_k19[[tmp]]     + labs(title = NULL, caption = NULL, subtitle = "Acclimated Thermal Response (K19)") +
  #   p_xi[[tmp]]      + labs(title = NULL, caption = NULL, subtitle = "Acclimated Stomatal Conductance (Xi)") +
  #   p_vj_phi[[tmp]]  + labs(title = NULL, caption = NULL, subtitle = "Acclimated Baserate (Vcmax25, Jmax25)") +
  #   plot_layout(guides = "collect", nrow = 1) &
  #   # plot_annotation(title = "Individual Processes") &
  #   theme(legend.position = "bottom")
  # 
  # ## 2/3 acclimated processes
  # p_comb <-
  #   p_k19_vj_phi[[tmp]] + labs(title = NULL, caption = NULL, subtitle = "Base rate + response") +
  #   p_xi_k19[[tmp]]     + labs(title = NULL, caption = NULL, subtitle = "Response + xi") +
  #   p_xi_vj_phi[[tmp]]  + labs(title = NULL, caption = NULL, subtitle = "Base rate + xi") +
  #   plot_layout(guides = "collect", nrow = 1) &
  #   # plot_annotation(title = "Combination of Processes") &
  #   theme(legend.position = "bottom")
  # 
  # ## 3/3 acclimated processes
  # p_matrix <-
  #   (p_vj_phi[[tmp]] + labs(subtitle = "base rate only", caption = NULL, title = NULL) + theme(legend.position = "none") +
  #      p_xi_vj_phi[[tmp]] + labs(subtitle = "base rate + stomatal sensitivity", caption = NULL, title = NULL, x = NULL, y = NULL) + theme(legend.position = "none") +
  #      p_k19_vj_phi[[tmp]] + labs(subtitle = "base rate + thermal response", caption = NULL, title = NULL, x = NULL, y = NULL) + theme(legend.position = "none")
  #   ) /
  #   (plot_spacer() +
  #      p_xi[[tmp]] + labs(subtitle = "stomatal sensitivity only", caption = NULL, title = NULL) + theme(legend.position = "none") +
  #      p_xi_k19[[tmp]] + labs(subtitle = "stomatal sensitivity + thermal response", caption = NULL, title = NULL, x = NULL, y = NULL) + theme(legend.position = "none")
  #   ) /
  #   (plot_spacer() +
  #      plot_spacer() +
  #      p_k19$modobs_topt + labs(subtitle = "thermal response only", caption = NULL, title = NULL) + theme(legend.position = "none")
  #   ) +
  #   plot_layout(guides = "collect", nrow = 3) &
  #   plot_annotation(title = "Effect of individual acclimation processes") &
  #   theme(legend.position = "bottom")
  
  ## Non, partial, full acclimated processes
  

  
  if (modobs_or_modtcair == "modobs") {
    if (trait == "topt")  {
      y_axis <- expression("Observed" ~ T[opt] ~ "[°C]")
      big_title <- expression(bold("Prediction of " ~ T[opt] ~ "Across Model Setups"))
    }
    if (trait == "aopt")  {
      y_axis <- expression("Observed" ~ A[opt] ~ "[µmol" ~ m ^-2 ~ s ^-1 ~ "]")
      big_title <- expression(bold("Prediction of " ~ A[opt] ~ "Across Model Setups"))
    }
    if (trait == "tspan") {
      y_axis <- expression("Observed" ~ T[span] ~ "[°C]")
      big_title <- expression(bold("Prediction of " ~ T[span] ~ "Across Model Setups"))
    }
  }
  
  if (modobs_or_modtcair == "modtcair") {
    if (trait == "topt")  {
      y_axis <- expression(T[opt] ~ "[°C]")
      big_title <- expression(bold("Thermal Acclimation of " ~ T[opt] ~ " Across Model Setups"))
    }
    if (trait == "aopt")  {
      y_axis <- expression(A[opt] ~ "[µmol" ~ m ^-2 ~ s ^-1 ~ "]")
      big_title <- expression(bold("Thermal Acclimation of " ~ A[opt] ~ " Across Model Setups"))
    }
    if (trait == "tspan") {
      y_axis <- expression(T[span] ~ "[°C]")
      big_title <- expression(bold("Thermal Acclimation of " ~ T[span] ~ " Across Model Setups"))
    }
  }
  
  if (modobs_or_modtcair == "modobs") {
    p_all <- 
      ## Top row
      p_noacc_scaled[[tmp]] + labs(x = NULL, y = y_axis, title = NULL, caption = NULL, subtitle = "No Acclimation") +
      p_fullacc[[tmp]]      + labs(x = NULL, y = NULL, title = NULL, caption = NULL, subtitle = "Full Acclimation") +
      plot_spacer() +
      
      ## Mid Row
      p_vj_phi[[tmp]]       + labs(x = NULL, y = y_axis, title = NULL, caption = NULL, subtitle = "PC") +
      p_xi[[tmp]]           + labs(x = NULL, y = NULL, title = NULL, caption = NULL, subtitle = "SS") +
      p_k19[[tmp]]          + labs(x = NULL, y = NULL, title = NULL, caption = NULL, subtitle = "ER") +
      
      ## Bottom Row
      p_k19_vj_phi[[tmp]]   + labs(title = NULL, y = y_axis, caption = NULL, subtitle = "PC + ER") +
      p_xi_vj_phi[[tmp]]    + labs(title = NULL, y = NULL, caption = NULL, subtitle = "PC + SS") +
      p_xi_k19[[tmp]]       + labs(title = NULL, y = NULL, caption = NULL, subtitle = "ER + SS") +
      
      ## Aesthetics
      plot_layout(guides = "collect", nrow = 3) +
      plot_annotation(title = big_title) &
      theme_classic() &
      theme(legend.position = "bottom",
            plot.title = element_text(size = 16, face = "bold", hjust = 0),
            plot.subtitle = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 12))
  }
  
  if (modobs_or_modtcair == "modtcair") {
    p_all <- 
      ## Top row
      p_noacc_scaled[[tmp]] + labs(x = NULL, y = y_axis, title = NULL, caption = NULL, subtitle = "No Acclimation") +
      p_fullacc[[tmp]]      + labs(x = NULL, y = NULL, title = NULL, caption = NULL, subtitle = "Full Acclimation") +
      p_fullacc[[obs_tcair]] + labs(x = NULL, y = NULL, title = NULL, caption = NULL, subtitle = "Observations") +
      # plot_spacer() +
      
      ## Mid Row
      p_vj_phi[[tmp]] + labs(x = NULL, y = y_axis, title = NULL, caption = NULL, subtitle = "PC") +
      p_xi[[tmp]]     + labs(x = NULL, y = NULL, title = NULL, caption = NULL, subtitle = "SS") +
      p_k19[[tmp]]    + labs(x = NULL, y = NULL, title = NULL, caption = NULL, subtitle = "ER") +
      
      ## Bottom Row
      p_k19_vj_phi[[tmp]] + labs(title = NULL, y = y_axis, caption = NULL, subtitle = "PC + ER") +
      p_xi_vj_phi[[tmp]]  + labs(title = NULL, y = NULL, caption = NULL, subtitle = "PC + SS") +
      p_xi_k19[[tmp]]     + labs(title = NULL, y = NULL, caption = NULL, subtitle = "ER + SS") +
      
      ## Aesthetics
      plot_layout(guides = "collect", nrow = 3) +
      plot_annotation(title = big_title) &
      theme_classic() &
      theme(legend.position = "bottom",
            plot.title = element_text(size = 16, face = "bold", hjust = 0),
            plot.subtitle = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 12))
  }
  
  ## Save
  path <- paste0(dir_figs, "/", modobs_or_modtcair, "_", trait)
  
  # ggsave(paste0(path, "_full-vs-no_acclimation.pdf"),
  #        p_fullno,
  #        height = 5,
  #        width = 10)
  # 
  # ggsave(paste0(path, "_1of3_effects.pdf"),
  #        p_single,
  #        height = 5,
  #        width = 15)
  # 
  # ggsave(paste0(path, "_2of3_effects.pdf"),
  #        p_comb,
  #        height = 5,
  #        width = 15)
  # 
  # ggsave(paste0(path, "_3of3_matrix_effects.pdf"),
  #        p_matrix,
  #        height = 12.5,
  #        width = 12.5)
  
  ggsave(paste0(path, "_all_effects.pdf"),
         p_all,
         height = 10,
         width = 10)
  
  return(p_all)
}

plot_capacities <- function(setup,
                            p_topt,
                            p_aopt,
                            p_tspan,
                            dir_figs) {
  
  p <- 
    (p_topt  + labs(subtitle = NULL, title = "Modelled", caption = NULL, x = NULL, y = bquote(T[opt] ~ " [°C]")) + theme(plot.title = element_text(hjust = 0))) + 
    p_aopt   + labs(subtitle = NULL,  title = NULL, caption = NULL, x = NULL, y = bquote(A[opt] ~ "[µmol" ~ m ^-2 ~ s ^-1 ~ "]")) + 
    p_tspan  + labs(subtitle = NULL,  title = NULL, caption = NULL, y = bquote(T[span] ~ " [°C]")) + 
    plot_layout(guides = "collect", 
                ncol = 1) &
    plot_annotation(title = setup) &
    theme(legend.position = "bottom") 
  
  ggsave(paste0(dir_figs, "thermal_acclimation_capacities", "/", setup, ".pdf"),
         p,
         height = 10,
         width = 5)
}

plot_iso_acc_capa <- function(p_noacc_scaled,
                              p_fullacc,
                              p_k19,
                              p_xi,
                              p_vj_phi,
                              p_k19_vj_phi,
                              p_xi_k19,
                              p_xi_vj_phi,
                              dir_figs) {
  
  plot_capacities(
    "noacc",
    p_noacc_scaled$p_topt_mod_tcair,
    p_noacc_scaled$p_aopt_mod_tcair,
    p_noacc_scaled$p_tspan_mod_tcair,
    dir_figs)
  
  plot_capacities(
    "noacc_scaled",
    p_noacc_scaled$p_topt_mod_tcair,
    p_noacc_scaled$p_aopt_mod_tcair,
    p_noacc_scaled$p_tspan_mod_tcair,
    dir_figs)
  
  plot_capacities(
    "fullacc",
    p_fullacc$p_topt_mod_tcair,
    p_fullacc$p_aopt_mod_tcair,
    p_fullacc$p_tspan_mod_tcair,
    dir_figs)
  
  plot_capacities(
    "k19",
    p_k19$p_topt_mod_tcair,
    p_k19$p_aopt_mod_tcair,
    p_k19$p_tspan_mod_tcair,
    dir_figs)
  
  plot_capacities(
    "xi",
    p_xi$p_topt_mod_tcair,
    p_xi$p_aopt_mod_tcair,
    p_xi$p_tspan_mod_tcair,
    dir_figs)
  
  plot_capacities(
    "vj_phi",
    p_vj_phi$p_topt_mod_tcair,
    p_vj_phi$p_aopt_mod_tcair,
    p_vj_phi$p_tspan_mod_tcair,
    dir_figs)
  
  plot_capacities(
    "k19_vj_phi",
    p_k19_vj_phi$p_topt_mod_tcair,
    p_k19_vj_phi$p_aopt_mod_tcair,
    p_k19_vj_phi$p_tspan_mod_tcair,
    dir_figs)
  
  plot_capacities(
    "xi_k19",
    p_xi_k19$p_topt_mod_tcair,
    p_xi_k19$p_aopt_mod_tcair,
    p_xi_k19$p_tspan_mod_tcair,
    dir_figs)
  
  plot_capacities(
    "xi_vj_phi",
    p_xi_vj_phi$p_topt_mod_tcair,
    p_xi_vj_phi$p_aopt_mod_tcair,
    p_xi_vj_phi$p_tspan_mod_tcair,
    dir_figs)
  
}

table_acclimation_contributions <- function(trait,
                                            p_noacc_scaled,
                                            p_fullacc,
                                            p_k19,
                                            p_xi,
                                            p_vj_phi,
                                            p_k19_vj_phi,
                                            p_xi_k19,
                                            p_xi_vj_phi,
                                            dir_tabs,
                                            comparison){
  if (comparison == "trait-temp") {
    tmp <- paste0("metrics_modtemp_", trait)
  } else if (comparison == "mod-obs") {
    tmp <- paste0("metrics_modobs_", trait)
  } else {
    stop("Error in selected comparison")
  }
  
  p_thresh <- 0.01
  
  p_noacc_scaled[[tmp]]$setup <- "p_noacc_scaled"
  p_fullacc[[tmp]]$setup <- "p_fullacc"
  p_k19[[tmp]]$setup <- "p_er"
  p_xi[[tmp]]$setup <- "p_sb"
  p_vj_phi[[tmp]]$setup <- "p_pc"
  p_k19_vj_phi[[tmp]]$setup <- "p_er_pc"
  p_xi_k19[[tmp]]$setup <- "p_er_sb"
  p_xi_vj_phi[[tmp]]$setup <- "p_sb_pc"
  
  p_noacc_scaled[[paste0("metrics_obstemp_", trait)]]$setup <- "obs"
  
  tab_acc_capa <- 
    rbind(
      p_noacc_scaled[[paste0("metrics_obstemp_", trait)]],
      p_fullacc[[tmp]],
      p_noacc_scaled[[tmp]],
      p_vj_phi[[tmp]],
      p_xi[[tmp]],
      p_k19[[tmp]],
      p_xi_vj_phi[[tmp]],
      p_k19_vj_phi[[tmp]],
      p_xi_k19[[tmp]]) %>% 
    relocate(setup)
  
  
  
  txt_top      <- paste0("setup, intercept, slope vs T_growth (* = sign. at 0.05), ip_xy, n_ac / n_aj")
  txt_obs      <- paste0(p_noacc_scaled[[paste0("metrics_obstemp_", trait)]]$setup, ", ", p_noacc_scaled[[paste0("metrics_obstemp_", trait)]]$b0, ", ", p_noacc_scaled[[paste0("metrics_obstemp_", trait)]]$b1, ifelse(p_noacc_scaled[[paste0("metrics_obstemp_", trait)]]$pval > p_thresh, "", "*"), ", ", p_noacc_scaled[[paste0("metrics_obstemp_", trait)]]$ip_xy, ", NA")
  txt_noacc    <- paste0(p_noacc_scaled[[tmp]]$setup, ", ", p_noacc_scaled[[tmp]]$b0, ", ", p_noacc_scaled[[tmp]]$b1, ifelse(p_noacc_scaled[[tmp]]$pval > p_thresh, "", "*"), ", ", p_noacc_scaled[[tmp]]$ip_xy, ", ", p_noacc_scaled[[tmp]]$n_ac, " / ", p_noacc_scaled[[tmp]]$n_aj)
  txt_fullacc  <- paste0(p_fullacc[[tmp]]$setup, ", ", p_fullacc[[tmp]]$b0, ", ", p_fullacc[[tmp]]$b1, ifelse(p_fullacc[[tmp]]$pval > p_thresh, "", "*"), ", ", p_fullacc[[tmp]]$ip_xy, ", ", p_fullacc[[tmp]]$n_ac, " / ", p_fullacc[[tmp]]$n_aj)
  txt_er       <- paste0(p_k19[[tmp]]$setup, ", ", p_k19[[tmp]]$b0, ", ", p_k19[[tmp]]$b1, ifelse(p_k19[[tmp]]$pval > p_thresh, "", "*"), ", ", p_k19[[tmp]]$ip_xy, ", ", p_k19[[tmp]]$n_ac, " / ", p_k19[[tmp]]$n_aj)
  txt_sb       <- paste0(p_xi[[tmp]]$setup, ", ", p_xi[[tmp]]$b0, ", ", p_xi[[tmp]]$b1, ifelse(p_xi[[tmp]]$pval > p_thresh, "", "*"), ", ", p_xi[[tmp]]$ip_xy, ", ", p_xi[[tmp]]$n_ac, " / ", p_xi[[tmp]]$n_aj)
  txt_pc       <- paste0(p_vj_phi[[tmp]]$setup, ", ", p_vj_phi[[tmp]]$b0, ", ", p_vj_phi[[tmp]]$b1, ifelse(p_vj_phi[[tmp]]$pval > p_thresh, "", "*"), ", ", p_vj_phi[[tmp]]$ip_xy, ", ", p_vj_phi[[tmp]]$n_ac, " / ", p_vj_phi[[tmp]]$n_aj)
  txt_er_pc    <- paste0(p_k19_vj_phi[[tmp]]$setup, ", ", p_k19_vj_phi[[tmp]]$b0, ", ", p_k19_vj_phi[[tmp]]$b1, ifelse(p_k19_vj_phi[[tmp]]$pval > p_thresh, "", "*"), ", ", p_k19_vj_phi[[tmp]]$ip_xy, ", ", p_k19_vj_phi[[tmp]]$n_ac, " / ", p_k19_vj_phi[[tmp]]$n_aj)
  txt_sb_pc    <- paste0(p_xi_k19[[tmp]]$setup, ", ", p_xi_k19[[tmp]]$b0, ", ", p_xi_k19[[tmp]]$b1, ifelse(p_xi_k19[[tmp]]$pval > p_thresh, "", "*"), ", ", p_xi_k19[[tmp]]$ip_xy, ", ", p_xi_k19[[tmp]]$n_ac, " / ", p_xi_k19[[tmp]]$n_aj)
  txt_er_sb    <- paste0(p_xi_vj_phi[[tmp]]$setup, ", ", p_xi_vj_phi[[tmp]]$b0, ", ", p_xi_vj_phi[[tmp]]$b1, ifelse(p_xi_vj_phi[[tmp]]$pval > p_thresh, "", "*"), ", ", p_xi_vj_phi[[tmp]]$ip_xy, ", ", p_xi_vj_phi[[tmp]]$n_ac, " / ", p_xi_vj_phi[[tmp]]$n_aj)

  tab_red <- 
    list(txt_top,
           txt_obs,
           txt_fullacc,
           txt_noacc,
           txt_pc,
           txt_sb,
           txt_er,
           txt_sb_pc,
           txt_er_pc,
           txt_er_sb)
  
  # Cleaning to avoid confusion of mod-obs and topt-tgrowth 
  if (comparison == "mod-obs") {
    tab_acc_capa <- tab_acc_capa %>% dplyr::filter(!str_detect(setup, "obs"))
    tab_red      <- tab_red[-2]
  }
    
  write_csv(tab_acc_capa, paste0(here(dir_tabs), "/", comparison, "_full_", trait, ".csv")) 
  write_lines(tab_red, paste0(here(dir_tabs), "/", comparison, "_clean_", trait, ".csv"))
  
  return(list(tab_red, tab_acc_capa))
}

### Seasonality of Traits ----
run_inst_for_365_days <- function(df_acc, 
                                  settings, 
                                  setup,
                                  dir_mods,
                                  set_add,
                                  step_size_tcair   = 0.5) {
  
  df_tmp  <- tibble()
  
  ## Take average ppfd across all measurements (is ca. 1500umol/m2/s which is light-saturated)
  ppfd_data <- 
    df_acc %>% 
    unnest(meas_cond) %>% 
    pull(PARi) %>% 
    mean(., na.rm = T) *
    1e-6
  
  if (set_add$high_light) ppfd_data <- 2000e-6
  
  for (i in 1:nrow(df_acc)) {
    # print(i)
    
    df_i <-  tibble(tc_leaf  = seq(step_size_tcair, 40, step_size_tcair), 
                    tc_air   = NA)
    
    df_i$sitename <- df_acc$sitename[[i]]
    df_i$date     <- df_acc$date[[i]]
    
    df_i$patm     <- df_acc$forcing_d[[i]]$patm
    df_i$co2      <- df_acc$forcing_d[[i]]$co2
    df_i$tc_air   <- df_acc$forcing_d[[i]]$temp
    df_i$ppfd     <- ppfd_data
    vpd           <- df_acc$forcing_d[[i]]$vpd
    
    if (set_add$vpd_scaled) {
      df_i$vpd  <- 
        VPDleafToAir(
          vpd / 1000,
          df_i$tc_air,
          df_i$tc_leaf,
          df_i$patm / 1000
        ) * 1000
      
      df_i$vpd <- ifelse(df_i$vpd > 10, df_i$vpd,  10) # If scaled larger than 100, take scaled, else take 100  
    } else {
      df_i$vpd <- vpd
    }
    
    df_i <- df_i %>% nest(forcing_inst = !any_of(c("sitename", "date")))
    df_tmp <- rbind(df_tmp, df_i)
  }
  
  ## Attach to df_acc again
  df_acc <-
    df_acc %>% 
    left_join(df_tmp)
  
  ## Extract relevant variables for instant response
  df_tmp <-
    df_acc %>%
    mutate(settings = list(as_tibble(settings))) %>% 
    mutate(tc_home           = purrr::map_dbl(siteinfo,       ~pull(., tc_home)),
           ppfd_measurement  = purrr::map_dbl(siteinfo,       ~pull(., ppfd_measurement)),
           tc_growth_air  = purrr::map_dbl(forcing_growth, ~pull(., temp)),
           tc_growth_leaf = purrr::map_dbl(rpm_acc,        ~pull(., tc_leaf)),
           kphio          = purrr::map_dbl(rpm_acc,        ~pull(., kphio  )),
           vcmax25        = purrr::map_dbl(rpm_acc,        ~pull(., vcmax25)),
           rd25           = purrr::map_dbl(rpm_acc,        ~pull(., rd25)),
           jmax25         = purrr::map_dbl(rpm_acc,        ~pull(., jmax25)),
           xi             = purrr::map_dbl(rpm_acc,        ~pull(., xi))) %>%
    dplyr::select(-any_of(c("meas_cond", "fit_opt", "rpm_acc", "forcing_growth",
                            "forcing_d", "data_org", "fit_opt", "siteinfo"))) %>% 
    unnest(c("settings", "forcing_inst"))  %>% 
    mutate(id = row_number(), # Have to add a row_id for 1-by-1 feeding of inputs
           method_eb = "off") # To ensure EB is *not* called here
  
  ## Remove sites where vcmax or jmax were predicted to be zero
  df_tmp <-
    df_tmp %>% 
    dplyr::filter(vcmax25 != 0, jmax25 != 0)
  
  ## Run instantaneous response
  df_rpm_inst <-
    df_tmp %>% 
    nest(inputs = !any_of(c("sitename", "date", "id"))) %>%
    rowwise() %>%
    mutate(rpm_inst = rpmodel_inst(inputs) %>%
             as_tibble() %>%
             list()
    ) %>%
    ungroup()
  
  df_inst <-
    df_rpm_inst %>% 
    mutate(
      tc_growth_air  = purrr::map_dbl(inputs, ~pull(., tc_growth_air )),
      tc_growth_leaf = purrr::map_dbl(inputs, ~pull(., tc_growth_leaf)),
      # Round to nearest 0.5 to extract fitting A_growth
      # tc_growth_air = round(tc_growth_air/0.5)*0.5,
      # tc_growth_leaf = round(tc_growth_leaf/0.5)*0.5,
    ) %>%
    # dplyr::select(-inputs, -id) %>% 
    unnest(rpm_inst) %>% 
    nest(rpm_inst = !all_of(c("sitename", "date"))) %>% 
    mutate(
      tc_opt        = purrr::map_dbl(rpm_inst, ~ slice_max(., anet)   %>% pull(tc_leaf)),
      tc_opt_agross = purrr::map_dbl(rpm_inst, ~ slice_max(., agross) %>% pull(tc_leaf)),
      tc_opt_ac     = purrr::map_dbl(rpm_inst, ~ slice_max(., ac)     %>% pull(tc_leaf)),
      tc_opt_aj     = purrr::map_dbl(rpm_inst, ~ slice_max(., aj)     %>% pull(tc_leaf)),
      
      min_a         = purrr::map_chr(rpm_inst, ~ slice_max(., anet)   %>% pull(min_a)),
      min_a         = as.factor(min_a),
      
      agross_opt    = purrr::map_dbl(rpm_inst, ~ slice_max(., agross) %>% pull(agross) * 1e6),
      anet_opt      = purrr::map_dbl(rpm_inst, ~ slice_max(., anet) %>% pull(anet) * 1e6),
      
      anet_growth   = 0, # purrr::map_dbl(rpm_inst, ~ dplyr::filter(., tc_growth_air == tc_leaf) %>% pull(anet) * 1e6),
      min_agrowth   = "ac", # purrr::map_chr(rpm_inst, ~ dplyr::filter(., tc_growth_air == tc_leaf) %>% pull(min_a)),
      min_agrowth   = "ac", # as.factor(min_agrowth),
      
      tspan_l       = purrr::map_dbl(rpm_inst, 
                                     ~ dplyr::filter(., anet > (0.9 * max(anet))) %>% 
                                       slice_min(tc_leaf) %>% 
                                       pull(tc_leaf)),
      tspan_h       = purrr::map_dbl(rpm_inst, 
                                     ~ dplyr::filter(., anet > (0.9 * max(anet))) %>% 
                                       slice_max(tc_leaf) %>% 
                                       pull(tc_leaf)),
      tspan          = tspan_h - tspan_l
    ) %>%
    nest(rpm_sim = !all_of(c("sitename", "date", "rpm_inst"))) %>%
    right_join(df_acc)
  
  tmp_dir <- here(dir_mods, setup)
  if (!dir.exists(here(tmp_dir))) dir.create(tmp_dir, recursive = T, showWarnings = F)
  
  saveRDS(df_inst, paste0(tmp_dir, "/df_inst.rds"))
  return(df_inst)
}


temporal_individual_plots <- function(df_inst,
                                      setup,
                                      min_topt_per_year = 3,
                                      dir_figs) {
  
  # Define variables
  vec_cols <- c(RColorBrewer::brewer.pal(3, "Dark2")[1],
                RColorBrewer::brewer.pal(3, "Dark2")[3],
                RColorBrewer::brewer.pal(3, "Dark2")[2])
  
  facet_cols <- ifelse(min_topt_per_year == 3, 1, 2) 
  p_height      <- ifelse(min_topt_per_year == 3, 8, 12)
  p_width       <- ifelse(min_topt_per_year == 3, 6,12)
  
  # Topt
  p_topt <- 
    df_inst %>% 
    rename(mod = rpm_sim, obs = fit_opt) %>% 
    unnest(c(mod, obs), names_sep = "_") %>% 
    unnest(c(forcing_growth)) %>% 
    pivot_longer(cols = c(mod_tc_opt, temp), names_to = "tc", values_to = "tc_value") %>% 
    ggplot() +
    # theme_linedraw() +
    aes(x = date) +
    geom_line(aes(y = tc_value, color = tc)) +
    geom_errorbar(aes(ymin = obs_tc_opt - obs_tc_opt_se, ymax = obs_tc_opt + obs_tc_opt_se, y = obs_tc_opt, color = "obs_tc_opt"), size = 0.75) +
    geom_point(aes(y = obs_tc_opt, color = "obs_tc_opt"), size = 1.5) +
    facet_wrap(~site_year, ncol = facet_cols , scales = "free") +
    scale_color_manual(
      name = NULL,
      breaks = c("mod_tc_opt", 
                 "temp",
                 "obs_tc_opt"),
      labels = c(bquote("Mod." ~ T[opt]), 
                 bquote(T[growth]),
                 bquote("Obs." ~ T[opt])),
      values = vec_cols,
      guide = guide_legend(
        # label.position = "bottom",
        # frame.colour = "black",
        override.aes = list(shape = c(NA, NA, 16),
                            size = c(1, 1, 1.5),
                            linetype = c(1, 1, 0),
                            # fill = c("black", "grey", "white"),
                            color = vec_cols))) +
    ylim(0, 40) +
    labs(y = bquote("Temperature [°C]"),
         x = "Date")
  
  ggsave(paste0(dir_figs, "/", setup, "_topt_fig.pdf"), p_topt, height = p_height, width = p_width)
  
  p_topt <- 
    df_inst %>% 
    rename(mod = rpm_sim, obs = fit_opt) %>% 
    unnest(c(mod, obs), names_sep = "_") %>% 
    unnest(c(forcing_growth)) %>% 
    pivot_longer(cols = c(mod_tc_opt, temp), names_to = "tc", values_to = "tc_value") %>% 
    ggplot() +
    # theme_linedraw() +
    aes(x = date) +
    geom_line(aes(y = tc_value, color = tc)) +
    geom_errorbar(aes(ymin = obs_tc_opt - obs_tc_opt_se, ymax = obs_tc_opt + obs_tc_opt_se, y = obs_tc_opt, color = "obs_tc_opt"), size = 0.75) +
    geom_point(aes(y = obs_tc_opt, color = "obs_tc_opt"), size = 1.5) +
    facet_wrap(~site_year, ncol = facet_cols , scales = "free") +
    scale_color_manual(
      name = NULL,
      breaks = c("mod_tc_opt", 
                 "temp",
                 "obs_tc_opt"),
      labels = c(bquote("Mod." ~ T[opt]), 
                 bquote(T[growth]),
                 bquote("Obs." ~ T[opt])),
      values = vec_cols,
      guide = guide_legend(
        # label.position = "bottom",
        # frame.colour = "black",
        override.aes = list(shape = c(NA, NA, 16),
                            size = c(1, 1, 1.5),
                            linetype = c(1, 1, 0),
                            # fill = c("black", "grey", "white"),
                            color = vec_cols))) +
    ylim(0, 40) +
    labs(y = bquote("Temperature [°C]"),
         x = "Date") +
    theme_linedraw()
  
  ggsave(paste0(dir_figs, "/", setup, "topt_fig.pdf"), p_topt, height = p_height, width = p_width)
  
  p_aopt <- 
    df_inst %>% 
    rename(mod = rpm_sim, obs = fit_opt) %>% 
    unnest(c(mod, obs), names_sep = "_") %>% 
    unnest(c(forcing_growth)) %>% 
    pivot_longer(cols = c(mod_anet_opt, temp), names_to = "tc", values_to = "tc_value") %>% 
    ggplot() +
    # theme_linedraw() +
    aes(x = date) +
    geom_line(aes(y = tc_value, color = tc)) +
    geom_errorbar(aes(ymin = obs_anet_opt - obs_anet_opt_se, ymax = obs_anet_opt + obs_anet_opt_se, y = obs_anet_opt, color = "obs_anet_opt"), size = 0.75) +
    geom_point(aes(y = obs_anet_opt, color = "obs_anet_opt"), size = 1.5) +
    facet_wrap(~site_year, ncol = facet_cols , scales = "free") +
    scale_color_manual(
      name = NULL,
      breaks = c("mod_anet_opt", 
                 "temp",
                 "obs_anet_opt"),
      labels = c(bquote("Mod." ~ A[opt] ~ "[µmol" ~ m ^-2 ~ s ^-1 ~ "]"), 
                 bquote(T[growth] ~ A[opt] ~ "[°C]"),
                 bquote("Obs." ~ A[opt] ~ "[µmol" ~ m ^-2 ~ s ^-1 ~ "]")),
      values = vec_cols,
      guide = guide_legend(
        override.aes = list(shape = c(NA, NA, 16),
                            size = c(1, 1, 1.5),
                            linetype = c(1, 1, 0),
                            # fill = c("black", "grey", "white"),
                            color = vec_cols))) +
    ylim(0, 40) +
    labs(y = bquote(A[opt] ~ "[µmol" ~ m ^-2 ~ s ^-1 ~ "]"),
         x = "Date") +
    theme_linedraw()
  
  ggsave(paste0(dir_figs, "/", setup, "aopt_fig.pdf"), p_aopt, height = p_height, width = p_width)
  
  p_tspan <- 
    df_inst %>% 
    rename(mod = rpm_sim, obs = fit_opt) %>% 
    unnest(c(mod, obs), names_sep = "_") %>% 
    unnest(c(forcing_growth)) %>% 
    pivot_longer(cols = c(mod_tspan, temp), names_to = "tc", values_to = "tc_value") %>% 
    ggplot() +
    # theme_linedraw() +
    aes(x = date) +
    geom_line(aes(y = tc_value, color = tc)) +
    geom_point(aes(y = obs_tspan, color = "obs_tspan"), size = 1.5) +
    facet_wrap(~site_year, ncol = facet_cols , scales = "free") +
    scale_color_manual(
      name = NULL,
      breaks = c("mod_tspan", 
                 "temp",
                 "obs_tspan"),
      labels = c(bquote("Mod." ~ T[span]), 
                 bquote(T[growth]),
                 bquote("Obs." ~ T[span])),
      values = vec_cols,
      guide = guide_legend(
        # label.position = "bottom",
        # frame.colour = "black",
        override.aes = list(shape = c(NA, NA, 16),
                            size = c(1, 1, 1.5),
                            linetype = c(1, 1, 0),
                            # fill = c("black", "grey", "white"),
                            color = vec_cols))) +
    ylim(0, 40) +
    labs(y = bquote("Temperature [°C]"),
         x = "Date") +
    theme_linedraw()
  
  ggsave(paste0(dir_figs, "/", setup, "aopt_fig.pdf"), p_tspan, height = p_height, width = p_width)
  
  # Old Plots without T_grwoth
  # # Aopt
  # p_aopt <- 
  #   df_inst %>% 
  #   rename(mod = rpm_sim, obs = fit_opt) %>% 
  #   unnest(c(mod, obs), names_sep = "_") %>% 
  #   unnest(c(forcing_growth, rpm_acc)) %>% 
  #   ggplot() +
  #   # theme_linedraw() +
  #   aes(x = date) +
  #   geom_line(aes(y = mod_anet_opt, color = "mod")) +
  #   geom_errorbar(aes(ymin = obs_anet_opt - obs_anet_opt_se, ymax = obs_anet_opt + obs_anet_opt_se, y = obs_anet_opt, color = "obs"), size = 0.75) +
  #   geom_point(aes(y = obs_anet_opt, color = "obs"), size = 1.5) +
  #   facet_wrap(~site_year, ncol = facet_cols , scales = "free") +
  #   scale_color_manual(
  #     name = NULL,
  #     breaks = c("mod",
  #                "obs"),
  #     labels = c(bquote("Mod." ~ A[opt]), 
  #                bquote("Obs." ~ A[opt])),
  #     values = c(vec_cols[1], vec_cols[3])) +
  #   ylim(0, 40) +
  #   labs(y = bquote(A[opt] ~ "[µmol" ~ m ^-2 ~ s ^-1 ~ "]"),
  #        x = "Date")
  # 
  # ggsave(paste0(dir_figs, "/", setup, "aopt_fig.pdf"), p_aopt, height = p_height, width = p_width)
  # 
  # # Tspan
  # p_tspan <- 
  #   df_inst %>% 
  #   rename(mod = rpm_sim, obs = fit_opt) %>% 
  #   unnest(c(mod, obs), names_sep = "_") %>% 
  #   unnest(c(forcing_growth)) %>% 
  #   ggplot() +
  #   # theme_linedraw() +
  #   aes(x = date) +
  #   geom_line(aes(y = mod_tspan, color = "mod")) +
  #   geom_point(aes(y = obs_tspan, color = "obs"), size = 1.5) +
  #   facet_wrap(~site_year, ncol = facet_cols , scales = "free") +
  #   scale_color_manual(
  #     name = NULL,
  #     breaks = c("mod",
  #                "obs"),
  #     labels = c(bquote("Mod." ~ T[span]), 
  #                bquote("Obs." ~ T[span])),
  #     values = c(vec_cols[1], vec_cols[3])) +
  #   ylim(0, 40) +
  #   labs(y = bquote("Temperature [°C]"),
  #        x = "Date")
  # 
  # ggsave(paste0(dir_figs, "/", setup, "tspan_fig.pdf"), p_tspan, height = p_height, width = p_width)
  
  ## All three
  p_all <- 
    (p_topt + p_aopt + p_tspan) +
    plot_layout(ncol = 3) &
    theme(legend.position = "bottom")
  
  ggsave(paste0(dir_figs, "/", setup, "all.pdf"), p_all, height = 14, width = 14)
  
  return(list(
    p_topt = p_topt,
    p_aopt = p_aopt,
    p_tspan = p_tspan,
    p_all = p_all)
  )
}

## Global Maps ----
add_climatezone_to_siteinfo <- function(df_in){
  ## Description:
  # Using parts of official R source code from Kottek et al. (2016) to extract climate zone data
  # Function with input dataframe with "sitename" and nested "siteinfo" that holds lat and lon data
  # Rasterdata has to be saved in "~/data/climate_zones/..."
  
  ## Code
  # Check if required packages are loaded from here: https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/
  packages <- c("rgdal", "raster", "tidyverse")
  package.check <- lapply(packages,
                          FUN = function(x) {if (!require(x, character.only = TRUE)) {
                            install.packages(x, dependencies = TRUE)
                            library(x, character.only = TRUE)}})
  
  # Load raster files
  path <- here("data/climate_zones/")
  period='1986-2010'
  r <- raster(paste(path, 'KG_', period, '.grd', sep=''))
  
  # Color palette for climate classification
  climate.colors=c("#960000", "#FF0000", "#FF6E6E", "#FFCCCC", "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64", "#007800", "#005000", "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00", "#550055", "#820082", "#C800C8", "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF", "#64FFFF", "#F5FFFF")
  
  # Legend must correspond to all climate classes, insert placeholders
  r0 <- r[1:32]; r[1:32] <- seq(1,32,1)
  
  # Converts raster field to categorical data
  r <- ratify(r); rat <- levels(r)[[1]]
  
  # Legend is always drawn in alphabetic order
  rat$climate <- c('Af', 'Am', 'As', 'Aw', 'BSh', 'BSk', 'BWh', 'BWk', 'Cfa', 'Cfb','Cfc', 'Csa', 'Csb', 'Csc', 'Cwa','Cwb', 'Cwc', 'Dfa', 'Dfb', 'Dfc','Dfd', 'Dsa', 'Dsb', 'Dsc', 'Dsd','Dwa', 'Dwb', 'Dwc', 'Dwd', 'EF','ET', 'Ocean')
  
  # Remove the placeholders
  r[1:32] <- r0; levels(r) <- rat
  
  # Find climate zone for each sample in siteinfo
  siteinfo    <- df_in %>% dplyr::select(sitename, siteinfo) %>% unnest(siteinfo)
  siteinfo$cl <- NA
  
  for (i in 1:nrow(siteinfo)) {
    # print(i)
    lon <- siteinfo$lon[i]
    lat <- siteinfo$lat[i]
    cl_id <- r[cellFromXY(r, c(lon, lat))]
    siteinfo$cl[i] <- rat$climate[which(rat$ID == cl_id)]
  }
  
  # Return df_in with updated siteinfo
  df_out <-
    siteinfo %>%
    nest(siteinfo = !all_of(c("sitename"))) %>% 
    left_join(df_in %>% dplyr::select(-siteinfo))
  
  return(df_out)
}

## .................................................................................................

get_kg_climate_plot <- function(){
  packages <- c("rgdal", "raster", "tidyverse")
  package.check <- lapply(packages,
                          FUN = function(x) {if (!require(x, character.only = TRUE)) {
                            install.packages(x, dependencies = TRUE)
                            library(x, character.only = TRUE)}})
  
  ## Global Maps Basis
  kg <- readOGR(dsn = here("data/climate_zones/1976-2000_GIS"), layer = "1976-2000")
  kg <- spTransform(kg, CRS("+proj=longlat +datum=WGS84"))
  kg_f <- fortify(kg, region = "GRIDCODE")
  key <- data.frame( id = c(11, 12, 13, 14, 21, 22, 26, 27, 31, 32, 33, 34, 35, 36, 37, 38, 39, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 61, 62),
                     cl = c('Af', 'Am', 'As', 'Aw', 'BSh', 'BSk', 'BWh', 'BWk', 'Cfa', 'Cfb','Cfc', 'Csa', 'Csb', 'Csc', 'Cwa','Cwb', 'Cwc', 'Dfa', 'Dfb', 'Dfc','Dfd', 'Dsa', 'Dsb', 'Dsc', 'Dsd','Dwa', 'Dwb', 'Dwc', 'Dwd', 'EF','ET'))
  kg_final <- merge(kg_f, key)
  
  p_base <-
    ggplot() +
    # theme_classic() +
    geom_polygon(data = kg_final, aes(x = long, y = lat, group = group, fill = cl), alpha = 0.8) +
    # geom_path(data = kg_final, aes(x = long, y = lat, group = group, fill = cl), size = 0.05, color = "black") +
    ylab("Latitude (Decimal Degree)") +
    xlab("Longitude (Decimal Degree)") +
    coord_cartesian(ylim=c(-80, 80)) +
    scale_fill_manual(values = c("#960000", "#FF0000", "#FF6E6E", "#FFCCCC", "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64", "#007800", "#005000", "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00", "#550055", "#820082", "#C800C8", "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF", "#64FFFF", "#F5FFFF")) +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(ncol = 12, direction = "horizontal", title.position = "top")) +
    labs(fill = "Köppen-Geiger Climate Zone")
  
  return(p_base)
} 


plot_global_samples <- function(df_in){
  
  df <-
    df_in %>% 
    unnest(sitedata) %>% 
    mutate(lat = purrr::map_dbl(siteinfo, ~pull(., lat) %>% unique()),
           lon = purrr::map_dbl(siteinfo, ~pull(., lon) %>% unique())) %>% 
    dplyr::select(sitename, agg_date, lat, lon) %>% 
    group_by(sitename, lat, lon) %>% 
    nest() %>% 
    mutate(n = purrr::map_dbl(data, ~nrow(.)))
    
  p_base <- get_kg_climate_plot()
  p_out <-
    p_base +
    geom_point(data = df, aes(x = lon, y = lat, size = n),
               color = "black",
               pch = 21, 
               fill = "white") +    
    ggtitle(paste0("Global Map of Samples")) +
    labs(size = expression(italic("N"))) +
    guides(size = guide_legend(ncol = 2, direction = "horizontal", title.position = "top")) 
  
  return(p_out) 
}

## .................................................................................................
