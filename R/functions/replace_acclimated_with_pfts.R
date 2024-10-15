# Replace acclimated for pft-specific Vcmax25 and Jmax25

# In: df_acc
# Out: df_acc
# Purpose: Replace vcmax25 and jmax25 from df_acc by with pft-specific values

replace_acclimated_with_pfts <- function(df_in,
                                         replace_pc = F,
                                         replace_xi = F,
                                         jvr_method = "kk07",
                                         vcmax25_source = "none") {

  # Check if input dataframe is of right shape
  if (!all(c("sitename", "date", "rpm_acc", "data_org", "fit_opt", "meas_cond",
             "siteinfo", "forcing_d", "forcing_growth") %in% names(df_in))) {
    stop("Input df has unintended variables.")
  }
  
  # Get list of vars to nest by 
  rpm_acc_vars <- names(df_in$rpm_acc[[1]])
  forcing_growth_vars <- names(df_in$forcing_growth[[1]])
  siteinfo_vars <- names(df_in$siteinfo[[1]])
  
  # Get list of pfts
  return_only_g1 <- ifelse(replace_xi & !replace_pc, TRUE, FALSE)
  list_of_pfts <- get_list_of_fixed_pfts(vcmax25_source, return_only_g1 = return_only_g1)
  # print(vcmax25_source)
  # print(list_of_pfts)
  
  tmp <- 
    df_in %>% 
    mutate(pft = purrr::map_chr(data_org, ~pull(., PFT) %>% unique())) %>% 
    unnest(rpm_acc) %>% 
    unnest(c(forcing_growth, siteinfo))
  
  # Replace values where needed
  for (i in 1:nrow(list_of_pfts)) {
    # cat("Working on pft:", list_of_pfts$pft[i] , "\n")
    
    # To set Medlyn and Prentice models equal, they must have the same units
    i_g1 <- list_of_pfts$g1[i] * (1000)^0.5 # g1 is in units kPa ^ 0.5, we need Pa^0.5
    
    i_co2   <- tmp %>% dplyr::filter(pft == list_of_pfts$pft[i]) %>% pull(co2) %>% mean()  # Keep ppm, needed for Medlyn model
    i_patm  <- tmp %>% dplyr::filter(pft == list_of_pfts$pft[i]) %>% pull(patm) %>% mean()
    i_vpd   <- tmp %>% dplyr::filter(pft == list_of_pfts$pft[i]) %>% pull(vpd) %>% mean()
    i_temp  <- tmp %>% dplyr::filter(pft == list_of_pfts$pft[i]) %>% pull(temp) %>% mean()
    i_thome <- tmp %>% dplyr::filter(pft == list_of_pfts$pft[i]) %>% pull(tc_home) %>% unique()
    
    i_ca    <- co2_to_ca(i_co2, i_patm) # In Pa, needed for Prentice model
    i_rstar <- calc_gammastar(i_temp, i_patm) # 
    
    for (j in 1:nrow(tmp)) {
      
      # Verbose Output
      # message("Row ", j, " for source ", vcmax25_source, ": Comparing: PFT ", list_of_pfts$pft[i], "\t to\t site-PFT ", tmp$pft[j])
      
      if (tmp$pft[j] == list_of_pfts$pft[i]) {
        
        # Verbose Output
        # message("Sitename -- PFT list | site -- old vcmax25 -- new vcmax25:\t ",
                # tmp$sitename[j], " -- ", list_of_pfts$pft[i], " | ", tmp$pft[j]," -- ", round(tmp$vcmax25[j]*1e6), " -- ", round(list_of_pfts$vcmax25[i]*1e6))

        # Replacing photosynthetic capacities
        if (replace_pc) {
          
          tmp$kphio[j]   <- 0.04479 
          tmp$vcmax25[j] <- list_of_pfts$vcmax25[i] 
          tmp$rd25[j]    <- list_of_pfts$vcmax25[i] * 0.015
          
          
          # Calculate jmax25
          if (!(jvr_method %in% c("leuning02", "kk07", "kumarathunge19", "kumarathunge19_fixed"))){
            stop("Input for jvr_method method scaling is not valid: ", jvr_method)
            
          } else if (jvr_method %in% c("kk07", "leuning02")) {
            tmp$jmax25[j] <- tmp$vcmax25[j] * 1.88 # Ratio from KK07
            
          } else if (jvr_method == "kumarathunge19_fixed") {
            # Method to remove acclimation effect from K19 parametrization
            
            tc_growth <- get_avg_tc_growth()
            # tc_home   <- get_avg_tc_home() # or i_thome to allow for adaptation effect
            tc_home <- i_thome
            
            tmp$jmax25[j] <- tmp$vcmax25[j] * ( 2.56 - 0.0375 * tc_home - 0.0202 * (tc_growth - tc_home))
            
          } else if (jvr_method == "kumarathunge19"){
            tc_growth <- i_temp
            tc_home   <- i_thome
            tmp$jmax25[j] <- tmp$vcmax25[j] * ( 2.56 - 0.0375 * tc_home - 0.0202 * (tc_growth - tc_home))
          }
        }
        
        # Replacing stomatal sensitivity
        if (replace_xi) {
          
          xi_from_g1 <- i_g1 - ( i_rstar * (sqrt(i_vpd)  + i_g1) / i_ca)
          
          # cat("\n", tmp$sitename[j], " ", list_of_pfts$pft[i] ,": xi from", tmp$xi[j], " to ", xi_from_g1)
          tmp$xi[j]  <- xi_from_g1
        }
      }
    }
  }
  
  # Nest vars again
  tmp <- tmp %>% 
    nest(rpm_acc = all_of(rpm_acc_vars)) %>% 
    nest(forcing_growth = all_of(forcing_growth_vars)) %>% 
    nest(siteinfo = all_of(siteinfo_vars)) %>% 
    dplyr::select(-pft)
  
  df_out <- tmp
  
  # Check if in and out are the same
  
  if (all((names(df_out) %in% names(df_in)))) {
    return(df_out)
  } else {
    stop("Input and output dataframe are not the same but should be.")
  }
  
}
