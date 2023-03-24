# Replace acclimated for pft-specific Vcmax25 and Jmax25

# In: df_acc
# Out: df_acc
# Purpose: Replace vcmax25 and jmax25 from df_acc by with pft-specific values

replace_acclimated_with_pfts <- function(df_in,
                                         replace_pc = F,
                                         replace_xi = F) {

  # Check if input dataframe is of right shape
  if (!all(c("sitename", "date", "rpm_acc", "data_org", "fit_opt", "meas_cond",
             "siteinfo", "forcing_d", "forcing_growth") %in% names(df_in))) {
    stop("Input df has unintended variables.")
  }
  
  # Get list of vars to nest by 
  rpm_acc_vars <- names(df_in$rpm_acc[[1]])
  forcing_growth_vars <- names(df_in$forcing_growth[[1]])
  
  # Get list of pfts
  list_of_pfts <- get_list_of_fixed_pfts()
  
  tmp <- 
    df_in %>% 
    mutate(pft = purrr::map_chr(data_org, ~pull(., PFT) %>% unique())) %>% 
    unnest(rpm_acc) %>% 
    unnest(forcing_growth)
  
  # Replace values where needed
  for (i in 1:nrow(list_of_pfts)) {
    
    # To set Medlyn and Prentice models equal, they must have the same units
    i_g1 <- list_of_pfts$g1[i] * (1000)^0.5 # g1 is in units kPa ^ 0.5, we need Pa^0.5
    
    i_co2  <- tmp %>% dplyr::filter(pft == list_of_pfts$pft[i]) %>% pull(co2) %>% mean()  # Keep ppm, needed for Medlyn model
    i_patm <- tmp %>% dplyr::filter(pft == list_of_pfts$pft[i]) %>% pull(patm) %>% mean()
    i_vpd  <- tmp %>% dplyr::filter(pft == list_of_pfts$pft[i]) %>% pull(vpd) %>% mean()
    i_temp <- tmp %>% dplyr::filter(pft == list_of_pfts$pft[i]) %>% pull(temp) %>% mean()
    
    i_ca    <- co2_to_ca(i_co2, i_patm) # In Pa, needed for Prentice model
    i_rstar <- calc_gammastar(i_temp, i_patm) # 
    
    for (j in 1:nrow(tmp)) {
      
      if (tmp$pft[j] == list_of_pfts$pft[i]) {
        
        # Verbose Output
        # message("Sitename -- PFT -- old vcmax25 -- new vcmax25 \n ",
        #         tmp$sitename[j], " -- ",tmp$pft[j], " -- ", tmp$vcmax25[j], " -- ",list_of_pfts$vcmax25[i])
        
        # Replacing photosynthetic capacities
        if (replace_pc) {
          tmp$vcmax25[j] <- list_of_pfts$vcmax25[i] 
          tmp$rd25[j]    <- list_of_pfts$vcmax25[i] * 0.015
          tmp$jmax25[j]  <- list_of_pfts$jmax25[i] 
          tmp$kphio[j]   <- 0.04479
          
        # Replacing stomatal sensitivity
        }
        
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
    dplyr::select(-pft)
  
  df_out <- tmp
  
  # Check if in and out are the same
  
  if (all((names(df_out) %in% names(df_in)))) {
    return(df_out)
  } else {
    stop("Input and output dataframe are not the same but should be.")
  }
  
}
