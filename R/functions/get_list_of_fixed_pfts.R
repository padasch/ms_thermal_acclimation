# Function to return list of fixed PFTs

# Source of vcmax25: Rogers et al. (2014), Table 1, Averages of clearly defined PFTs
# Source of ratio jmax/vcmax25: Kattge & Knorr (2007)
# Source of g1: CLM5 Tech Note

get_list_of_fixed_pfts <- function(vcmax25_source = "none", return_only_g1 = FALSE) {
  
  
  if (return_only_g1){
    # Creat pfts dataframe  
    # Attention: This is a lazy bugfix! When changing g1 parameters, make sure to change them below too.
    pfts <- 
      data.frame(
        pft = c("Arctic", "BDT-Te", "BET-Te", "BET-Tr", "NET-Bo", "NET-Te"),
        g1  = c(    2.22,     4.45,     4.12,   4.12,      2.35,     2.35) # in kPa
      )
    
    return(pfts)
  }
  
  # Get model-specific pfts
  model_vcmax25 <-  tibble(
    pfts       = c("Arctic", "BDT-Te", "BET-Te", "BET-Tr", "NET-Bo", "NET-Te"),
    avim       = c(NA,       60,       68,       64,       58,       60),
    bethy      = c(20,       54,       58,       NA,       58,       NA),
    clm        = c(52,       52,       72,       72,       54,       61),
    # hybrid     = c(37,       NA,       NA,       NA,       NA,       NA),
    ibis       = c(NA,       75,       100,      163,      63,       75),
    kumar19    = c(78.3,     39,       82.9,     39.4,     80.4,     42.8),
    ocn        = c(NA,       NA,       59,       24,       17,       17),
    orchidee   = c(NA,       NA,       21,       18,       38,       32),
    mean       = c(46.8,     56.0,     65.8,     63.4,     52.6,     48.0),
    mean_nok19 = c(36.3,     60.3,     63.0,     68.2,     48.0,	   49.0)
  )
  
  # Check if input is valid
  avl_models = names(model_vcmax25 |> select(-pfts))
  if (!(vcmax25_source %in% avl_models)) {
    stop("Invalid vcmax25_source input. Pick one of:\n\t ", paste0(avl_models, sep = " | "))
  }
  
  # Replace entries that are NA with the overall mean value
  for (irow in 1:nrow(model_vcmax25)) {
    model_mean = model_vcmax25$mean[irow]
    for (icol in 1:ncol(model_vcmax25)) {
      if (is.na(model_vcmax25[irow, icol])) {
        model_vcmax25[irow, icol] <- model_mean
      }
    }
  }
  
  # Creat pfts dataframe  
  pfts <- 
    data.frame(
      pft = c("Arctic", "BDT-Te", "BET-Te", "BET-Tr", "NET-Bo", "NET-Te"),
      g1  = c(    2.22,     4.45,     4.12,   4.12,      2.35,     2.35), # in kPa
      vcmax25 = model_vcmax25[[vcmax25_source]] * 1e-6
    )
  
  return(pfts)
}

