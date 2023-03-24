# Function to return list of fixed PFTs

# Source of vcmax25: Rogers et al. (2014), Table 1, Averages of clearly defined PFTs
# Source of ratio jmax/vcmax25: Kattge & Knorr (2007)
# Source of g1: CLM5 Tech Note

get_list_of_fixed_pfts <- function() {
  
  pfts <- 
    data.frame(
    pft     = c("ARCTIC", "NET_TE", "BET_TE", "NET_B", "BDT_TE", "BET_Tr"),
    vcmax25 = c(      36,       49,       63,     48,        60,       68), # umol/m2/s
    g1      = c(    2.22,     2.35,     4.12,   2.35,      4.45,     4.12)  # kPa
    )
  
  # Turn umol into mol
  pfts$vcmax25 <- pfts$vcmax25 * 1e-6
  
  # Calcualte non-acclimated jmax25
  pfts$jmax25 <- pfts$vcmax25 * 1.88
  
  return(pfts)
}


