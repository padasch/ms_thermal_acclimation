# Get forcing data ----

## Disclaimer ----
#' To run this file, {ingestr} must be installed. However, this packages is under
#' continual development without official releases. Therefore, the final output
#' of this file is provided in the manuscripts repository.

### Finalize forcing dataframe ----
df_forc_tmp <- readRDS("data/final/k19_df_forc_from cluster.rds")

## Attach CO2 on annual basis
df_forc_tmp2 <-
  df_forc_tmp %>% 
  unnest(forcing) %>% 
  left_join(ingestr::ingest(siteinfo = df_forc_tmp %>% 
                              unnest(siteinfo) %>% 
                              dplyr::select(-forcing),
                            source  = "co2_mlo",
                            verbose = FALSE,
                            timescale = "h") %>% 
              unnest(data) %>% 
              dplyr::select(-year)
  ) %>% 
  nest(forcing = !any_of(c("sitename", "siteinfo")))

## Calculate MAT and MAP per site and calculate leaf size
df_forc_tmp3 <-
  df_forc_tmp2 %>% 
  unnest(forcing) %>% 
  mutate(year = year(date),
         month = month(date))  %>% 
  nest(forcing_month = !any_of(c("sitename", "siteinfo", "year", "month"))) %>%
  mutate(tc_mean_month = purrr::map_dbl(forcing_month, ~pull(., temp) %>% mean(na.rm = T))) %>% 
  unnest(forcing_month) %>% 
  nest(forcing_year = !any_of(c("sitename", "siteinfo", "year"))) %>% 
  mutate(map   = purrr::map_dbl(forcing_year, ~pull(., prec) %>% mean(na.rm = T)),
         map   = map * 60*60*24,
         tc_wm = purrr::map_dbl(forcing_year, ~pull(., tc_mean_month) %>% max())) %>% 
  rowwise() %>% 
  mutate(leaf_width = calc_leaf_size(tc_wm, map),
         leaf_width = sqrt(leaf_width)/100) %>% 
  unnest(forcing_year) %>% 
  dplyr::select(-year, -month) %>% 
  nest(forcing_hh = !any_of(c("sitename", "siteinfo")))

## Save final forcing df for further analysis
saveRDS(df_forc_tmp3, here("data", "final", "k19_sitename_siteinfo_forcing_hh.rds"))
