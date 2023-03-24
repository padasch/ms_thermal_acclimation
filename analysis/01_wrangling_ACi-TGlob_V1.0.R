# Wrangling of ACi-TGlob V1.0 Dataset ----
# TODO: Add explanation what QC according to K19 means, maybe also where date was found when linking to K19


## 1. Source Packages ----
library(here)
source(here("R", "source.R"))

## 2. Download Data ----
# TODO: Add code to download dataset and add it to data/raw/...

## 3. Load Raw Data ----
df_raw <- read_csv(here("data", "raw", "ACi-TGlob", "ACi-TGlob_V1.0.csv"))

## 4. Checking validity of data ----
# TODO: Change lat-lon information for sites landing in water
df_raw <- 
    df_raw %>% 
    ## Correcting wronlgy spelled names
    mutate(Data_contributor = ifelse(str_detect(Data_contributor, "Cater"), "Kelsey Carter", Data_contributor))

## 5. Filter Raw Data ----
# TODO: Should I filter only for mature plants and hold-out saplings? 
#       Information is in docx in ACi-TGlob but not in the csv itself

df_tmp <- df_raw %>%
  dplyr::filter(
      ## Take only field measurements of native species
      Growth_condition %in% c("Field (NE)", "Field")) %>%
  
  ## Add unique sitename for identification using contributor name and sampling coordinates
  mutate(sitename = tolower(
    paste0(
      sub(".* ", "", Data_contributor), "_",
      abs(round(seed_source_latitude)),
      abs(round(seed_source_longitude))
    )
  )) %>%
    
  ## Nest by sitename for further wrangling
  nest(data = !any_of("sitename")) %>%
  arrange(sitename)

#_______________________________________________________________________________
## Unifying sampling dates ----
df_battaglia_43147 <-
  df_tmp[which(df_tmp$sitename == "battaglia_43147"), ]     %>%
  unnest(data) %>%
  mutate(sample_date = dmy(Date)) %>% 
  nest(data = !all_of(c("sitename")))  

df_carter_1866 <-
  df_tmp[which(df_tmp$sitename == "carter_1866"), ]     %>%
  unnest(data) %>%
  mutate(sample_date = dmy(Date)) %>% 
  nest(data = !all_of(c("sitename")))  

df_cavaleri_1866 <- 
  df_tmp[which(df_tmp$sitename == "cavaleri_1866"), ] %>%
  unnest(data) %>% 
  ## Add sampling dates from source data (see repository mentioned by K19)
  mutate(Date = ifelse(Curve_Id %in% c(1),      "09/08/2014", Date),
         Date = ifelse(Curve_Id %in% c(2:6),    "06/08/2014", Date),
         Date = ifelse(Curve_Id %in% c(7:11),   "27/08/2014", Date),
         Date = ifelse(Curve_Id %in% c(12),     "28/08/2014", Date),
         Date = ifelse(Curve_Id %in% c(13:14),  "11/08/2014", Date),
         Date = ifelse(Curve_Id %in% c(15:17),  "12/08/2015", Date),
         Date = ifelse(Curve_Id %in% c(18:25),  "10/03/2015", Date),
         Date = ifelse(Curve_Id %in% c(26:31),  "13/08/2014", Date),
         sample_date = dmy(Date)) %>% 
  nest(data = !all_of(c("sitename")))  

df_cernusak_14131 <-
  df_tmp[which(df_tmp$sitename == "cernusak_14131"), ] %>%
  unnest(data) %>%
  mutate(sample_date = dmy(Date)) %>% 
  nest(data = !all_of(c("sitename")))

# TODO: Dates missing and not findable in any doucment on photom repository or referenced publication
df_crous_34151 <- 
  df_tmp[which(df_tmp$sitename == "crous_34151"), ] %>%
    unnest(data) %>%
    mutate(sample_date = dmy(Date)) %>% 
  nest(data = !all_of(c("sitename")))

df_ellsworth_3679 <- 
  df_tmp[which(df_tmp$sitename == "ellsworth_3679"), ]  %>%
  unnest(data) %>% 
  mutate(dataset = "ellsworth",
         sample_date = dmy(paste("15-", Date))) %>% 
  nest(data = !all_of(c("sitename")))

df_han_35139 <- 
  df_tmp[which(df_tmp$sitename == "han_35139"), ] %>%
  unnest(data) %>%
  mutate(
    sample_date = Date,
    sample_date = ifelse(Curve_Id %in% 1:3, "2002-05-15", sample_date),
    sample_date = ifelse(Curve_Id %in% 4:6, "2001-07-15", sample_date),
    sample_date = ifelse(Curve_Id %in% 7:10, "2001-11-15", sample_date),
    sample_date = ymd(sample_date)) %>% 
  nest(data = !all_of(c("sitename")))

df_han_36140 <- 
  df_tmp[which(df_tmp$sitename == "han_36140"), ] %>% 
  unnest(data) %>% 
  # Assuming mid of month for reported months of sampling
  mutate(sample_date = dmy(paste0("15-", Date))) %>% 
  nest(data = !all_of(c("sitename")))

df_hikosaka_43142  <-
  df_tmp[which(df_tmp$sitename == "hikosaka_43142"), ] %>%
  unnest(data) %>% 
  # Dates reconstructed using "hikosaka_aci_data.csv" file in photom repository
  mutate(Date = ifelse(Date == "2001-06-01", "2001-06-11", Date),
         Date = ifelse(Date == "2001-08-01", "2001-08-10", Date),
         Date = ifelse(Date == "2001-09-01", "2001-09-20", Date),
         Date = ifelse(Date == "2002-06-01", "2002-06-10", Date),
         Date = ifelse(Date == "2002-07-01", "2002-07-29", Date),
         Date = ifelse(Date == "2002-09-01", "2002-09-28", Date),
         sample_date = dmy(Date)) %>% 
  nest(data = !all_of(c("sitename")))

df_jensen_4893 <-
  df_tmp[which(df_tmp$sitename == "jensen_4893"), ]     %>% 
  unnest(data) %>% 
  
  # Converting DOY and year information into day, -1 needed because as_date takes 0 as 1. January)
  left_join(read_csv(here("data", "raw", "files_from_photom_repo", "SPRUCE_3_cohort_ACi_data.csv")) %>% 
              dplyr::select(DOY, Year, Photo, Cond)) %>% 
  
  # To reconstruct date, file from photom repository is used "SPRUCE_3_cohort_ACi_data.csv"
  mutate(Date = as.character(as_date(DOY-1, origin = paste0(Year, "-01-01")))) %>% 
  dplyr::select(-DOY, -Year) %>% 
  mutate(sample_date = ymd(Date)) %>% 
  nest(data = !all_of(c("sitename")))

df_kelly_16145  <-
  df_tmp[which(df_tmp$sitename == "kelly_16145"), ]     %>%
  unnest(data) %>% 
  # Date reconstruction based on file "Daintree_ACidata_processed.csv" in photom repository
  mutate(season = ifelse(season == "summer", "2011-04-07", season),
         season = ifelse(season == "winter", "2010-07-22", season),
         sample_date = ymd(season)) %>% 
  nest(data = !all_of(c("sitename")))

df_medlyn_355 <-
  df_tmp[which(df_tmp$sitename == "medlyn_355"), ]      %>% 
  unnest(data) %>% 
  mutate(sample_date = dmy(Date)) %>% 
  nest(data = !all_of(c("sitename")))

df_medlyn_36148 <-
  df_tmp[which(df_tmp$sitename == "medlyn_36148"), ]      %>% 
  unnest(data) %>% 
  mutate(sample_date = dmy(Date)) %>% 
  nest(data = !all_of(c("sitename")))

df_medlyn_441 <-
  df_tmp[which(df_tmp$sitename == "medlyn_441"), ]      %>% 
  unnest(data) %>% 
  mutate(sample_date = dmy(Date)) %>% 
  nest(data = !all_of(c("sitename")))

df_onoda_41141     <-
  df_tmp[which(df_tmp$sitename == "onoda_41141"), ]     %>% 
  unnest(data) %>% 
  # Date reconstruction using reported sampling period in Onoda et al. (2005)
  mutate(Date = ifelse(str_detect(Date, "May"), paste0(mean.Date(c(dmy("08-05-2002"), dmy("12-05-2002")))),  Date),
         Date = ifelse(str_detect(Date, "Aug"), paste0(mean.Date(c(dmy("29-07-2002"), dmy("06-08-2002")))), Date),
         Date = ifelse(str_detect(Date, "Oct"), paste0(mean.Date(c(dmy("14-10-2002"), dmy("25-10-2002")))), Date),
         sample_date = ymd(Date)) %>% 
  nest(data = !all_of(c("sitename")))

df_rogers_71157     <-
  df_tmp[which(df_tmp$sitename == "rogers_71157"), ]     %>% 
  unnest(data) %>% 
  
  # To reconstruct date, file from photom repository is used "Arctic_A-Ci_curves_2012-2015_V2.csv"
  # Cannot bind df's by "Photo" or "Cond" due to differences in decimals, thus using row id arranged by "Photo"
  arrange(Photo) %>% 
  mutate(row_id = row_number()) %>% 
  dplyr::select(-Date) %>% 
  left_join(read_csv(here("data", "raw", "files_from_photom_repo", "Arctic_A-Ci_curves_2012-2015_V2.csv")) %>% 
              arrange(Photo) %>% 
              mutate(row_id = row_number(),
                     Date = as.character(Sample_Date)) %>% 
              dplyr::select(Date, row_id)) %>% 
  dplyr::select(-row_id) %>% 
  mutate(sample_date = ymd(Date)) %>% 
  nest(data = !all_of(c("sitename")))

df_slot_980 <-
  df_tmp[which(df_tmp$sitename == "slot_980"), ] %>%
  unnest(data) %>%
  mutate(sample_date = dmy(Date)) %>% 
  nest(data = !all_of(c("sitename")))

df_tarvainen_5812  <- 
  df_tmp[which(df_tmp$sitename == "tarvainen_5812"), ]  %>% 
  unnest(data) %>% 
  # TODO: Date reconstruction very vague, only reported as period in K19 data
  mutate(sample_date = Date,
         sample_date = ifelse(str_detect(Date, "Jun-10"),       paste0("2010-06-15"), sample_date),
         sample_date = ifelse(str_detect(Date, "Jun-sep-2009"), paste0("2009-08-01"), sample_date),
         sample_date = ymd(sample_date)) %>% 
  nest(data = !all_of(c("sitename")))

df_tarvainen_6420  <- 
  df_tmp[which(df_tmp$sitename == "tarvainen_6420"), ]  %>% 
  unnest(data) %>% 
  mutate(sample_date = dmy(Date)) %>% 
  nest(data = !all_of(c("sitename")))

df_togashi_30121  <- 
  df_tmp[which(df_tmp$sitename == "togashi_30121"), ]  %>% 
  unnest(data) %>% 
  # TODO: Could not identify any clear date information (not in original data, not in publication)
  mutate(sample_date = dmy(Date)) %>% 
  nest(data = !all_of(c("sitename")))

df_tribuzy_360  <- 
  df_tmp[which(df_tmp$sitename == "tribuzy_360"), ]  %>% 
  unnest(data) %>% 
  mutate(sample_date = dmy(Date)) %>% 
  nest(data = !all_of(c("sitename")))

df_wang_6331  <- 
  df_tmp[which(df_tmp$sitename == "wang_6331"), ]  %>% 
  unnest(data) %>% 
  mutate(sample_date = dmy(ifelse(Date == "Jul-98", "15-07-1994", Date))) %>% 
  nest(data = !all_of(c("sitename")))

# Bind all dataframes together
df_tmp_date <-
  rbind(
    df_battaglia_43147, 
    df_cavaleri_1866,
    df_carter_1866,
    df_cernusak_14131,  
    df_crous_34151,     
    df_ellsworth_3679,  
    df_han_35139,       
    df_han_36140,       
    df_hikosaka_43142,  
    df_jensen_4893,     
    df_kelly_16145,     
    df_medlyn_355,      
    df_medlyn_36148,    
    df_medlyn_441,      
    df_onoda_41141,     
    df_rogers_71157,    
    df_slot_980,        
    df_tarvainen_5812,  
    df_tarvainen_6420,  
    df_togashi_30121,   
    df_tribuzy_360,     
    df_wang_6331)

# Summary
per_dates_pre <-
  round((df_tmp %>%
           unnest(data) %>% 
           dplyr::filter(
             is.na(Date)) %>% nrow()) / (nrow(df_tmp %>% unnest(data))) * 100, 0)

per_dates_post <-
  round((df_tmp_date %>%
           unnest(data) %>% 
           dplyr::filter(
             is.na(Date)) %>% nrow()) / (nrow(df_tmp %>% unnest(data))) * 100, 0)

cat("\n Data points without given date BEFORE unifying dates: ", per_dates_pre, " %",
    "\n Data points without given date AFTER  unifying dates: ", per_dates_post, " %")

#__________________________________________________________________________________________________#
## 6. T_opt data ----

### General Filtering ----
# Remove data outside reasonable CO2 range (250-450ppm)

# Plots for inspection
# df_tmp_date %>%
#   dplyr::filter(!is.na(CO2S)) %>% 
#   ggplot() +
#   geom_histogram(aes(x=CO2S)) +
#   facet_wrap(~sitename, scales = "free")
# 
# df_tmp_date %>%
#   dplyr::filter(!is.na(CO2R)) %>% 
#   ggplot() +
#   geom_histogram(aes(x=CO2R)) +
#   facet_wrap(~sitename, scales = "free")
# 
# plot(df_tmp_date$CO2R, df_tmp_date$CO2S)
# abline(0, 1)

# Filter
df_tmp_date_flt <-
  df_tmp_date %>% 
  unnest(data) %>% 
  dplyr::filter(
    between(Ci, 150, 450) &
      (between(CO2S, 250, 450) | between(CO2R, 250, 450)))
  
# Plots for inspection
# visdat::vis_miss(df_tmp_date_flt %>% 
#                    dplyr::select(Ci, CO2R, CO2S, VpdL, Tleaf, Photo), 
#                  show_perc_col = T,
#                  sort_miss = T)

### Extraction ----
#### battaglia_43147 ----
df_battaglia_43147     <-
    df_tmp_date_flt[which(df_tmp$sitename == "battaglia_43147"), ]     %>%
    mutate(agg_date = sample_date) %>%
    nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
    mutate(
      fit_opt = purrr::map(data_org,  ~fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
      fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "daily"))
      )

## Check output
# df_battaglia_43147 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_battaglia_43147)

#### carter_1866 ----
df_carter_1866 <-
  df_tmp_date_flt[which(df_tmp_date$sitename == "carter_1866"), ] %>%
  mutate(agg_date = sample_date) %>%
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org,~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "daily"))
    )

# df_carter_1866 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_carter_1866)

#### cavaleri_1866 ----
df_cavaleri_1866 <-
  df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "cavaleri_1866")) %>%
  # TODO: Removed data from 13/08/2014 because it is has a flat rate of GPP over all temperatures (could be done in the end?)
  # filter(!str_detect(date, "13/08/2014")) %>%
  mutate( # agg_date = sample_date) %>% # Take daily values
    agg_date = floor_date(sample_date, unit = "week") # Flooring to week
  ) %>% 
  # agg_date = floor_date(sample_date, unit = "month")) %>% # Flooring to month
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "weekly"))
    )

## Check output
# df_cavaleri_1866 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_cavaleri_1866)

#### cernusak_14131 ----
df_cernusak_14131 <-
  df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "cernusak_14131")) %>%
  dplyr::filter(
    between(CO2R, 380, 440), # Only use data taken at ambient CO2
    PARi > 1800
  ) %>% # Only use data at light-saturation
  mutate(agg_date = sample_date) %>%
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "daily"))
    )

## Check output
# df_cernusak_14131 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_cernusak_14131)

#### crous_34151 ----
df_crous_34151 <-
  df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "crous_34151")) %>%
  dplyr::filter(between(CO2R, 380, 440)) %>% # Only use data taken at ambient CO2
  mutate(agg_date = sample_date) %>%
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "daily"))
    )

## Check output
# df_crous_34151 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_crous_34151)


#### ellsworth_3679 ----
df_ellsworth_3679 <-
  df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "ellsworth_3679")) %>%
  dplyr::filter(between(CO2S, 340, 370)) %>% # QC according to code in repository by K19
  mutate(agg_date = sample_date) %>%
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "daily"))
    )

## Check output
# df_ellsworth_3679 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_ellsworth_3679)

#### han_35139 ----
df_han_35139 <- df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "han_35139")) %>%
  # QC according to code by K19 in photom repository
  dplyr::filter(between(CO2S, 330, 350)) %>%
  mutate(agg_date = sample_date) %>%
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "daily"))
    )

## Check output
# df_han_35139 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_han_35139)

#### han_36140 ----
df_han_36140 <-
  df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "han_36140")) %>%
  # QC according to code by K19 in photom repository
  dplyr::filter(between(Ci, 175, 230)) %>%
  mutate(agg_date = sample_date) %>%
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "daily"))
    )

## Check output
# df_han_36140 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_han_36140)

#### hikosaka_43142 ----
df_hikosaka_43142 <-
  df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "hikosaka_43142")) %>%
  # QC according to code by K19 in photom repository
  dplyr::filter(
    between(CO2S, 360, 375),
    !(Curve_Id %in% c(108:110, 115:117))
  ) %>%
  mutate(agg_date = sample_date) %>%
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "daily"))
    )

## Check output
# df_hikosaka_43142 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_hikosaka_43142)

#### jensen_4893 ----
df_jensen_4893 <-
  df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "jensen_4893")) %>%
  # QC according to code by K19 in photom repository
  dplyr::filter(between(CO2R, 395, 410)) %>%
  mutate(agg_date = sample_date) %>% # agg_date = floor_date(sample_date, unit = "month")) %>% # Flooring to week
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "daily"))
    )

## Check output
# df_jensen_4893 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_jensen_4893)

#### kelly_16145 ----
df_kelly_16145 <-
  df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "kelly_16145")) %>%
  # QC according to code by K19 in photom repository
  dplyr::filter(between(CO2R, 375, 425), Ci > 175) %>%
  mutate(agg_date = sample_date) %>%
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "daily"))
    )

## Check output
# df_kelly_16145 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_kelly_16145)

#### medlyn_355 ----
df_medlyn_355 <-
  df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "medlyn_355")) %>%
  # QC according to code by K19 in photom repository
  dplyr::filter(between(Ci, 175, 350), Curve_Id != 62) %>%
  mutate(agg_date = sample_date) %>%
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "daily"))
    )

## Check output
# df_medlyn_355 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_medlyn_355)

#### medlyn_36148 ----
df_medlyn_36148 <-
  df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "medlyn_36148")) %>%
  # QC according to code by K19 in photom repository
  dplyr::filter(PARi > 1200, Cond > 0.09) %>%
  mutate(agg_date = sample_date) %>%
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "daily"))
    )

## Check output
# df_medlyn_36148 %>% dplyr::select(-data_org) %>% unnest(fit_opt) %>% drop_na(aopt)
# plot_topt_extraction(df_medlyn_36148)

#### medlyn_441 ----
df_medlyn_441 <-
  df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "medlyn_441")) %>%
  # QC according to code by K19 in photom repository
  dplyr::filter(between(Ci, 175, 350), Curve_Id != 62) %>%
  mutate(agg_date = sample_date) %>%
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "daily"))
    )

## Check output
# df_medlyn_441 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_medlyn_441)

#### onoda_41141 ----
# Comment: Sampling data shows no clear pattern of A_net - T_leaf
df_onoda_41141 <-
  df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "onoda_41141")) %>%
  # QC according to code by K19 in photom repository
  dplyr::filter(Photo < 8.2, between(Ci, 250, 350)) %>%
  mutate(agg_date = sample_date) %>%
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "daily"))
    )

## Check output
# df_onoda_41141 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_onoda_41141)

#### rogers_71157 ----
df_rogers_71157 <-
  df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "rogers_71157")) %>%
  # QC according to code by K19 in photom repository
  dplyr::filter(between(CO2S, 370, 410)) %>%
  mutate(agg_date = round_date(sample_date, unit = "week")) %>%
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "weekly"))
    )

## Check output
# df_rogers_71157 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_rogers_71157)

#### slot_980 ----
df_slot_980 <-
  df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "slot_980")) %>%
  mutate(agg_date = sample_date) %>%
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "daily"))
    )

## Check output
# df_slot_980 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_slot_980)

#### tarvainen_5812 ----
df_tarvainen_5812 <-
  df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "tarvainen_5812")) %>%
  # QC according to code by K19 in photom repository
  dplyr::filter(between(CO2R, 395, 405)) %>%
  mutate(agg_date = sample_date) %>%
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "daily"))
    )

## Check output
# df_tarvainen_5812 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_tarvainen_5812)

#### tarvainen_6420 ----
df_tarvainen_6420 <-
  df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "tarvainen_6420")) %>%
  # QC according to code by K19 in photom repository
  dplyr::filter(between(CO2R, 395, 405)) %>%
  mutate(agg_date = floor_date(sample_date, unit = "week")) %>% # Flooring to week
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "weekly"))
    )

## Check output
# df_tarvainen_6420 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_tarvainen_6420)

#### togashi_30121 ----
df_togashi_30121 <-
  df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "togashi_30121")) %>%
  # QC according to code by K19 in photom repository
  dplyr::filter(between(CO2R, 365, 435), !(Curve_Id %in% c(19, 23, 27, 44))) %>%
  mutate(agg_date = sample_date) %>%
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "daily"))
    )

## Check output
# df_togashi_30121 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_togashi_30121)

#### tribuzy_360 ----
df_tribuzy_360 <-
  df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "tribuzy_360")) %>%
  # QC according to code by K19 in photom repository
  dplyr::filter(between(CO2R, 365, 435), !(Curve_Id %in% c(19, 23, 27, 44))) %>%
  mutate(agg_date = sample_date) %>%
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "daily"))
    )

## Check output
# df_tribuzy_360 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_tribuzy_360)

#### wang_6331 ----
df_wang_6331 <-
  df_tmp_date_flt %>% 
  dplyr::filter(str_detect(sitename, "wang_6331")) %>%
  # QC according to code by K19 in photom repository
  dplyr::filter(between(Ci, 240, 400)) %>%
  mutate(agg_date = sample_date) %>%
  nest(data_org = !any_of(c("sitename", "agg_date"))) %>%
  mutate(
    fit_opt = purrr::map(data_org, ~ fit_nonlinear_topt(dat = .x, x = "Tleaf", y = "Photo", random = "Species")),
    fit_opt = purrr::map(fit_opt, ~mutate(., agg_scale = "daily"))
    )

## Check output
# df_wang_6331 %>% dplyr::select(-data_org) %>% unnest(fit_opt)
# plot_topt_extraction(df_wang_6331)

### Bind all datasets together ----
df_pre_qc <-
  rbind(
    df_battaglia_43147,
    df_cavaleri_1866,
    df_carter_1866,
    df_cernusak_14131,
    df_crous_34151,
    df_ellsworth_3679,
    df_han_35139,
    df_han_36140,
    df_hikosaka_43142,
    df_jensen_4893,
    df_kelly_16145,
    df_medlyn_355,
    df_medlyn_36148,
    df_medlyn_441,
    df_onoda_41141,
    df_rogers_71157,
    df_slot_980,
    df_tarvainen_5812,
    df_tarvainen_6420,
    df_togashi_30121,
    df_tribuzy_360,
    df_wang_6331
  ) %>%
  mutate(id = paste0(sitename, "_", agg_date)) 
  

#__________________________________________________________________________________________________#
## 7. Quality Control ----
## Filter datasets based on goodness-of-fit
df_post_qc <- 
    df_pre_qc %>% 
    mutate(n_qc = purrr::map_dbl(data_org, ~nrow(.))) %>% 
    unnest(fit_opt) %>% 
    drop_na(topt, agg_date) %>%      # Remove missing fits or missing agg_dates
    dplyr::filter(n_qc > 5,          # Take out all fits with less than 5 data points
                  aopt.se < 5,       # Take out all fits where standard error of A_opt is bigger than 5 µmol/m2/s
                  b > 0,             # Take out all fits with inverted parabola
                  topt.se < 5) %>%   # Take out all fits where standard error of T_opt is bigger than 5 degC
    dplyr::select(-n_qc) %>% 
    nest(fit_opt = !all_of(c("sitename", "agg_date", "data_org", "id")))

## Get sites that were dropped in automatic routine
ids_of_kept_sites <- df_post_qc$id
df_dropped <- df_pre_qc %>% dplyr::filter(!(id %in% ids_of_kept_sites)) 

## Adding/Removing based on visual inspection of plots
add_to_final <-
  c("medlyn_36148_2002-05-10" # Reasonable fit but with SE = 5.94 was included nonetheless
    )

remove_from_final <-
  c(
    # Data from two leaf temperatures that are spread far away which gives
    # an unreasonably good parabola fit
    "ellsworth_3679_1998-12-15",
    "ellsworth_3679_1999-08-15",
    
    # Data from 35 different curves that have visibly different T_opt but
    # no additional information given to extract proper fit
    "jensen_4893_2013-08-02",
    
    ## Unreasonably low assimilation rates, potentially a unit error in original data
    "jensen_4893_2013-04-23",
    "jensen_4893_2013-04-24",
    "jensen_4893_2013-04-25"
    )

df_post_qc <-
  df_dropped %>% 
  filter(id %in% add_to_final) %>%
  bind_rows(df_post_qc) %>% 
  filter(!(id %in% remove_from_final))

ids_of_kept_sites <- df_post_qc$id
df_dropped <- df_pre_qc %>% dplyr::filter(!(id %in% ids_of_kept_sites)) 

# Check if no sites is in both dfs:
if (T %in% c(df_post_qc$id %in% df_dropped$id, df_dropped$id %in% df_post_qc$id) |
    nrow(df_dropped) + nrow(df_post_qc) != nrow(df_pre_qc)) {
  stop("Careful, exchanging sites between final and dropped df did not work properly!")
}

cat("Number of points removed through QC: ", nrow(df_dropped), "out of total: ", nrow(df_pre_qc))

## Make plots of final thermal response curves
dir_1 <- here("output", "figures", "tc_opt_extraction", "final_sites")
dir_2 <- here("output", "figures", "tc_opt_extraction", "dropped_sites")

big_plot_qc      <- plot_topt_extraction(df_post_qc, make_one_plot = T, dir = dir_1)
big_plot_dropped <- plot_topt_extraction(df_dropped, make_one_plot = T, dir = dir_2)

#_______________________________________________________________________________
## Adding Siteinfo
df_final <-
    df_post_qc %>% 
    mutate(
       
         ## Data from nested data
        lat = purrr::map_dbl(data_org, ~pull(., seed_source_latitude)  %>% unique()),
        lon = purrr::map_dbl(data_org, ~pull(., seed_source_longitude) %>% unique()),
        pft = purrr::map_chr(data_org, ~pull(., PFT) %>% unique())) %>% 
        
    ## Extract start and end of agg_date to get sampling data
    nest(nest_tmp = !any_of("sitename")) %>% 
    mutate(year_start = purrr::map_dbl(nest_tmp,  ~pull(., agg_date) %>% as.character() %>% min() %>% year()),
           year_end   = purrr::map_dbl(nest_tmp,  ~pull(., agg_date) %>% as.character() %>% max() %>% year())) %>% 
    unnest(nest_tmp) %>%
    mutate(
        
        ## Data from explanatory file "data_descriptor_V1.0.docx" of ACi-T Database
        tc_growth_air_k19 = NA,
        tc_growth_air_k19 = ifelse(str_detect(sitename, "battaglia_43147"), 18.4, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "carter_1866"), 24.0, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "cavaleri_1866"), 24.0, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "cernusak_14131"), 37.8, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "crous_34151"), 17.1, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "ellsworth_3679"), 31.4, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "han_35139"), 26.0, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "han_36140"), 30.1, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "hikosaka_43142"), 23.8, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "jensen_4893"), 25.8, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "kelly_16145"), 31.3, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "medlyn_355"), 26.7, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "medlyn_36148"), 22.9, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "medlyn_441"), 26.7, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "onoda_41141"), mean(15.9, 25.3), tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "rogers_71157"), 7.7, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "slot_980"), 26.9, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "tarvainen_5812"), 20.1, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "tarvainen_6420"), 20.1, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "togashi_30121"), 34.5, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "tribuzy_360"), 32.7, tc_growth_air_k19),
        tc_growth_air_k19 = ifelse(str_detect(sitename, "wang_6331"), 20.4, tc_growth_air_k19),
    
        ## Data from explanatory file "data_descriptor_V1.0.docx" of ACi-T Database
        tc_home_k19 = NA,
        tc_home_k19 = ifelse(str_detect(sitename, "battaglia_43147"), 8.5, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "carter_1866"), 29.8, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "cavaleri_1866"), 29.8, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "cernusak_14131"), 27.4, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "crous_34151"), 28.9, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "ellsworth_3679"), 14.6, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "han_35139"), 11.8, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "han_36140"), 13.7, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "hikosaka_43142"), 12.3, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "jensen_4893"), 12.2, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "kelly_16145"), 24.4, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "medlyn_355"), 12.6, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "medlyn_36148"), 8.5, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "medlyn_441"), 12.6, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "onoda_41141"), NA, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "rogers_71157"), 3.1, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "slot_980"), 32.6, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "tarvainen_5812"), 9.4, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "tarvainen_6420"), 8.5, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "togashi_30121"), 18.6, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "tribuzy_360"), 27.2, tc_home_k19),
        tc_home_k19 = ifelse(str_detect(sitename, "wang_6331"), 8.4, tc_home_k19),
        
        ## Experimental ppfd as reported in SI of Kumarathunge et al. (2019)
        ppfd_measurement = NA,
        ppfd_measurement = ifelse(str_detect(sitename, "battaglia_43147"), 1500, ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "carter_1866"), 800, ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "cavaleri_1866"), mean(600, 800), ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "cernusak_14131"), 2000, ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "crous_34151"), 1800, ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "ellsworth_3679"), 1800, ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "han_35139"), 1100, ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "han_36140"), 1100, ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "hikosaka_43142"), 1000, ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "jensen_4893"), 1700, ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "kelly_16145"), 1000, ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "medlyn_355"), 1400, ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "medlyn_36148"), 1500, ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "medlyn_441"), 1400, ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "onoda_41141"), mean(1000, 2000), ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "rogers_71157"), 2000, ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "slot_980"), 1500, ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "tarvainen_5812"), 1047, ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "tarvainen_6420"), 1500, ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "togashi_30121"), 1800, ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "tribuzy_360"), NA, ppfd_measurement),
        ppfd_measurement = ifelse(str_detect(sitename, "wang_6331"), 1200, ppfd_measurement),
        
        ## Experimental conditions from original data
        meas_cond   = purrr::map(data_org, ~dplyr::select(., CO2S, CO2R, VpdL, Ci, PARi, Tleaf)),) %>% 
    
    nest(sitedata = all_of(c("data_org", "fit_opt", "meas_cond"))) %>% 
    nest(siteinfo = !any_of(c("sitename", "sitedata", "agg_date")))

## Add information on climate zone to siteinfo
df_final <-  add_climatezone_to_siteinfo(df_final)

### Save final T_opt df  ----
filn <- here("data", "final", "k19_sitename_siteinfo_sitedata.rds")
saveRDS(df_final, filn)


## 8. Analysis of cleaned data ----
### Load df with processed data ----
df_final <- readRDS(here("data", "final", "k19_sitename_siteinfo_sitedata.rds"))

### Global map ----
p_map <- plot_global_samples(df_final)
# p_map + dark_theme_classic() + theme(legend.position = "bottom") + ggtitle(NULL)
dir_tmp <- here("output/global_maps")
if (!dir.exists(dir_tmp)) dir.create(dir_tmp, recursive = T, showWarnings = F)
ggsave(paste0(dir_tmp, "/global_map.pdf"), p_map, width = 9, height = 7)

### Check data availability and distributions ----
## Ridges to see distributions of data
# df_final %>%
#   unnest(sitedata) %>% 
#   unnest(meas_cond) %>% 
#   pivot_longer(cols = c("CO2S", "CO2R", "Ci", "VpdL", "Tleaf", "PARi")) %>% 
#   mutate(name = as.factor(name)) %>% 
#   ggplot() +
#   ggridges::geom_density_ridges(aes(y = sitename, x = value),
#                                 alpha = 0.5,
#                                 jittered_points = TRUE,
#                                 point_alpha=0.25,
#                                 # point_shape=1,
#                                 scale = 0.9) +
#   facet_wrap(~name, scales = "free_x", nrow = 1) +
#   theme_classic()
# 
# visdat::vis_miss(
#   df_final %>%
#   unnest(sitedata) %>% 
#   unnest(data_org) %>% 
#   dplyr::select(Ci, CO2R, CO2S, VpdL, Tleaf, Photo), 
#   show_perc_col = T,
#   sort_miss = T)

### Meta Info Table ----
tab_meta <- 
  df_final %>% 
  unnest(c(sitedata, siteinfo)) %>% 
  unnest(c(fit_opt, data_org)) %>% 
  group_by(sitename) %>% 
  nest() %>% 
  mutate(n = purrr::map_chr(data, ~pull(., agg_date) %>% unique() %>% length()),
         agg_scale = purrr::map_chr(data, ~pull(., agg_scale) %>% unique() %>% str_flatten(collapse = ", ")),
         agg_date = purrr::map_chr(data, ~pull(., agg_date) %>% unique() %>% str_flatten(collapse = ", ")),
         species = purrr::map_chr(data, ~pull(., Species) %>% unique() %>% str_flatten(collapse = ", ")),
         species = str_replace(species, "�", replacement = " "),
         # pft = purrr::map_chr(data, ~pull(., PFT) %>% unique() %>% str_flatten(collapse = ", ")),
         org_ref = purrr::map_chr(data, ~pull(., Reference) %>% unique() %>% str_flatten(collapse = ", ")),
         org_ref = str_remove(org_ref, "\xa0"),
         lat = purrr::map_chr(data, ~pull(., lat) %>% unique()),
         lon = purrr::map_chr(data, ~pull(., lon) %>% unique()),
         cl = purrr::map_chr(data, ~pull(., cl) %>% unique()),
         pft = purrr::map_chr(data, ~pull(., PFT) %>% unique()),
         agg_scale = purrr::map_chr(data, ~pull(., agg_scale) %>% unique())) %>%
  dplyr::select(-data)

dir_tmp <- here("output/metadata")
if (!dir.exists(dir_tmp)) dir.create(dir_tmp, recursive = T, showWarnings = F)

write.csv(tab_meta, paste0(dir_tmp, "/si_site-metadata.csv"), row.names = FALSE)

## 9. Get forcing data ----
#__________________________________________________________________________________________________#
### Euler Code to Extract Forcing ----
# 2022-07-01: This code in in file `./R/rscript_ingest_siteinfo `
# Save final df_forc to data/final!
#__________________________________________________________________________________________________#

### Finalize forcing dataframe ----
df_forc_tmp <- readRDS("data/final/k19_df_forc_from cluster.rds")

## Attach CO2 on annual basis
df_forc_tmp2 <-
    df_forc_tmp %>% 
    unnest(forcing) %>% 
    left_join(ingest(siteinfo = df_forc_tmp %>% 
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

#_______________________________________________________________________________
# Extraction of Vcmax ----
#' I have to extract all the values
#' I have to take an average based on sampling dates
#' Then I can predict these vcmax values based on sampling dates and compare them
#' For this I need:
#'   - Vcmax instant from data
#'   - Tleaf instant from data
#'   - Vcmax acc from rpmodel
#'   - Tgrowth (and Thome) from climate data
#'   - I do not even have to look at the sampling-date average because we
#'     are predicting instant vcmax depending on the leaf temperature 
#'     they were measured on! I only need the sampling date for calculating the
#'     acclimated Vcmax values. But from these values I can then predict inst
#'     vcmax using the parametrization by Kumarathunge et al. directly from 
#'     T_growth to instant T_leaf.
#'     OR I could downscale either to vcmax25 using the parametrization BUT
#'     I should not use the plantecophys::fitaci() option that returns Vcmax25,
#'     I should rather use Vcmax at the given leaf temperature directly.
#_______________________________________________________________________________

## Nest dataframe by curves
df_curves <-
  df_tmp_date %>%
  unnest(data) %>% 
  ## Nest original data by A-Ci curve
  nest(curve_data = !any_of(c("sitename", "sample_date", "Curve_number"))) %>% 
  ## Drop sites without sampling dates
  drop_na(sample_date) %>% 
  arrange(Curve_number)

## Loop through curves and extract aci fits
df_acifits <- tibble()

for (i in 1:nrow(df_curves)) {
  
  df_i <- df_curves %>% slice(i)
  message("> i = ", i, "/", nrow(df_curves), " | ", df_i$sitename, " ", df_i$sample_date, " ", df_i$Curve_number)

  try_out <- try(
    df_i <-
      df_i %>%
      mutate(fitaci = purrr::map(curve_data, ~ fitaci(.,
        fitmethod = "bilinear",
        Tcorrect = F
      ))),
    silent = T
  )

  if (inherits(try_out, "try-error")) {
    df_i$fitaci_qc <- F
    df_i$fitaci <- NA
  } else {
    df_i$fitaci_qc <- T
  }
  
  df_acifits <- rbind(df_acifits, df_i)
}

cat("\n Percentage of failed acifits: ",
    round((df_acifits %>% dplyr::filter(fitaci_qc == F) %>% nrow())/
      (df_acifits %>% nrow()), 2)*100, " %")

# Extract and summarize traits to site and sample_date level
df_traits <-
  df_acifits %>% 
  ## Remove curves where fitaci-function failed:
  dplyr::filter(fitaci_qc == T) %>% 
  ## Extract values for instant vcmax, jmax and mean leaf temperature
  mutate(
    fit_tleaf = purrr::map_dbl(fitaci, ~ mean(.$df$Tleaf, na.rm = T)),
    fit_vcmax = purrr::map_dbl(fitaci, ~ .$pars[[1]] * 1e-6),
    fit_jmax  = purrr::map_dbl(fitaci, ~ .$pars[[2]] * 1e-6),
    fit_rmse  = purrr::map_dbl(fitaci, ~ .$RMSE * 1e-6)
    # fit_vcmax_se    = purrr::map_dbl(fitaci, ~.$pars[[4]] * 1e-6),
    # fit_jmax_se     = purrr::map_dbl(fitaci, ~.$pars[[5]] * 1e-6),
  ) %>%
  ## Aggregate data to date of aggregation and attach mean values for
  ## photosynthetic traits to fit_opt
  unnest(curve_data) %>% 
  nest(site_date_data = !all_of(c("sitename", "sample_date"))) %>% 
  mutate(tleaf = purrr::map_dbl(site_date_data, ~mean(.$fit_tleaf, na.rm = T)),
         tleaf_se = purrr::map_dbl(site_date_data, ~standard_error(.$fit_tleaf)),
         vcmax = purrr::map_dbl(site_date_data, ~mean(.$fit_vcmax, na.rm = T)),
         vcmax_se = purrr::map_dbl(site_date_data, ~standard_error(.$fit_vcmax)),
         jmax = purrr::map_dbl(site_date_data, ~mean(.$fit_jmax, na.rm = T)),
         jmax_se = purrr::map_dbl(site_date_data, ~standard_error(.$fit_jmax)),
         rmse = purrr::map_dbl(site_date_data, ~mean(.$fit_rmse, na.rm = T))) %>% 
  nest(traits_data = !all_of(c("sitename", "sample_date", "site_date_data")))

# Attach siteinfo from above
df_traits <-
  df_traits %>% 
  inner_join(df_final %>% dplyr::select(sitename, siteinfo) %>% distinct())


# Vcmax Analyis -----------------------------------------------------------
df_traits %>% 
  unnest(traits_data) %>% 
  ggplot() +
  geom_boxplot(aes(x = vcmax, y = sitename)) +
  geom_point(aes(x = vcmax, y = sitename, fill = rmse), shape = 21)


# Save rd files
saveRDS(df_traits %>% dplyr::select(-site_date_data), here("data", "final", "traits_without_raw_data.rds"))
# saveRDS(df_traits, here("data", "final", "traits_with_raw_data.rds")) # Outcommented because takes too much time

## End of Script ___________________________________________________________________________________
