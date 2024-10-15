# Sensitivity analysis ----
## FUNCTIONS ----
sensana_inst_to_plot <- function(df_inst){
  do_not_nest <- names(df_inst)
  
  df_tmp <-
    df_inst %>% 
    ## Remove rpm_inst where agross values were equal (e.g. 0 because no light)
    mutate( inst_qc   = purrr::map_dbl(rpm_inst, ~ pull(., agross) %>% unique() %>% length())) %>% 
    mutate( tspan_qc  = purrr::map_lgl(rpm_sim, ~ pull(., tspan) %>% is.na())) %>% 
    dplyr::filter(inst_qc > 1, tspan_qc != TRUE) %>% 
    dplyr::select(-inst_qc, -tspan_qc) %>% 
    
    ## Extract relevant variables
    mutate(
      tc_opt        = purrr::map_dbl(rpm_inst, ~ slice_max(., anet)   %>% pull(tc_leaf)),
      tc_opt_agross = purrr::map_dbl(rpm_inst, ~ slice_max(., agross) %>% pull(tc_leaf)),
      tc_opt_ac     = purrr::map_dbl(rpm_inst, ~ slice_max(., ac_net) %>% pull(tc_leaf)),
      tc_opt_aj     = purrr::map_dbl(rpm_inst, ~ slice_max(., aj_net) %>% pull(tc_leaf)),
      
      lim_anet      = purrr::map_chr(rpm_inst, ~ slice_max(., anet)   %>% pull(min_a)),
      lim_anet      = as.factor(lim_anet),
      
      agross_opt    = purrr::map_dbl(rpm_inst, ~ slice_max(., agross) %>% pull(agross) * 1e6),
      ac_net_opt    = purrr::map_dbl(rpm_inst, ~ slice_max(., ac_net) %>% pull(ac_net) * 1e6),
      aj_net_opt    = purrr::map_dbl(rpm_inst, ~ slice_max(., aj_net) %>% pull(aj_net) * 1e6),
      anet_opt      = purrr::map_dbl(rpm_inst, ~ slice_max(., anet)   %>% pull(anet)   * 1e6),
      
      # anet_growth   = purrr::map_dbl(rpm_inst, ~ dplyr::filter(., tc_growth_air == tc_leaf) %>% pull(anet) * 1e6),
      # min_agrowth   = purrr::map_chr(rpm_inst, ~ dplyr::filter(., tc_growth_air == tc_leaf) %>% pull(min_a)),
      # min_agrowth   = as.factor(min_agrowth),
      
      tspan_l       = purrr::map_dbl(rpm_inst, 
                                     ~ dplyr::filter(., anet > (0.9 * max(anet))) %>% 
                                       slice_min(tc_leaf) %>% 
                                       pull(tc_leaf)),
      tspan_h       = purrr::map_dbl(rpm_inst, 
                                     ~ dplyr::filter(., anet > (0.9 * max(anet))) %>% 
                                       slice_max(tc_leaf) %>% 
                                       pull(tc_leaf)),
      
      tspan_anet    = tspan_h - tspan_l,
      
      tspan_l       = purrr::map_dbl(rpm_inst, 
                                     ~ dplyr::filter(., ac_net > (0.9 * max(ac_net))) %>% 
                                       slice_min(tc_leaf) %>% 
                                       pull(tc_leaf)),
      tspan_h       = purrr::map_dbl(rpm_inst, 
                                     ~ dplyr::filter(., ac_net > (0.9 * max(ac_net))) %>% 
                                       slice_max(tc_leaf) %>% 
                                       pull(tc_leaf)),
      
      tspan_ac    = tspan_h - tspan_l,
      
      tspan_l       = purrr::map_dbl(rpm_inst, 
                                     ~ dplyr::filter(., aj_net > (0.9 * max(aj_net))) %>% 
                                       slice_min(tc_leaf) %>% 
                                       pull(tc_leaf)),
      tspan_h       = purrr::map_dbl(rpm_inst, 
                                     ~ dplyr::filter(., aj_net > (0.9 * max(aj_net))) %>% 
                                       slice_max(tc_leaf) %>% 
                                       pull(tc_leaf)),
      
      tspan_aj    = tspan_h - tspan_l
    ) %>%
    nest(traits_inst  = !all_of(do_not_nest))
  
  return(df_tmp)
}

# REFERENCE DF ----
# Start with working dataframe
df_forcing <- 
  k19_create_input_for_acc(settings = get_settings()) %>% 
  dplyr::select(-data_org, 
                -fit_opt,
                -meas_cond) %>%  
  slice(1) %>% 
  mutate(sitename = "dummy-site",
         date     = "dummy-date")

# Overwrite current forcings with standard values
df_std <- get_df_ref()

df_forcing$forcing_growth[[1]]$co2     <- df_std$co2
df_forcing$forcing_growth[[1]]$patm    <- df_std$patm
df_forcing$forcing_growth[[1]]$ppfd    <- df_std$ppfd
df_forcing$forcing_growth[[1]]$temp    <- df_std$temp
df_forcing$forcing_growth[[1]]$vpd     <- df_std$vpd

df_forcing$forcing_d[[1]]$co2     <- df_std$co2
df_forcing$forcing_d[[1]]$patm    <- df_std$patm
df_forcing$forcing_d[[1]]$ppfd    <- df_std$ppfd
df_forcing$forcing_d[[1]]$temp    <- df_std$temp
df_forcing$forcing_d[[1]]$vpd     <- df_std$vpd

df_forcing$siteinfo[[1]]$tc_home  <- df_std$temp
df_forcing$siteinfo[[1]]$ppfd_measurement  <- df_std$ppfd

# Attaching measuring conditions
steps_inst <- 0.5

df_forcing <-
  left_join(df_forcing,
            tibble(tc_leaf = seq(0, 40, steps_inst)) %>%
              mutate(sitename = "dummy-site",
                     date     = "dummy-date",
                     tc_air   = 25,
                     patm     = 101325,
                     vpd      = 1000,
                     co2      = 400,
                     ppfd     = 1500e-6) %>%
              nest(forcing_inst = !c(sitename, date))
            )


# GROWTH FORCING ----
## Replicate dataframes for sensitivity analysis ----
# Get sensitivity df
n_steps <- 20
df_sens0 <- tibble()

for (i in 1:n_steps) {
  df_sens0 <- rbind(df_sens0, df_forcing)
}

df_for_acc <- tibble()
vars <- c("co2", "patm", "ppfd", "temp_k19", "temp_l02", "vpd", "tc_home")

for (v in vars) {
  
  # Get fixed range of current variable
  if (str_detect(v, "temp") | 
      v == "tc_home")     r <- seq(0,40, length.out = n_steps)
  if (v == "ppfd")        r <- seq(0, 2000e-6, length.out = n_steps)
  if (v == "vpd")         r <- seq(10, 6000, length.out = n_steps)
  if (v == "patm")        r <- seq(60000, 110000, length.out = n_steps)
  if (v == "co2")         r <- seq(280, 560, length.out = n_steps)
  
  
  # Loop over each row in df_sens1
  df_sens1 <- df_sens0 # Get fresh df for every variable
  
  for (i in 1:nrow(df_sens1)) {
    
    df_sens1$sitename[[i]] <- paste0(v)
    df_sens1$date[[i]]     <- paste0(round(r[i], 6))
    
    if (v == "tc_home") {
      df_sens1$siteinfo[[i]]$tc_home <- r[i]
      
    } else if (str_detect(v, "temp")) {
      df_sens1$forcing_growth[[i]]$temp <- r[i]
      
    } else {
      df_sens1$forcing_growth[[i]][[v]] <- r[i]
    }
  }
  
  df_for_acc <- rbind(df_for_acc, df_sens1)
}

## Run acc to inst ----
## K19 part
settings <- get_settings()
df_acc_k19   <- k19_run_acc(df_for_acc %>% dplyr::filter(!str_detect(sitename, "l02")), 
                            settings)
df_inst_k19  <- k19_run_inst(df_acc_k19, settings, sensana = TRUE)

## L02 part
settings <- get_settings()
settings$method_ftemp <- "kumarathunge19_fixed"
df_acc_l02   <- k19_run_acc(df_for_acc %>% dplyr::filter(str_detect(sitename, "l02")), 
                            settings)
df_inst_l02  <- k19_run_inst(df_acc_l02, settings, sensana = TRUE)

# Concatenate
df_sensana_growth <- sensana_inst_to_plot(rbind(df_inst_l02, df_inst_k19))

# Remove ER sensitivity analysis because Leuning2002 not used anymore
df_sensana_growth <- df_sensana_growth |> filter(sitename != "temp_l02")

## Plot it ----
labels <- 
  c(
    'tc_opt' = "T_opt [°C]",
    't_span' = "T_span [°C]",
    'a_opt' = "A_opt [mu mol/m2/s]",
    'temp_l02' = "Tgrowth fixed [°C]",
    'temp_k19' = "Tgrowth [°C]",
    'tc_home' = "T_home [°C]",
    'co2' = "CO2 [ppm]",
    'patm' = "P_atm [kPa]",
    'vpd' = "VPD_air [kPa]",
    'ppfd' = "I_abs [mu mol/m2/s]")

p_growth <- 
  df_sensana_growth %>% 
  dplyr::select(sitename, date, traits_inst) %>% 
  unnest(traits_inst) %>% 
  pivot_longer(cols = c("tc_opt_ac", "tc_opt_aj",
                        "ac_net_opt", "aj_net_opt",
                        "tspan_ac", "tspan_aj"),
               values_to = "trait_value",
               names_to = "trait") %>% 
  mutate(rate = ifelse(str_detect(trait, "ac"), "Ac", "Aj"),
         rate = as.factor(rate),
         trait = ifelse(str_detect(trait, "tc"), "tc_opt", trait),
         trait = ifelse(str_detect(trait, "net"), "a_opt", trait),
         trait = ifelse(str_detect(trait, "tspan"), "t_span", trait),
         var = as.factor(sitename),
         val = as.double(date),
         val = ifelse(var == "patm", val/1000, val),
         val = ifelse(var == "vpd", val/1000, val),
         val = ifelse(var == "ppfd", val*1000, val)) %>% 
  ggplot() +
  aes(val, 
                trait_value, 
                color = rate, 
                group = rate) +
  # geom_point() +
  geom_smooth(linewidth = 1.5, se = F) +
  facet_grid(vars(trait), 
             vars(var),
             labeller = as_labeller(labels),
             switch = "both",
             scales = "free") +
  ylim(0, 35)  +
  labs(x = "Trait Value",
       y = "Parameter Value",
       color = "Assimil.\nRate",
       title = "Sensitivity of traits against growth forcing") +
  theme(plot.title = element_text(face = "bold"))


# INSTANT FORCING ----
## Replicate dataframes for sensitivity analysis ----
# Get sensitivity df
n_steps <- 20
df_sens0 <- tibble()

for (i in 1:n_steps) {
  df_sens0 <- rbind(df_sens0, df_forcing)
}

df_for_acc <- tibble()
vars <- c("co2", "patm", "ppfd", "temp", "vpd")

for (v in vars) {
  
  # Get fixed range of current variable
  if (v == "temp")     r <- seq(0,40, length.out = n_steps)
  if (v == "ppfd")     r <- seq(0, 2000e-6, length.out = n_steps)
  if (v == "vpd")      r <- seq(10, 6000, length.out = n_steps)
  if (v == "patm")     r <- seq(60000, 110000, length.out = n_steps)
  if (v == "co2")      r <- seq(280, 560, length.out = n_steps)
  
  # Loop over each row in df_sens1
  df_sens1 <- df_sens0 # Get fresh df for every variable
  
  for (i in 1:nrow(df_sens1)) {
    
    df_sens1$sitename[[i]] <- paste0(v)
    df_sens1$date[[i]]     <- paste0(round(r[i], 6))

    df_sens1$forcing_inst[[i]][[v]] <- rep(r[i], nrow(df_sens1$forcing_inst[[i]]))

  }
  
  df_for_acc <- rbind(df_for_acc, df_sens1)
}

## Run acc to inst ----
settings <- get_settings()
df_acc   <- k19_run_acc(df_for_acc, settings)
df_inst  <- k19_run_inst(df_acc, settings, sensana = TRUE)

df_sensana_inst <- sensana_inst_to_plot(df_inst) 

## Plot it ----
labels <- 
  c(
    'tc_opt' = "T_opt [°C]",
    't_span' = "T_span [°C]",
    'a_opt' = "A_opt [mu mol/m2/s]",
    'temp' = "T_air [°C]",
    'co2' = "CO2 [ppm]",
    'patm' = "P_atm [kPa]",
    'vpd' = "VPD_air [kPa]",
    'ppfd' = "I_abs [mu mol/m2/s]")

p_inst <- 
  df_sensana_inst %>% 
  dplyr::select(sitename, date, traits_inst) %>% 
  unnest(traits_inst) %>% 
  pivot_longer(cols = c("tc_opt_ac", "tc_opt_aj",
                        "ac_net_opt", "aj_net_opt",
                        "tspan_ac", "tspan_aj"),
               values_to = "trait_value",
               names_to = "trait") %>% 
  mutate(rate = ifelse(str_detect(trait, "ac"), "Ac", "Aj"),
         rate = as.factor(rate),
         trait = ifelse(str_detect(trait, "tc"), "tc_opt", trait),
         trait = ifelse(str_detect(trait, "net"), "a_opt", trait),
         trait = ifelse(str_detect(trait, "tspan"), "t_span", trait),
         var = as.factor(sitename),
         val = as.double(date),
         val = ifelse(var == "patm", val/1000, val),
         val = ifelse(var == "vpd", val/1000, val),
         val = ifelse(var == "ppfd", val*1000, val)) %>% 
  ggplot() +
  aes(val, 
                trait_value, 
                color = rate, 
                group = rate) +
  # geom_point() +
  geom_smooth(linewidth = 1.5, se = F) +
  facet_grid(vars(trait), 
             vars(var),
             labeller = as_labeller(labels),
             switch = "both",
             scales = "free") +
  ylim(0, 35)  +
  labs(x = "Trait Value",
       y = "Parameter Value",
       color = "Assimil.\nRate",
       title = "Sensitivity of traits against instantaneous forcing") +
  theme(plot.title = element_text(face = "bold"))

# ACCLIMATION ----
## Get reference acclimated dataframe ----
df_ref_acc <-
 left_join(df_forcing,
           tibble(
             sitename = "dummy-site",
             date     = "dummy-date",
             vcmax25 = 50e-6,
             jmax25  = 100e-6,
             rd25    = 50e-6 * 1.5/100,
             xi      = 60,
             kphio   = 0.07,
             tc_leaf = NA
           ) %>%
             nest(rpm_acc = !c(sitename, date))
 )

## Leuning02 ----
### Get sensitivity range
n_steps <- 20
df_sens0 <- tibble()

for (i in 1:n_steps) {
  df_sens0 <- rbind(df_sens0, df_ref_acc)
}

df_for_acc <- tibble()
vars <- c("vcmax25", "jmax25", "rd25", "xi", "kphio", "temp_leuning")

for (v in vars) {
  
  # Get fixed range of current variable
  if (v == "vcmax25")    r <- seq(5,100, length.out = n_steps) / 1e6 
  if (v == "jmax25")     r <- seq(5,100, length.out = n_steps) / 1e6 *2
  if (v == "rd25")       r <- seq(5,100, length.out = n_steps) / 1e6 * 1.5 / 100
  if (v == "xi")         r <- seq(10,140, length.out = n_steps)
  if (v == "kphio")      r <- seq(0.01, 0.125, length.out = n_steps)
  if (v == "temp_leuning") r <- seq(0,40, length.out = n_steps)

  
  # Loop over each row in df_sens1
  df_sens1 <- df_sens0 # Get fresh df for every variable
  
  for (i in 1:nrow(df_sens1)) {
    
    df_sens1$sitename[[i]] <- paste0(v)
    df_sens1$date[[i]]     <- paste0(round(r[i], 6))
    
    if (str_detect(v, "25")) df_sens1$date[[i]] <- paste0(round(r[i]*1e6, 2))

    
    
    if (v == "temp_leuning") {
      df_sens1$forcing_growth[[i]]$temp <- r[i]
    } else {
      df_sens1$rpm_acc[[i]][[v]] <- r[i]
    }
  }
  
  df_for_acc <- rbind(df_for_acc, df_sens1)
}

### Run model
settings <- get_settings()
settings$method_ftemp <- "kumarathunge19_fixed"
df_inst  <- k19_run_inst(df_for_acc, settings, sensana = TRUE)
df_sensana_with_leuning <- sensana_inst_to_plot(df_inst) 

## Kumarathunge2019
### Get 
n_steps <- 20
df_sens0 <- tibble()

for (i in 1:n_steps) {
  df_sens0 <- rbind(df_sens0, df_ref_acc)
}

df_for_acc <- tibble()
vars <- c("temp_kumarathunge")

for (v in vars) {
  
  # Get fixed range of current variable
  if (v == "temp_kumarathunge") r <- seq(0,40, length.out = n_steps)
  
  # Loop over each row in df_sens1
  df_sens1 <- df_sens0 # Get fresh df for every variable
  
  for (i in 1:nrow(df_sens1)) {
    
    df_sens1$sitename[[i]] <- paste0(v)
    df_sens1$date[[i]]     <- paste0(round(r[i], 6))
    df_sens1$forcing_growth[[i]]$temp <- r[i]
  }
  
  df_for_acc <- rbind(df_for_acc, df_sens1)
}

### Run model
settings <- get_settings()
df_inst  <- k19_run_inst(df_for_acc, settings, sensana = TRUE)
df_sensana_with_kumara <- sensana_inst_to_plot(df_inst) 

df_sensana_acc <- 
  rbind(df_sensana_with_leuning,
        df_sensana_with_kumara)

# Remove ER sensitivity, not needed anymore after removing Leuning2002
df_sensana_acc <- df_sensana_acc |> filter(
  sitename != "temp_leuning",
  sitename != "temp_kumarathunge"
)

## Plot it ----
### Labeller
labels <- 
  c(
    'tc_opt' = "T_opt [°C]",
    't_span' = "T_span [°C]",
    'a_opt' = "A_opt [mu mol/m2/s]",
    'vcmax25' = "vcmax25 [mu mol/m2/s]",
    'jmax25' = "jmax25 [mu mol/m2/s]",
    'rd25' = "rd25 [mu mol/m2/s]",
    'xi' = "xi [Pa^0.5]",
    'kphio' = "kphio [-]",
    'temp_leuning' = "Tgrowth fixed [°C]",
    'temp_kumarathunge' = "Tgrowth flexible [°C]")

p_acc <- 
  df_sensana_acc %>% 
  dplyr::select(sitename, date, traits_inst) %>% 
  unnest(traits_inst) %>% 
  pivot_longer(cols = c("tc_opt_ac", "tc_opt_aj",
                        "ac_net_opt", "aj_net_opt",
                        "tspan_ac", "tspan_aj"),
               values_to = "trait_value",
               names_to = "trait") %>% 
  mutate(rate = ifelse(str_detect(trait, "ac"), "Ac", "Aj"),
         rate = as.factor(rate),
         trait = ifelse(str_detect(trait, "tc"), "tc_opt", trait),
         trait = ifelse(str_detect(trait, "net"), "a_opt", trait),
         trait = ifelse(str_detect(trait, "tspan"), "t_span", trait),
         var = as.factor(sitename),
         val = as.double(date),
         val = ifelse(var == "patm", val/1000, val),
         val = ifelse(var == "vpd", val/1000, val),
         val = ifelse(var == "ppfd", val*1000, val)) %>% 
  ggplot() +
  aes(val, 
      trait_value, 
      color = rate, 
      group = rate) +
  # geom_point() +
  geom_smooth(linewidth = 1.5, se = F) +
  facet_grid(vars(trait), 
             vars(var),
             labeller = as_labeller(labels),
             switch = "both",
             scales = "free") +
  ylim(0, 35)  +
  labs(x = "Trait Value",
       y = "Parameter Value",
       color = "Assimil.\nRate",
       title = "Sensitivity of traits against acclimation-related parameters") +
  theme(plot.title = element_text(face = "bold"))


# SAVE PLOTS ----
dir.create(here("output/sensitivity_analysis/", today()), recursive = T, showWarnings = F)

ggsave(paste0(here("output/sensitivity_analysis/", today()), "/sensana_growth_forcing.pdf"),
       p_growth,
       height = 6,
       width  = 15)

ggsave(paste0(here("output/sensitivity_analysis/", today()), "/sensana_inst_forcing.pdf"),
       p_inst,
       height = 6,
       width  = 15)

ggsave(paste0(here("output/sensitivity_analysis/", today()), "/sensana_acc.pdf"),
       p_acc,
       height = 6,
       width  = 15)
