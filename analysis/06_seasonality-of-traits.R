# SEASONALITY OF TRAITS ----
source("R/source.R")

# ├ Meta toggles ----
debug_slice <- F
step_size_tcair   <- 1
set_add <- list()
set_add$min_topt_per_year <- 3
set_add$without_noise <- T
set_add$high_light    <- T
set_add$vpd_scaled    <- F
vcmax25_source <- "kumar19"

dir_analysis <- "seasonality"
dir_tmp      <- "final_stepsize-0.5"


## Define directories ----
# First check how many runs have been made to add counter
dir_mods <- make_new_dynamic_dir(dir_analysis, "model_rds", dir_tmp)
dir_figs <- make_new_dynamic_dir(dir_analysis, "figures", dir_tmp)
dir_tabs <- make_new_dynamic_dir(dir_analysis, "tables", dir_tmp)

# ├ FORCING ----------------------------------
# ├— df_ref - Get relevant sites (minimum of N samplings per year) -----------------------------------------------------
df_ref_org <- k19_create_input_for_acc(get_settings()) # Settings should be irrelevant here

df_ref_org$site_year <- NA
for (i in 1:nrow(df_ref_org)) {
  
  # cat("\n", i, "/", nrow(df_ref_org))
  
  df_ref_org$site_year[[i]] =  paste0(
    str_extract_all(df_ref_org$sitename[[i]], "[a-z]+"),
    year(df_ref_org$date[[i]]))
}

df_ref <- 
  df_ref_org %>% 
  nest(data = !any_of(c("site_year"))) %>%
  mutate(n = purrr::map_dbl(data, ~nrow(.))) %>% 
  dplyr::filter(n >= set_add$min_topt_per_year) %>% 
  dplyr::select(-n) %>% 
  unnest(data)

df_ref$tmp_id <- NA
for (i in 1:nrow(df_ref)) {
  df_ref$tmp_id[[i]] <- paste0(df_ref$sitename[[i]], df_ref$date[[i]])
}

# Update df_ref to hold fitted response curves
df_tmp <- df_ref %>% mutate(fit_thermal_response = list(tibble()))

for (i in 1:nrow(df_tmp)) {
  df_tmp$fit_thermal_response[[i]] <- 
    tibble(
      tleaf = seq(0, 50, 0.1),
      anet = df_tmp$fit_opt[[i]]$aopt - ( df_tmp$fit_opt[[i]]$b * ( tleaf - df_tmp$fit_opt[[i]]$aopt) ^ 2),
      tc_growth_air = round(df_tmp$forcing_growth[[i]]$temp/0.1)*0.1)
}

df_ref <-
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
  dplyr::select(-fit_opt) %>% 
  nest(fit_opt = any_of(c("tc_opt", "tc_opt_se", "anet_opt", "anet_opt_se", "agrowth", "tspan_l", "tspan_h", "tspan"))) %>%
  dplyr::select(sitename, date, fit_opt) %>%
  right_join(df_ref %>% dplyr::select(-fit_opt)) %>% 
  dplyr::select(-c(forcing_d, forcing_growth))

# ├— B15 ----
settings <- get_settings()
settings$save_plots <- F
settings$daily_conditions <- "B"
settings$tau <- 15

# Get forcing
df_forc_org <- k19_create_df_forcing(settings)
nest_vars <- c("date", names(df_forc_org)[1:4])
df_forc_b15 <-
  df_forc_org %>% 
  dplyr::select(!c(forcing_hh, forcing_d)) %>% 
  unnest(forcing_growth) %>% 
  nest(forcing_growth = !any_of(nest_vars)) %>% 
  left_join(df_forc_org %>% 
              dplyr::select(sitename, forcing_d) %>% 
              unnest(forcing_d)) %>% 
  nest(forcing_d = !any_of(c("sitename", "date", "siteinfo", "forcing_growth"))) 

# Expand to per-day-level
df_forc_b15$site_year <- NA
for (i in 1:nrow(df_forc_b15)) {
  
  # cat("\n", i, "/", nrow(df_forc_b15))
  
  df_forc_b15$site_year[[i]] =  paste0(
    str_extract_all(df_forc_b15$sitename[[i]], "[a-z]+"),
    year(df_forc_b15$date[[i]]))
}

# Attach tmp_id
df_forc_b15$tmp_id <- NA
for (i in 1:nrow(df_forc_b15)) {
  df_forc_b15$tmp_id[[i]] <- paste0(df_forc_b15$sitename[[i]], df_forc_b15$date[[i]])
}

# Reduce to available site_year
df_forc_b15 <-
  df_forc_b15 %>% 
  dplyr::filter(site_year %in% df_ref$site_year)

# Attach measurement data and overwrite original data to get PFTs when replacing traits 
df_forc_b15$meas_cond <- list(NULL)
df_forc_b15$data_org  <- list(NULL)
df_forc_b15$fit_opt   <- list(NULL)

for (i in 1:nrow(df_ref)) {
  
  # Get reference data
  i_sy <- df_ref %>% slice(i) %>% pull(site_year)
  id   <- df_ref %>% slice(i) %>% pull(tmp_id)
  df_ref_i <- df_ref[which(df_ref$tmp_id == id), ]
  
  # Overwrite in df_forc_b15
  df_forc_b15[which(df_forc_b15$tmp_id == id), ]$meas_cond <- df_ref_i$meas_cond
  df_forc_b15[which(df_forc_b15$tmp_id == id), ]$fit_opt   <- df_ref_i$fit_opt  
  df_forc_b15[which(df_forc_b15$site_year == i_sy), ]$data_org  <- df_ref_i$data_org 
}

# ├ Replace daily forcing with growth forcing to remove day-to-day noise ----
if (set_add$without_noise) {
  for (i in 1:nrow(df_forc_b15)) {
    
      # Replace daily by growth conditions
      df_forc_b15$forcing_d[[i]]    <- df_forc_b15$forcing_growth[[i]]
      # df_forc_noacc$forcing_d[[i]]  <- df_forc_b15$forcing_growth[[i]]
  }
}

# ├ Run Setups ----
# ├— Full Acc. ----------------------------------------------------------------------
settings <- get_settings()
settings$save_plots <- F
df_acc_fullacc   <- k19_run_acc(df_forc_b15, settings)
df_inst_fullacc  <- run_inst_for_365_days(df_acc_fullacc, 
                                          settings, 
                                          setup = "fullacc", 
                                          dir_mods,
                                          set_add)

# ├— No Acc. ----------------------------------------------------------------------
settings <- get_settings()
settings$save_plots <- F
settings$method_vcmax25 = "kumar19"
settings$method_ftemp <- "kumarathunge19_fixed"

df_acc_noacc   <- k19_run_acc(df_forc_b15, settings)
df_acc_noacc   <- replace_acclimated_with_pfts(df_acc_noacc, TRUE, TRUE, jvr_method = settings$method_ftemp, vcmax25_source)
df_inst_noacc  <- run_inst_for_365_days(df_acc_noacc, 
                                        settings, 
                                        setup = "noacc", 
                                        dir_mods,
                                        set_add)

saveRDS(df_inst_noacc, here(dir_mods, "/noacc/df_inst.rds"))


# # ├— ER - Enzymatic Reactions (ftemp) ------------------------------------------
# 
# ## Update Settings
# settings <- get_settings()
# settings$save_plots <- F
# settings$method_vcmax25 = "kumar19"
# settings$method_ftemp <- "kumarathunge19"
# 
# df_acc_er   <- k19_run_acc(df_forc_b15, settings)
# df_acc_er   <- replace_acclimated_with_pfts(df_acc_er, TRUE, TRUE, jvr_method = settings$method_ftemp, vcmax25_source)
# df_inst_er  <- run_inst_for_365_days(df_acc_er, 
#                                         settings, 
#                                         setup = "er", 
#                                         dir_mods,
#                                         set_add)
# 
# saveRDS(df_inst_er, here(dir_mods, "/er/df_inst.rds"))
# 
# # ├— SB - Stomatal Behavior (Xi) -----------------------------------------------
# settings <- get_settings()
# settings$save_plots <- F
# settings$method_vcmax25 = "kumar19"
# settings$method_ftemp <- "kumarathunge19_fixed"
# 
# df_acc_sb   <- k19_run_acc(df_forc_b15, settings)
# df_acc_sb   <- replace_acclimated_with_pfts(df_acc_sb, replace_pc = TRUE, replace_xi = FALSE, jvr_method = settings$method_ftemp, vcmax25_source)
# df_inst_sb  <- run_inst_for_365_days(df_acc_sb, 
#                                      settings, 
#                                      setup = "sb", 
#                                      dir_mods,
#                                      set_add)
# 
# saveRDS(df_inst_sb, here(dir_mods, "/sb/df_inst.rds"))
# 
# 
# # ├— PC - Photosynthetic Capacity (Vcmax_acc, Jmax_acc) -------------------------
# settings <- get_settings()
# settings$save_plots <- F
# settings$method_ftemp <- "kumarathunge19_fixed"
# 
# df_acc_pc   <- k19_run_acc(df_forc_b15, settings)
# df_acc_pc   <- replace_acclimated_with_pfts(df_acc_pc, replace_pc = FALSE, replace_xi = TRUE, jvr_method = settings$method_ftemp, vcmax25_source)
# df_inst_pc  <- run_inst_for_365_days(df_acc_pc, 
#                                      settings, 
#                                      setup = "pc", 
#                                      dir_mods,
#                                      set_add)
# 
# saveRDS(df_inst_pc, here(dir_mods, "/pc/df_inst.rds"))

# ├ Analysis ------------------------------------------
# ├— Load Data ------------------------------------------

df_inst_fullacc <- readRDS(here(dir_mods, "/fullacc/df_inst.rds"))
df_inst_noacc <- readRDS(here(dir_mods, "/noacc/df_inst.rds"))
# df_inst_er <- readRDS(here(dir_mods, "/er/df_inst.rds"))
# df_inst_sb <- readRDS(here(dir_mods, "/sb/df_inst.rds"))
# df_inst_pc <- readRDS(here(dir_mods, "/pc/df_inst.rds"))

# ├— Plots (individual) ------------------------------------------
p_fullacc <- temporal_individual_plots(df_inst_fullacc, "fullacc", dir_figs = dir_figs)
p_noacc <- temporal_individual_plots(df_inst_noacc, "noacc", dir_figs = dir_figs)
# p_er <- temporal_individual_plots(df_inst_er, "er", dir_figs = dir_figs)
# p_sb <- temporal_individual_plots(df_inst_sb, "sb", dir_figs = dir_figs)
# p_pc <- temporal_individual_plots(df_inst_pc, "pc", dir_figs = dir_figs)

# Reduce dfs and combine them per setup
df_fullacc <- 
  df_inst_fullacc %>% 
  dplyr::select(sitename, date, fit_opt, rpm_sim, site_year) %>% 
  rename(mod = rpm_sim, obs = fit_opt) %>% 
  unnest(c(mod, obs), names_sep = "_") %>% 
  mutate(setup = "fullacc")

df_noacc   <- 
  df_inst_noacc %>% 
  dplyr::select(sitename, date, fit_opt, rpm_sim, site_year) %>% 
  rename(mod = rpm_sim, obs = fit_opt) %>% 
  unnest(c(mod, obs), names_sep = "_") %>% 
  mutate(setup = "noacc")

# df_er      <- 
#   df_inst_er %>% 
#   dplyr::select(sitename, date, fit_opt, rpm_sim, site_year) %>% 
#   rename(mod = rpm_sim, obs = fit_opt) %>% 
#   unnest(c(mod, obs), names_sep = "_") %>% 
#   mutate(setup = "er")
# 
# df_sb      <- 
#   df_inst_sb %>% 
#   dplyr::select(sitename, date, fit_opt, rpm_sim, site_year) %>% 
#   rename(mod = rpm_sim, obs = fit_opt) %>% 
#   unnest(c(mod, obs), names_sep = "_") %>% 
#   mutate(setup = "sb")
# 
# df_pc      <- 
#   df_inst_pc %>% 
#   dplyr::select(sitename, date, fit_opt, rpm_sim, site_year) %>% 
#   rename(mod = rpm_sim, obs = fit_opt) %>% 
#   unnest(c(mod, obs), names_sep = "_") %>% 
#   mutate(setup = "pc")

# ├— Plots (multiple) ------------------------------------------
# ├—— noacc - fullacc - observations ----

# Aesthetics
vec_cols <- RColorBrewer::brewer.pal(6, "Dark2")

sh <- c(rep(NA, 2), 21)
sz <- c(rep(1, 2), 2)
lt <- c(rep(1, 2), 0)
fl <- c(rep(NA, 2), "black")
cl <- c(vec_cols[1:2], "black")
br <- c("noacc", "fullacc", "obs")
lb <- c("No Acc.", "Full Acc.", "Observ.")


# df_all <- 
#   rbind(df_fullacc, 
#         df_noacc) %>%
#   left_join(df_inst_fullacc %>% 
#               dplyr::select(sitename, date, forcing_growth, forcing_d, site_year))
# 
# # T_opt
# p_topt_season <- 
#   df_all %>% 
#   ggplot() +
#   aes(x = date) +
#   geom_line(aes(y = mod_tc_opt, group = setup, color = setup)) +
#   geom_errorbar(aes(y = obs_tc_opt, 
#                       ymin = obs_tc_opt - obs_tc_opt_se, 
#                       ymax = obs_tc_opt + obs_tc_opt_se,
#                       color = "obs")) +
#   geom_point(aes(y = obs_tc_opt, 
#                  color = "obs",
#                  ),
#              fill = "black",
#              size = 2,
#              shape = 21) +
#   scale_color_manual(name = "",
#                      breaks = br,
#                      labels = lb,
#                      values = cl,
#                      guide = guide_legend(override.aes = list(shape  = sh,
#                                                               size   = sz,
#                                                               linetype = lt,
#                                                               fill  = fl,
#                                                               color = cl))) +
#   guides(fill = "none") +
#   theme_linedraw() +
#   facet_wrap(~site_year, scales = "free_x", ncol = 1) +
#   ylim(0, 40) +
#   labs(y = bquote("Temperature [°C]"),
#        x = "Date") 
# 
# ggsave(paste0(dir_figs, "/noacc-fullacc_topt.png"), p_topt_season, width = 12, height = 12)
# 
# ## A_opt
# p_aopt_season <- 
#   df_all %>% 
#   ggplot() +
#   aes(x = date) +
#   geom_line(aes(y = mod_anet_opt, group = setup, color = setup)) +
#   geom_errorbar(aes(y = obs_tc_opt, 
#                     ymin = obs_anet_opt - obs_anet_opt_se, 
#                     ymax = obs_anet_opt + obs_anet_opt_se,
#                     color = "obs")) +
#   geom_point(aes(y = obs_anet_opt, 
#                  color = "obs",
#   ),
#   fill = "black",
#   size = 2,
#   shape = 21) +
#   scale_color_manual(name = "",
#                      breaks = br,
#                      labels = lb,
#                      values = cl,
#                      guide = guide_legend(override.aes = list(shape  = sh,
#                                                               size   = sz,
#                                                               linetype = lt,
#                                                               fill  = fl,
#                                                               color = cl))) +
#   guides(fill = "none") +
#   theme_linedraw() +
#   facet_wrap(~site_year, scales = "free_x", ncol = 1) +
#   ylim(0, 40) +
#   labs(y = bquote(A[opt] ~ "[µmol m¯² s¯¹]"),
#        x = "Date") 
# 
# ggsave(paste0(dir_figs, "/noacc-fullacc_aopt.png"), p_aopt_season, width = 12, height = 12)
# 
# ## T_span
# p_tspan_season <- 
#   df_all %>% 
#   ggplot() +
#   aes(x = date) +
#   geom_line(aes(y = mod_tspan, group = setup, color = setup)) +
#   geom_point(aes(y = obs_tspan, 
#                  color = "obs",
#   ),
#   fill = "black",
#   size = 2,
#   shape = 21) +
#   scale_color_manual(name = "",
#                      breaks = br,
#                      labels = lb,
#                      values = cl,
#                      guide = guide_legend(override.aes = list(shape  = sh,
#                                                               size   = sz,
#                                                               linetype = lt,
#                                                               fill  = fl,
#                                                               color = cl))) +
#   guides(fill = "none") +
#   theme_linedraw() +
#   facet_wrap(~site_year, scales = "free_x", ncol = 1) +
#   ylim(0, 40) +
#   labs(y = bquote("Temperature [°C]"),
#        x = "Date") 
# 
# ggsave(paste0(dir_figs, "/noacc-fullacc_tspan.png"), p_tspan_season, width = 12, height = 12)
# 
# p_all <- 
#   (p_topt_season + p_aopt_season + p_tspan_season) +
#   plot_layout(ncol = 3,
#               guides = "collect") &
#   theme(legend.position = "bottom")
# 
# ggsave(paste0(dir_figs, "/noacc-fullacc_ALL.png"), p_all, width = 12, height = 12)
# 
# # ├—— all setups at once ----
# df_all <- 
#   rbind(df_er,
#         df_fullacc,
#         df_noacc,
#         df_sb,
#         df_pc) %>%
#   left_join(df_inst_fullacc %>% 
#               dplyr::select(sitename, date, forcing_growth, forcing_d, site_year))
# 
# # Aesthetics
# # br <- c("er", "sb", "pc", "obs")
# # lb <- c("ER", "SB", "PC", "Observ.")
# br <- c("noacc", "fullacc", "er", "sb", "pc", "obs")
# lb <- c("No Acc.", "Full Acc.", "ER", "SB", "PC", "Observ.")
# 
# n_tmp <- length(br) - 1
# 
# sh <- c(rep(NA, n_tmp), 21)
# sz <- c(rep(1, n_tmp), 2)
# lt <- c(rep(1, n_tmp), 0)
# fl <- c(rep(NA, n_tmp), "black")
# cl <- c(vec_cols[1:n_tmp], "black")
# 
# # Try to get a facet_grid (needs to much effort right now)
# 
# # T_opt
# p_topt_season <- 
#   df_all %>% 
#   ggplot() +
#   aes(x = date) +
#   geom_line(aes(y = mod_tc_opt, group = setup, color = setup)) +
#   geom_errorbar(aes(y = obs_tc_opt, 
#                       ymin = obs_tc_opt - obs_tc_opt_se, 
#                       ymax = obs_tc_opt + obs_tc_opt_se,
#                       color = "obs")) +
#   geom_point(aes(y = obs_tc_opt, 
#                  color = "obs",
#                  ),
#              fill = "black",
#              size = 2,
#              shape = 21) +
#   scale_color_manual(name = "",
#                      breaks = br,
#                      labels = lb,
#                      values = cl,
#                      guide = guide_legend(override.aes = list(shape  = sh,
#                                                               size   = sz,
#                                                               linetype = lt,
#                                                               fill  = fl,
#                                                               color = cl))) +
#   guides(fill = "none") +
#   theme_linedraw() +
#   facet_wrap(~site_year, scales = "free_x", ncol = 1) +
#   ylim(0, 40) +
#   labs(y = bquote("Temperature [°C]"),
#        x = "Date") 
# 
# ggsave(paste0(dir_figs, "/all_topt.png"), p_topt_season, width = 12, height = 12)
# 
# ## A_opt
# p_aopt_season <- 
#   df_all %>% 
#   ggplot() +
#   aes(x = date) +
#   geom_line(aes(y = mod_anet_opt, group = setup, color = setup)) +
#   geom_errorbar(aes(y = obs_tc_opt, 
#                     ymin = obs_anet_opt - obs_anet_opt_se, 
#                     ymax = obs_anet_opt + obs_anet_opt_se,
#                     color = "obs")) +
#   geom_point(aes(y = obs_anet_opt, 
#                  color = "obs",
#   ),
#   fill = "black",
#   size = 2,
#   shape = 21) +
#   scale_color_manual(name = "",
#                      breaks = br,
#                      labels = lb,
#                      values = cl,
#                      guide = guide_legend(override.aes = list(shape  = sh,
#                                                               size   = sz,
#                                                               linetype = lt,
#                                                               fill  = fl,
#                                                               color = cl))) +
#   guides(fill = "none") +
#   theme_linedraw() +
#   facet_wrap(~site_year, scales = "free_x", ncol = 1) +
#   ylim(0, 40) +
#   labs(y = bquote(A[opt] ~ "[µmol m¯² s¯¹]"),
#        x = "Date") 
# 
# ggsave(paste0(dir_figs, "/all_aopt.png"), p_aopt_season, width = 12, height = 12)
# 
# ## T_span
# p_tspan_season <- 
#   df_all %>% 
#   ggplot() +
#   aes(x = date) +
#   geom_line(aes(y = mod_tspan, group = setup, color = setup)) +
#   geom_point(aes(y = obs_tspan, 
#                  color = "obs",
#   ),
#   fill = "black",
#   size = 2,
#   shape = 21) +
#   scale_color_manual(name = "",
#                      breaks = br,
#                      labels = lb,
#                      values = cl,
#                      guide = guide_legend(override.aes = list(shape  = sh,
#                                                               size   = sz,
#                                                               linetype = lt,
#                                                               fill  = fl,
#                                                               color = cl))) +
#   guides(fill = "none") +
#   theme_linedraw() +
#   facet_wrap(~site_year, scales = "free_x", ncol = 1) +
#   ylim(0, 40) +
#   labs(y = bquote("Temperature [°C]"),
#        x = "Date") 
# 
# ggsave(paste0(dir_figs, "/all_tspan.png"), p_tspan_season, width = 12, height = 12)
# 
# p_all <- 
#   (p_topt_season + p_aopt_season + p_tspan_season) +
#   plot_layout(ncol = 3,
#               guides = "collect") &
#   theme(legend.position = "bottom")
# 
# ggsave(paste0(dir_figs, "/all_ALL.png"), p_all, width = 12, height = 12)

# ├—— Full Acc and Growth Climate ----

# var_selection <- c("mod_tc_opt", "mod_anet_opt", "mod_tspan", "vpd", 'temp')
var_selection <- c("mod_tc_opt", "vpd", 'temp')

df_tmp <- 
  df_fullacc %>% 
  left_join(df_inst_fullacc %>% 
              dplyr::select(sitename, date, forcing_growth, forcing_d, site_year)) %>% 
  unnest(forcing_growth) %>% 
  mutate(mod_tc_opt = mod_tc_opt/max(mod_tc_opt),
         mod_anet_opt = mod_anet_opt/max(mod_anet_opt),
         mod_tspan = mod_tspan/max(mod_tspan),
         obs_tc_opt = obs_tc_opt/max(obs_tc_opt),
         obs_tc_opt_se = obs_tc_opt_se/max(obs_tc_opt),
         vpd = vpd/max(vpd),
         temp = temp/max(temp)) %>% 
  pivot_longer(cols = var_selection,
             values_to = "val",
             names_to  = "nam")
  
# p_topt_climate <-
  df_tmp %>% 
  ggplot() +
  aes(x = date) +
  geom_line(aes(y = val, group = nam, color = nam)) +
  # geom_errorbar(aes(y = obs_tc_opt, 
  #                   ymin = obs_tc_opt - obs_tc_opt_se, 
  #                   ymax = obs_tc_opt + obs_tc_opt_se,
  #                   color = "obs")) +
  # geom_point(aes(y = obs_tc_opt, 
  #                color = "obs",
  # ),
  # fill = "black",
  # size = 2,
  # shape = 21) +
  # scale_color_manual(name = "",
  #                    breaks = br,
  #                    labels = lb,
  #                    values = cl,
  #                    guide = guide_legend(override.aes = list(shape  = sh,
  #                                                             size   = sz,
  #                                                             linetype = lt,
  #                                                             fill  = fl,
  #                                                             color = cl))) +
  guides(fill = "none") +
  theme_linedraw() +
  facet_wrap(~site_year, scales = "free_x", ncol = 1) +
  ylim(0, 1) +
  labs(y = bquote("Scaled Variable"),
       x = "Date") 

df_tmp %>% 
  ggplot() +
  aes(x = date) +
  geom_line(aes(y = val, group = nam, color = nam)) +
  # scale_color_manual(name = "",
  #                    breaks = br,
  #                    labels = lb,
  #                    values = cl,
  #                    guide = guide_legend(override.aes = list(shape  = sh,
  #                                                             size   = sz,
  #                                                             linetype = lt,
  #                                                             fill  = fl,
  #                                                             color = cl))) +
  guides(fill = "none") +
  theme_linedraw() +
  facet_wrap(~site_year, scales = "free_x", ncol = 1) +
  ylim(0, 1) +
  labs(y = bquote("Scaled Variable"),
       x = "Date") 


# ├—— Facet grid (Full-No, all Traits) ----
df_obs_1 <-
  df_fullacc %>% 
  drop_na(obs_tc_opt) %>% 
  mutate(setup = "obs") %>% 
  pivot_longer(c(obs_tc_opt, obs_anet_opt, obs_tspan),
               names_to = "trait",
               values_to = "value") %>% 
  mutate(trait = str_remove_all(trait, "obs_"),
         trait = ifelse(str_detect(trait, "anet"),
                        paste0("aopt"),
                        trait),
         trait = ifelse(str_detect(trait, "tc_opt"),
                        paste0("topt"),
                        trait)) %>% 
  dplyr::select(sitename, date, trait, value, setup, site_year)
  
df_obs_2 <-
  df_fullacc %>% 
  drop_na(obs_tc_opt) %>% 
  mutate(setup = "obs",
         obs_tspan_se = 0) %>% 
  pivot_longer(c(obs_tc_opt_se, obs_anet_opt_se, obs_tspan_se),
               names_to = "trait",
               values_to = "se") %>% 
  mutate(trait = str_remove_all(trait, "obs_|_se"),
         trait = ifelse(str_detect(trait, "anet"),
                        paste0("aopt"),
                        trait),
         trait = ifelse(str_detect(trait, "tc_opt"),
                        paste0("topt"),
                        trait)) %>% 
  dplyr::select(sitename, date, trait, se, setup, site_year)


df_obs <- left_join(df_obs_1, df_obs_2)


df_mod <- 
  df_fullacc %>% 
  rbind(df_fullacc, 
        df_noacc) %>%  
  pivot_longer(c(mod_tc_opt, mod_anet_opt, mod_tspan),
               names_to = "trait",
               values_to = "value") %>% 
  mutate(trait = str_remove_all(trait, "mod_"),
         trait = ifelse(str_detect(trait, "anet"),
                        paste0("aopt"),
                        trait),
         trait = ifelse(str_detect(trait, "tc_opt"),
                        paste0("topt"),
                        trait),
         se    = 0) %>% 
  dplyr::select(sitename, date, trait, value, se, setup, site_year)

df_long <- 
  rbind(df_mod, df_obs) %>% 
  mutate(doy = lubridate::yday(date))

p_facet <- 
  ggplot() +
  geom_line(data = df_long %>% dplyr::filter(setup != "obs"),
            aes(x = doy,
                y = value,
                group = setup,
                color = setup)) +
  geom_point(data = df_long %>% dplyr::filter(setup == "obs"),
             aes(x = doy,
                 y = value,
                 group = setup,
                 color = setup)) +
  geom_errorbar(data = df_long %>% dplyr::filter(setup == "obs"),
           aes(x = doy,
               ymin = value - se,
               ymax = value + se,
               group = setup,
               color = setup)) +
  # facet_wrap(~site_year + trait, 
  #            scales = "free_x",
  #            ncol = 3,
  #           labeller = as_labeller(
  #             c(
  #               'topt'  = "T[opt] ~ (degC)",
  #               'aopt'  = "A[opt] ~ (µmol/m^2/s)",
  #               'tspan' = "T[span] ~ (degC)",
  #               'jensen2011' = "Jensen~et~al.~(2015)",
  #               'hikosaka2001' = "Hikosaka~et~al.~(2007)",
  #               'hikosaka2002' = "Hikosaka~et~al.~(2007)",
  #               'slot2016' = "Slot~et~al.~(2017)",
  #               'medlyn2002' = "Medlyn~et~al.~(2007)"),
  #             label_parsed)
  #           ) +
  facet_grid(vars(site_year),
             vars(trait),
             scales = "free",
             as.table = T,
             # switch = "x",
             labeller = as_labeller(
               c(
                 'topt'  = "T[opt] ~ (degree*C)",
                 'aopt'  = "A[opt] ~ (µmol/m^2/s)",
                 'tspan' = "T[span] ~ (degree*C)",
                 'cavaleri2014' = "Cavaleri~et~al.~(2014~unpub.)",
                 'rogers2012' = "Rogers~et~al.~(2012)",
                 'jensen2011' = "Jensen~et~al.~(2015)",
                 'hikosaka2001' = "Hikosaka~et~al.~(2007)",
                 'hikosaka2002' = "Hikosaka~et~al.~(2007)",
                 'slot2016' = "Slot~et~al.~(2017)",
                 'medlyn2000' = "Medlyn~et~al.~(2002)"),
                 # 'medlyn2000' = "Medlyn~et~al.~(2002)~(2)",
                 # 'medlyn1999' = "Medlyn~et~al.~(2002)~(1)"),
               label_parsed)
             ) +
  ylim(0, 40) +
  # coord_flip() +
  scale_color_brewer(palette = "Dark2",
                     name = "Source:",
                     labels = c("Model (Full Acc.)",
                                "Model (No Acc.)",
                                "Obs. (with standard error)")) +
  theme_linedraw() +
  theme(legend.position = "bottom",
        # axis.title = element_text(size = 12),
        # strip.text = element_text(size = 10),
        text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(title = NULL, #expression(bold("Seasonal Acclimation of Traits")),
       x = "Day of Year",
       y = "Trait Value") 


ggsave(paste0(dir_figs, "seasonality_noacc-fullacc-facet.pdf"), p_facet, width = 8, height = 11)

beepr::beep()
stop("End Of Script")

# ├— Plots (other checks) ------------------------------------------

# ├—— Forcings ----
df_inst_fullacc %>%
  unnest(forcing_d) %>% 
  ggplot() +
  aes(x = date, y = vpd) +
  geom_line() +
  facet_wrap(~site_year, scales = "free_x", ncol = 1) +
  theme_linedraw() 

df_inst_fullacc %>%
  unnest(forcing_d) %>% 
  ggplot() +
  aes(x = date, y = patm) +
  geom_line() +
  facet_wrap(~site_year, scales = "free_x", ncol = 1) +
  theme_linedraw() 

df_inst_fullacc %>%
  rename(fd = forcing_d, fg = forcing_growth) %>% 
  unnest(c(fd, fg), names_sep = "_") %>% 
  pivot_longer(cols = c("fd_vpd", "fg_temp"), values_to = "climate_val", names_to = "climate_var") %>% 
  ggplot() +
  aes(x = date, y = climate_val) +
  geom_line() +
  # facet_wrap(~site_year, scales = "free_x", ncol = 1) +
  facet_grid(vars(climate_var), vars(site_year), scales = "free") +
  theme_linedraw() 

climate_vars <- c("co2", "patm", "temp", 'vpd')

df_inst_fullacc %>% 
  unnest(forcing_growth) %>% 
  pivot_longer(cols = climate_vars, 
               values_to = "climate_val", 
               names_to = "climate_var") %>% 
  ggplot() +
  aes(x = date, y = climate_val, color = climate_var) +
  geom_line() +
  # facet_wrap(~site_year, scales = "free", ncol = 1) +
  facet_grid(vars(climate_var), vars(site_year), scales = "free") +
  theme_linedraw() 
    
df_inst_fullacc %>% 
  unnest(forcing_d) %>% 
  pivot_longer(cols = climate_vars, 
               values_to = "climate_val", 
               names_to = "climate_var") %>% 
  ggplot() +
  aes(x = date, y = climate_val, color = climate_var) +
  geom_line() +
  # facet_wrap(~site_year, scales = "free", ncol = 1) +
  facet_grid(vars(climate_var), vars(site_year), scales = "free") +
  theme_linedraw() 


# ├—— Acclimated Variables ----------------------------------------------------------------------

acc_vars <- c("vcmax", "vcmax25", 
              "jmax", "jmax25", 
              "rd25", "xi", "kphio")

## Unscaled 
df_inst_fullacc %>% 
  unnest(rpm_acc) %>% 
  pivot_longer(cols = acc_vars, 
               values_to = "val", 
               names_to = "var") %>% 
  rowwise() %>% 
  mutate(val = ifelse(var %in% c("xi", "kphio"), val, val*1e6)) %>% 
  ungroup() %>% 
  ggplot() +
  aes(x = date, y = val, color = var) +
  geom_line() +
  # facet_wrap(~site_year, scales = "free", ncol = 1) +
  facet_grid(vars(var), vars(site_year), scales = "free", switch = "y") +
  theme_linedraw() 

## Focus on kphio
acc_vars <- c("anet_opt", "kphio", "temp", "vcmax25")

df_inst_fullacc %>%
  unnest(rpm_acc, rpm_sim, forcing_growth) %>% 
  mutate(anet_opt = anet_opt / max(anet_opt),
         kphio = kphio / max(kphio),
         temp  = temp / max(temp),
         vcmax25 = vcmax25 / max(vcmax25)) %>% 
  pivot_longer(cols = acc_vars, 
               values_to = "val", 
               names_to = "var") %>% 
  ggplot() +
  aes(x = date, y = val, color = var) +
  geom_line() +
  facet_wrap(~site_year, scales = "free", ncol = 1) +
  # facet_grid(vars(var), vars(site_year), scales = "free", switch = "y") +
  theme_linedraw() 


# ├—— Simulated Variables ----
sim_vars <- c("tc_opt_ac", "tc_opt_aj", "tc_opt", "tspan")

df_inst_fullacc %>% 
  unnest(rpm_sim) %>% 
  pivot_longer(cols = sim_vars, 
               values_to = "val", 
               names_to = "var") %>% 
  rowwise() %>% 
  mutate(val = ifelse(var == "anet_opt", val*1e6, val)) %>% 
  ungroup() %>% 
  ggplot() +
  aes(x = date, y = val, color = var) +
  geom_line() +
  facet_wrap(~site_year, scales = "free", ncol = 1) +
  # facet_grid(vars(var), vars(site_year), scales = "free", switch = "y") +
  theme_linedraw() 

# ├—— Instantaneous Variables ----
inst_vars <- c("ac", "aj", "rd", "chi", "gs_c")

df_inst_fullacc %>% 
  unnest(c(rpm_inst, forcing_d)) %>% 
  mutate(chi = ci/co2_to_ca(co2, patm)) %>% 
  pivot_longer(cols = inst_vars, 
               values_to = "val", 
               names_to = "var") %>%
  dplyr::filter(tc_leaf == 25) %>% 
  rowwise() %>% 
  mutate(val = ifelse(var == "chi", val, val*1e6)) %>% 
  ungroup() %>% 
  ggplot() +
  aes(x = date, y = val, color = var) +
  geom_line() +
  # facet_wrap(~site_year, scales = "free", ncol = 1) +
  facet_grid(vars(var), vars(site_year), scales = "free", switch = "y") +
  theme_linedraw() 
