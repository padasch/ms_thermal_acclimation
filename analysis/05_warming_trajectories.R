# TRAJECTORIES BETWEEN DIFFERENT SITES ----
## Settings ----
dir_analysis <- "trajectories"
dir_tmp <- NA
use_pfts <- T

## Define directories ----
# First check how many runs have been made to add counter
dir_mods <- make_new_dynamic_dir(dir_analysis, "model_rds", dir_tmp)
dir_figs <- make_new_dynamic_dir(dir_analysis, "figures", dir_tmp)
dir_tabs <- make_new_dynamic_dir(dir_analysis, "tables", dir_tmp)

## Get Forcing Data ----
### Reference ----

# Get site-date df
settings <- get_settings()
settings$save_plots <- F
settings$daily_conditions <- "B"
settings$tau <- 15
settings$rpmodel_exp <- 2000e-6

df_b15           <- k19_create_df_forcing(settings)
df_site_date_b15 <- k19_create_input_for_acc(settings, df_b15)

# Select sites of interest, one from each climate zone
df_tmp_0 <- 
  df_site_date_b15 %>% 
  mutate(cl_ = purrr::map_chr(siteinfo, ~pull(., cl) %>% substr(., 1, 1))) %>% 
  group_split(cl_)
  
df_tmp_1 <- tibble()
set.seed(2023)

for (i in 1:length(df_tmp_0)) {
  df_tmp_1 <- 
    rbind(
      df_tmp_1,
      df_tmp_0[[i]] %>% 
        slice_sample(n = 1) %>% 
        dplyr::select(-cl_)
      )
}

### Warming ----
# Add warming to tc_growth and vpd_growth
max_warming   <- 5 
steps_warming <- 5
df_tmp_3      <- tibble()
nest_forcing  <- names(df_tmp_1$forcing_growth[[1]]) 
nest_siteinfo <- c(names(df_tmp_1$siteinfo[[1]]), "warming", "site_id")

for (i in 1:nrow(df_tmp_1)) {
  cat("\n Site: ", df_tmp_1$sitename[[i]])
  
  for (add_degc in seq(0, max_warming, max_warming/steps_warming)) {
    
    # unnest forcing and siteinfo
    df_tmp_2 <- df_tmp_1 %>% slice(i) %>% unnest(forcing_growth, siteinfo)
    
    # add sitename and warming to siteinfo
    df_tmp_2$warming <- add_degc
    df_tmp_2$site_id <- df_tmp_2$sitename

    # Make sitename unique
    df_tmp_2$sitename <- paste0(df_tmp_2$sitename, "_+", add_degc)
    
    # Add warming to temp and vpd
    old_tmp <- df_tmp_2$temp
    df_tmp_2$temp <- old_tmp + add_degc
    
    old_vpd <- df_tmp_2$vpd
    df_tmp_2$vpd <- VPDairToLeaf(old_vpd/1000, old_tmp, df_tmp_2$temp) * 1000 # Toggle on/off
    
    # Verbose output
    cat("\n   Temp: from ", old_tmp, "to", df_tmp_2$temp, 
        " --- VPD: from", old_vpd, " to ", df_tmp_2$vpd)
    
    # Nest dataframes
    df_tmp_2 <- df_tmp_2 %>% nest(forcing_growth = all_of(nest_forcing))
    df_tmp_2 <- df_tmp_2 %>% nest(siteinfo = all_of(nest_siteinfo))
    
    # Overwrite daily forcing with growth forcing, to keep consistent warming
    df_tmp_2$forincg_d <- df_tmp_2$forcing_growth
    
    # Bind dataframes
    df_tmp_3 <- rbind(df_tmp_3, df_tmp_2)
  }
}

## Running Models ----
### Full Accl. ----
df_acc   <- k19_run_acc(df_tmp_3, settings)
df_inst_fullacc  <- k19_run_inst(df_acc, settings)
df_final_fullacc <- k19_create_final_df(df_inst_fullacc, settings)
saveRDS(list(df = df_final_fullacc, set = settings), here(dir_mods, "df_final_fullacc.rds"))

### No Acclimation ----
settings$method_ftemp <- "leuning02"
settings <- k19_from_forc_to_plot(settings, returndf = "update_settings")
df_acc   <- k19_run_acc(df_tmp_3, settings)

if (use_pfts) df_acc <- replace_acclimated_with_pfts(df_acc, 
                                                     replace_pc = T, 
                                                     replace_xi = T)

df_inst_noacc    <- k19_run_inst(df_acc, settings)
df_final_noacc <- k19_create_final_df(df_inst_noacc, settings)
saveRDS(list(df = df_final_noacc, set = settings), here(dir_mods, "df_final_noacc.rds"))

### ER ----
settings$method_ftemp <- "kumarathunge19"
settings <- k19_from_forc_to_plot(settings, returndf = "update_settings")

df_acc   <- k19_run_acc(df_tmp_3, settings)

if (use_pfts) df_acc <- replace_acclimated_with_pfts(df_acc, 
                                                     replace_pc = T, 
                                                     replace_xi = T)

df_inst_er    <- k19_run_inst(df_acc, settings)
df_final_er <- k19_create_final_df(df_inst_er, settings)
saveRDS(list(df = df_final_er, set = settings), here(dir_mods, "df_final_er.rds"))

## Load RDS ----
all_files  <- list.files(here(dir_mods), pattern = ".rds")
plot_names <- paste0("p_", str_remove_all(all_files, ".rds|df_final_"))
out_names  <- paste0("out_", str_remove_all(all_files, ".rds|df_final_"))

for (i in 1:length(all_files)) {
  tmp    <- readRDS(here(dir_mods, all_files[i]))
  assign(out_names[i],  tmp)
  # assign(plot_names[i], plot_all_final_plots_from_df_plot(tmp$df, tmp$set))
}

## Wrangling dataframes ----
# We want to have a fullacc model that omits temporal acclimation
# Replace the predicted response curves from the non-warmed df with the warmed dfs

df_fullacc_spatial <- 
  out_fullacc$df %>%
  nest(data = !any_of(c("sitename", "temp", "warming")))

overwrite <- TRUE

for (i in 1:nrow(df_fullacc_spatial)) {
  
  # Get df to overwrite subsequent ones with
  if (overwrite) {
    keep_df <- df_fullacc_spatial$data[[i]]
    overwrite <- FALSE
  }
  
  # Overwrite df with old one
  df_fullacc_spatial$data[[i]] <- keep_df
  
  # Check if next site is up, get new df
  if (i %% (steps_warming + 1) == 0) overwrite <- TRUE
}

df_fullacc_spatial <- 
  df_fullacc_spatial %>% 
  unnest(c("data"))

# Repeat the same for the enzymatic response df
df_er_spatial <- 
  out_er$df %>%
  nest(data = !any_of(c("sitename", "temp", "warming")))

overwrite <- TRUE

for (i in 1:nrow(df_er_spatial)) {
  
  # Get df to overwrite subsequent ones with
  if (overwrite) {
    keep_df <- df_er_spatial$data[[i]]
    overwrite <- FALSE
  }
  
  # Overwrite df with old one
  df_er_spatial$data[[i]] <- keep_df
  
  # Check if next site is up, get new df
  if (i %% (steps_warming + 1) == 0) overwrite <- TRUE
}

df_er_spatial <- 
  df_er_spatial %>% 
  unnest(c("data"))

## Make dfs longer
df_obs     <- turn_dfplot_into_dfdelta(out_fullacc$df, setup = "obs")
df_noacc   <- turn_dfplot_into_dfdelta(out_noacc$df, setup = "noacc")
df_er_spatial      <- turn_dfplot_into_dfdelta(df_er_spatial, setup = "er_spatial")
df_fullacc_spatial <- turn_dfplot_into_dfdelta(df_fullacc_spatial, setup = "fullacc_spatial")
df_fullacc_temporal <- turn_dfplot_into_dfdelta(out_fullacc$df, setup = "fullacc_temporal")

## Get one long dataframe
df_long <- 
  rbind(df_fullacc_temporal,
        df_fullacc_spatial,
        df_noacc,
        df_er_spatial)

tab_traj_agrowth <- 
  df_long %>% 
  dplyr::filter(warming %in% c(0, 5)) %>% 
  group_by(site_id, setup) %>% 
  nest() %>% 
  mutate(anetperc_0 = purrr::map_dbl(
             data,
             ~dplyr::filter(., warming == 0) %>% pull(rel_red_aopt)
           ),
         anetperc_5 = purrr::map_dbl(
           data,
           ~dplyr::filter(., warming == 5) %>% pull(rel_red_aopt)
           ),
         change_realized_anet = round((anetperc_5 - anetperc_0) * 100),
         setup = as.character(setup),
         setup = ifelse(str_detect(setup, "fullacc_temporal"), "Acclimation", setup),
         setup = ifelse(str_detect(setup, "fullacc_spatial"), "Adaptation", setup),
         setup = ifelse(str_detect(setup, "noacc"), "Fixed", setup),
         setup = ifelse(str_detect(setup, "er_spatial"), "Adaptation (ER only)", setup),
         site_id = ifelse(str_detect(site_id, "slot"), "Tropical", site_id),
         site_id = ifelse(str_detect(site_id, "tarv"), "Temperate", site_id),
         site_id = ifelse(str_detect(site_id, "jens"), "Boreal", site_id),
         site_id = ifelse(str_detect(site_id, "roge"), "Arctic", site_id),
         ecosystem = site_id) %>% 
  ungroup() %>% 
  dplyr::select(ecosystem, setup, change_realized_anet) %>% 
  arrange(ecosystem, setup)

tab_traj_agrowth %>% knitr::kable()

readr::write_csv(tab_traj_agrowth, paste(dir_tabs, "changes_in_realized_anet.csv"))

## Plotting ----

## Get one long dataframe
df_long <- 
  rbind(df_fullacc_temporal,
        df_fullacc_spatial,
        df_noacc) %>% 
  mutate(
    across(
      setup, factor, 
      levels=c("noacc", "fullacc_spatial", "fullacc_temporal"))) %>% 
  dplyr::filter(
    # warming %in% c(0, 2, 4)
    # warming %in% c(0, 5)
  )

### Fixed - Spatial - Temporal ----
## Presave scales
p_facet <- 
  facet_wrap(~setup,
             nrow = 1,
             labeller = as_labeller(c(
               'noacc'  = "Fixed",
               'fullacc_spatial'   = "Adaptation",
               'fullacc_temporal'  = "Acclimation"
             )))

p_fill <- 
  scale_fill_viridis_d(
    "Warming [°C]: ",
    option = "viridis",
    # labels = c("+0", "+5"),
    labels = c("+0", "+1", "+2", "+3", "+4", "+5"),
    guide = guide_legend(
      override.aes = list(
        # shape = c(21, 21)
        shape = c(21, 21, 21, 21, 21, 21)
      )
    )
  )

p_shape <- 
  scale_shape_manual(
    "Ecosystem: ",
    labels = c("Arctic", "Temperate", "Boreal", "Tropical"),
    values = c(22, 21, 25, 24),
    guide = guide_legend(
      nrow   = 2,
      byrow  = TRUE,
      override.aes = list(
        shape = c(21, 24, 22, 25),
        fill  = c("white", "white", "white", "white")))
  )

p_theme <- 
  theme_linedraw(base_size = 14) +
  theme(
    # legend.position = c(0.9, 0.05),
    # legend.justification = c(1, 0),
    # legend.position = "bottom",
    # strip.background =element_rect(fill="white"),
    # strip.text = element_text(colour = 'black', face = "bold"),
    panel.grid = element_blank(),
    plot.subtitle = element_text(hjust = 0, face = "bold"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14))
  
## Absolute Change
p_absolute <-
  df_long %>% 
  ggplot() +
  
  # > Geoms
  aes(
    # y = rel_red_aopt*100,
    y = agrowth,
    # x = delta_t
    x = temp
  ) +
  geom_path(
    aes(
      group = site_id
    )
  ) +
  geom_point(
    aes(
      fill = as.factor(warming),
      shape = site_id,
    ),
    size = 2.5) +
  
  # > Legend and colors
  p_shape +
  p_fill  +
  p_facet +
  p_theme +
  
  # > Layout
  # ylim(40, 100) +
  # xlim(-15, 15) +
  labs(y = expression(A[net] ~ "at" ~ T[growth] ~ "[µmol" ~ m ^-2 ~ s ^-1 ~ "]"),
       # x = expression(T[opt] ~ "-" ~ T[growth] ~ "[°C]")
       x = expression(T[growth] ~ "[°C]"),
       subtitle = "A.)"
       # subtitle = expression(bold("Trajectories of" ~ A[net] ~ "under warming"))
       )


## Relative change to Aopt
p_relative <-
  df_long %>% 
  ggplot() +
  
  # > Geoms
  aes(
    y = rel_red_aopt*100,
    # y = agrowth,
    # x = delta_t
    x = temp
  ) +
  geom_path(
    aes(
      group = site_id
    )
  ) +
  geom_point(
    aes(
      fill = as.factor(warming),
      shape = site_id,
    ),
    size = 2.5) +

  # > Legend and colors
  p_shape +
  p_fill  +
  p_facet +
  p_theme +
  # ylim(40, 100) +
  # xlim(-15, 15) +
  labs(y = expression(A[net] ~ "at" ~ T[growth] ~ "[% of"~ A[opt] ~ "]"),
       # x = expression(T[opt] ~ "-" ~ T[growth] ~ "[°C]")
       x = expression(T[growth] ~ "[°C]"),
       subtitle = "B.)"
       # subtitle = expression(bold("Trajectories of" ~ A[net] ~ "under warming"))
       )


## Saving
p_save <- 
  p_absolute / p_relative +
    plot_layout(
      guides = "collect"
    ) +
    plot_annotation(
      title = expression(bold("Trajectories of" ~ A[net] ~ "under warming"))
    ) &
    theme(legend.position = "bottom")
# 
# ggsave(paste0(dir_figs, "/trajectory-absolute-and-relative.pdf"),
#        p_save,
#        height = 10,
#        width = 10)

# Relative scale only
p_save <-
  p_relative + 
  labs(title = expression(bold("Trajectories of" ~ A[net] ~ "under warming")),
       subtitle = NULL) +
  theme(plot.title = element_text(hjust = 0)) + 
  # Add tags
  geom_text(
    data = 
      data.frame(
        setup   = c('noacc', 'fullacc_spatial', 'fullacc_temporal'),
        label   = c('(a)', '(b)', '(c)')
      ) %>% 
      mutate(across(setup, factor,levels=c("noacc", "fullacc_spatial", "fullacc_temporal"))),
    aes(y = 32, 
        x = 10, 
        label = label),
    fontface = 'bold', 
    size = 4) +
  ylim(30, 100)

ggsave(paste0(dir_figs, "/trajectory-relative-only.pdf"),
       p_save,
       height = 4.25,
       width = 10)

# Absolute scale only
p_save <- 
  p_absolute + 
  labs(title = expression(bold("Trajectories of" ~ A[net] ~ "under warming")),
       subtitle = NULL) +
  theme(plot.title = element_text(hjust = 0))


# ggsave(paste0(dir_figs, "/trajectory-absolute-only.pdf"),
#        p_save,
#        height = 4.25,
#        width = 10)

### Noacc - ER - Full ----

## Bind dataframes to one
# df_long <- 
#   rbind(df_fullacc_temporal, 
#         df_noacc,
#         df_er) %>% 
#   mutate(
#     across(
#       setup, factor, 
#       levels=c("noacc", "er", "fullacc_temporal")))
# 
# p <-
#   df_long %>% 
#   ggplot() +
#   
#   # > Geoms
#   aes(
#     # y = rel_red_aopt*100,
#     y = agrowth,
#     # x = delta_t
#     x = temp
#   ) +
#   geom_path(
#     aes(
#       group = site_id
#     )
#   ) +
#   geom_point(
#     aes(
#       fill = as.factor(warming),
#       shape = site_id,
#     ),
#     size = 2.5) +
#   facet_wrap(~setup,
#              nrow = 1,
#              labeller = as_labeller(c(
#                'fullacc_temporal'  = "Full Acclimation",
#                'noacc'  = "No Acclimation",
#                'er' = "ER"
#              ))) +
# 
#   # > Legend and colors
#   theme_linedraw() +
#   theme(
#     # legend.position = c(0.9, 0.05),
#     # legend.justification = c(1, 0),
#     # legend.position = "bottom",
#     # strip.background =element_rect(fill="white"),
#     # strip.text = element_text(colour = 'black', face = "bold"),
#     panel.grid = element_blank(),
#     strip.text = element_text(face = "bold"),
#     plot.title = element_text(hjust = 0.5, size = 14)) +
#   scale_fill_viridis_d(
#     "Warming [°C]: ",
#     option = "viridis",
#     labels = c("+0", "+1", "+2", "+3", "+4", "+5"),
#     guide = guide_legend(
#       override.aes = list(
#         shape = c(21, 21, 21, 21, 21, 21)
#       )
#     )
#   ) +
#   scale_shape_manual(
#     "Climate Zone: ",
#     labels = c("Arctic", "Temperate", "Boreal", "Tropical"),
#     values = c(22, 21, 25, 24),
#     guide = guide_legend(
#       override.aes = list(
#         shape = c(21, 24, 22, 25),
#         fill  = c("black", "black", "black", "black")))
#     ) +
#   
#   # > Layout
#   # ylim(40, 100) +
#   # xlim(-15, 15) +
#   labs(y = expression(A[net] ~ "at" ~ T[growth] ~ " as percentage of " ~ A[opt] ~ "[%]"),
#        x = expression(T[opt] ~ "-" ~ T[growth] ~ "[°C]"),
#        # x = expression(T[growth] ~ "[°C]"),
#        subtitle = expression(bold("Trajectories of" ~ A[net] ~ "under warming")))
# 
# ggsave(paste0(dir_figs, "/trajectory-agrowth_vs_deltat.pdf"),
#        p,
#        height = 4,
#        width = 12)
# 
# df_tmp_4 <- df_long %>% group_split(setup)
# 
# for (i in 1:length(df_tmp_4)) {
#   
#   cat("\n-----------------\n",
#       "Setup: ", as.character(df_tmp_4[[i]]$setup)[1])
#   
#   df_tmp_5 <- df_tmp_4[[i]] %>% group_split(site_id)
#   
#   for (j in 1:length(df_tmp_5)) {
#     
#     cat("\n Site: ", df_tmp_5[[j]]$site_id[[1]],
#         "    d Anet = ", round(abs(max(df_tmp_5[[j]]$rel_red_aopt - min(df_tmp_5[[j]]$rel_red_aopt))), 2))
#   }
# }
