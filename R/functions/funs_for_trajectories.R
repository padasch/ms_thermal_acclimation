# Final Script ----

turn_dfplot_into_dfdelta <- function(df_plot,
                                     setup) {
  
  df_out <- 
    df_plot %>% 
    unnest(rpm_inst) %>% 
    mutate(temp_rnd = round(temp),
           id       = paste(sitename, date)) %>%
    dplyr::filter(temp_rnd == tc_leaf) %>% 
    mutate(delta_t     = temp - rpm_sim_tc_opt,
           delta_t_abs = abs(delta_t),
           anet        = anet*1e6,
           delta_a     = anet - rpm_sim_anet_opt,
           delta_a_abs = abs(delta_a),
           agrowth     = anet,
           rel_red_aopt= 1 - (rpm_sim_anet_opt - anet) / rpm_sim_anet_opt,
           aopt        = rpm_sim_anet_opt,
           setup       = as.factor(setup))
  
  if (setup == "obs") {
    df_out <- 
      df_out %>% 
      mutate(delta_t = abs(eva_data_tc_opt - temp))
  }
  
  df_out <- 
    df_out %>% 
    dplyr::select(
      any_of(
        c(
          "id", 
          "delta_t",
          "delta_t_abs",
          "delta_a",
          "delta_a_abs",
          "rel_red_aopt",
          "aopt",
          "setup",
          "temp",
          "warming", 
          "site_id",
          "agrowth")
      )
      )
  
  return(df_out)
}

plot_agrowth_deltat <- function(out_fullacc,
                                out_noacc,
                                out_er) {

  df_obs     <- turn_dfplot_into_dfdelta(out_fullacc$df, setup = "obs")
  df_fullacc <- turn_dfplot_into_dfdelta(out_fullacc$df, setup = "fullacc")
  df_noacc   <- turn_dfplot_into_dfdelta(out_noacc$df, setup = "noacc")
  df_er      <- turn_dfplot_into_dfdelta(out_er$df, setup = "er")
  
  
  ## Bind dataframes to one
  df_long <- 
    rbind(df_fullacc, 
          df_noacc,
          df_er)
  
  ## Get coordinates for illustrative points
  id_t0_dt <- df_long %>% arrange(delta_t) %>% slice(1) %>% pull(id)
  df_t0_dt <- df_long %>% dplyr::filter(id == id_t0_dt)
  
  # Adjust t0 coordinate to fall on trajectory
  # noacc 
  df_t0_dt$rel_red_aopt[2] <- 0.555 # original 0.574
  df_t0_dt$delta_t[2]      <- -14.6  # original -14.9
  
  # fullacc
  df_t0_dt$rel_red_aopt[1] <- 0.905 # original 0.912
  df_t0_dt$delta_t[1]      <- -6.65  # original -6.44
  
  # er 
  # no adjustment needed, original values fit well
  
  # Adjust t1 coordinate to fall on trajectory
  t1_warming_degc <- 5
  df_t1_dt <- df_t0_dt
  
  # noacc 
  df_t1_dt$rel_red_aopt[2] <- 0.795 # original 0.574
  df_t1_dt$delta_t[2]      <- df_t0_dt$delta_t[2] + t1_warming_degc  
  
  # fullacc
  df_t1_dt$rel_red_aopt[1] <- 0.995 # original 0.912
  df_t1_dt$delta_t[1]      <- df_t0_dt$delta_t[1] + t1_warming_degc  
  
  # er 
  df_t1_dt$rel_red_aopt[3] <- 0.96 # original 0.725
  df_t1_dt$delta_t[3]      <- df_t0_dt$delta_t[3] + t1_warming_degc   
  
  ## Merging illustrative points into one df
  df_dt <- 
    rbind(
      df_t0_dt %>% dplyr::mutate(t0_or_t1 = "0"),
      df_t1_dt %>% dplyr::mutate(t0_or_t1 = "1")
    )
  
  ## Saving trajectory information on illustrative points
  df_traj <- data.frame(matrix(ncol = 3, nrow = 3))
  colnames(df_traj) <- c("setup", "change_in_tgrowth", "change_in_anet")
  
  df_traj$setup <- c("fullacc", "noacc", "er")
  df_traj$change_in_tgrowth <- 
    c(
      df_t1_dt$delta_t[1] - df_t0_dt$delta_t[1],
      df_t1_dt$delta_t[2] - df_t0_dt$delta_t[2],
      df_t1_dt$delta_t[3] - df_t0_dt$delta_t[3]
    )
  df_traj$change_in_anet <- 
    c(
      df_t1_dt$rel_red_aopt[1] - df_t0_dt$rel_red_aopt[1],
      df_t1_dt$rel_red_aopt[2] - df_t0_dt$rel_red_aopt[2],
      df_t1_dt$rel_red_aopt[3] - df_t0_dt$rel_red_aopt[3]
    )
  
  df_traj$da_dt <- df_traj$change_in_anet / df_traj$change_in_tgrowth
  
  ## Make Figure
  p_dt_tgrowth <-
    df_long %>% 
    ggplot() + 
    aes(delta_t, rel_red_aopt*100, color = setup) +
    scale_color_brewer("Model Setup: ", 
                       labels = c("Full Acclimation", "No Acclimation", "ER"), 
                       palette = "Dark2") +
    scale_fill_brewer("Model Setup: ",
                      labels = c("Full Acclimation", "No Acclimation", "ER"),
                      palette  = "Dark2",
                      guide = guide_legend(
                        override.aes = list(shape = c(NA, NA, NA))
                      )) +
    # geom_vline(xintercept = 0, linetype = "dotted") +
    geom_point(alpha = 0.25, shape = 4) +
    geom_smooth(se = FALSE) +
    geom_point(
      data = df_dt,
      aes(delta_t, 
          rel_red_aopt*100, 
          fill = setup,
          shape = t0_or_t1),
      color = "black",
      size  = 3
    ) +
    scale_shape_manual(
      expression("Ecosystem at" ~ T[growth] ~ ": "),
      labels = c(expression(T[0]), expression(T[0] + 5 ~ "°C")),
      values = c(21, 24),
      guide = guide_legend(
        override.aes = list(shape = c(16, 17))
        ),
    ) +
    # annotate("text", label = "Arctic",  x = -0.8*15,  y = 50 #, size = 8, colour = "red"
    # ) +
    # annotate("text", label = "Boreal",  x = -0.20*15,  y = 50 #, size = 8, colour = "red"
    # ) +
    # annotate("text", label = "Temperate",  x = 0.20*15,  y = 50 #, size = 8, colour = "red"
    # ) +
    # annotate("text", label = "Tropic",  x = 0.8*15,  y = 50 #, size = 8, colour = "red"
    # ) +
    # theme_bw() + # to include minor grid lines
    theme_classic() +
    theme(legend.position = "right",
          legend.text = element_text(hjust = 0),
          # legend.box = "vertical",
          plot.subtitle = element_text(face="bold")) +
    ylim(50, 100) +
    xlim(-15, 15) +
    labs(y = expression(A[net] ~ "(" ~ T[growth] ~ ") as percentage of " ~ A[opt] ~ "[%]"),
         x = expression(T[opt] ~ "-" ~ T[growth] ~ "[°C]"),
         subtitle = expression(bold("Trajectories of" ~ A[net] ~ "under warming")))
  
  out <- list(
    p = p_dt_tgrowth,
    df_traj = df_traj,
    df_dt = df_dt
  )
    
  return(out)
}

# Working Zone Below ----

# stop("Script sourced.")
# 
# settings <- get_settings()
# settings$daily_conditions <- "B"
# settings$tau <- 15
# out <- k19_from_forc_to_plot(settings)
# 
# df_plot <-
#   out$df %>%
#   unnest(rpm_inst) %>%
#   mutate(temp_rnd = round(temp)) %>%
#   dplyr::filter(temp_rnd == tc_leaf) %>%
#   mutate(delta_t_obs = abs(eva_data_tc_opt - temp),
#          delta_t_mod = abs(rpm_sim_tc_opt - temp),
#          delta_a     = abs(anet*1e6 - rpm_sim_anet_opt))
# 
# df_plot %>%
#   ggplot() +
#   aes(delta_t_mod, delta_a) +
#   geom_point() +
#   geom_smooth()
# 
# df_plot %>%
#   # mutate(delta_t_obs = abs(eva_data_tc_opt - temp)) %>%
#   # mutate(delta_t_obs = eva_data_tc_opt - temp) %>%
#   ggplot(aes(temp, delta_t_obs)) +
#   # geom_point() +
#   geom_smooth(se = FALSE) +
#   geom_hline(yintercept = 0, linetype = "dotted") +
#   theme_classic() +
#   ylim(0, 10) +
#   xlim(0, 35)
# 
# df_plot %>%
#   ggplot(aes(temp, delta_t_mod)) +
#   # geom_point() +
#   geom_smooth() +
#   geom_hline(yintercept = 0, linetype = "dotted") +
#   theme_classic() +
#   ylim(-10, 10)
# 
# df_tmp <-
# df_tmp <-
# 
# df_noacc <-
#   df_tmp$df %>%
#   unnest(rpm_inst) %>%
#   mutate(temp_rnd = round(temp)) %>%
#   dplyr::filter(temp_rnd == tc_leaf) %>%
#   mutate(delta_t_obs = abs(eva_data_tc_opt - temp),
#          delta_t_mod = abs(rpm_sim_tc_opt - temp),
#          delta_a     = abs(anet - rpm_sim_anet_opt))
# 
# df_noacc %>%
#   ggplot(aes(temp, delta_t_mod)) +
#   # geom_point() +
#   geom_smooth() +
#   geom_hline(yintercept = 0, linetype = "dotted") +
#   theme_classic() +
#   ylim(-10, 10)
# 
# df_noacc %>%
#   # mutate(delta_t_obs = abs(rpm_sim_tc_opt - temp)) %>%
#   mutate(delta_t_obs = rpm_sim_tc_opt - temp) %>%
#   ggplot(aes(temp, delta_t_obs)) +
#   geom_point() +
#   geom_smooth() +
#   geom_hline(yintercept = 0, linetype = "dotted") +
#   theme_classic() +
#   ylim(-10, 10) +
#   xlim(0, 30)
# 
# 
# df_plot %>%
#   dplyr::mutate(obs = eva_data_tc_opt - temp,
#                 fullacc = rpm_sim_tc_opt - temp) %>%
#   dplyr::select(sitename, date, obs, fullacc, temp) %>%
#   left_join(df_plot_noacc$df %>%
#               dplyr::mutate(noacc = rpm_sim_tc_opt - temp) %>%
#               dplyr::select(sitename, date, noacc, temp)) %>%
#   pivot_longer(cols = c(obs, fullacc, noacc)) %>%
#   mutate(value = abs(value)) %>%
#   ggplot() +
#   aes(temp, value, color = name) %>%
#   geom_smooth(se = F) +
#   theme_classic() +
#   geom_hline(yintercept = 0, linetype = "dotted")  +
#   ylim(-10, 10) +
#   xlim(0, 35)
# 
# ## ANALYSIS OF DELTA_T AND DELTA_A
# 
# turn_dfplot_into_dfdelta <- function(df_plot,
#                                      setup) {
# 
#   df_out <-
#     df_plot %>%
#     unnest(rpm_inst) %>%
#     mutate(temp_rnd = round(temp),
#            id       = paste(sitename, date)) %>%
#     dplyr::filter(temp_rnd == tc_leaf) %>%
#     mutate(delta_t     = temp - rpm_sim_tc_opt,
#            delta_t_abs = abs(delta_t),
#            anet        = anet*1e6,
#            delta_a     = anet - rpm_sim_anet_opt,
#            delta_a_abs = abs(delta_a),
#            rel_red_aopt= 1 - (rpm_sim_anet_opt - anet) / rpm_sim_anet_opt,
#            setup       = as.factor(setup))
# 
#   if (setup == "obs") {
#     df_out <-
#       df_out %>%
#       mutate(delta_t = abs(eva_data_tc_opt - temp))
#   }
# 
#   df_out <-
#     df_out %>%
#     dplyr::select(id,
#                   delta_t,
#                   delta_t_abs,
#                   delta_a,
#                   delta_a_abs,
#                   rel_red_aopt,
#                   setup,
#                   temp)
# 
#   return(df_out)
# }
# 
# out_fullacc <- readRDS("~/Insync/96p.schneider@gmail.com/repos/padasch/pivate_paper_therm_acclimation/output/model_runs/isolated_acclimation_processes/2023-02-10/run_01/df_final_fullac.rds")
# out_noacc   <- readRDS("~/Insync/96p.schneider@gmail.com/repos/padasch/pivate_paper_therm_acclimation/output/model_runs/isolated_acclimation_processes/2023-02-10/run_01/df_final_noacc_scaled.rds")
# out_er      <- readRDS("~/Insync/96p.schneider@gmail.com/repos/padasch/pivate_paper_therm_acclimation/output/model_runs/isolated_acclimation_processes/2023-02-10/run_01/df_final_er.rds")
# 
# df_obs <- turn_dfplot_into_dfdelta(out_fullacc$df, setup = "obs")
# df_fullacc <- turn_dfplot_into_dfdelta(out_fullacc$df, setup = "fullacc")
# df_noacc   <- turn_dfplot_into_dfdelta(out_noacc$df, setup = "noacc")
# df_er   <- turn_dfplot_into_dfdelta(out_er$df, setup = "er")
# 
# 
# ## Final Figure
# # Bind dataframes to one
# df_long <-
#   rbind(df_fullacc,
#         df_noacc,
#         df_er)
# 
# # Get coordinates for illustrative points
# 
# # From df_long
# id_minimum_dt <- df_long %>% arrange(delta_t) %>% slice(1) %>% pull(id)
# df_minimum_dt <- df_long %>% dplyr::filter(id == id_minimum_dt)
# 
# df_minimum_dt
# 
# # From hand
# # noacc
# df_minimum_dt$rel_red_aopt[2] <- 0.555 # original 0.574
# df_minimum_dt$delta_t[2]      <- -14.6  # original -14.9
# 
# # fullacc
# df_minimum_dt$rel_red_aopt[1] <- 0.905 # original
# df_minimum_dt$delta_t[1]      <- -6.65  # original -6.44
# 
# # er
# # no adjustment needed, original values fit well
# 
# # p_dt_tgrowth <-
#   df_long %>%
#   ggplot() +
#   aes(delta_t, rel_red_aopt*100, color = setup) +
#   scale_color_brewer("Model Setup: ",
#                      labels = c("Full Acclimation", "No Acclimation", "ER"),
#                      palette = "Set1") +
#   scale_fill_brewer("Model Setup: ",
#                     labels = c("Full Acclimation", "No Acclimation", "ER"),
#                     palette  = "Set1") +
#   # geom_point(alpha = 0.1) +
#   geom_smooth(se = FALSE) +
#   geom_point(
#     data = df_minimum_dt,
#     aes(delta_t, rel_red_aopt*100, fill = setup),
#     color = "black",
#     shape = 21,
#     size  = 3
#   ) +
#   annotate("text", label = "Arctic",  x = -12,  y = 50 #, size = 8, colour = "red"
#   ) +
#   annotate("text", label = "Boreal",  x = -1.5,  y = 50 #, size = 8, colour = "red"
#   ) +
#   annotate("text", label = "Temperate",  x = 1.75,  y = 50 #, size = 8, colour = "red"
#   ) +
#   annotate("text", label = "Tropic",  x = 10,  y = 50 #, size = 8, colour = "red"
#   ) +
#   geom_vline(xintercept = 0, linetype = "dotted") +
#   theme_classic() +
#   theme(legend.position = "bottom") +
#   ylim(50, 100) +
#   xlim(-15, 15) +
#   labs(y = expression(A[net] ~ "at " ~ T[growth] ~ "as percentage of " ~ A[opt] ~ "[%]"),
#        x = expression(T[opt] ~ "-" ~ T[growth] ~ "[°C]"))
# 
# # p_dt_tgrowth
# 
# 
# ## Figures playground
# # Delta T versus Tgrowth
# # Absolute delta T
# rbind(df_fullacc,
#       # df_obs,
#       df_noacc,
#       df_er) %>%
#   ggplot() +
#   aes(temp, delta_t_abs, color = setup) +
#   # geom_point() +
#   geom_smooth(se = FALSE) +
#   geom_hline(yintercept = 0, linetype = "dotted") +
#   theme_classic() +
#   ylim(0, 10) +
#   xlim(0, 35) +
#   labs(y = "|T_opt - T_growth| [degC]",
#        x = "T_growth [degC]")
# 
# # Normal delta T
# p_dt_tgrowth <-
#   rbind(df_fullacc,
#         # df_obs,
#         df_noacc,
#         df_er) %>%
#   ggplot() +
#   aes(temp, delta_t, color = setup) +
#   # geom_point() +
#   geom_smooth(se = FALSE) +
#   geom_hline(yintercept = 0, linetype = "dotted") +
#   theme_classic() +
#   ylim(-10, 10) +
#   xlim(0, 35) +
#   labs(y = "T_opt - T_growth [degC]",
#        x = "T_growth [degC]")
# 
# p_dt_tgrowth
# 
# # Delta A versus Delta T
# # Relative to 100% Aopt
# p_rel_red <-
#   rbind(df_fullacc,
#         df_noacc,
#         df_er) %>%
#   ggplot() +
#   aes(delta_t_abs, rel_red_aopt*100, color = setup) +
#   # geom_point() +
#   geom_smooth(se = FALSE) +
#   geom_hline(yintercept = 0, linetype = "dotted") +
#   theme_classic() +
#   ylim(50, 100) +
#   xlim(0, 15) +
#   labs(y = "Change in A_net at T_growth [% of A_opt]",
#        x = "|T_opt - T_growth| [degC]")
# 
# p_rel_red
# 
# # Absolute delta A
# rbind(df_fullacc,
#       df_noacc,
#       df_er) %>%
#   ggplot() +
#   aes(delta_t_abs, delta_a_abs, color = setup) +
#   # geom_point() +
#   geom_smooth(se = FALSE) +
#   geom_hline(yintercept = 0, linetype = "dotted") +
#   theme_classic() +
#   ylim(0, 3) +
#   xlim(0, 15) +
#   labs(y = "|A_opt - A_growth| [umol/m2/s]",
#        x = "|T_opt - T_growth| [degC]")
# 
# # Normal delta A
# rbind(df_fullacc,
#       df_noacc,
#       df_er) %>%
#   ggplot() +
#   aes(delta_t, delta_a, color = setup) +
#   # geom_point() +
#   geom_smooth(se = FALSE) +
#   geom_hline(yintercept = 0, linetype = "dotted") +
#   theme_classic() +
#   ylim(-3, 3) +
#   xlim(0, 15) +
#   labs(y = "A_opt - A_growth [umol/m2/s]",
#        x = "|T_opt - T_growth| [degC]")
# 
# ## Patchwork
# (p_dt_tgrowth + p_rel_red) +
#   plot_layout(guides = "collect") &
#   theme(legend.position = "bottom")
# 
# 
# 
