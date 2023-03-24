# Script for running smaller analyses ----
dir_tmp <- here("output/additional_analyses/")
if (!dir.exists(dir_tmp)) dir.create(dir_tmp, recursive = T, showWarnings = F)

# SINLGE RUNS ----

settings <- get_settings()
settings$daily_conditions <- "B"
settings$tau <- 15
# settings$method_ftemp <- "leuning02" 
# settings$rpmodel_exp <- 2000
settings$load_forcing <- F
out <- k19_from_forc_to_plot(settings)
plots <- plot_all_final_plots_from_df_plot(out$df, out$set)

# BEST MODEL FIGURES ----
settings <- get_settings()
settings$daily_conditions <- "B"
settings$tau <- 15
# settings$load_forcing <- F
settings$save_plots <- T
out <- k19_from_forc_to_plot(settings)

# ├ Observed vs. Modelled Patterns ----
p_mod <- 
  (out$p$p_topt_mod_tcair  + labs(subtitle = NULL, title = "Modelled", caption = NULL, x = NULL, y = bquote(T[opt] ~ " [°C]")) + theme(plot.title = element_text(hjust = 0)) + 
  out$p$p_aopt_mod_tcair  + labs(subtitle = NULL,  title = NULL, caption = NULL, x = NULL, y = bquote(A[opt] ~ "[µmol" ~ m ^-2 ~ s ^-1 ~ "]")) + 
  out$p$p_tspan_mod_tcair + labs(subtitle = NULL,  title = NULL, caption = NULL, y = bquote(T[span] ~ " [°C]"))) +
  plot_layout(guides = "collect", ncol = 1) &
  theme(legend.position = "bottom")

p_obs <- 
  (out$p$p_topt_obs_tcair  + labs(subtitle = NULL, title = "Observed", caption = NULL, x = NULL, y = NULL) + theme(plot.title = element_text(hjust = 0))) + 
  out$p$p_aopt_obs_tcair  + labs(subtitle = NULL, title = NULL, caption = NULL, x = NULL, y = NULL) + 
  out$p$p_tspan_obs_tcair + labs(subtitle = NULL, title = NULL, caption = NULL, y = NULL) +
  plot_layout(guides = "collect", ncol = 1) &
  theme(legend.position = "bottom")
  
(p <- p_obs | p_mod) +
  plot_layout(guides = "collect", ncol = 2) &
  theme(legend.position = "bottom")

# ggsave(paste0(dir_bm, "/modelled-observed-patterns.pdf"),
#        p,
#        height = 12,
#        width = 8)

# ├ Modobs ----
(p <-
   out$p$modobs_topt  + labs(caption = NULL, title = NULL) +
   out$p$modobs_aopt  + labs(caption = NULL, title = NULL) +
   out$p$modobs_tspan + labs(caption = NULL, title = NULL) + 
   plot_layout(guides = "collect") +
   plot_annotation(title = expression(underline(bold("Observed vs. Modelled Traits (Full Model Setup)")))) &
   theme(legend.position = "bottom",
         plot.title = element_text(hjust = 0, size = 14)) 
 )

dir_save <- paste0(dir_tmp, "/best_model_setup")
if (!dir.exists(dir_save)) dir.create(dir_save, recursive = T, showWarnings = F)

ggsave(paste0(dir_save, "/modobs-best-model.pdf"),
       p,
       height = 5,
       width = 10)

# ├ Curve shapes per climate zone ----
ggsave(paste0(dir_save, "/curve-shape-per-climate.pdf"),
       out$p$p_scaled_cl_2,
       height = 5,
       width = 10)

# ├ Curve shapes per climate zone ----
ggsave(paste0(dir_save, "/curve-shape-per-climate-parab.pdf"),
       out$p$p_scaled_cl_2_parab,
       height = 5,
       width = 10)

# ├ P-Model Acclimated Traits ---------------------------------------------------------------
ggsave(paste0(dir_save, "/acclimated-traits.pdf"),
       out$p$p_traits,
       height = 6,
       width  = 9)

# ├ Modelled curve for each site ---------------------------------------------------------------
ggsave(paste0(dir_save, "/response-of-each-site.pdf"),
       out$p$p_big_plot,
       height = 15,
       width  = 26)

#__________________________________________________________________________________________________
# LEUNING02 RUN ----
settings <- get_settings()
settings$daily_conditions <- "B"
settings$tau <- 15
settings$method_ftemp <- "leuning02"
out <- k19_from_forc_to_plot(settings)

dir_save <- paste0(dir_tmp, "replacing_k19_with_l02")
if (!dir.exists(dir_save)) dir.create(dir_save, recursive = T, showWarnings = F)

ggsave(paste0(dir_save, "/curve-shape-per-climate-l02.pdf"),
       out$p$p_scaled_cl_2 + 
         labs(subtitle = "Model Setup: PC + SS (using Leuning 2002 for ER)") + 
         theme(plot.subtitle = element_text(hjust = 0)),
       height = 4,
       width = 8)


#.----
# BEST MODEL SEARCH ------------------------------------------------------------

# ├ Effect of phi_0  -----------------------------------------------------------
# ├— Setup: Create vector of parameter to investigate ----
par_vec <- seq(0.025, 0.125, length.out=10) # 0.05 based on Smith2019, 0.09 based on calibrations in Stocker2020, theoretical 1.25 Mengoli2022
par_vec <- sort(c(par_vec, 0.04479))
par_vec <- round(par_vec, 5)

# Define model settings that are not associated with parameter
settings <- get_settings()
settings$vpd_inst_method <- "climate_scaled"
settings$dir_prefix <- "effect_of_kphio_on_modobs"
settings$daily_conditions <- "B"
settings$tau <- 15
settings$save_plots <- F

# Define empty list for plots and tibble for metrics
metrics_topt   <- tibble()
metrics_aopt   <- tibble()
metrics_agro   <- tibble()
metrics_tspan  <- tibble()
plots_topt     <- list()
plots_aopt     <- list()
plots_agro     <- list()
plots_tspan    <- list()

# ├— Run Model Loop ----
for (i in 1:length(par_vec)) {
  
  # Preparations
  message("\014 [>] Working on ", i, "/", length(par_vec), ": ", par_vec[i])
  
  # Overwrite settings with looped variables
  settings$kphio_calib <- par_vec[i]
  
  # Run Model
  out         <- k19_from_forc_to_plot(settings, returndf = "df_plot") 
  plots       <- plot_all_for_one_model(out$df)
  
  # T_opt
  plots_topt[[paste0("kphiocalib_", par_vec[i])]] <- plots$modobs_topt + labs(caption = NULL)
  plots$metrics_modobs_topt$kphio_calib  <- par_vec[i]
  metrics_topt <- rbind(metrics_topt, plots$metrics_modobs_topt)
  
  # A_opt
  plots_aopt[[paste0("kphiocalib_", par_vec[i])]] <- plots$modobs_aopt + labs(caption = NULL)
  plots$metrics_modobs_aopt$kphio_calib  <- par_vec[i]
  metrics_aopt <- rbind(metrics_aopt, plots$metrics_modobs_aopt)
  
  # # A_growth
  # plots_agro[[paste0("kphiocalib_", par_vec[i])]] <- plots$modobs_agrowth + labs(caption = NULL)
  # plots$metrics_modobs_agrowth$kphio_calib  <- par_vec[i]
  # metrics_agro <- rbind(metrics_agro, plots$metrics_modobs_agrowth)
  
  # T_span
  plots_tspan[[paste0("kphiocalib_", par_vec[i])]] <- plots$modobs_tspan + labs(caption = NULL)
  plots$metrics_modobs_tspan$kphio_calib  <- par_vec[i]
  metrics_tspan <- rbind(metrics_tspan, plots$metrics_modobs_tspan)
}

# ├— Attaching model performance label ----
metrics_topt$model_perf <- NA
metrics_aopt$model_perf <- NA
metrics_tspan$model_perf <- NA

for (i in 1:nrow(metrics_topt)) {
  ## Topt
  prf <- "Good"
  if ((metrics_topt$r2[i] < 0.5) | (metrics_topt$b0_sign_not_0[i] | metrics_topt$b1_sign_not_one[i])) prf <- "Poor"
  metrics_topt$model_perf[i] <- prf
  
  ## Aopt
  prf <- "Good"
  if ((metrics_aopt$r2[i] < 0.5) | (metrics_aopt$b0_sign_not_0[i] | metrics_aopt$b1_sign_not_one[i])) prf <- "Poor"
  metrics_aopt$model_perf[i] <- prf
  
  ## Tspan
  prf <- "Good"
  if ((metrics_tspan$r2[i] < 0.5) | (metrics_tspan$b0_sign_not_0[i] | metrics_tspan$b1_sign_not_one[i])) prf <- "Poor"
  metrics_tspan$model_perf[i] <- prf
}

## Attaching NA for missing model performance factor for easier ggplotting
if (length(unique(metrics_topt$model_perf)) == 1) {
  newrow <- metrics_topt[1, ]
  newrow[1, ] <- NA
  newrow$model_perf <- str_remove("GoodPoor", unique(metrics_topt$model_perf)) # Keep missing model performance factor
  metrics_topt <- rbind(metrics_topt, newrow)
}

if (length(unique(metrics_aopt$model_perf)) == 1) {
  newrow <- metrics_aopt[1, ]
  newrow[1, ] <- NA
  newrow$model_perf <- str_remove("GoodPoor", unique(metrics_aopt$model_perf)) # Keep missing model performance factor
  metrics_aopt <- rbind(metrics_aopt, newrow)
}

if (length(unique(metrics_tspan$model_perf)) == 1) {
  newrow <- metrics_tspan[1, ]
  newrow[1, ] <- NA
  newrow$model_perf <- str_remove("GoodPoor", unique(metrics_tspan$model_perf)) # Keep missing model performance factor
  metrics_tspan <- rbind(metrics_tspan, newrow)
}


# ├— Creating metrics plot ----
# Make labeller function
vnames <-list("r2" = bquote(R^2 ~ "[-]"), "rmse"  = bquote("RMSE [-]"), "bias" = bquote("Bias [-]"))
vlabeller <- function(variable,value){return(vnames[value])}

# Create plot metric vs. parameter
p_metrics_topt <-
  metrics_topt %>%
  pivot_longer(cols = c(r2, rmse, bias), names_to = "metric", values_to = "value") %>%
  ggplot() +
  aes(x = kphio_calib, y = value) +
  geom_line() +
  geom_point(aes(shape = model_perf), size = 1.5, alpha = 0.8) +
  scale_shape_manual(breaks = c("Good", "Poor"),
                     values = c(19, 1)) +
  scale_color_brewer(palette = "Dark2") +
  labs(color = "Condition:", shape = "Model Fit: ") +
  theme_classic() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  facet_wrap(~metric, ncol = 1, scales = "free_y",  labeller = vlabeller) +
  ylab("Metric Value") +
  xlab(bquote(varphi[0] ~ "[-]")) +
  ggtitle(bquote(T[opt]))

p_metrics_aopt <-
  metrics_aopt %>% 
  pivot_longer(cols = c(r2, rmse, bias), names_to = "metric", values_to = "value") %>%
  ggplot() +
  aes(x = kphio_calib, y = value) +
  geom_line() +
  geom_point(aes(shape = model_perf), size = 1.5, alpha = 0.8) +
  scale_shape_manual(breaks = c("Good", "Poor"),
                     values=c(19, 1)) +
  scale_color_brewer(palette = "Dark2") +
  labs(color = "Condition:", shape = "Model Fit: ") +
  theme_classic() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  facet_wrap(~metric, ncol = 1, scales = "free_y",  labeller = vlabeller) +
  ylab(NULL) +
  xlab(bquote(varphi[0] ~ "[-]")) +
  ggtitle(bquote(A[opt]))

p_metrics_tspan <-
  metrics_tspan %>% 
  pivot_longer(cols = c(r2, rmse, bias), names_to = "metric", values_to = "value") %>%
  ggplot() +
  aes(x = kphio_calib, y = value) +
  geom_line() +
  geom_point(aes(shape = model_perf), size = 1.5, alpha = 0.8) +
  scale_shape_manual(breaks = c("Good", "Poor"),
                     values=c(19, 1)) +
  scale_color_brewer(palette = "Dark2") +
  labs(color = "Condition:", shape = "Model Fit: ") +
  theme_classic() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  facet_wrap(~metric, ncol = 1, scales = "free_y",  labeller = vlabeller) +
  ylab(NULL) +
  xlab(bquote(varphi[0] ~ "[-]")) +
  ggtitle(bquote(T[span]))

## ├— Make and save final plot ----
p_metrics <- p_metrics_topt + p_metrics_aopt + p_metrics_tspan + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.direction = "horizontal") 

dir_save <- paste0(dir_tmp, "calibration/kphio/B15_narrower_range/")
if (!dir.exists(dir_save)) dir.create(dir_save, recursive = T, showWarnings = F)

ggsave(paste0(dir_save, "000_metrics_comparison.pdf"), p_metrics, height = 6, width = 10)

## ├— Create overview plot of each modobs plot: ----
# Create first plot to which others are added
kphio_i        <- par_vec[[1]]
p_modobs_topt  <- plots_topt[[1]] + labs(x = NULL, title = NULL, subtitle = bquote(varphi[0] ~ "=" ~ .(kphio_i)))
p_modobs_aopt  <- plots_aopt[[1]] + labs(x = NULL, title = NULL, subtitle = bquote(varphi[0] ~ "=" ~ .(kphio_i)))
# p_modobs_agro  <- plots_agro[[1]] + labs(x = NULL, title = NULL, subtitle = bquote(varphi[0] ~ "=" ~ .(kphio_i)))
p_modobs_tspan <- plots_tspan[[1]] + labs(x = NULL, title = NULL, subtitle = bquote(varphi[0] ~ "=" ~ .(kphio_i)))

for (i in c(1, 3, 6, 11)) {
  # Skip first plot because provided in outer loop
  if(i == 1) next
  
  kphio_i <- par_vec[i]
  
  p_tmp_topt <- plots_topt[[paste0("kphiocalib_", kphio_i)]] + labs(title = NULL, subtitle = bquote(varphi[0] ~ "=" ~ .(kphio_i)))
  p_tmp_aopt <- plots_aopt[[paste0("kphiocalib_", kphio_i)]] + labs(title = NULL, subtitle = bquote(varphi[0] ~ "=" ~ .(kphio_i)))
  # p_tmp_agro <- plots_agro[[paste0("kphiocalib_", kphio_i)]] + labs(title = NULL, subtitle = paste0("kphio_calib = ", .(kphio_i)))
  p_tmp_tspan <- plots_tspan[[paste0("kphiocalib_", kphio_i)]] + labs(title = NULL, subtitle = bquote(varphi[0] ~ "=" ~ .(kphio_i)))

  ## Top right plot does not need xy axis
  if (i == 3) {
    p_tmp_topt <- p_tmp_topt + labs(x = NULL, y = NULL)
    p_tmp_aopt <- p_tmp_aopt + labs(x = NULL, y = NULL)
    # p_tmp_agro <- p_tmp_agro + labs(x = NULL, y = NULL)
    p_tmp_tspan <- p_tmp_tspan + labs(x = NULL, y = NULL)
  }  

  ## Bottom right plot does not need y axis
  if (i == 11) {
    p_tmp_topt <- p_tmp_topt + labs(y = NULL)
    p_tmp_aopt <- p_tmp_aopt + labs(y = NULL) + xlim(0, max(p_tmp_aopt$data$rpm_sim_anet_opt))
    # p_tmp_agro <- p_tmp_agro + labs(y = NULL)
    p_tmp_tspan <- p_tmp_tspan + labs(y = NULL)
  }  
  
  p_modobs_topt <- p_modobs_topt + p_tmp_topt
  p_modobs_aopt <- p_modobs_aopt + p_tmp_aopt
  # p_modobs_agro <- p_modobs_agro + p_tmp_agro
  p_modobs_tspan <- p_modobs_tspan + p_tmp_tspan
}

## ├— Save plot for given condition ----
p_final_topt <- 
  p_modobs_topt +
  plot_layout(guides = "collect", nrow = 2) & 
  theme(legend.position = "bottom") &
  plot_annotation(subtitle = bquote("Effect of" ~ varphi[0] ~ "on" ~ T[opt]))

p_final_aopt <- 
  p_modobs_aopt +
  plot_layout(guides = "collect", nrow = 2) & 
  theme(legend.position = "bottom") &
  plot_annotation(subtitle = bquote("Effect of" ~ varphi[0] ~ "on" ~ A[opt]))

# p_final_agro <- 
#   p_modobs_agro +
#   plot_layout(guides = "collect", nrow = 1) & 
#   theme(legend.position = "bottom") &
#   plot_annotation(title = paste0("Effect of A_growth on kphio_calib: "))

p_final_tspan <- 
  p_modobs_tspan +
  plot_layout(guides = "collect", nrow = 2) & 
  theme(legend.position = "bottom") &
  plot_annotation(subtitle = bquote("Effect of" ~ varphi[0] ~ "on" ~ T[span]))

sub_dir <- here("output/additional_analyses/calibration/kphio")

ggsave(paste0(sub_dir, "/000_topt.pdf"), p_final_topt, height = 4, width = 4)
ggsave(paste0(sub_dir, "/000_aopt.pdf"), p_final_aopt, height = 4, width = 4)
# ggsave(paste0(sub_dir, "/000_agro.pdf"), p_final_agro, height = 4, width = 4)
ggsave(paste0(sub_dir, "/000_tspan.pdf"), p_final_tspan, height = 4, width = 4)

# ├— Final SI Fig. ----
p <- (p_metrics / p_final_aopt) & 
  plot_annotation(tag_levels = "A", tag_suffix = ".)") &
  theme(text = element_text(size = 10),
        plot.tag = element_text(face = "bold"))

ggsave(here("output/additional_analyses/calibration/kphio/si-effect-kphio.pdf"), p, height = 12, width = 8)

# ├ Effect of acclimation assumption and timescale -----------------------------
# This section is to test the effect of different acclimation timescales
# on the prediction of t_opt. It runs different timescales and compares their 
# effect on model metrics and saves plots.
# The longer a timescale, the more "averaged" is the input 

# Debug:
vec_cond <- c("A", "B")
vec_time <- c(15, 30, 90)

# ├— Setup: Create vector of timescales to investigate ----
vec_cond <- c("A", "B", "D")
vec_time <- c(1, 3, 5, 7, 10, 15, 30, 45, 60, 90, 180)

# Define model settings that are not associated with timescale and condition
settings <- get_settings()
settings$vpd_inst_method <- "climate_scaled"
settings$save_plots <- F
# settings$dir_prefix <- "effect-of-acclimation-assumptions"

# Define empty list for plots and tibble for metrics
metrics_topt   <- tibble()
metrics_aopt   <- tibble()
metrics_tspan  <- tibble()
plots_topt     <- list()
plots_aopt     <- list()
plots_tspan    <- list()

# ______________________________________________________________________________
# ├— Run Model Loop ----
for (cond in vec_cond) {

    for (time in vec_time) {
      
      ## Short-cut to add conditions if others were run before already (just run loop again after break)
      if ("cond" %in% names(metrics_topt)) {
        if (paste0(cond, time) %in% paste0(metrics_topt$cond, metrics_topt$tau)) {
          next
        }
      }
      
      message("\014 [>] Working on: ", cond, " ", time)

      # Overwrite settings with looped variables
      settings$daily_conditions <- cond
      settings$tau <- time
      
      # Run Model
      out         <- k19_from_forc_to_plot(settings)
      
      # T_opt
      plots_topt[[cond]][[paste0("tau_", time)]] <- out$p$modobs_topt
      out$p$metrics_modobs_topt$tau  <- time
      out$p$metrics_modobs_topt$cond <- cond
      metrics_topt <- rbind(metrics_topt, out$p$metrics_modobs_topt)
      
      # A_opt
      plots_aopt[[cond]][[paste0("tau_", time)]] <- out$p$modobs_aopt
      out$p$metrics_modobs_aopt$tau  <- time
      out$p$metrics_modobs_aopt$cond <- cond
      metrics_aopt <- rbind(metrics_aopt, out$p$metrics_modobs_aopt)
      
      # T_opt
      plots_tspan[[cond]][[paste0("tau_", time)]] <- out$p$modobs_tspan
      out$p$metrics_modobs_tspan$tau  <- time
      out$p$metrics_modobs_tspan$cond <- cond
      metrics_tspan <- rbind(metrics_tspan, out$p$metrics_modobs_tspan)
  }
}

# ______________________________________________________________________________
# ├— Attach model performance ----
# Note: Definition of good and poor model performance
# b0_sign_not_0:   If true, then poor performance
# b1_sign_not_one: If true, then poor performance
# If b0_sign_not_0 OR b1_sign_not_one TRUE, then poor performance

metrics_topt$model_perf <- NA
metrics_aopt$model_perf <- NA
metrics_tspan$model_perf <- NA

for (i in 1:nrow(metrics_topt)) {
  ## Topt
  prf <- "Good"
  if ((metrics_topt$pval[i] > 0.05) | metrics_topt$b0_sign_not_0[i] | metrics_topt$b1_sign_not_one[i]) prf <- "Poor"
  metrics_topt$model_perf[i] <- prf
  
  ## Aopt
  prf <- "Good"
  if ((metrics_aopt$pval[i] > 0.05) | metrics_aopt$b0_sign_not_0[i] | metrics_aopt$b1_sign_not_one[i]) prf <- "Poor"
  metrics_aopt$model_perf[i] <- prf
  
  ## Tspan
  prf <- "Good"
  if ((metrics_tspan$pval[i] > 0.05) | metrics_tspan$b0_sign_not_0[i] | metrics_tspan$b1_sign_not_one[i]) prf <- "Poor"
  metrics_tspan$model_perf[i] <- prf
}

# ├— Plotting ----
# Make labeller function
vnames <-list("r2" = bquote(R^2 ~ "[-]"),
              "rmse"  = bquote("RMSE [-]"),
              "bias" = bquote("Bias [-]"))

vlabeller <- function(variable,value){
  return(vnames[value])
}

# Day-Conditions to filter for
filt_cond <- c("A_night", "A", "B", "C","D")


## ├—— Metrics vs. Tau ----
p_metrics_topt <-
  metrics_topt %>%
  filter(cond %in% filt_cond) %>%
  mutate(
    cond = ifelse(cond == "A_night", paste0("Daily Mean (Light without night values)"), cond),
    cond = ifelse(cond == "A", paste0("Daily Mean"), cond),
    cond = ifelse(cond == "B", paste0("Midday"), cond),
    cond = ifelse(cond == "C", paste0("Midday Light, Daily Mean Rest"), cond),
    cond = ifelse(cond == "D", paste0("Daily Max."), cond),
    cond = as.factor(cond)
  ) %>%
  pivot_longer(cols = c(r2, rmse, bias), names_to = "metric", values_to = "value") %>%
  ggplot() +
  aes(x = tau, y = value, color = cond) +
  geom_line() +
  geom_point(aes(shape = model_perf), size = 1.5, alpha = 0.8) +
  scale_shape_manual(values=c(19, 1)) +
  scale_color_brewer(palette = "Dark2") +
  labs(color = "Condition:", shape = "Model Fit: ") +
  theme_classic() +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  facet_wrap(~metric, ncol = 1, scales = "free_y",  labeller = vlabeller) +
  ylab("Metric Value") +
  xlab(bquote(tau ~ " [days]")) +
  ggtitle(bquote(T[opt]))

p_metrics_aopt <-
  metrics_aopt %>% 
  filter(cond %in% filt_cond) %>%
  mutate(
    cond = ifelse(cond == "A_night", paste0("Daily Mean (Light without night values)"), cond),
    cond = ifelse(cond == "A", paste0("Daily Mean"), cond),
    cond = ifelse(cond == "B", paste0("Midday"), cond),
    cond = ifelse(cond == "C", paste0("Midday Light, Daily Mean Rest"), cond),
    cond = ifelse(cond == "D", paste0("Daily Max."), cond),
    cond = as.factor(cond)
  ) %>%
  pivot_longer(cols = c(r2, rmse, bias), names_to = "metric", values_to = "value") %>%
  ggplot() +
  aes(x = tau, y = value, color = cond) +
  geom_line() +
  geom_point(aes(shape = model_perf), size = 1.5, alpha = 0.8) +
  scale_shape_manual(breaks = c("Good", "Poor"),
                     values=c(19, 1)) +
  scale_color_brewer(palette = "Dark2") +
  labs(color = "Condition:", shape = "Model Fit: ") +
  theme_classic() +
  theme(legend.position = "none", legend.direction = "vertical") +
  facet_wrap(~metric, ncol = 1, scales = "free_y",  labeller = vlabeller) +
  ylab(NULL) +
  xlab(bquote(tau ~ " [days]")) +
  ggtitle(bquote(A[opt]))

p_metrics_tspan <-
  metrics_tspan %>% 
  filter(cond %in% filt_cond) %>%
  mutate(
    cond = ifelse(cond == "A_night", paste0("Daily Mean (Light without night values)"), cond),
    cond = ifelse(cond == "A", paste0("Daily Mean"), cond),
    cond = ifelse(cond == "B", paste0("Midday"), cond),
    cond = ifelse(cond == "C", paste0("Midday Light, Daily Mean Rest"), cond),
    cond = ifelse(cond == "D", paste0("Daily Max."), cond),
    cond = as.factor(cond)
  ) %>%
  pivot_longer(cols = c(r2, rmse, bias), names_to = "metric", values_to = "value") %>%
  ggplot() +
  aes(x = tau, y = value, color = cond) +
  geom_line() +
  geom_point(aes(shape = model_perf), size = 1.5, alpha = 0.8) +
  scale_shape_manual(breaks = c("Good", "Poor"),
                     values=c(19, 1)) +
  scale_color_brewer(palette = "Dark2") +
  labs(color = "Condition:", shape = "Model Fit: ") +
  theme_classic() +
  theme(legend.position = "none", legend.direction = "vertical") +
  facet_wrap(~metric, ncol = 1, scales = "free_y",  labeller = vlabeller) +
  ylab(NULL) +
  xlab(bquote(tau ~ " [days]")) +
  ggtitle(bquote(T[span]))

p_metrics <- p_metrics_topt + p_metrics_aopt + p_metrics_tspan + plot_layout(guides = "collect")

p_full <- 
  p_metrics / (plots_topt$B$tau_1    + labs(title = expression(paste(tau, " = 1")), caption = NULL) +
               plots_topt$B$tau_15 + labs(title = expression(paste(tau, " = 15")), caption = NULL) +
               plots_topt$B$tau_60 + labs(title = expression(paste(tau, " = 60")), caption = NULL)) +
  plot_layout(guides = "collect") &
  plot_annotation(tag_levels = "A", tag_suffix = ".)") &
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        plot.tag = element_text(face = "bold"))

dir_save <- paste0(dir_tmp, "calibration/acclimation-assumption/")
if (!dir.exists(dir_save)) dir.create(dir_save, recursive = T, showWarnings = F)
ggsave(paste0(dir_save, "/acc-assumption.pdf"), p_full, height = 10, width = 10)

## ├—— Metrics vs. Tau with example plots below ----

plots_topt$B$tau_1  + plots_topt$B$tau_15 + plots_topt$B$tau_60
plots_aopt$B$tau_1  + plots_aopt$B$tau_15 + plots_aopt$B$tau_60
plots_tspan$B$tau_1 + plots_tspan$B$tau_15 + plots_tspan$B$tau_60

## ├— Create overview plot of each modobs plot ----
for (cond in vec_cond) {
  # Create first plot to which others are added
  time <- vec_time[[1]]
  p_modobs_topt <- plots_topt[[paste0(cond)]][[1]] + ggtitle(paste0("tau = ", time))
  p_modobs_aopt <- plots_aopt[[paste0(cond)]][[1]] + ggtitle(paste0("tau = ", time))
  p_modobs_tspan <- plots_tspan[[paste0(cond)]][[1]] + ggtitle(paste0("tau = ", time))
  
  for (time in vec_time) {
    # Skip first plot because provided in outer loop
    if(time == vec_time[[1]]) next
    # Add all other plots
    p_modobs_topt <- p_modobs_topt + plots_topt[[cond]][[paste0("tau_", time)]] + ggtitle(paste0("tau = ", time))
    p_modobs_aopt <- p_modobs_aopt + plots_aopt[[cond]][[paste0("tau_", time)]] + ggtitle(paste0("tau = ", time))
    p_modobs_tspan <- p_modobs_tspan + plots_tspan[[cond]][[paste0("tau_", time)]] + ggtitle(paste0("tau = ", time))
  }
  
  ## Save plot for given condition
  p_final_topt <- 
    p_modobs_topt +
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom") &
    plot_annotation(title = paste0("Condition: ", cond))
  
  p_final_aopt <- 
    p_modobs_aopt +
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom") &
    plot_annotation(title = paste0("Condition: ", cond))
  
  p_final_tspan <- 
    p_modobs_tspan +
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom") &
    plot_annotation(title = paste0("Condition: ", cond))
  
  
  dir.create(here(dir_save, cond), recursive = T, showWarnings = F)
  
  ggsave(here(dir_save, cond, "/topt.pdf"), p_final_topt, height = 14, width = 14)
  ggsave(here(dir_save, cond, "/aopt.pdf"), p_final_aopt, height = 14, width = 14)
  ggsave(here(dir_save, cond, "/tspan.pdf"), p_final_tspan, height = 14, width = 14)
}


# . ----
# OTHER TESTS ----
# ├ High Light Test ----
dir_save <- here("output/additional_analyses/high-light")
if (!dir.exists(dir_save)) dir.create(dir_save, recursive = T, showWarnings = F)
settings <- get_settings()
settings$daily_conditions <- "B"
settings$tau <- 15
settings$save_plots <- F
settings$load_forcing <- F
settings$rpmodel_exp <- 2000
out <- k19_from_forc_to_plot(settings)

p <-
  out$p$modobs_topt  + labs(caption = NULL, title = bquote(T[opt])) +
  out$p$modobs_aopt  + labs(caption = NULL, title = bquote(A[opt])) +
  out$p$modobs_tspan + labs(caption = NULL, title = bquote(T[span])) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(paste0(dir_save, "/modobs-all.pdf"),
       p,
       height = 5,
       width  = 12)

# ├ VPD Approach ----

dir_tmp <- here("output/additional_analyses/vpd-approach/")
if (!dir.exists(dir_tmp)) dir.create(dir_tmp, recursive = T, showWarnings = F)

## Best Model (Climate Scaled)
settings     <- get_settings()
settings$save_plots <- F
settings$daily_conditions <- "B"
settings$tau <- 15
out_bm     <- k19_from_forc_to_plot(settings)

## VPD Const
settings$vpd_inst_method <- "climate_const"
out_climate_const     <- k19_from_forc_to_plot(settings)

## Meas. Scaled
settings$vpd_inst_method <- "meas_interp"
out_meas_interp     <- k19_from_forc_to_plot(settings)

## Reduce dfs to sites where VpdL was available for comparison
out_climate_const$df <- 
  out_climate_const$df %>% 
  mutate(qc = map_dbl(data_org, ~pull(., VpdL) %>% mean(.))) %>% 
  drop_na(qc)

out_bm$df <- 
  out_bm$df %>% 
  mutate(qc = map_dbl(data_org, ~pull(., VpdL) %>% mean(.))) %>% 
  drop_na(qc)

out_meas_interp$df <- 
  out_meas_interp$df %>% 
  mutate(qc = map_dbl(data_org, ~pull(., VpdL) %>% mean(.))) %>% 
  drop_na(qc)

out_bm$p            <- plot_all_final_plots_from_df_plot(out_bm$df, out_bm$set)
out_climate_const$p <- plot_all_final_plots_from_df_plot( out_climate_const$df, out_climate_const$set)
out_meas_interp$p   <- plot_all_final_plots_from_df_plot(out_meas_interp$df, out_meas_interp$set)

## Curve Plots
ggsave(paste0(dir_tmp, "/curve-shape-per-climate-vpdconst.pdf"),
       out_climate_const$p$p_scaled_cl_2,
       height = 4,
       width = 8)

## Modobs Plots
p_topt <-
  out_bm$p$modobs_topt             + labs(caption = NULL, title = NULL, subtitle = expression(paste("Inst. ", italic("VPD"), " scaled, from climate"))) +
  out_climate_const$p$modobs_topt  + labs(y = NULL, caption = NULL, title = NULL, subtitle = expression(paste("Inst. ", italic("VPD"), " constant, from climate"))) +
  out_meas_interp$p$modobs_topt    + labs(y = NULL, caption = NULL, title = NULL, subtitle = expression(paste("Inst. ", italic("VPD"), " scaled, from measurements")))+ 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(paste0(dir_tmp, "/topt_comparison.pdf"),
       p_topt,
       height = 5,
       width  = 12)

p_aopt <-
  out_bm$p$modobs_aopt             + labs(caption = NULL, title = NULL, subtitle = expression(paste("Inst. ", italic("VPD"), " scaled, from climate"))) +
  out_climate_const$p$modobs_aopt  + labs(y = NULL, caption = NULL, title = NULL, subtitle = expression(paste("Inst. ", italic("VPD"), " constant, from climate"))) +
  out_meas_interp$p$modobs_aopt    + labs(y = NULL, caption = NULL, title = NULL, subtitle = expression(paste("Inst. ", italic("VPD"), " scaled, from measurements")))+ 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(paste0(dir_tmp, "/aopt_comparison.pdf"),
       p_aopt,
       height = 5,
       width  = 12)

p_tspan <-
  out_bm$p$modobs_tspan             + labs(caption = NULL, title = NULL, subtitle = expression(paste("Inst. ", italic("VPD"), " scaled, from climate"))) +
  out_climate_const$p$modobs_tspan  + labs(y = NULL, caption = NULL, title = NULL, subtitle = expression(paste("Inst. ", italic("VPD"), " constant, from climate"))) +
  out_meas_interp$p$modobs_tspan    + labs(y = NULL, caption = NULL, title = NULL, subtitle = expression(paste("Inst. ", italic("VPD"), " scaled, from measurements")))+ 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(paste0(dir_tmp, "/tspan_comparison.pdf"),
       p_tspan,
       height = 5,
       width  = 12)

# p_all <- (p_topt / p_aopt / p_tspan) + plot_layout(guides = "collect")

p_all <-
  out_bm$p$modobs_topt               + labs(caption = NULL, title = NULL, subtitle = expression(bold(paste("Inst. ", italic("VPD"), " scaled, from climate")))) +
  out_climate_const$p$modobs_topt  + labs(y = NULL, caption = NULL, title = NULL, subtitle = expression(bold(paste("Inst. ", italic("VPD"), " constant, from climate")))) +
  out_meas_interp$p$modobs_topt    + labs(y = NULL, caption = NULL, title = NULL, subtitle = expression(bold(paste("Inst. ", italic("VPD"), " scaled, from measurements")))) + 
  
  out_bm$p$modobs_aopt             + labs(caption = NULL, title = NULL, subtitle = NULL) +
  out_climate_const$p$modobs_aopt  + labs(y = NULL, caption = NULL, title = NULL, subtitle = NULL) +
  out_meas_interp$p$modobs_aopt    + labs(y = NULL, caption = NULL, title = NULL, subtitle = NULL)+ 
  
  out_bm$p$modobs_tspan             + labs(caption = NULL, title = NULL, subtitle = NULL) +
  out_climate_const$p$modobs_tspan  + labs(y = NULL, caption = NULL, title = NULL, subtitle = NULL) +
  out_meas_interp$p$modobs_tspan    + labs(y = NULL, caption = NULL, title = NULL, subtitle = NULL)+ 
  plot_layout(guides = "collect", ncol = 3) &
  theme(legend.position = "bottom",
        plot.subtitle = element_text(hjust = 0))

ggsave(paste0(dir_tmp, "/vpd_approach_comparison.pdf"),
       p_all,
       height = 12,
       width  = 12)

# ├ Growth Conditions as Instant Forcing ----
dir_tmp <- here("output/additional_analyses/growth_forcing_as_instant/")
if (!dir.exists(dir_tmp)) dir.create(dir_tmp, recursive = T, showWarnings = F)

settings <- get_settings()
settings$daily_conditions <- "B"
settings$tau <- 15
settings$save_plots <- F

df_site_date <- k19_create_input_for_acc(settings)
df_acc       <- k19_run_acc(df_site_date, settings) 

for (i in 1:nrow(df_acc)) {
  df_acc$forcing_d[[i]] <- df_acc$forcing_growth[[i]] 
}

df_inst      <- k19_run_inst(df_acc, settings)
df_final     <- k19_create_final_df(df_inst, settings)
p_all        <- plot_all_for_one_model(df_final)

p_modobs <-
  p_all$modobs_topt  + labs(caption = NULL, title = bquote(T[opt])) +
  p_all$modobs_aopt  + labs(caption = NULL, title = bquote(A[opt])) +
  p_all$modobs_tspan + labs(caption = NULL, title = bquote(T[span])) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

p_mod <- 
  (out$p$p_topt_mod_tcair  + labs(subtitle = NULL, title = "Modelled", caption = NULL, x = NULL, y = bquote(T[opt] ~ " [°C]")) + theme(plot.title = element_text(hjust = 0)) + 
     out$p$p_aopt_mod_tcair  + labs(subtitle = NULL,  title = NULL, caption = NULL, x = NULL, y = bquote(A[opt] ~ "[µmol" ~ m ^-2 ~ s ^-1 ~ "]")) + 
     out$p$p_tspan_mod_tcair + labs(subtitle = NULL,  title = NULL, caption = NULL, y = bquote(T[span] ~ " [°C]"))) +
  plot_layout(guides = "collect", ncol = 1) &
  theme(legend.position = "bottom")

p_obs <- 
  (out$p$p_topt_obs_tcair  + labs(subtitle = NULL, title = "Observed", caption = NULL, x = NULL, y = NULL) + theme(plot.title = element_text(hjust = 0))) + 
  out$p$p_aopt_obs_tcair  + labs(subtitle = NULL, title = NULL, caption = NULL, x = NULL, y = NULL) + 
  out$p$p_tspan_obs_tcair + labs(subtitle = NULL, title = NULL, caption = NULL, y = NULL) +
  plot_layout(guides = "collect", ncol = 1) &
  theme(legend.position = "bottom")

p_patterns <- 
  (p_obs | p_mod) +
  plot_layout(guides = "collect", ncol = 2) &
  theme(legend.position = "bottom")

ggsave(paste0(dir_tmp, "/modobs-all.pdf"),
       p_modobs,
       height = 5,
       width  = 12)

ggsave(paste0(dir_tmp, "/modelled-observed-patterns.pdf"),
       p_patterns,
       height = 12,
       width = 8)

ggsave(paste0(dir_tmp, "/curve-shape-per-climate.pdf"),
       p_all$p_scaled_cl_2,
       height = 4,
       width = 8)
