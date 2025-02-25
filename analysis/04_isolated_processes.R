# ISOLATED ACCLIMATION PROCESSES ----
source("R/source.R")

# ├ Settings ----
dir_analysis <- "isolated_processes"
dir_tmp      <- "final"
k19_or_l02   <- "kumarathunge19_fixed"
use_pfts     <- T
vcmax25_source <- "kumar19"

# ├ Define directories ----
dir_mods <- make_new_dynamic_dir(dir_analysis, "model_rds", dir_tmp)
dir_figs <- make_new_dynamic_dir(dir_analysis, "figures", dir_tmp)
dir_tabs <- make_new_dynamic_dir(dir_analysis, "tables", dir_tmp)

dir_mods <- make_new_dynamic_dir(dir_analysis, "model_rds", dir_tmp, return_latest = TRUE)
dir_figs <- make_new_dynamic_dir(dir_analysis, "figures", dir_tmp, return_latest = TRUE)
dir_tabs <- make_new_dynamic_dir(dir_analysis, "tables", dir_tmp, return_latest = TRUE)

## Setups ----
### Full Accl. ----
settings <- get_settings()
settings$save_plots <- F
df_forc  <- k19_create_df_forcing(settings)
df_sd    <- k19_create_input_for_acc(settings, df_forc)
df_acc   <- k19_run_acc(df_sd, settings)
df_inst_fullacc  <- k19_run_inst(df_acc, settings)
df_final_fullacc <- k19_create_final_df(df_inst_fullacc, settings)
saveRDS(list(df = df_final_fullacc, set = settings), here(dir_mods, "df_final_fullacc.rds"))

### No Acclimation ----
settings <- get_settings()
settings$save_plots <- F
settings$method_ftemp <- k19_or_l02
df_acc   <- k19_run_acc(df_sd, settings)

if (use_pfts) df_acc <- replace_acclimated_with_pfts(df_acc, 
                                                     replace_pc = T, 
                                                     replace_xi = T,
                                                     jvr_method = settings$method_ftemp,
                                                     vcmax25_source = vcmax25_source)

df_inst_noacc  <- k19_run_inst(df_acc, settings)
df_final_noacc <- k19_create_final_df(df_inst_noacc, settings)
saveRDS(list(df = df_final_noacc, set = settings), here(dir_mods, "df_final_noacc_scaled.rds"))

### PC ----
settings <- get_settings()
settings$save_plots <- F
settings$method_ftemp <- k19_or_l02
settings <- k19_from_forc_to_plot(settings, returndf = "update_settings")

df_acc   <- k19_run_acc(df_sd, settings)

if (use_pfts) df_acc <- replace_acclimated_with_pfts(df_acc, 
                                                     replace_pc = F, 
                                                     replace_xi = T,
                                                     jvr_method = settings$method_ftemp,
                                                     vcmax25_source = vcmax25_source)

df_inst_pc    <- k19_run_inst(df_acc, settings)
df_final_pc <- k19_create_final_df(df_inst_pc, settings)
saveRDS(list(df = df_final_pc, set = settings), here(dir_mods, "df_final_pc.rds"))


### SB ----
settings <- get_settings()
settings$save_plots <- F
settings$method_ftemp <- k19_or_l02
settings <- k19_from_forc_to_plot(settings, returndf = "update_settings")

df_acc   <- k19_run_acc(df_sd, settings)

if (use_pfts) df_acc <- replace_acclimated_with_pfts(df_acc, 
                                                     replace_pc = T, 
                                                     replace_xi = F,
                                                     jvr_method = settings$method_ftemp,
                                                     vcmax25_source = vcmax25_source)

df_inst_sb    <- k19_run_inst(df_acc, settings)
df_final_sb <- k19_create_final_df(df_inst_sb, settings)
saveRDS(list(df = df_final_sb, set = settings), here(dir_mods, "df_final_sb.rds"))

### ER ----
settings <- get_settings()
settings$save_plots <- F
settings <- k19_from_forc_to_plot(settings, returndf = "update_settings")

df_acc   <- k19_run_acc(df_sd, settings)

if (use_pfts) df_acc <- replace_acclimated_with_pfts(df_acc, 
                                                     replace_pc = T, 
                                                     replace_xi = T,
                                                     jvr_method = "kumarathunge19_fixed",
                                                     vcmax25_source = vcmax25_source)

df_inst_er    <- k19_run_inst(df_acc, settings)
df_final_er <- k19_create_final_df(df_inst_er, settings)
saveRDS(list(df = df_final_er, set = settings), here(dir_mods, "df_final_er.rds"))


### ER + SB ----
settings <- get_settings()
settings$save_plots <- F

df_acc   <- k19_run_acc(df_sd, settings)

if (use_pfts) df_acc <- replace_acclimated_with_pfts(df_acc, 
                                                     replace_pc = T, 
                                                     replace_xi = F,
                                                     jvr_method = "kumarathunge19_fixed",
                                                     vcmax25_source = vcmax25_source)

df_inst_er_sb    <- k19_run_inst(df_acc, settings)
df_final_er_sb <- k19_create_final_df(df_inst_er_sb, settings)
saveRDS(list(df = df_final_er_sb, set = settings), here(dir_mods, "df_final_er_sb.rds"))

### SB + PC ----
settings <- get_settings()
settings$save_plots <- F
settings$method_ftemp <- k19_or_l02

df_acc   <- k19_run_acc(df_sd, settings)
df_inst_sb_pc    <- k19_run_inst(df_acc, settings)
df_final_sb_pc <- k19_create_final_df(df_inst_sb_pc, settings)
saveRDS(list(df = df_final_sb_pc, set = settings), here(dir_mods, "df_final_sb_pc.rds"))

### ER + PC ----
settings <- get_settings()
settings$save_plots <- F

df_acc   <- k19_run_acc(df_sd, settings)

if (use_pfts) df_acc <- replace_acclimated_with_pfts(df_acc,
                                                     replace_pc = F,
                                                     replace_xi = T,
                                                     jvr_method = "kumarathunge19_fixed",
                                                     vcmax25_source = vcmax25_source)

df_inst_er_pc    <- k19_run_inst(df_acc, settings)
df_final_er_pc <- k19_create_final_df(df_inst_er_pc, settings)
saveRDS(list(df = df_final_er_pc, set = settings), here(dir_mods, "df_final_er_pc.rds"))

message("✅ All models successfully run and saved!")

# Load model runs ----
all_files <- list.files(here(dir_mods), pattern = ".rds")
plot_names <- paste0("p_", str_remove_all(all_files, ".rds|df_final_"))
out_names  <- paste0("out_", str_remove_all(all_files, ".rds|df_final_"))

for (i in 1:length(all_files)) {
  tmp    <- readRDS(here(dir_mods, all_files[i]))
  assign(out_names[i],  tmp)
  tmp$set$dir_now <- "automatic" # Needs fix later
  assign(plot_names[i], plot_all_final_plots_from_df_plot(tmp$df, tmp$set))
}

# Make Output ----
## Generic ----

# Isolated mod-obs
p_modobs <- list()

for (t in c("topt", "aopt", "tspan")) {
  message("Working on ", t)
  p_modobs[[t]] <- 
    plot_iso_acc_modobs(t,
                        p_noacc_scaled = p_noacc_scaled,
                        p_fullacc = p_fullacc,
                        p_k19 = p_er,
                        p_xi = p_sb,
                        p_vj_phi = p_pc,
                        p_k19_vj_phi = p_er_pc,
                        p_xi_k19 = p_er_sb,
                        p_xi_vj_phi = p_sb_pc,
                        dir_figs,
                        modobs_or_modtcair = "modobs")
}

# Isolated trait-tgrowth
p_modtcair <- list()
for (t in c("topt", "aopt", "tspan")) {
  p_modtcair[[t]] <- 
    plot_iso_acc_modobs(t,
                        p_noacc_scaled = p_noacc_scaled,
                        p_fullacc = p_fullacc,
                        p_k19 = p_er,
                        p_xi = p_sb,
                        p_vj_phi = p_pc,
                        p_k19_vj_phi = p_er_pc,
                        p_xi_k19 = p_er_sb,
                        p_xi_vj_phi = p_sb_pc,
                        dir_figs,
                        modobs_or_modtcair = "modtcair")
}


# Individual plots from isolated dfs
if (!dir.exists(here(dir_figs, "response_curves_per_climate"))) dir.create(here(dir_figs, "response_curves_per_climate"), recursive = T, showWarnings = F)
if (!dir.exists(here(dir_figs, "traits_vs_tgrowth"))) dir.create(here(dir_figs, "traits_vs_tgrowth"), recursive = T, showWarnings = F)
if (!dir.exists(here(dir_figs, "jvr_patterns"))) dir.create(here(dir_figs, "jvr_patterns"), recursive = T, showWarnings = F)

# Response Curves
# Add labels
add_labels <- 
  geom_text(
    data = 
      data.frame(
        cl   = c('A', 'C', 'D', 'ET'),
        label   = c('(a)', '(b)', '(c)', '(d)')
      ) %>% 
      mutate(across(cl, factor, levels=c('ET', 'D', 'C', 'A'))),
    aes(y = 135, 
        x = -18, 
        label = label),
    fontface = 'bold', 
    size = 4)

# Curves
ggsave(here(dir_figs, "response_curves_per_climate/curve-shape-per-climate-noacc_scaled.pdf"), p_noacc_scaled$p_scaled_cl_2 + add_labels, height = 5, width = 10)
ggsave(here(dir_figs, "response_curves_per_climate/curve-shape-per-climate-fullacc.pdf"), p_fullacc$p_scaled_cl_2 + add_labels, height = 4, width = 10)
ggsave(here(dir_figs, "response_curves_per_climate/curve-shape-per-climate-er.pdf"), p_er$p_scaled_cl_2 + add_labels, height = 5, width = 10)
ggsave(here(dir_figs, "response_curves_per_climate/curve-shape-per-climate-sb.pdf"), p_sb$p_scaled_cl_2 + add_labels, height = 5, width = 10)
ggsave(here(dir_figs, "response_curves_per_climate/curve-shape-per-climate-pc.pdf"), p_pc$p_scaled_cl_2 + add_labels, height = 5, width = 10)
ggsave(here(dir_figs, "response_curves_per_climate/curve-shape-per-climate-er_pc.pdf"), p_er_pc$p_scaled_cl_2 + add_labels, height = 5, width = 10)
ggsave(here(dir_figs, "response_curves_per_climate/curve-shape-per-climate-er_sb.pdf"), p_er_sb$p_scaled_cl_2 + add_labels, height = 5, width = 10)
ggsave(here(dir_figs, "response_curves_per_climate/curve-shape-per-climate-sb_pc.pdf"), p_sb_pc$p_scaled_cl_2 + add_labels, height = 5, width = 10)

ggsave(here(dir_figs, "response_curves_per_climate/parabola_curve-shape-per-climate-noacc_scaled.pdf"), p_noacc_scaled$p_scaled_cl_2_parab + add_labels, height = 5, width = 10)
ggsave(here(dir_figs, "response_curves_per_climate/parabola_curve-shape-per-climate-fullacc.pdf"), p_fullacc$p_scaled_cl_2_parab + add_labels, height = 5, width = 10)
ggsave(here(dir_figs, "response_curves_per_climate/parabola_curve-shape-per-climate-er.pdf"), p_er$p_scaled_cl_2_parab + add_labels, height = 5, width = 10)
ggsave(here(dir_figs, "response_curves_per_climate/parabola_curve-shape-per-climate-sb.pdf"), p_sb$p_scaled_cl_2_parab + add_labels, height = 5, width = 10)
ggsave(here(dir_figs, "response_curves_per_climate/parabola_curve-shape-per-climate-pc.pdf"), p_pc$p_scaled_cl_2_parab + add_labels, height = 5, width = 10)
ggsave(here(dir_figs, "response_curves_per_climate/parabola_curve-shape-per-climate-er_pc.pdf"), p_er_pc$p_scaled_cl_2_parab + add_labels, height = 5, width = 10)
ggsave(here(dir_figs, "response_curves_per_climate/parabola_curve-shape-per-climate-er_sb.pdf"), p_er_sb$p_scaled_cl_2_parab + add_labels, height = 5, width = 10)
ggsave(here(dir_figs, "response_curves_per_climate/parabola_curve-shape-per-climate-sb_pc.pdf"), p_sb_pc$p_scaled_cl_2_parab + add_labels, height = 5, width = 10)

# Traits
ggsave(here(dir_figs, "traits_vs_tgrowth/traits-tgrowth_noacc_scaled.pdf"), p_noacc_scaled$p_traits + plot_annotation(subtitle = "No acc."), height = 4, width = 8)
ggsave(here(dir_figs, "traits_vs_tgrowth/traits-tgrowth_fullacc.pdf"), p_fullacc$p_traits + plot_annotation(subtitle = "Full acc."), height = 4, width = 8)
ggsave(here(dir_figs, "traits_vs_tgrowth/traits-tgrowth_er.pdf"), p_er$p_traits + plot_annotation(subtitle = "ER"), height = 4, width = 8)
ggsave(here(dir_figs, "traits_vs_tgrowth/traits-tgrowth_sb.pdf"), p_sb$p_traits + plot_annotation(subtitle = "SB"), height = 4, width = 8)
ggsave(here(dir_figs, "traits_vs_tgrowth/traits-tgrowth_pc.pdf"), p_pc$p_traits + plot_annotation(subtitle = "PC"), height = 4, width = 8)
ggsave(here(dir_figs, "traits_vs_tgrowth/traits-tgrowth_er_pc.pdf"), p_er_pc$p_traits + plot_annotation(subtitle = "ER + PC"), height = 4, width = 8)
ggsave(here(dir_figs, "traits_vs_tgrowth/traits-tgrowth_er_sb.pdf"), p_er_sb$p_traits + plot_annotation(subtitle = "ER + SB"), height = 4, width = 8)
ggsave(here(dir_figs, "traits_vs_tgrowth/traits-tgrowth_sb_pc.pdf"), p_sb_pc$p_traits + plot_annotation(subtitle = "SB + PC"), height = 4, width = 8)

# JVr
ggsave(here(dir_figs, "jvr_patterns/jvr_noacc_scaled.pdf"), p_noacc_scaled$p_jvr_all + plot_annotation(subtitle = "No acc."), height = 4, width = 8)
ggsave(here(dir_figs, "jvr_patterns/jvr_fullacc.pdf"), p_fullacc$p_jvr_all + plot_annotation(subtitle = "Full acc."), height = 4, width = 8)
ggsave(here(dir_figs, "jvr_patterns/jvr_er.pdf"), p_er$p_jvr_all + plot_annotation(subtitle = "ER"), height = 4, width = 8)
ggsave(here(dir_figs, "jvr_patterns/jvr_sb.pdf"), p_sb$p_jvr_all + plot_annotation(subtitle = "SB"), height = 4, width = 8)
ggsave(here(dir_figs, "jvr_patterns/jvr_pc.pdf"), p_pc$p_jvr_all + plot_annotation(subtitle = "PC"), height = 4, width = 8)
ggsave(here(dir_figs, "jvr_patterns/jvr_er_pc.pdf"), p_er_pc$p_jvr_all + plot_annotation(subtitle = "ER + PC"), height = 4, width = 8)
ggsave(here(dir_figs, "jvr_patterns/jvr_er_sb.pdf"), p_er_sb$p_jvr_all + plot_annotation(subtitle = "ER + SB"), height = 4, width = 8)
ggsave(here(dir_figs, "jvr_patterns/jvr_sb_pc.pdf"), p_sb_pc$p_jvr_all + plot_annotation(subtitle = "SB + PC"), height = 4, width = 8)

## Tables ----
# Topt - Tgrowth
tables_acc <- list()
for (t in c("topt", "aopt", "tspan")) {
  tables_acc[[t]] <- 
    table_acclimation_contributions(t,
                                    p_noacc_scaled = p_noacc_scaled,
                                    p_fullacc = p_fullacc,
                                    p_k19 = p_er,
                                    p_xi = p_sb,
                                    p_vj_phi = p_pc,
                                    p_k19_vj_phi = p_er_pc,
                                    p_xi_k19 = p_er_sb,
                                    p_xi_vj_phi = p_sb_pc,
                                    dir_tabs,
                                    comparison = "trait-temp")
}

tables_modobs <- list()
for (t in c("topt", "aopt", "tspan")) {
  tables_modobs[[t]] <- 
    table_acclimation_contributions(t,
                                    p_noacc_scaled = p_noacc_scaled,
                                    p_fullacc = p_fullacc,
                                    p_k19 = p_er,
                                    p_xi = p_sb,
                                    p_vj_phi = p_pc,
                                    p_k19_vj_phi = p_er_pc,
                                    p_xi_k19 = p_er_sb,
                                    p_xi_vj_phi = p_sb_pc,
                                    dir_tabs,
                                    comparison = "mod-obs")
}

## Specific ----
### 3x3 Prediced - Observed ----
no_x_axis  <- theme_classic() + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
no_y_axis  <- theme_classic() + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
no_xy_axis <- theme_classic() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank())    

p <-
  ## Top row
  p_noacc_scaled$p_topt_mod_tcair   + no_x_axis  + annotate("text", x = 38, y = 2, label = "(a)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "Modelled (no acclimation)", y = bquote(T[opt] ~ " [°C]")) +  # "Modelled (No acclimation)"
  p_fullacc$p_topt_mod_tcair        + no_xy_axis + annotate("text", x = 38, y = 2, label = "(b)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "Modelled (full acclimation)") +  # "Modelled (Full acclimation)"
  p_fullacc$p_topt_obs_tcair        + no_xy_axis + annotate("text", x = 38, y = 2, label = "(c)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "Observations") +  # "Observed"
  ## Mid Row
  p_noacc_scaled$p_aopt_mod_tcair   + no_x_axis  + ylim(0, 22) + annotate("text", x = 38, y = 2, label = "(d)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, y = bquote(A[opt]  ~ "[µmol" ~ m ^-2 ~ s ^-1 ~ "]")) +
  p_fullacc$p_aopt_mod_tcair        + no_xy_axis + ylim(0, 22) + annotate("text", x = 38, y = 2, label = "(e)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL) +
  p_fullacc$p_aopt_obs_tcair        + no_xy_axis + ylim(0, 22) + annotate("text", x = 38, y = 2, label = "(f)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL) +
  ## Bottom Row
  p_noacc_scaled$p_tspan_mod_tcair  + theme_classic() + ylim(0, 22) + annotate("text", x = 38, y = 2, label = "(g)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, y = bquote(T[span] ~ " [°C]")) +
  p_fullacc$p_tspan_mod_tcair       + no_y_axis  + ylim(0, 22) + annotate("text", x = 38, y = 2, label = "(h)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL) +
  p_fullacc$p_tspan_obs_tcair       + no_y_axis  + ylim(0, 22) + annotate("text", x = 38, y = 2, label = "(i)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL) +
  
  ## Aesthetics
  plot_layout(guides = "collect",
              nrow = 3) &
  # plot_annotation(
    # title = expression(bold("Modelled and observed patterns of thermal acclimation")),
    # subtitle = ""
    # ) &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0, size=16),
        plot.subtitle = element_text(hjust = 0.5, face="bold", size=12),
        axis.title = element_text(size = 12))

ggsave(here(dir_figs, "modelled-observed-patterns_3x3.pdf"),
       p,
       height = 7.75,
       width = 7.75)

ggsave(here(dir_figs, "png-for-github_modelled-observed-patterns_3x3.png"),
       p,
       height = 7.75,
       width = 7.75)


### 3x3 Mod-Obs Topt ----
p <-
  ## Top row
  p_noacc_scaled$modobs_topt   + no_x_axis  + annotate("text", x = 38, y = 2, label = "(a)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "No acclimation") +  # "Modelled (No acclimation)"
  p_fullacc$modobs_topt        + no_xy_axis + annotate("text", x = 38, y = 2, label = "(b)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "Full acclimation") +  # "Modelled (Full acclimation)"
  plot_spacer() +
  ## Mid Row
  p_pc$modobs_topt             + no_x_axis   + annotate("text", x = 38, y = 2, label = "(c)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "PC") +
  p_sb$modobs_topt             + no_xy_axis  + annotate("text", x = 38, y = 2, label = "(d)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "SS") +
  p_er$modobs_topt             + no_xy_axis  + annotate("text", x = 38, y = 2, label = "(e)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "ER") +
  ## Bottom Row
  p_sb_pc$modobs_topt          + theme_classic() + annotate("text", x = 38, y = 2, label = "(f)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "PC + SS") +
  p_er_pc$modobs_topt          + no_y_axis + annotate("text", x = 38, y = 2, label = "(g)", fontface = "bold", size = 10/.pt)       + labs(title = "", caption = NULL, subtitle = "PC + ER") +
  p_er_sb$modobs_topt          + no_y_axis + annotate("text", x = 38, y = 2, label = "(h)", fontface = "bold", size = 10/.pt)       + labs(title = "", caption = NULL, subtitle = "SS + ER") +
  
  ## Aesthetics
  plot_layout(guides = "collect",
              nrow = 3) &
  # plot_annotation(title = expression(bold("Performance of model setups to predict" ~ T[opt])),
                  # subtitle = "") &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0, size=16),
        plot.subtitle = element_text(hjust = 0, face="bold", size=12),
        axis.title = element_text(size = 12))

ggsave(here(dir_figs, "mod-obs_topt_3x3.pdf"),
       p,
       height = 7.75,
       width = 7.75)

### 3x3 Aopt-Tgrowth ----
p <-
  ## Top row
  p_noacc_scaled$p_aopt_mod_tcair   + no_x_axis +  annotate("text", x = 38, y = 2, label = "(a)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "No acclimation",   y = bquote(A[opt]  ~ "[µmol" ~ m ^-2 ~ s ^-1 ~ "]")) +  # "Modelled (No acclimation)"
  p_fullacc$p_aopt_mod_tcair        + no_xy_axis + annotate("text", x = 38, y = 2, label = "(b)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "Full acclimation") +  # "Modelled (Full acclimation)"
  p_fullacc$p_aopt_obs_tcair        + no_xy_axis + annotate("text", x = 38, y = 2, label = "(c)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "Observations") +  # "Modelled (Full acclimation)"
  ## Mid Row
  p_pc$p_aopt_mod_tcair             + no_x_axis +  annotate("text", x = 38, y = 2, label = "(c)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "PC", y = bquote(A[opt]  ~ "[µmol" ~ m ^-2 ~ s ^-1 ~ "]")) +
  p_sb$p_aopt_mod_tcair             + no_xy_axis + annotate("text", x = 38, y = 2, label = "(d)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "SS") +
  p_er$p_aopt_mod_tcair             + no_xy_axis + annotate("text", x = 38, y = 2, label = "(e)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "ER") +
  ## Bottom Row
  p_sb_pc$p_aopt_mod_tcair          + theme_classic() + annotate("text", x = 38, y = 2, label = "(f)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "PC + SS", y = bquote(A[opt]  ~ "[µmol" ~ m ^-2 ~ s ^-1 ~ "]")) +
  p_er_pc$p_aopt_mod_tcair          + no_y_axis + annotate("text", x = 38, y = 2, label = "(g)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "PC + ER") +
  p_er_sb$p_aopt_mod_tcair          + no_y_axis + annotate("text", x = 38, y = 2, label = "(h)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "SS + ER") +
  
  ## Aesthetics
  plot_layout(guides = "collect",
              nrow = 3) &
  # plot_annotation(title = expression(bold("Predicted and observed thermal acclimation of" ~ A[opt])),
  #                 subtitle = "") &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0, size=16),
        plot.subtitle = element_text(hjust = 0, face="bold", size=12),
        axis.title = element_text(size = 12))

ggsave(here(dir_figs, "aopt-tgrowth_3x3.pdf"),
       p,
       height = 7.75,
       width = 7.75)


### Thermacc of traits (vcmax25, jmax25, xi) ----
# PFT Fixed
p_pft <- 
  (
    p_noacc_scaled$p_traits_list$vcmax25 + labs(subtitle = NULL, y = bquote(V[cmax]^25 ~ "[µmol" ~ m^-2 ~ s ^-1 ~ "]")) +
    p_noacc_scaled$p_traits_list$jmax25  + labs(subtitle = NULL, y = bquote(J[max]^25 ~ "[µmol" ~ m^-2 ~ s ^-1 ~ "]")) +
    p_noacc_scaled$p_traits_list$xi      + labs(subtitle = NULL, y = bquote(xi ~ "[" ~ Pa^-0.5 ~ "]"))
    ) +
  plot_layout(guides = "collect") &
  # plot_annotation(title = expression("PFT-specific values for " ~ V[cmax]^25 ~ "," ~ J[max]^25 ~ ", and" ~ xi)) &
  theme_classic() &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0, size=16),
        plot.subtitle = element_text(hjust = 0, face="bold", size=12),
        axis.title = element_text(size = 12))

p_pft

ggsave(here(dir_figs, "traits-pft-specific.pdf"),
       p_pft,
       height = 4,
       width = 10)

# Optimality based
#  --- Bugfix for separate optimality plot with abc.)
p_opt <- 
  (
    p_fullacc$p_traits_list$vcmax25   + annotate("text", x = 35, y = 2, label = "(a)", fontface = "bold", size = 10/.pt) + labs(subtitle = NULL, y = bquote(V[cmax]^25 ~ "[µmol" ~ m^-2 ~ s ^-1 ~ "]")) +
      p_fullacc$p_traits_list$jmax25  + annotate("text", x = 35, y = 4, label = "(b)", fontface = "bold", size = 10/.pt) + labs(subtitle = NULL, y = bquote(J[max]^25 ~ "[µmol" ~ m^-2 ~ s ^-1 ~ "]")) +
      p_fullacc$p_traits_list$xi      + annotate("text", x = 35, y = 2, label = "(c)", fontface = "bold", size = 10/.pt) + labs(subtitle = NULL, y = bquote(xi ~ "[" ~ Pa^-0.5 ~ "]"))
  ) +
  plot_layout(guides = "collect") &
  theme_classic() &
  labs(fill = "Plant Functional Type") &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0, size=16),
        plot.subtitle = element_text(hjust = 0, face="bold", size=12),
        axis.title = element_text(size = 12))

ggsave(here(dir_figs, "traits-optimality-predicted.pdf"),
       p_opt,
       height = 4,
       width = 10)

#  --- Original code again
p_opt <- 
  (
    p_fullacc$p_traits_list$vcmax25   + annotate("text", x = 2, y = 150, label = "(d)", fontface = "bold", size = 10/.pt) + labs(subtitle = NULL, y = bquote(V[cmax]^25 ~ "[µmol" ~ m^-2 ~ s ^-1 ~ "]")) +
      p_fullacc$p_traits_list$jmax25  + annotate("text", x = 2, y = 300, label = "(e)", fontface = "bold", size = 10/.pt) + labs(subtitle = NULL, y = bquote(J[max]^25 ~ "[µmol" ~ m^-2 ~ s ^-1 ~ "]")) +
      p_fullacc$p_traits_list$xi      + annotate("text", x = 2, y = 150, label = "(f)", fontface = "bold", size = 10/.pt) + labs(subtitle = NULL, y = bquote(xi ~ "[" ~ Pa^-0.5 ~ "]"))
  ) +
  plot_layout(guides = "collect") &
  # plot_annotation(title = expression("Optimaliy-prediced values for " ~ V[cmax]^25 ~ "," ~ J[max]^25 ~ ", and" ~ xi)) &
  theme_classic() &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0, size=16),
        plot.subtitle = element_text(hjust = 0, face="bold", size=12),
        axis.title = element_text(size = 12))

# p_opt

ggsave(here(dir_figs, "traits-optimality-predicted_no-labels.pdf"),
       p_opt,
       height = 4,
       width = 10)

# Both together
p_pft <- 
  (
    p_er$p_traits_list$vcmax25   + annotate("text", x = 2, y = 150, label = "(a)", fontface = "bold", size = 10/.pt) + labs(subtitle = NULL, y = bquote(V[cmax]^25 ~ "[µmol" ~ m^-2 ~ s ^-1 ~ "]")) +
      p_er$p_traits_list$jmax25  + annotate("text", x = 2, y = 300, label = "(b)", fontface = "bold", size = 10/.pt) + labs(subtitle = NULL, y = bquote(J[max]^25 ~ "[µmol" ~ m^-2 ~ s ^-1 ~ "]")) +
      p_er$p_traits_list$xi      + annotate("text", x = 2, y = 150, label = "(c)", fontface = "bold", size = 10/.pt) + labs(subtitle = NULL, y = bquote(xi ~ "[" ~ Pa^-0.5 ~ "]"))
  ) +
  plot_layout(guides = "collect") &
  # plot_annotation(title = expression("PFT-specific values for " ~ V[cmax]^25 ~ "," ~ J[max]^25 ~ ", and" ~ xi)) &
  # theme_classic(base_size = 16) &
  theme_classic(base_size = 12) &
  theme(legend.position = "none")

p_both <- 
  p_pft / p_opt &
  # plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  labs(fill ='Plant Functional Type') &
  theme(plot.tag = element_text(face = "bold"))

p_both

ggsave(here(dir_figs, "traits-both.pdf"),
       p_both,
       height = 6,
       width = 10)
 
### 3x3 JVr ~ Tgrowth ----
no_x_axis  <- theme_classic() + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
no_y_axis  <- theme_classic() + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
no_xy_axis <- theme_classic() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank())    

p <-
  ## Top row
  p_noacc_scaled$p_jvr_tgrowth   + no_x_axis  + annotate("text", x = 1, y = 3, label = "", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "No acclimation", y = bquote(JV[r])) + ylim(1.25, 3) + xlim(0, 40) +  # "Modelled (No acclimation)"
  p_fullacc$p_jvr_tgrowth        + no_xy_axis + annotate("text", x = 1, y = 3, label = "", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "Full acclimation") + ylim(1.25, 3) + xlim(0, 40) +  # "Modelled (Full acclimation)"
  # p_fullacc$p_jvr_tgrowth        + no_xy_axis + annotate("text", x = 1, y = 3, label = "", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "Observations") +  # "Observed"
  plot_spacer() +
  ## Mid Row
  p_pc$p_jvr_tgrowth   + no_x_axis  + annotate("text", x = 1, y = 3, label = "", fontface = "bold", size = 10/.pt) + labs(subtitle = "PC", caption = NULL, y = bquote(JV[r])) + ylim(1.25, 3) + xlim(0, 40) +
  p_sb$p_jvr_tgrowth   + no_xy_axis + annotate("text", x = 1, y = 3, label = "", fontface = "bold", size = 10/.pt) + labs(subtitle = "SS", caption = NULL) + ylim(1.25, 3) + xlim(0, 40) +
  p_er$p_jvr_tgrowth   + no_xy_axis + annotate("text", x = 1, y = 3, label = "", fontface = "bold", size = 10/.pt) + labs(subtitle = "ER", caption = NULL) + ylim(1.25, 3) + xlim(0, 40) +
  ## Bottom Row
  p_sb_pc$p_jvr_tgrowth  + theme_classic() + annotate("text", x = 1, y = 3, label = "", fontface = "bold", size = 10/.pt) + labs(subtitle = "PC + SS", caption = NULL, y = bquote(JV[r])) + ylim(1.25, 3) + xlim(0, 40) +
  p_er_pc$p_jvr_tgrowth  + no_y_axis  + annotate("text", x = 1, y = 3, label = "", fontface = "bold", size = 10/.pt)      + labs(subtitle = "PC + ER", caption = NULL) + ylim(1.25, 3) + xlim(0, 40) +
  p_er_sb$p_jvr_tgrowth  + no_y_axis  + annotate("text", x = 1, y = 3, label = "", fontface = "bold", size = 10/.pt)      + labs(subtitle = "SS + ER", caption = NULL) + ylim(1.25, 3) + xlim(0, 40) +
  
  ## Aesthetics
  plot_layout(guides = "collect",
              nrow = 3) &
  # plot_annotation(title = expression(bold("Modelled Response of " ~ JV[r] ~ " to " ~ T[growth])),
  #                 subtitle = "") &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0, size=16),
        plot.subtitle = element_text(hjust = 0, face="bold", size=12),
        axis.title = element_text(size = 12))

ggsave(here(dir_figs, "jvr-tgrowth_3x3.pdf"),
       p,
       height = 7.75,
       width = 7.75)

### 3x3 Topt ~ JVr ----
no_x_axis  <- theme_classic() + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
no_y_axis  <- theme_classic() + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
no_xy_axis <- theme_classic() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank())    

p <-
  ## Top row
  p_noacc_scaled$p_topt_jvr   + no_x_axis  + annotate("text", x = 2.85, y = 2.5, label = "(a)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "No acclimation", x = bquote(JV[r])) + xlim(1.25, 3) + ylim(0, 40) +  # "Modelled (No acclimation)"
  p_fullacc$p_topt_jvr        + no_xy_axis + annotate("text", x = 2.85, y = 2.5, label = "(b)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "Full acclimation") + xlim(1.25, 3) + ylim(0, 40) +  # "Modelled (Full acclimation)"
  # p_fullacc$p_topt_jvr        + no_xy_axis + annotate("text", x = 2.85, y = 2.5, label = "(c)", fontface = "bold", size = 10/.pt) + labs(title = "", caption = NULL, subtitle = "Observations") +  # "Observed"
  plot_spacer() +
  ## Mid Row
  p_pc$p_topt_jvr   + no_x_axis  + annotate("text", x = 2.85, y = 2.5, label = "(c)", fontface = "bold", size = 10/.pt) + labs(subtitle = "PC", caption = NULL, x = bquote(JV[r])) + xlim(1.25, 3) + ylim(0, 40) +
  p_sb$p_topt_jvr   + no_xy_axis + annotate("text", x = 2.85, y = 2.5, label = "(d)", fontface = "bold", size = 10/.pt) + labs(subtitle = "SS", caption = NULL) + xlim(1.25, 3) + ylim(0, 40) +
  p_er$p_topt_jvr   + no_xy_axis + annotate("text", x = 2.85, y = 2.5, label = "(e)", fontface = "bold", size = 10/.pt) + labs(subtitle = "ER", caption = NULL) + xlim(1.25, 3) + ylim(0, 40) +
  ## Bottom Row
  p_sb_pc$p_topt_jvr  + theme_classic() + annotate("text", x = 2.85, y = 2.5, label = "(f)", fontface = "bold", size = 10/.pt) + labs(subtitle = "PC + SS", caption = NULL, x = bquote(JV[r])) + xlim(1.25, 3) + ylim(0, 40) +
  p_er_pc$p_topt_jvr  + no_y_axis  + annotate("text", x = 2.85, y = 2.5, label = "(g)", fontface = "bold", size = 10/.pt)      + labs(subtitle = "PC + ER", caption = NULL) + xlim(1.25, 3) + ylim(0, 40) +
  p_er_sb$p_topt_jvr  + no_y_axis  + annotate("text", x = 2.85, y = 2.5, label = "(h)", fontface = "bold", size = 10/.pt)      + labs(subtitle = "SS + ER", caption = NULL) + xlim(1.25, 3) + ylim(0, 40) +
  
  ## Aesthetics
  plot_layout(guides = "collect",
              nrow = 3) +
  plot_annotation(
    # title = expression(bold("Modelled Response of " ~ T[opt] ~ " to " ~ JV[r])),
    subtitle = ""
    ) &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0, size=16),
        plot.subtitle = element_text(hjust = 0, face="bold", size=12),
        axis.title = element_text(size = 12))

ggsave(here(dir_figs, "jvr-topt_3x3.pdf"),
       p,
       height = 7.75,
       width = 7.75)

beepr::beep()
