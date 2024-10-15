# Testing Aopt Prediction Across PFT Parametriations ----
# ISOLATED ACCLIMATION PROCESSES ----

source("R/source.R")

# ├ Settings ----
dir_analysis <- "aopt_pfts_comparison"
dir_tmp      <- "final"
use_pfts     <- T
vec_vcmax25_source <- c("avim", "bethy", "clm", "ibis", "kumar19", "ocn", "orchidee", "mean", "mean_nok19")
pval_threshold <- 0.01

# ├ Define directories ----
dir_mods <- make_new_dynamic_dir(dir_analysis, "model_rds", dir_tmp)
dir_figs <- make_new_dynamic_dir(dir_analysis, "figures", dir_tmp)
dir_tabs <- make_new_dynamic_dir(dir_analysis, "tables", dir_tmp)

dir_mods <- make_new_dynamic_dir(dir_analysis, "model_rds", dir_tmp, return_latest = TRUE)
dir_figs <- make_new_dynamic_dir(dir_analysis, "figures", dir_tmp, return_latest = TRUE)
dir_tabs <- make_new_dynamic_dir(dir_analysis, "tables", dir_tmp, return_latest = TRUE)

## Run Setups ----
### PC+SS+ER Setup ----
settings <- get_settings()
settings$save_plots <- F
settings$method_ftemp <- "kumarathunge19"
df_forc  <- k19_create_df_forcing(settings)
df_sd    <- k19_create_input_for_acc(settings, df_forc)
df_acc   <- k19_run_acc(df_sd, settings)
df_inst_fullacc  <- k19_run_inst(df_acc, settings)
df_final_fullacc <- k19_create_final_df(df_inst_fullacc, settings)

# Check if change in Aopt is significant
pvalue = summary(lm(df_final_fullacc$rpm_sim_anet_opt ~ df_final_fullacc$temp))$coef[2, 4] 
df_final_fullacc[["pvalue_aopt"]] <- pvalue

pvalue = summary(lm(df_final_fullacc$rpm_sim_tc_opt ~ df_final_fullacc$temp))$coef[2, 4] 
df_final_fullacc[["pvalue_topt"]] <- pvalue

pvalue = summary(lm(df_final_fullacc$rpm_sim_tspan ~ df_final_fullacc$temp))$coef[2, 4] 
df_final_fullacc[["pvalue_tspan"]] <- pvalue

# Add information
df_final_fullacc[["vcmax25_source"]] = "optimality"
df_final_fullacc <- df_final_fullacc |> mutate(pft = purrr::map_chr(data_org, ~pull(., PFT) %>% unique()))

### ER setup for different vcmax25 sources ----
settings <- get_settings()
settings$save_plots <- F
settings$method_ftemp <- "kumarathunge19"
df_forc  <- k19_create_df_forcing(settings)
df_sd    <- k19_create_input_for_acc(settings, df_forc)

df_out <- df_final_fullacc
for (vcmax25_source_in in vec_vcmax25_source) {
  
  # Verbose
  message("\n - Working on \t", vcmax25_source_in, "\n")
  
  df_acc   <- k19_run_acc(df_sd, settings)
  if (use_pfts) df_acc <- replace_acclimated_with_pfts(df_acc, 
                                                       replace_pc = T, 
                                                       replace_xi = F,
                                                       jvr_method = settings$method_ftemp,
                                                       vcmax25_source = vcmax25_source_in)
  
  df_inst_noacc  <- k19_run_inst(df_acc, settings)
  df_final_noacc <- k19_create_final_df(df_inst_noacc, settings)
  # Add source and pft info
  df_final_noacc[["vcmax25_source"]] = vcmax25_source_in
  # Extract PFT information:
  df_final_noacc <- df_final_noacc |> mutate(pft = purrr::map_chr(data_org, ~pull(., PFT) %>% unique()))
  # Remove PFTs where model had NA values
  avim_nas     <- c("Arctic")
  bethy_nas    <- c("BET-Tr", "NET-Te")
  ibis_nas     <- c("Arctic")
  ocn_nas      <- c("Arctic", "BDT-Te")
  orchidee_nas <- c("Arctic", "BDT-Te")
  df_final_noacc <- df_final_noacc |>
    filter(
     !(vcmax25_source == "avim"     & pft %in% avim_nas),
     !(vcmax25_source == "bethy"    & pft %in% bethy_nas),
     !(vcmax25_source == "ibis"     & pft %in% ibis_nas),
     !(vcmax25_source == "ocn"      & pft %in% ocn_nas),
     !(vcmax25_source == "orchidee" & pft %in% orchidee_nas),
    )

  # Check if change in Aopt is significant
  pvalue = summary(lm(df_final_noacc$rpm_sim_anet_opt ~ df_final_noacc$temp))$coef[2, 4] 
  df_final_noacc[["pvalue_aopt"]] <- pvalue
  
  pvalue = summary(lm(df_final_noacc$rpm_sim_tc_opt ~ df_final_noacc$temp))$coef[2, 4] 
  df_final_noacc[["pvalue_topt"]] <- pvalue
  
  pvalue = summary(lm(df_final_noacc$rpm_sim_tspan ~ df_final_noacc$temp))$coef[2, 4] 
  df_final_noacc[["pvalue_tspan"]] <- pvalue
  
  # Concatenate dataframe
  df_out = rbind(df_out, df_final_noacc)
}

# Attach pvalue info
df_out[["aopt_temp_sign"]] <- NA
df_out[["topt_temp_sign"]] <- NA
df_out[["tspan_temp_sign"]] <- NA

for (i in 1:nrow(df_out)) {
  # Aopt
  if (df_out$pvalue_aopt[i] > pval_threshold) {
    df_out$"aopt_temp_sign"[i] = FALSE
  } else {
    df_out$"aopt_temp_sign"[i] = TRUE
  }
  # Topt
  if (df_out$pvalue_topt[i] > pval_threshold) {
    df_out$"topt_temp_sign"[i] = FALSE
  } else {
    df_out$"topt_temp_sign"[i] = TRUE
  }
  # Tspan
  if (df_out$pvalue_tspan[i] > pval_threshold) {
    df_out$"tspan_temp_sign"[i] = FALSE
  } else {
    df_out$"tspan_temp_sign"[i] = TRUE
  }
}

df_out$aopt_temp_sign <- as.factor(df_out$aopt_temp_sign)
df_out$topt_temp_sign <- as.factor(df_out$topt_temp_sign)
df_out$tspan_temp_sign <- as.factor(df_out$tspan_temp_sign)

# Number of unique vcmax25_source categories and Generate palette
num_categories <- length(unique(df_out$vcmax25_source))
all_colors <- pals::kelly()
my_pal <- all_colors[-1]
my_pal <- my_pal[1:num_categories]

# Aopt ----
## Plot ----

paopt <- 
  ggplot(df_out) + 
  aes(x = temp,
      y = rpm_sim_anet_opt, 
      color = vcmax25_source, 
      # shape = 4,
      linetype = aopt_temp_sign) +
  # geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +
  # scale_color_brewer(palette = "Paired") +
  scale_color_manual(values = my_pal,
                     name = bquote("Source of "~V[cmax]^25)) +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dotdash"),
                        labels = c("TRUE" = "Yes", "FALSE" = "No"),
                        name = paste("Significant at ", pval_threshold)) +
  guides(linetype = guide_legend(override.aes = list(color = "black"))) +
  xlab(expression(T[growth] ~ "[°C]")) + 
  ylab(expression("Modelled" ~ A[opt] ~ "[µmol" ~ m ^-2 ~ s ^-1 ~ "]")) +
  xlim(0, 40) +
  ylim(0, 40) +
  theme_classic(base_size = 16)

paopt

## Save ----
ggsave(
  here(dir_figs, "/aopt_versus_tgrowth_for_all_vcmax25_sources.pdf"), 
  paopt,
  height = 5, 
  width = 7
  )

# Topt ----
## Plot ----

ptopt <- 
  ggplot(df_out) + 
  aes(x = temp,
      y = rpm_sim_tc_opt, 
      # linetype = tc_opt_temp_sign,
      color = vcmax25_source
      # shape = 4,
  ) +
  # geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +
  # scale_color_brewer(palette = "Paired") +
  scale_color_manual(values = my_pal,
                     name = bquote("Source of "~V[cmax]^25)) +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dotdash"),
                        labels = c("TRUE" = "Yes", "FALSE" = "No"),
                        name = paste("Significant at ", pval_threshold)) +
  guides(linetype = guide_legend(override.aes = list(color = "black"))) +
  xlab(expression(T[growth] ~ "[°C]")) + 
  ylab(expression("Modelled" ~ T[opt] ~ "[°C]")) +
  xlim(0, 40) +
  ylim(0, 40) +
  theme_classic(base_size = 16)

ptopt

## Save ----
ggsave(
  here(dir_figs, "/topt_versus_tgrowth_for_all_vcmax25_sources.pdf"), 
  ptopt,
  height = 5, 
  width = 7
)


# Tspan ----
## Plot ----
ptspan <- 
  ggplot(df_out) + 
  aes(x = temp,
      y = rpm_sim_tspan, 
      # linetype = tspan_temp_sign,
      color = vcmax25_source
      # shape = 4,
      ) +
  # geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +
  # scale_color_brewer(palette = "Paired") +
  scale_color_manual(values = my_pal,
                     name = bquote("Source of "~V[cmax]^25)) +
  # scale_linetype_manual(
  #   values = c("TRUE" = "solid", "FALSE" = "dotted"),  # Include additional level
  #   breaks = c(levels(df_out$tspan_temp_sign), FALSE),  # Include additional level
  #   labels = c("TRUE" = "Yes", "FALSE" = "No"),
  #   name = paste("Significant at ", pval_threshold)
  # ) +
  # guides(linetype = guide_legend(override.aes = list(color = "black"))) +
  guides(linetype = NULL) +
  xlab(expression(T[growth] ~ "[°C]")) + 
  ylab(expression("Modelled" ~ T[span] ~ "[°C]")) +
  xlim(0, 40) +
  ylim(0, 40) +
  theme_classic(base_size = 16)

ptspan

## Save ----
ggsave(
  here(dir_figs, "/tspan_versus_tgrowth_for_all_vcmax25_sources.pdf"), 
  ptspan,
  height = 5, 
  width = 7
)

# Patchwork ----
pp <- 
  (paopt +  annotate("text", x = 3, y = 39, label = "(a)", fontface = "bold", size = 10/.pt)) +   
  (ptopt +  annotate("text", x = 3, y = 39, label = "(b)", fontface = "bold", size = 10/.pt)) +  
  (ptspan + annotate("text", x = 3, y = 39, label = "(c)", fontface = "bold", size = 10/.pt)) +    
  plot_layout(guides = "collect", nrow = 1)

ggsave(
  here(dir_figs, "/all_traits_per_vcmax25source.pdf"), 
  pp,
  height = 5, 
  width = 12
)
  

