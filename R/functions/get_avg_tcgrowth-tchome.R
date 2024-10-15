get_avg_tc_growth <- function(){
  # tc_growth was calculated from:
   # settings <- get_settings()
   # settings$daily_conditions <- "A"
   # k19_create_df_forcing(settings) %>% 
   #   k19_create_input_for_acc(settings, .) |>
   #   unnest(forcing_d) |> 
   #   filter(temp > 5) |> 
   #   pull(temp) |> 
   #   mean()
  
  tcgrowth <- 20.63868
  return(tcgrowth)
}

get_avg_tc_home <- function(){
  # tc_home was calculated from (standard settings):
  # settings <- get_settings()
  # k19_create_df_forcing(settings) %>% k19_create_input_for_acc(settings, .) |> 
  # unnest(siteinfo) |> group_by(sitename) |> summarise(tc_home = mean(tc_home)) |> 
  # pull(tc_home) |> mean()
  
  tchome <- 25.97692
  return(tchome)
}
