make_new_dynamic_dir <- function(analysis  = NA,
                                 file_type = NA,
                                 appendix  = NA,
                                 return_latest  = FALSE) {
  
  if (is.na(analysis)) stop("Specify name of analysis in make_new_dir")
  if (is.na(file_type)) stop("Specify file type in make_new_dir")
  
  if (return_latest) {
    
    dir_top <- here("output", analysis, file_type)
    dir_mid <- tail(list.dirs(dir_top, recursive = F), 1) # Get latest date
    dir_tmp <- tail(list.dirs(dir_mid, recursive = F), 1) # Get latest run
    
  } else {
  
    # Get basic dir structure
    dir_tmp <- here("output", analysis, file_type, Sys.Date())
    
    # print(dir_tmp)
    
    # Check if multiple runs were done for this day and label dynamically
    if (!dir.exists(dir_tmp)) {
      n_runs <- 0
      dir.create(dir_tmp, recursive = T, showWarnings = F)
    } else {
      
      n_runs <- length(list.dirs(dir_tmp, recursive = F))
    }
  
    if (n_runs < 10) {
      n_runs <- paste0(0, n_runs)
    }
    
    # Make naming nice if with or without suffix
    if (is.na(appendix)) {
      dir_tmp <- paste0(dir_tmp, "/", n_runs, "/")
    } else {
      dir_tmp <- paste0(dir_tmp, "/", n_runs, "_", appendix, "/")
    }
    
    # Make final dir and return it
    dir.create(dir_tmp, recursive = T, showWarnings = F)
    
  }
  
  return(dir_tmp)
}
