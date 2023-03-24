# Function to fit Anet-Tleaf parabola to observations ----
## Disclaimer ----
#' This code was adapted from the conde by Kumarathunge et al. 2019 (DOI: 10.1111/nph.15668)
#' Original code available online: 
#' https://bitbucket.org/Kumarathunge/photom/src/master/R/Functions/functions_for_analysis.R

fit_nonlinear_topt <- function(dat, x = "", y = "", random = ""){
  
  ## Checks ####
  # Catching data aggregations with too little entries
  if (nrow(dat) <= 3) {
    fit_method <- "NA"
    out <- tibble(aopt = NA,
                  aopt.se = NA,
                  topt = NA,
                  topt.se = NA,
                  b = NA,
                  b.se = NA,
                  r2 = NA,
                  s = NA,
                  pvalue = NA,
                  fit_method = fit_method)
    return(out)
  }
  
  # Define b as negative for if-checks later due to nlme() algorithm
  # being unable to constrain estimates to be positive
  b <- -1
  
  ## Generalize ####
  dat$random <- dat[[random]]
  dat$y <- dat[[y]]
  dat$x <- dat[[x]]
  
  fit_method <- "NLME"
  
  ## Try NLME fit ####
  fit <- tryCatch(
    {fit <- nlme::nlme(y~Aopt-(b*(x-Topt)^2),fixed=list(Aopt + Topt + b ~ 1), random = Aopt+Topt ~ 1 | random,
                       start=list(fixed = c(Aopt = max(dat$y), Topt = 25, b = 0.05)),
                       data=dat,
                       # Stable:
                       # control = list(msMaxIter=1000))
                       
                       # In development: Worked on 2021-06-23, Problem: parameter cannot be constrained
                       # control = nlmeControl(opt = "nlminb", maxiter=1000, lower = c(-Inf, -Inf, 0.000001)))
                       
                       # Working update from 2021-06-23:
                       control = nlmeControl(opt = "nlminb", maxiter=10000))
    },
    
    warning = function(cond){
      # print("There was a warning")
      return(NA)
    },
    error = function(cond){
      # print("This message will not be printed.")
      return(NA)
    },
    finally = {
      #pass
    })
  
  ## Extract and return fit information
  if (length(fit) != 1) {
    aopt<-summary(fit)$tTable[[1]]
    topt<-summary(fit)$tTable[[2]]
    b<-summary(fit)$tTable[[3]]
    aopt.se<-summary(fit)$tTable[[4]]
    topt.se<-summary(fit)$tTable[[5]]
    b.se<-summary(fit)$tTable[[6]]
    
    #to get R2 between fitted and observed photosynthesis
    r<-cor(fitted(fit),dat$Photo)
    r2<-r*r
    
    #test for normality of residuals
    rest<-residuals(fit)
    norm<-shapiro.test(rest)
    s<-norm$statistic
    pvalue<-norm$p.value
    
    param<-cbind(aopt,topt,b,aopt.se,topt.se,b.se)
    
    names(param)[1:6] <- c("aopt","topt","b","aopt.se","topt.se","b.se")
    
    #return(param)
    #return(c(param,r2,s,pvalue))
    
    out <- tibble(aopt = aopt,
                  aopt.se = aopt.se,
                  topt = topt,
                  topt.se = topt.se,
                  b = b,
                  b.se = b.se,
                  r2 = r2,
                  s = s,
                  pvalue = summary(fit)$tTable[[14]],
                  fit_method = fit_method)
  }
  
  ## Try NLS fit, if NLME failed ####
  if (length(fit) == 1 | b < 0) {
    fit_method <- "NLS"
    
    fit <- tryCatch(
      {fit <- stats::nls(Photo ~ Aopt - (b * (x - Topt)^2),
                         start = list(Aopt = max(dat$Photo), Topt = 25, b = 0.05),
                         data = dat,
                         # Stable:
                         # control = list(msMaxIter=1000))
                         
                         # In development: Worked on 2021-06-23
                         control = list(maxiter = 10000),
                         lower = c(-Inf, -Inf, 0.00001),
                         algorithm = "port")
      },
      
      warning = function(cond){
        # print("There was a warning")
        return(NA)
      },
      error = function(cond){
        # print("This message will not be printed.")
        return(NA)
      },
      finally = {
        #pass
      })
    
    ## Extract and return fit information
    if (fit_method == "NLS" & length(fit) != 1) {
      # In dev.:
      results <- summary(fit)$coefficients
      
      
      # Old and working:
      # A.1<-summary(fit)
      # results<-(A.1$coefficients[1:6])
      # names(results)[1:6] <- c("aopt","topt","b","aopt.se","topt.se","b.se")
      
      # Catching non-convergence errors that give NA for the fit's s.e.
      if ((is.nan(results[5]))) {
        fit_method <- "NA"
        out <- tibble(aopt = NA,
                      aopt.se = NA,
                      topt = NA,
                      topt.se = NA,
                      b = NA,
                      b.se = NA,
                      r2 = NA,
                      s = NA,
                      pvalue = NA,
                      fit_method = fit_method)
        return(out)
      }
      
      #TT.i <- seq(min(dat$x),max(dat$x),length=51)
      #predicts <- predictNLS(fit, newdata=data.frame(x = TT.i),interval="confidence",level=0.95)
      #predicts.df <- data.frame(predicts$summary)
      #predicts.df$x <- TT.i
      
      r<-cor(fitted(fit), dat$Photo)
      r2<-r*r
      
      #test for normality of residuals
      rest<-residuals(fit)
      norm<-shapiro.test(rest)
      s<-norm$statistic
      pvalue<-norm$p.value
      
      out <- tibble(aopt = results[1],
                    aopt.se = results[4],
                    topt = results[2],
                    topt.se = results[5],
                    b = results[3],
                    b.se = results[6],
                    r2 = r2,
                    s = s,
                    pvalue = results[11],
                    fit_method = fit_method)
    }
  }
  
  ## Return NA, if all fits failed ####
  if (fit_method == "NLS" & length(fit) == 1) {
    
    fit_method <- "NA"
    out <- tibble(aopt = NA,
                  aopt.se = NA,
                  topt = NA,
                  topt.se = NA,
                  b = NA,
                  b.se = NA,
                  r2 = NA,
                  s = NA,
                  pvalue = NA,
                  fit_method = fit_method)
  }
  return(out)
}
