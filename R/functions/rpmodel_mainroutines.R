# Main Routines of P-Model ----
## Disclaimer ----
#' This model code was taken from a developmental version of the P-Model as
#' accessible under https://github.com/geco-bern/rpmodel. See that repository
#' for the latest version. 
#' For detailed description, see: Stocker et al. 2020 https://doi.org/10/gjjxmd 

# Acclimated Response --------------------------------------------------

rpmodel <- function(inputs){
  
    #__________________________________________________________________________________________________#
    # 0. Check input ----
    inp_vars <- names(inputs)
    
    ## Necessary variables that have to be provided
    nec_vars <- c("tc_air", "tc_home","ppfd", "vpd", "co2", "kphio_calib")
    
    ## Additional variables that can be provided
    add_vars <- c(
        "fapar", "patm", "tc_leaf", ## Environment
        "beta", "soilm", "meanalpha", "apar_soilm", "bpar_soilm", "c4", "soilmstress", ## P-Model Parameters
        "wind", "leaf_width", "stomatal_ratio", "leaf_abs") ## Energy Balance Parameters
    
    ## Settings that can be provided
    set_vars <- c("method_optim", "method_ci", "method_jmaxlim", "method_ftemp", "method_eb",
                  "num_jmax_via_ratio", "do_ftemp_kphio", "do_ftemp_theta", "do_soilmstress",
                  "returnvar", "verbose", "method_rd25")
    
    nec_missing <- setdiff(nec_vars, inp_vars)
    add_missing <- setdiff(add_vars, inp_vars)
    set_missing <- setdiff(set_vars, inp_vars)
    
    if (length(nec_missing) > 0) {
        stop("[!] Missing necessary variables: ", paste0(nec_missing, sep = ", "))
    }
    
    # if (length(add_missing) > 0) {
    #     stop("[!] Missing additional variables: ", paste0(add_missing, sep = ", "))
    # }
    # 
    # if (length(set_missing) > 0) {
    #     stop("[!] Missing setting variables: ", paste0(set_missing, sep = ", "))
    # }
    
    ## Definition of additional variables if not provided:
    if ("tc_leaf" %in% add_missing)            inputs$tc_leaf <- NA
    if ("fapar" %in% add_missing)              inputs$fapar <- 1
    if ("beta" %in% add_missing)               inputs$beta <- 146
    if ("soilm" %in% add_missing)              inputs$soilm <- 1
    if ("meanalpha" %in% add_missing)          inputs$meanalpha <- 1
    if ("apar_soilm" %in% add_missing)         inputs$apar_soilm <- 0.0
    if ("bpar_soilm" %in% add_missing)         inputs$bpar_soilm <- 0.73
    if ("c4" %in% add_missing)                 inputs$c4 <- FALSE
    if ("soilmstress" %in% add_missing)        inputs$soilmstress <- 1

    if ("method_optim" %in% set_missing)       inputs$method_optim <- "analytical"
    if ("method_ci" %in% set_missing)       inputs$method_ci <- "prentice14"
    if ("method_jmaxlim" %in% set_missing)     inputs$method_jmaxlim <- "wang17"
    if ("method_ftemp" %in% set_missing)       inputs$method_ftemp <- "kumarathunge19"
    if ("method_eb" %in% set_missing)          inputs$method_eb <- "off"
    if ("method_rd25" %in% set_missing)          inputs$method_rd25 <- "heskel16"
    if ("num_jmax_via_ratio" %in% set_missing) inputs$num_jmax_via_ratio <- TRUE
    if ("do_ftemp_kphio" %in% set_missing)     inputs$do_ftemp_kphio <- TRUE
    if ("do_ftemp_theta" %in% set_missing)     inputs$do_ftemp_theta <- FALSE
    if ("do_soilmstress" %in% set_missing)     inputs$do_soilmstress <- FALSE
    if ("verbose" %in% set_missing)            inputs$verbose <- FALSE
    
    if (inputs$method_eb == "plantecophys") {
        ## Standard values for energy balance are defined in that function itself
        if ("wind" %in% add_missing)           inputs["wind"]           <- NA
        if ("leaf_width" %in% add_missing)     inputs["leaf_width"]     <- NA
        if ("stomatal_ratio" %in% add_missing) inputs["stomatal_ratio"] <- NA
        if ("leaf_abs" %in% add_missing)       inputs["leaf_abs"]       <- NA
    }
    
    ## Check of local variables
    if (!("patm" %in% inp_vars)) {
        if ("elv" %in% inp_vars) {
            inputs$patm <- calc_patm(inputs$elv)
            if (inputs$verbose) warning("Atmospheric pressure (patm) not provided. Calculating it as a function of elevation (elv), assuming standard atmosphere (101325 Pa at sea level).")
        } else { 
            stop("Aborted. Provide either elevation (arugment elv) or atmospheric pressure (argument patm).")
        }
    }
    
    ## Definition of local variables
    # Reduce to necessary variables for higher computation first
    inputs <- inputs %>% dplyr::select(any_of(c(nec_vars, add_vars, set_vars)))
    
    # Absorbed Irradiance
    inputs$iabs <- inputs$fapar * inputs$ppfd
    
    # Make local variables
    for(i in 1:length(inputs)) assign(names(inputs)[i], inputs[[i]])
    
    # Constants
    c_molmass <- 12.0107  # molecular mass of carbon (g)
    kPo <- 101325.0       # standard atmosphere, Pa (Allen, 1973)
    kTo <- 25.0           # base temperature, deg C (Prentice, unpublished)
    rd_to_vcmax <- 0.015  # Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous 
    
    # Optimization flags
    opt_convergence <- NA
    eb_convergence  <- NA
    
    #__________________________________________________________________________________________________#
    # 1. Analytical Model ----
      
    ## No energy balance accounting, assume tc_leaf = tc_air
    if ("tc_leaf" %in% add_missing) {tc_leaf <- tc_air}
        
    ## Rescale VPD to tc_leaf 
    vpd <- VPDairToLeaf(vpd/1000, tc_air, tc_leaf, patm/1000)*1000
    
    ## 1.2 Photosynthesis Variables ----
    ftemp_kphio <- 1
    
    if (do_ftemp_kphio) ftemp_kphio <- calc_ftemp_kphio( tc_leaf, c4 ) 
    if (do_soilmstress) soilmstress <- soilmstress( soilm, meanalpha, apar_soilm, bpar_soilm )
    
    ca         <- co2_to_ca( co2, patm )
    gammastar  <- calc_gammastar( tc_leaf, patm )
    kmm        <- calc_kmm( tc_leaf, patm )
    ns_star    <- calc_viscosity_h2o( tc_leaf, patm ) / calc_viscosity_h2o( kTo, kPo )  # (unitless)
    xi         <- sqrt( (beta * ( kmm + gammastar ) ) / ( 1.6 * ns_star ) )
    kphio_accl <- kphio_calib * ftemp_kphio
    iabs       <- fapar * ppfd
    
    if (do_ftemp_theta) {
        theta  <- 0.344375 + 0.012134 * tc_leaf
        inputs$theta <- theta
    } else {
        theta <- 0.85
        inputs$theta <- theta
    }
        
      ## The heart of the P-model: calculate ci:ca ratio (chi) and additional terms
      ## 1.3 Get chi, vcmax, jmax ----
      
      if (c4) {
        out_optchi <- calc_chi_c4()
      } else if (method_ci == "prentice14") {
        out_optchi <- calc_optimal_chi(kmm, gammastar, ns_star, ca, vpd, beta)
      } else {
        stop("rpmodel(): argument method_ci not idetified.")
      }
      
      ### Get vcmax and jmax
      if (c4) {
        out_lue_vcmax <- calc_lue_vcmax_c4(
          kphio_calib,
          ftemp_kphio,
          c_molmass,
          soilmstress
        )
      } else if (method_jmaxlim == "wang17" | method_jmaxlim == "smith37") {
      
        ## apply correction by Jmax limitation
        out_lue_vcmax <- calc_lue_vcmax_wang17(
          out_optchi,
          kphio_calib,
          ftemp_kphio,
          c_molmass,
          soilmstress
        )
      } else if (method_jmaxlim == "smith19" | method_jmaxlim == "farquhar89") {
        out_lue_vcmax <- calc_lue_vcmax_smith19(
          out_optchi = out_optchi,
          kphio = kphio_calib,
          ftemp_kphio = ftemp_kphio,
          c_molmass = c_molmass,
          soilmstress = soilmstress,
          theta = theta
        )
      } else if (method_jmaxlim == "none") {
        out_lue_vcmax <- calc_lue_vcmax_none(
          out_optchi,
          kphio_calib,
          ftemp_kphio,
          c_molmass,
          soilmstress
        )
      } else {
        stop("rpmodel(): argument method_jmaxlim not idetified.")
      }
      
      ## Corrolaries
      chi        <- out_optchi$chi
      xi         <- out_optchi$xi
      ci         <- chi * ca # leaf-internal CO2 partial pressure (Pa)
      iwue       <- ( ca - ci ) / 1.6  # intrinsic water use efficiency (in Pa)
      
      ## 1.4 Get GPP ----         
      # Gross Primary Productivity
      gpp <- iabs * out_lue_vcmax$lue   # in g C m-2 s-1
      
      # Carboxylation rate
      ftemp_vcmax       <- calc_ftemp_inst_vcmax(tcleaf = tc_leaf, tcgrowth = tc_air, method_ftemp = method_ftemp)
      vcmax             <- iabs * out_lue_vcmax$vcmax_unitiabs
      vcmax25           <- vcmax / ftemp_vcmax
      ac                <- calc_ac(ci, ca, gammastar, kmm, vcmax, model = method_optim)$ac
      
      ## Electron transport rate
      jmax              <- calc_jmax(kphio_accl, iabs, ci, gammastar, method = method_jmaxlim)
      ftemp_jmax        <- calc_ftemp_inst_jmax(tcleaf = tc_leaf, tcgrowth = tc_air, tchome = tc_home, method_ftemp = method_ftemp)
      jmax25            <- jmax / ftemp_jmax
      aj                <- calc_aj(kphio_accl, ppfd, jmax, gammastar, ci, ca, fapar, j_method = method_jmaxlim, model = method_optim, theta = theta)$aj
      
      ## Dark respiration
      ftemp_inst_rd      <- calc_ftemp_inst_rd( tc_leaf )
      rd_unitiabs        <- rd_to_vcmax * (ftemp_inst_rd / ftemp_vcmax) * out_lue_vcmax$vcmax_unitiabs
      rd                 <- iabs * rd_unitiabs
      
      ## Gross assimilation
      a_gross <- ifelse(aj < ac , aj, ac)
      
      ## Average stomatal conductance
      gs_c <- a_gross / (ca - ci)

      ## Calculate Rd25
      rd25 <- calc_rd25(vcmax25, method_rd25)
    
    
      ## 3. Output definition ----
      out <- tibble(
          # gpp             = gpp,   # remove this again later
          ca              = ca,
          gammastar       = gammastar,
          # kmm             = kmm,
          ns_star         = ns_star,
          chi               = chi,
          xi                = xi,
          # mj              = out_optchi$mj,
          # mc              = out_optchi$mc,
          ci                = ci,
          # iwue              = iwue,
          gs_c              = gs_c,
          vcmax             = vcmax,
          vcmax25           = vcmax25,
          jmax              = jmax,
          jmax25            = jmax25,
          kphio             = kphio_accl,
          tc_leaf           = tc_leaf,
          rd25              = rd25,
          # rd              = rd,
          opt_convergence   = opt_convergence,
          eb_convergence    = eb_convergence,
          a_gross           = a_gross,
          # inputs            = inputs
      )
    
  return( out )
}

#__________________________________________________________________________________________________#
# Instantaneous Response ----
rpmodel_inst <- function(inputs, make_checks = T) {
  
    ## 0. Check Input ----
    if (make_checks) {
        inp_vars <- names(inputs)
        
        ## Necessary variables that have to be provided
        nec_vars <- c("tc_growth_air", "tc_home", "kphio", "vcmax25", "jmax25", "rd25", "xi", ## Acclimated Conditions
                      "tc_air", "ppfd", "vpd", "co2") ## Instant Conditions
        
        ## Additional variables that can be provided
        add_vars <- c("fapar", "patm", "tc_leaf", ## Environment
                      "beta", "soilm", "meanalpha", "apar_soilm", "bpar_soilm", "c4", "soilmstress", "theta", ## P-Model Parameters
                      "wind", "leaf_width", "stomatal_ratio", "leaf_abs") ## Energy Balance Parameters
        
        ## Settings that can be provided
        set_vars <- c("method_optim", "method_ci", "method_jmaxlim", "method_ftemp", "method_eb",
                      "num_jmax_via_ratio", "do_ftemp_kphio", "do_ftemp_theta", "do_soilmstress",
                      "method_rd25", "method_rd", "q10", ## Instant Parameters
                      "returnvar", "verbose")
        
        nec_missing <- setdiff(nec_vars, inp_vars)
        add_missing <- setdiff(add_vars, inp_vars)
        set_missing <- setdiff(set_vars, inp_vars)
        
        if (length(nec_missing) > 0) {
            stop("[!] Missing necessary variables: ", paste0(nec_missing, sep = ", "))
        }
        
        ## Definition of global variables if not provided:
        if ("tc_leaf" %in% add_missing)            inputs$tc_leaf <- NA
        if ("fapar" %in% add_missing)              inputs$fapar <- 1
        if ("beta" %in% add_missing)               inputs$beta <- 146
        if ("theta" %in% add_missing)              inputs$theta <- 0.85
        if ("soilm" %in% add_missing)              inputs$soilm <- 1
        if ("meanalpha" %in% add_missing)          inputs$meanalpha <- 1
        if ("apar_soilm" %in% add_missing)         inputs$apar_soilm <- 0.0
        if ("bpar_soilm" %in% add_missing)         inputs$bpar_soilm <- 0.73
        if ("c4" %in% add_missing)                 inputs$c4 <- FALSE
        if ("soilmstress" %in% add_missing)        inputs$soilmstress <- 1
        
        if ("method_optim" %in% set_missing)       inputs$method_optim <- "analytical"
        if ("method_ci" %in% set_missing)          inputs$method_ci <- "prentice14"
        if ("method_jmaxlim" %in% set_missing)     inputs$method_jmaxlim <- "wang17"
        if ("method_ftemp" %in% set_missing)       inputs$method_ftemp <- "kumarathunge19"
        if ("method_eb" %in% set_missing)          inputs$method_eb <- "off"
        if ("num_jmax_via_ratio" %in% set_missing) inputs$num_jmax_via_ratio <- TRUE
        if ("do_ftemp_kphio" %in% set_missing)     inputs$do_ftemp_kphio <- TRUE
        if ("do_ftemp_theta" %in% set_missing)     inputs$do_ftemp_theta <- FALSE
        if ("do_soilmstress" %in% set_missing)     inputs$do_soilmstress <- FALSE
        if ("verbose" %in% set_missing)            inputs$verbose <- FALSE
        
        if ("method_rd" %in% set_missing)          inputs$method_rd <- "heskel16"
        if ("q10" %in% set_missing)                inputs$q10 <- 2
        
        if (inputs$method_eb == "plantecophys") {
            ## Standard values for energy balance are defined in that function itself
            if ("wind" %in% add_missing)           inputs["wind"]           <- NA
            if ("leaf_width" %in% add_missing)     inputs["leaf_width"]     <- NA
            if ("stomatal_ratio" %in% add_missing) inputs["stomatal_ratio"] <- NA
            if ("leaf_abs" %in% add_missing)       inputs["leaf_abs"]       <- NA
        }
        
        ## Check of local variables
        if (!("patm" %in% inp_vars)) {
            if ("elv" %in% inp_vars) {
                inputs$patm <- calc_patm(inputs$elv)
                if (inputs$verbose) warning("Atmospheric pressure (patm) not provided. Calculating it as a function of elevation (elv), assuming standard atmosphere (101325 Pa at sea level).")
            } else { 
                stop("Aborted. Provide either elevation (arugment elv) or atmospheric pressure (argument patm).")
            }
        }
        
        # if ( !(identical(NA, inputs$tc_leaf)) & inputs$method_eb != "off") {
        #     warning("Leaf temperature given in input but energy balance model called. Skipping energy balance!")
        # }
    }
    
    ## Definition of local variables
    # Reduce to necessary variables for higher computation first
    # message("Reducing")
    # inputs <- inputs %>% select(any_of(c(nec_vars, add_vars, set_vars)))
    
    # message("Making local vars")
    # Make local variables
    # for(i in 1:length(inputs)) assign(names(inputs)[i], inputs[[i]])
    
    # Constants
    c_molmass <- 12.0107  # molecular mass of carbon (g)
    kPo <- 101325.0       # standard atmosphere, Pa (Allen, 1973)
    kTo <- 25.0           # base temperature, deg C (Prentice, unpublished)
    rd_to_vcmax <- 0.015  # Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous 
    
    # Optimization flags
    opt_convergence <- NA
    eb_convergence  <- NA
    
    ## Prepare output
    agross <- NA
    anet <- NA
    vcmax <- NA
    jmax <- NA
    rd <- NA
    ci <- NA
    gs <- NA
    aj <- NA
    ac <- NA
    min_a <- NA
    gs_c <- NA
    eb_convergence <- NA
    dummy <- NA
    dummy_var <- NA

    # __________________________________________________________________________________________________#
    ## 1. Analytical Model ----
    ## Without energy balance accounting, assume tc_leaf = tc_air
    if ("tc_leaf" %in% add_missing) {inputs$tc_leaf <- inputs$tc_air}
  
    ## Photosynthetic parameters
    if (inputs$do_ftemp_theta) {
      theta  <- 0.344375 + 0.012134 * inputs$tc_growth_air
    } else {
      theta <- 0.85
    }
  
    ca <- co2_to_ca(inputs$co2, inputs$patm)
    gammastar <- calc_gammastar(inputs$tc_leaf, inputs$patm)
    kmm <- calc_kmm(inputs$tc_leaf, inputs$patm)
    ns_star <- calc_viscosity_h2o(inputs$tc_leaf, inputs$patm) / calc_viscosity_h2o(tc = 25.0, p = 101325.0)
    vcmax <- inputs$vcmax25 * calc_ftemp_inst_vcmax(inputs$tc_leaf, inputs$tc_growth_air, tcref = 25.0, method_ftemp = inputs$method_ftemp)
    jmax  <- inputs$jmax25  * calc_ftemp_inst_jmax(inputs$tc_leaf, inputs$tc_growth_air, inputs$tc_home, tcref = 25.0, method_ftemp = inputs$method_ftemp)
    
    # Debug
    # cat(" Tleaf: ", round(inputs$tc_leaf), " Tgrowth: ", inputs$tc_growth_air, "  Scalar: ", calc_ftemp_inst_vcmax(inputs$tc_leaf, inputs$tc_growth_air, method_ftemp = inputs$method_ftemp), " Vcmax25: ", inputs$vcmax25, " Vcmax: ",vcmax, "\n")
    
    # Instantaneous leaf-internal CO2: ci
    ci <- calc_ci(ca, gammastar, inputs$xi, inputs$vpd, inputs$patm, inputs$method_ci)
    
    # Electron Transport Rate: Aj
    aj <- calc_aj(inputs$kphio, inputs$ppfd, jmax, gammastar, ci, ca, inputs$fapar, theta, j_method = inputs$method_jmaxlim, model = "analytical", gs = NA)$aj
    
    # Carboxylation Rate: Ac
    ac <- calc_ac(ci = ci, ca, gammastar = gammastar, kmm = kmm, vcmax = vcmax, model = "analytical")$ac
    
    # Dark Respiration: Rd
    rd <- calc_rd(inputs$tc_leaf, inputs$rd25, inputs$tc_growth_air, q10 = inputs$q10, method_rd_scale = inputs$method_rd)
    
    # Gross assimilation rate: A_gross
    agross <- min(aj, ac)
    
    if (aj == ac) {
    min_a <- "colimit"
    } else if (agross == ac) {
    min_a <- "ac"
    } else if (agross == aj) {
    min_a <- "aj"
    }
    
    # Net assimilation rate: A_net
    anet <- agross - rd
    
    # Stomatal conductance: gs
    gs_c <- anet / (ca - ci)
    
    # __________________________________________________________________________________________________#
    ## 2. Definition of output ----
    out <- list(
      agross = agross,
      anet = anet,
      vcmax = vcmax,
      jmax = jmax,
      rd = rd,
      ci = ci,
      aj = aj,
      aj_net = aj - rd,
      ac = ac,
      ac_net = ac - rd,
      min_a = min_a,
      tc_leaf = inputs$tc_leaf,
      eb_convergence = eb_convergence,
      gs_c = gs_c,
      dummy = dummy_var
    )
    
  return(out)
}
