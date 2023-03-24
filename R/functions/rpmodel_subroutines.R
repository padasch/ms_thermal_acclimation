# Sub Routines of P-Model ----
## Disclaimer ----
#' This model code was partially taken from a developmental version of the P-Model as
#' accessible under https://github.com/geco-bern/rpmodel. See that repository
#' for the latest version. 
#' For detailed description, see: Stocker et al. 2020 https://doi.org/10/gjjxmd

# FvCB-related Functions ----
# .............................................................................
calc_ci <- function(ca, gammastar, xi, vpd, patm, ci_fixed_at_ppm = NA){
    
    if (is.numeric(ci_fixed_at_ppm)) {
        ci <- co2_to_ca(ci_fixed_at_ppm, patm)
        
    } else if (is.na(ci_fixed_at_ppm) | ci_fixed_at_ppm == "prentice14") {
        ci <- ca - (ca - gammastar) / (1 + xi / sqrt(vpd))
    }
    
    return(ci)
}

# .............................................................................
calc_aj <- function(kphio, ppfd, jmax, gammastar, ci = NA, ca, fapar, theta = 0.85, j_method = "smith37", model = "analytical", gs = NA) {
    
    ## Get iabs ....................................................................................
    iabs <- ppfd * fapar
    
    ## Get J based on Jmax-Limitation formulation ..................................................
    ## smith37:    Jmax Limitation following Smith et al. (1937):
    if (j_method %in% c("smith37", "wang17")) {j <- (4 * kphio * iabs) / sqrt(1 + ((4 * kphio * iabs)/(jmax))^2)}
    
    ## farquhar89: Jmax Limitation following Farquhar et al. (1989):
    if (j_method %in% c("farquhar89", "smith19")) {j <- (kphio * iabs + jmax - sqrt(( kphio * iabs + jmax)^2 - (4*kphio*theta*iabs*jmax))) / (2*theta)}
    
    
    ## Find ci if not given
    if (model == "numerical" | is.na(ci)) {
        ## Solve Eq. system:
        ## A = gs (ca - ci)
        ## A = j/4 * (ci-gammastar)/ci+2*gammastar)
        ## This leads to a quadratic equation:
        ## A * ci^2 + B * ci + C  = 0
        ## 0 = a + b*x + c*x^2
        ## with
        A  <- -gs
        B  <- gs * ca - 2 * gammastar * gs - j/4
        C  <- 2 * gammastar * gs * ca + gammastar * j/4
        ci <- QUADM(A, B, C)
    }
    
    aj   <- j/4 * (ci - gammastar)/(ci + 2 * gammastar)
    
    out = list(aj   = aj,
               j    =  j,
               ci   = ci)
    
    return(out)
    
    # Test Case:
    calc_aj(kphio = 0.1, ppfd = 130, jmax=8, gs = 0.5,
            gammastar=4, ci = NA, ca=40, fapar=1, theta = 0.85, 
            j_method = "smith37", model = "analytical")
}

# .............................................................................
calc_ac <- function(ci = NA, ca, gammastar, kmm, vcmax, model = "analytical", gs) {
    
    ## Find ci if not given: 
    if (model == "numerical" | is.na(ci)) {
        ## Solve Eq. system:
        ## A = gs (ca- ci)
        ## A = Vcmax * (ci - gammastar)/(ci + Kmm)
        ## This leads to a quadratic equation:
        ## A * ci^2 + B * ci + C  = 0
        ## 0 = a + b*x + c*x^2
        ## with
        A  <- -1.0 * gs
        B  <- gs * ca - gs * kmm - vcmax
        C  <- gs * ca * kmm + vcmax * gammastar
        ci <- QUADM(A, B, C)
    }
    
    ac  <- vcmax * (ci - gammastar)/(ci + kmm)
    out <- list(ac = ac, ci = ci)
    
    return(out)
    
    
}

# .............................................................................
calc_ac <- function(ci, ca, gammastar, kmm, vcmax, model = "analytical", gs) {
    
    if (model == "analytical") {
        ac <- vcmax * (ci - gammastar)/(ci + kmm)
        
        out <- list(ac = ac)
        
    } else if (model == "numerical") {
        ## Rubisco is limiting
        ## Solve Eq. system
        ## A = gs (ca- ci)
        ## A = Vcmax * (ci - gammastar)/(ci + Kmm)
        
        ## This leads to a quadratic equation:
        ## A * ci^2 + B * ci + C  = 0
        ## 0 = a + b*x + c*x^2
        
        ## with
        A <- -1.0 * gs
        B <- gs * ca - gs * kmm - vcmax
        C <- gs * ca * kmm + vcmax * gammastar
        
        ci_c <- QUADM(A, B, C)
        ac <- vcmax * (ci_c - gammastar) / (ci_c + kmm)
        
        out <- list(ac = ac,
                    ci = ci_c)
    }
    
    return(out)
}

# .............................................................................
calc_rd <- function(tc_leaf, rd25, tc_growth, q10 = 2, method_rd_scale = "heskel16") {
    
    ## Get temperature scaling for Rd:
    if (method_rd_scale == "heskel16") {
        apar <- 0.1012
        bpar <- 0.0005
        f <- exp( apar * (tc_leaf - 25.0) - bpar * (tc_leaf^2 - 25.0^2) )
    }
    
    if (method_rd_scale == "arrhenius") {
        f <- calc_ftemp_arrh(tc_leaf + 273.15, dha = 20700) # dha: Activation energy taken from Kumarathunge et al. (2019), Table 1, Mature Natural Environment
    }
    
    if (method_rd_scale == "q10") {
        f <- (3.22 - 0.046 * tc_leaf)^(tc_leaf - 25.0)/10 # Q10 changes with tc_leaf acc. to Tjoelker et al. (2001)
    }
    
    rd <- rd25 * f
    
    return(rd)
}

# .............................................................................
calc_rd25 <- function(vcmax25, method_rd25 = "atkin15"){
  ## Get base rate of Rd at 25 degC:
  if (method_rd25 == "atkin15") {
    rd_to_vcmax <- 0.015 # Ratio of Rdark to Vcmax25, Atkin et al., 2015 for C3 herbaceous
  }
  
  if (method_rd25 == "kumarathunge19") {
    rd_to_vcmax <- 0.0360 - 0.0010 * tc_growth # Acclimated rd_to_vcmax taken from Kumarathunge et al. (2019), Table 1, Mature Natural Environment
  }
  
  rd25 <- vcmax25 * rd_to_vcmax
  
  return(rd25)
}

# .............................................................................
calc_jmax <- function(kphio, iabs, ci, gammastar, method, theta = NA, c = NA){ 
    
    if (method == "smith37" | method == "wang17") {
        c      <- 0.103 # Estimated by Wang et al. (2017)
        c_star <- 4*c
        jmax   <- 4*kphio*iabs / (sqrt(1 / (1 - (c_star*(ci + 2*gammastar)/(ci-gammastar))^(2/3)) - 1))  
    }
    
    if (method == "farquhar89" | method == "smith19") {
        theta <- 0.85
        c     <- 0.053 # Estimated by Smith et al. (2019)
        m     <- (ci - gammastar) / (ci + 2 * gammastar)
        omega <- -(1 - 2 * theta) + sqrt((1-theta) * (1 / (4*c/m * (1-theta*4*c/m)) - 4*theta))
        jmax  <- kphio * iabs * omega
    }
    
    return(jmax)
}


# .............................................................................
calc_ftemp_arrh <- function(
    tk,
    dha,
    tkref = 298.15
){
    
    # Note that the following forms are equivalent:
    # ftemp = exp( dha * (tk - 298.15) / (298.15 * kR * tk) )
    # ftemp = exp( dha * (tc - 25.0)/(298.15 * kR * (tc + 273.15)) )
    # ftemp = exp( (dha/kR) * (1/298.15 - 1/tk) )
    
    kR   <- 8.3145     # Universal gas constant, J/mol/K
    ftemp <- exp( dha * (tk - tkref) / (tkref * kR * tk) )
    
    return(ftemp)
}

# .............................................................................

calc_ftemp_kphio <- function( tc, c4 = FALSE ){
    
    if (c4){
        ftemp = -0.008 + 0.00375 * tc - 0.58e-4 * tc^2   # Based on calibrated values by Shirley
    } else {
        ftemp <- 0.352 + 0.022 * tc - 3.4e-4 * tc^2
    }
    
    ## avoid negative values
    ftemp <- ifelse(ftemp < 0.0, 0.0, ftemp)
    
    return(ftemp)
}

# .............................................................................

calc_ftemp_inst_vcmax <- function( tcleaf, tcgrowth = NA, tcref = 25.0, method_ftemp = "kumarathunge19" ){
    
  # message("Using Method ", method_ftemp)
  
  # local parameters
  Rgas   <- 8.3145 # universal gas constant (J/mol/K)
  tkref  <- tcref + 273.15  # to Kelvin
  tkleaf <- tcleaf + 273.15  # conversion of temperature to Kelvin, tcleaf is the instantaneous leaf temperature in degrees C.
  
  if (method_ftemp %in% c("kattge07", "kumarathunge19")) {
  
    # Kattge2007 Parametrization
    Hd     <- 200000 # deactivation energy (J/mol)
    Ha    <- 71513  # activation energy (J/mol)
    a_ent <- 668.39 # offset of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K)
    b_ent <- 1.07   # slope of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K^2)
    
    if (method_ftemp == "kumarathunge19"){
        # Kumarathunge2019 Implementation:
        # local parameters
        a_ent = 645.13 # offset of entropy vs. temperature relationship (J/mol/K)
        b_ent = 0.38   # slope  of entropy vs. temperature relationship (J/mol/K^2)
        
        # local variables
        Ha = 42600 + (1140 * tcgrowth) # Acclimation for vcmax
    }
    
    # calculate entropy following Kattge & Knorr (2007), negative slope and y-axis intersect is when expressed as a function of temperature in degrees Celsius, not Kelvin !!!
    dent <- a_ent - (b_ent * tcgrowth)  # 'tcgrowth' corresponds to 'tmean' in Nicks, 'tc25' is 'to' in Nick's
    
    fva <- calc_ftemp_arrh( tkleaf, Ha, tkref = tkref )
    fvb <- (1 + exp( (tkref * dent - Hd)/(Rgas * tkref) ) ) / (1 + exp( (tkleaf * dent - Hd)/(Rgas * tkleaf) ) )
    fv  <- fva * fvb
  } else if (method_ftemp == "leuning02") {
    # Ref: Leuning, R. (2002). Temperature dependence of two parameters in a photosynthesis model. Plant, Cell & Environment, 25(9), 1205–1210. https://doi.org/10.1046/j.1365-3040.2002.00898.x
    # Table 2:
    Ha <- 73637
    Hd <- 149252
    Sv <- 486
    
    term_1 <- 1 + exp( (Sv * tkref  - Hd) / (Rgas * tkref) ) 
    term_3 <- 1 + exp( (Sv * tkleaf - Hd) / (Rgas * tkleaf) )
    term_2 <- exp((Ha / (Rgas * tkref)) * (1 - tkref/tkleaf)) # Careful: In Eq. (1) in Leuning et al. (1992), there is a bracket missing in this term!

    fv <- term_1 * term_2 / term_3
  } else {
    stop("Invalid method_ftemp.")
  }
    
    return( fv )
}

# .............................................................................

calc_ftemp_inst_jmax <- function( tcleaf, tcgrowth, tchome = NA, tcref = 25.0, method_ftemp = "kumarathunge19" ){
    
  # local parameters
  Rgas  <- 8.3145 # universal gas constant (J/mol/K)
  tkref <- tcref + 273.15  # to Kelvin
  tkleaf <- tcleaf + 273.15  # conversion of temperature to Kelvin, tcleaf is the instantaneous leaf temperature in degrees C.
  
  if (method_ftemp %in% c("kattge07", "kumarathunge19")) {
      
    # Kattge2007 Parametrization
    Hd    <- 200000 # deactivation energy (J/mol)
    Ha    <- 49884  # activation energy (J/mol)
    a_ent <- 659.70 # offset of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K)
    b_ent <- 0.75   # slope of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K^2)
    
    # calculate entropy following Kattge & Knorr (2007), negative slope and y-axis intersect is when expressed as a function of temperature in degrees Celsius, not Kelvin !!!
    dent <- a_ent - b_ent * tcgrowth   # 'tcgrowth' corresponds to 'tmean' in Nicks, 'tc25' is 'to' in Nick's
    
    if(method_ftemp == "kumarathunge19"){
        # Kumarathunge2019 Implementation:
        # local parameters
        Ha    = 40710  # activation energy (J/mol)
        a_ent = 658.77 # offset of entropy vs. temperature relationship (J/mol/K)
        b_ent = 0.84   # slope of entropy vs. temperature relationship (J/mol/K^2)
        c_ent = 0.52   # 2nd slope of entropy vs. temperature (J/mol/K^2)
        
        # Entropy calculation, equations given in Celsius, not in Kelvin
        dent = a_ent - (b_ent * tchome) - c_ent * (tcgrowth - tchome)
    }
    
    fva <- calc_ftemp_arrh( tkleaf, Ha, tkref = tkref )
    fvb <- (1 + exp( (tkref * dent - Hd)/(Rgas * tkref) ) ) / (1 + exp( (tkleaf * dent - Hd)/(Rgas * tkleaf) ) )
    fv  <- fva * fvb
    
  } else if (method_ftemp == "leuning02") {
    # Ref: Leuning, R. (2002). Temperature dependence of two parameters in a photosynthesis model. Plant, Cell & Environment, 25(9), 1205–1210. https://doi.org/10.1046/j.1365-3040.2002.00898.x
    # Table 2:
    Ha <- 50300
    Hd <- 152044
    Sv <- 495
    
    term_1 <- 1 + exp( (Sv * tkref  - Hd) / (Rgas * tkref) ) 
    term_3 <- 1 + exp( (Sv * tkleaf - Hd) / (Rgas * tkleaf) )
    term_2 <- exp((Ha / (Rgas * tkref)) * (1 - tkref/tkleaf)) # Careful: In Eq. (1) in Leuning et al. (1992), there is a bracket missing in this term!
    
    fv <- term_1 * term_2 / term_3
  } else {
    stop("Invalid method_ftemp.")
  }
    
    return( fv )
}

calc_gammastar <- function( tc, patm ) {
    
    # (J/mol) Activation energy, Bernacchi et al. (2001)
    dha    <- 37830
    
    # Pa, value based on Bernacchi et al. (2001), 
    # converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.
    gs25_0 <- 4.332
    
    gammastar <- gs25_0 * patm / calc_patm(0.0) * calc_ftemp_arrh( (tc + 273.15), dha=dha )
    
    return( gammastar )
}

calc_kmm <- function( tc, patm ) {
    
    dhac   <- 79430      # (J/mol) Activation energy, Bernacchi et al. (2001)
    dhao   <- 36380      # (J/mol) Activation energy, Bernacchi et al. (2001)
    kco    <- 2.09476e5  # (ppm) O2 partial pressure, Standard Atmosphere
    
    ## k25 parameters are not dependent on atmospheric pressure
    kc25 <- 39.97   # Pa, value based on Bernacchi et al. (2001), converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.
    ko25 <- 27480   # Pa, value based on Bernacchi et al. (2001), converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.
    
    ## conversion to Kelvin
    tk <- tc + 273.15
    
    kc <- kc25 * calc_ftemp_arrh( tk, dha=dhac )
    ko <- ko25 * calc_ftemp_arrh( tk, dha=dhao )
    
    po  <- kco * (1e-6) * patm         # O2 partial pressure
    kmm <- kc * (1.0 + po/ko)
    
    return(kmm)
}

density_h2o <- function(tc, p){
  
  # Calculate lambda, (bar cm^3)/g:
  my_lambda <- 1788.316 +
    21.55053*tc +
    -0.4695911*tc*tc +
    (3.096363e-3)*tc*tc*tc +
    -(7.341182e-6)*tc*tc*tc*tc
  
  # Calculate po, bar
  po <- 5918.499 +
    58.05267*tc +
    -1.1253317*tc*tc +
    (6.6123869e-3)*tc*tc*tc +
    -(1.4661625e-5)*tc*tc*tc*tc
  
  # Calculate vinf, cm^3/g
  vinf <- 0.6980547 +
    -(7.435626e-4)*tc +
    (3.704258e-5)*tc*tc +
    -(6.315724e-7)*tc*tc*tc +
    (9.829576e-9)*tc*tc*tc*tc +
    -(1.197269e-10)*tc*tc*tc*tc*tc +
    (1.005461e-12)*tc*tc*tc*tc*tc*tc +
    -(5.437898e-15)*tc*tc*tc*tc*tc*tc*tc +
    (1.69946e-17)*tc*tc*tc*tc*tc*tc*tc*tc +
    -(2.295063e-20)*tc*tc*tc*tc*tc*tc*tc*tc*tc
  
  # Convert pressure to bars (1 bar <- 100000 Pa)
  pbar <- (1e-5)*p
  
  # Calculate the specific volume (cm^3 g^-1):
  v <- vinf + my_lambda/(po + pbar)
  
  # Convert to density (g cm^-3) -> 1000 g/kg; 1000000 cm^3/m^3 -> kg/m^3:
  rho <- (1e3/v)
  
  return(rho)
}

calc_viscosity_h2o <- function(tc, p) {
    
    # Define reference temperature, density, and pressure values:
    tk_ast  <- 647.096    # Kelvin
    rho_ast <- 322.0      # kg/m^3
    mu_ast  <- 1e-6       # Pa s
    
    # Get the density of water, kg/m^3
    rho <- density_h2o(tc, p)
    
    # Calculate dimensionless parameters:
    tbar  <- (tc + 273.15)/tk_ast
    tbarx <- tbar^(0.5)
    tbar2 <- tbar^2
    tbar3 <- tbar^3
    rbar  <- rho/rho_ast
    
    # Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
    mu0 <- 1.67752 + 2.20462/tbar + 0.6366564/tbar2 - 0.241605/tbar3
    mu0 <- 1e2*tbarx/mu0
    
    # Create Table 3, Huber et al. (2009):
    h_array <- array(0.0, dim=c(7,6))
    h_array[1,] <- c(0.520094, 0.0850895, -1.08374, -0.289555, 0.0, 0.0)  # hj0
    h_array[2,] <- c(0.222531, 0.999115, 1.88797, 1.26613, 0.0, 0.120573) # hj1
    h_array[3,] <- c(-0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.0) # hj2
    h_array[4,] <- c(0.161913,  0.257399, 0.0, 0.0, 0.0, 0.0) # hj3
    h_array[5,] <- c(-0.0325372, 0.0, 0.0, 0.0698452, 0.0, 0.0) # hj4
    h_array[6,] <- c(0.0, 0.0, 0.0, 0.0, 0.00872102, 0.0) # hj5
    h_array[7,] <- c(0.0, 0.0, 0.0, -0.00435673, 0.0, -0.000593264) # hj6
    
    # Calculate mu1 (Eq. 12 & Table 3, Huber et al., 2009):
    mu1 <- 0.0
    ctbar <- (1.0/tbar) - 1.0
    # print(paste("ctbar",ctbar))
    # for i in xrange(6):
    for (i in 1:6){
        coef1 <- ctbar^(i-1)
        # print(paste("i, coef1", i, coef1))
        coef2 <- 0.0
        for (j in 1:7){
            coef2 <- coef2 + h_array[j,i] * (rbar - 1.0)^(j-1)
        }
        mu1 <- mu1 + coef1 * coef2
    }
    mu1 <- exp( rbar * mu1 )
    # print(paste("mu1",mu1))
    
    # Calculate mu_bar (Eq. 2, Huber et al., 2009)
    #   assumes mu2 = 1
    mu_bar <- mu0 * mu1
    
    # Calculate mu (Eq. 1, Huber et al., 2009)
    mu <- mu_bar * mu_ast    # Pa s
    
    return( mu )
}

calc_ftemp_inst_rd <- function( tc ){
    
    # local parameters
    apar <- 0.1012
    bpar <- 0.0005
    
    fr <- exp( apar * (tc - 25.0) - bpar * (tc^2 - 25.0^2) )
    
    return(fr)
}

calc_patm <- function(elv, patm0 = 101325) 
{
  kTo <- 298.15
  kL <- 0.0065
  kG <- 9.80665
  kR <- 8.3145
  kMa <- 0.028963
  patm <- patm0 * (1 - kL * elv/kTo)^(kG * kMa/(kR * kL))
  return(patm)
}

co2_to_ca <- function( co2, patm ){
  # Input:    - float, annual atm. CO2, ppm (co2)
  #           - float, monthly atm. pressure, Pa (patm)
  # Output:   - ca in units of Pa
  # Features: Converts ca (ambient CO2) from ppm to Pa.
  ca   <- ( 1.0e-6 ) * co2 * patm         # Pa, atms. CO2
  return( ca )
}

# P-Model related functions ----
calc_optimal_chi <- function(kmm, gammastar, ns_star, ca, vpd, beta ){
    
    # Input:    - float, 'kmm' : Pa, Michaelis-Menten coeff.
    #           - float, 'ns_star'  : (unitless) viscosity correction factor for water
    #           - float, 'vpd' : Pa, vapor pressure deficit
    # Output:   float, ratio of ci/ca (chi)
    # Features: Returns an estimate of leaf internal to ambient CO2
    #           partial pressure following the "simple formulation".
    # Depends:  - kc
    #           - ns
    #           - vpd
    
    ## Avoid negative VPD (dew conditions), resolves issue #2 (https://github.com/stineb/rpmodel/issues/2)
    vpd <- ifelse(vpd < 0, 0, vpd)
    
    ## leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio
    xi  <- sqrt( (beta * ( kmm + gammastar ) ) / ( 1.6 * ns_star ) )
    chi <- gammastar / ca + ( 1.0 - gammastar / ca ) * xi / ( xi + sqrt(vpd) )
    
    ## more sensible to use chi for calculating mj - equivalent to code below
    # # Define variable substitutes:
    # vdcg <- ca - gammastar
    # vacg <- ca + 2.0 * gammastar
    # vbkg <- beta * (kmm + gammastar)
    #
    # # Check for negatives, vectorized
    # mj <- ifelse(ns_star>0 & vpd>0 & vbkg>0,
    #              mj(ns_star, vpd, vacg, vbkg, vdcg, gammastar),
    #              rep(NA, max(length(vpd), length(ca)))
    #              )
    
    ## alternative variables
    gamma <- gammastar / ca
    kappa <- kmm / ca
    
    ## use chi for calculating mj
    mj <- (chi - gamma) / (chi + 2 * gamma)
    
    ## mc
    mc <- (chi - gamma) / (chi + kappa)
    
    ## mj:mv
    mjoc <- (chi + kappa) / (chi + 2 * gamma)
    
    # format output list
    out <- list(
        xi = xi,
        chi = chi,
        mc = mc,
        mj = mj,
        mjoc = mjoc
    )
    return(out)
}


calc_mprime <- function( mc ){
    # Input:  mc   (unitless): factor determining LUE
    # Output: mpi (unitless): modified m accounting for the co-limitation
    #                         hypothesis after Prentice et al. (2014)
    
    kc <- 0.41          # Jmax cost coefficient
    
    mpi <- mc^2 - kc^(2.0/3.0) * (mc^(4.0/3.0))
    
    # Check for negatives:
    mpi <- ifelse(mpi>0, sqrt(mpi), NA)
    
    return(mpi)
}


calc_lue_vcmax_wang17 <- function(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress){
    
    ## Include effect of Jmax limitation
    len <- length(out_optchi[[1]])
    mprime <- calc_mprime( out_optchi$mj )
    
    out <- list(
        
        mprime = mprime,
        
        ## Light use efficiency (gpp per unit absorbed light)
        lue = kphio * ftemp_kphio * mprime * c_molmass * soilmstress,
        
        ## Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
        vcmax_unitiabs = kphio * ftemp_kphio * out_optchi$mjoc * mprime / out_optchi$mj * soilmstress,
        
        ## complement for non-smith19
        omega      = rep(NA, len),
        omega_star = rep(NA, len)
        
    )
    
    return(out)
}



calc_lue_vcmax_smith19 <- function(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress, theta = 0.85){
    
    len <- length(out_optchi[[1]])
    
    # Adopted from Nick Smith's code:
    # Calculate omega, see Smith et al., 2019 Ecology Letters
    omega <- function( theta, c_cost, m ){
        
        cm <- 4 * c_cost / m                        # simplification term for omega calculation
        v  <- 1/(cm * (1 - theta * cm)) - 4 * theta # simplification term for omega calculation
        
        # account for non-linearities at low m values
        capP <- (((1/1.4) - 0.7)^2 / (1-theta)) + 3.4
        aquad <- -1
        bquad <- capP
        cquad <- -(capP * theta)
        m_star <- (4 * c_cost) / polyroot(c(aquad, bquad, cquad))
        
        omega <- ifelse(  m < Re(m_star[1]),
                          -( 1 - (2 * theta) ) - sqrt( (1 - theta) * v),
                          -( 1 - (2 * theta))  + sqrt( (1 - theta) * v)
        )
        return(omega)
    }
    
    ## constants
    c_cost <- 0.05336251
    
    ## factors derived as in Smith et al., 2019
    omega <- omega( theta = theta, c_cost = c_cost, m = out_optchi$mj )          # Eq. S4
    omega_star <- 1.0 + omega - sqrt( (1.0 + omega)^2 - (4.0 * theta * omega) )       # Eq. 18
    
    ## Effect of Jmax limitation
    mprime <- out_optchi$mj * omega_star / (8.0 * theta)
    
    ## Light use efficiency (gpp per unit absorbed light)
    lue <- kphio * ftemp_kphio * mprime * c_molmass * soilmstress
    
    # calculate Vcmax per unit aborbed light
    vcmax_unitiabs  <- kphio * ftemp_kphio * out_optchi$mjoc * omega_star / (8.0 * theta) * soilmstress   # Eq. 19
    
    out <- list(
        lue            = lue,
        vcmax_unitiabs = vcmax_unitiabs,
        omega          = omega,
        omega_star     = omega_star
    )
    
    return(out)
}



calc_lue_vcmax_none <- function(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress){
    ## Do not include effect of Jmax limitation
    len <- length(out_optchi[[1]])
    
    out <- list(
        
        ## Light use efficiency (gpp per unit absorbed light)
        lue = kphio * ftemp_kphio * out_optchi$mj * c_molmass * soilmstress,
        
        ## Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
        vcmax_unitiabs = kphio * ftemp_kphio * out_optchi$mjoc * soilmstress,
        
        ## complement for non-smith19
        omega               = rep(NA, len),
        omega_star          = rep(NA, len)
    )
    
    return(out)
}



calc_lue_vcmax_c4 <- function( kphio, ftemp_kphio, c_molmass, soilmstress ){
    
    len <- length(kphio)
    out <- list(
        ## Light use efficiency (gpp per unit absorbed light)
        lue = kphio * ftemp_kphio * c_molmass * soilmstress,
        
        ## Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
        vcmax_unitiabs = kphio * ftemp_kphio * soilmstress,
        
        ## complement for non-smith19
        omega               = rep(NA, len),
        omega_star          = rep(NA, len)
    )
    
    return(out)
}


calc_chi_c4 <- function(){
    
    # (Dummy-) ci:ca for C4 photosynthesis
    out <- list( chi=9999, mc=1, mj=1, mjoc=1 )
    return(out)
}

# Additional Functions ----
VPDairToLeaf <- function(VPD, Tair, Tleaf, Pa=101){
  
  # This function was adopted from the R package 'plantecophys'  
  # Duursma (2015) https://doi.org/10/bkmj.
  
  e <- esat(Tair, Pa) - VPD*1000
  vpd <- esat(Tleaf, Pa) - e
  
  return(vpd/1000)
}

VPDtoRH <- function(VPD, TdegC, Pa=101){
  # This function was adoted from the R package 'plantecophys'  
  # Duursma (2015) https://doi.org/10/bkmj.
  
  ## VPD and Pa in kPa
  esatval <- esat(TdegC, Pa)
  e <- pmax(0, esatval - VPD*1000)
  RH <- 100 * e/esatval
  return(RH)
}

# Calculation of vapor pressure in Pa
esat <- function(TdegC, Pa=101){  
  
  # This function was adopted from the R package 'plantecophys'  
  # Duursma (2015) https://doi.org/10/bkmj.
  
  ## Pa in kPa
  a <- 611.21
  b <- 17.502
  c <- 240.97
  f <- 1.0007 + 3.46 * 10^-8 * Pa * 1000
  esatval <- f * a * (exp(b * TdegC/(c + TdegC)))
  return(esatval)
}
