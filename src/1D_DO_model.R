#
# Authors: Jake Diamond, Robert Hensley, Matt Cohen
# Purpose: To model dissolved oxygen spatiotemporal evolution in streams/rivers
# Date: October 5, 2019
# 

# Load libraries
library(ReacTran)

# Overall model function ----------------------------------------------------------
model <- function(time, state, parms){
  with(as.list(parms), {
    
    ############################ Set up 1D grid and states ####################
    grid <- grid_fun(L, dx) # framework
    D.grid <- dgrid_fun(D, grid) # dispersion
    v.grid <- vgrid_fun(Q, w, d, grid) # velocity
    
    # Total number of boxes
    N <- L / dx
    
    # Unpack states
    O2 <- state[1:N]
    O2stor <- state[(N+1):(2*N)]

    ############################ Internal parameters ###########################
    # Calculate parameters in correct units (time unit is hour)
    k <- K * d / 24 # gas exchange velocity (m/h)
    A <- w * d # rectangular river cross-sectional area (m2)
    A_s <- A * A_stor_frac # transient storage area (m2)
    
    # mean PAR at the time step (umol/m2/s)
    par <- par_fun(time, day = doy, latitude = latitude)
    
    # mean daily PAR
    meanpar <- mean(par_fun(0:23, day = doy, latitude = latitude))
    
    # Calculate saturation value for O2
    O2_sat <- O2_sat(temp, p_atm)
    
    # Boundary conditions for O2 concentration
    Cdown <- ifelse(is.na(C_down), O2_sat, C_down)
    Cup   <- ifelse(is.na(C_up),   O2_sat, C_up)
    
    ############################ Calculate fluxes ###########################
    # GPP is switching function depending on user choice (g O2/m2/h)
    # First calculate GPP as linear function of mean PAR 
    gpp_base <- (gpp_mean / 24) * (par / meanpar)
    
    # Constant hourly ecosystem respiration (mol O2/m2/hr)
    er_const <- (er_mean / 24)
    
    # Gas exchange (mol/m2/h) from river perspective (+ is into river)
    reaeration <- k * (O2_sat - O2)
    
    ############################ Spatial variability ##########################
    gpp_choice_chr <- case_when(
      gpp_choice == 0 ~ "none",
      gpp_choice == 1 ~ "constant",
      gpp_choice == 2 ~ "ramp",
      gpp_choice == 3 ~ "sine"
    )
    
    # GPP varies spatially according to user choice
    gpp <- switch(gpp_choice_chr,
                  none = rep(0, N),
                  constant = rep(gpp_base, N),
                  ramp = seq(gpp_base * (1 - ramp_rate),
                             gpp_base * (1 + ramp_rate),
                             length.out = N),
                  sine = gpp_base * sin((2 * pi * f) * 
                                          (seq(0, 1, length.out = N) + phase)
                                        ) + gpp_base)

    # advection-dispersion in each box (g O2/m3/h)
    adv_dis <- tran.1D(C = O2,
                       C.up = Cup, # upstream boundary concentration
                       C.down = Cdown, # downstream boundary concentration
                       D = D.grid, 
                       v = v.grid, 
                       dx = grid)$dC 
    
    # Ecosystem respiration in the grid array (g O2/m2/h)
    er <- rep(er_const, N)
    
    # Change in storage concentration exchange with river (g O2/m3/h)
    storage <- alpha * (A / A_s) * (O2 - O2stor)
    
    # Change in stream O2 concentration due to transfer to transient storage (g O2/m3/h)
    storage_str <- alpha * (O2stor - O2)
    
    ################### Difference equations ################################
    # Rates of change (g/m3/hr)
    dO2 <- (adv_dis + storage_str + (gpp - er + reaeration) / d)
    dO2stor <- storage
    
    ################### Model output ################################
    diffs <- c(dO2 = dO2, dO2stor = dO2stor)
    
    # Fluxes
    fluxes <- c(GPP = gpp, ER = er)
    # Other outputs
    O2sat  <- O2 / O2_sat * 100
    
    # Output of model
    return(list(diffs, fluxes, PAR = rep(par, N)))
  })
}  # end of model
