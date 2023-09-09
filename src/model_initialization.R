#
# Authors: Jake Diamond, Robert Hensley, Matt Cohen
# Purpose: To define model paramters and initial conditions
# Date: November 19, 2021
# 

# Get simulation times ----------------------------------------------------
del_t    <- 0.5             # time step (h)
days     <- 5               # number of days to simulate
sim_time <- (24 * days ) - 1 # simulation time (h)
times    <- seq(0, sim_time, # sequence of in-model times to simulate
                by = del_t)  

# Default model parameters --------------------------------------------------------
parms <- c(
  # Site info
  doy = 180, # day of year, used for modeling light
  latitude = 29, # latitude, used for modeling light
  temp = 20, # temperature of water, degrees C
  p_atm = 1, # atmospheric pressure (atm)
  
  # Grid parameters and boundary conditions
  dx = 200, # segment length (m)
  L = 20000, # reach length (m)
  C_up = NA_real_, # upstream concentration (mg/L); for Rainbow river = 5.5
  C_down = NA_real_, # downstream concentration (mg/L), if NA_real_, then set to saturation
  
  # Hydraulic parameters
  Q = 21.6, # discharge (m3/s); for Rainbow river choose 20.4 for mean u=0.18 m/s
  d = 2, # average depth (m)
  w = 60, # average channel width (m)
  D = 4 * 3600, # dispersivity (m2/h); 12240 for Rainbow River (Hensley and Cohen 2012)
  
  # Reactive parameters
  gpp_mean = 12, # mean daily GPP rate (g O2/m2/d)
  er_mean = 12, # mean daily ER rate (g O2/m2/d)
  K = 2, # gas exchange coefficient (1/d)
  
  # Spatial variability parameters
  gpp_choice = 1, # none, constant, ramp, or sine
  ramp_rate = 1, # 100%rate GPP ramp; +1 = increase from 0-200%; -1 = decrease from 200-0% (GPP%/box) 
  f = 1, # frequency of longitudinal GPP sine wave, 1 = 1 wave per reach (1/reach_length)
  phase = 0, #phase of longitudinal GPP sine wave (-)
  
  # Storage parameters
  A_stor_frac = 0.2, # fraction of channel area that is storage area (-); for Rainbow River (Hensley and Cohen 2012)
  alpha = 0.252 # exchange coefficient (1/h); 0.252 for Rainbow River (Hensley and Cohen 2012)
)

# Initial conditions ------------------------------------------------------
# Two vectors of stream and transient storage DO concentrations, start at sat.
DO_ini <- O2_sat(with(as.list(parms), temp)) # mg/L
DO_stor_ini <- DO_ini # mg/L
yini <- c(DO = rep(DO_ini, with(as.list(parms), L / dx)),
          DO_stor = rep(DO_stor_ini, with(as.list(parms), L / dx)))