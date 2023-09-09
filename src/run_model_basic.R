#
# Authors: Jake Diamond
# Purpose: Run the spatiotemporal metabolism model with O2 output
# Date: 2023 February 17
# 

# Load libraries
library(tidyverse)
library(deSolve)

# Get model functions -----------------------------------------------------
source(file.path("src", "functions_for_model.R"))
source(file.path("src", "1D_DO_model.R"))
source(file.path("src", "model_treatment_functions.R"))

# Get simulation times ----------------------------------------------------
del_t    <- 0.1             # time step (h)
days     <- 2               # number of days to simulate
sim_time <- (24 * days ) - 1 # simulation time (h)
times    <- seq(0, sim_time, # sequence of in-model times to simulate
                by = del_t)  

# Define model parameters --------------------------------------------------------
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
  Q = 21.6, # discharge (m3/s); for Rainbow river choose 21.6 for mean u=0.18m3/s
  d = 2, # average depth (m)
  w = 60, # average channel width (m)
  D = 4 * 3600, # dispersivity (m2/h)
  
  # Reactive parameters
  gpp_choice = 3, # 0,1,2, or 3 = none, constant, ramp, or sine
  gpp_mean = 12, # mean daily GPP rate (g O2/m2/d)
  er_mean = 12, # mean daily ER rate (g O2/m2/d)
  K = 2, # gas exchange coefficient (1/d)
  
  # Spatial variability parameters
  ramp_rate = 1, # 100%rate GPP ramp; +1 = increase from 0-200%; -1 = decrease from 200-0% (GPP%/box) 
  f = 1, # frequency of longitudinal GPP sine wave, 1 = 1 wave per reach (1/reach_length)
  phase = 0, #phase of longitudinal GPP sine wave (-)
  
  # Storage parameters
  A_stor_frac = 2, # fraction of channel area that is storage area (-); for Rainbow River (Hensley and Cohen 2012)
  alpha = 0.252 # exchange coefficient (1/h); 0.252 for Rainbow River (Hensley and Cohen 2012)
)

# Initial conditions ------------------------------------------------------
# Two vectors of stream and transient storage DO concentrations, start at sat.
DO_ini <- O2_sat(with(as.list(parms), temp)) # mg/L
DOstor_ini <- DO_ini # mg/L
yini <- c(DO = rep(DO_ini, with(as.list(parms), L / dx)),
          DOstor = rep(DOstor_ini, with(as.list(parms), L / dx)))

# Run the model -----------------------------------------------------------
# uses the R deSolve function (lsoda method)
out <- ode.1D(y = yini,
              times = times,
              func = model,
              parms = parms,
              nspec = 2,
              dimens = with(as.list(parms), L / dx))

# Examine the model output ------------------------------------------------
# Reorganize the data into long-form
df <- df_fun(out)

# Plot the data over time, can pick multiple reaches
ggplot(data = filter(df, reach %in% c(50)),
       aes(x = time,
           y = value,
           color = dist,
           group = reach)) +
  geom_line() +
  facet_grid(rows = vars(type_plot), scales = "free_y",
             labeller = label_parsed) +
  scale_color_viridis_c(name = "Distance (m)") +
  scale_x_continuous(breaks = seq(0, sim_time, 24)) +
  theme_bw() +
  xlab("time (h)") +
  ylab("")

# Plot the data over space, can pick multiple times
ggplot(data = filter(df, time %in% c(24, 28, 32, 36, 40, 44),
                     type == "DO"),
       aes(x = dist,
           y = as.numeric(value),
           color = as.factor(as.numeric(time)))) +
  geom_line() +
  scale_color_viridis_d(name = "time") +
  theme_bw() +
  xlab("distance (m)") +
  ylab("")

