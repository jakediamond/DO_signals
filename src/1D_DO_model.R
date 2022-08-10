#
# Authors: Jake Diamond, Robert Hensley, Matt Cohen
# Purpose: To model dissolved oxygen spatiotemporal evolution in streams/rivers
# Date: October 5, 2019
# 

# Load libraries
library(deSolve)
library(tidyverse)
library(ReacTran)
library(lubridate)

# Define model parameters --------------------------------------------------------
parms <- c(
  # Grid and time parameters
  dx = 200, # segment length (m)
  L = 20000, # reach length (m)
  C_down = 9, # downstream boundary concentration (mg/L)
  C_up = 9, # upstream boundary concentration (mg/L); for Rainbow river = 5.5
  
  # Hydraulic parameters
  Q = 21.6, # discharge (m3/s); for Rainbow river choose 21.6 for mean u=0.18m3/s
  d = 2, # average depth (m)
  w = 60, # average channel width (m)
  S = 0.001, # average channel slope (m/m)
  D = 3600, # dispersivity (m/h)
  
  # Reactive parameters
  GPP_max = 2, # maximum GPP rate (g O2/m2/h); multiply by 6 to estimate daily GPP
  k_gpp = 500, # PAR at which half of GPP_max is produced (umol/m2/s)
  do_sat = 9, # saturation for DO, should be T dependent (mg/L)
  K = 1.7, # gas exchange coefficient (1/d)
  er_const = -0.5, # constant ER rate (g O2/m2/h)
  ramp_rate = 2/3, # 100%rate GPP ramp; +1 = increase from 0-200%; -1 = decrease from 200-0% (GPP%/box) 
  f = 2, # frequency of longitudinal GPP sine wave, 1 = 1 wave per reach (1/reach_length)
  phase = 0, #phase of longitudinal GPP sine wave, fractions of pi (-)
  amp = 1, # amplitude of GPP sine wave, needs to be less than GPP_max to avoid negative GPP
  
  # Storage parameters
  A_stor_frac = 0.25, # fraction of channel area that is storage area (-); for Rainbow River (Hensley and Cohen 2012)
  stor_depth = 1, # arbitrary depth of storage to convert from area to volume (m)
  alpha = 0.252 # exchange coefficient (1/h); 0.252 for Rainbow River (Hensley and Cohen 2012)
 )

# Set up model grid -------------------------------------------------------
# Grid of reaches with total length L and reach length L / dx
grid <- with(as.list(parms),
             setup.grid.1D(L = L, N = L / dx)
             )
# Dispersion grid
D.grid <- with(as.list(parms),
                 setup.prop.1D(value = D, # (m2/h)
                               grid = grid)
               ) 
# Velocity grid, all the same based on discharge and area (assumed rectangular)
v.grid <- with(as.list(parms),
                 setup.prop.1D(value = Q * 3600 / (w * d), # (m/h)
                               grid = grid)
               )

# Overall model function ----------------------------------------------------------
model <- function(time, state, parms, gpp_choice){
  with(as.list(parms), {
    # Need to set the seed to avoid crazy computation issues, and to allow rep.
    set.seed(42)
    # Total number of boxes
    N = L / dx
    
    # Unpack states
    DO <- state[1:N]
    DO_stor <- state[(N+1):(2*N)]

    # Calculate parameters in correct units (time unit is hour)
    k = K * d / 24 # gas exchange velocity (m/h)
    A = w * d # rectangular river cross-sectional area (m2)
    A_s = A * A_stor_frac # transient storage area (m2)
    
    # Calculate in-stream reaction in each box (gpp and gas exchange)
    # Photosynthetically reactive radiation (umol/m2/s)
    par = ifelse(-500+2000*sin((2*pi*(time-(6/del_t)) / (24/del_t))) < 0,
                 0,
                 -500+2000*sin((2*pi*(time-(6/del_t)) / (24/del_t))))
    
    # GPP is switching function depending on user choice (g O2/m2/h)
    # First calculate base rate of GPP as Michaelis-Menten
    gpp_base = GPP_max * par / (k_gpp + par)

    # Then modulate this based on user choice
    gpp = switch(gpp_choice,
                 none = rep(0, N),
                 constant = rep(gpp_base, N),
                 sine = amp*gpp_base*sin(2*pi*f*seq(0,1,length.out=N)+phase*pi)+gpp_base,
                 ramp = seq(gpp_base*(1-ramp_rate),
                            gpp_base*(1+ramp_rate),
                            length.out = N))
                 
    # Reaeration in each box (g O2/m2/h)
    reaeration = k * (do_sat - DO)

    # advection-dispersion in each box (g O2/m3/h)
    adv_dis = tran.1D(C = DO,
                      C.up = C_up, # upstream boundary concentration
                      C.down = C_down, # downstream boundary concentration
                      D = D.grid, 
                      v = v.grid, 
                      dx = grid)$dC 
    
    # Ecosystem respiration amount (g O2/m2/h)
    er = rep(er_const, N)
    
    # Ecosystem respiration location changes where DO is consumed
    # Change in storage concentration due to respiration (g O2/m3/h)
    storage = alpha * (A / A_s) * (DO - DO_stor)
    
    # Change in stream DO concentration due to transfer to transient storage (g O2/m3/h)
    storage_str = alpha * (DO_stor - DO)
    
    # Rates of change (g O2/m3)
    dDO = (adv_dis + storage_str + (gpp + er + reaeration) / d) * del_t
    dDO_stor = storage * del_t
    diffs = c(dDO = dDO, dDO_stor = dDO_stor)
    
    # Fluxes
    fluxes = c(GPP = gpp * del_t, ER = er * del_t)
    
    # Output of model
    return(list(diffs, fluxes, PAR = rep(par, N)))
  })
}  # end of model

# Initial conditions ------------------------------------------------------
# Two vectors of stream and transient storage DO concentrations
DO_ini = 9 # mg/L
DO_stor_ini = 9 # mg/L
yini <- c(DO = rep(DO_ini, with(as.list(parms), L / dx)),
          DO_stor = rep(DO_stor_ini, with(as.list(parms), L / dx)))

# Run the model -----------------------------------------------------------
# Get simulation times
del_t <- 0.5 # time step (h)
days <- 5 # number of days to simulate
simulation_time <- 24 * days / del_t # simulation time (h)
times <- seq(0, simulation_time, by = del_t)

# GPP choice (choose between: "none", "constant", "sine,"ramp")
gpp_choice <- "sine"

# uses the R deSolve function (lsoda method)
out <- ode.1D(y = yini,
              times = times,
              func = model,
              parms = parms,
              nspec = 2,
              dimens = with(as.list(parms), L / dx),
              gpp_choice = gpp_choice)

# Examine the model output ------------------------------------------------
# Reorganize the data into long-form
df <- as_tibble(out) %>%
  tidyfast::dt_pivot_longer(cols = -time, names_to = "key") %>%
  tidyfast::dt_separate(key, c("type", "reach"), "(?<=[A-Za-z])(?=[0-9])", fixed = FALSE, perl = T) %>%
  mutate(time_hr = time * del_t,
         dist = as.numeric(reach) * with(as.list(parms), dx)) %>%
  mutate(type_plot = recode(type,
         `DO` = "DO~(mg~L^{-1})",
         `DO_stor` = "DO[storage]~(mg~L^{-1})",
         `ER.DO_stor` = "ER~(g~O^2~m^{-2}~h^{-1})",
         `ER` = "ER~(g~O^2~m^{-2}~h^{-1})",
         `GPP` = "GPP~(g~O^2~m^{-2}~h^{-1})",
         `PAR` = "PAR~({`mu`}*mol~m^{-2}~s^{-1})"))

# save base case data to compare later
# saveRDS(df, file.path("data", "DO_model_base_case.RDS"))

# Plot the data
ggplot(data = filter(df, reach %in% c(98)),
       aes(x = time_hr,
           y = value)) +
  geom_line() +
  # geom_point(alpha = 0.4) +
  facet_grid(rows = vars(type_plot), scales = "free_y",
             labeller = label_parsed) +
  scale_color_viridis_c(name = "Distance (m)") +
  scale_x_continuous(breaks = seq(0, simulation_time * del_t, 24)) +
  theme_bw() +
  xlab("Time (h)") +
  ylab("")

# Plot the data
ggplot(data = filter(df, time_hr %in% c(48, 52, 56, 60, 64, 68),
                     type == "DO"),
       aes(x = dist,
           y = as.numeric(value),
           color = as.factor(as.numeric(time_hr)))) +
  geom_line() +
  # geom_point(alpha = 0.4) +
  # facet_grid(rows = vars(type_plot), scales = "free_y",
  #            labeller = label_parsed) +
  scale_color_viridis_d(name = "time") +
  scale_y_continuous(limits = c(7, 13)) +
  # scale_x_continuous(breaks = seq(0, simulation_time * del_t, 24)) +
  theme_bw() +
  xlab("distance (m)") +
  ylab("")



ggplot(data = filter(df,
                     type == "GPP",
                     time_hr ==12),
       aes(x = dist,
           y = as.numeric(value))) +
  geom_line() +
  # geom_point(alpha = 0.4) +
  # facet_grid(rows = vars(type_plot), scales = "free_y",
  #            labeller = label_parsed) +
  # scale_color_viridis_d(name = "time") +
  # scale_y_continuous(limits = c(7, 13)) +
  # scale_x_continuous(breaks = seq(0, simulation_time * del_t, 24)) +
  theme_bw() +
  xlab("distance (m)") +
  ylab("")
