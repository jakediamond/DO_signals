# -------------------------------------
# Author: Jake Diamond
# Purpose: Functions for use in model
# Date: 2023-04-10
# -------------------------------------

# Helper functions --------------------------------------------------------
# O2 saturation function
O2_sat <- function (temp = 25, press = 1) {
  # Conversions
  mlL_mgL   <- 1.42905 # O2 mL/L to mg/L, per USGS memo 2011.03
  mmHg_atm  <- 760 # mmHg to atm
  
  # Calculate saturated concentrations of O2 (mg/L) using Garcia-Benson
  # vapor pressure of water (mmHg)
  u <- 10^(8.10765 - 1750.286 / (235 + temp)) 
  #press correction (=1 at 1atm) USGS memos 81.11 and 81.15
  press_corr <- (press * mmHg_atm - u) / (760 - u) 
  Ts <- log((298.15 - temp)/(273.15 + temp)) # scaled temperature
  lnC <- 2.00907 + 3.22014 * Ts + 4.0501 * Ts^2 + 4.94457 * 
    Ts^3 + -0.256847 * Ts^4 + 3.88767 * Ts^5
  o2.sat <- exp(lnC) #O2 saturation (mL/L)
  o2_sat <- o2.sat * mlL_mgL * press_corr # (mg/L)
  return(o2_sat)
}

# Photosynthetically active radiation function (umol/m2/s)
par_fun <- function (hour, day = 180, latitude = 0, max.insolation = 2326) {
  declin <- (23.45 * sin((2*pi / 365) * (284 + day))) * pi / 180
  hour.angle <- (360/24) * (hour - 12) * pi / 180
  lat <- latitude * pi / 180
  zenith <- acos(sin(lat) * sin(declin) + 
                   cos(lat) * cos(declin) * cos(hour.angle))
  insolation <- max.insolation * cos(zenith)
  insolation <- pmax(insolation, 0)
}

# Set up model grid 
grid_fun <- function(L, dx) {
  # Grid of reaches with total length L and reach length L / dx
  grid <- setup.grid.1D(L = L, N = L / dx)
}

dgrid_fun <- function(D, grid) {
  # Dispersion grid
  grid <- setup.prop.1D(value = D, # (m2/h)
                               grid = grid)
}

vgrid_fun <- function(Q, w, d, grid) {
  # Velocity grid, all the same based on discharge and area (assumed rectangular)
  grid <- setup.prop.1D(value = Q * 3600 / (w * d), # (m/h)
                               grid = grid)
}
