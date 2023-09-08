#
# Authors: Jake Diamond, Robert Hensley, Matt Cohen
# Purpose: To compare 1D-modeled dissolved oxygen with different parameter sets
# Date: November 19, 2021
# 

# Load libraries
library(furrr) # for running in parallel
library(tidyverse)

# Set up parallelization
plan(multisession, workers = 4)

# Load the model
source(file.path("src", "1D_DO_model.R"))
source(file.path("src", "functions_for_model.R"))
# Load the functions
source(file.path("src", "model_comparison_functions.R"))

# Get simulation times ----------------------------------------------------
del_t    <- 0.1             # time step (h)
days     <- 2               # number of days to simulate
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
  gpp_choice = 1, #constant
  gpp_mean = 12, # mean daily GPP rate (g O2/m2/d)
  er_mean = 12, # mean daily ER rate (g O2/m2/d)
  K = 2, # gas exchange coefficient (1/d)
  
  # Spatial variability parameters
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

# Now get the different treatments ----------------------------------------
# Create dataframe of treatments
# 12 treatments of varying K, D, and storage area, final is plug flow
trts <- tibble(
  K = c(1, 5, 20, rep(2, 9)),
  D = c(4, 4, 4, 2, 8, 32, 64, 4, 4, 4, 4, 0) * 3600,
  A_stor_frac = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.01, 0.1, 0.5, 2, 0.2),
  alpha = c(rep(0.252, 11), 0)
  )

# Apply 1D DO model to the treatments
mods <- trts |>
  mutate(out = future_pmap(list(
    K = K,
    D = D,
    A_stor_frac = A_stor_frac),
    func_mod))

# Create dataframe for with input data for different reaches
res <- mods |>
  mutate(data = future_map(out, df_fun))

df_res <- select(res, -out) |>
  tidytable::mutate(data = tidytable::map(data, ~ .x |> 
                      mutate_all(as.character))) |> 
  tidytable::unnest(data) |>
  filter(reach == 50) |>
  type_convert() |>
  filter(type == "DO")

# changes in K
p_k <- ggplot() +
  geom_line(data = filter(df_res, D == 4*3600, A_stor_frac == 0.2),
            aes(x = time, y = value, color = as.factor(K), group = K),
            size = 1.5) +
  scale_color_viridis_d(expression(K~"("*d^{-1}*")")) +
  geom_line(data = filter(df_res, D == 0),
            aes(x = time, y = value), color = "black",
            size = 1) +
  scale_x_continuous(breaks = seq(0, 48, 24)) +
  theme_classic(base_size = 10) +
  labs(x = "hr",
       y = expression(O[2]~"("*mg~L^{-1}*")")) +
  theme(legend.position = c(0.12, 0.7),
        legend.key.size = unit(0.25, "cm"),
        legend.background = element_rect(fill = "transparent"))
p_k

# changes in D
p_d <- ggplot() +
  geom_line(data = filter(df_res, K == 2, A_stor_frac == 0.2, D!=0) |>
              mutate(D = D / 3600),
            aes(x = time, y = value, color = as.factor(D), group = D),
            size = 1.5) +
  scale_color_viridis_d(expression(D~"("*m^{2}~s^{-1}*")")) +
  geom_line(data = filter(df_res, D == 0),
            aes(x = time, y = value), color = "black",
            size = 1) +
  scale_x_continuous(breaks = seq(0, 48, 24)) +
  theme_classic(base_size = 10) +
  labs(x = "hr",
       y = expression(O[2]~"("*mg~L^{-1}*")")) +
  theme(legend.position = c(0.12, 0.7),
        legend.key.size = unit(0.25, "cm"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_rect(fill = "transparent"))
p_d

# Changes in storage
p_s <- ggplot() +
  geom_line(data = filter(df_res, K == 2, D == 4*3600),
            aes(x = time, y = value, color = as.factor(A_stor_frac), group = A_stor_frac),
            size = 1.5) +
  scale_color_viridis_d(expression(A[S]~"(-)")) +
  geom_line(data = filter(df_res, D == 0),
            aes(x = time, y = value), color = "black",
            size = 1) +
  scale_x_continuous(breaks = seq(0, 48, 24)) +
  theme_classic(base_size = 10) +
  labs(x = "hr",
       y = expression(O[2]~"("*mg~L^{-1}*")")) +
  theme(legend.position = c(0.12, 0.7),
        legend.key.size = unit(0.25, "cm"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_rect(fill = "transparent"))
p_s

p <- p_k | p_d | p_s
p

ggsave(plot = p,
       filename = file.path("results", "replicate_fig2.png"),
       dpi = 300,
       units = "cm",
       height = 9,
       width = 18.4)
