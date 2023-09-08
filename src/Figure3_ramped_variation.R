#
# Authors: Jake Diamond, Robert Hensley, Matt Cohen
# Purpose: To compare 1D-modeled dissolved oxygen with different parameter sets
# Date: November 19, 2021
# 

# Load libraries
library(furrr) # for running in parallel
library(tidyverse)
library(patchwork)

# Set up parallelization
plan(multisession, workers = 4)

# Load the model and necessary functions
source(file.path("src", "1D_DO_model.R"))
source(file.path("src", "functions_for_model.R"))
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
  gpp_mean = 12, # mean daily GPP rate (g O2/m2/d)
  er_mean = 12, # mean daily ER rate (g O2/m2/d)
  K = 2, # gas exchange coefficient (1/d)
  
  # Spatial variability parameters
  gpp_choice = 2, # none, constant, ramp, or sine
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
# 7 treatments of ramping of GPP
trts <- tibble(
  gpp_choice = 2, #ramped
  ramp_rate = c(-1, -2/3, -1/3, 0, 1/3, 2/3, 1)
)

# Apply 1D DO model to the treatments
mods <- trts |>
  mutate(out = future_pmap(list(
    ramp_rate = ramp_rate),
    func_mod))

# Create dataframe for with input data for different reaches
res <- mods |>
  mutate(data = future_map(out, df_fun))

# Clean up data
df_res <- select(res, -out) |>
  tidytable::mutate(data = tidytable::map(data, ~ .x |> 
                                            mutate_all(as.character))) |> 
  tidytable::unnest(data) |>
  type_convert()

# Plot over time
# Get legend info for plotting lines
line_info <- tibble(rate = unique(df_res$ramp_rate),
                    lab = c("2-0", "1.67-0.33", "1.33-0.67",
                            "none", "0.67-1.33", "0.33-1.67", "0-2"),
                    color = c("pink", "red", "darkred", "black", "darkblue",
                              "blue", "lightblue"))

# Look at variation in GPP at noon
p_gpp <- ggplot(data = filter(df_res, type == "GPP", time == 36),
                aes(x = dist,
                    y = value,
                    color = as.factor(ramp_rate))) +
  geom_line() +
  scale_color_manual(name = "",
                     values = line_info$color,
                     labels = line_info$lab) +
  theme_classic() +
  scale_x_continuous(breaks = c(0, 10000, 20000)) +
  labs(x = "distance (m)",
       y = expression(GPP~"("*g~O[2]~m^{-2}~d^{-1}*")")) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill='transparent', color=NA))
p_gpp

# Plot over time
p_t <- ggplot() +
  geom_line(data = filter(df_res, type == "DO", reach == 90),
            aes(x = time, y = value, 
                color = as.factor(ramp_rate)),
                # color = as.factor(round(ramp_rate,2)), group = ramp_rate),
            size = 1.5) +
  scale_color_manual(name = "ramp rate",
                     values = line_info$color,
                     labels = line_info$lab) +
  # scale_color_viridis_d(expression(K~"("*d^{-1}*")")) +
  scale_x_continuous(breaks = seq(0, 48, 24)) +
  scale_y_continuous(limits = c(6, 14)) +
  theme_classic(base_size = 10) +
  labs(x = "hr",
       y = expression(O[2]~"("*mg~L^{-1}*")")) +
  theme(legend.position = c(0.5, 0.1),
        legend.direction = "horizontal")

# Plot the data over space, can pick multiple times
p_d <- ggplot(data = filter(df_res, time %in% c(24, 28, 32, 36, 40, 44),
                     type == "DO",
                     ramp_rate %in% c(-1, 0, 1)),
       aes(x = dist,
           y = as.numeric(value),
           color = as.factor(paste0(as.numeric(time %% 24), ":00")),
           group = time)) +
  geom_line() +
  facet_grid(rows = vars(ramp_rate)) +
  scale_color_viridis_d(name = "time") +
  theme_classic() +
  xlab("distance (m)") +
  ylab(expression(O[2]~"("*mg~L^{-1}*")")) +
  theme(legend.position = c(0.45, 0.65),
        legend.background = element_rect(fill = "transparent"),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.2, "cm"))

p_d

layout <- c(
  area(t = 1, l = 1, b = 8, r = 12),
  area(t = 1, l = 6, b = 2, r = 9),
  area(t = 1, l = 13, b = 8, r = 18)
)

p <- p_t + p_gpp + p_d + plot_layout(design = layout)
p

ggsave(plot = p,
       filename = file.path("results", "replicate_fig3.png"),
       dpi = 300,
       units = "cm",
       height = 16,
       width = 18.4)

