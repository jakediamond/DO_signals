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
  gpp_choice = 3, # none, constant, ramp, or sine
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
# 6 treatments of varying GPP frequency and phase, last is normal
trts <- tibble(
  gpp_choice = c(rep(3, 6), 1), #sine then constant
  f = c(1, 1, 1, 2, 4, 8, 0),
  phase = c(0, 0.25, -0.25, 0, 0, 0, 0)
)

# Apply 1D DO model to the treatments
mods <- trts |>
  mutate(out = future_pmap(list(
    f = f,
    phase = phase),
    func_mod))

# Create dataframe for with input data for different reaches
res <- mods |>
  mutate(data = future_map(out, df_fun))

# Clean the results up for plotting
df_res <- select(res, -out) |>
  tidytable::mutate(data = tidytable::map(data, ~ .x |> 
                                            mutate_all(as.character))) |> 
  tidytable::unnest(data) |>
  type_convert() |>
  mutate(line = if_else(phase != 0,
                        paste0("freq ", f, ", ", "phase ", phase),
                        paste("freq", f)),
         line = if_else(f == 0, "uniform", line))


# Plot over time
# Get legend info for plotting lines
line_info <- tibble(line = unique(df_res$line),
                    color = c(rep("darkred", 3), "darkblue", "blue", "lightblue", "black"),
                    linetype = c(1, 2, 3, rep(1,4)))

# Look at variation in GPP at noon
p_gpp <- ggplot(data = filter(df_res, type == "GPP", time == 36,
                     line != "uniform")|>
         mutate(fac = if_else(f == 1, "b", "a")),
       aes(x = dist,
           y = value,
           color = line,
           linetype = line)) +
  geom_line() +
  scale_color_manual(name = "",
                     values = line_info$color,
                     labels = line_info$line) +
  scale_linetype_manual(name = "",
                        values = line_info$linetype,
                        labels = line_info$line) +
  facet_grid(rows = vars(fac)) +
  theme_classic() +
  labs(x = "distance (m)",
       y = expression(GPP~"("*g~O[2]~m^{-2}~d^{-1}*")")) +
  theme(legend.position = "none",
        strip.text = element_blank())

# Plot over time
p_t <- ggplot() +
  geom_line(data = filter(df_res, type == "DO", reach == 95),
            aes(x = time, y = value, 
                color = line, 
                linetype = line,
                group = line),
            size = 1) +
  scale_color_manual(name = "",
                     values = line_info$color,
                     labels = line_info$line) +
  scale_linetype_manual(name = "",
                        values = line_info$linetype,
                        labels = line_info$line) +
  scale_x_continuous(breaks = seq(0, 48, 24)) +
  theme_classic(base_size = 10) +
  labs(x = "time (hr)",
       y = expression(O[2]~"("*mg~L^{-1}*")")) +
  theme(legend.position = c(0.15, 0.8),
        legend.background = element_rect(color = "transparent",
                                         fill = "transparent"))

# Plot the data over space, can pick multiple times
p_d <- ggplot(data = filter(df_res, time %in% c(24, 28, 32, 36, 40, 44),
                     type == "DO", line != "uniform"),
       aes(x = dist,
           y = as.numeric(value),
           color = as.factor(paste0(as.numeric(time %% 24), ":00")),
           group = time)) +
  geom_line() +
  facet_grid(rows = vars(line)) +
  scale_color_viridis_d(name = "") +
  # scale_y_continuous()
  theme_classic() +
  xlab("distance (m)") +
  ylab(expression(O[2]~"("*mg~L^{-1}*")")) +
  theme(legend.position = c(0.5, 0.15),
        legend.background = element_rect(fill = "transparent"),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.2, "cm"))
  # guides(colour = guide_legend(nrow = 1))

layout <- c(
  area(t = 1, l = 1, b = 2, r = 4),
  area(t = 3, l = 1, b = 7, r = 4),
  area(t = 1, l = 5, b = 7, r = 6)
)

p <- p_gpp + p_t + p_d + plot_layout(design = layout)
p

ggsave(plot = p,
       filename = file.path("results", "replicate_fig4.png"),
       dpi = 300,
       units = "cm",
       height = 16,
       width = 18.4)
