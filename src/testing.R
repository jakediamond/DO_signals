del_t    <- 1/60             # time step (h)
days     <- 2               # number of days to simulate
sim_time <- (24 * days ) - 1 # simulation time (h)
times    <- seq(0, sim_time, # sequence of in-model times to simulate
                by = del_t)  


parms_pf <- c(
  # Site info
  doy = 180, # day of year, used for modeling light
  latitude = 29, # latitude, used for modeling light
  temp = 20, # temperature of water, degrees C
  p_atm = 1, # atmospheric pressure (atm)
  
  # Grid parameters and boundary conditions
  dx = 800, # segment length (m)
  L = 80000, # reach length (m)
  C_up = NA_real_, # upstream concentration (mg/L); for Rainbow river = 5.5
  C_down = NA_real_, # downstream concentration (mg/L), if NA_real_, then set to saturation
  
  # Hydraulic parameters
  Q = 21.6, # discharge (m3/s); for Rainbow river choose 20.4 for mean u=0.18 m/s
  d = 2, # average depth (m)
  w = 60, # average channel width (m)
  D = 0 * 3600, # dispersivity (m2/h); 12240 for Rainbow River (Hensley and Cohen 2012)
  
  # Reactive parameters
  gpp_mean = 0, # mean daily GPP rate (g O2/m2/d)
  er_mean = 0, # mean daily ER rate (g O2/m2/d)
  K = 0, # gas exchange coefficient (1/d)
  
  # Spatial variability parameters
  ramp_rate = 1, # 100%rate GPP ramp; +1 = increase from 0-200%; -1 = decrease from 200-0% (GPP%/box) 
  f = 1, # frequency of longitudinal GPP sine wave, 1 = 1 wave per reach (1/reach_length)
  phase = 0, #phase of longitudinal GPP sine wave (-)
  
  # Storage parameters
  A_stor_frac = 0.2, # fraction of channel area that is storage area (-); for Rainbow River (Hensley and Cohen 2012)
  alpha = 0 # exchange coefficient (1/h); 0.252 for Rainbow River (Hensley and Cohen 2012)
)



# Try with events ---------------------------------------------------------
DOpulse <- data.frame(var = 1,
                      time = seq(0,15, by = del_t),
                      value = 15,
                      method = "replace")
out_event_pf <- ode.1D(y = yini,
              times = times,
              func = model,
              parms = parms_pf,
              nspec = 2,
              dimens = with(as.list(parms_pf), L / dx),
              gpp_choice = gpp_choice,
              events = list(data = DOpulse))

df_pulse_pf <- df_fun(out_event_pf, pars = parms_pf)
# distinct(df_pulse_pf, reach)

ggplot(data = filter(df_pulse_pf, reach %in% c(2, 5, 25, 45, 49),
                     type == "DO"),
       aes(x = time,
           y = value,
           group = reach,
           color = dist)) +
  geom_line() +
  # geom_point(alpha = 0.4) +
  # facet_grid(rows = vars(type_plot), scales = "free_y",
  #            labeller = label_parsed) +
  scale_color_viridis_c(name = "Distance (m)") +
  scale_x_continuous(breaks = seq(0, sim_time, 24)) +
  theme_bw() +
  xlab("time (h)") +
  ylab("")

locs <- tibble(time = round(c(648, 2000, 5000, 10000, 15000, 20000) / 648),
            type = "DO") |>
  mutate(reach = round(time * 648/ with(as.list(parms_pf), dx))) #c(6, 10, 25, 50, 75, 99),

event_d <- semi_join(mutate(df_pulse, time = as.numeric(time),
                           reach = as.numeric(reach)),
                    locs) |>
  mutate(drop = (value- DO_ini) / DO_ini)

event_pf <- semi_join(mutate(df_pulse_pf, time = as.numeric(time),
                              reach = as.numeric(reach)),
                       locs) |>
  mutate(drop = (value- DO_ini) / DO_ini)
