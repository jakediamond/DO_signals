#
# Authors: Jake Diamond, Robert Hensley, Matt Cohen
# Purpose: To compare 1D-modeled dissolved oxygen with estimates from streammetabolizer
# Date: November 19, 2021
# 

# Load libraries
# library(furrr)
library(streamMetabolizer)
library(tidyverse)

# Load the model
source(file.path("src", "1D_DO_model.R"))

# Load the functions
source(file.path("src", "model_streammetabolizer_comparison_functions.R"))

# Create dataframe of treatments
trts <- tibble(#A_stor_frac = c(0, 0.01, 0.5, 2),
               #K = c(0, 1, 5, 20), #prescribe the O2 gas exchange constant (1/d)
               #D = c(0, 8*3600, 32*3600, 64*3600) #prescribe dispersion (m2/h), delete line if you want it to be calculated internally
               # gpp_choice = "ramp",
               # GPP_max = runif(20, min = 0, max = 2),
               f = c(1, 2 ,4, 8)
               # phase = ,
               # ramp_rate = c(-1, -2/3, -1/3, 0, 1/3, 2/3, 1)
               ) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "treatment")

# Apply 1D DO model to the treatments
mods <- trts %>%
  mutate(out = map2(parameter, treatment, func_mod, gpp_choice = "sine"))

# Create dataframe for with input data for different reaches
res <- mods %>%
  mutate(data = map(out, df_fun)) %>%
  crossing(reach_no = c(98)) # right now, just choose one reach in the middle as representative, but can choose multiple

# Now add a column for the streammetabolizer input data
res <- res %>%
  mutate(sm_data = map2(data, reach_no, sm_fun))

# Get data in nice format for easy plotting
res_p <- select(res, -sm_data, -out, -reach_no) %>%
  mutate(data = map(data, ~ .x %>% 
                      mutate_all(as.character))) %>% 
  unnest(data) %>%
  filter(reach == 98) %>%
  type_convert()


df_fig2 <- filter(res_p, type == "DO")  %>%
  mutate(value = as.numeric(value)) %>%
  filter(between(time_hr, 24, 72))

p_fig2 <- ggplot(data = df_fig2,
                 aes(x = time_hr - 24, y = value, 
                     # color = (time_hr - 7)%%24,
                     linetype = as.factor(treatment),
                     group = treatment)) +
  geom_line(size = 1.5) +
  theme_bw(base_size = 10) +
  scale_x_continuous(breaks = seq(0,71,24), expand = c(0, 0)) +
  facet_wrap(~parameter) +
  # scale_color_viridis_c(guide = "none") +
  labs(x = "time (hour)", y = expression("dissolved oxygen (mg "*L^{-1}*")")) + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())
p_fig2

# Plot the data
res_p_d <- select(res, -sm_data, -out, -reach_no) %>%
  mutate(data = map(data, ~ .x %>% 
                      mutate_all(as.character))) %>% 
  unnest(data) %>%
  filter(type == "DO",
         time_hr %in% c(48, 52, 56, 60, 64, 68),
         treatment %in% c(-1, 0, 1)) %>%
  type_convert() %>%
  mutate(value = as.numeric(value))

p_fig3_dist <- ggplot(data = res_p_d,
       aes(x = dist,
           y = as.numeric(value),
           color = as.factor(as.numeric(time_hr %% 24)))) +
  geom_line(size = 1.2) +
  facet_grid(rows = vars(treatment)) +
  scale_color_viridis_d(name = "time") +
  scale_y_continuous(limits = c(7, 14)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "distance (m)", 
       y = expression("dissolved oxygen (mg "*L^{-1}*")"))
p_fig3_dist











# Get streammetabolizer estimates based on 1D DO model inputs
results_sm <- res %>%
  mutate(sm = map(sm_data, sm_res, 
                  modtype = "mle", # can prescribe mle or bayesian ("bayes")
                  gpp_fun = "satlight",  #can have gpp be a saturating function of light or linear ("linlight)
                  er_fun = "constant", #er is constant, can be a function of temp ("q10temp")
                  k = "none" #prescribe initial K600 value if you want (1/d), leave as 'none" you want the model to choose
                  )
         ) 

# Dataframe of comparisons between streammetabolizer and DO model
df_com <- results_sm %>%
  mutate(temp = map_dbl(sm_data, ~mean(.$temp.water)), #stream temperature from model
         Sc = 1801 - 120.1 * temp + 3.782 * temp^2 -0.0476 * temp^3, #schmidt number based on temp
         k600_mod = k_O2 /sqrt(600/Sc), #get k600 from DO model instead of KO2
         modinfo = map2(data, reach_no, do_mod_res)) %>%
  select(alpha, gpp_choice, er_loc_choice, k_O2, reach_no, D, sm, temp, k600_mod, modinfo) %>%
  unnest(cols = c(sm, modinfo))

# Save these results so you don't have to run again
# saveRDS(df_com, "C:/Users/jake.diamond/Dropbox/nonlinear_programming/mle_results_disp_alpha_kconst")

# Get into good format for plotting
df_com_l <- df_com %>%
  filter(!str_detect(msgs.fit, c("W", "E"))) %>% #no days with errors in streammetabolizer
  select(-GPP.lower, -GPP.upper, -ER.lower, -ER.upper, -date_mod, -msgs.fit,
         -warnings, -errors, -temp, -date) %>%
  rename(GPP_sm = GPP, ER_sm = ER, k600_sm = K600.daily) %>%
  pivot_longer(cols = GPP_sm:GPP_mod,
               names_to = c("type", "model"),
               names_pattern = "(.*)_(.*)") %>%
  pivot_wider(names_from = model, values_from = value, values_fn = mean) %>%
  mutate(pctdiff = (sm - mod) / mod)

# Plot the data
ggplot(data = df_com_l,
       aes(x = D,
           y = alpha,
           fill = pctdiff *100)) +
  facet_wrap(~type) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  geom_tile() +
  scale_fill_viridis_c(name = "% diff.", direction = -1) +
  theme_bw() +
  labs(x = expression("dispersion ("*m^2~h^{-1}*")"),
       y = expression("storage exchange ("*h^{-1}*")"))


