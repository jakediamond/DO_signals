#
# Authors: Jake Diamond, Robert Hensley, Matt Cohen
# Purpose: To define functions for comparison of DO model results with streammetabolizer
# Date: November 19, 2021
# 

# Load libraries
library(deSolve)
library(tidyverse)
library(ReacTran)
library(lubridate)

# Function to run the 1D model with choices of treatments
func_mod <- function(...){
  # Get the changes to parameters
  arguments = list(...)
  
  gpp_use = arguments[["gpp"]]
  er_use = arguments[["er"]]
  D_use = arguments[["D"]]
  gpp = ifelse(is.character(gpp_use), gpp_use, "constant")
  er = ifelse(is.character(er_use), er_use, "channel")
  D_mod = ifelse(is.numeric(D_use), D_use, 3600)
  
  # Prescribe dispersion
  # D_use = if(is.numeric(D)){D} else{d * sqrt(9.8 * d * S) * 3600} # (m2/h)
  
  # Dispersion grid (need to make a global variabel with <<-)
  D.grid <<- with(as.list(parms),
                 setup.prop.1D(value = D_mod,
                               grid = grid)
  )
  # Model with user specifications
  ode.1D(y = yini,
         times = times,
         func = model,
         parms = replace(parms,
                         names(arguments),
                         unlist(arguments)),
         nspec = 2,
         dimens = with(as.list(parms), L / dx),
         gpp_choice = gpp,
         er_loc_choice = er,
         er_choice = "constant")
}

# Turn model outputs into dataframes to compare with streammetabolizer
# function to do so
df_fun <- function(mod, reach_choice){
  as_tibble(mod) %>%
    pivot_longer(cols = -time, names_to = "key", values_to = "value") %>%
    separate(key, c("type", "reach"), "(?<=[A-Za-z])(?=[0-9])") %>%
    mutate(time_hr = as.numeric(time * del_t),
           reach = as.numeric(reach),
           dist = as.numeric(reach) * with(as.list(parms), dx)) %>%
    mutate(type_plot = recode(type,
                              `DO` = "DO~(mg~L^{-1})",
                              `DO_stor` = "DO[storage]~(mg~L^{-1})",
                              `ER.DO_stor` = "ER~(g~O^2~m^{-2}~h^{-1})",
                              `ER` = "ER~(g~O^2~m^{-2}~h^{-1})",
                              `GPP` = "GPP~(g~O^2~m^{-2}~h^{-1})",
                              `PAR` = "PAR~({`mu`}*mol~m^{-2}~s^{-1})"))
}

# Function to get data ready for streammetabolizer, pick a reach number
sm_fun <- function(data, reach_no){
  data %>%
    filter(type %in% c("PAR", "DO"),
           reach == reach_no) %>%
    select(time_hr, type, value) %>%
    pivot_wider(names_from = type, values_from = value) %>%
    rename(light = PAR,
           DO.obs = DO) %>%
    mutate(temp.water = 20,
           DO.sat = with(as.list(parms), do_sat),
           depth = with(as.list(parms), d),
           day_no = floor(time_hr / 24) + 1,
           time = paste(floor(time_hr %% 24), 
                        round((time_hr %% 24 - floor(time_hr %% 24)) * 60,
                              digits = 2), sep=":"),
           datetime = ymd_hm(paste0("2019-07-0", day_no, " ", time))) %>%
    select(-time_hr, -time, -day_no) 
  
}

# Run streamMetabolizer function
sm_res <- function(df_sm, modtype = "mle", gpp_fun = "satlight", er_fun = "constant", k = 7.5, ksig = 0.01){
  # Get the data and get the streamMetabolizer model name
  x = df_sm
  mod_name = mm_name(type=modtype, GPP_fun=gpp_fun, ER_fun=er_fun, pool_K600 = "none")
  # Get the SM model specs
  if(modtype=="mle"){
    mod_specs = specs(mod_name)
    if(is.numeric(k)){
      data.daily = mutate(x, date = date(datetime)) %>%
        distinct(date) %>%
        mutate(K600.daily = k)
    }
  } else {
    mod_specs = specs(mod_name)
    if(is.numeric(k)){
      mod_specs = specs(mod_name, K600_daily_meanlog = log(k), K600_daily_sdlog = ksig)
    }
  }
  # Get the data in good form to sm, time needs to be in solar time
  x$solar.time = lubridate::force_tz(x$datetime, 'Etc/GMT+0')
  x$solar.time = streamMetabolizer::calc_solar_time(x$solar.time, 
                                                    longitude=0)
  # Remove first and last four hours because sm starts at 4am
  x$datetime = NULL
  x = x[-c(1:4, 174:193), ]
  
  # do the metabolism model
  if(modtype == "mle"){
    mm = metab(mod_specs, data=x)
    if(is.numeric(k)){
      mm = metab(mod_specs, data=x, data_daily = data.daily)
    }
  } else {
    mm = metab(mod_specs, data=x)
  }
  mm_res = left_join(predict_metab(mm), get_params(mm) %>%
                       select(date, K600.daily))
}

# Extract GPP, ER, and K600 from 1D DO model
do_mod_res <- function(data, reach_no){
  filter(data, reach == reach_no, type %in% c("GPP", "ER")) %>%
    mutate(day_no = floor(time_hr / 24) + 1,
           date = ymd(paste0("2019-07-0", day_no))) %>%
    group_by(type, date) %>%
    filter(n() == 24) %>% #only full days
    summarize(mod = sum(value)) %>%
    pivot_wider(names_from = type, values_from = mod) %>%
    rename(GPP_mod = GPP, ER_mod = ER, date_mod = date)
  
}