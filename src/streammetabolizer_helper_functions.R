#
# Authors: Jake Diamond, Robert Hensley, Matt Cohen
# Purpose: To define functions for comparison of DO model results with streammetabolizer
# Date: November 19, 2021
# 

# Load libraries
library(tidyverse)
library(streamMetabolizer)

# Function to get data ready for streammetabolizer, pick a reach number
sm_input_fun <- function(data, reach_no){
  data |>
    filter(type %in% c("PAR", "DO"),
           reach == reach_no) |>
    select(time, type, value) |>
    tidytable::pivot_wider(names_from = type, values_from = value) |>
    rename(light = PAR,
           DO.obs = DO) |>
    mutate(temp.water = 20,
           DO.sat = with(as.list(parms), O2_sat(temp)),
           depth = with(as.list(parms), d),
           day_no = floor(time / 24) + 1,
           time = paste(floor(time %% 24), 
                        round((time %% 24 - floor(time %% 24)) * 60,
                              digits = 2), sep=":"),
           datetime = ymd_hm(paste0("2019-07-0", day_no, " ", time))) |>
    select(-time, -day_no) 
}

# Run streamMetabolizer function
sm_fun <- function(input, modtype = "mle", gpp_fun = "linlight", 
                   er_fun = "constant", k = 7.5, ksig = 0.01) {
  
  # Get the delta t in hours
  dt <- as.numeric(difftime(input$datetime[2], 
                 input$datetime[1], 
                 units = "hours"))
  
  # The model specs for streammetabolizer
  mod_name <- mm_name(type = modtype, GPP_fun = gpp_fun, 
                      ER_fun = er_fun, pool_K600 = "none")
  
  # Get the SM model specs
  if(modtype=="mle"){
    mod_specs = specs(mod_name)
    if(is.numeric(k)){
      data.daily = mutate(input, date = date(datetime)) |>
        distinct(date) |>
        mutate(K600.daily = k)
    }
  } else {
    mod_specs = specs(mod_name)
    if(is.numeric(k)){
      mod_specs = specs(mod_name, K600_daily_meanlog = log(k), K600_daily_sdlog = ksig)
    }
  }
  
  # Get the data in good form to sm, time needs to be in solar time
  input$solar.time <- force_tz(input$datetime, 'Etc/GMT+0')
  input$solar.time <- streamMetabolizer::calc_solar_time(input$solar.time, 
                                                    longitude = 0)
  
  # Remove first and last four hours because sm starts at 4am
  input$datetime <- NULL
  input <- input[-c(1:(4/dt), (nrow(input) - (4/dt)):nrow(input)), ]
  
  # do the metabolism model
  if(modtype == "mle"){
    mm = metab(mod_specs, data = input)
    if(is.numeric(k)){
      mm = metab(mod_specs, data = input, data_daily = data.daily)
    }
  } else {
    mm = metab(mod_specs, data = input)
  }
  mm_res = left_join(predict_metab(mm), get_params(mm) |>
                       select(date, K600.daily))
}

# Extract GPP, ER, and K600 from 1D DO model
do_mod_res <- function(mod, reach_no){
  mod |>
  filter(reach == reach_no, type %in% c("GPP", "ER")) |>
    mutate(day_no = floor(time / 24) + 1,
           date = ymd(paste0("2019-07-0", day_no))) |>
    group_by(type, date) |>
    # filter(n() == 24) |> #only full days
    summarize(mod = mean(value*  24)) |>
    pivot_wider(names_from = type, values_from = mod) |>
    rename(GPP_mod = GPP, ER_mod = ER, date_mod = date)
  }
