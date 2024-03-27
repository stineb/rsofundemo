library(dplyr)
library(tidyr)
library(rsofun)
library(here)
library(ncdf4)
library(cwd)

# function used to create dataframes
get_annual_aet_pet <- function(df){
  df |>
    mutate(year = lubridate::year(date)) |>
    group_by(year) |>
    summarise(aet = sum(aet),
              pet = sum(pet)) |>
    ungroup() |>
    summarise(aet = mean(aet),
              pet = mean(pet))
}

get_annual_prec_cond <- function(df){
  df |>
    mutate(year = lubridate::year(date)) |>
    group_by(year) |>
    summarise(prec_cond = sum(prec_cond)) |>
    ungroup() |>
    summarise(prec_cond = mean(prec_cond))
}

get_annual_gpp_netrad <-  function(df){
  df |>
    mutate(year = lubridate::year(date)) |>
    group_by(year) |>
    summarise(gpp = sum(gpp),
               netrad = sum(netrad)) |>
    ungroup() |>
    summarise(gpp = mean(gpp),
              netrad = mean(netrad))
}

transfrom_le_ET <-  function(df){
  df |>
    mutate(ET = convert_et(le,temp,patm))  |>
    mutate(year = lubridate::year(date)) |>
    group_by(year) |>
    summarise(ET = sum(ET)) |>
    ungroup() |>
    summarise(ET = mean(ET))
}


# paramter for p model
params_modl <- list(
  kphio              = 0.04998,    # setup ORG in Stocker et al. 2020 GMD
  kphio_par_a        = 0.0,        # set to zero to disable temperature-dependence of kphio
  kphio_par_b        = 1.0,
  soilm_thetastar    = 0.6 * 240,  # to recover old setup with soil moisture stress
  soilm_betao        = 0.0,
  beta_unitcostratio = 146.0,
  rd_to_vcmax        = 0.014,      # value from Atkin et al. 2015 for C3 herbaceous
  tau_acclim         = 30.0,
  kc_jmax            = 0.41
)

cost_whc_driver <- var_whc_driver <- readRDS(here("data","rsofun_driver_data_v3.rds"))

# change whc to previous result

nc_file <- nc_open(here("data","whc_2m.nc"))

whc = ncvar_get(nc_file, "whc_2m")
lons = ncvar_get(nc_file, "lon")
lats = ncvar_get(nc_file, "lat")

geo <- cost_whc_driver |>
  unnest(site_info) |>
  select(lon  , lat)

geo$sitename <- cost_whc_driver$sitename

n <- 1 # parameter to select size of slice to average

old_whc <- lapply(geo$sitename, function(x){
  tmp <- geo[geo$sitename == x,]
  lonid <- which(lons > tmp$lon)[1]
  latid <- which(lats > tmp$lat)[1]
  whc_grid <- whc[(lonid-n):(lonid+n), (latid-n):(latid+n)]
  whc_site <- mean(as.numeric(whc_grid, na.rm=T))
  return(whc_site)
})

old_whc = unlist(old_whc)

for(i in 1:dim(cost_whc_driver)[1]){
  cost_whc_driver$site_info[i][[1]][4] <- old_whc[i]
}
rm(whc) # for memory

# run p model

cost_whc_output <- rsofun::runread_pmodel_f(
  cost_whc_driver,
  par = params_modl
)

# some simulations failed due to missing values in the forcing. If that happens,
# the number of rows in the output data is 1.
cost_whc_output <- cost_whc_output |>
  mutate(len = purrr::map_int(data, ~nrow(.))) |>
  filter(len != 1) |>
  select(-len)

# create dataframe
cost_whc_adf <- cost_whc_output |>
  mutate(cost_whc_adf = purrr::map(data, ~get_annual_aet_pet(.))) |>
  unnest(cost_whc_adf) |>
  select(sitename, aet,pet)

cost_whc_adf <- cost_whc_driver |>
  unnest(forcing) |>
  left_join(
    cost_whc_output |>
      unnest(data) |>
      select(-snow, -netrad, -fapar, gpp_pmodel = gpp),
    by = c("sitename", "date")
  ) |>
  mutate(prec = (rain + snow) * 60 * 60 * 24) |>
  mutate(prec_cond = prec + cond) |>
  group_by(sitename) |>
  nest() |>
  mutate(df = purrr::map(data, ~get_annual_prec_cond(.))) |>
  unnest(df) |>
  select(sitename, prec_cond) |>
  right_join(
    cost_whc_adf,
    by = "sitename"
  )

whc <-  cost_whc_output |>
  unnest(site_info) |>
  select(whc)

cost_whc_adf <- cost_whc_output |>
  mutate(df = purrr::map(data, ~get_annual_gpp_netrad(.))) |>
  unnest(df) |>
  select(sitename,gpp,netrad) |>
  right_join(
    cost_whc_adf,
    by = "sitename"
  )

cost_whc_adf <- cost_whc_driver |>
  mutate(df = purrr::map(forcing, ~transfrom_le_ET(.))) |>
  unnest(df) |>
  select(sitename,ET) |>
  right_join(
    cost_whc_adf,
    by = "sitename"
  )

whc <-  cost_whc_output |>
  unnest(site_info) |>
  select(whc)

cost_whc_adf$whc <- whc[[1]]

# run p model
var_whc_output <- rsofun::runread_pmodel_f(
  var_whc_driver,
  par = params_modl
)

# some simulations failed due to missing values in the forcing. If that happens,
# the number of rows in the output data is 1.
var_whc_output <- var_whc_output |>
  mutate(len = purrr::map_int(data, ~nrow(.))) |>
  filter(len != 1) |>
  select(-len)

# create dataframe
var_whc_adf <- var_whc_output |>
  mutate(var_whc_adf = purrr::map(data, ~get_annual_aet_pet(.))) |>
  unnest(var_whc_adf) |>
  select(sitename, aet, pet)

var_whc_adf <- var_whc_driver |>
  unnest(forcing) |>
  left_join(
    var_whc_output |>
      unnest(data) |>
      select(-snow, -netrad, -fapar, gpp_pmodel = gpp),
    by = c("sitename", "date")
  ) |>
  mutate(prec = (rain + snow) * 60 * 60 * 24) |>
  mutate(prec_cond = prec + cond) |>
  group_by(sitename) |>
  nest() |>
  mutate(df = purrr::map(data, ~get_annual_prec_cond(.))) |>
  unnest(df) |>
  select(sitename, prec_cond) |>
  right_join(
    var_whc_adf,
    by = "sitename"
  )

whc <-  var_whc_output |>
  unnest(site_info) |>
  select(whc)

var_whc_adf$whc <- whc[[1]]

var_whc_adf <- var_whc_output |>
  mutate(df = purrr::map(data, ~get_annual_gpp_netrad(.))) |>
  unnest(df) |>
  select(sitename,gpp,netrad) |>
  right_join(
    var_whc_adf,
    by = "sitename"
  )

var_whc_adf <- var_whc_driver |>
  mutate(df = purrr::map(forcing, ~transfrom_le_ET(.))) |>
  unnest(df) |>
  select(sitename,ET) |>
  right_join(
    var_whc_adf,
    by = "sitename"
  )

saveRDS(cost_whc_adf,here("data","costant_whc.rds"))
saveRDS(var_whc_adf,here("data","variable_whc.rds"))
