# create validation file used in p_model
library(dplyr)
library(FluxDataKit)
library(rsofun)

# insert all the files in the data-raw folder need also valid_years.csv (found in rsofun repo)

sites <- c("ES-Amo","FR-Pue")
# to use all sites 
# sites <- FluxDataKit::fdk_site_info

# the data-raw are the csv files and the lsm files obtained from flux data kit repo

files_csv <- list.files(here::here("data-raw//"))
lsm_path <- here::here("data-raw//")

valid_years <- read.csv(paste0(here::here("data-raw//"),"/valid_years_final.csv"), header = T, fileEncoding = "UTF-16")

validation <- lapply(sites,function(site){

  # Get filename for HH data for matching site
  file_csv <- files_csv[intersect(grep(site, files_csv),
                                 grep("HH", files_csv))]

  # get half-hourly data  --------------------------------------------------------
  # message("- convert to FLUXNET standard CSV file")
  hhdf <- suppressWarnings(
    try(
      fdk_convert_lsm(
        site = site,
        fluxnet_format = TRUE,
        path = lsm_path)))

  if(inherits(hhdf, "try-error")){
    message("!!! conversion to FLUXNET failed  !!!")
    return(NULL)}

  # Add date and time columns to hhdf for easier further processing.
  # ---------------------------------------------------------
  hhdf <- hhdf |>
    mutate(time = lubridate::as_datetime(as.character(TIMESTAMP_START), tz = "GMT", format="%Y%m%d%H%M")) |>
    mutate(date = lubridate::as_date(time))
  
  # Aggregate to daily 24-hr means  ----------------------------------------------------------
  message("- downsampling FLUXNET format - 24 hr means")
  ddf_24hr_mean <-
    try(
      hhdf |>
        group_by(date) |>
        select(-TIMESTAMP_START, -TIMESTAMP_END,-time) |>
        summarize_all(.funs = mean))

  # Get valid data years
  ystart <- valid_years %>% filter(Site==site) %>% pull(start_year)
  yend <- valid_years %>% filter(Site==site) %>% pull(end_year)

  # Write validation data
  message("- compiling validation")

  p_hydro_validation <- p_model_validation
  p_hydro_validation$sitename[[1]] = site

  p_hydro_validation$data <- ddf_24hr_mean |>
    dplyr::filter(lubridate::year(date) %in% ystart:yend) |>
    dplyr::filter(!(lubridate::mday(date) == 29 & lubridate::month(date) == 2)) |>
    group_by(date) |>
    summarise(
      date = date,
      gpp = GPP_DT_VUT_REF,
      gpp_unc = GPP_DT_VUT_SE
    ) |>
    mutate(gpp = gpp*86400/1e6*12)  |> # convert [umol m-2 s-1] to [gC m-2 s-1]
    mutate(gpp_unc = gpp_unc*86400/1e6*12)  |> # convert [umol m-2 s-1] to [gC m-2 s-1]
    list()
  return(p_hydro_validation)
})

validation <- dplyr::bind_rows(validation)

# save validation file, used in vignette
saveRDS(
  validation, paste0(here::here("data//"),"rsofun_validation_data.rds"),
  compress = "xz")
