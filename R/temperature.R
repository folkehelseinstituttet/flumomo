#' download_weather
#' @param country_code a
#' @export
download_weather <- function(country_code) {
  tmp <- tempfile()
  url <- paste0("http://euromomo.eu/methods/weather/wdata_", country_code, ".txt")
  utils::download.file(url, tmp)
  retval <- fread(tmp)
  retval <- retval[, c("date", "pop3", "NUTS3", "temp")]

  return(retval)
}

# Temperature data from WeatherData
# PARM ET data containing then variables: date, pop3, NUTS3, temp
# PARM start_year first year to be included. Must be >= 2000
# PARM start_week first week in start_year to be included
# PARM end_year last year to be included
# PARM end_week last week in end_year to be included
TemperatureData <- function(data_weather, start_year, start_week, end_year, end_week) {
  retval <- copy(data_weather)
  setDT(retval)
  retval[, date := as.Date(date)]
  retval[is.na(pop3), pop3 <- 1]

  retval <- retval[, .(
    temp = mean(temp, na.rm = T),
    pop3 = mean(pop3, na.rm = T)
  ), keyby = .(NUTS3, date)]

  retval[, pop3_sum := sum(pop3), by = date]

  retval <- retval[, .(
    temp = sum(temp * pop3) / sum(pop3_sum)
  ), keyby = .(date)]

  # restrict
  isoyrwk_start <- paste0(start_year, "-W", formatC(start_week, width = 2, flag = "0"))
  isoyrwk_end <- paste0(end_year, "-W", formatC(end_week, width = 2, flag = "0"))

  retval <- retval[date %in% seq.Date(as.Date(paste0(start_year, "/01/01")), as.Date(paste0(end_year + 1, "/01/01")), by = "day")]
  retval[, ISOweek := ISOweek::ISOweek(date)]

  retval <- retval[isoyrwk_start <= ISOweek & ISOweek <= isoyrwk_end]

  # aggregate over week
  retval <- retval[, .(
    temp = mean(temp, na.rm = T),
    tmin = min(temp, na.rm = T),
    tmax = min(temp, na.rm = T)
  ), keyby = .(ISOweek)]

  retval[, wk := 1:.N]
  retval[, sin52 := sin((2 * pi / (365.25 / 7)) * wk)]
  retval[, cos52 := cos((2 * pi / (365.25 / 7)) * wk)]

  # imputing
  retval$ptemp <- stats::predict(stats::glm(temp ~ sin52 + cos52, data = retval[!(is.na(temp) | is.infinite(temp)), ]), retval)
  retval$ptmin <- stats::predict(stats::glm(tmin ~ sin52 + cos52, data = retval[!(is.na(retval$tmin) | is.infinite(retval$tmin)), ]), retval)
  retval$ptmax <- stats::predict(stats::glm(tmax ~ sin52 + cos52, data = retval[!(is.na(retval$tmax) | is.infinite(retval$tmax)), ]), retval)
  retval[, c("wk", "cos52", "sin52") := NULL]

  retval[, t := temp]
  retval[is.na(temp) | is.infinite(temp), t := ptemp]

  retval[, ET := (t - ptmax) * (t > ptmax) + (t - ptmin) * (t < ptmin)]

  retval[, year := as.numeric(substr(ISOweek, 1, 4))]
  retval[, week := as.numeric(substr(ISOweek, 7, 8))]

  retval[, c("ISOweek", "t") := NULL]
  return(retval)
}
