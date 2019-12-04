
#' run
#' @param country a
#' @param country_code a
#' @param start_year a
#' @param start_week a
#' @param end_year a
#' @param end_week a
#' @param data_deaths a
#' @param data_ia a
#' @param data_weather a
#' @param data_population a
#' @param IArest a
#' @param IAlags a
#' @param ETlags a
#' @export
run <- function(
                country,
                country_code,
                start_year = 2013,
                start_week = 40,
                end_year,
                end_week = 40,
                data_deaths,
                data_ia,
                data_weather = download_weather(country_code),
                data_population = NULL,
                IArest = TRUE,
                IAlags = 2,
                ETlags = 2) {
  stopifnot(sum(!c("group", "YoDi", "WoDi", "nb", "nbc") %in% names(data_deaths)) == 0)
  stopifnot(sum(!c("date", "pop3", "NUTS3", "temp") %in% names(data_weather)) == 0)
  if (!is.null(data_population)) stopifnot(sum(!c("group", "year", "N") %in% names(data_population)) == 0)
  stopifnot(sum(!c("group", "year", "week", "IA") %in% names(data_ia)) == 0)

  setDT(data_deaths)
  setDT(data_weather)
  if (!is.null(data_population)) setDT(data_population)
  setDT(data_ia)

  data_deaths <- data_deaths[, c("group", "YoDi", "WoDi", "nb", "nbc")]
  data_weather <- data_weather[, c("date", "pop3", "NUTS3", "temp")]
  data_population <- data_population[, c("group", "year", "N")]
  data_ia <- data_ia[, c("group", "year", "week", "IA")]

  res <- estimate(
    country = country,
    country_code = country_code,
    start_year = start_year,
    start_week = start_week,
    end_year = end_year,
    end_week = end_week,
    data_deaths = data_deaths,
    data_ia = data_ia,
    data_weather = data_weather,
    data_population = data_population,
    IArest = IArest,
    IAlags = IAlags,
    ETlags = ETlags
  )

  return(res)
}
