##################################
### Version 4.2 - October 2018 ###
##################################
### Any questions or problems please contact: Jens Nielsen, Phone: +45 32683965 E-mail: nls@ssi.dk

# Install required packages
# install.packages("ISOweeks")

estimate <- function(
                     country,
                     country_code,
                     start_year,
                     start_week,
                     end_year,
                     end_week,
                     data_deaths,
                     data_ia,
                     data_weather,
                     data_population,
                     IArest,
                     IAlags,
                     ETlags) {
  ptrend <- 0.05
  p26 <- 0.05
  p52 <- 0.10

  data_deaths <- copy(data_deaths)
  data_ia <- copy(data_ia)
  data_weather <- copy(data_weather)
  if (!is.null(data_population)) data_population <- copy(data_population)

  ### START: deaths data ###
  data_deaths$deaths <- pmax(data_deaths$nb, data_deaths$nbc)
  data_deaths$agegrp <- match(data_deaths$group, c("0to4", "5to14", "15to64", "65P", "Total")) - 1 # Convert group to integer (0-5)
  data_deaths <- data_deaths[!is.na(agegrp), c("agegrp", "YoDi", "WoDi", "deaths")]
  names(data_deaths)[2:3] <- c("year", "week")

  data_deaths <- data_deaths[(start_year * 100 + start_week <= year * 100 + week) & (year * 100 + week <= end_year * 100 + end_week)]
  data_deaths[is.na(deaths), deaths := 0]
  setorder(data_deaths, agegrp, year, week)
  ### END: deaths data ###

  ### Start: Temperature data ###
  ET <- TemperatureData(data_weather, start_year, start_week, end_year, end_week)
  ### END: Temperature data ###

  ### START: data_population
  if (!is.null(data_population)) {
    data_population[, agegrp := match(group, c("0to4", "5to14", "15to64", "65P", "Total")) - 1] # Convert group to integer (0-5)
    data_population[, group := NULL]
  }

  ### START: data_ia ###
  data_ia[, agegrp := match(group, c("0to4", "5to14", "15to64", "65P", "Total")) - 1] # Convert group to integer (0-5)
  data_ia[, group := NULL]

  ### START: Merge data ###
  data_big <- merge(data_deaths, data_ia, by = c("agegrp", "year", "week"), all = FALSE)
  data_big <- merge(data_big, ET[, c("year", "week", "ET")], by = c("year", "week"), all = FALSE)
  if (!is.null(data_population)) {
    data_big <- merge(data_big, data_population, by = c("year", "agegrp"), all = FALSE)
  } else {
    data_big[, N := 1]
  }

  ### END: Merge data ###

  ### START: Prepare data ###


  data_big <- data_big[(year * 100 + week >= start_year * 100 + start_week) & (year * 100 + week <= end_year * 100 + end_week)]
  setorder(data_big, agegrp, year, week)
  data_big[, season := year - (week < 27)]
  data_big[, summer := week >= 21 & week <= 39]
  data_big[, winter := !summer]

  results4 <- as.data.frame(data_big)

  # lags IA
  for (a in 0:4) {
    for (s in sort(unique(results4[results4$agegrp == a, "season"]))) {
      for (d in 0:IAlags) {
        results4[results4$agegrp == a, paste("d", d, "_IA", s, sep = "")] <- vecshift(results4[results4$agegrp == a, "IA"], -d) * as.numeric(results4[results4$agegrp == a, "season"] == s)
      }
    }
  }
  # warm/cold summer/winter and lags
  for (s in c("summer", "winter")) {
    results4[, paste0("cold_", s)] <- with(results4, -((ET < 0) & get(s)) * ET)
    results4[, paste0("warm_", s)] <- with(results4, ((ET > 0) & get(s)) * ET)
    for (t in c("cold", "warm")) {
      for (d in 0:ETlags) {
        results4[, paste0("d", d, "_", t, "_", s)] <- NA
        for (a in 0:4) {
          results4[results4$agegrp == a, paste0("d", d, "_", t, "_", s)] <- vecshift(results4[results4$agegrp == a, paste0(t, "_", s)], -d)
        }
      }
    }
    results4[, paste0("cold_", s)] <- NULL
    results4[, paste0("warm_", s)] <- NULL
  }
  results4$summer <- NULL
  results4$winter <- NULL

  # Create wk
  results4 <- results4[with(results4, order(agegrp, year, week)), ]
  results4$wk <- unlist(lapply(table(results4$agegrp), seq))

  results4$sin52 <- sin((2 * pi / (365.25 / 7)) * results4$wk)
  results4$cos52 <- cos((2 * pi / (365.25 / 7)) * results4$wk)
  results4$sin26 <- sin((4 * pi / (365.25 / 7)) * results4$wk)
  results4$cos26 <- cos((4 * pi / (365.25 / 7)) * results4$wk)

  results4[is.na(results4)] <- 0
  ### END: Prepare data ###

  ### START: Estimation ###
  results4$EB <- NA
  results4$VlogB <- NA
  results4$EIA <- NA
  results4$VlogIA <- NA
  results4$EET <- NA
  results4$VlogET <- NA
  results4$Vdeaths <- NA

  for (a in 0:4) {
    print(paste("### Age group ", a, "###"))
    wk <- "wk"
    f <- paste(c(
      "deaths ~ ", wk, " sin52 + cos52 + sin26 + cos26",
      grep("^d[0-9]", names(results4), value = TRUE)
    ), collapse = " + ")
    m <- try(stats::glm(f, stats::quasipoisson, offset = log(N), data = results4[results4$agegrp == a, ]))
    if (!inherits(m, "try-error") & m$converge & (stats::median(results4[results4$agegrp == a, "deaths"]) > 0)) {
      fa <- paste(c(
        "deaths ~ sin52 + cos52 + sin26 + cos26",
        grep("^d[0-9]", names(results4), value = TRUE)
      ), collapse = " + ")
      ma <- stats::glm(fa, stats::quasipoisson, offset = log(N), data = results4[results4$agegrp == a, ])
      if (stats::anova(m, ma, dispersion = max(1, sum(stats::residuals(m, type = "deviance")^2) / stats::df.residual(m)), test = "LRT")$`Pr(>Chi)`[2] > ptrend) {
        wk <- ""
        m <- ma
      }
      fa <- paste(c("deaths ~ ", wk, grep("^d[0-9]", names(results4), value = TRUE)), collapse = " + ")
      ma <- stats::glm(fa, stats::quasipoisson, offset = log(N), data = results4[results4$agegrp == a, ])
      if (stats::anova(m, ma, dispersion = max(1, sum(stats::residuals(m, type = "deviance")^2) / stats::df.residual(m)), test = "LRT")$`Pr(>Chi)`[2] > max(p26, p52)) {
        m <- ma
      } else {
        fa <- paste(c("deaths ~ ", wk, " + cos52 + sin52", grep("^d[0-9]", names(results4), value = TRUE)), collapse = " + ")
        ma <- stats::glm(fa, stats::quasipoisson, offset = log(N), data = results4[results4$agegrp == a, ])
        if (stats::anova(m, ma, dispersion = max(1, sum(stats::residuals(m, type = "deviance")^2) / stats::df.residual(m)), test = "LRT")$`Pr(>Chi)`[2] > p26) {
          m <- ma
        } else {
          fa <- paste(c("deaths ~ ", wk, grep("^d[0-9]", names(results4), value = TRUE)), collapse = " + ")
          ma <- stats::glm(fa, stats::quasipoisson, offset = log(N), data = results4[results4$agegrp == a, ])
          if (stats::anova(m, ma, dispersion = max(1, sum(stats::residuals(m, type = "deviance")^2) / stats::df.residual(m)), test = "LRT")$`Pr(>Chi)`[2] > p52) {
            m <- ma
          }
        }
      }
    } else {
      if (inherits(m, "try-error")) print("### Could not fit model ###")
      if (!m$converge) print("### Model did not converge ###")
      if (stats::median(results4[results4$agegrp == a, "deaths"]) == 0) print("### Zero inflated ###")
      print("### Simple model with only trend used ###")
      f <- paste(c("deaths ~ wk"))
      m <- stats::glm(f, stats::quasipoisson, offset = log(N), data = results4[results4$agegrp == a, ])
    }
    print(summary(m, dispersion = max(1, sum(stats::residuals(m, type = "deviance")^2) / stats::df.residual(m))))
    # Baseline
    results4.B <- results4[results4$agegrp == a, ]
    for (d in 0:IAlags) {
      for (s in min(results4.B$season):max(results4.B$season)) {
        results4.B[, paste("d", d, "_IA", s, sep = "")] <- 0
      }
      for (s in c("summer", "winter")) {
        results4.B[, paste("d", d, "_warm_", s, sep = "")] <- 0
        results4.B[, paste("d", d, "_cold_", s, sep = "")] <- 0
      }
    }
    results4[results4$agegrp == a, ]$EB <- exp(stats::predict.glm(m, newdata = results4.B, se.fit = TRUE)$fit)
    results4[results4$agegrp == a, ]$VlogB <- stats::predict.glm(m, newdata = results4.B, dispersion = max(1, sum(stats::residuals(m, type = "deviance")^2) / stats::df.residual(m)), se.fit = TRUE)$se.fit^2
    # IA
    results4.IA <- results4[results4$agegrp == a, ]
    for (s in c("summer", "winter")) {
      for (d in 0:ETlags) {
        results4.IA[, paste("d", d, "_warm_", s, sep = "")] <- 0
        results4.IA[, paste("d", d, "_cold_", s, sep = "")] <- 0
      }
    }
    results4[results4$agegrp == a, ]$EIA <- exp(stats::predict.glm(m, newdata = results4.IA, se.fit = TRUE)$fit)
    results4[results4$agegrp == a, ]$VlogIA <- stats::predict.glm(m, newdata = results4.IA, dispersion = max(1, sum(stats::residuals(m, type = "deviance")^2) / stats::df.residual(m)), se.fit = TRUE)$se.fit^2
    # ET
    results4.ET <- results4[results4$agegrp == a, ]
    for (d in 0:IAlags) {
      for (s in sort(unique(results4.ET$season))) {
        results4.ET[, paste("d", d, "_IA", s, sep = "")] <- 0
      }
    }
    results4[results4$agegrp == a, ]$EET <- exp(stats::predict.glm(m, newdata = results4.ET, se.fit = T)$fit)
    results4[results4$agegrp == a, ]$VlogET <- stats::predict.glm(m, newdata = results4.ET, dispersion = max(1, sum(stats::residuals(m, type = "deviance")^2) / stats::df.residual(m)), se.fit = T)$se.fit^2

    results4[results4$agegrp == a, ]$Vdeaths <- with(results4[results4$agegrp == a, ], EB * max(1, sum(stats::residuals(m, type = "deviance")^2) / stats::df.residual(m)))
  }

  # Delete external objects
  rm(a, f, fa, s, d, m, ma, results4.B, results4.ET, results4.IA, ptrend, p26, p52, wk)
  ### END: Estimation ###

  ### START: Post estimation ###
  # Keep relevant variables #
  results4 <- results4[, c("agegrp", "year", "week", "deaths", "N", "IA", "ET", "Vdeaths", "EB", "VlogB", "EIA", "VlogIA", "EET", "VlogET")]

  # Baseline, EIA and EET estimation variances - not on log scale
  results4$VB <- with(results4, (exp(VlogB) - 1) * exp(2 * log(EB) + VlogB))
  results4$VIA <- with(results4, (exp(VlogIA) - 1) * exp(2 * log(EIA) + VlogIA))
  results4$VET <- with(results4, (exp(VlogET) - 1) * exp(2 * log(EET) + VlogET))

  # Effects of IA and ET
  results4$EdIA <- with(results4, EIA - EB)
  results4[is.na(results4$EdIA), "EdIA"] <- 0
  results4$EdET <- with(results4, EET - EB)
  results4[is.na(results4$EdET), "EdET"] <- 0
  # Excess relative to baseline
  results4$excess <- with(results4, deaths - EB)
  # Unexplained excess
  results4$uexcess <- with(results4, deaths - (EB + EdIA + EdET))
  # Exclude negative IA effects
  if (IArest) results4$uexcess <- with(results4, deaths - (EB + pmax(0, EdIA) + EdET))

  # Baseline 2/3 residual confidence intervals
  results4$RVB <- with(results4, ((2 / 3) * (EB^(2 / 3 - 1))^2) * Vdeaths + ((2 / 3) * (EB^(2 / 3 - 1))^2) * VB)
  results4[with(results4, is.na(RVB) | is.infinite(RVB)), "RVB"] <- 0
  results4$EB_95L <- with(results4, pmax(0, sign((sign(EB) * abs(EB)^(2 / 3)) - 1.96 * sqrt(RVB)) * abs((sign(EB) * abs(EB)^(2 / 3)) - 1.96 * sqrt(RVB))^(3 / 2)))
  results4$EB_95U <- with(results4, sign((sign(EB) * abs(EB)^(2 / 3)) + 1.96 * sqrt(RVB)) * abs((sign(EB) * abs(EB)^(2 / 3)) + 1.96 * sqrt(RVB))^(3 / 2))

  # EdIA 2/3 residual confidence intervals
  results4$RVdIA <- with(results4, ((2 / 3) * (EB^(2 / 3 - 1))^2) * VB + ((2 / 3) * (EIA^(2 / 3 - 1))^2) * VIA)
  results4[with(results4, is.na(RVdIA) | is.infinite(RVdIA)), "RVdIA"] <- 0
  results4$EdIA_95L <- with(results4, sign((sign(EdIA) * abs(EdIA)^(2 / 3)) - 1.96 * sqrt(RVdIA)) * abs((sign(EdIA) * abs(EdIA)^(2 / 3)) - 1.96 * sqrt(RVdIA))^(3 / 2))
  results4$EdIA_95U <- with(results4, pmax(EdIA_95L, sign((sign(EdIA) * abs(EdIA)^(2 / 3)) + 1.96 * sqrt(RVdIA)) * abs((sign(EdIA) * abs(EdIA)^(2 / 3)) + 1.96 * sqrt(RVdIA))^(3 / 2)))

  # EdET 2/3 residual confidence intervals
  results4$RVdET <- with(results4, ((2 / 3) * (EB^(2 / 3 - 1))^2) * VB + ((2 / 3) * (EET^(2 / 3 - 1))^2) * VET)
  results4[with(results4, is.na(RVdET) | is.infinite(RVdET)), "RVdET"] <- 0
  results4$EdET_95L <- with(results4, sign((sign(EdET) * abs(EdET)^(2 / 3)) - 1.96 * sqrt(RVdET)) * abs((sign(EdET) * abs(EdET)^(2 / 3)) - 1.96 * sqrt(RVdET))^(3 / 2))
  results4$EdET_95U <- with(results4, pmax(EdET_95L, sign((sign(EdET) * abs(EdET)^(2 / 3)) + 1.96 * sqrt(RVdET)) * abs((sign(EdET) * abs(EdET)^(2 / 3)) + 1.96 * sqrt(RVdET))^(3 / 2)))

  results4$summer <- with(results4, ifelse(week >= 21 & week <= 39, year, NA))
  results4$winter <- with(results4, ifelse(week <= 20 | week >= 40, year - (week <= 20), NA))

  # Order
  results4 <- results4[with(results4, order(agegrp, year, week)), ]


  # EdIA and EdET residual variances
  for (s in c("summer", "winter", "year")) {
    nn <- !is.na(results4[, s]) # logical vector indicating season
    for (a in 0:4) {
      ag <- (results4$agegrp == a)
      for (v in c("excess", "uexcess")) {
        results4[nn & ag, paste0("c", v, "_", s)] <- NA
        for (sy in sort(unique(results4[, s]))) {
          results4[(results4[, s] == sy) & nn & ag, paste0("c", v, "_", s)] <-
            cummulate(results4[(results4[, s] == sy) & nn & ag, c(v, "EdIA")], FALSE)
        }
      }
      for (v in c("B", "IA")) {
        results4[nn & ag, paste0("cE", v, "_", s)] <- NA
        results4[nn & ag, paste0("cV", v, "_", s)] <- NA
        for (sy in sort(unique(results4[, s]))) {
          results4[(results4[, s] == sy) & nn & ag, paste0("cE", v, "_", s)] <-
            cummulate(results4[(results4[, s] == sy) & nn & ag, c(paste0("E", v), "EdIA")], IArest)
          results4[(results4[, s] == sy) & nn & ag, paste0("cV", v, "_", s)] <-
            cummulate(results4[(results4[, s] == sy) & nn & ag, c(paste0("V", v), "EdIA")], IArest)
        }
      }
      results4[nn & ag, paste0("cEd", v, "_", s)] <- NA
      results4[nn & ag, paste0("cEd", v, "_", s, "_95L")] <- NA
      results4[nn & ag, paste0("cEd", v, "_", s, "_95U")] <- NA
      results4[nn & ag, c(paste0("cEd", v, "_", s), paste0("cEd", v, "_", s, "_95L"), paste0("cEd", v, "_", s, "_95U"))] <-
        cCI(results4[nn & ag, c(paste0("cEB_", s), paste0("cVB_", s), paste0("cE", v, "_", s), paste0("cV", v, "_", s))])
      for (v in c("B", "ET")) {
        results4[nn & ag, paste0("cE", v, "_", s)] <- NA
        results4[nn & ag, paste0("cV", v, "_", s)] <- NA
        for (sy in sort(unique(results4[, s]))) {
          results4[(results4[, s] == sy) & nn & ag, paste0("cE", v, "_", s)] <-
            cummulate(results4[(results4[, s] == sy) & nn & ag, c(paste0("E", v), "EdIA")], FALSE)
          results4[(results4[, s] == sy) & nn & ag, paste0("cV", v, "_", s)] <-
            cummulate(results4[(results4[, s] == sy) & nn & ag, c(paste0("V", v), "EdIA")], FALSE)
        }
      }
      results4[nn & ag, paste0("cEd", v, "_", s)] <- NA
      results4[nn & ag, paste0("cEd", v, "_", s, "_95L")] <- NA
      results4[nn & ag, paste0("cEd", v, "_", s, "_95U")] <- NA
      results4[nn & ag, c(paste0("cEd", v, "_", s), paste0("cEd", v, "_", s, "_95L"), paste0("cEd", v, "_", s, "_95U"))] <-
        cCI(results4[nn & ag, c(paste0("cEB_", s), paste0("cVB_", s), paste0("cE", v, "_", s), paste0("cV", v, "_", s))])
    }
  }

  results4$country <- country
  results4$IArestricted <- as.integer(IArest)
  results4 <- results4[, c(
    "country", "IArestricted", "agegrp", "year", "week", "deaths", "Vdeaths", "N", "IA", "ET",
    "EB", "EB_95L", "EB_95U", "VB",
    "EIA", "VIA", "EET", "VET",
    "EdIA", "EdIA_95L", "EdIA_95U",
    "EdET", "EdET_95L", "EdET_95U",
    "cexcess_year", "cuexcess_year",
    "cEdIA_year", "cEdIA_year_95L", "cEdIA_year_95U",
    "cEdET_year", "cEdET_year_95L", "cEdET_year_95U",
    "summer", "cexcess_summer", "cuexcess_summer",
    "cEdIA_summer", "cEdIA_summer_95L", "cEdIA_summer_95U",
    "cEdET_summer", "cEdET_summer_95L", "cEdET_summer_95U",
    "winter", "cexcess_winter", "cuexcess_winter",
    "cEdIA_winter", "cEdIA_winter_95L", "cEdIA_winter_95U",
    "cEdET_winter", "cEdET_winter_95L", "cEdET_winter_95U"
  )]

  return(results4)
}
### END: Post estimation ###

# Function to get lag, forward and backward
vecshift <- function(x, shift = 0) {
  if (shift == 0) {
    return(x)
  }
  if (shift > 0) {
    return(c(x[-(1:shift)], rep(NA, shift)))
  }
  if (shift < 0) {
    return(c(rep(NA, -shift), x[1:(length(x) + shift)]))
  }
}

cummulate <- function(x, pos = FALSE) {
  # x = variable, EdIA
  y <- x[, 1]
  if (pos) y <- x[, 1] * (x[, 2] > 0)
  return(cumsum(y))
}
cCI <- function(x) {
  # x = cEB, cVB, cE*, cV*
  colnames(x) <- c("cEB", "cVB", "cE", "cV")
  x$cEd <- x$cE - x$cEB
  x$cRVd <- ((2 / 3) * (x$cEB^(2 / 3 - 1))^2) * x$cVB + ((2 / 3) * (x$cE^(2 / 3 - 1))^2) * x$cV
  x[with(x, is.na(cRVd) | is.infinite(cRVd)), "cRVd"] <- 0
  x$cCI_95L <- with(x, sign((sign(cEd) * abs(cEd)^(2 / 3)) - 1.96 * sqrt(cRVd)) * abs((sign(cEd) * abs(cEd)^(2 / 3)) - 1.96 * sqrt(cRVd))^(3 / 2))
  x$cCI_95U <- with(x, sign((sign(cEd) * abs(cEd)^(2 / 3)) + 1.96 * sqrt(cRVd)) * abs((sign(cEd) * abs(cEd)^(2 / 3)) + 1.96 * sqrt(cRVd))^(3 / 2))
  return(x[, c("cEd", "cCI_95L", "cCI_95U")])
}
