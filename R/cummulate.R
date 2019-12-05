csvout <- function(x, s) {
  y <- do.call("rbind", lapply(sort(unique(x[, s])), function(i) {
    cbind(
      unique(x[, "agegrp"]), min(x[(x[, s] == i), "yw"]), max(x[(x[, s] == i), "yw"]),
      round(tail(
        x[
          (x[, s] == i),
          c(
            paste0("cEdIA_", s), paste0("cEdIA_", s, "_95L"), paste0("cEdIA_", s, "_95U"),
            paste0("cEdET_", s), paste0("cEdET_", s, "_95L"), paste0("cEdET_", s, "_95U"),
            paste0("cexcess_", s), paste0("cuexcess_", s), "N"
          )
        ],
        n = 1
      ))
    )
  }))
  colnames(y) <- c(
    "AgeGroup", "Start", "End", "cIA", "cIA_95L", "cIA_95U",
    "cET", "cET_95L", "cET_95U", "cexcess", "cuexcess", "N"
  )
  return(y)
}

cummulate <- function(results) {
  res <- as.data.frame(results)
  res$yw <- sprintf("%04iw%02i", res$year, res$week)

  retval <- list()
  for (s in c("summer", "winter", "year")) {
    out <- do.call("rbind", lapply(0:4, function(a) {
      csvout(res[(!is.na(res[, s])) & (res$agegrp == a), ], s)
    }))
    out$AgeGroup <- factor(
      out$AgeGroup,
      levels = c(0, 1, 2, 3, 4),
      labels = c("0-4 years", "5-14 years", "15-64 years", "Aged 65", "Total")
    )
    retval[[s]] <- out[, c(
      "AgeGroup", "Start", "End", "cIA", "cIA_95L", "cIA_95U",
      "cET", "cET_95L", "cET_95U", "cexcess", "cuexcess"
    )]
  }
  retval
}
