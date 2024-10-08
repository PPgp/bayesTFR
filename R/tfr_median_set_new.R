tfr.shift.estimation.to.wpp <- function(sim.dir, ..., verbose = TRUE){
  mcmc.set <- get.tfr.mcmc(sim.dir)
  country_code <- NULL # for CRAN check not to complain
  if(!has.est.uncertainty(mcmc.set)) stop("Function can be only applied to estimation with uncertainty.")
  meta <- mcmc.set$meta
  wpp.year <- meta$wpp.year
  countries <- get.countries.table(mcmc.set)$code
  if(meta$annual.simulation && wpp.year < 2022) stop("No annual WPP data available.")

  wppdata <- data.table::data.table(load.from.wpp("tfr", wpp.year = wpp.year, annual = meta$annual.simulation))
  wppdatal <- melt(wppdata, id.vars = c("country_code", "name"), value.name = "wpp", variable.name = "period", variable.factor = FALSE)
  wppdatal$year <- if(meta$annual.simulation) as.integer(wppdatal$period) else as.integer(substr(wppdatal$period, 1, 4)) + 3
  meta$median.shift.estimation <- NULL
  if(verbose) cat("\n")
  for(icntry in seq_along(countries)) {
    if (verbose) {
      if(interactive()) cat("\rAdjusting countries' estimation to wpp", wpp.year, " ... ", round(icntry/length(countries) * 100), ' %')
      else {
        if (icntry == 1)
          cat("Adjusting countries' estimation to wpp", wpp.year, " ... ")
        cat(icntry, ", ")
      }
    }
    cntry <- countries[icntry]
    est.obj <- get.tfr.estimation(mcmc.set, country = cntry, probs=0.5, adjust =  FALSE, ...) # extract median estimates
    est.obj$tfr_quantile$period <- if(!meta$annual.simulation) paste(est.obj$tfr_quantile$year - 3, est.obj$tfr_quantile$year + 2, sep = "-") else as.character(est.obj$tfr_quantile$year)
    merged <- merge(data.table::data.table(est.obj$tfr_quantile), wppdatal[country_code == cntry], by = "period")
    merged$shift <- merged$wpp - merged$V1
    if(sum(merged$shift) != 0)
      meta$median.shift.estimation[[as.character(cntry)]] <- merged$shift
  }
  store.bayesTFR.meta.object(meta, meta$output.dir)
  predictions <- available.tfr.predictions(sim.dir = sim.dir)
  for(pred.subdir in predictions) { # need to save this also in the collapsed chain of the prediction
    pred <- get.tfr.prediction(sim.dir = sim.dir, subdir = pred.subdir)
    pred$mcmc.set$meta$median.shift.estimation <- meta$median.shift.estimation
    store.bayesTFR.prediction(pred)
    store.bayesTFR.meta.object(pred$mcmc.set$meta, pred$mcmc.set$meta$output.dir)
  }
  if(verbose) cat("\n")
  invisible(get.tfr.mcmc(sim.dir))
}

get.tfr.shift.estimation <- function(country.code, meta)
{
  if(is.null(meta$median.shift.estimation)) return(NULL)
  return(meta$median.shift.estimation[[as.character(country.code)]])
}

tfr.median.set.all <- function(sim.dir, country, values, years=NULL, burnin = 0, thin = 1, 
                               subdir = "predictions")
{
  mcmc.set <- get.tfr.mcmc(sim.dir)
  meta <- mcmc.set$meta
  country.obj <- get.country.object(country, meta=meta)
  has.predictions <- has.tfr.prediction(sim.dir = sim.dir, subdir = subdir)
  # browser()
  bdem.shift.estimation <- do.call('get.tfr.shift.estimation', list(country.obj$code, meta))
  estimation.years <- as.numeric(dimnames(meta$tfr_matrix)[[1]])
  nr.estimation <- length(estimation.years)
  bdem.shift <- NULL
  nr.proj <- 0
  has.uncertainty <- !is.null(mcmc.set$mcmc.list[[1]]$uncertainty) && mcmc.set$mcmc.list[[1]]$uncertainty
  years.missing <- is.null(years)
  if (has.predictions)
  {
    pred <- get.tfr.prediction(sim.dir = sim.dir, subdir = subdir)
    pred.years <- as.numeric(dimnames(pred$quantiles)[[3]])
    bdem.shift <- do.call('get.tfr.shift', list(country.obj$code, pred))
    nr.proj <- pred$nr.projections + 1
  }
  else if(!has.uncertainty) stop('No estimation and no prediction to adjust.')
  if(is.null(country.obj$name)) stop('Country not found.')
  if(is.null(years)) 
  {
    if (has.uncertainty) years <- c(years, estimation.years)
    if (has.predictions) years <- c(years, pred.years[2:nr.proj])
  }
  mid.years <- years
  lvalues <- length(values)
  if (lvalues != length(years)){
      if(!years.missing || lvalues > length(years))
            stop("Mismatch between values and years.")
      years <- years[1:lvalues]       
  }
  if(is.null(bdem.shift.estimation)) bdem.shift.estimation <- rep(0, nr.estimation)
  if(is.null(bdem.shift)) bdem.shift <- rep(0, nr.proj)
  which.estimation <- is.element(years, estimation.years)
  if (any(which.estimation))
  {
    if (!has.uncertainty)
      stop("Provided years for estimation but uncertainty is not considered in the modeling step")
    tfr.object <- get.tfr.estimation(mcmc.list = mcmc.set, country = country.obj$code, probs = 0.5, burnin = burnin, thin = thin)
    estimation.medians <- as.numeric(unlist(tfr.object$tfr_quantile[,1]))
    bdem.shift.estimation[is.element(estimation.years, years)] <- values[which.estimation] - estimation.medians[is.element(estimation.years, years)]
  }
  if (has.predictions) 
  {
    which.pred <- is.element(years, pred.years)
    pred.medians <- pred$quantiles[country.obj$index, '0.5',]
    if (any(which.pred))
    {
      bdem.shift[is.element(pred.years, years)] <- values[which.pred] - pred.medians[is.element(pred.years, years)]
    }
  }
  if (sum(bdem.shift) == 0) bdem.shift <- NULL
  if (sum(bdem.shift.estimation) == 0) bdem.shift.estimation <- NULL
  output <- list()
  if (has.uncertainty)
  {
    meta$median.shift.estimation[[as.character(country.obj$code)]] <- bdem.shift.estimation
    store.bayesTFR.meta.object(meta, meta$output.dir)
    output[['meta']] <- meta
    if(has.predictions) { # need to save this also in the collapsed chain of the prediction
        pred$mcmc.set$meta$median.shift.estimation <- meta$median.shift.estimation
        store.bayesTFR.prediction(pred)
        store.bayesTFR.meta.object(pred$mcmc.set$meta, pred$mcmc.set$meta$output.dir)
    }
  }
  if (has.predictions)
  {
    pred$median.shift[[as.character(country.obj$code)]] <- bdem.shift
    store.bayesTFR.prediction(pred)
    output[['pred']] <- pred
  }
  
  invisible(output)
}

tfr.median.reset.estimation <- function(sim.dir, countries = NULL)
{
  mcmc.set <- get.tfr.mcmc(sim.dir)
  meta <- mcmc.set$meta
  predictions <- available.tfr.predictions(sim.dir = sim.dir)
  has.estimation <- (!is.null(mcmc.set$mcmc.list[[1]]$uncertainty) && mcmc.set$mcmc.list[[1]]$uncertainty)
  output <- list()
  if (has.estimation && !is.null(meta$median.shift.estimation)) 
  {
    if(is.null(countries)) meta$median.shift.estimation <- NULL # reset all countries
    else {
      for (country in countries) 
      {
        country.obj <- get.country.object(country, meta=meta)
        meta$median.shift.estimation[[as.character(country.obj$code)]] <- NULL
      }
    }
    output[['meta']] <- meta
    store.bayesTFR.meta.object(meta, meta$output.dir)
    for(pred.subdir in predictions) { # need to save this also in the prediction mcmc object
        pred <- get.tfr.prediction(sim.dir = sim.dir, subdir = pred.subdir)
        pred$mcmc.set$meta$median.shift.estimation <- meta$median.shift.estimation
        store.bayesTFR.prediction(pred)
    }
  }
    for(pred.subdir in predictions) {
        pred <- get.tfr.prediction(sim.dir = sim.dir, subdir = pred.subdir)
        if(is.null(pred$median.shift)) next
        if(is.null(countries)) pred$median.shift <- NULL # reset all countries
        else {
            for (country in countries) 
            {
                country.obj <- get.country.object(country, meta=meta)
                pred$median.shift[[as.character(country.obj$code)]] <- NULL
            }
        }
        output[['pred']] <- pred # this just returns the last prediction object which is not correct if there are more than one prediction directories
        store.bayesTFR.prediction(pred)
  }
  invisible(output)
}

.do.shift.prediction.to.wpp <- function(pred, wpp.dataset, stat = "median", wpp.year = NULL, verbose = TRUE){
  country_code <- NULL # for CRAN check not to complain
  meta <- pred$mcmc.set$meta
  wpp.year <- if(!is.null(wpp.year)) wpp.year else meta$wpp.year
  countries <- get.countries.table(pred$mcmc.set)$code
  if(meta$annual.simulation && wpp.year < 2022) stop("No annual WPP data available.")
  
  wppdata <- data.table::data.table(load.from.wpp(wpp.dataset, wpp.year = wpp.year, annual = meta$annual.simulation))
  wppdatal <- melt(wppdata, id.vars = c("country_code", "name"), value.name = "wpp", variable.name = "period", variable.factor = FALSE)
  wppdatal$year <- if(meta$annual.simulation) as.integer(wppdatal$period) else as.integer(substr(wppdatal$period, 1, 4)) + 3
  
  pred.years <- as.numeric(dimnames(pred$quantiles)[[3]])
  pred$median.shift <- NULL
  
  if(verbose) cat("\n")
  for(icntry in seq_along(countries)) {
    if (verbose) {
      if(interactive()) cat("\rAdjusting countries' prediction to wpp", wpp.year, " ... ", round(icntry/length(countries) * 100), ' %')
      else {
        if (icntry == 1)
          cat("Adjusting countries' prediction to wpp", wpp.year, " ... ")
        cat(icntry, ", ")
      }
    }
    cntry <- countries[icntry]
    if(!stat %in% c("median", "mean")) stop("Argument 'stat' must be 'median' or 'mean', but is ", stat, ".")
    to.match <- merge(data.table::data.table(year = pred.years, median = pred$quantiles[icntry, "0.5", ]), 
                          wppdatal[country_code == cntry], by = "year", all.x = TRUE)
    if(stat == "mean") to.match[["median"]] <- to.match[["median"]] + pred$traj.mean.sd[icntry, 1, ] - to.match[["median"]] # difference between the mean and median
    to.match$wpp[is.na(to.match$wpp)] <- to.match$median[is.na(to.match$wpp)] # no shift for years that don't match
    to.match$shift <- to.match$wpp - to.match$median
    if(sum(to.match$shift) != 0)
      pred$median.shift[[as.character(cntry)]] <- to.match$shift
  }
  if(verbose) cat("\n")
  return(pred)
}


tfr.shift.prediction.to.wpp <- function(sim.dir, subdir = "predictions", ...){
  pred <- get.tfr.prediction(sim.dir, subdir = subdir)
  new.pred <- .do.shift.prediction.to.wpp(pred, wpp.dataset = "tfrprojMed", ...)
  store.bayesTFR.prediction(new.pred)
  invisible(new.pred)
}
