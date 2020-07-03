require(mvtnorm)

get.obs.estimate.diff <- function(mcmc)
{
  mcmc$eps_unc <- list()
  for(country in 1:min(mcmc$meta$nr_countries, length(mcmc$meta$raw.data)))
  {
    if (is.null(mcmc$meta$raw.data[[country]])) next
    tmp.year <- mcmc$meta$raw.data[[country]]$Year
    tmp.tfr <- mcmc$meta$raw.data[[country]]$DataValue
    if (mcmc$meta$annual.simulation)
    {
      tmp.ind <- (tmp.year - mcmc$meta$start.year + 1)
      ref.tfr <- mcmc$meta$tfr_matrix_all[, country]
      interpolated.tfr <- ref.tfr[round(tmp.ind + 0.01)]
    }
    else
    {
      tmp.ind <- (tmp.year - mcmc$meta$start.year + 2) / 5
      left.distance <- tmp.ind - floor(tmp.ind)
      ref.tfr <- mcmc$meta$tfr_matrix_all[, country]
      interpolated.tfr <- ref.tfr[pmax(floor(tmp.ind), 1)] * (5 - left.distance) / 5 + 
        ref.tfr[pmin(pmax(floor(tmp.ind), 1) + 1, length(ref.tfr))] * left.distance / 5
    }
    delta <- tmp.tfr - interpolated.tfr
    mcmc$eps_unc[[country]] <- delta
  }
  return(mcmc)
}

get.obs.estimate.diff.original <- function(mcmc)
{
  mcmc$meta$raw_data.original$eps <- NA
  if (mcmc$meta$annual.simulation)
  {
    for(year in mcmc$meta$start.year:mcmc$meta$present.year)
    {
      year.ind <- year - mcmc$meta$start.year + 1
      if (length(mcmc$meta$ind.by.year[[year.ind]]) == 0) next
      ref.tfr <- mcmc$meta$tfr_matrix_all[year.ind, mcmc$meta$country.ind.by.year[[year.ind]]]
      idx <- mcmc$meta$ind.by.year[[year.ind]]
      mcmc$meta$raw_data.original$eps[idx] <- mcmc$meta$raw_data.original$DataValue[idx] - ref.tfr
    }
  }
  else 
  {
    years <- as.numeric(rownames(mcmc$meta$tfr_matrix_all))
    count <- 0
    for (year in c(years, rev(years)[1]+5))
    {
      count <- count + 1
      if (length(mcmc$meta$ind.by.year[[count]]) == 0) next
      idx <- mcmc$meta$ind.by.year[[count]]
      tfr.left <- mcmc$meta$tfr_matrix_all[max(count-1, 1), mcmc$meta$country.ind.by.year[[count]]]
      tfr.right <- mcmc$meta$tfr_matrix_all[min(count, nrow(mcmc$meta$tfr_matrix_all)), mcmc$meta$country.ind.by.year[[count]]]
      left.distance <- mcmc$meta$left.distance[[count]]
      ref.tfr <- (tfr.left * (5 - left.distance) + tfr.right * left.distance) / 5
      mcmc$meta$raw_data.original$eps[idx] <- mcmc$meta$raw_data.original$DataValue[idx] - ref.tfr
    }
  }
  return(mcmc)
}

update.obs.estimate.diff <- function(mcmc, country, year)
{
  if (mcmc$meta$annual.simulation)
  {
    indices <- which(round(mcmc$meta$raw_data_extra[[country]]$Year - mcmc$meta$start.year + 1.01) == year)
  }
  else
  {
    indices <- which((mcmc$meta$raw_data_extra[[country]]$Year < (year * 5 + mcmc$meta$start.year + 3)) & 
                       (mcmc$meta$raw_data_extra[[country]]$Year > (year * 5 + mcmc$meta$start.year -7)))
  }
  if (length(indices) > 0)
  {
    tmp.tfr <- mcmc$meta$raw_data_extra[[country]]$DataValue[indices]
    tmp.year <- mcmc$meta$raw_data_extra[[country]]$Year[indices]
    if (mcmc$meta$annual.simulation)
    {
      tmp.ind <- (tmp.year - mcmc$meta$start.year + 1)
      ref.tfr <- mcmc$meta$tfr_all[, country]
      interpolated.tfr <- ref.tfr[round(tmp.ind +0.01)]
    }
    else
    {
      tmp.ind <- (tmp.year - mcmc$meta$start.year + 2) / 5
      left.distance <- tmp.ind - floor(tmp.ind)
      ref.tfr <- mcmc$meta$tfr_all[, country]
      interpolated.tfr <- ref.tfr[pmax(floor(tmp.ind), 1)] * (5 - left.distance) / 5 + 
        ref.tfr[pmin(floor(tmp.ind) + 1, length(ref.tfr))] * left.distance / 5
    }
    return (tmp.tfr - interpolated.tfr)
  }
  else
  {
    return (numeric(0))
  }
}

estimate.bias.sd.raw <- function(mcmc)
{
  for(country in 1:min(mcmc$meta$nr_countries, length(mcmc$meta$raw.data)))
  {
    if (is.null(mcmc$meta$raw.data[[country]])) next
    source <- factor(as.character(mcmc$meta$raw.data[[country]]$DataProcess))
    estimation.method <- factor(as.character(mcmc$meta$raw.data[[country]]$Estimating.Methods))
    diff <- mcmc$eps_unc[[country]]
    regressor <- '1'
    if (length(levels(source)) > 1)
    {
      regressor <- paste0(regressor, ' + source')
    }
    if (length(levels(estimation.method)) > 1)
    {
      regressor <- paste0(regressor, ' + estimation.method')
    }
    m1 <- lm(as.formula(paste0('diff ~ ', regressor)))
    bias <- predict(m1)
    abs.residual <- abs(residuals(m1))
    m2 <- lm(as.formula(paste0('abs.residual ~ ', regressor)))
    std <- predict(m2) * sqrt(pi/2)
    std[std < 1e-6] <- 0.1
    std[std < (abs(bias) / 2)] <- abs(bias[std < (abs(bias) / 2)])
    mcmc$meta$raw.data[[country]]$bias <- bias
    mcmc$meta$raw.data[[country]]$std <- std
  }
  
  return(mcmc)
}

estimate.bias.sd.original <- function(mcmc, iso.unbiased=NULL, covariates=c('DataProcess', 'Estimating.Methods'), 
                                      cont_covariates=NULL)
{
  mcmc$meta$raw_data.original$bias <- NA
  mcmc$meta$raw_data.original$std <- NA
  if (is.null(mcmc$meta$bias_model))
  {
    mcmc$meta$bias_model <- list()
    mcmc$meta$std_model <- list()
  }
  for(country in 1:mcmc$meta$nr_countries)
  {
    ISO.code <- get.country.object(country, meta=mcmc$meta, index=TRUE)
    if (is.na(ISO.code$code)) next
    index.by.country <- which(mcmc$meta$raw_data.original$ISO.code == ISO.code$code)
    if (length(index.by.country) == 0) next
    regressor <- '1'
    diff <- mcmc$meta$raw_data.original$eps[index.by.country]
    if (length(covariates) > 0)
    {
      for (i in 1:length(covariates))
      {
        covariate <- covariates[i]
        assign(paste0('covariate_', i), factor(as.character(mcmc$meta$raw_data.original[index.by.country, covariate])))
        if (length(levels(get(paste0('covariate_', i)))) > 1)
          regressor <- paste0(regressor, ' + covariate_', i)
      }
    }
    if (length(cont_covariates) > 0)
    {
      for (i in 1:length(cont_covariates))
      {
        covariate <- cont_covariates[i]
        assign(paste0('cont_covariate_', i), mcmc$meta$raw_data.original[index.by.country, covariate])
        if (max(get(paste0('cont_covariate_', i))) - min(get(paste0('cont_covariate_', i))) > 1e-6)
        {
          regressor <- paste0(regressor, ' + cont_covariate_', i)
        }
      }
    }
    
    m1 <- lm(as.formula(paste0('diff ~ ', regressor)))
    bias <- predict(m1)
    abs.residual <- abs(residuals(m1))
    m2 <- lm(as.formula(paste0('abs.residual ~ ', regressor)))
    std <- predict(m2) * sqrt(pi/2)
    std[std < 1e-6] <- 0.1
    std[std < (abs(bias) / 2)] <- abs(bias[std < (abs(bias) / 2)]) 
    mcmc$meta$raw_data.original$bias[index.by.country] <- bias
    mcmc$meta$raw_data.original$std[index.by.country] <- std
    ## Optional
    if (ISO.code$code %in% iso.unbiased)
    {
      index.by.country.vr.estimate <- which((mcmc$meta$raw_data.original$ISO.code == ISO.code$code) & 
                                              (mcmc$meta$raw_data.original[, covariates[1]] %in% c("VR", 'Estimate')))
      mcmc$meta$raw_data.original$bias[index.by.country.vr.estimate] <- 0
      mcmc$meta$raw_data.original$std[index.by.country.vr.estimate] <- 0.016
    }
    mcmc$meta$bias_model[[country]] <- m1
    mcmc$meta$std_model[[country]] <- m2
  }
  
  
  return(mcmc)
}

get.eps.all.phases <- function(Dlpar, mcmc, country)
{
  eps_return <- numeric(length = mcmc$meta$T_end-1)
  id2 <- mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country] - 1)
  dl <- - Dlpar[5]/(1 + exp(- 2*log(9)/Dlpar[1] *(mcmc$meta$tfr_all[id2, country] - Dlpar[1] - Dlpar[2] -Dlpar[3] - Dlpar[4] + 0.5*Dlpar[1]))) + 
          Dlpar[5]/(1 + exp(- 2*log(9)/Dlpar[3] *(mcmc$meta$tfr_all[id2, country] - Dlpar[4] - 0.5*Dlpar[3])))
  dl <- dl * ifelse(mcmc$meta$annual.simulation, 1, 5)
  dl[mcmc$meta$tfr_all[id2, country] < 1] <- 0
  dl[dl<0] <- 0
  eps_return[id2] <- mcmc$meta$tfr_all[id2 + 1, country] - mcmc$meta$tfr_all[id2, country] + dl
  if (mcmc$meta$start_c[country] > 1)
  {
    eps_return[1:(mcmc$meta$start_c[country] - 1)] <- mcmc$meta$tfr_all[2:mcmc$meta$start_c[country], country] - 
      mcmc$meta$tfr_all[1:(mcmc$meta$start_c[country]-1), country]
  }
  if (mcmc$meta$lambda_c[country] < mcmc$meta$T_end)
  {
    id3 <- which(mcmc$meta$id_phase3 == country)
    eps_return[mcmc$meta$lambda_c[country]:(mcmc$meta$T_end - 1)] <- mcmc$meta$tfr_all[(mcmc$meta$lambda_c[country]+1):mcmc$meta$T_end, country] - 
      mcmc$meta$tfr_all[mcmc$meta$lambda_c[country]:(mcmc$meta$T_end-1), country] * mcmc$rho.c[id3] - (1-mcmc$rho.c[id3]) * mcmc$mu.c[id3]
  }
  return (eps_return)
}

get.sd.all.phases <- function(mcmc, country)
{
  sd_return <- mcmc$sd_Tc[, country]
  if (mcmc$meta$start_c[country] > 1)
  {
    sd_return[1:(mcmc$meta$start_c[country] - 1)] <- mcmc$sd_eps_tau
  }
  if (mcmc$meta$lambda_c[country] < mcmc$meta$T_end)
  {
    sd_return[mcmc$meta$lambda_c[country]:(mcmc$meta$T_end - 1)] <- mcmc$sigma.eps
  }
  return (sd_return)
}

mcmc.update.tfr <- function(country, mcmc)
{
  Dlpar <- c((mcmc$U_c[country]-mcmc$Triangle_c4[country])*
               exp(mcmc$gamma_ci[country,])/
               sum(exp(mcmc$gamma_ci[country,])), 
             mcmc$Triangle_c4[country], 
             mcmc$d_c[country])
  epsT.idx <- mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country]-1)
  eps_tfr_prev <- get.eps.all.phases(Dlpar, mcmc, country)
  sd_tfr_prev <- get.sd.all.phases(mcmc, country)
  eps_tfr_prop <- eps_tfr_prev
  sd_tfr_prop <- sd_tfr_prev
  if (country %in% mcmc$meta$id_phase3) {id3 <- which(mcmc$meta$id_phase3 == country)}
  for (year in 1:mcmc$meta$T_end)
  {
    tmp.years <- max(1, year-1):min(year, mcmc$meta$T_end -1)
    if (!is.null(mcmc$meta$ar.phase2) && mcmc$meta$ar.phase2)
    {
      tmp.extra.year <- max(2, year-1):min(year+1, mcmc$meta$T_end -1)
      tmp.extra.year <- tmp.extra.year[tmp.extra.year > mcmc$meta$start_c[country] & tmp.extra.year < mcmc$meta$lambda_c[country]]
      tmp.more.years <- tmp.years[!(tmp.years %in% tmp.extra.year)]
    }
    sd_f <- ifelse(mcmc$finished.iter>=100, 2.4 * max(mcmc$meta$tfr_sd_all[year, country], 0.01), 0.2)
    prop_tfr <- rnorm(1, mcmc$meta$tfr_all[year, country], sd_f)
    while (prop_tfr < 0.5) prop_tfr <- rnorm(1, mcmc$meta$tfr_all[year, country], sd_f)
    if (mcmc$meta$annual.simulation)
    {
      indices <- which(round(mcmc$meta$raw_data_extra[[country]]$Year - mcmc$meta$start.year + 1.01) == year)
    }
    else
    {
      indices <- which((mcmc$meta$raw_data_extra[[country]]$Year < (year * 5 + mcmc$meta$start.year + 3)) & 
                         (mcmc$meta$raw_data_extra[[country]]$Year > (year * 5 + mcmc$meta$start.year - 7)))
    }
    tfr_old <- mcmc$meta$tfr_all[year, country]
    eps_prev <- mcmc$eps_unc[[country]][indices]
    loglik_prev <- sum(dnorm(eps_prev, mean=mcmc$meta$raw_data_extra[[country]]$bias[indices], sd=mcmc$meta$raw_data_extra[[country]]$std[indices], log=TRUE))
    if (!is.null(mcmc$meta$ar.phase2) && mcmc$meta$ar.phase2)
    {
      loglik_prev <- loglik_prev + sum(dnorm(eps_tfr_prev[tmp.extra.year], mean = mcmc$rho.phase2 * eps_tfr_prev[tmp.extra.year - 1], 
                                             sd = sd_tfr_prev[tmp.extra.year], log=TRUE))
      loglik_prev <- loglik_prev + sum(dnorm(eps_tfr_prev[tmp.more.years], mean = 0, sd = sd_tfr_prev[tmp.more.years], log=TRUE))
    }
    else loglik_prev <- loglik_prev + sum(dnorm(eps_tfr_prev[tmp.years], mean = 0, sd = sd_tfr_prev[tmp.years], log=TRUE))
    mcmc$meta$tfr_all[year, country] <- prop_tfr
    eps_new <- update.obs.estimate.diff(mcmc, country, year)
    loglik_new <- sum(dnorm(eps_new, mean=mcmc$meta$raw_data_extra[[country]]$bias[indices], sd=mcmc$meta$raw_data_extra[[country]]$std[indices], log=TRUE))
    if (year < mcmc$meta$start_c[country])
    {
      eps_tfr_prop[tmp.years] <- mcmc$meta$tfr_all[tmp.years + 1, country] - mcmc$meta$tfr_all[tmp.years, country]
    }
    else if (year > mcmc$meta$lambda_c[country])
    {
      eps_tfr_prop[tmp.years] <- mcmc$meta$tfr_all[tmp.years + 1, country] - mcmc$meta$tfr_all[tmp.years, country] * 
        mcmc$rho.c[id3] - (1-mcmc$rho.c[id3]) * mcmc$mu.c[id3]
    }
    else 
      eps_tfr_prop <- get.eps.all.phases(Dlpar, mcmc, country)
    
    if ((year >= mcmc$meta$start_c[country]) && (year < mcmc$meta$lambda_c[country]))
    {
      tmp <- (prop_tfr - mcmc$S_sd) * ifelse(prop_tfr > mcmc$S_sd, -mcmc$a_sd, mcmc$b_sd)
      sd_tfr_prop[year] <- ifelse(mcmc$const_sd_dummie_Tc[year, country] == 1, mcmc$const_sd, 1) * 
        ifelse((mcmc$sigma0 + tmp > mcmc$meta$sigma0.min), mcmc$sigma0 + tmp, mcmc$meta$sigma0.min)
    }
    if (!is.null(mcmc$meta$ar.phase2) && mcmc$meta$ar.phase2)
    {
      loglik_new <- loglik_new + sum(dnorm(eps_tfr_prop[tmp.extra.year], mean = mcmc$rho.phase2 * eps_tfr_prop[tmp.extra.year - 1], 
                                             sd = sd_tfr_prop[tmp.extra.year], log=TRUE))
      loglik_new <- loglik_new + sum(dnorm(eps_tfr_prop[tmp.more.years], mean = 0, sd = sd_tfr_prop[tmp.more.years], log=TRUE))
    }
    else loglik_new <- loglik_new + sum(dnorm(eps_tfr_prop[tmp.years], mean = 0, sd = sd_tfr_prop[tmp.years], log=TRUE))
    
    # loglik_new <- loglik_new + sum(dnorm(eps_tfr_prop[tmp.years], mean=0, sd=sd_tfr_prop[tmp.years], log=TRUE))
    accept.prob <- exp(loglik_new - loglik_prev)
    if (runif(1) < accept.prob)
    {
      if (year %in% epsT.idx)
      {
        if (!is.null(mcmc$meta$ar.phase2) && mcmc$meta$ar.phase2 && year > mcmc$meta$start_c[country]) 
          mcmc$eps_Tc[year, country] <- eps_tfr_prop[year] - mcmc$rho.phase2 * eps_tfr_prop[year - 1]
        else mcmc$eps_Tc[year, country] <- eps_tfr_prop[year]
        mcmc$sd_Tc[year, country] <- sd_tfr_prop[year]
        mcmc$add_to_sd_Tc[year, country] <- (prop_tfr - mcmc$S_sd) * ifelse(prop_tfr > mcmc$S_sd, -mcmc$a_sd, mcmc$b_sd)
      }
      if ((year - 1) %in% epsT.idx)
      {
        if (!is.null(mcmc$meta$ar.phase2) && mcmc$meta$ar.phase2 && year > mcmc$meta$start_c[country] + 1) 
          mcmc$eps_Tc[year - 1, country] <- eps_tfr_prop[year - 1] - mcmc$rho.phase2 * eps_tfr_prop[year - 2]
        else mcmc$eps_Tc[year - 1, country] <- eps_tfr_prop[year - 1]
        mcmc$sd_Tc[year - 1, country] <- sd_tfr_prop[year - 1]
      }
      if (!is.null(mcmc$meta$ar.phase2) && mcmc$meta$ar.phase2 && (year + 1) %in% epsT.idx && year > mcmc$meta$start_c[country] - 1)
        mcmc$eps_Tc[year + 1, country] <- eps_tfr_prop[year + 1] - mcmc$rho.phase2 * eps_tfr_prop[year]
      
      mcmc$eps_unc[[country]][indices] <- eps_new
      eps_tfr_prev <- eps_tfr_prop
      sd_tfr_prev <- sd_tfr_prop
      
    }
    else
    {
      mcmc$meta$tfr_all[year, country] <- tfr_old
      sd_tfr_prop <- sd_tfr_prev
    }
    mu0 <- mcmc$meta$tfr_mu_all[year, country]
    sigma0 <- mcmc$meta$tfr_sd_all[year, country]
    mcmc$meta$tfr_mu_all[year, country] <- (mcmc$finished.iter * mu0 + mcmc$meta$tfr_all[year, country]) / (mcmc$finished.iter + 1)
    mu <- mcmc$meta$tfr_mu_all[year, country]
    mcmc$meta$tfr_sd_all[year, country] <- sqrt(((mcmc$finished.iter-1) * sigma0 ** 2 + 
                                                  mcmc$finished.iter * (mu0 - mu) ** 2 + 
                                                  (mcmc$meta$tfr_all[year, country] - mu) ** 2)/mcmc$finished.iter)
    
  }
  if (country %in% mcmc$meta3$id_phase3) mcmc$observations[[id3]] <- mcmc$meta$tfr_all[mcmc$meta$lambda_c[country]:mcmc$meta$T_end, country]
  mcmc$data.list[[country]][epsT.idx] <- mcmc$meta$tfr_all[epsT.idx, country]
}

# get.log.lik.year <- function(year.ind, mcmc, Dlpar, phase3par, id_phase1, id_phase2, id_phase3, tfr=NULL, prev=TRUE)
# {
# ## Work for computing log-likelihood from year.ind -> year.ind+1
#   if (is.null(tfr))
#   {
#     tfr_0 <- mcmc$meta$tfr_all[year.ind, ]
#     tfr_1 <- mcmc$meta$tfr_all[year.ind + 1, ]
#   }
#   else if (prev)
#   {
#     tfr_0 <- tfr
#     tfr_1 <- mcmc$meta$tfr_all[year.ind + 1, ]
#   }
#   else
#   {
#     tfr_0 <- mcmc$meta$tfr_all[year.ind, ]
#     tfr_1 <- tfr
#   }
#   
#   eps <- numeric(mcmc$meta$nr_countries)
#   eps[id_phase1] <- tfr_1[id_phase1] - tfr_0[id_phase1]
#   if (!is.null(tfr))
#   {
#     dl <- - Dlpar[id_phase2,5]/(1 + exp(- 2*log(9)/Dlpar[id_phase2, 1] * 
#                                           (tfr_0[id_phase2] - 0.5 * Dlpar[id_phase2, 1] - Dlpar[id_phase2, 2] -
#                                              Dlpar[id_phase2, 3] - Dlpar[id_phase2, 4]))) + 
#       Dlpar[id_phase2, 5]/(1 + exp(- 2*log(9)/Dlpar[id_phase2, 3] *
#                                      (tfr_0[id_phase2] - Dlpar[id_phase2, 4] - 0.5*Dlpar[id_phase2,3])))
#     dl <- dl * ifelse(mcmc$meta$annual.simulation, 1, 5)
#     dl[tfr_0[id_phase2] < 1] <- 0
#     dl[dl<0] <- 0
#     eps[id_phase2] <- tfr_1[id_phase2] - (tfr_0[id_phase2] - dl)
#   }
#   else
#   {
#     eps[id_phase2] <- mcmc$eps_Tc[year.ind, id_phase2]
#   }
#   eps[id_phase3] <- tfr_1[id_phase3] - phase3par[id_phase3, 1] - 
#     phase3par[id_phase3, 2] * (tfr_0[id_phase3] - phase3par[id_phase3, 1])
#   
#   std <- numeric(mcmc$meta$nr_countries)
#   std[id_phase1] <- mcmc$sd_eps_tau
#   std[id_phase3] <- mcmc$sigma.eps
#   if (!is.null(tfr))
#   {
#     const_dummy <- ifelse(mcmc$const_sd_dummie_Tc[year.ind, 1], mcmc$const_sd, 1)
#     tmp_add <- tfr_0[id_phase2] - mcmc$S_sd
#     tmp_add[tfr_0[id_phase2] > mcmc$S_sd] <- - tmp_add[tfr_0[id_phase2] > mcmc$S_sd] * mcmc$a_sd
#     tmp_add[tfr_0[id_phase2] <= mcmc$S_sd] <- tmp_add[tfr_0[id_phase2] <= mcmc$S_sd] * mcmc$b_sd
#     std[id_phase2] <- tmp_add + mcmc$sigma0
#     std[std <= mcmc$meta$sigma0.min] <- mcmc$meta$sigma0.min
#   }
#   else std[id_phase2] <- mcmc$sd_Tc[year.ind, id_phase2]
#   
#   log.lik <- dnorm(eps, mean=0, sd=std, log=TRUE)
#   output <- list(log.lik = log.lik)
#   if (!is.null(tfr))
#   {
#     output$eps <- eps
#     output$std <- std
#   }
#   return (output)
# }

get.log.lik.year <- function(year.ind, mcmc, Dlpar, phase3par, tfr=NULL, prev=TRUE, prev_two=FALSE)
{
  ## Work for computing log-likelihood from year.ind -> year.ind+1
  id_phase1 <- mcmc$meta$id_phase1_by_year[[year.ind]]
  id_phase2 <- mcmc$meta$id_phase2_by_year[[year.ind]]
  id_phase3 <- mcmc$meta$id_phase3_by_year[[year.ind]]
  if (year.ind > 1) id_phase2_prev <- mcmc$meta$id_phase2_by_year[[year.ind-1]]
  
  if (is.null(tfr) || prev_two)
  {
    tfr_0 <- mcmc$meta$tfr_all[year.ind, ]
    tfr_1 <- mcmc$meta$tfr_all[year.ind + 1, ]
  }
  else if (prev)
  {
    tfr_0 <- tfr
    tfr_1 <- mcmc$meta$tfr_all[year.ind + 1, ]
  }
  else
  {
    tfr_0 <- mcmc$meta$tfr_all[year.ind, ]
    tfr_1 <- tfr
  }
  
  eps <- numeric(mcmc$meta$nr_countries)
  eps[id_phase1] <- tfr_1[id_phase1] - tfr_0[id_phase1]
  if (!is.null(tfr))
  {
    dl <- - Dlpar[id_phase2,5]/(1 + exp(- 2*log(9)/Dlpar[id_phase2, 1] * 
                                          (tfr_0[id_phase2] - 0.5 * Dlpar[id_phase2, 1] - Dlpar[id_phase2, 2] -
                                             Dlpar[id_phase2, 3] - Dlpar[id_phase2, 4]))) + 
      Dlpar[id_phase2, 5]/(1 + exp(- 2*log(9)/Dlpar[id_phase2, 3] *
                                     (tfr_0[id_phase2] - Dlpar[id_phase2, 4] - 0.5*Dlpar[id_phase2,3])))
    dl <- dl * ifelse(mcmc$meta$annual.simulation, 1, 5)
    dl[tfr_0[id_phase2] < 1] <- 0
    dl[dl<0] <- 0
    eps[id_phase2] <- tfr_1[id_phase2] - (tfr_0[id_phase2] - dl)
    if ((year.ind > 1) && !(is.null(mcmc$meta$ar.phase2)) && mcmc$meta$ar.phase2)
    {
      joint_phase2 <- intersect(id_phase2, id_phase2_prev)
      if (prev_two) tfr_minus_1 <- tfr
      else tfr_minus_1 <- mcmc$meta$tfr_all[year.ind - 1, ]
      dl_prev <- - Dlpar[joint_phase2,5]/(1 + exp(- 2*log(9)/Dlpar[joint_phase2, 1] * 
                                            (tfr_minus_1[joint_phase2] - 0.5 * Dlpar[joint_phase2, 1] - Dlpar[joint_phase2, 2] -
                                               Dlpar[joint_phase2, 3] - Dlpar[joint_phase2, 4]))) + 
        Dlpar[joint_phase2, 5]/(1 + exp(- 2*log(9)/Dlpar[joint_phase2, 3] *
                                       (tfr_minus_1[joint_phase2] - Dlpar[joint_phase2, 4] - 0.5*Dlpar[joint_phase2,3])))
      dl_prev <- dl_prev * ifelse(mcmc$meta$annual.simulation, 1, 5)
      dl_prev[tfr_minus_1[joint_phase2] < 1] <- 0
      dl_prev[dl_prev<0] <- 0
      eps_prev <- tfr_0[joint_phase2] - (tfr_minus_1[joint_phase2] - dl_prev)
      eps[joint_phase2] <- eps[joint_phase2] - mcmc$rho.phase2 * eps_prev
    }
  }
  else
  {
    eps[id_phase2] <- mcmc$eps_Tc[year.ind, id_phase2]
  }
  eps[id_phase3] <- tfr_1[id_phase3] - phase3par[id_phase3, 1] - 
    phase3par[id_phase3, 2] * (tfr_0[id_phase3] - phase3par[id_phase3, 1])
  
  std <- numeric(mcmc$meta$nr_countries)
  std[id_phase1] <- mcmc$sd_eps_tau
  std[id_phase3] <- mcmc$sigma.eps
  if (!is.null(tfr))
  {
    const_dummy <- ifelse(mcmc$const_sd_dummie_Tc[year.ind, 1], mcmc$const_sd, 1)
    tmp_add <- tfr_0[id_phase2] - mcmc$S_sd
    tmp_add[tfr_0[id_phase2] > mcmc$S_sd] <- - tmp_add[tfr_0[id_phase2] > mcmc$S_sd] * mcmc$a_sd
    tmp_add[tfr_0[id_phase2] <= mcmc$S_sd] <- tmp_add[tfr_0[id_phase2] <= mcmc$S_sd] * mcmc$b_sd
    std[id_phase2] <- tmp_add + mcmc$sigma0
    std[std <= mcmc$meta$sigma0.min] <- mcmc$meta$sigma0.min
  }
  else std[id_phase2] <- mcmc$sd_Tc[year.ind, id_phase2]
  
  log.lik <- dnorm(eps, mean=0, sd=std, log=TRUE)
  output <- list(log.lik = log.lik)
  if (!is.null(tfr))
  {
    output$eps <- eps
    output$std <- std
  }
  return (output)
}

get.log.lik.raw <- function(year.ind, mcmc, tfr)
{
  ## Work for computing log-likelihood from year.ind -> year.ind+1
  if ((length(mcmc$meta$ind.by.year[[year.ind]]) == 0) && mcmc$meta$annual.simulation) return (numeric(mcmc$meta$nr_countries))
  if (!mcmc$meta$annual.simulation && (length(mcmc$meta$ind.by.year[[year.ind]]) == 0 && length(mcmc$meta$ind.by.year[[year.ind + 1]]) == 0)) 
    return (numeric(mcmc$meta$nr_countries))
  output.loglik <- numeric(mcmc$meta$nr_countries)
  
  if (length(mcmc$meta$ind.by.year[[year.ind]]) > 0)
  {
    bias <- mcmc$meta$raw_data.original$bias[mcmc$meta$ind.by.year[[year.ind]]]
    std <- mcmc$meta$raw_data.original$std[mcmc$meta$ind.by.year[[year.ind]]]
    if (mcmc$meta$annual.simulation) tfr.ref <- tfr[mcmc$meta$country.ind.by.year[[year.ind]]]
    else
    {
      tfr.ref <- tfr[mcmc$meta$country.ind.by.year[[year.ind]]] * mcmc$meta$left.distance[[year.ind]] / 5
      tfr.ref <- tfr.ref + mcmc$meta$tfr_matrix_all[max(1, year.ind - 1), mcmc$meta$country.ind.by.year[[year.ind]]] * 
        (5 - mcmc$meta$left.distance[[year.ind]]) / 5
    }
    delta <- mcmc$meta$raw_data.original$DataValue[mcmc$meta$ind.by.year[[year.ind]]] - tfr.ref
    loglik.delta <- dnorm(delta, mean=bias, sd=std, log=TRUE)
    df <- data.frame(idx = mcmc$meta$country.ind.by.year[[year.ind]], loglik = loglik.delta)
    res <- aggregate(df$loglik, by=list(ind=df$idx), sum)
    output.loglik[res$ind] <- res$x
  }
  if (!mcmc$meta$annual.simulation && length(mcmc$meta$ind.by.year[[year.ind + 1]]) > 0)
  {
    bias <- mcmc$meta$raw_data.original$bias[mcmc$meta$ind.by.year[[year.ind + 1]]]
    std <- mcmc$meta$raw_data.original$std[mcmc$meta$ind.by.year[[year.ind + 1]]]
    tfr.ref <- tfr[mcmc$meta$country.ind.by.year[[year.ind + 1]]] * (5 - mcmc$meta$left.distance[[year.ind + 1]]) / 5
    tfr.ref <- tfr.ref + mcmc$meta$tfr_matrix_all[min(mcmc$meta$T_end, year.ind + 1), mcmc$meta$country.ind.by.year[[year.ind + 1]]] * 
      mcmc$meta$left.distance[[year.ind + 1]] / 5
    delta <- mcmc$meta$raw_data.original$DataValue[mcmc$meta$ind.by.year[[year.ind + 1]]] - tfr.ref
    loglik.delta <- dnorm(delta, mean=bias, sd=std, log=TRUE)
    df <- data.frame(idx = mcmc$meta$country.ind.by.year[[year.ind + 1]], loglik = loglik.delta)
    res <- aggregate(df$loglik, by=list(ind=df$idx), sum)
    output.loglik[res$ind] <- output.loglik[res$ind] + res$x
  }
  
  return (output.loglik)
}

mcmc.update.tfr.year <- function(mcmc, countries = NULL)
{
  Dlpar <- cbind((mcmc$U_c - mcmc$Triangle_c4) * 
                   exp(mcmc$gamma_ci)/ apply(exp(mcmc$gamma_ci), 1, sum), 
                 mcmc$Triangle_c4, mcmc$d_c)
  nr_countries <- mcmc$meta$nr_countries
  phase3par <- matrix(nrow = nr_countries, ncol=2)
  phase3par[mcmc$meta$id_phase3,] <- cbind(mcmc$mu.c, mcmc$rho.c)
  ## _r represents phases from year -> year + 1
  id_phase1_r <- mcmc$meta$id_phase1_by_year[[1]]
  id_phase3_r <- mcmc$meta$id_phase3_by_year[[1]]
  id_phase2_r <- mcmc$meta$id_phase2_by_year[[1]]
  if (!is.null(countries))
  {
    id_phase1_r <- intersect(id_phase1_r, countries)
    id_phase3_r <- intersect(id_phase3_r, countries)
    id_phase2_r <- intersect(id_phase2_r, countries)
  }
  for (year in 1:mcmc$meta$T_end)
  {
    if (year > 1)
    {
      id_phase1_l <- id_phase1_r
      id_phase2_l <- id_phase2_r
      id_phase3_l <- id_phase3_r
      id_phase1_r <- mcmc$meta$id_phase1_by_year[[year]]
      id_phase3_r <- mcmc$meta$id_phase3_by_year[[year]]
      id_phase2_r <- mcmc$meta$id_phase2_by_year[[year]]
      if (!is.null(countries))
      {
        id_phase1_r <- intersect(id_phase1_r, countries)
        id_phase3_r <- intersect(id_phase3_r, countries)
        id_phase2_r <- intersect(id_phase2_r, countries)
      }
    }
    loglik_orig <- numeric(nr_countries)
    loglik_proposed <- numeric(nr_countries)
    tfr_orig <- mcmc$meta$tfr_all[year, ]
    if (mcmc$finished.iter <= 100) sd_f <- rep(0.2, nr_countries) else sd_f <- pmax(2.4 * mcmc$meta$tfr_sd_all[year, ], 0.01)
    tfr_proposed <- rep(1, nr_countries)
    if (is.null(countries))
      tfr_proposed <- rnorm(nr_countries, tfr_orig, sd_f)
    else
      tfr_proposed[countries] <- rnorm(length(countries), tfr_orig[countries], sd_f[countries])
    while(any(tfr_proposed < 0.5)) 
    {
      idx <- which(tfr_proposed < 0.5)
      if (mcmc$finished.iter <= 100) sd_f_std <- 0.2 else sd_f_std <- sd_f[idx]
      tfr_proposed[idx] <- rnorm(length(idx), tfr_orig[idx], sd_f_std)
    }
    
    ## Compute Log likelihood
    if ((year < mcmc$meta$T_end - 1) && !is.null(mcmc$meta$ar.phase2) && mcmc$meta$ar.phase2)
    {
      loglik_next <- get.log.lik.year(year + 1, mcmc, Dlpar, phase3par)$log.lik
      loglik_orig <- loglik_orig + loglik_next
      loglik_next_prop <- get.log.lik.year(year + 1, mcmc, Dlpar, phase3par, tfr_proposed, prev_two = TRUE)
      loglik_proposed <- loglik_proposed + loglik_next_prop$log.lik
    }
    if (year < mcmc$meta$T_end)
    {
      loglik_mid <- get.log.lik.year(year, mcmc, Dlpar, phase3par)$log.lik
      loglik_orig <- loglik_orig + loglik_mid
      loglik_mid_prop <- get.log.lik.year(year, mcmc, Dlpar, phase3par, tfr_proposed)
      loglik_proposed <- loglik_proposed + loglik_mid_prop$log.lik
    }
    if (year > 1)
    {
      loglik_orig <- loglik_orig + loglik_prev
      loglik_prev_prop <- get.log.lik.year(year-1, mcmc, Dlpar, phase3par, tfr_proposed, prev = FALSE)
      loglik_proposed <- loglik_proposed + loglik_prev_prop$log.lik
    }
    loglik_orig <- loglik_orig + get.log.lik.raw(year, mcmc, mcmc$meta$tfr_all[year, ])
    loglik_proposed <- loglik_proposed + get.log.lik.raw(year, mcmc, tfr_proposed)
    
    # Accept those with likelihood not dropping so much
    accept.prob <- exp(loglik_proposed - loglik_orig)
    rv_unif <- runif(nr_countries)
    idx_accept <- which(rv_unif < accept.prob)
    mcmc$meta$tfr_all[year, idx_accept] <- tfr_proposed[idx_accept]
    if (year < mcmc$meta$T_end)
    {
      idx_update <- intersect(idx_accept, id_phase2_r)
      mcmc$eps_Tc[year, idx_update] <- loglik_mid_prop$eps[idx_update]
      mcmc$sd_Tc[year, idx_update] <- loglik_mid_prop$std[idx_update]
      tmp <- tfr_proposed[idx_update] - mcmc$S_sd
      idx <- which(tmp < 0)
      tmp[tmp > 0] <- tmp[tmp > 0] * (-mcmc$a_sd)
      tmp[idx] <- tmp[idx] * mcmc$b_sd
      mcmc$add_to_sd_Tc[year, idx_update] <- tmp
    }
    if (year > 1)
    {
      idx_update <- intersect(idx_accept, id_phase2_l)
      mcmc$eps_Tc[year - 1, idx_update] <- loglik_prev_prop$eps[idx_update]
    }
    if ((year < mcmc$meta$T_end - 1) && !is.null(mcmc$meta$ar.phase2) && mcmc$meta$ar.phase2)
    {
      idx_update <- intersect(idx_accept, mcmc$meta$id_phase2_by_year[[year + 1]])
      mcmc$eps_Tc[year + 1, idx_update] <- loglik_next_prop$eps[idx_update]
    }
    
    loglik_prev <- loglik_mid
    loglik_prev[idx_accept] <- loglik_mid_prop$log.lik[idx_accept]
    
    mu0 <- mcmc$meta$tfr_mu_all[year, ]
    sigma0 <- mcmc$meta$tfr_sd_all[year, ]
    mcmc$meta$tfr_mu_all[year, ] <- (mcmc$finished.iter * mu0 + mcmc$meta$tfr_all[year, ]) / (mcmc$finished.iter + 1)
    mu <- mcmc$meta$tfr_mu_all[year, ]
    mcmc$meta$tfr_sd_all[year, ] <- sqrt(((mcmc$finished.iter-1) * sigma0 ** 2 + 
                                            mcmc$finished.iter * (mu0 - mu) ** 2 + 
                                            (mcmc$meta$tfr_all[year, ] - mu) ** 2)/mcmc$finished.iter)
    
  }
  for (country in 1:nr_countries)
  {
    if (!is.null(countries) && !(country %in% countries)) next
    idx <- mcmc$meta$lambda_c[country]:(mcmc$meta$T_end-1)
    if (country %in% mcmc$meta3$id_phase3) 
    {
      id3 <- which(mcmc$meta$id_phase3 == country)  
      mcmc$observations[[id3]] <- mcmc$meta$tfr_all[c(idx, mcmc$meta$T_end), country]
    }
    idx2 <- mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country] - 1)
    mcmc$data.list[[country]][idx2] <- mcmc$meta$tfr_all[idx2, country]
  }
}


one.step.mcmc3.sampling <- function(mcmc)
{
  meta <- mcmc$meta
  niter <- mcmc$iter
  nr.countries <- meta$nr.countries
  countries.index <- meta$id_phase3
  ardata <- mcmc$observations
  Ts <- mcmc$length_obs
  gamma.mu.low <- 1/(meta$sigma.mu.prior.range[2])^2
  gamma.mu.up <- if (meta$sigma.mu.prior.range[1] == 0) NA else 1/(meta$sigma.mu.prior.range[1])^2
  gamma.rho.low <- 1/(meta$sigma.rho.prior.range[2])^2
  gamma.rho.up <- if (meta$sigma.rho.prior.range[1] == 0) NA else 1/(meta$sigma.rho.prior.range[1])^2
  gamma.eps.low <- 1/(meta$sigma.eps.prior.range[2])^2
  gamma.eps.up <- if (meta$sigma.eps.prior.range[1] == 0) NA else 1/(meta$sigma.eps.prior.range[1])^2

  #### Start MCMC

  mu.integral.to.mC <- (mu.rho.integral(mcmc[['mu']], mcmc[['sigma.mu']], low=0))^(-nr.countries)
  prop.mu <- proposal.mu.rho(mcmc[['mu.c']], mcmc[['sigma.mu']], nr.countries,
                             meta[['mu.prior.range']][1], meta[['mu.prior.range']][2])
  accept.prob <- min(((mu.rho.integral(prop.mu, mcmc[['sigma.mu']], low=0))^(-nr.countries))/mu.integral.to.mC, 1)
  if (runif(1) < accept.prob) {
    mcmc[['mu']] <- prop.mu
    mcmc[['recompute.mu.integral']] <- TRUE
  } else mcmc[['recompute.mu.integral']] <- FALSE

  # Metropolis-Hastings for sigma_mu=1/sqrt(lambda_mu)
  S <- sum((mcmc[['mu.c']]-mcmc[['mu']])^2)
  if(mcmc[['recompute.mu.integral']]) mu.integral.to.mC <- mu.rho.integral(mcmc[['mu']], mcmc[['sigma.mu']], low=0)^(-nr.countries)
  prop.lambda.mu <- rgamma.trunc((nr.countries-1)/2, S/2, low=gamma.mu.low, high=gamma.mu.up)
  accept.prob <- min(((mu.rho.integral(mcmc[['mu']], 1/prop.lambda.mu, low=0))^(-nr.countries))/mu.integral.to.mC, 1)
  if (runif(1) < accept.prob) {
    mcmc[['sigma.mu']] <- 1/sqrt(prop.lambda.mu)
    recompute.mu.integral <- TRUE
  } else recompute.mu.integral <- FALSE

  # Slice sampling for rho
  mcmc[['rho']] <- slice.sampling(mcmc[['rho']], logdensity.mu.rho, 1,
                                  low=meta[['rho.prior.range']][1], up=meta[['rho.prior.range']][2],
                                  par.c=mcmc[['rho.c']], sd=mcmc[['sigma.rho']],
                                  c.low=0, c.up=1)
  # Metropolis-Hastings for rho
  # rho.integral.to.mC <- (mu.rho.integral(mcmc[['rho']], mcmc[['sigma.rho']], low=0, up=1))^(-nr.countries)
  # prop.rho <- proposal.mu.rho(mcmc[['rho.c']], mcmc[['sigma.rho']], nr.countries,
  # meta[['rho.prior.range']][1], meta[['rho.prior.range']][2])
  # accept.prob <- min(((mu.rho.integral(prop.rho, mcmc[['sigma.rho']], low=0, up=1))^(-nr.countries))/rho.integral.to.mC, 1)
  # if (runif(1) < accept.prob) {
  # mcmc[['rho']] <- prop.rho
  # recompute.rho.integral <- TRUE
  # } else recompute.rho.integral <- FALSE

  mcmc[['sigma.rho']] <- slice.sampling(mcmc[['sigma.rho']], logdensity.sigma.mu.rho, 1,
                                        low=meta$sigma.rho.prior.range[1], up=meta$sigma.rho.prior.range[2],
                                        par.c=mcmc[['rho.c']], mean=mcmc[['rho']],
                                        c.low=0, c.up=1)

  # Metropolis-Hastings for sigma_rho=1/sqrt(lambda_rho)
  # S <- sum((mcmc[['rho.c']]-mcmc[['rho']])^2)
  # if(recompute.rho.integral)
  # rho.integral.to.mC <- (mu.rho.integral(mcmc[['rho']], mcmc[['sigma.rho']], low=0, up=1))^(-nr.countries)
  # prop.lambda.rho <- rgamma.trunc((nr.countries-1)/2, S/2, low=gamma.rho.low, high=gamma.rho.up)
  # accept.prob <- min(((mu.rho.integral(mcmc[['rho']], 1/prop.lambda.rho, low=0, up=1))^(-nr.countries))/rho.integral.to.mC, 1)
  # if (runif(1) < accept.prob) {
  # mcmc[['sigma.rho']] <- 1/sqrt(prop.lambda.rho)
  # recompute.rho.integral <- TRUE
  # } else recompute.rho.integral <- FALSE

  sigma.eps.sq <- mcmc$sigma.eps^2
  sigma.mu.sq <- mcmc$sigma.mu^2
  sigma.mu.sq.inv <- 1/sigma.mu.sq
  mu.over.sigma.sq <- mcmc$mu/sigma.mu.sq
  sigma.rho.sq <- mcmc$sigma.rho^2
  sigma.rho.sq.inv <- 1/sigma.rho.sq
  rho.over.sigma.sq <- mcmc$rho/sigma.rho.sq
  S.eps <- STn <- 0
  one.minus.rho <- 1-mcmc$rho.c
  one.minus.rho.sq <- one.minus.rho^2
  W <- one.minus.rho.sq/sigma.eps.sq

  # country-specific parameters - Gibbs sampler
  for(country in 1:nr.countries) {
    f.ct <- ardata[[country]][2:Ts[country]]
    f.ctm1 <- ardata[[country]][1:(Ts[country]-1)]
    # mu.c
    s <- sum((f.ct - mcmc$rho.c[country]*f.ctm1)/one.minus.rho[country])
    nomin <- W[country] * s + mu.over.sigma.sq
    denom <- (Ts[country]-1) * W[country] + sigma.mu.sq.inv
    mcmc$mu.c[country] <- rnorm.trunc(mean=nomin/denom, sd=1/sqrt(denom), low=0)

    # rho.c
    d1 <- f.ctm1 - mcmc$mu.c[country]
    a <- sum(d1^2)/sigma.eps.sq
    b <- sum(d1*(f.ct - mcmc$mu.c[country]))/sigma.eps.sq
    nomin <- b + rho.over.sigma.sq
    denom <- a + sigma.rho.sq.inv
    mcmc$rho.c[country] <- rnorm.trunc(mean=nomin/denom, sd=1/sqrt(denom), low=0, #high=1-10*.Machine$double.xmin
                                       high=0.999999)
    S.eps <- S.eps + sum((f.ct - mcmc$mu.c[country] - mcmc$rho.c[country]*d1)^2)
    STn <- STn + Ts[country]-1
  }
  # Gibbs for sigma.eps
  mcmc$sigma.eps <- 1/sqrt(rgamma.trunc((STn-1)/2, S.eps/2, low=gamma.eps.low, high=gamma.eps.up))

}

