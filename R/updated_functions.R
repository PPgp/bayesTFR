require(mvtnorm)

get.obs.estimate.diff <- function(mcmc)
{
  mcmc$eps_unc <- list()
  for(country in 1:mcmc$meta$nr_countries)
  {
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

update.obs.estimate.diff <- function(mcmc, country, year)
{
  if (mcmc$meta$annual.simulation)
  {
    indices <- which(round(mcmc$meta$raw.data[[country]]$Year - mcmc$meta$start.year + 1 + 0.01) == year)
  }
  else
  {
    indices <- which((mcmc$meta$raw.data[[country]]$Year < (year * 5 + mcmc$meta$start.year + 3)) & 
                       (mcmc$meta$raw.data[[country]]$Year > (year * 5 + mcmc$meta$start.year -7)))
  }
  if (length(indices) > 0)
  {
    tmp.tfr <- mcmc$meta$raw.data[[country]]$DataValue[indices]
    tmp.year <- mcmc$meta$raw.data[[country]]$Year[indices]
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
  for(country in 1:mcmc$meta$nr_countries)
  {
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
    std[std < 1e-6] <- 0.5
    mcmc$meta$raw.data[[country]]$bias <- bias
    mcmc$meta$raw.data[[country]]$std <- std
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
  
  for (year in 1:mcmc$meta$T_end)
  {
    # if ((mcmc$finished.iter == 146) && (country == 113) && (year == 6)) {browser()}
    tmp.years <- max(1, year-1):min(year, mcmc$meta$T_end -1)
    sd_f <- ifelse(mcmc$finished.iter>=100, 2.4 * max(mcmc$meta$tfr_sd_all[year, country], 0.01), 1)
    prop_tfr <- rnorm(1, mcmc$meta$tfr_all[year, country], sd_f)
    while (prop_tfr < 0.5) prop_tfr <- rnorm(1, mcmc$meta$tfr_all[year, country], sd_f)
    if (mcmc$meta$annual.simulation)
    {
      indices <- which(round(mcmc$meta$raw.data[[country]]$Year - mcmc$meta$start.year + 1.01) == year)
      # indices <- which((mcmc$meta$raw.data[[country]]$Year < (year + mcmc$meta$start.year)) & 
      #                    (mcmc$meta$raw.data[[country]]$Year > (year + mcmc$meta$start.year -2)))
    }
    else
    {
      indices <- which((mcmc$meta$raw.data[[country]]$Year < (year * 5 + mcmc$meta$start.year + 3)) & 
                         (mcmc$meta$raw.data[[country]]$Year > (year * 5 + mcmc$meta$start.year - 7)))
    }
    tfr_old <- mcmc$meta$tfr_all[year, country]
    eps_prev <- mcmc$eps_unc[[country]][indices]
    loglik_prev <- sum(dnorm(eps_prev, mean=mcmc$meta$raw.data[[country]]$bias[indices], sd=mcmc$meta$raw.data[[country]]$std[indices], log=TRUE))
    loglik_prev <- loglik_prev + sum(dnorm(eps_tfr_prev[tmp.years], mean =0, sd = sd_tfr_prev[tmp.years], log=TRUE))
    mcmc$meta$tfr_all[year, country] <- prop_tfr
    eps_new <- update.obs.estimate.diff(mcmc, country, year)
    loglik_new <- sum(dnorm(eps_new, mean=mcmc$meta$raw.data[[country]]$bias[indices], sd=mcmc$meta$raw.data[[country]]$std[indices], log=TRUE))
    if (year < mcmc$meta$start_c[country])
    {
      eps_tfr_prop[tmp.years] <- mcmc$meta$tfr_all[tmp.years + 1, country] - mcmc$meta$tfr_all[tmp.years, country]
    }
    else if (year > mcmc$meta$lambda_c[country])
    {
      id3 <- which(mcmc$meta$id_phase3 == country)
      eps_tfr_prop[tmp.years] <- mcmc$meta$tfr_all[tmp.years + 1, country] - mcmc$meta$tfr_all[tmp.years, country] * 
        mcmc$rho.c[id3] - (1-mcmc$rho.c[id3]) * mcmc$mu.c[id3]
    }
    else 
      eps_tfr_prop <- get.eps.all.phases(Dlpar, mcmc, country)
    
    if ((year >= mcmc$meta$start_c[country]) && (year < mcmc$meta$lambda_c[country]))
    {
      tmp <- (prop_tfr - mcmc$S_sd) * ifelse(prop_tfr > mcmc$S_sd, -mcmc$a_sd, mcmc$b_sd)
      sd_tfr_prop[year] <- ifelse(mcmc$const_sd_dummie_Tc[year] == 1, mcmc$const_sd, 1) * ifelse((mcmc$sigma0 + tmp > 0), mcmc$sigma0 + tmp, mcmc$meta$sigma0.min)
    }
      
    loglik_new <- loglik_new + sum(dnorm(eps_tfr_prop[tmp.years], mean=0, sd=sd_tfr_prop[tmp.years], log=TRUE))
    accept.prob <- exp(loglik_new - loglik_prev)
    if (runif(1) < accept.prob)
    {
      if (year %in% epsT.idx)
      {
        mcmc$eps_Tc[year, country] <- eps_tfr_prop[year]
        mcmc$sd_Tc[year, country] <- sd_tfr_prop[year]
        mcmc$add_to_sd_Tc[year, country] <- (prop_tfr - mcmc$S_sd) * ifelse(prop_tfr > mcmc$S_sd, -mcmc$a_sd, mcmc$b_sd)
      }
      if ((year - 1) %in% epsT.idx)
      {
        mcmc$eps_Tc[year - 1, country] <- eps_tfr_prop[year - 1]
        mcmc$sd_Tc[year - 1, country] <- sd_tfr_prop[year - 1]
      }
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
  mcmc$data.list[[country]][epsT.idx] <- mcmc$meta$tfr_all[epsT.idx, country]
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

