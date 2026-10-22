#########################################################
# MCMC sampling for Phase III
#########################################################

tfr3.mcmc.sampling <- function(mcmc, thin=1, start.iter=2, verbose=FALSE, verbose.iter=10) {
  if (!is.null(mcmc$rng.state)) .Random.seed <- mcmc$rng.state
	meta <- mcmc$meta
  niter <- mcmc$iter
  nr.countries <- meta$nr.countries
  countries.index <- meta$id_phase3
  ardata <- list()
  Ts <- rep(0, nr.countries)
  for(country in 1:nr.countries) {
    data <- get.observed.tfr(countries.index[country], meta$parent, 'tfr_matrix_all')
		ardata[[country]] <- data[meta$parent$lambda_c[countries.index[country]]:meta$parent$T_end_c[countries.index[country]]]
		Ts[country] <- length(ardata[[country]])
  }
  mcmc$observations <- ardata
  mcmc$length_obs <- Ts
  mcmc$recompute.mu.integral <- TRUE
  mcmc$recompute.rho.integral <- TRUE
  mcenv <- as.environment(mcmc)
    
  # Start MCMC
	############
  for (iter in start.iter:niter) {
    if(verbose.iter > 0 && (iter %% verbose.iter == 0))
      cat('\nIteration:', iter, '--', date())
    unblock.gtk('bDem.TFRmcmc')
		# Metropolis-Hastings for mu
		one.step.mcmc3.sampling(mcenv)
		mcenv$finished.iter <- mcenv$finished.iter + 1
		mcenv$rng.state <- .Random.seed
    if (iter %% thin == 0){
      mcenv$length <- mcenv$length + 1
      flush.buffer <- FALSE
      if (iter + thin > niter) flush.buffer <- TRUE                
      store.mcmc3(mcenv, append=TRUE, flush.buffer=flush.buffer, verbose=verbose)
    }
  }
  .cleanup.mcmc(mcenv)
  resmc <- as.list(mcenv)
  class(resmc) <- class(mcmc)
  return(resmc)
}

one.step.mcmc3.sampling <- function(mcmc)
{
  meta <- mcmc$meta
  nr.countries <- meta$nr.countries
  gamma.mu.low <- 1/(meta$sigma.mu.prior.range[2])^2
  gamma.mu.up <- if (meta$sigma.mu.prior.range[1] == 0) NA else 1/(meta$sigma.mu.prior.range[1])^2
  gamma.eps.low <- 1/(meta$sigma.eps.prior.range[2])^2
  gamma.eps.up <- if (meta$sigma.eps.prior.range[1] == 0) NA else 1/(meta$sigma.eps.prior.range[1])^2

  #### Start MCMC

  # Metropolis-Hastings for mu
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
  }

  # Slice sampling for the world-level AR parameters and their standard deviations:
  # rho for AR(1); partial autocorrelations alpha1 in [0,1) and alpha2 in [-1,1) for AR(2)
  ar.pars <- if(is.ar2.phase3(meta)) c('alpha1', 'alpha2') else 'rho'
  ar.pars.low <- if(is.ar2.phase3(meta)) c(0, -1) else 0
  for(i in seq_along(ar.pars)) {
    par <- ar.pars[i]
    sigma.par <- paste0('sigma.', par)
    par.c <- paste0(par, '.c')
    mcmc[[par]] <- slice.sampling(mcmc[[par]], logdensity.mu.rho, 1,
                                  low=meta[[paste0(par, '.prior.range')]][1], up=meta[[paste0(par, '.prior.range')]][2],
                                  par.c=mcmc[[par.c]], sd=mcmc[[sigma.par]],
                                  c.low=ar.pars.low[i], c.up=1)
    mcmc[[sigma.par]] <- slice.sampling(mcmc[[sigma.par]], logdensity.sigma.mu.rho, 1,
                                        low=meta[[paste0(sigma.par, '.prior.range')]][1], up=meta[[paste0(sigma.par, '.prior.range')]][2],
                                        par.c=mcmc[[par.c]], mean=mcmc[[par]],
                                        c.low=ar.pars.low[i], c.up=1)
  }

  # country-specific parameters - Gibbs sampler
  res <- mcmc3.update.cs.pars(mcmc)

  # Gibbs for sigma.eps
  mcmc$sigma.eps <- 1/sqrt(rgamma.trunc((res$STn-1)/2, res$S.eps/2, low=gamma.eps.low, high=gamma.eps.up))
}

mcmc3.update.cs.pars <- function(mcmc) {
  # Gibbs sampler for the country-specific Phase III parameters, i.e. mu.c and rho.c for AR(1),
  # or mu.c, alpha1.c, alpha2.c for AR(2). mcmc is expected to be an environment.
  # Returns the sum of squared residuals and their count, needed for updating sigma.eps.
  ardata <- mcmc$observations
  Ts <- mcmc$length_obs
  sigma.eps.sq <- mcmc$sigma.eps^2
  sigma.mu.sq <- mcmc$sigma.mu^2
  sigma.mu.sq.inv <- 1/sigma.mu.sq
  mu.over.sigma.sq <- mcmc$mu/sigma.mu.sq
  S.eps <- STn <- 0

  if(!is.ar2.phase3(mcmc$meta)) {
    sigma.rho.sq <- mcmc$sigma.rho^2
    sigma.rho.sq.inv <- 1/sigma.rho.sq
    rho.over.sigma.sq <- mcmc$rho/sigma.rho.sq
    one.minus.rho <- 1-mcmc$rho.c
    one.minus.rho.sq <- one.minus.rho^2
    W <- one.minus.rho.sq/sigma.eps.sq
    for(country in seq_len(mcmc$meta$nr.countries)) {
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
    return(list(S.eps=S.eps, STn=STn))
  }
  # AR(2): f_t - mu = rho1 (f_{t-1} - mu) + rho2 (f_{t-2} - mu) + eps_t,
  # with rho1 = alpha1 (1 - alpha2), rho2 = alpha2; conditioned on the first two observations
  sigma.alpha1.sq <- mcmc$sigma.alpha1^2
  sigma.alpha1.sq.inv <- 1/sigma.alpha1.sq
  alpha1.over.sigma.sq <- mcmc$alpha1/sigma.alpha1.sq
  sigma.alpha2.sq <- mcmc$sigma.alpha2^2
  sigma.alpha2.sq.inv <- 1/sigma.alpha2.sq
  alpha2.over.sigma.sq <- mcmc$alpha2/sigma.alpha2.sq
  for(country in seq_len(mcmc$meta$nr.countries)) {
    f.ct <- ardata[[country]][3:Ts[country]]
    f.ctm1 <- ardata[[country]][2:(Ts[country]-1)]
    f.ctm2 <- ardata[[country]][1:(Ts[country]-2)]
    alpha1 <- mcmc$alpha1.c[country]
    alpha2 <- mcmc$alpha2.c[country]
    rho1 <- alpha1 * (1 - alpha2)
    # mu.c: f_t - rho1 f_{t-1} - rho2 f_{t-2} = (1 - rho1 - rho2) mu + eps_t
    one.minus.rho <- 1 - rho1 - alpha2
    nomin <- one.minus.rho * sum(f.ct - rho1*f.ctm1 - alpha2*f.ctm2)/sigma.eps.sq + mu.over.sigma.sq
    denom <- (Ts[country]-2) * one.minus.rho^2/sigma.eps.sq + sigma.mu.sq.inv
    mcmc$mu.c[country] <- rnorm.trunc(mean=nomin/denom, sd=1/sqrt(denom), low=0)

    d0 <- f.ct - mcmc$mu.c[country]
    d1 <- f.ctm1 - mcmc$mu.c[country]
    d2 <- f.ctm2 - mcmc$mu.c[country]
    # alpha1.c: d0 - alpha2 d2 = alpha1 (1 - alpha2) d1 + eps
    u <- (1 - alpha2) * d1
    nomin <- sum(u * (d0 - alpha2*d2))/sigma.eps.sq + alpha1.over.sigma.sq
    denom <- sum(u^2)/sigma.eps.sq + sigma.alpha1.sq.inv
    alpha1 <- mcmc$alpha1.c[country] <- rnorm.trunc(mean=nomin/denom, sd=1/sqrt(denom), low=0, high=0.999999)

    # alpha2.c: d0 - alpha1 d1 = alpha2 (d2 - alpha1 d1) + eps
    v <- d2 - alpha1 * d1
    nomin <- sum(v * (d0 - alpha1*d1))/sigma.eps.sq + alpha2.over.sigma.sq
    denom <- sum(v^2)/sigma.eps.sq + sigma.alpha2.sq.inv
    alpha2 <- mcmc$alpha2.c[country] <- rnorm.trunc(mean=nomin/denom, sd=1/sqrt(denom), low=-0.999999, high=0.999999)

    S.eps <- S.eps + sum((d0 - alpha1*(1 - alpha2)*d1 - alpha2*d2)^2)
    STn <- STn + Ts[country]-2
  }
  return(list(S.eps=S.eps, STn=STn))
}

proposal.mu.rho <- function(par.c, sigma, C, low, up=NA)
 return(rnorm.trunc(mean=sum(par.c)/C, sd=sigma/sqrt(C), low=low, high=up))
 
mu.rho.integral <- function(par, sigma, low=0, up=NA) {
	term1 <- 1
	if(!is.na(up)) term1 <- pnorm((up-par)/sigma)
	return(term1 - pnorm((low-par)/sigma))
}

logdensity.mu.rho <- function(par, low, up, par.c, sd, c.low, c.up) {
	return(log(max(dunif(par, low, up), .Machine$double.xmin)) + sum(log(pmax(dnorm.trunc(par.c, par, sd, c.low, c.up),
				.Machine$double.xmin))))
}

logdensity.sigma.mu.rho <- function(par, low, up, par.c, mean, c.low, c.up) {
	return(log(max(dunif(par, low, up), .Machine$double.xmin)) + sum(log(pmax(dnorm.trunc(par.c, mean, par, c.low, c.up),
				.Machine$double.xmin))))
}



.get.trunc.condition <- function(temp, i, maxit, low, high)
	return((temp<low || temp>high) && i <= maxit)

.get.ltrunc.condition <- function(temp, i, maxit, low, ...)
	return(temp<low && i <= maxit)
	
rnorm.trunc<-function(mean,sd,low,high=NA){
	temp<--999
	maxit <- 10
	i <- 1
	cond.fct <- if(is.na(high)) '.get.ltrunc.condition' else '.get.trunc.condition'
	while(do.call(cond.fct, list(temp, i, maxit, low, high))) {
		temp<-rnorm(1,mean=mean,sd=sd)
		i <- i+1
	}
	if (i > maxit) {
		temp <- if(temp<low) low else {if(!is.na(high)) high else temp}
		#warning(paste('Maximum iterations reached in rnorm.trunc(', 
		#			mean, ',', sd, '). Value truncated to ', temp, '.', sep=''), immediate.=TRUE)
	}
	return(temp)
}

dnorm.trunc<-function(x,mean,sd,low,high){
  out<-dnorm(x,mean=mean,sd=sd)/(pnorm(high,mean=mean,sd=sd)-pnorm(low,mean=mean,sd=sd))
  out[x<low]<-0
  out[x>high]<-0
  return(out)
}

rgamma.trunc<-function(shape,rate,low,high=NA){
  temp<--999
  maxit <- 10
  i <- 1
  cond.fct <- if(is.na(high)) '.get.ltrunc.condition' else '.get.trunc.condition'
  while(do.call(cond.fct, list(temp, i, maxit, low, high))) {
     temp<-rgamma(1,shape=shape,rate=rate)
     i <- i+1
  }
  if (i > maxit) {
  	temp <- if(temp<low) low else {if(!is.na(high)) high else temp}
  	#warning(paste('Maximum iterations reached in rgamma.trunc(', shape, ',', rate, ').', sep=''), immediate.=TRUE)
  }
  return(temp)
}

slice.sampling <- function(x0, fun, width,  ..., low, up, maxit=50) {
	# Slightly modified version of 
	# http://www.cs.toronto.edu/~radford/ftp/slice-R-prog (Radford M. Neal, 17 March 2008)
	gx0 <- fun(x0, ..., low=low, up=up)
	z <- gx0 - rexp(1) # slice S={x: z < gx0}
	L <- x0 - runif(1, 0, width)
	R <- L + width # should guarantee that x0 is in [L,R], even with roundoff
	#print(c(L,R,z))
	# Expand the interval until its ends are outside the slice, or until
	# the limit on steps is reached.
	J <- floor(runif(1,0,maxit))
    K <- (maxit-1) - J
    while (J>0 && L > low && fun(L,  ..., low=low, up=up)>z) {
      L <- L - width
      J <- J - 1
    }
    while (K>0 && R < up && fun(R,  ..., low=low, up=up)>z) {
      R <- R + width
      K <- K - 1
    }
    #print(c(maxit - K - J, z, L, R))
	# Shrink interval to lower and upper bounds.
	if (L<low) L <- low
  	if (R>up) R <- up
 	#if(debug) print(c('Slice sampling begin:', L, R, z, x0))
	# Sample from the interval, shrinking it on each rejection.
	i<-1
	while(i<=maxit) {
		x1 <- runif(1,L,R)
		if(z <= fun(x1,  ..., low=low, up=up)) {
			#if(debug) print(c('Slice sampling end:', L, R, x1))
			return(x1)
		}
		if (x1 < x0) L <- x1
		else R <- x1
		i <- i+1
	}
	stop('Problem in slice sampling')
}
