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
