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
    gamma.mu.low <- 1/(meta$sigma.mu.prior.range[2])^2
    gamma.mu.up <- if (meta$sigma.mu.prior.range[1] == 0) NA else 1/(meta$sigma.mu.prior.range[1])^2
    gamma.rho.low <- 1/(meta$sigma.rho.prior.range[2])^2
    gamma.rho.up <- if (meta$sigma.rho.prior.range[1] == 0) NA else 1/(meta$sigma.rho.prior.range[1])^2
    gamma.eps.low <- 1/(meta$sigma.eps.prior.range[2])^2
    gamma.eps.up <- if (meta$sigma.eps.prior.range[1] == 0) NA else 1/(meta$sigma.eps.prior.range[1])2
    recompute.mu.integral <- TRUE
    recompute.rho.integral <- TRUE
    
    # Start MCMC
	############
    for (iter in start.iter:niter) {
    	if(verbose.iter > 0 && (iter %% verbose.iter == 0))
        	cat('\nIteration:', iter, '--', date())
        unblock.gtk('bDem.TFRmcmc')
		# Metropolis-Hastings for mu
		mu.integral.to.mC <- (mu.rho.integral(mcmc[['mu']], mcmc[['sigma.mu']], low=0))^(-nr.countries)
		prop.mu <- proposal.mu.rho(mcmc[['mu.c']], mcmc[['sigma.mu']], nr.countries, 
						meta[['mu.prior.range']][1], meta[['mu.prior.range']][2])
		accept.prob <- min(((mu.rho.integral(prop.mu, mcmc[['sigma.mu']], low=0))^(-nr.countries))/mu.integral.to.mC, 1)
		if (runif(1) < accept.prob) {
			mcmc[['mu']] <- prop.mu
			recompute.mu.integral <- TRUE
		} else recompute.mu.integral <- FALSE
		
		# Metropolis-Hastings for sigma_mu=1/sqrt(lambda_mu)
		S <- sum((mcmc[['mu.c']]-mcmc[['mu']])^2)
		if(recompute.mu.integral) mu.integral.to.mC <- mu.rho.integral(mcmc[['mu']], mcmc[['sigma.mu']], low=0))^(-nr.countries)
		prop.lambda.mu <- rgamma.trunc((nr.countries-1)/2, S/2, low=gamma.mu.low, up=gamma.mu.up)
		accept.prob <- min(((mu.rho.integral(mcmc[['mu']], 1/prop.lambda.mu, low=0))^(-nr.countries))/mu.integral.to.mC, 1)
		if (runif(1) < accept.prob) {
			mcmc[['sigma.mu']] <- 1/sqrt(prop.lambda.mu)
			recompute.mu.integral <- TRUE
		} else recompute.mu.integral <- FALSE	
		
		# Metropolis-Hastings for rho
		rho.integral.to.mC <- (mu.rho.integral(mcmc[['rho']], mcmc[['sigma.rho']], low=0, up=1))^(-nr.countries)
		prop.rho <- proposal.mu.rho(mcmc[['rho.c']], mcmc[['sigma.rho']], nr.countries, 
						meta[['rho.prior.range']][1], meta[['rho.prior.range']][2])
		accept.prob <- min(((mu.rho.integral(prop.rho, mcmc[['sigma.rho']], low=0, up=1))^(-nr.countries))/rho.integral.to.mC, 1)
		if (runif(1) < accept.prob) {
			mcmc[['rho']] <- prop.rho
			recompute.rho.integral <- TRUE
		} else recompute.rho.integral <- FALSE
		
		# Metropolis-Hastings for sigma_rho=1/sqrt(lambda_rho)
		S <- sum((mcmc[['rho.c']]-mcmc[['rho']])^2)
		if(recompute.rho.integral) 
			rho.integral.to.mC <- (mu.rho.integral(mcmc[['rho']], mcmc[['sigma.rho']], low=0, up=1))^(-nr.countries)
		prop.lambda.rho <- rgamma.trunc((nr.countries-1)/2, S/2, low=gamma.rho.low, up=gamma.rho.up)
		accept.prob <- min(((mu.rho.integral(mcmc[['rho']], 1/prop.lambda.rho, low=0, up=1))^(-nr.countries))/rho.integral.to.mC, 1)
		if (runif(1) < accept.prob) {
			mcmc[['sigma.rho']] <- 1/sqrt(prop.lambda.rho)
			recompute.rho.integral <- TRUE
		} else recompute.rho.integral <- FALSE
		
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
			mcmc$rho.c[country] <- rnorm.trunc(mean=nomin/denom, sd=1/sqrt(denom), low=0, high=1-.Machine$double.xmin)
			S.eps <- S.eps + sum((f.ct - mcmc$mu.c[country] - mcmc$rho.c[country]*d1)^2)
			STn <- STn + Ts[country]-1
		}
		# Gibbs for sigma.eps
		mcmc$sigma.eps <- 1/sqrt(rgamma.trunc((STn-1)/2, S.eps/2, low=gamma.eps.low, high=gamma.eps.up))
		
		mcmc$finished.iter <- mcmc$finished.iter+1
        mcmc$rng.state <- .Random.seed  
        if (iter %% thin == 0){
        	mcmc$length <- mcmc$length + 1
        	flush.buffer <- FALSE
            if (iter + 1 > niter) flush.buffer <- TRUE                
            store.mcmc3(mcmc, append=TRUE, flush.buffer=flush.buffer, verbose=verbose)
         }
	}
	return(mcmc)
}

proposal.mu.rho <- function(par.c, sigma, C, low, up=NA)
 return(rnorm.trunc(mean=sum(par.c)/C, sd=sigma/sqrt(C), low=low, high=up))
 
mu.rho.integral <- function(par, sigma, low=0, up=NA) {
	term1 <- 1
	if(!is.na(up)) term1 <- pnorm((up-par)/sigma)
	return(term1 - pnorm((low-par)/sigma))
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

