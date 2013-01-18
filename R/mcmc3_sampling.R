#########################################################
# MCMC sampling for Phase III
#########################################################

tfr3.mcmc.sampling <- function(mcmc, thin=1, start.iter=2, verbose=FALSE, verbose.iter=10) {
	library(msm)
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
    # par.c <- par.range <- par.sig <- par.sig.range <- list()
    # for(par in c('mu', 'rho')) {
    	# par.c[[par]] <- paste(par, 'c', sep='.')
    	# par.range[[par]] <- paste(par, 'prior', 'range', sep='.')
    	# par.sig[[par]] <- paste('sigma', par, sep='.')
    	# par.sig.range[[par]] <- paste(par.sig[[par]], 'prior', 'range', sep='.')
    # }
    mcmc$observations <- ardata
      
    # Start MCMC
	############
    for (iter in start.iter:niter) {
    	if(verbose.iter > 0 && (iter %% verbose.iter == 0))
        	cat('\nIteration:', iter, '--', date())
        unblock.gtk('bDem.TFRmcmc')

		# world parameters mu, rho, sigma.mu, sigma.rho
		# for(par in c('mu', 'rho')) {
			# sum.par <- sum(mcmc[[par.c[[par]]]])
			# par.mean <- sum.par/nr.countries
			# par.sd <- mcmc[[par.sig[[par]]]]/sqrt(nr.countries)
			# mcmc[[par]] <- rtnorm(1, mean=par.mean, sd=par.sd, lower=meta[[par.range[[par]]]][1], upper=meta[[par.range[[par]]]][2])
			# S <- sum((mcmc[[par.c[[par]]]]-mcmc[[par]])^2)
			# mcmc[[par.sig[[par]]]] <- 1/sqrt(rgamma.ltrunc((nr.countries-1)/2, S/2, low=1/meta[[par.sig.range[[par]]]][2]^2))
			# #mcmc[[par.sig]] <- slice.sampling(mcmc[[par.sig]], logdensity.sigma.mu.rho, 0.5, 
			# #					par=mcmc[[par]], par.c=mcmc[[par.c]], 
			# #					low=meta[[par.sig.range]][1], up=meta[[par.sig.range]][2])
		# }
		prop.mu <- proposal.mu.rho(mcmc[['mu.c']], mcmc[['sigma.mu']], nr.countries, 
						meta[['mu.prior.range']][1], meta[['mu.prior.range']][2])
		accept.prob <- min(((mu.rho.integral(prop.mu, mcmc[['sigma.mu']], low=0))^(-nr.countries))/(
							(mu.rho.integral(mcmc[['mu']], mcmc[['sigma.mu']], low=0))^(-nr.countries)), 1)
		if (runif(1) < accept.prob) mcmc[['mu']] <- prop.mu
		
		S <- sum((mcmc[['mu.c']]-mcmc[['mu']])^2)
		prop.lambda.mu <- rgamma.ltrunc((nr.countries-1)/2, S/2, low=1/meta[['sigma.mu.prior.range']][2]^2)
		accept.prob <- min(((mu.rho.integral(mcmc[['mu']], 1/prop.lambda.mu, low=0))^(-nr.countries))/(
						(mu.rho.integral(mcmc[['mu']], mcmc[['sigma.mu']], low=0))^(-nr.countries)), 1)
		if (runif(1) < accept.prob) mcmc[['sigma.mu']] <- 1/sqrt(prop.lambda.mu)
		
		prop.rho <- proposal.mu.rho(mcmc[['rho.c']], mcmc[['sigma.rho']], nr.countries, 
						meta[['rho.prior.range']][1], meta[['rho.prior.range']][2])
		accept.prob <- min(((mu.rho.integral(prop.rho, mcmc[['sigma.rho']], low=0, up=1))^(-nr.countries))/(
							(mu.rho.integral(mcmc[['rho']], mcmc[['sigma.rho']], low=0, up=1))^(-nr.countries)), 1)
		if (runif(1) < accept.prob) mcmc[['rho']] <- prop.rho
		
		S <- sum((mcmc[['rho.c']]-mcmc[['rho']])^2)
		prop.lambda.rho <- rgamma.ltrunc((nr.countries-1)/2, S/2, low=1/meta[['sigma.rho.prior.range']][2]^2)
		accept.prob <- min(((mu.rho.integral(mcmc[['rho']], 1/prop.lambda.rho, low=0, up=1))^(-nr.countries))/(
							(mu.rho.integral(mcmc[['rho']], mcmc[['sigma.rho']], low=0, up=1))^(-nr.countries)), 1)
		if (runif(1) < accept.prob) mcmc[['sigma.rho']] <- 1/sqrt(prop.lambda.rho)
		
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
		# country-specific parameters
		for(country in 1:nr.countries) {
			f.ct <- ardata[[country]][2:Ts[country]]
			f.ctm1 <- ardata[[country]][1:(Ts[country]-1)]
			# mu.c
			s <- sum((f.ct - mcmc$rho.c[country]*f.ctm1)/one.minus.rho[country])
			nomin <- W[country] * s + mu.over.sigma.sq
			denom <- (Ts[country]-1) * W[country] + sigma.mu.sq.inv
			mcmc$mu.c[country] <- rtnorm(1, mean=nomin/denom, sd=1/sqrt(denom), lower=0)	

			# rho.c
			d1 <- f.ctm1 - mcmc$mu.c[country]
			a <- sum(d1^2)/sigma.eps.sq
			b <- sum(d1*(f.ct - mcmc$mu.c[country]))/sigma.eps.sq
			nomin <- b + rho.over.sigma.sq
			denom <- a + sigma.rho.sq.inv
			mcmc$rho.c[country] <- rtnorm(1, mean=nomin/denom, sd=1/sqrt(denom), lower=0, upper=1-.Machine$double.xmin)

			# mcmc$mu.c[country] <- slice.sampling(mcmc$mu.c[country],
							# logdensity.mu.rho.c, 1, mean=mcmc$mu, sd=mcmc$sigma.mu, sigma.eps=mcmc$sigma.eps, 
							# rho.c=mcmc$rho.c[country], data=ardata[[country]], low=0, up=999999
							# )
			#mcmc$rho.c[country] <- slice.sampling(mcmc$rho.c[country],
			#				 logdensity.rho.c, 1, rho=mcmc$rho, sigma.rho=mcmc$sigma.rho, sigma.eps=mcmc$sigma.eps, 
			#				 mu.c=mcmc$mu.c[country], data=ardata[[country]], low=0, up=1, Ts=Ts[country]
			#				 )
			S.eps <- S.eps + sum((f.ct - mcmc$mu.c[country] - mcmc$rho.c[country]*d1)^2)
			STn <- STn + Ts[country]-1
		}
		# mcmc$sigma.eps <- slice.sampling(mcmc$sigma.eps, logdensity.sigma.eps, 0.3, 
								# data=ardata, mu.c=mcmc$mu.c, rho.c=mcmc$rho.c, 
								# low=meta$sigma.eps.prior.range[1], up=meta$sigma.eps.prior.range[2],
								# C=nr.countries, Ts=Ts)
		mcmc$sigma.eps <- 1/sqrt(rgamma.ltrunc((STn-1)/2, S.eps/2, low=1/meta$sigma.eps.prior.range[2]^2))
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
 
logdensity.rho.c <- function(x, rho, sigma.rho, mu.c, sigma.eps, data, low, up, Ts) {
	res <- sum(log(pmax(dnorm((data[2:Ts]-mu.c)/(data[1:(Ts-1)]-mu.c), x, sd=sigma.eps/abs(data[1:(Ts-1)]-mu.c)),.Machine$double.xmin)))
	return(log(dnorm.trunc(x, rho, sigma.rho, 0, 1)) + res)
}

logdensity.mu.rho.c <- function(x, mean, sd, sigma.eps, rho.c, data, low, up) {
	logdens <- 0.0
	res <- .C("dologdensityAR1", as.double(x), as.double(mean), as.double(sd), as.double(sigma.eps), as.double(rho.c), 
				as.double(data), as.double(length(data)), as.double(low), as.double(up), logdens=logdens)
	return(res$logdens)	
}

logdensity.sigma.eps <- function(x, data, mu.c, rho.c, low, up, C, Ts) {
	log.d.cDens <- 0
	for(country in 1:C) {
		log.d.cDens <- log.d.cDens + sum(log(pmax(dnorm(data[[country]][2:Ts[country]], 
										mu.c[country]+rho.c[country]*(data[[country]][1:(Ts[country]-1)]-mu.c[country]), 
										sd=x),.Machine$double.xmin)))
	}
	return(log(dunif(x, low, up)) + log.d.cDens)
}

logdensity.sigma.mu.rho <- function(x, par, par.c, low, up) {
	log.d.cDens <- sum(log(pmax(dnorm(par.c, par, sd=x),.Machine$double.xmin)))
	return(log(dunif(x, low, up)) + log.d.cDens)
}

rnorm.trunc<-function(mean,sd,low,high=NA){
	if(is.na(high)) return(rnorm.ltrunc(mean, sd, low))
  temp<--999
  maxit <- 10
  i <- 1
  while((temp<low || temp>high) && i <= maxit) {
     temp<-rnorm(1,mean=mean,sd=sd)
     i <- i+1
  }
  if (i > maxit) {
  	temp <- if(temp<low) low else high
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

rnorm.ltrunc<-function(mean,sd,low){
  temp<--999
  maxit <- 10
  i <- 1
  while((temp<low) && i <= maxit) {
     temp<-rnorm(1,mean=mean,sd=sd)
     i <- i+1
  }
  if (i > maxit) {
  	temp <- low
  	#warning(paste('Maximum iterations reached in rnorm.ltrunc(', 
  	#			mean, ',', sd, '). Value truncated to ', temp, '.', sep=''), immediate.=TRUE)
  }
  return(temp)
}

slice.sampling <- function(x0, fun, width,  ..., low, up) {
	gx0 <- fun(x0, ..., low=low, up=up)
	z <- gx0 - rexp(1) # slice S={x: z < gx0}
	u <- runif(2)
	L <- max(x0 - width*u[1], low)
	R <- min(L + width, up)
	i<-1
	maxit <- 50
	while (TRUE) {
		u <- runif(1)
		x1 <- L + u*(R-L)
		if(z < fun(x1,  ..., low=low, up=up)) return(x1)
		if (x1 < x0) L <- x1
		else R <- x1
		i <- i+1
		if(i>maxit) stop('Problem in slice sampling')
	}
}

rgamma.ltrunc<-function(shape,rate,low){
  temp<--999
  maxit <- 10
  i <- 1
  while(temp<low && i <= maxit) {
     temp<-rgamma(1,shape=shape,rate=rate)
     i <- i+1
  }
  if (i > maxit) {
  	temp <- low
  	#warning(paste('Maximum iterations reached in rgamma.ltrunc(', shape, ',', rate, ').', sep=''), immediate.=TRUE)
  }
  return(temp)
}

