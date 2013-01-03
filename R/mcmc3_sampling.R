#########################################################
# MCMC sampling for Phase III
#########################################################

tfr.mcmc3.sampling <- function(mcmc, thin=1, start.iter=2, verbose=FALSE, verbose.iter=10) {

	if (!is.null(mcmc$rng.state)) .Random.seed <- mcmc$rng.state
	meta <- mcmc$meta
    niter <- mcmc$iter
    nr.countries <- meta$nr.countries
    countries.index <- meta$id_phase3
    
    # Start MCMC
	############
    for (iter in start.iter:niter) {
    	if(verbose.iter > 0 && (iter %% verbose.iter == 0))
        	cat('\nIteration:', iter, '--', date())
        unblock.gtk('bDem.TFRmcmc')

		for(par in c('mu', 'rho')) {
			par.c <- paste(par, 'c', sep='.')			
			par.range <- paste(par, 'prior', 'range', sep='.')
			sum.par <- sum(mcmc[[par.c]])
			sig.sq <- mcmc[[par.sig]]^2
			par.mean <- sum.par/nr.countries
			par.sd <- sqrt(sig.sq/nr.countries)
			#stop('')
			mcmc[[par]] <- rnorm.trunc(mean=par.mean, sd=par.sd, low=meta[[par.range]][1], high=meta[[par.range]][2])
			par.sig <- paste('sigma', par, sep='.')
			par.sig.range <- paste(par.sig, 'prior', 'range', sep='.')
			S <- sum((mcmc[[par.c]]-mcmc[[par]])^2)
			mcmc[[par.sig]] <- 1/sqrt(rgamma.ltrunc(nr.countries/2 - 3, S/2, low=1/sqrt(meta[[par.sig.range]][2])))
		}
		for(country in 1:nr.countries) {
			data <- get.observed.tfr(countries.index[country], meta$parent, 'tfr_matrix_all')
			ardata <- data[meta$parent$lambda_c[countries.index[country]]:meta$parent$T_end_c[countries.index[country]]]
			mcmc$mu.c[country] <- slice.sampling(mcmc$mu.c[country],
							logdensity.mu.rho.c, 1, mean=mcmc$mu, sd=mcmc$sigma.mu, sigma.eps=mcmc$sigma.eps, 
							rho.c=mcmc$rho.c[country], data=ardata, low=0, up=999999
							)
			mcmc$rho.c[country] <- slice.sampling(mcmc$rho.c[country],
							logdensity.mu.rho.c, 1, mean=mcmc$rho, sd=mcmc$sigma.rho, sigma.eps=mcmc$sigma.eps, 
							rho.c=mcmc$rho.c[country], data=ardata, low=0, up=1
							)
		}
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

logdensity.mu.rho.c <- function(x, mean, sd, sigma.eps, rho.c, data, low, up) {
	logdens <- 0.0
	res <- .C("dologdensityAR1", as.double(x), as.double(mean), as.double(sd), as.double(sigma.eps), as.double(rho.c), 
				as.double(data), as.double(length(data)), as.double(low), as.double(up), logdens=logdens)
	return(res$logdens)	
}

rnorm.trunc<-function(mean,sd,low,high){
  temp<--999
  maxit <- 10
  i <- 1
  while((temp<low || temp>high) && i <= maxit) {
     temp<-rnorm(1,mean=mean,sd=sd)
     i <- i+1
  }
  if (i > maxit) {
  	temp <- if(temp<low) low else high
  	warning(paste('Maximum iterations reached in rnorm.trunc(', 
  				mean, ',', sd, '). Value truncated to ', temp, '.', sep=''), immediate.=TRUE)
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
  	warning(paste('Maximum iterations reached in rgamma.ltrunc(', shape, ',', rate, ').', sep=''), immediate.=TRUE)
  }
  return(temp)
}