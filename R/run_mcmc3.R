run.tfr3.mcmc <- function(sim.dir, nr.chains=3, iter=5000,
						thin=1, replace.output=FALSE,
						# meta parameters
						my.tfr.file = NULL, buffer.size=100,
						mu.prior.range=c(0, 2.1), 
						 rho.prior.range=c(0,1),
						 sigma.mu.prior.range=c(0,0.318),
						 sigma.rho.prior.range=c(0,0.289),
						# starting values (length of 1 or nr.chains)
						mu.ini = NULL, mu.ini.range=mu.prior.range, 
						rho.ini=NULL, rho.ini.range=rho.prior.range, 
						sigma.mu.ini=NULL, sigma.mu.ini.range=sigma.mu.prior.range,
						sigma.rho.ini=NULL, sigma.rho.ini.range=sigma.rho.prior.range,
						sigma.eps.ini=0.01,
					 	seed = NULL, parallel=FALSE, nr.nodes=nr.chains, 
					 	compression.type='None',
					 	auto.conf = list(max.loops=5, iter=5000, iter.incr=2000, nr.chains=3, thin=5, burnin=500),
						verbose=FALSE, verbose.iter = 100, ...) {
	get.init.values <- function(range)
		ifelse(rep(nr.chains==1, nr.chains), sum(range)/2, seq(range[1], to=range[2], length=nr.chains))

	mc <- get.tfr.mcmc(sim.dir)
	output.dir <- file.path(sim.dir, 'phaseIII')
	if(file.exists(output.dir)) {
		if(!replace.output) stop('MCMCs for Phase III already exist in ', sim.dir, 
                        '.\nSet replace.output=TRUE if you want to overwrite existing results.')
        unlink(output.dir, recursive=TRUE)
	}
	dir.create(output.dir)
	default.auto.conf <- formals(run.tfr3.mcmc)$auto.conf
	for (par in names(default.auto.conf))
		if(is.null(auto.conf[[par]])) auto.conf[[par]] <- default.auto.conf[[par]]
	auto.run <- FALSE
	if(iter == 'auto') { # defaults for auto-run (includes convergence diagnostics)
		iter <- auto.conf$iter
		nr.chains <- auto.conf$nr.chains
		auto.run <- TRUE		
	}
	if (verbose) {
		cat('\nStarting Bayesian Hierarchical Model for TFR - Phase III.\n')
		cat('=========================================================\n')
		cat('Initialize simulation -', nr.chains, 'chain(s) in total.\n')
	}
	if(!is.null(seed)) set.seed(seed)
	
	# starting values
	#==============================
	for(varname in c('mu.ini', 'rho.ini', 'sigma.mu.ini', 'sigma.rho.ini')) {
		if(is.null(get(varname)))
			assign(varname, get.init.values(get(paste(varname, '.range', sep=''))))
	}
	bayesTFR.mcmc.meta <- structure(list(nr.chains=nr.chains,
								my.tfr.file=my.tfr.file, output.dir=output.dir,
								phase=3, id_phase3 = which(mc$meta$lambda_c < mc$meta$T_end_c),
								nr.countries=sum(mc$meta$lambda_c < mc$meta$T_end_c),
								mu.prior.range=mu.prior.range, rho.prior.range=rho.prior.range,
						 		sigma.mu.prior.range=sigma.mu.prior.range, 
						 		sigma.rho.prior.range=sigma.rho.prior.range,
								mu.ini = mu.ini, mu.ini.range=mu.ini.range, 
								rho.ini=rho.ini, rho.ini.range=rho.ini.range, 
								sigma.mu.ini=sigma.mu.ini, sigma.mu.ini.range=sigma.mu.ini.range,
								sigma.rho.ini=sigma.rho.ini, sigma.rho.ini.range=sigma.rho.ini.range,
								sigma.eps.ini=sigma.eps.ini,
								compression.type=compression.type, buffer.size=buffer.size, auto.conf=auto.conf
								))	
	store.bayesTFR.meta.object(bayesTFR.mcmc.meta, output.dir)
	bayesTFR.mcmc.meta$parent <- mc$meta
	# propagate initial values for all chains if needed
    starting.values <- list()
    for (var in c('mu.ini', 'rho.ini', 'sigma.mu.ini', 'sigma.rho.ini', 'iter')) {
    	if (length(get(var)) < nr.chains) 
            assign(var, rep(get(var), nr.chains)[1:nr.chains])
        if (var != 'iter') starting.values[[var]] <- get(var)
    }

	if (parallel) { # run chains in parallel
		require(snowFT)
		chain.set <- performParallel(nr.nodes, 1:nr.chains, mcmc3.run.chain, 
                                     initfun=init.nodes, meta=bayesTFR.mcmc.meta, 
                                     thin=thin, iter=iter, 
                                     starting.values=starting.values,                                     
                                     verbose=verbose, verbose.iter=verbose.iter, ...)
	} else { # run chains sequentially
		chain.set <- list()
		for (chain in 1:nr.chains) {
			chain.set[[chain]] <- mcmc3.run.chain(chain, bayesTFR.mcmc.meta, thin=thin, 
                                                iter=iter, starting.values=starting.values, 
                                                verbose=verbose, verbose.iter=verbose.iter)
		}
	}
	names(chain.set) <- 1:nr.chains
	mcmc.set <- structure(list(meta=bayesTFR.mcmc.meta, mcmc.list=chain.set), class='bayesTFR.mcmc.set')
    cat('\nResults stored in', sim.dir,'\n')
    
    # if(auto.run) {
		# diag <- try(tfr.diagnose(sim.dir=output.dir, keep.thin.mcmc=TRUE, 
						# thin=auto.conf$thin, burnin=auto.conf$burnin,
						# verbose=verbose))
		# if(auto.conf$max.loops>1) {
			# for(loop in 2:auto.conf$max.loops) {
				# if(!inherits(diag, "try-error") && has.mcmc.converged(diag)) break
				# mcmc.set <- continue.e0.mcmc(iter=auto.conf$iter.incr, output.dir=output.dir, nr.nodes=nr.nodes,
										  # parallel=parallel, verbose=verbose, verbose.iter=verbose.iter)
				# diag <- try(e0.diagnose(sim.dir=output.dir, keep.thin.mcmc=TRUE, 
							# thin=auto.conf$thin, burnin=auto.conf$burnin,
							# verbose=verbose))
			# }
		# }
	# }
    if (verbose) 
		cat('\nSimulation successfully finished!!!\n')
    invisible(mcmc.set)
				
}

mcmc3.run.chain <- function(chain.id, meta, thin=1, iter=100, starting.values=NULL,
							verbose=FALSE, verbose.iter=10) {
								
	cat('\n\nChain nr.', chain.id, '\n')
    if (verbose) 
    	cat('************\n')
    	
    this.sv <- list()
	for(var in names(starting.values)) {
		this.sv[[var]] <- starting.values[[var]][chain.id]
	}
	mcmc <- do.call('mcmc3.ini', c(list(chain.id, meta, iter=iter[chain.id], thin=thin, starting.values=this.sv) ))
	if (verbose) {
    	cat('Starting values:\n')
		print(unlist(mcmc[names(starting.values)]))
        cat('Store initial values into ', mcmc$output.dir, '\n')
    }
	store.mcmc3(mcmc, append=FALSE, flush.buffer=TRUE, verbose=verbose)
	
	if (verbose) 
		cat('Start sampling -', mcmc$iter, 'iterations in total.\n')
	mcmc <- tfr.mcmc3.sampling(mcmc, thin=thin, verbose=verbose, verbose.iter=verbose.iter)
	return(mcmc)
}

mcmc3.ini <- function(chain.id, mcmc.meta, iter=100, thin=1, starting.values=NULL,
				     verbose=FALSE) {
                                                        
    if (!exists(".Random.seed")) runif(1)
	mcmc <- structure(c(starting.values, list(
        				output.dir=paste('mc', chain.id, sep=''), 
        				thin=thin, finished.iter=1, length = 1,
        				iter=iter, id=chain.id, traces=0,
        				traces.burnin=0, rng.state = .Random.seed,
        				compression.type=mcmc.meta$compression.type,
        				meta = mcmc.meta)), class='bayesTFR.mcmc')
    # country-specific initial values
    for(varname in c('mu', 'rho')) {
    	var <- paste(varname, 'c', sep='.')
    	range.var <- paste(varname,'ini.range', sep='.')
    	mcmc[[var]] <- runif(mcmc.meta$nr.countries, mcmc.meta[[range.var]][1], mcmc.meta[[range.var]][2])
    	mcmc[[varname]] <- starting.values[[paste(varname,'ini', sep='.')]]
    }
    for(varname in c('sigma.mu', 'sigma.rho', 'sigma.eps')) {
    	mcmc[[varname]] <- starting.values[[paste(varname,'ini', sep='.')]]
    }
    return(mcmc) 
}

