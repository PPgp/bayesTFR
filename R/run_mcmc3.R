run.tfr3.mcmc <- function(sim.dir, nr.chains=3, iter=50000,
						thin=10, replace.output=FALSE,
						# meta parameters
						my.tfr.file = NULL, buffer.size=100,
						use.extra.countries=FALSE,
						mu.prior.range=c(0, 2.1), 
						rho.prior.range=c(0,1-.Machine$double.xmin),
						 sigma.mu.prior.range=c(1e-5,0.318),
						 sigma.rho.prior.range=c(1e-5,0.289),
						 sigma.eps.prior.range=c(1e-5, 0.5),
						# starting values (length of 1 or nr.chains)
						mu.ini = NULL, mu.ini.range=mu.prior.range, 
						rho.ini=NULL, rho.ini.range=rho.prior.range, 
						sigma.mu.ini=NULL, sigma.mu.ini.range=sigma.mu.prior.range,
						sigma.rho.ini=NULL, sigma.rho.ini.range=sigma.rho.prior.range,
						sigma.eps.ini=NULL, sigma.eps.ini.range=sigma.eps.prior.range,
					 	seed = NULL, parallel=FALSE, nr.nodes=nr.chains, 
					 	compression.type='None',
					 	auto.conf = list(max.loops=5, iter=50000, iter.incr=20000, nr.chains=3, thin=60, burnin=10000),
						verbose=FALSE, verbose.iter = 1000, ...) {
	get.init.values <- function(range) {
		ifelse(rep(nr.chains==1, nr.chains), sum(range)/2, 
				#seq(range[1], to=range[2], length=nr.chains)
				runif(nr.chains, range[1], range[2])
				)
	}
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
	for(varname in c('mu.ini', 'rho.ini', 'sigma.mu.ini', 'sigma.rho.ini', 'sigma.eps.ini')) {
		if(is.null(get(varname)))
			assign(varname, get.init.values(get(paste(varname, '.range', sep=''))))
	}
	c.index <- 1: (if(use.extra.countries) get.nr.countries(mc$meta) else get.nr.countries.est(mc$meta))
	bayesTFR.mcmc.meta <- structure(list(nr.chains=nr.chains,
								my.tfr.file=my.tfr.file, output.dir=output.dir,
								phase=3, id_phase3 = which(mc$meta$lambda_c[c.index] < mc$meta$T_end_c[c.index]),
								nr.countries=sum(mc$meta$lambda_c[c.index] < mc$meta$T_end_c[c.index]),
								mu.prior.range=mu.prior.range, rho.prior.range=rho.prior.range,
						 		sigma.mu.prior.range=sigma.mu.prior.range, 
						 		sigma.rho.prior.range=sigma.rho.prior.range,
						 		sigma.eps.prior.range=sigma.eps.prior.range,
								mu.ini = mu.ini, mu.ini.range=mu.ini.range, 
								rho.ini=rho.ini, rho.ini.range=rho.ini.range, 
								sigma.mu.ini=sigma.mu.ini, sigma.mu.ini.range=sigma.mu.ini.range,
								sigma.rho.ini=sigma.rho.ini, sigma.rho.ini.range=sigma.rho.ini.range,
								sigma.eps.ini=sigma.eps.ini, sigma.eps.ini.range=sigma.eps.ini.range,
								compression.type=compression.type, buffer.size=buffer.size, auto.conf=auto.conf
								), class='bayesTFR.mcmc.meta')	
	store.bayesTFR.meta.object(bayesTFR.mcmc.meta, output.dir)
	meta <- bayesTFR.mcmc.meta
	if(meta$nr.countries <= 0) return(NULL)
	meta$parent <- mc$meta
	meta$regions <- mc$meta$regions
	# propagate initial values for all chains if needed
    starting.values <- list()
    for (var in c('mu.ini', 'rho.ini', 'sigma.mu.ini', 'sigma.rho.ini', 'sigma.eps.ini', 'iter')) {
    	if (length(get(var)) < nr.chains) 
            assign(var, rep(get(var), nr.chains)[1:nr.chains])
        if (var != 'iter') starting.values[[var]] <- get(var)
    }

	if (parallel) { # run chains in parallel
		chain.set <- bDem.performParallel(nr.nodes, 1:nr.chains, mcmc3.run.chain, 
                                     initfun=init.nodes, meta=meta, 
                                     thin=thin, iter=iter, 
                                     starting.values=starting.values,                                     
                                     verbose=verbose, verbose.iter=verbose.iter, ...)
	} else { # run chains sequentially
		chain.set <- list()
		for (chain in 1:nr.chains) {
			chain.set[[chain]] <- mcmc3.run.chain(chain, meta, thin=thin, 
                                                iter=iter, starting.values=starting.values, 
                                                verbose=verbose, verbose.iter=verbose.iter)
		}
	}
	names(chain.set) <- 1:nr.chains
	mcmc.set <- structure(list(meta=meta, mcmc.list=chain.set), class='bayesTFR.mcmc.set')
    cat('\nResults stored in', sim.dir,'\n')
    
    if(auto.run) {
		diag <- try(tfr3.diagnose(sim.dir=sim.dir, express=TRUE, 
						thin=auto.conf$thin, burnin=auto.conf$burnin,
						verbose=verbose))
		if(auto.conf$max.loops>1) {
			for(loop in 2:auto.conf$max.loops) {
				if(!inherits(diag, "try-error") && has.mcmc.converged(diag)) break
				mcmc.set <- continue.tfr3.mcmc(sim.dir=sim.dir, iter=auto.conf$iter.incr, nr.nodes=nr.nodes,
										  parallel=parallel, verbose=verbose, verbose.iter=verbose.iter, ...)
				diag <- try(tfr3.diagnose(sim.dir=sim.dir, express=TRUE, 
							thin=auto.conf$thin, burnin=auto.conf$burnin,
							verbose=verbose))
			}
		}
	}
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
	mcmc <- tfr3.mcmc.sampling(mcmc, thin=thin, verbose=verbose, verbose.iter=verbose.iter)
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

continue.tfr3.mcmc <- function(sim.dir, iter, chain.ids=NULL, parallel=FALSE, nr.nodes=NULL, auto.conf = NULL, 
								verbose=FALSE, verbose.iter=1000, ...) {
	mcmc.set <- get.tfr3.mcmc(sim.dir)

	auto.run <- FALSE
	if(iter == 'auto') { # defaults for auto-run (includes convergence diagnostics)
		default.auto.conf <- mcmc.set$meta$auto.conf
		if(is.null(auto.conf)) auto.conf <- list()
		for (par in names(default.auto.conf))
			if(is.null(auto.conf[[par]])) auto.conf[[par]] <- default.auto.conf[[par]]
		iter <- auto.conf$iter.incr
		auto.run <- TRUE
		fiter <- sapply(mcmc.set$mcmc.list, function(x) x$finished.iter)
		if (!all(fiter== fiter[1])) stop('All chains must be of the same length if the "auto" option is used.')
	}
	if (is.null(chain.ids) || auto.run) {
		chain.ids <- names(mcmc.set$mcmc.list)
	}
	if (parallel) { # run chains in parallel
		if(is.null(nr.nodes)) nr.nodes<-length(chain.ids)
		chain.list <- bDem.performParallel(nr.nodes, chain.ids, mcmc3.continue.chain, 
						initfun=init.nodes, mcmc.list=mcmc.set$mcmc.list, iter=iter, verbose=verbose, 
						verbose.iter=verbose.iter, ...)
		for (i in 1:length(chain.ids))
			mcmc.set$mcmc.list[[chain.ids[i]]] <- chain.list[[i]]
	} else { # run chains sequentially
		for (chain.id in chain.ids) {
			mcmc.set$mcmc.list[[chain.id]] <- mcmc3.continue.chain(chain.id, mcmc.set$mcmc.list, 
												iter=iter, verbose=verbose, verbose.iter=verbose.iter)
		}
	}
	cat('\n')
	if(auto.run) {
		diag <- try(tfr3.diagnose(sim.dir=sim.dir, express=TRUE, 
						thin=auto.conf$thin, burnin=auto.conf$burnin,
						verbose=verbose))
		if(auto.conf$max.loops>1) {
			for(loop in 2:auto.conf$max.loops) {
				if(!inherits(diag, "try-error") && has.mcmc.converged(diag)) break
				mcmc.set <- continue.tfr3.mcmc(sim.dir=sim.dir, iter=auto.conf$iter.incr, nr.nodes=nr.nodes,
										  parallel=parallel, verbose=verbose, verbose.iter=verbose.iter)
				diag <- try(tfr3.diagnose(sim.dir=sim.dir, express=TRUE, 
							thin=auto.conf$thin, burnin=auto.conf$burnin,
							verbose=verbose))
			}
		}
	}
	invisible(mcmc.set)
}
	
mcmc3.continue.chain <- function(chain.id, mcmc.list, iter, verbose=FALSE, verbose.iter=10) {
	cat('\n\nChain nr.', chain.id, '\n')
	if (verbose)
		cat('************\n')
	mcmc <- mcmc.list[[chain.id]]
	mcmc$iter <- mcmc$finished.iter + iter
	if (verbose) 
		cat('Continue sampling -', iter, 'additional iterations,', mcmc$iter, 'iterations in total.\n')

	mcmc <- tfr3.mcmc.sampling(mcmc, thin=mcmc$thin, start.iter=mcmc$finished.iter+1, verbose=verbose, verbose.iter=verbose.iter)
	return(mcmc)
}

run.tfr3.mcmc.subn1d <- function(countries, sim.dir=file.path(getwd(), 'bayesTFR.output'),
										iter=10000, world.sim.dir = sim.dir, post.burnin = 10000, 
										buffer.size=100, seed = NULL, parallel=FALSE, nr.nodes=NULL, 
										verbose=FALSE, verbose.iter=100, ...) {
	if(is.null(world.sim.dir)) world.sim.dir <- sim.dir
	if(!has.tfr3.mcmc(world.sim.dir)) stop('No phase III available in ', world.sim.dir)
	output.dir <- file.path(sim.dir, 'subnat')
	mcmc2 <- get.tfr.mcmc(world.sim.dir)	
	mcmc3 <- get.tfr3.mcmc(world.sim.dir)
	world.use.tfr3 <- FALSE
	post.thin <- NULL
	if(has.tfr.prediction(sim.dir=world.sim.dir)) {
		p <- get.tfr.prediction(sim.dir=world.sim.dir)
		world.use.tfr3 <- p$use.tfr3
	}
	if(world.use.tfr3) {
		post.burnin <- p$burnin3
		total.iter <- p$nr.traj
		post.thin <- p$thin3
	} else {
		l <- mcmc3$mcmc.list[[1]]$finished.iter
		if(post.burnin > l) {
			post.burnin <- as.integer(l/2)
			warning('post.burnin larger than MCMC length. Adjusted to ', post.burnin)
		}
		total.iter <- mcmc3$mcmc.list[[1]]$length - get.thinned.burnin(mcmc3$mcmc.list[[1]], post.burnin)
	}
	if(is.null(iter)) iter <- total.iter
	if(iter > total.iter) post.idx <- sample(1:total.iter, iter, replace=TRUE)
	else {
		post.idx <- get.thinning.index(iter, total.iter)
		iter <- post.idx$nr.points
		post.idx <- post.idx$index
	}
	sampling.idx <- 1:iter			
	mcmc3$meta$buffer.size <- buffer.size
	if(!is.null(seed)) set.seed(seed)
	result <- list()
	for (country in countries) {
		country.obj <- get.country.object(country, mcmc2$meta)	
		if(verbose) 
			cat('\nStoring Phase III for sub-regions of ', country.obj$name, '\n')
		m2.sim.dir <- file.path(output.dir, paste0('c', country.obj$code))
		if(!has.tfr.mcmc(m2.sim.dir)) {
			warnings('No Phase II MCMC for ', country.obj$name)
			next
		}
		m2 <- get.tfr.mcmc(m2.sim.dir)
		this.output.dir <- file.path(m2.sim.dir, 'phaseIII')		
		m3 <- mcmc3
		m3$meta$output.dir <- this.output.dir
		m3$meta$id_phase3 <- which(m2$meta$lambda_c < m2$meta$T_end_c)
		m3$meta$parent <- m2$meta
		m3$meta$regions <- m2$meta$regions
		m3$meta$nr.countries <- length(m3$meta$id_phase3)
		bayesTFR.mcmc.meta <- m3$meta
		if(file.exists(this.output.dir)) unlink(this.output.dir, recursive=TRUE)
		dir.create(this.output.dir)
		store.bayesTFR.meta.object(bayesTFR.mcmc.meta, this.output.dir)
		m3$mcmc.list <- hyperpars3.from.country.pars(mcmc3$mcmc.list, mcmc3$meta, m3$meta, world.use.tfr3, this.output.dir, country.obj, post.idx,  
								burnin=post.burnin, thin=post.thin)
		nr.chains <- length(m3$mcmc.list)
		if(length(m3$meta$id_phase3) > 0) {
			nrc <- get.nr.countries(m3$meta)
			if (parallel  && nr.chains > 1) { # run chains in parallel
				if(is.null(nr.nodes)) nr.nodes <- nr.chains
				chain.set <- bDem.performParallel(nr.nodes, 1:nr.chains, mcmc3.run.chain.subnat.extra, 
						initfun=init.nodes, mcmc.set=m3, region.idx=1:nrc, 
						posterior.sample=sampling.idx, iter=iter, burnin=0, verbose=verbose, 
                        verbose.iter=verbose.iter, ...)
			} else { # run chains sequentially
				chain.set <- list()
				for (chain in 1:nr.chains) {
					chain.set[[chain]] <- mcmc3.run.chain.subnat.extra(chain, m3, region.idx=1:nrc, posterior.sample=sampling.idx,
												iter=iter, burnin=0, verbose=verbose, verbose.iter=verbose.iter)
				}
			}
			names(chain.set) <- 1:nr.chains
			m3$mcmc.list <- chain.set
		}
		result[[as.character(country.obj$code)]] <- structure(list(meta=bayesTFR.mcmc.meta, mcmc.list=m3$mcmc.list), 
																			class='bayesTFR.mcmc.set')
	}
	invisible(result)				
}

mcmc3.run.chain.subnat.extra <- function(chain.id, mcmc.set, region.idx, posterior.sample, 
												iter=NULL, burnin=0, verbose=FALSE, verbose.iter=100) {
	cat('\n\nChain nr.', chain.id, '\n')
	if (verbose)
		cat('************\n')
	mcmc <- mcmc.set$mcmc.list[[chain.id]]
	mcmc$length <- mcmc$finished.iter <- iter
	mcmc$rho.c <- rep(mcmc$rho.ini, get.nr.countries(mcmc.set$meta))
	mcmc$mu.c <- rep(mcmc$mu.ini, get.nr.countries(mcmc.set$meta))
	mcmc <- tfr3.mcmc.sampling.extra(mcmc, mcmc.list=mcmc.set$mcmc.list, countries=region.idx, 
									posterior.sample=posterior.sample, do.sample.posterior=FALSE,
									iter=iter, burnin=burnin, verbose=verbose, verbose.iter=verbose.iter)
	return(mcmc)
}

hyperpars3.from.country.pars <- function(mcmc.list, world.meta, meta, thin.and.collapse, output.dir, country.obj, posterior.idx, burnin, thin=NULL) {
	if(thin.and.collapse) { # thin and collapse so that it matches predictions
		new.mcmc.list <- mcmc.list[1]
	} else {
		new.mcmc.list <- mcmc.list
	}
	#if(country.obj$code == 203) stop('')
	for (i in 1:length(new.mcmc.list)) {
		mc <- if(thin.and.collapse) mcmc.list else mcmc.list[[i]]
		if(country.obj$index %in% world.meta$id_phase3) { # country has phase 3 pars
			traces <- data.frame(get.tfr3.parameter.traces.cs(mc, country.obj,
									par.names=c('mu.c', 'rho.c'), burnin=burnin, thin=thin))[posterior.idx,]
			colnames(traces) <- c('mu', 'rho')
		} else {
			traces <- data.frame(get.tfr3.parameter.traces(mc, par.names=c('mu', 'rho'), 
									burnin=burnin, thin=thin))[posterior.idx,]
		}
		traces <- cbind(traces, data.frame(get.tfr3.parameter.traces(mc, 
									par.names=c('sigma.eps', 'sigma.rho', 'sigma.mu'), burnin=burnin, thin=thin))[posterior.idx,]) 
		result <- list()
		for(par in colnames(traces)) result[[par]] <- traces[,par]
		m <- new.mcmc.list[[i]]
		m$length <- nrow(traces)
		m$finished.iter <- m$length
		m$thin <- 1
		m$meta <- meta
		write.list.into.file.cindep(m, result)
		store.bayesTFR.object(m, file.path(output.dir, m$output.dir))
		new.mcmc.list[[i]] <- m
	}
	names(new.mcmc.list) <- 1:length(new.mcmc.list)
	return(new.mcmc.list)
}

