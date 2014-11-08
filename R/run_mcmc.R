
run.tfr.mcmc <- function(nr.chains=3, iter=62000, output.dir=file.path(getwd(), 'bayesTFR.output'), 
						thin=1, replace.output=FALSE,
						# meta parameters
						start.year=1750, present.year=2010, wpp.year=2012,
						my.tfr.file = NULL, buffer.size=100,
					 	U.c.low=5.5, U.up=8.8, U.width=3,
					 	mean.eps.tau0 = -0.25, sd.eps.tau0 = 0.4, nu.tau0 = 2,                                                
        				Triangle_c4.low = 1, Triangle_c4.up = 2.5,
        				Triangle_c4.trans.width=2, Triangle4.0 = 0.3, 
        				delta4.0 = 0.8, nu4 = 2,
					 	S.low=3.5, S.up=6.5, S.width=0.5,
					 	a.low=0, a.up=0.2, a.width=0.02,
					 	b.low=a.low, b.up=a.up, b.width=0.02,
					 	sigma0.low=0.01, sigma0.up=0.6, sigma0.width=0.1,
					 	sigma0.min=0.001, 
					 	const.low=0.8, const.up=2, const.width=0.3,
					 	d.low=0.05, d.up=0.5, d.trans.width=1,
					 	chi0=-1.5, psi0=0.6, nu.psi0=2,
					 	alpha0.p=c(-1, 0.5, 1.5), delta0=1, nu.delta0=2,
					 	dl.p1=9, dl.p2=9,
						# starting values (length of 1 or nr.chains)
						S.ini=NULL, a.ini=NULL, b.ini=NULL, 
					 	sigma0.ini=NULL, Triangle_c4.ini=NULL, const.ini=NULL, gamma.ini=1, 
					 	proposal_cov_gammas = NULL, # should be a list with elements 'values' and 'country_codes'
					 	seed = NULL, parallel=FALSE, nr.nodes=nr.chains, 
					 	save.all.parameters = FALSE, compression.type='None',
					 	auto.conf = list(max.loops=5, iter=62000, iter.incr=10000, nr.chains=3, thin=80, burnin=2000),
						verbose=FALSE, verbose.iter = 10, ...) {

	if(file.exists(output.dir)) {
		if(length(list.files(output.dir)) > 0 & !replace.output)
			stop('Non-empty directory ', output.dir, 
			' already exists.\nSet replace.output=TRUE if you want to overwrite existing results.')
		unlink(output.dir, recursive=TRUE)
	}
	dir.create(output.dir)
	
	default.auto.conf <- formals(run.tfr.mcmc)$auto.conf
	for (par in names(default.auto.conf))
		if(is.null(auto.conf[[par]])) auto.conf[[par]] <- default.auto.conf[[par]]
	auto.run <- FALSE
	if(iter == 'auto') { # defaults for auto-run (includes convergence diagnostics)
		iter <- auto.conf$iter
		nr.chains <- auto.conf$nr.chains
		auto.run <- TRUE		
	}
	
	if (verbose) {
		cat('\nStarting Bayesian Hierarchical Model for TFR - Phase II.\n')
		cat('========================================================\n')
		cat('Initialize simulation -', nr.chains, 'chain(s) in total.\n')
	}
	if(!is.null(seed)) set.seed(seed)
	
	# starting values (length of 1 or nr.chains)
	if (missing(S.ini) || is.null(S.ini)) 
		S.ini <- ifelse(rep(nr.chains==1, nr.chains), 
					 		(S.low+S.up)/2, seq(S.low, to=S.up, length=nr.chains))
	if (missing(a.ini) || is.null(a.ini)) 
		a.ini <- ifelse(rep(nr.chains==1, nr.chains), 
					 		(a.low+a.up)/2, seq(a.low, to=a.up, length=nr.chains))
	if (missing(b.ini) || is.null(b.ini)) 
		b.ini <- ifelse(rep(nr.chains==1, nr.chains), 
					 		(b.low+b.up)/2, seq(b.low, to=b.up, length=nr.chains))
	if (missing(sigma0.ini) || is.null(sigma0.ini)) 
		sigma0.ini <- ifelse(rep(nr.chains==1, nr.chains), 
					 		(sigma0.low+sigma0.up)/2, 
					 		seq(sigma0.low, to=sigma0.up, length=nr.chains))
	if (missing(Triangle_c4.ini) || is.null(Triangle_c4.ini)) 
		Triangle_c4.ini <- ifelse(rep(nr.chains==1, nr.chains), 
					 		(Triangle_c4.low+Triangle_c4.up)/2, 
					 		seq(Triangle_c4.low+0.0001, to=Triangle_c4.up-0.0001, length=nr.chains))
	if (missing(const.ini) || is.null(const.ini))  
		const.ini <- ifelse(rep(nr.chains==1, nr.chains), 
					 		(const.low+const.up)/2, 
					 		seq(const.low, to=const.up, length=nr.chains))
					 		
	bayesTFR.mcmc.meta <- mcmc.meta.ini(
						nr.chains=nr.chains,
						start.year=start.year, present.year=present.year, 
						wpp.year=wpp.year, my.tfr.file = my.tfr.file,
						output.dir=output.dir, phase=2,
					 	U.c.low=U.c.low, U.up=U.up, U.width=U.width,
					 	mean.eps.tau0=mean.eps.tau0, sd.eps.tau0 = sd.eps.tau0, nu.tau0 = nu.tau0,                                            
        				Triangle4.0 = Triangle4.0,  
        				Triangle_c4.low = Triangle_c4.low , Triangle_c4.up = Triangle_c4.up,
        				Triangle_c4.trans.width=Triangle_c4.trans.width,
        				delta4.0 = delta4.0, nu4=nu4,
					 	S.low=S.low, S.up=S.up, S.width=S.width,
					 	a.low=a.low, a.up=a.up, a.width=a.width,
					 	b.low=b.low, b.up=b.up, b.width=b.width,
					 	sigma0.low=sigma0.low, sigma0.up=sigma0.up, sigma0.width=sigma0.width,
					 	sigma0.min=sigma0.min, 
					 	const.low=const.low, const.up=const.up, const.width=const.width,
					 	d.low=d.low, d.up=d.up, d.trans.width=d.trans.width,
					 	chi0=chi0, psi0=psi0, nu.psi0=nu.psi0,
					 	alpha0.p = alpha0.p, delta0=delta0, nu.delta0=nu.delta0,
					 	dl.p1=dl.p1, dl.p2=dl.p2, 
					 	proposal_cov_gammas = proposal_cov_gammas,
					 	buffer.size=buffer.size, compression.type=compression.type, 
					 	auto.conf=auto.conf, verbose=verbose)
	store.bayesTFR.meta.object(bayesTFR.mcmc.meta, output.dir)
			
	# propagate initial values for all chains if needed
	for (var in c('S.ini', 'a.ini', 'b.ini', 'sigma0.ini', 'const.ini', 'gamma.ini', 'Triangle_c4.ini', 'iter')) {
		if (length(get(var)) < nr.chains) {
			if (length(get(var)) == 1) {
				assign(var, rep(get(var), nr.chains))
				} else {
				warning(var, ' has the wrong length. Either 1 or ', nr.chains, 
				' is allowed.\nValue set to ', get(var)[1], ' for all chains.')
				assign(var, rep(get(var)[1], nr.chains))
				}
			}
		}
	if (parallel) { # run chains in parallel
		chain.set <- bDem.performParallel(nr.nodes, 1:nr.chains, mcmc.run.chain, 
						initfun=init.nodes, meta=bayesTFR.mcmc.meta, 
						thin=thin, iter=iter, S.ini=S.ini, a.ini=a.ini,
                        b.ini=b.ini, sigma0.ini=sigma0.ini, Triangle_c4.ini=Triangle_c4.ini, const.ini=const.ini,
                        gamma.ini=gamma.ini, save.all.parameters=save.all.parameters, verbose=verbose, 
                        verbose.iter=verbose.iter, ...)
	} else { # run chains sequentially
		chain.set <- list()
		for (chain in 1:nr.chains) {
			chain.set[[chain]] <- mcmc.run.chain(chain, bayesTFR.mcmc.meta, thin=thin, 
					 	iter=iter, S.ini=S.ini, a.ini=a.ini, b.ini=b.ini, 
					 	sigma0.ini=sigma0.ini, Triangle_c4.ini=Triangle_c4.ini, const.ini=const.ini, 
					 	gamma.ini=gamma.ini, save.all.parameters=save.all.parameters,
					 	verbose=verbose, verbose.iter=verbose.iter)
		}
	}
	names(chain.set) <- 1:nr.chains
	mcmc.set <- structure(list(meta=bayesTFR.mcmc.meta, mcmc.list=chain.set), class='bayesTFR.mcmc.set')
	cat('\nResults stored in', output.dir,'\n')
	
	if(auto.run) {
		diag <- try(tfr.diagnose(sim.dir=output.dir, keep.thin.mcmc=TRUE, 
						thin=auto.conf$thin, burnin=auto.conf$burnin,
						verbose=verbose))
		if(auto.conf$max.loops>1) {
			for(loop in 2:auto.conf$max.loops) {
				if(!inherits(diag, "try-error") && has.mcmc.converged(diag)) break
				mcmc.set <- continue.tfr.mcmc(iter=auto.conf$iter.incr, output.dir=output.dir, nr.nodes=nr.nodes,
										  parallel=parallel, verbose=verbose, verbose.iter=verbose.iter, ...)
				diag <- try(tfr.diagnose(sim.dir=output.dir, keep.thin.mcmc=TRUE, 
							thin=auto.conf$thin, burnin=auto.conf$burnin,
							verbose=verbose))
			}
		}
	}
	if (verbose)
		cat('\nSimulation successfully finished!!!\n')
	invisible(mcmc.set)
}


mcmc.run.chain <- function(chain.id, meta, thin=1, iter=100, 
							S.ini, a.ini, b.ini, sigma0.ini, Triangle_c4.ini, const.ini, gamma.ini=1,
							save.all.parameters=FALSE,
							verbose=FALSE, verbose.iter=10) {
								
	cat('\n\nChain nr.', chain.id, '\n')
    if (verbose) {
    	cat('************\n')
    	cat('Starting values:\n')
    	sv <- c(S.ini[chain.id], a.ini[chain.id], b.ini[chain.id], sigma0.ini[chain.id], Triangle_c4.ini[chain.id],
    			const.ini[chain.id], gamma.ini[chain.id])
    	names(sv) <- c('S', 'a', 'b', 'sigma0', 'Triangle_c4', 'const', 'gamma')
    	print(sv)
    }

	mcmc <- mcmc.ini(chain.id, meta, iter=iter[chain.id],
                                     S.ini=S.ini[chain.id],
                                     a.ini=a.ini[chain.id],
                                     b.ini=b.ini[chain.id],
                                     sigma0.ini=sigma0.ini[chain.id],
                                     Triangle_c4.ini=Triangle_c4.ini[chain.id],
                                     const.ini=const.ini[chain.id],
                                     gamma.ini=gamma.ini[chain.id],
                                     save.all.parameters=save.all.parameters,
                                     verbose=verbose)
	
	if (verbose) 
		cat('Store initial values into ', mcmc$output.dir, '\n')
	store.mcmc(mcmc, append=FALSE, flush.buffer=TRUE, verbose=verbose)
	
	if (verbose) 
		cat('Start sampling -', mcmc$iter, 'iterations in total.\n')
	mcmc <- tfr.mcmc.sampling(mcmc, thin=thin, verbose=verbose, verbose.iter=verbose.iter)
	return(mcmc)
}
	
continue.tfr.mcmc <- function(iter, chain.ids=NULL, output.dir=file.path(getwd(), 'bayesTFR.output'), 
								parallel=FALSE, nr.nodes=NULL, auto.conf = NULL, verbose=FALSE, verbose.iter=10, ...) {
	mcmc.set <- get.tfr.mcmc(output.dir)

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
		chain.list <- bDem.performParallel(nr.nodes, chain.ids, mcmc.continue.chain, 
						initfun=init.nodes, mcmc.list=mcmc.set$mcmc.list, iter=iter, verbose=verbose, 
						verbose.iter=verbose.iter, ...)
		for (i in 1:length(chain.ids))
			mcmc.set$mcmc.list[[chain.ids[i]]] <- chain.list[[i]]
	} else { # run chains sequentially
		for (chain.id in chain.ids) {
			mcmc.set$mcmc.list[[chain.id]] <- mcmc.continue.chain(chain.id, mcmc.set$mcmc.list, 
												iter=iter, verbose=verbose, verbose.iter=verbose.iter)
		}
	}
	cat('\n')
	if(auto.run) {
		diag <- try(tfr.diagnose(sim.dir=output.dir, keep.thin.mcmc=TRUE, 
						thin=auto.conf$thin, burnin=auto.conf$burnin,
						verbose=verbose))
		if(auto.conf$max.loops>1) {
			for(loop in 2:auto.conf$max.loops) {
				if(!inherits(diag, "try-error") && has.mcmc.converged(diag)) break
				mcmc.set <- continue.tfr.mcmc(iter=auto.conf$iter.incr, output.dir=output.dir, nr.nodes=nr.nodes,
										  parallel=parallel, verbose=verbose, verbose.iter=verbose.iter, ...)
				diag <- try(tfr.diagnose(sim.dir=output.dir, keep.thin.mcmc=TRUE, 
							thin=auto.conf$thin, burnin=auto.conf$burnin,
							verbose=verbose))
			}
		}
	}
	invisible(mcmc.set)
}
	
mcmc.continue.chain <- function(chain.id, mcmc.list, iter, verbose=FALSE, verbose.iter=10) {
	cat('\n\nChain nr.', chain.id, '\n')
	if (verbose)
		cat('************\n')
	mcmc <- mcmc.list[[chain.id]]
	mcmc$iter <- mcmc$finished.iter + iter
	if (verbose) 
		cat('Continue sampling -', iter, 'additional iterations,', mcmc$iter, 'iterations in total.\n')

	mcmc <- tfr.mcmc.sampling(mcmc, thin=mcmc$thin, start.iter=mcmc$finished.iter+1, verbose=verbose, verbose.iter=verbose.iter)
	return(mcmc)
}

run.tfr.mcmc.extra <- function(sim.dir=file.path(getwd(), 'bayesTFR.output'), 
								countries = NULL, my.tfr.file = NULL, iter = NULL,
								thin=1, burnin=2000, parallel=FALSE, nr.nodes=NULL, 
								verbose=FALSE, verbose.iter=100, ...) {
									
	mcmc.set <- get.tfr.mcmc(sim.dir)
	Eini <- mcmc.meta.ini.extra(mcmc.set, countries=countries, my.tfr.file=my.tfr.file, 
												burnin=burnin, verbose=verbose)
	meta <- Eini$meta
	if(length(Eini$index) <= 0) {
		cat('\nNothing to be done.\n')
		return(invisible(mcmc.set))
	}
	chain.ids <- names(mcmc.set$mcmc.list)
	mcthin <- 1
	if(verbose) cat('\n')
	for (chain in chain.ids) { # update meta in each chain
		if(verbose) cat('Updating meta in chain', chain, '\n')
		mcmc.set$mcmc.list[[chain]]$meta <- meta
		mcmc.set$mcmc.list[[chain]] <- mcmc.ini.extra(mcmc.set$mcmc.list[[chain]], countries=Eini$index,
												index.replace=Eini$index.replace)
		mcthin <- max(mcthin, mcmc.set$mcmc.list[[chain]]$thin)
	}
	if(length(Eini$index_DL) <= 0) {
		cat('\nNo DL countries or regions. Nothing to be done.\n')
		store.bayesTFR.meta.object(meta, meta$output.dir)
		mcmc.set$meta <- meta
		return(invisible(mcmc.set))
	}
	mcthin <- mcmc.set$mcmc.list[[1]]$thin
	total.iter <- mcmc.set$mcmc.list[[1]]$length - get.thinned.burnin(mcmc.set$mcmc.list[[1]], burnin)
	thin <- max(thin, mcthin)
	post.idx <- if (thin > mcthin) unique(round(seq(thin, total.iter, by=thin/mcthin)))
				else 1:total.iter
	if (!is.null(mcmc.set$mcmc.list[[1]]$rng.state)) .Random.seed <- mcmc.set$mcmc.list[[1]]$rng.state
	
	if (parallel) { # run chains in parallel
		if(is.null(nr.nodes)) nr.nodes<-length(chain.ids)
		chain.list <- bDem.performParallel(nr.nodes, chain.ids, mcmc.run.chain.extra, 
						initfun=init.nodes, mcmc.list=mcmc.set$mcmc.list, countries=Eini$index_DL, 
						posterior.sample=post.idx, iter=iter, burnin=burnin, verbose=verbose, verbose.iter=verbose.iter, ...)
		for (i in 1:length(chain.ids))
			mcmc.set$mcmc.list[[chain.ids[i]]] <- chain.list[[i]]
	} else { # run chains sequentially
		for (chain.id in chain.ids) {
			mcmc.set$mcmc.list[[chain.id]] <- mcmc.run.chain.extra(chain.id, mcmc.set$mcmc.list, 
												countries=Eini$index_DL, posterior.sample=post.idx, iter=iter,  
												burnin=burnin, verbose=verbose, verbose.iter=verbose.iter)
		}
	}
	store.bayesTFR.meta.object(meta, meta$output.dir)
	mcmc.set$meta <- meta
	cat('\n')
	invisible(mcmc.set)
}
	
mcmc.run.chain.extra <- function(chain.id, mcmc.list, countries, posterior.sample, 
												iter=NULL, burnin=2000, verbose=FALSE, verbose.iter=100) {
	cat('\n\nChain nr.', chain.id, '\n')
	if (verbose)
		cat('************\n')
	mcmc <- mcmc.list[[chain.id]]
		
	if (verbose) 
		cat('MCMC sampling for additional countries and regions.\n')

	mcmc <- tfr.mcmc.sampling.extra(mcmc, mcmc.list=mcmc.list, countries=countries, 
									posterior.sample=posterior.sample, 
									iter=iter, burnin=burnin, verbose=verbose, verbose.iter=verbose.iter)
	return(mcmc)
}

init.nodes <- function() {
	library(bayesTFR)
}

run.tfr.mcmc.subnat <- function(countries, my.tfr.file, sim.dir=file.path(getwd(), 'bayesTFR.output'), 
								nr.chains=3, iter=30000, thin=10,
								start.year=1950, present.year=2010, 								
								use.world.posterior = FALSE, world.sim.dir = sim.dir,
								post.burnin = 2000, buffer.size=100,
								seed = NULL, parallel=FALSE, nr.nodes=nr.chains, compression.type='None',
								auto.conf = list(max.loops=5, iter=30000, iter.incr=10000, nr.chains=3, thin=40, burnin=2000),
								verbose=FALSE, verbose.iter=100, ...) {
	if(is.null(world.sim.dir)) world.sim.dir <- sim.dir
	mcmc.set <- get.tfr.mcmc(world.sim.dir)
	meta <- mcmc.set$meta
	output.dir <- file.path(sim.dir, 'subnat')
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	# modify priors of hyper-paramters using posterior of a world-country simulation
	if(!use.world.posterior) {
		#meta$S.low <- 3.503
		#meta$S.up <- 5.474
		#meta$a.low <- 0.001
		#meta$a.up <- 0.057
		#meta$b.low <- 0.014
		#meta$b.up <- 0.062
		#meta$sigma0.low <- 0.221
		#meta$sigma0.up <- 0.301
		#meta$const.low <- 1.188
		#meta$const.up <- 1.59
		meta$chi0 <- -1.571
		meta$psi0 <- 0.126
		meta$alpha0.p <- c(-0.853, 0.459, 1.596)
		meta$delta0 <- 0.654
		meta$Triangle4.0 <- 0.51
		meta$delta4.0 <- 0.284
		meta$mean.eps.tau0 <- -0.272
		meta$sd.eps.tau0 <- 0.03
		meta$psi0r <- 0.869
		meta$nu.psi0 <- 57.441
		meta$delta0r <- 0.727
		meta$nu.delta0 <- 6.084
		meta$delta4.0r <- 0.869
		meta$nu4 <- 4.266
		meta$sd.eps.tau0r <- 0.364
		meta$nu.tau0 <- 142.914
	} else {
		library(MASS)
		if(has.tfr.prediction(sim.dir=world.sim.dir))
			m <- get.tfr.prediction(world.sim.dir)$mcmc.set
		else {
			l <- mcmc.set$mcmc.list[[1]]$finished.iter
			if(post.burnin > l) {
				post.burnin <- as.integer(l/2)
				warning('post.burnin larger than MCMC length. Adjusted to ', post.burnin)
			}
			m <- get.tfr.mcmc(world.sim.dir, burnin=post.burnin)
		}
		ml <- get.mcmc.list(m)
		
		# Uniform priors
		tr <- as.vector(get.tfr.parameter.traces(ml, 'S_sd'))
		meta$S.low <- quantile(tr, 0.001)
		meta$S.up <- quantile(tr, 0.995)
		tr <- as.vector(get.tfr.parameter.traces(ml, 'a_sd'))
		meta$a.low <- quantile(tr, 0.001)
		meta$a.up <- quantile(tr, 0.995)
		tr <- as.vector(get.tfr.parameter.traces(ml, 'b_sd'))
		meta$b.low <- quantile(tr, 0.001)
		meta$b.up <- quantile(tr, 0.999)
		tr <- as.vector(get.tfr.parameter.traces(ml, 'sigma0'))
		meta$sigma0.low <- quantile(tr, 0.001)
		meta$sigma0.up <- quantile(tr, 0.999)
		tr <- as.vector(get.tfr.parameter.traces(ml, 'const_sd'))
		meta$const.low <- quantile(tr, 0.001)
		meta$const.up <- quantile(tr, 0.999)
		
		# Normal priors
		tr <- as.vector(get.tfr.parameter.traces(ml, 'chi'))
		meta$chi0 <- mean(tr)
		meta$psi0 <- sd(tr)
		tr <- get.tfr.parameter.traces(ml, 'alpha')
		meta$alpha0.p <- apply(tr, 2, mean)
		meta$delta0 <- max(apply(tr, 2, sd))
		tr <- as.vector(get.tfr.parameter.traces(ml, 'Triangle4'))
		meta$Triangle4.0 <- mean(tr)
		meta$delta4.0 <- sd(tr)
		tr <- as.vector(get.tfr.parameter.traces(ml, 'mean_eps_tau'))
		meta$mean.eps.tau0 <- mean(tr)
		meta$sd.eps.tau0 <- sd(tr)
		
		# Gamma priors
		tr <- as.vector(get.tfr.parameter.traces(ml, 'psi'))
		gfit <- fitdistr(1/(tr)^2, 'gamma')
		meta$psi0r <- sqrt(gfit$estimate['rate']/gfit$estimate['shape'])
		meta$nu.psi0 <- gfit$estimate['shape']*2
		tr <- get.tfr.parameter.traces(ml, 'delta')
		dr <- nur <- c()
		for(i in 1:3) {
			gfit <- fitdistr(1/(tr[,i])^2, 'gamma')
			dr <- c(dr, sqrt(gfit$estimate['rate']/gfit$estimate['shape']))
			nur <- c(nur, gfit$estimate['shape']*2)
		}
		meta$delta0r <- min(dr)
		meta$nu.delta0 <- min(nur)
		tr <- as.vector(get.tfr.parameter.traces(ml, 'delta4'))
		gfit <- fitdistr(1/(tr)^2, 'gamma')
		meta$delta4.0r <- sqrt(gfit$estimate['rate']/gfit$estimate['shape'])
		meta$nu4 <- gfit$estimate['shape']*2
		tr <- as.vector(get.tfr.parameter.traces(ml, 'sd_eps_tau'))
		gfit <- fitdistr(1/(tr)^2, 'gamma')
		meta$sd.eps.tau0r <- sqrt(gfit$estimate['rate']/gfit$estimate['shape'])
		meta$nu.tau0 <- gfit$estimate['shape']*2
	}
	if(verbose)
		print(tfr.priors(meta))
	
	meta$buffer.size <- buffer.size
	meta$compression.type <- compression.type

	default.auto.conf <- formals(run.tfr.mcmc.subnat)$auto.conf
	for (par in names(default.auto.conf))
		if(is.null(auto.conf[[par]])) auto.conf[[par]] <- default.auto.conf[[par]]
	auto.run <- FALSE
	if(iter == 'auto') { # defaults for auto-run (includes convergence diagnostics)
		iter <- auto.conf$iter
		nr.chains <- auto.conf$nr.chains
		auto.run <- TRUE		
	}
	if(!is.null(seed)) set.seed(seed)
	# starting values (length of 1 or nr.chains)
	S.ini <- with(meta, ifelse(rep(nr.chains==1, nr.chains), 
					 		(S.low+S.up)/2, seq(S.low, to=S.up, length=nr.chains))) 
	a.ini <- with(meta, ifelse(rep(nr.chains==1, nr.chains), 
					 		(a.low+a.up)/2, seq(a.low, to=a.up, length=nr.chains))) 
	b.ini <- with(meta, ifelse(rep(nr.chains==1, nr.chains), 
					 		(b.low+b.up)/2, seq(b.low, to=b.up, length=nr.chains)))
	sigma0.ini <- with(meta, ifelse(rep(nr.chains==1, nr.chains), 
					 		(sigma0.low+sigma0.up)/2, 
					 		seq(sigma0.low, to=sigma0.up, length=nr.chains)))
	Triangle_c4.ini <- with(meta, ifelse(rep(nr.chains==1, nr.chains), 
					 		(Triangle_c4.low+Triangle_c4.up)/2, 
					 		seq(Triangle_c4.low+0.0001, to=Triangle_c4.up-0.0001, length=nr.chains))) 
	const.ini <- with(meta, ifelse(rep(nr.chains==1, nr.chains), 
					 		(const.low+const.up)/2, 
					 		seq(const.low, to=const.up, length=nr.chains)))
	gamma.ini <- 1
	# propagate initial values for all chains if needed
	for (var in c('S.ini', 'a.ini', 'b.ini', 'sigma0.ini', 'const.ini', 'gamma.ini', 'Triangle_c4.ini', 'iter')) {
		if (length(get(var)) < nr.chains) {
			if (length(get(var)) == 1) {
				assign(var, rep(get(var), nr.chains))
			} else {
			warning(var, ' has the wrong length. Either 1 or ', nr.chains, 
				' is allowed.\nValue set to ', get(var)[1], ' for all chains.')
			assign(var, rep(get(var)[1], nr.chains))
			}
		}
	}
	results <- list()
	for (country in countries) {
		country.obj <- get.country.object(country, mcmc.set$meta)	
		if(verbose) 
			cat('\nSimulating Phase II for', country.obj$name, '\n')							
		ini <- mcmc.meta.ini.subnat(meta, country=country.obj$code, my.tfr.file=my.tfr.file, 
									start.year=start.year, present.year=present.year, verbose=verbose)
		this.output.dir <- file.path(output.dir, paste0('c', country.obj$code))
		ini$output.dir <- this.output.dir
		if(file.exists(this.output.dir)) unlink(this.output.dir, recursive=TRUE)
		dir.create(this.output.dir)
		bayesTFR.mcmc.meta <- ini
		store.bayesTFR.meta.object(bayesTFR.mcmc.meta, this.output.dir)
		if (parallel) { # run chains in parallel
			chain.set <- bDem.performParallel(nr.nodes, 1:nr.chains, mcmc.run.chain, 
						initfun=init.nodes, meta=bayesTFR.mcmc.meta, 
						thin=thin, iter=iter, S.ini=S.ini, a.ini=a.ini,
                        b.ini=b.ini, sigma0.ini=sigma0.ini, Triangle_c4.ini=Triangle_c4.ini, const.ini=const.ini,
                        gamma.ini=gamma.ini, verbose=verbose, 
                        verbose.iter=verbose.iter, ...)
		} else { # run chains sequentially
			chain.set <- list()
			for (chain in 1:nr.chains) {
				chain.set[[chain]] <- mcmc.run.chain(chain, bayesTFR.mcmc.meta, thin=thin, 
					 	iter=iter, S.ini=S.ini, a.ini=a.ini, b.ini=b.ini, 
					 	sigma0.ini=sigma0.ini, Triangle_c4.ini=Triangle_c4.ini, const.ini=const.ini, 
					 	gamma.ini=gamma.ini, 
					 	verbose=verbose, verbose.iter=verbose.iter)
			}
		}
		names(chain.set) <- 1:nr.chains
		results[[as.character(country.obj$code)]] <- structure(list(meta=bayesTFR.mcmc.meta, mcmc.list=chain.set), 
																			class='bayesTFR.mcmc.set')
		cat('\nResults stored in', this.output.dir,'\n')
		if(auto.run) {
			diag <- try(tfr.diagnose(sim.dir=this.output.dir, keep.thin.mcmc=TRUE, 
						thin=auto.conf$thin, burnin=auto.conf$burnin,
						verbose=verbose))
			if(auto.conf$max.loops>1) {
				for(loop in 2:auto.conf$max.loops) {
					if(!inherits(diag, "try-error") && has.mcmc.converged(diag)) break
					results[[as.character(country.obj$code)]] <- continue.tfr.mcmc(
								iter=auto.conf$iter.incr, output.dir=this.output.dir, nr.nodes=nr.nodes,
								parallel=parallel, verbose=verbose, verbose.iter=verbose.iter)
					diag <- try(tfr.diagnose(sim.dir=this.output.dir, keep.thin.mcmc=TRUE, 
								thin=auto.conf$thin, burnin=auto.conf$burnin,
								verbose=verbose))
				}
			}
		}
	}
	if (verbose)
		cat('\nSimulation successfully finished!!!\n')
	invisible(results)
		
}

set.default.cltype <- function() {
	if(!is.element(snow::getClusterOption("type"), c("MPI", "SOCK"))) 
		snow::setDefaultClusterOptions(type="SOCK")
}

bDem.performParallel <- function(...) {
	set.default.cltype()
	snowFT::performParallel(...)
}
