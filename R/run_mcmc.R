run.tfr.mcmc <- function(nr.chains=3, iter=62000, output.dir=file.path(getwd(), 'bayesTFR.output'), 
						thin=1, replace.output=FALSE, annual = FALSE,
						# meta parameters
						start.year=1950, present.year=2020, wpp.year=2019,
						my.tfr.file = NULL, my.locations.file = NULL, buffer.size=100,
					 	U.c.low=5.5, U.up=8.8, U.width=3,
					 	mean.eps.tau0 = -0.25, sd.eps.tau0 = 0.4, nu.tau0 = 2,                                                
        				Triangle_c4.low = 1, Triangle_c4.up = 2.5,
        				Triangle_c4.trans.width=2, Triangle4.0 = 0.3, 
        				delta4.0 = 0.8, nu4 = 2,
					 	S.low=3.5, S.up=6.5, S.width=0.5,
					 	a.low=0, a.up=0.2, a.width=0.02,
					 	b.low=a.low, b.up=a.up, b.width=0.05,
					 	sigma0.low=if(annual) 0.0045 else 0.01, sigma0.up=0.6, sigma0.width=0.1,
					 	sigma0.min=0.04, 
					 	const.low=0.8, const.up=2, const.width=0.3,
					 	d.low=0.05, d.up=0.5, d.trans.width=1,
					 	chi0=-1.5, psi0=0.6, nu.psi0=2,
					 	alpha0.p=c(-1, 0.5, 1.5), delta0=1, nu.delta0=2,
					 	dl.p1=9, dl.p2=9,
						phase3.parameter=NULL,
						# starting values (length of 1 or nr.chains)
						S.ini=NULL, a.ini=NULL, b.ini=NULL, 
					 	sigma0.ini=NULL, Triangle_c4.ini=NULL, const.ini=NULL, gamma.ini=1, 
						phase3.starting.values=NULL,
					 	proposal_cov_gammas = NULL, # should be a list with elements 'values' and 'country_codes'
					 	seed = NULL, parallel=FALSE, nr.nodes=nr.chains, 
					 	save.all.parameters = FALSE, compression.type='None',
					 	auto.conf = list(max.loops=5, iter=62000, iter.incr=10000, nr.chains=3, thin=80, burnin=2000),
						verbose=FALSE, verbose.iter = 10, uncertainty = FALSE, 
						my.tfr.raw.file=NULL, 
						iso.unbiased=NULL, covariates=c('source', 'method'), cont_covariates=NULL, ar.phase2=FALSE, ...) 
{
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
					 		seq(max(sigma0.low, 0.1), to=sigma0.up, length=nr.chains))
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
						wpp.year=wpp.year, annual.simulation = annual,
						my.tfr.file = my.tfr.file, my.locations.file=my.locations.file,
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
					 	auto.conf=auto.conf, package.version = packageVersion("bayesTFR"),
						verbose=verbose, uncertainty=uncertainty, my.tfr.raw.file=my.tfr.raw.file, ar.phase2=ar.phase2)
	if (uncertainty)
	{
	  bayesTFR.mcmc.meta[["covariates"]] <- covariates
	  bayesTFR.mcmc.meta[["cont_covariates"]] <- cont_covariates
	}
	store.bayesTFR.meta.object(bayesTFR.mcmc.meta, output.dir)
	
	starting.values <- NULL
	if (uncertainty)
	{
	  get.init.values <- function(range) {
	    ifelse(rep(nr.chains==1, nr.chains), sum(range)/2, 
	           #seq(range[1], to=range[2], length=nr.chains)
	           runif(nr.chains, range[1], range[2])
	    )
	  }
	  
	  dir.create(file.path(output.dir, 'phaseIII'))
	  if ('mu.prior.range' %in% names(phase3.parameter)) {mu.prior.range <- phase3.parameter[['mu.prior.range']]}
	  else {mu.prior.range <- c(0, 2.1)}
	  if ('rho.prior.range' %in% names(phase3.parameter)) {rho.prior.range <- phase3.parameter[['rho.prior.range']]}
	  else {rho.prior.range <- c(0,1-.Machine$double.xmin)}
	  if ('sigma.mu.prior.range' %in% names(phase3.parameter)) {sigma.mu.prior.range <- phase3.parameter[['sigma.mu.prior.range']]}
	  else {sigma.mu.prior.range <- c(1e-5,0.318)}
	  if ('sigma.rho.prior.range' %in% names(phase3.parameter)) {sigma.rho.prior.range <- phase3.parameter[['sigma.rho.prior.range']]}
	  else {sigma.rho.prior.range <- c(1e-5,0.289)}
	  if ('sigma.eps.prior.range' %in% names(phase3.parameter)) {sigma.eps.prior.range <- phase3.parameter[['sigma.eps.prior.range']]}
	  else {sigma.eps.prior.range <- c(1e-5, 0.5)}
	  
	  for (varname in c('mu', 'rho', 'sigma.mu', 'sigma.rho', 'sigma.eps'))
	  {
	    assign(paste0(varname, '.ini.range'), get(paste0(varname, '.prior.range')))
	    if (paste0(varname, '.ini') %in% names(phase3.starting.values))
	    {
	      assign(paste0(varname, '.ini'), phase3.starting.values[[paste0(varname, '.ini')]])
	    }
	    else
	    {
	      assign(paste0(varname, '.ini'), get.init.values(get(paste0(varname, '.prior.range'))))
	    }
	  }
	  c.index <- 1:get.nr.countries(bayesTFR.mcmc.meta)
	  meta <- structure(list(nr.chains=nr.chains,
	                         my.tfr.file=my.tfr.file, output.dir=output.dir,
	                         phase=3, id_phase3 = which(bayesTFR.mcmc.meta$lambda_c[c.index] < bayesTFR.mcmc.meta$T_end_c[c.index]),
	                         nr.countries=sum(bayesTFR.mcmc.meta$lambda_c[c.index] < bayesTFR.mcmc.meta$T_end_c[c.index]),
	                         mu.prior.range=mu.prior.range, rho.prior.range=rho.prior.range,
	                         sigma.mu.prior.range=sigma.mu.prior.range, 
	                         sigma.rho.prior.range=sigma.rho.prior.range,
	                         sigma.eps.prior.range=sigma.eps.prior.range,
	                         mu.ini = get("mu.ini"), mu.ini.range=get("mu.ini.range"), 
	                         rho.ini=get("rho.ini"), rho.ini.range=get("rho.ini.range"), 
	                         sigma.mu.ini=get("sigma.mu.ini"), sigma.mu.ini.range=get("sigma.mu.ini.range"),
	                         sigma.rho.ini=get("sigma.rho.ini"), sigma.rho.ini.range=get("sigma.rho.ini.range"),
	                         sigma.eps.ini=get("sigma.eps.ini"), sigma.eps.ini.range=get("sigma.eps.ini.range"),
	                         compression.type=compression.type, buffer.size=buffer.size, auto.conf=auto.conf
	  ), class='bayesTFR.mcmc.meta')	
	  store.bayesTFR.meta.object(meta, file.path(output.dir, 'phaseIII'))
	  if(meta$nr.countries <= 0) run.phase3 <- FALSE else run.phase3 <- TRUE
	  for (name in names(meta))
	  {
	    if (!(name %in% names(bayesTFR.mcmc.meta)))
	    {
	      bayesTFR.mcmc.meta[[name]] <- meta[[name]]
	    }
	  }
	  bayesTFR.mcmc.meta$run.phase3 <- run.phase3
	  bayesTFR.mcmc.meta$parent <- bayesTFR.mcmc.meta
	  for (var in c('mu.ini', 'rho.ini', 'sigma.mu.ini', 'sigma.rho.ini', 'sigma.eps.ini', 'iter')) {
	    if (length(get(var)) < nr.chains) 
	      assign(var, rep(get(var), nr.chains)[1:nr.chains])
	    if (var != 'iter') starting.values[[var]] <- get(var)
	  }
	}
	
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
						initfun=init.nodes, seed = seed, meta=bayesTFR.mcmc.meta, 
						thin=thin, starting.values=starting.values, iter=iter, S.ini=S.ini, a.ini=a.ini,
                        b.ini=b.ini, sigma0.ini=sigma0.ini, Triangle_c4.ini=Triangle_c4.ini, const.ini=const.ini,
                        gamma.ini=gamma.ini, save.all.parameters=save.all.parameters, verbose=verbose, 
                        verbose.iter=verbose.iter, uncertainty=uncertainty, iso.unbiased=iso.unbiased, 
						            covariates=covariates, cont_covariates=cont_covariates, ...)
	} else { # run chains sequentially
		chain.set <- list()
		for (chain in 1:nr.chains) {
			chain.set[[chain]] <- mcmc.run.chain(chain, bayesTFR.mcmc.meta, thin=thin, starting.values=starting.values, 
					 	iter=iter, S.ini=S.ini, a.ini=a.ini, b.ini=b.ini, 
					 	sigma0.ini=sigma0.ini, Triangle_c4.ini=Triangle_c4.ini, const.ini=const.ini, 
					 	gamma.ini=gamma.ini, save.all.parameters=save.all.parameters,
					 	verbose=verbose, verbose.iter=verbose.iter, uncertainty=uncertainty, iso.unbiased=iso.unbiased, 
					 	covariates=covariates, cont_covariates=cont_covariates)
		}
	}
	names(chain.set) <- 1:nr.chains
	mcmc.set <- structure(list(meta=chain.set[[1]]$meta, mcmc.list=chain.set), class='bayesTFR.mcmc.set')
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


mcmc.run.chain <- function(chain.id, meta, thin=1, iter=100, starting.values=NULL, 
							S.ini, a.ini, b.ini, sigma0.ini, Triangle_c4.ini, const.ini, gamma.ini=1,
							save.all.parameters=FALSE,
							verbose=FALSE, verbose.iter=10, uncertainty=FALSE, iso.unbiased=NULL, 
							covariates=c('source', 'method'), cont_covariates=NULL) {
								
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
	                 verbose=verbose, uncertainty=uncertainty, iso.unbiased=iso.unbiased, 
	                 covariates=covariates, cont_covariates=NULL)
	if (uncertainty)
	{
	  this.sv <- list()
	  for(var in names(starting.values)) {
	    this.sv[[var]] <- starting.values[[var]][chain.id]
	  }
	  mcmc3 <- do.call('mcmc3.ini', c(list(chain.id, meta, iter=iter[chain.id], thin=thin, starting.values=this.sv) ))
	  for (name in names(mcmc3))
	  {
	    if (!(name %in% names(mcmc)))
	    {
	      mcmc[[name]] <- mcmc3[[name]]
	    }
	    else if (name == 'meta')
	    {
	      mcmc[['meta3']] <- mcmc3[['meta']]
	    }
	  }
	}
	
	if (verbose) 
		cat('Store initial values into ', mcmc$output.dir, '\n')
	
	store.mcmc(mcmc, append=FALSE, flush.buffer=TRUE, verbose=verbose)
	if (uncertainty) store.mcmc3(mcmc, append=FALSE, flush.buffer=TRUE, verbose=verbose)
	
	if (verbose) 
		cat('Start sampling -', mcmc$iter, 'iterations in total.\n')
	mcmc <- tfr.mcmc.sampling(mcmc, thin=thin, verbose=verbose, verbose.iter=verbose.iter, uncertainty=uncertainty)
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
	uncertainty <- mcmc$uncertainty
	if (verbose) 
		cat('Continue sampling -', iter, 'additional iterations,', mcmc$iter, 'iterations in total.\n')

	mcmc <- tfr.mcmc.sampling(mcmc, thin=mcmc$thin, start.iter=mcmc$finished.iter+1, verbose=verbose, verbose.iter=verbose.iter, uncertainty=uncertainty)
	return(mcmc)
}

run.tfr.mcmc.extra <- function(sim.dir=file.path(getwd(), 'bayesTFR.output'), 
								countries = NULL, my.tfr.file = NULL, iter = NULL,
								thin=1, thin.extra=1, burnin=2000, parallel=FALSE, nr.nodes=NULL, 
								my.locations.file = NULL,
								verbose=FALSE, verbose.iter=100, uncertainty=FALSE, 
								my.tfr.raw.file=NULL, 
								iso.unbiased=NULL, covariates=c('source', 'method'), cont_covariates=NULL, ...) {
  mcmc.set <- get.tfr.mcmc(sim.dir)
  meta.old <- mcmc.set$meta
  
  if(is.null(covariates) && is.null(cont_covariates) && uncertainty)
  {
    covariates <- meta.old[['covariates']]
    cont_covariates <- meta.old[['cont_covariates']]
  }
  
  Eini <- mcmc.meta.ini.extra(mcmc.set, countries=countries, my.tfr.file=my.tfr.file, 
												my.locations.file=my.locations.file, burnin=burnin, verbose=verbose, uncertainty=uncertainty, 
												my.tfr.raw.file=my.tfr.raw.file)
	if(length(Eini$index) <= 0) {
		cat('\nNothing to be done.\n')
		return(invisible(mcmc.set))
	}
	if (uncertainty && has.tfr3.mcmc(sim.dir))
	{
	  mcmc3.set <- get.tfr3.mcmc(sim.dir)
	  Eini$meta[['id_phase3']] <- intersect(mcmc3.set$meta$id_phase3, which(mcmc.set$meta$regions$country_code %in% countries))
	  for (par.name in tfr3.parameter.names())
	  {
	    for (suffix in c('prior.range', 'ini', 'ini.range'))
	    {
	      Eini$meta[[paste0(par.name, '.', suffix)]] <- mcmc3.set$meta[[paste0(par.name, '.', suffix)]]
	    }
	  }
	  Eini$meta$parent <- mcmc3.set$meta$parent
	}
	
	meta <- Eini$meta
	chain.ids <- names(mcmc.set$mcmc.list)
	mcthin <- 1
	
	if(verbose) cat('\n')
	for (chain in chain.ids) { # update meta in each chain
		if(verbose) cat('Updating meta in chain', chain, '\n')
		mcmc.set$mcmc.list[[chain]]$meta <- meta
		mcmc.set$mcmc.list[[chain]] <- mcmc.ini.extra(mcmc.set$mcmc.list[[chain]], countries=Eini$index,
												index.replace=Eini$index.replace)
		if (uncertainty)
		{
		  mcmc.set$mcmc.list[[chain]]$meta$nr.countries <- length(Eini$meta[['id_phase3']])
		  for(varname in c('mu', 'rho')) {
		    var <- paste(varname, 'c', sep='.')
		    range.var <- paste(varname,'ini.range', sep='.')
		    mcmc.set$mcmc.list[[chain]][[var]] <- runif(length(Eini$meta[['id_phase3']]), 
		                                                mcmc.set$mcmc.list[[chain]]$meta[[range.var]][1], 
		                                                mcmc.set$mcmc.list[[chain]]$meta[[range.var]][2])
		  }
		  if (length(Eini$meta[['id_phase3']]) > 0)
		  {
		    for (country in 1:length(Eini$meta[['id_phase3']]))
		      mcmc.set$mcmc.list[[chain]]$observations[[country]] <- 
		        mcmc.set$mcmc.list[[chain]]$meta$tfr_all[meta$lambda_c[meta$id_phase3[country]]:meta$T_end, meta$id_phase3[country]]
		  }
		}
		mcthin <- max(mcthin, mcmc.set$mcmc.list[[chain]]$thin)
		if(uncertainty)
		{
		  mcmc.set$mcmc.list[[chain]] <- get.obs.estimate.diff.original(mcmc.set$mcmc.list[[chain]])
		  mcmc.set$mcmc.list[[chain]] <- estimate.bias.sd.original(mcmc.set$mcmc.list[[chain]], iso.unbiased, covariates, cont_covariates)
		  mcmc.set$mcmc.list[[chain]]$eps_unc <- list()
		  if (is.null(mcmc.set$mcmc.list[[chain]]$meta$raw_data_extra)) mcmc.set$mcmc.list[[chain]]$meta$raw_data_extra <- list()
		  for (country in countries)
		  {
		    df_country <- mcmc.set$mcmc.list[[chain]]$meta$raw_data.original
		    country.obj <- get.country.object(country, meta.old)
		    if (!is.null(country.obj$index))
		    {
		      df_country <- df_country[df_country$country_code == country,]
		      mcmc.set$mcmc.list[[chain]]$meta$raw_data_extra[[country.obj$index]] <- df_country
		      mcmc.set$mcmc.list[[chain]]$eps_unc[[country.obj$index]] <- df_country$eps
		    }
		  }
		}
	}
	
	meta <- mcmc.set$mcmc.list[[1]]$meta
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
						posterior.sample=post.idx, iter=iter, thin = thin.extra, burnin=burnin, verbose=verbose, verbose.iter=verbose.iter, 
						uncertainty=uncertainty, ...)
		for (i in 1:length(chain.ids))
			mcmc.set$mcmc.list[[chain.ids[i]]] <- chain.list[[i]]
	} else { # run chains sequentially
		for (chain.id in chain.ids) {
			mcmc.set$mcmc.list[[chain.id]] <- mcmc.run.chain.extra(chain.id, mcmc.set$mcmc.list, 
												countries=Eini$index_DL, posterior.sample=post.idx, iter=iter, thin = thin.extra,
												burnin=burnin, verbose=verbose, verbose.iter=verbose.iter, uncertainty=uncertainty)
		}
	}
	if(uncertainty)
	{
	  if (!dir.exists(file.path(meta$output.dir, 'extra.meta'))) dir.create(file.path(meta$output.dir, 'extra.meta'))
	  if (!dir.exists(file.path(meta$output.dir, 'extra.meta', countries[1]))) dir.create(file.path(meta$output.dir, 'extra.meta', countries[1]))
	  store.bayesTFR.meta.object(meta, file.path(meta$output.dir, 'extra.meta', countries[1]))
	  for (name in c("country.ind.by.year", "ind.by.year", "id_phase1_by_year", "id_phase2_by_year", "id_phase3_by_year", 
	                 "id_phase3", "nr.countries"))
	  {
	    meta[[name]] <- meta.old[[name]]
	  }
	  if (is.null(meta[['extra']])) meta[['extra']] <- c()
	  if (is.null(meta[['extra_iter']])) meta[['extra_iter']] <- numeric(get.nrest.countries(meta.old))
	  if (is.null(meta[['extra_thin']])) meta[['extra_thin']] <- numeric(get.nrest.countries(meta.old))
	  if (is.null(meta[['extra_covariates']])) meta[['extra_covariates']] <- list()
	  if (is.null(meta[['extra_cont_covariates']])) meta[['extra_cont_covariates']] <- list()
	  for (country in countries)
	  {
	    country.idx <- get.country.object(country, meta.old)$index
	    if (!is.null(country.idx)) 
	    {
	      meta[['extra']] <- c(meta[['extra']], country.idx)
	      meta[['extra_iter']][country.idx] <- iter
	      meta[['extra_thin']][country.idx] <- thin.extra
	      meta[['extra_covariates']][[country.idx]] <- covariates
	      meta[['extra_cont_covariates']][[country.idx]] <- cont_covariates
	    }
	  }
	  meta[['extra']] <- sort(unique(meta[['extra']]))
	}
	store.bayesTFR.meta.object(meta, file.path(meta$output.dir))
	mcmc.set$meta <- meta
	cat('\n')
	invisible(mcmc.set)
}
	
mcmc.run.chain.extra <- function(chain.id, mcmc.list, countries, posterior.sample, 
												iter=NULL, burnin=2000, thin=1, verbose=FALSE, verbose.iter=100, uncertainty=FALSE) {
	cat('\n\nChain nr.', chain.id, '\n')
	if (verbose)
		cat('************\n')
	mcmc <- mcmc.list[[chain.id]]
	mcmc$uncertainty <- uncertainty
	
	if (verbose) 
		cat('MCMC sampling for additional countries and regions.\n')

	mcmc <- tfr.mcmc.sampling.extra(mcmc, mcmc.list=mcmc.list, countries=countries, 
									posterior.sample=posterior.sample, 
									iter=iter, burnin=burnin, thin=thin, verbose=verbose, verbose.iter=verbose.iter, uncertainty=uncertainty)
	return(mcmc)
}


init.nodes <- function() {
	library(bayesTFR)
}

set.default.cltype <- function() {
#	if(!is.element(snow::getClusterOption("type"), c("MPI", "SOCK"))) 
#		snow::setDefaultClusterOptions(type="SOCK")
}

bDem.performParallel <- function(..., cltype='SOCK') {
	set.default.cltype()
	snowFT::performParallel(..., cltype=cltype)
}
