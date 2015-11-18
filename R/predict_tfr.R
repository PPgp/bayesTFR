tfr.predict <- function(mcmc.set=NULL, end.year=2100,
						sim.dir=file.path(getwd(), 'bayesTFR.output'),
						replace.output=FALSE,
						start.year=NULL, nr.traj = NULL, thin = NULL, burnin=2000, 
						use.diagnostics=FALSE,
						use.tfr3=TRUE, mcmc3.set=NULL, burnin3=10000,
						mu=2.1, rho=0.8859, sigmaAR1=0.1016,
						use.phase3.from = NULL,
						use.correlation=FALSE,
						save.as.ascii=1000, output.dir = NULL,
						low.memory=TRUE,
						seed=NULL, verbose=TRUE, ...) {
	if(!is.null(mcmc.set)) {
		if (class(mcmc.set) != 'bayesTFR.mcmc.set') {
			stop('Wrong type of mcmc.set. Must be of type bayesTFR.mcmc.set.')
			}
		mcmc.set.passed <- TRUE
	} else {		
		mcmc.set <- get.tfr.mcmc(sim.dir, low.memory=low.memory, verbose=verbose)
		mcmc.set.passed <- FALSE
	}
	has.phase3 <- FALSE
	if(use.tfr3) {
		has.phase3 <- !is.null(mcmc3.set)
		if(!has.phase3) has.phase3 <- has.tfr3.mcmc(mcmc.set$meta$output.dir)
		if(!has.phase3 && !is.null(use.phase3.from) && has.tfr3.mcmc(use.phase3.from)) has.phase3 <- TRUE
		if(!has.phase3)
			warning('No Phase III MCMCs available. Switching to constant AR(1) parameters.', immediate. = TRUE)
	}
	if(!has.phase3) {
		if (is.null(rho) || is.na(rho) || is.null(sigmaAR1) || is.na(sigmaAR1)) {
			res <- get.ar1.parameters(mu = mu, mcmc.set$meta)
			if (is.null(rho) || is.na(rho)) 
				rho <- res$rho
			if(is.null(sigmaAR1) || is.na(sigmaAR1))
				sigmaAR1 <- res$sigmaAR1
		}
	}
	if(verbose) {
		if(has.phase3) cat('\nAR(1) simulated using phase III MCMCs.\n')
		else cat('\nAR(1) parameters for all countries: mu=', mu, ', rho=', rho, ', sigma=', sigmaAR1, '\n')
	}
	if(!is.null(seed)) set.seed(seed)
	
	# Get argument settings from existing convergence diagnostics
	if(use.diagnostics) {
		nrtraj.burnin <- get.burnin.nrtraj.from.diagnostics(mcmc.set)
		nr.traj <- nrtraj.burnin[1]
		burnin <- nrtraj.burnin[2]
		if(verbose)
			cat('\nUsing convergence settings: nr.traj=', nr.traj, ', burnin=', burnin, '\n')
	}
	invisible(make.tfr.prediction(mcmc.set, end.year=end.year, replace.output=replace.output,  
					start.year=start.year, nr.traj=nr.traj, burnin=burnin, thin=thin, use.tfr3=has.phase3, 
					mcmc3.set=mcmc3.set, burnin3=burnin3,
					mu=mu, rho=rho,  sigmaAR1 = sigmaAR1, use.phase3.from=use.phase3.from, 
					use.correlation=use.correlation,
					save.as.ascii=save.as.ascii,
					output.dir=output.dir, verbose=verbose, ...))			
}

get.burnin.nrtraj.from.diagnostics <- function(mcmc.set) {
	diag.list <- get.tfr.convergence.all(mcmc.set$meta$output.dir)
	ldiag <- length(diag.list)
	if (ldiag == 0) stop('There is no diagnostics available. Use manual settings of "nr.traj" or "thin".')
	use.nr.traj <- use.burnin <- rep(NA, ldiag)
	for(idiag in 1:ldiag) {
		if (has.mcmc.converged(diag.list[[idiag]])) {
			use.nr.traj[idiag] <- diag.list[[idiag]]$use.nr.traj
			use.burnin[idiag] <- diag.list[[idiag]]$burnin
		}
	}
	if(all(is.na(use.nr.traj)))
		stop('There is no diagnostics indicating convergence of the MCMCs. Use manual settings of "nr.traj" or "thin".')
	# Try to select those that suggest nr.traj >= 2000 (take the minimum of those)
	traj.is.notna <- !is.na(use.nr.traj)
	larger2T <- traj.is.notna & use.nr.traj>=2000
	nr.traj.idx <- if(sum(larger2T)>0) (1:ldiag)[larger2T][which.min(use.nr.traj[larger2T])] 
					else (1:ldiag)[traj.is.notna][which.max(use.nr.traj[traj.is.notna])]
	return(c(use.nr.traj[nr.traj.idx], use.burnin[nr.traj.idx]))
}

tfr.predict.extra <- function(sim.dir=file.path(getwd(), 'bayesTFR.output'), 
					prediction.dir=sim.dir, 
					countries = NULL, save.as.ascii=1000, verbose=TRUE) {
	# Run prediction for given countries/regions (as codes). If they are not given it will be set to countries 
	# for which there are MCMC results but no prediction.
	# It is to be used after running run.tfr.mcmc.extra
	
	mcmc.set <- get.tfr.mcmc(sim.dir)
	if(is.null(mcmc.set))
		stop('Error in "sim.dir" argument.')
	pred <- get.tfr.prediction(sim.dir=prediction.dir)
	if(is.null(pred))
		stop('Error in "prediction.dir" argument.')
	if(length(setdiff(pred$mcmc.set$meta$regions$country_code, mcmc.set$meta$regions$country_code)) > 0)
		stop('Prediction is inconsistent with the mcmc results. Use tfr.predict.')
	if(is.null(countries)) {
		countries.idx <- (1:mcmc.set$meta$nr_countries)[!is.element(mcmc.set$meta$regions$country_code, 
												pred$mcmc.set$meta$regions$country_code)]
	} else {
		countries.idx <- (1:mcmc.set$meta$nr_countries)[is.element(mcmc.set$meta$regions$country_code,
												countries)]
	}
	if(length(countries.idx) == 0) {
		cat('\nNothing to be done.\n')
		return(invisible(pred))	
	}
	use.tfr3 <- pred$use.tfr3
	if(pred$use.tfr3 && !has.tfr3.mcmc(sim.dir)) {
		warning('Prediction used BHM for phase III TFR but the MCMCs are not longer available. Switching to constant AR(1) parameters.')
		use.tfr3 <- FALSE
	}
	new.pred <- make.tfr.prediction(mcmc.set, start.year=pred$start.year, end.year=pred$end.year, replace.output=FALSE,
									nr.traj=pred$nr.traj, burnin=pred$burnin, use.tfr3=use.tfr3, burnin3=pred$burnin3,
									use.correlation=pred$use.correlation, mu=pred$mu, rho=pred$rho, sigmaAR1=pred$sigmaAR1, 
									countries=countries.idx, save.as.ascii=0, output.dir=prediction.dir,
									force.creating.thinned.mcmc=TRUE,
									write.summary.files=FALSE, write.trajectories=TRUE, verbose=verbose)
									
	# merge the two predictions
	code.other.countries <- setdiff(pred$mcmc.set$meta$regions$country_code, 
									mcmc.set$meta$regions$country_code[countries.idx])
	idx.pred.others <- (1:pred$mcmc.set$meta$nr_countries)[is.element(pred$mcmc.set$meta$regions$country_code, 
												code.other.countries)]
	idx.other.countries <- (1:mcmc.set$meta$nr_countries)[is.element(mcmc.set$meta$regions$country_code,
												code.other.countries)]
												
	prev.pred <- pred
	pred$quantiles <- new.pred$quantiles
	pred$quantiles[idx.other.countries,,] <- prev.pred$quantiles[idx.pred.others,,]
	
	pred$traj.mean.sd <- new.pred$traj.mean.sd
	pred$traj.mean.sd[idx.other.countries,,] <- prev.pred$traj.mean.sd[idx.pred.others,,]
	
	pred$tfr_matrix_reconstructed <- new.pred$tfr_matrix_reconstructed
	pred$tfr_matrix_reconstructed[,idx.other.countries] <- prev.pred$tfr_matrix_reconstructed[,idx.pred.others]
	
	pred$mcmc.set <- new.pred$mcmc.set
	
	# save updated prediction, convert trajectories and create summary files
	bayesTFR.prediction <- pred
	prediction.file <- file.path(pred$output.dir, 'prediction.rda')
	save(bayesTFR.prediction, file=prediction.file)
	
	do.convert.trajectories(pred=bayesTFR.prediction, n=save.as.ascii, output.dir=pred$output.dir, 
							verbose=verbose)
	tfr.write.projection.summary.and.parameters(pred=bayesTFR.prediction, output.dir=pred$output.dir)
	
	cat('\nPrediction stored into', pred$output.dir, '\n')
	invisible(bayesTFR.prediction)
}

estimate.scale.and.sd <- function(wmeta, meta, country.index, reg.index, intercept=FALSE) {
	wtfr <- get.observed.tfr(country.index, wmeta, 'tfr_matrix_all')
	regtfr <- get.observed.tfr(reg.index, meta, 'tfr_matrix_all')
	# get only common years that are not NA
	regtfr <- regtfr[!is.na(regtfr)]
	regtfr <- regtfr[names(regtfr) %in%names(wtfr)]
    wtfr <- wtfr[names(wtfr) %in%names(regtfr)]
    f <- "regtfr ~ wtfr"
    if (!intercept) f <- paste0(f, ' - 1')
    lmfit <- lm(formula(f))
    if(intercept) {
    	alpha <- lmfit$coef[1]
    	scale <- lmfit$coef[2]
    } else {
    	alpha <- 0
    	scale <- lmfit$coef[1]
    }
    return(list(scale=scale, alpha=alpha, sd=sqrt(mean(lmfit$residuals^2)), last.tfr=regtfr[length(regtfr)]))
}

estimate.scale <- function(wmeta, meta, country.index, reg.index) {
	wtfr <- get.observed.tfr(country.index, wmeta, 'tfr_matrix_all')
	regtfr <- get.observed.tfr(reg.index, meta, 'tfr_matrix_all')
	# get only common years that are not NA
	regtfr <- regtfr[!is.na(regtfr)]
	regtfr <- regtfr[names(regtfr) %in%names(wtfr)]
    wtfr <- wtfr[names(wtfr) %in%names(regtfr)]
    lmfit <- lm(regtfr ~ -1 + wtfr)
    scale <- lmfit$coef[1]
   
    return(list(scale=scale, last.tfr=regtfr[length(regtfr)], fit=lmfit))
}

tfr.predict.subnat <- function(countries, my.tfr.file, sim.dir=file.path(getwd(), 'bayesTFR.output'),
								end.year=2100, output.dir = NULL, nr.traj=NULL, seed = NULL,  
								sigma3=0.115, sigma3.default=0.115, use.lmfit=c(FALSE, TRUE), verbose=TRUE) {
	wpred <- get.tfr.prediction(sim.dir)
	wmeta <- wpred$mcmc.set$meta
	if(!is.null(seed)) set.seed(seed)
	if (is.null(output.dir)) output.dir <- wmeta$output.dir
	quantiles.to.keep <- as.numeric(dimnames(wpred$quantiles)[[2]])
	nr.project <- min(length(seq(wmeta$start.year+5, end.year, by=5)), wpred$nr.projections)
	
	result <- list()
	BHMp2 <- get.tfr.parameter.traces(wpred$mcmc.set$mcmc.list, 
							c('sigma0', 'a_sd', 'b_sd', 'S_sd', 'const_sd'), burnin=0)
	if(is.null(nr.traj)) nr.traj <- nrow(BHMp2)
	thinning.index <- unique(round(seq(1, nrow(BHMp2), length=nr.traj)))
	nr.trajs <- length(thinning.index)
	BHMp2 <- BHMp2[thinning.index,]
	if(length(sigma3)==1 && length(countries) > 1) {
		sigma3 <- rep(sigma3, length(countries))
		names(sigma3) <- countries
	}
	for (country in countries) {
		country.obj <- get.country.object(country, wmeta)
		if(verbose) 
			cat('\nPredicting TFR for ', country.obj$name, '\n')
		start.year <- if(is.null(wpred$start.year)) as.integer(dimnames(wpred$tfr_matrix_reconstructed)[[1]][wpred$present.year.index])+5
						else wpred$start.year
		meta <- mcmc.meta.ini.subnat(wmeta, country=country.obj$code, my.tfr.file=my.tfr.file, 
									start.year=1750, present.year=start.year-5, verbose=verbose)
		this.output.dir <- file.path(output.dir, 'subnat', paste0('c', country.obj$code))
		outdir <- file.path(this.output.dir, 'predictions')
		meta$output.dir <- this.output.dir
		
		if(file.exists(this.output.dir)) unlink(this.output.dir, recursive=TRUE)
		dir.create(outdir, recursive=TRUE)
		bayesTFR.mcmc.meta <- meta
		store.bayesTFR.meta.object(bayesTFR.mcmc.meta, this.output.dir)
		wtrajs <- get.tfr.trajectories(wpred, country.obj$code)		
		wtrajs <- wtrajs[1:(nr.project+1),thinning.index]
		nr.reg <- get.nr.countries(meta)
		PIs_cqp <- array(NA, c(nr.reg, length(quantiles.to.keep), nrow(wtrajs)),
						dimnames=list(meta$regions$country_code, dimnames(wpred$quantiles)[[2]], dimnames(wtrajs)[[1]]))
		mean_sd <- array(NA, c(nr.reg, 2, nrow(wtrajs)))
		meta$T_end_c <- rep(nrow(meta$tfr_matrix_all), nr.reg)
		this.sigma3 <- if(as.character(country) %in% names(sigma3)) sigma3[as.character(country)] else sigma3.default
		if(verbose) cat('\nPhase III SD = ', this.sigma3)
		for(region in 1:nr.reg) {
			reg.obj <- get.country.object(region, meta, index=TRUE)			
			scale.res <- estimate.scale(wmeta, meta, country.obj$index, reg.obj$index)
			scale <- scale.res$scale
			tfr.pred <- matrix(NA, nrow=nr.project+1, ncol=ncol(wtrajs))
			tfr.pred[1,] <- scale.res$last.tfr
			is.reg.phase3 <- !(reg.obj$index %in% meta$id_Tistau) & (length(meta$lambda_c[region]:meta$T_end_c[region]) > 1)
			for(s in 1:ncol(tfr.pred)) { # iterate over trajectories
				is.in.phase3 <- is.reg.phase3
				for(year in 2:(nr.project+1)) {
					if(!is.in.phase3 & year > 2) { # check if now in phase 3
						end.phase2 <- find.lambda.for.one.country(tfr.pred[1:(year-1),s], year-1)	 				
						is.in.phase3 <- end.phase2 < year-1
					}
					if(!is.in.phase3) {
						ulmfit <- use.lmfit[1]
						sigma.eps <- max(BHMp2[s,'sigma0'] + (tfr.pred[year-1,s] - BHMp2[s,'S_sd'])*
		  									ifelse(tfr.pred[year-1,s] > BHMp2[s,'S_sd'], 
		  											-BHMp2[s,'a_sd'], BHMp2[s,'b_sd']), wmeta$sigma0.min)
		  			} else { # phase 3
		  				ulmfit <- use.lmfit[2]
		  				sigma.eps <- this.sigma3
		  			}								
					tfr.pred[year, s] <- if(ulmfit) predict(scale.res$fit, data.frame(wtfr=wtrajs[year,s])) 
											else wtrajs[year,s] * scale + rnorm(1, 0, sd=sigma.eps)
				}
			}
			trajectories <- tfr.pred
			save(trajectories, file = file.path(outdir, paste0('traj_country', reg.obj$code, '.rda')))
			# compute quantiles
			PIs_cqp[region,,] <- apply(trajectories, 1, quantile, quantiles.to.keep, na.rm = TRUE)
 			mean_sd[region,1,] <- apply(trajectories, 1, mean, na.rm = TRUE)
 			mean_sd[region,2,] <- apply(trajectories, 1, sd, na.rm = TRUE) 	
  		}
  		present.year.index <- which(rownames(meta$tfr_matrix_all) == rownames(wmeta$tfr_matrix_all)[wpred$present.year.index])
		bayesTFR.prediction <- structure(list(
				quantiles = PIs_cqp,
				traj.mean.sd = mean_sd,
				nr.traj=nr.traj,
				tfr_matrix_reconstructed = meta$tfr_matrix_all,
				output.directory=outdir,
				na.index=wpred$na.index,
				mcmc.set=list(meta=meta, mcmc.list=list()),
				nr.projections=wpred$nr.projection,
				burnin=NA, thin=NA,
				end.year=wpred$end.year, use.tfr3=NA, burnin3=NA, thin3=NA,
				mu=NA, rho=NA, sigmaAR1 = NA, mu.c=NA, rho.c=NA,
				use.correlation=NA, start.year=wpred$start.year,
				present.year.index=present.year.index,
				present.year.index.all=present.year.index),
				class='bayesTFR.prediction')
			
		store.bayesTFR.prediction(bayesTFR.prediction, outdir)
		#do.convert.trajectories(pred=bayesTFR.prediction, n=save.as.ascii, output.dir=outdir, verbose=verbose)
		#if(write.summary.files)
		#	tfr.write.projection.summary.and.parameters(pred=bayesTFR.prediction, output.dir=outdir)
		result[[as.character(country.obj$code)]] <- bayesTFR.prediction
	}
	cat('\nPrediction stored into', output.dir, '\n')
	invisible(result)
}

tfr.predict.subnat.cor <- function(mcmc.set.list=NULL, sim.dir=file.path(getwd(), 'bayesTFR.output'), 
								countries = NULL, nr.traj = NULL, thin = NULL, burnin=2000,  use.diagnostics=FALSE,
								use.tfr3=TRUE, burnin3=2000,
								mu=2.1, rho=0.8859,
								use.phase3.from = NULL, use.correlation=FALSE, cor.method=c('bayes', 'median', 'mean'),
								cor.start.year=NULL,
								output.dir = NULL, ..., verbose=TRUE) {
	if(is.null(mcmc.set.list)) {
		cdirs <- list.files(file.path(sim.dir, 'subnat'), pattern='^c[0-9]+', full.names=TRUE)
		for(i in 1:length(cdirs)) {
			cdir.short <- strsplit(cdirs[i], split=.Platform$file.sep)
			cdir.short <- cdir.short[[1]][length(cdir.short[[1]])]
			country.code <- substring(cdir.short, 2)
			if(is.null(countries) || is.element(country.code, countries))
				mcmc.set.list[[country.code]] <- get.tfr.mcmc(cdirs[i])
		}
	}
	if(is.null(nr.traj) && is.null(thin) && !use.diagnostics)
		stop('Either nr.traj or thin must be given or use.diagnostics must be TRUE.')
	result <- cor.mat <- NULL
	outdir <- NULL
	for(country in names(mcmc.set.list)){
		if(use.correlation) {
			if(verbose) cat('\nCompute correlation for country', country, '...')
			snr.traj <- nr.traj
			sburnin <- burnin
			if(use.diagnostics) {
				nrtraj.burnin <- get.burnin.nrtraj.from.diagnostics(mcmc.set.list[[country]])
				sburnin <- nrtraj.burnin[2]
			}
			has.phase3 <- FALSE
			ar1pars <- NULL
			use.external.phase3 <- FALSE
			if(use.tfr3) {
				has.phase3 <- has.tfr3.mcmc(mcmc.set.list[[country]]$meta$output.dir)
				if(!has.phase3 && !is.null(use.phase3.from) && has.tfr3.mcmc(use.phase3.from)) {
					has.phase3 <- TRUE
					use.external.phase3 <- TRUE
				}
			}
			if(!has.phase3) {
				ar1pars <- if (is.null(rho) || is.na(rho)) get.ar1.parameters(mu = mu, mcmc.set.list[[country]]$meta)
							else list(mu=mu, rho=rho)
			}
			cor.mat <- tfr.correlation.subnat(mcmc.set.list[[country]], method=cor.method, burnin=sburnin, 
							burnin3=burnin3, use.external.phase3=use.external.phase3, use.phase3.from = use.phase3.from, 
							ar1pars=ar1pars, cor.start.year=cor.start.year, verbose=verbose)
			if(verbose) cat(' done.\n')
		}
		if(!is.null(output.dir)) {
			outdir <- file.path(output.dir, 'subnat', paste0('c', country))
		}		
		result[[country]] <- tfr.predict(mcmc.set=mcmc.set.list[[country]], nr.traj=nr.traj, thin=thin, burnin=burnin, 
										use.tfr3=use.tfr3, burnin3=burnin3, use.diagnostics=use.diagnostics,
										mu=mu, rho=rho, use.phase3.from = use.phase3.from,
										use.correlation=use.correlation, 
										correlation.matrices=cor.mat, output.dir=outdir, verbose=verbose, ...)
	}
	return(result)
}

tfr.predict.subnat1BHM <- function(countries, my.tfr.file, sim.dir=file.path(getwd(), 'bayesTFR.output'), 
								 nr.traj = NULL, #thin = NULL, burnin=2000,  use.diagnostics=FALSE,
								 start.year=NULL,
								burnin3=2000, mu=2.1, rho=0.8859,
								use.correlation=FALSE, cor.method=c('meth12'),
								cor.start.year=NULL,
								output.dir = NULL, ..., verbose=TRUE) {
	wpred <- get.tfr.prediction(sim.dir)
	wmcmc.set <- wpred$mcmc.set
	wmeta <- wmcmc.set$meta
	if (is.null(output.dir)) output.dir <- wmeta$output.dir
	stored.iter <- get.stored.mcmc.length(wmcmc.set$mcmc.list, burnin=0)		
	if(is.null(nr.traj)) nr.traj <- min(stored.iter, 2000)
	nr.traj <- min(nr.traj, stored.iter)
	do.not.thin <- nr.traj >= stored.iter-1
	if(!do.not.thin) {
		mcthin <- max(sapply(wmcmc.set$mcmc.list, function(x) x$thin))
		thin <- floor(stored.iter/nr.traj * mcthin)
		if(thin > 1) {
			thinned.mcmc <- get.thinned.tfr.mcmc(wmcmc.set, thin=thin, burnin=0)
			if(is.null(thinned.mcmc)) {
				wmcmc.set <- create.thinned.tfr.mcmc(wmcmc.set, thin=thin, burnin=0, 
				 							output.dir=output.dir, verbose=verbose)
				stored.iter <- get.stored.mcmc.length(wmcmc.set$mcmc.list, burnin=0)
			}
		}
		if(stored.iter != nr.traj) cat('\nNr. of trajectories changed to ', stored.iter, '\n')
	}
	has.p3 <- has.tfr3.mcmc(sim.dir)
	mcmc3 <- if(has.p3) get.tfr3.mcmc(sim.dir) else NULL
	result <- cor.mat <- NULL
	for(country in countries){
		country.obj <- get.country.object(country, wmeta)
		if(verbose) 
			cat('\nPredicting TFR for ', country.obj$name, '\n')
		if(is.null(start.year)) 
			start.year <- if(is.null(wpred$start.year)) as.integer(dimnames(wpred$tfr_matrix_reconstructed)[[1]][wpred$present.year.index])+5
						else wpred$start.year

		meta <- mcmc.meta.ini.subnat(wmeta, country=country.obj$code, my.tfr.file=my.tfr.file, 
									start.year=min(if(is.null(cor.start.year)) wmeta$start.year else cor.start.year, wmeta$start.year),
									present.year=wmeta$present.year, verbose=verbose)
		this.output.dir <- file.path(output.dir, 'subnat', paste0('c', country.obj$code))
		meta$output.dir <- this.output.dir		
		if(file.exists(this.output.dir)) unlink(this.output.dir, recursive=TRUE)
		dir.create(this.output.dir, recursive=TRUE)
		bayesTFR.mcmc.meta <- meta
		store.bayesTFR.meta.object(bayesTFR.mcmc.meta, this.output.dir)
		nr.reg <- get.nr.countries(meta)
		ar1pars <- NULL
		if(use.correlation) {
			cor.mat <- list()
			#for(w in c('low', 'high')) {
			#	cor.mat[[w]] <- matrix(0.99, nrow=nr.reg, ncol=nr.reg)
			#	diag(cor.mat[[w]]) <- 1
			#}
			if(!has.p3) {
				ar1pars <- if (is.null(rho) || is.na(rho)) get.ar1.parameters(mu = mu, wmeta)
							else list(mu=mu, rho=rho)
			}
			cor.mat <- tfr.correlation.subnat2(wmcmc.set, meta=meta, mcmc3.set=mcmc3, country=country, method=cor.method, burnin=0, 
							burnin3=burnin3, 
							ar1pars=ar1pars, cor.start.year=cor.start.year, verbose=verbose)
		}
		meta[['use.mcmc.from']] <- country.obj$code
		#meta$lambda_c[] <- wmeta$lambda_c[country.obj$index]
		#meta$id_Tistau <- if(country.obj$index %in% wmeta$id_Tistau) 1:nr.reg else c()
		data3.meta <- NULL	
		if(has.p3) {
			data3.meta <- mcmc3$meta
			# All sub-regions set as phase3 in order for all of them to draw from country-specific PIII parameters
			data3.meta$id_phase3 <- if(country.obj$index %in% mcmc3$meta$id_phase3) 1:nr.reg  else c()
			data3.meta$nr.countries <- length(mcmc3$meta$id_phase3)
		}
		result[[as.character(country)]] <- make.tfr.prediction(mcmc.set=wmcmc.set, start.year=start.year,
										use.correlation=use.correlation, use.tfr3=has.p3, mcmc3.set=mcmc3, burnin3=wpred$burnin3,
										correlation.matrices=cor.mat, output.dir=this.output.dir, 
										is.mcmc.set.thinned=TRUE, data.meta=meta, data3.meta=data3.meta,
										verbose=verbose, ...)
	}
	invisible(result)
}


make.tfr.prediction <- function(mcmc.set, start.year=NULL, end.year=2100, replace.output=FALSE,
								nr.traj = NULL, burnin=0, thin = NULL, 
								use.tfr3=TRUE, mcmc3.set=NULL, burnin3=0,
								mu=2.1, rho=0.9057, sigmaAR1 = 0.0922, 
								use.correlation=FALSE, correlation.matrices=NULL,
								countries = NULL,
								adj.factor1=NA, adj.factor2=0, forceAR1=FALSE,
								use.phase3.from = NULL,
								boost.first.period.in.phase2=TRUE,
							    save.as.ascii=1000, output.dir = NULL, write.summary.files=TRUE, 
							    is.mcmc.set.thinned=FALSE, force.creating.thinned.mcmc=FALSE,
							    write.trajectories=TRUE, data.meta=NULL, data3.meta=NULL,
							    verbose=verbose){
	# if 'countries' is given, it is an index
	# sigmaAR1 can be a vector. The last element will be repeated up to nr.projections
	meta <- if(is.null(data.meta)) mcmc.set$meta else data.meta
	#meta <- mcmc.set$meta
	#mc.meta <- if(!is.null(meta$use.mcmc.meta)) meta$use.mcmc.meta else meta
	mc.meta <- mcmc.set$meta
		
	present.year <- if(is.null(start.year)) meta$present.year else start.year - 5
	nr_project <- length(seq(present.year+5, end.year, by=5))
	suppl.T <- if(!is.null(meta$suppl.data$regions)) meta$suppl.data$T_end else 0
#	if (verbose)
		cat('\nPrediction from', present.year+5, 'until', end.year, '(i.e.', nr_project, 'projections)\n\n')
	l.sigmaAR1 <- length(sigmaAR1)
	sigma.end <- rep(sigmaAR1[l.sigmaAR1], nr_project + meta$T_end-l.sigmaAR1)
	sigmas_all <- c(sigmaAR1, sigma.end) 
	
	sigma0_s <- a_sd_s <- b_sd_s <- f_sd_s <- const_sd_s <- NULL
	burn <- if(is.mcmc.set.thinned) 0 else burnin
	total.iter <- get.total.iterations(mcmc.set$mcmc.list, burn)
	stored.iter <- get.stored.mcmc.length(mcmc.set$mcmc.list, burn)
	mcthin <- max(sapply(mcmc.set$mcmc.list, function(x) x$thin))
	if(!is.null(nr.traj) && !is.null(thin)) {
		warning('Both nr.traj and thin are given. Argument thin will be ignored.')
		thin <- NULL
	}
	if(is.null(nr.traj)) nr.traj <- min(stored.iter, 2000)
	else {
		if (nr.traj > stored.iter) 
			warning('nr.traj is larger than the available MCMC sample. Only ', stored.iter, ' trajectories will be generated.')
		nr.traj <- min(nr.traj, stored.iter)	
	}
	if(is.null(thin)) thin <- floor(stored.iter/nr.traj * mcthin)
	if(stored.iter <= 0 || thin == 0)
		stop('The number of simulations is 0. Burnin might be larger than the number of simulated values, or # trajectories is too big.')
	
	#setup output directory
	if (is.null(output.dir)) output.dir <- meta$output.dir
	outdir <- file.path(output.dir, 'predictions')
	
	if(is.null(countries)) {
		if(!replace.output && has.tfr.prediction(sim.dir=output.dir))
			stop('Prediction in ', outdir,
			' already exists.\nSet replace.output=TRUE if you want to overwrite existing projections.')
		unlink(outdir, recursive=TRUE)
		write.to.disk <- TRUE
		if(!file.exists(outdir)) 
			dir.create(outdir, recursive=TRUE)
	} else write.to.disk <- FALSE
	
	if(is.mcmc.set.thinned) { 
		thinned.mcmc <- mcmc.set
		has.thinned.mcmc <- TRUE
	} else {
		thinned.mcmc <- get.thinned.tfr.mcmc(mcmc.set, thin=thin, burnin=burnin)
		has.thinned.mcmc <- !is.null(thinned.mcmc) && thinned.mcmc$meta$parent.iter == total.iter
	}
	unblock.gtk('bDem.TFRpred')
	load.mcmc.set <- if(has.thinned.mcmc && !force.creating.thinned.mcmc) thinned.mcmc
					 else create.thinned.tfr.mcmc(mcmc.set, thin=thin, burnin=burnin, 
					 							output.dir=output.dir, verbose=verbose)
	nr_simu <- load.mcmc.set$mcmc.list[[1]]$finished.iter

	if (verbose) cat('Load variance parameters.\n')
	var.par.names <- c('sigma0', 'a_sd', 'b_sd', 'S_sd', 'const_sd')
	
	var.par.values <- get.tfr.parameter.traces(load.mcmc.set$mcmc.list, 
											var.par.names, burnin=0)
											
	prediction.countries <- if(is.null(countries)) 1:meta$nr_countries else countries
	nr_countries <- meta$nr_countries
	nr_countries_real <- length(prediction.countries)
	tfr_matrix_reconstructed <- get.tfr.reconstructed(meta$tfr_matrix_observed, meta)
	#ltfr_matrix <- dim(tfr_matrix_reconstructed)[1]
	#ltfr_matrix.all <- ltfr_matrix + suppl.T
	present.year.index <- get.estimation.year.index(meta, present.year)
	ltfr_matrix <- present.year.index
	ltfr_matrix.all <- present.year.index + suppl.T
	
	#keep these defaults for checking the out-of-sample projections
    quantiles.to.keep <- c(0,0.025,0.05,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,0.95,0.975,1)
	PIs_cqp <- array(NA, c(nr_countries, length(quantiles.to.keep), nr_project+1))
	dimnames(PIs_cqp)[[2]] <- quantiles.to.keep
	proj.middleyears <- get.prediction.years(meta, nr_project+1, present.year.index)
	dimnames(PIs_cqp)[[3]] <- proj.middleyears
	mean_sd <- array(NA, c(nr_countries, 2, nr_project+1))
	hasNAs <- rep(FALSE, nr_simu)
	adjust.true <- !is.na(adj.factor1)
	if (verbose) cat('Load parameters mean_eps_tau and sd_eps_tau.\n')
    tau.par.names <- c('mean_eps_tau', 'sd_eps_tau')
    tau.par.values <- get.tfr.parameter.traces(load.mcmc.set$mcmc.list, tau.par.names, burnin=0)

	if (verbose) cat('Load hierarchical parameters.\n')
	alpha.vars <- paste('alpha_',1:3, sep='')
	delta.vars <- paste('delta_',1:3, sep='')
	other.vars <- c('chi', 'psi', 'Triangle4', 'delta4')
	cs.par.values_hier <- newPointer(get.tfr.parameter.traces(load.mcmc.set$mcmc.list, 
										c(alpha.vars, delta.vars, other.vars), burnin=0))
										
	mid.years <- as.integer(c(if(suppl.T > 0) rownames(meta$suppl.data$tfr_matrix) else c(), rownames(tfr_matrix_reconstructed)))
	thin3 <- NA
	has.phase3 <- use.tfr3
	mu.c.mean <- rho.c.mean <- meta3 <- mc.meta3.pointer <- mcmc3.list.pointer <- NULL
	if(has.phase3) {
		mcmc3 <- NULL
		if(is.null(mcmc3.set)) {
			if(has.tfr3.mcmc(mc.meta$output.dir))
				mcmc3 <- get.tfr3.mcmc(mc.meta$output.dir) 
		} else mcmc3 <- mcmc3.set
		if((is.null(mcmc3) || length(mcmc3$meta$id_phase3)<=0) && !is.null(use.phase3.from)) {
			mcmc3 <- get.tfr3.mcmc(use.phase3.from) # use external phase3 parameters
			mcmc3$meta$id_phase3 <- c()             # only hierarchical pars. to be used
			if (has.tfr.prediction(sim.dir=use.phase3.from)) 
				burnin3 <- get.tfr.prediction(use.phase3.from)$burnin3
			if (verbose) cat('Using Phase III model (with burnin=', burnin3, ') from\n', use.phase3.from, '\n')
		}
					
		total.iter <- get.stored.mcmc.length(mcmc3$mcmc.list, burnin3)
		thinning.index <- unique(round(seq(1, total.iter, length=nr_simu)))
		if(length(thinning.index) < nr_simu) 
			stop('Length of MCMCs for phase 2 and 3 cannot be matched with these settings. Check arguments burnin, burnin3, nr.traj and thin.')
		m3.par.values <- get.tfr3.parameter.traces(mcmc3$mcmc.list, par.names=c('mu', 'rho', 'sigma.eps'), 
								thinning.index=thinning.index, burnin=burnin3)
		thin3 <- (thinning.index[2]-thinning.index[1])*mcmc3$mcmc.list[[1]]$thin
		if(dim(m3.par.values)[1] != nr_simu) stop('Mismatch in length of MCMCs for phase 2 and 3.')
		# used for the resulting prediction object
		mu <- mean(m3.par.values[,'mu'])
		rho <- mean(m3.par.values[,'rho'])
		sigmaAR1 <- mean(m3.par.values[,'sigma.eps'])
		mu.c.mean <- rho.c.mean <- rep(NA, nr_countries)
		meta3 <- if(is.null(data3.meta)) mcmc3$meta else data3.meta
		mc.meta3.pointer <- newPointer(mcmc3$meta)
		mcmc3.list.pointer <- newPointer(mcmc3$mcmc.list)
	}
	max.nr.project <- nr_project
	all.T_end.min <- ltfr_matrix.all
	# load country-specific parameters
	if (verbose) cat('Load country-specific parameters.\n')
	cs.par.values.list <- cs.var.names <- theta_si.list <- country.objects <- all.tfr.list <- m3.par.values.cs.list <- list()
	nmissing <- list()
	mc.meta.pointer <- newPointer(mc.meta)
	meta.pointer <- newPointer(meta)
	mcmc.list.pointer <- newPointer(mcmc.set$mcmc.list)
	load.mcmc.list.pointer <- newPointer(load.mcmc.set$mcmc.list)
	load.meta.pointer <- newPointer(load.mcmc.set$meta)
	
	# country loop for preparing data for projections
	if(!is.null(meta$use.mcmc.from)) {
		country.obj <- get.country.object(meta$use.mcmc.from, mc.meta)
		dprep <- .prepare.country.spec.pars.for.predictions(country.obj$index, country.obj, mc.meta.pointer, mcmc.list.pointer, 
														mc.meta.pointer, load.mcmc.list.pointer, 
														nr_simu, burnin, alpha.vars, delta.vars, has.phase3, 
														mc.meta3.pointer, mcmc3.list.pointer, burnin3, thinning.index, cs.par.values_hier)
	}
	for (country in prediction.countries){
		country.obj <- get.country.object(country, meta, index=TRUE)
		if(is.null(meta$use.mcmc.from)) {
			dprep <- .prepare.country.spec.pars.for.predictions(country, country.obj, meta.pointer, mcmc.list.pointer, 
														load.meta.pointer, load.mcmc.list.pointer, nr_simu, burnin,
														alpha.vars, delta.vars, has.phase3, mc.meta3.pointer, mcmc3.list.pointer, burnin3, thinning.index,
														cs.par.values_hier)
		}
		cs.var.names[[country]] <- list(U=dprep$U.var, Triangle_c4=dprep$Triangle_c4.var)
		cs.par.values.list[[country]] <- dprep$cs.par.values[,c(dprep$Triangle_c4.var, dprep$U.var)]
		theta_si.list[[country]] <- dprep$theta_si
		country.objects[[country]] <- country.obj
		# get the whole TFR time series including supplemental historical data
		all.tfr <- get.observed.tfr(country, meta, 'tfr_matrix_observed', 'tfr_matrix_all')
		all.tfr.list[[country]] <- all.tfr
		this.T_end.min <- min(meta$T_end_c[country], ltfr_matrix.all)
		all.T_end.min <- min(all.T_end.min, this.T_end.min)
		nmissing[[country]] <- ltfr_matrix.all - this.T_end.min # positive only if this.T_end < ltfr_matrix.all
		max.nr.project <- max(max.nr.project, nr_project + nmissing[[country]])
		
		# load phase3 country-specific parameter traces
		if(has.phase3 && is.element(country, meta3$id_phase3)) {
			m3.par.values.cs.list[[country]] <- dprep$m3.par.values.cs
			mu.c.mean[country] <- mean(m3.par.values.cs.list[[country]][,1])
			rho.c.mean[country] <- mean(m3.par.values.cs.list[[country]][,2])
		}
	} # end country prep loop
	
	if(use.correlation) {
		# prepare AR1 eps for joint sampling
		eps.correlation <- correlation.matrices
		if(is.null(eps.correlation))	
			eps.correlation <- tfr.correlation(meta) 
		cor.mat.na <- which(apply(is.na(eps.correlation$low), 2, sum) > dim(eps.correlation$low)[1]-2)
		nr.countries.no.na <- nr_countries - length(cor.mat.na)
		epsilons <- rep(NA, nr_countries)
		kappa<-eps.correlation$kappa
	}
	# array for results - includes also historical data for periods with missing data
	all.f_ps <- array(NA, dim=c(nr_countries_real, max.nr.project+1, nr_simu))
	# vector W with the weight for the first two periods:
	W <- matrix(0, nrow=nr_countries_real, ncol=max.nr.project+1)
	which.Wsecond <- rep(NA, nr_countries_real)
	# index of the last period within all.f_ps that is observed
	fps.end.obs.index <- dim(tfr_matrix_reconstructed)[1] - all.T_end.min + suppl.T + 1
	first.projection <- rep(1, nr_countries_real)
	for (icountry in 1:nr_countries_real) {
		# fill the result array with observed data 
		for(year in 1:fps.end.obs.index) 
			all.f_ps[icountry,year,] <- all.tfr.list[[prediction.countries[icountry]]][all.T_end.min+year-1]
		first.two.na <- which(is.na(all.f_ps[icountry,,1]))[1:2]
		which.Wsecond[icountry] <- first.two.na[2]
		W[icountry,first.two.na] <- c(adj.factor1, adj.factor2)
		first.projection[icountry] <- first.two.na[1]
	}
	W[is.na(W)] <- 0
	mu.c <- rho.c <- rep(NA, nr_countries)
	sigma.epsAR1 <- list()
	if(length(sigmas_all) < max.nr.project) {
		sigmas_all <- c(rep(sigmas_all[1], max.nr.project-length(sigmas_all)), sigmas_all)
	}
	if(!has.phase3) {
		mu.c <- rep(mu, nr_countries)
		rho.c <- rep(rho, nr_countries)
		sigma.epsAR1 <- rep(list(sigmas_all), nr_countries)
	}
	traj.counter <- 0
	country.loop.max <- 20
	status.for.gui <- paste('out of', nr_simu, 'trajectories.')
	gui.options <- list()
	if (verbose) {
		verbose.iter <- as.integer(max(1, nr_simu/100))
		if(interactive()) cat('\n')
	}
	#########################################
	for (s in 1:nr_simu){ # Iterate over trajectories
	#########################################
		if(getOption('bDem.TFRpred', default=FALSE)) {
			# This is to unblock the GUI, if the run is invoked from bayesDem
			# and pass info about its status
			# In such a case the gtk libraries are already loaded
			traj.counter <- traj.counter + 1
			gui.options$bDem.TFRpred.status <- paste('finished', traj.counter, status.for.gui)
			unblock.gtk('bDem.TFRpred', gui.options)
		}
		if (verbose) {
			if(interactive()) cat('\rProjecting TFR trajectories ... ', round(s/nr_simu * 100), ' %')
			else {
				if (s %% verbose.iter == 0) 
					cat('TFR projection trajectory ', s, '\n')
				}
		}
		if(has.phase3) { # set country-spec parameters for phase 3 - time-invariant
			mu.c[] <- m3.par.values[s,'mu']
			rho.c[] <- m3.par.values[s,'rho']
			sigma.epsAR1 <- rep(list(rep(m3.par.values[s,'sigma.eps'], max.nr.project)), nr_countries)
			for (country in prediction.countries[prediction.countries %in% meta3$id_phase3]){	
					mu.c[country] <- m3.par.values.cs.list[[country]][s,1]
					rho.c[country] <- m3.par.values.cs.list[[country]][s,2]
			}
		}
		is.in.phase3 <- rep(forceAR1, nr_countries_real)
		S11 <- rep(0, nr_countries_real)
		#########################################
		for (year in 2:(max.nr.project+1)) { # Iterate over time
		#########################################
			ALLtfr.prev <- all.f_ps[,year-1,s]
			if(use.correlation) {
				cor.mat <- eps.correlation$low
				hiTFR <- which(ALLtfr.prev >= kappa)
				# this leaves low on spots where both countries are low
				cor.mat[hiTFR,] <- eps.correlation$high[hiTFR,]
        		cor.mat[,hiTFR] <- eps.correlation$high[,hiTFR]		
        		cor.mat.no.na <- if(length(cor.mat.na)==0) cor.mat else cor.mat[-cor.mat.na, -cor.mat.na]
        		if(!is.cor.positive.definite(cor.mat.no.na))
        			cor.mat.no.na <- zero.neg.evals(cor.mat.no.na)
        	} 
        	stop.country.loop <- FALSE
			country.loop <- 1  
			# loop for resampling if tfr is outside of the bounds
			# if use.correlation is FALSE, it goes through only once
        	while(!stop.country.loop && (country.loop<=country.loop.max)) { 
				country.loop <- country.loop + 1
				stop.country.loop <- TRUE
				if(use.correlation) {
        			epsilons.no.na <- mvrnorm(1,rep(0,nr.countries.no.na), cor.mat.no.na)
        			if(length(cor.mat.na)>0) epsilons[-cor.mat.na] <- epsilons.no.na
        			else epsilons[] <- epsilons.no.na
				}
				tfr.c <- all.f_ps[, year,s]
				#stop('')
				#########################################
				for (icountry in 1:nr_countries_real){ # Iterate over countries
				#########################################
					if(!is.na(all.f_ps[icountry, year,s])) next
					country <- prediction.countries[icountry]	# index within meta (icountry is index within countries for which this is run)
					this.T_end <- meta$T_end_c[country]
					all.tfr <- all.tfr.list[[country]]
		 			if(!is.in.phase3[icountry]) {
		 				# check if now in phase 3
						if(year == first.projection[icountry]) { # first projection period
							if(!is.element(country, meta$id_Tistau)) {
								is.in.phase3[icountry] <- ((meta$lambda_c[country] < this.T_end) || 
	                							((min(all.tfr[1:this.T_end], na.rm=TRUE) <= 
	                									cs.par.values.list[[country]][s, cs.var.names[[country]]$Triangle_c4]) && 
	                 								(all.tfr[this.T_end] > all.tfr[this.T_end-1])))
	                 		}
	                 	} else is.in.phase3[icountry] <- ((min(all.f_ps[icountry, 1:(year-1),s]) <= 
	                 										cs.par.values.list[[country]][s, cs.var.names[[country]]$Triangle_c4]) && 
	                 									(all.f_ps[icountry, year-1,s] > all.f_ps[icountry,year-2,s]))
	                }
					if(adjust.true) {
						if(year == first.projection[icountry]) { # first projection period
							D11 <- (all.tfr[this.T_end-1] - all.tfr[this.T_end])
				 			if(!is.in.phase3[icountry]) { # country in Phase II				
		           				d11 <- DLcurve(theta_si.list[[country]][s,], all.tfr[this.T_end-1], meta$dl.p1, meta$dl.p2)
			 					S11[icountry] <- D11 - d11
			  				} else { # country in Phase III	
								S11[icountry] <- D11 - (all.tfr[this.T_end-1] - 
													(mu.c[country] + rho.c[country]*(all.tfr[this.T_end-1]-mu.c[country])))
							}
						}
						if(is.in.phase3[icountry] && year == which.Wsecond[icountry]) {
							# if a country with adjustment enters Phase III in second proj. step (in normal case corresponds to year ==3)
							# then the adjustment needs to be changed, based on observed diff in last proj step and AR(1) decrement
			  				S11[icountry] <- ( all.f_ps[icountry,year-2,s] - all.f_ps[icountry, year-1,s]) - ( 
			  							all.f_ps[icountry, year-2,s] - (mu.c[country] + rho.c[country]*( 
			  									all.f_ps[icountry,year-2,s]-mu.c[country])))
				 		}
		  			}
		  			#if(s==5) {
		  			#	cobj <- get.country.object(country, meta, index=TRUE)
		  			#	print(c(year, cobj$name, c('II', 'III')[is.in.phase3[icountry]+1]))
		  			#}
		  			# Simulate projection
					if (!is.in.phase3[icountry]){ # Phase II						
						old.tfr <- all.f_ps[icountry,year-1,s]
						# for subnational TFR get decrements based on U_c as maximum (for TFR larger than U_c the decrements would be close to zero)
						if(!is.null(meta$use.mcmc.from) && all.f_ps[icountry,year-1,s] > sum(theta_si.list[[country]][s,1:4]))
							old.tfr <- sum(theta_si.list[[country]][s,1:4])
						new.tfr <- (all.f_ps[icountry,year-1,s]- DLcurve(theta_si.list[[country]][s,], old.tfr, meta$dl.p1, meta$dl.p2) - 
										W[icountry,year]*S11[icountry])
						# get errors
						if(boost.first.period.in.phase2 && is.element(country, meta$id_Tistau) && (year == first.projection[icountry])) {
							eps.mean <- tau.par.values[s, 'mean_eps_tau']
							sigma_eps <- tau.par.values[s, 'sd_eps_tau']
							if(use.correlation && !is.na(epsilons[country])) sigma_eps <- rnorm(1, eps.mean, sigma_eps)
						} else {
							eps.mean <- 0
							sigma_eps <- max(var.par.values[s,'sigma0'] + (all.f_ps[icountry, year -1,s] - var.par.values[s,'S_sd'])*
		  									ifelse(all.f_ps[icountry, year -1,s] > var.par.values[s,'S_sd'], 
		  											-var.par.values[s,'a_sd'], var.par.values[s,'b_sd']), meta$sigma0.min)
						}
						if(!use.correlation || is.na(epsilons[country])) {
							passed <- FALSE
		       				for(i in 1:50) {
		       					err <- rnorm(1, eps.mean, sigma_eps)
		                    	if( (new.tfr + err) > 0.5 && 
		                    		(new.tfr + err) <= cs.par.values.list[[country]][s,cs.var.names[[country]]$U] ) {passed <- TRUE; break}
		                	}
		                	if(!passed) err <- min(max(err, 0.5-new.tfr), cs.par.values.list[[country]][s,cs.var.names[[country]]$U]-new.tfr)
		                } else { # joint predictions
		                	err <- sigma_eps*epsilons[country]
		                	if(err < 0.5 - new.tfr || err > cs.par.values.list[[country]][s,cs.var.names[[country]]$U]-new.tfr) {# TFR outside of bounds
		                		stop.country.loop <- FALSE
		                		if(country.loop < country.loop.max) break
		                		else err <- min(max(err, 0.5-new.tfr), cs.par.values.list[[country]][s,cs.var.names[[country]]$U]-new.tfr)
		                	}
		                }
					} else { # Phase III
						new.tfr <- (mu.c[country] + rho.c[country]*(all.f_ps[icountry,year-1,s] - mu.c[country]) 
										- W[icountry,year]*S11[icountry])						
						if(!use.correlation || is.na(epsilons[country])) {
							passed <- FALSE
	 						for(i in 1:50){
	 							err <- rnorm(1, 0, sigma.epsAR1[[country]][year-1])
	 							if (new.tfr + err > 0.5 )   {passed <- TRUE; break}
							}
							if(!passed) err <- 0.5 - new.tfr
						} else { # joint predictions
							err <- sigma.epsAR1[[country]][year-1]*epsilons[country]
							if(err < 0.5 - new.tfr) {
		                		stop.country.loop <- FALSE
		                		if(country.loop < country.loop.max) break
		                		else err <- 0.5 - new.tfr
							}
						}		
					}
					tfr.c[icountry] <- new.tfr + err
				} # end countries loop
			} # end while loop
			all.f_ps[,year,s] <- tfr.c
		} # end time loop
	} # end simu loop
	if(verbose && interactive()) cat('\n')
	##############
	# Impute missing values if any and compute quantiles
	for (icountry in 1:nr_countries_real){
		country <- prediction.countries[icountry]
		# Ignore trajectories that go below 0.5 
		isnotNA <- apply(1-(all.f_ps[icountry,,] < 0.5), 2, prod) 
		isnotNA <- ifelse(is.na(isnotNA),0,isnotNA)
		all.f_ps[icountry,,isnotNA==0] <- NA
		# extract the future trajectories (including the present period)
		f_ps_future <- all.f_ps[icountry,(dim(all.f_ps)[2]-nr_project):dim(all.f_ps)[2],]
		if (nmissing[[country]] > 0) { # data imputation
			f_ps_future[1,] <- quantile(f_ps_future[1,], 0.5, na.rm = TRUE) # set all trajectories in the first time period to the median
			tfr_matrix_reconstructed[(ltfr_matrix-fps.end.obs.index+2):ltfr_matrix,country] <- apply(
											all.f_ps[icountry,2:fps.end.obs.index,,drop=FALSE],
												2, quantile, 0.5, na.rm = TRUE)
			if (verbose) 
				cat('\t', nmissing[[country]], 'data points reconstructed for', country.objects[[country]]$name,'\n')
		}
		this.hasNAs <- apply(is.na(f_ps_future), 2, any)
		hasNAs[this.hasNAs] <- TRUE
		if(write.trajectories) {
			trajectories <- f_ps_future # save only trajectories simulated for the future time
  			save(trajectories, file = file.path(outdir, paste('traj_country', country.objects[[country]]$code, '.rda', sep='')))
  		}
  		# compute quantiles
 		PIs_cqp[country,,] <- apply(f_ps_future, 1, quantile, quantiles.to.keep, na.rm = TRUE)
 		mean_sd[country,1,] <- apply(f_ps_future, 1, mean, na.rm = TRUE)
 		mean_sd[country,2,] <- apply(f_ps_future, 1, sd, na.rm = TRUE) 		
 	}
	load.mcmc.set <- remove.tfr.traces(load.mcmc.set)
	if(!is.null(data.meta)) load.mcmc.set$meta <- data.meta
	bayesTFR.prediction <- structure(list(
				quantiles = PIs_cqp,
				traj.mean.sd = mean_sd,
				nr.traj=nr_simu,
				tfr_matrix_reconstructed = tfr_matrix_reconstructed,
				output.directory=outdir,
				na.index=(1:nr_simu)[hasNAs],
				mcmc.set=load.mcmc.set,
				nr.projections=nr_project,
				burnin=burnin, thin=thin,
				end.year=end.year, use.tfr3=has.phase3, burnin3=burnin3, thin3=thin3,
				mu=mu, rho=rho, sigmaAR1 = sigmaAR1, mu.c=mu.c.mean, rho.c=rho.c.mean,
				use.correlation=use.correlation, start.year=start.year,
				present.year.index=present.year.index,
				present.year.index.all=ltfr_matrix.all),
				class='bayesTFR.prediction')
			
	if(write.to.disk) {
		store.bayesTFR.prediction(bayesTFR.prediction, outdir)
		do.convert.trajectories(pred=bayesTFR.prediction, n=save.as.ascii, output.dir=outdir, verbose=verbose)
		if(write.summary.files)
			tfr.write.projection.summary.and.parameters(pred=bayesTFR.prediction, output.dir=outdir)
		cat('\nPrediction stored into', outdir, '\n')
	}
	invisible(bayesTFR.prediction)
}
newPointer <- function(inputValue){
	object <- new.env(parent=globalenv())
	object$value <- inputValue
	class(object) <- 'pointer'
	return(object)
}

getValue <- function(pointer)
	return(pointer$value)


.prepare.country.spec.pars.for.predictions <- function(country, country.obj, meta, mcmc.list, load.meta, load.mcmc.list, nr_simu, burnin,
														alpha.vars, delta.vars, has.phase3, meta3, mcmc3.list, burnin3, thinning.index,
														cs.par.values_hier) {
	if (is.element(country,getValue(meta)$id_DL)){
		U.var <- paste0('U_c', country.obj$code)
		d.var <- paste0('d_c', country.obj$code)
		Triangle_c4.var <- paste0("Triangle_c4_c", country.obj$code)
		gamma.vars <- paste0('gamma_',1:3,'_c', country.obj$code)
		if(country <= getValue(meta)$nr_countries_estimation) {
			cs.par.values <- get.tfr.parameter.traces.cs(getValue(load.mcmc.list), country.obj, 
									tfr.parameter.names.cs(trans=FALSE), burnin=0)
		} else {
			# if it's country that was not included in the estimation, determine the posterior index
			# again, since the MCMC might be of different length
			if (is.element(country.obj$code, getValue(load.meta)$regions$country_code)) {
				cs.par.values <- get.tfr.parameter.traces.cs(getValue(load.mcmc.list), country.obj, 
									tfr.parameter.names.cs(trans=FALSE), burnin=0)
			} else { # there are no thinned traces for this country, use the full traces 
				cs.par.values <- get.tfr.parameter.traces.cs(getValue(mcmc.list), country.obj, 
									tfr.parameter.names.cs(trans=FALSE), burnin=burnin)
				selected.simu <- get.thinning.index(nr_simu, dim(cs.par.values)[1])
				if (length(selected.simu$index) < nr_simu)
					selected.simu$index <- sample(selected.simu$index, nr_simu, replace=TRUE)
				cs.par.values <- cs.par.values[selected.simu$index,]
			}
		}
		pc_si <- matrix(NA, nr_simu, 3)
		for (i in 1:3)
	 		pc_si[,i] <- exp(cs.par.values[,gamma.vars[i]])/
	 								apply(exp(cs.par.values[,gamma.vars]), 1,sum)
		theta_si <- cbind(
					(cs.par.values[,U.var] - cs.par.values[, Triangle_c4.var])*pc_si[,1],
					(cs.par.values[,U.var] - cs.par.values[, Triangle_c4.var])*pc_si[,2],
					(cs.par.values[,U.var] - cs.par.values[, Triangle_c4.var])*pc_si[,3],
					cs.par.values[, Triangle_c4.var], 
					cs.par.values[,d.var])
	} else { #Tistau countries
    	# sample decline parameters from the hier distributions                 
		Triangle4_tr_s <- rnorm(nr_simu, mean = getValue(cs.par.values_hier)[,'Triangle4'], 
										sd = getValue(cs.par.values_hier)[, 'delta4'])
		Triangle_c4_s <- (getValue(meta)$Triangle_c4.up*exp(Triangle4_tr_s) + getValue(meta)$Triangle_c4.low)/(1+exp(Triangle4_tr_s))
	
		# need U and Triangle_c4 in cs... later in loop for start of phase III and prior on f_t
		cs.par.values = rep(get.observed.tfr(country, getValue(meta), 'tfr_matrix_all')[getValue(meta)$tau_c[country]], nr_simu)
		Triangle_c4.var <- 'Triangle_c4'
		U.var <- 'U'
		cs.par.values <- cbind(cs.par.values, Triangle_c4_s)
		colnames(cs.par.values) <- c(U.var, Triangle_c4.var)
		d_tr_s <- rnorm(nr_simu, mean = getValue(cs.par.values_hier)[,'chi'], sd = getValue(cs.par.values_hier)[,'psi'])
		d_s <-  (getValue(meta)$d.up*(exp(d_tr_s) + getValue(meta)$d.low)/(1+exp(d_tr_s)))	
		gamma_si <- matrix(NA, nr_simu, 3)
		for (i in 1:3)
	 		gamma_si[,i] <- rnorm(nr_simu, mean = getValue(cs.par.values_hier)[,alpha.vars[i]], 
	 											sd = getValue(cs.par.values_hier)[,delta.vars[i]])
		pc_si <- matrix(NA, nr_simu, 3)
		for (i in 1:3) pc_si[,i] <- exp(gamma_si[,i])/apply(exp(gamma_si), 1,sum)
		theta_si <- cbind((cs.par.values[,U.var] - Triangle_c4_s)*pc_si[,1],
	                      (cs.par.values[,U.var] - Triangle_c4_s)*pc_si[,2],
	                      (cs.par.values[,U.var] - Triangle_c4_s)*pc_si[,3],
	                          Triangle_c4_s, d_s) 
	}
	m3.par.values.cs <- NULL
	if(has.phase3 && !is.null(getValue(meta3)) && is.element(country, getValue(meta3)$id_phase3))
			m3.par.values.cs <- get.tfr3.parameter.traces.cs(getValue(mcmc3.list), country.obj=country.obj,
											par.names=c('mu.c', 'rho.c'), burnin=burnin3, thinning.index=thinning.index)
	return(list(U.var=U.var, Triangle_c4.var=Triangle_c4.var, theta_si=theta_si, cs.par.values=cs.par.values,
			m3.par.values.cs=m3.par.values.cs))
}


is.cor.positive.definite <- function(m, tol = 1e-06) {
	E <- eigen(m, symmetric=TRUE)
	ev <- E$values
    return(all(ev >= -tol * abs(ev[1L])))
}

zero.neg.evals <- function(eps.cor)
{
	zero <- 1e-10 
	E <- eigen(eps.cor,symmetric=TRUE)

    # compute eigenvectors/-values
	V   <- E$vectors
	D   <- E$values

	# replace negative eigenvalues by zero
	D   <- pmax(D,zero)

	# reconstruct correlation matrix
	new.cor.mat  <- V %*% diag(D) %*% t(V)

	# rescale correlation matrix
	T   <- 1/sqrt(diag(new.cor.mat))
	TT  <- outer(T,T)
	return(new.cor.mat * TT)
}


get.ar1.countries <- function(meta) {
	index <- get.ar1.countries.index(meta)
	return(data.frame(country_name=meta$regions$country_name[index], country_code=meta$regions$country_code[index]))
}

get.ar1.countries.index <- function(meta) {
	nr.countries <- get.nr.countries.est(meta)
	return(seq(1, nr.countries)[meta$lambda_c[1:nr.countries]!=meta$T_end_c[1:nr.countries]])
}

get.ar1.data <- function(meta) {
	tfr_prev <- tfr_now <- NULL
    for (country in get.ar1.countries.index(meta)){  
		tfr <- get.observed.tfr(country, meta, 'tfr_matrix_all')
		tfr_prev <- c(tfr_prev, tfr[meta$lambda_c[country]:(meta$T_end_c[country]-1)])
		tfr_now <- c(tfr_now, tfr[(meta$lambda_c[country]+1):meta$T_end_c[country]] )
	}
	return(list(tfr_prev=tfr_prev, tfr_now=tfr_now, countries=get.ar1.countries(meta)))
}

get.ar1.parameters <- function(mu = 2.1, meta){
	ar1data <- get.ar1.data(meta)
	yt <- ar1data$tfr_now - mu
	ytm1 <- ar1data$tfr_prev - mu
	mod = lm(yt~-1 +ytm1)
	rho = mod$coeff[1]
	sigmaAR1 = sqrt(sum(mod$residuals^2)/(length(ar1data$tfr_now)-1))
	#tfr = ifelse(meta$tfr_matrix_all[,1:nr.countries]<=mu,meta$tfr_matrix_all[,1:nr.countries], NA)
	#sd_tot = sd(c(tfr, 2*mu-tfr), na.rm = TRUE)
	#sigma.end = sd_tot*sqrt(1-rho^2)
	return( #list(rho = round(rho,2), sigmaAR1 = round(sigmaAR1,2))
			list(rho = rho, sigmaAR1 = sigmaAR1, mu=mu, data=ar1data)
				)
}

remove.tfr.traces <- function(mcmc.set) {
	for (i in 1:length(mcmc.set$mcmc.list)) {
		mcmc.set$mcmc.list[[i]]$traces <- 0
		mcmc.set$mcmc.list[[i]]$burnin <- 0
	}
	invisible(mcmc.set)
}

"get.traj.ascii.header" <- function(meta, ...) UseMethod("get.traj.ascii.header")
get.traj.ascii.header.bayesTFR.mcmc.meta <- function(meta, ...) 
	return (list(country_code='LocID', period='Period', year='Year', trajectory='Trajectory', tfr='TF'))
		
store.traj.ascii <- function(trajectories, n, output.dir, country.code, meta, index, append=FALSE, present.index=NULL) {
	# Store trajectories into ASCII files of a specific UN format 
	#header <- list(country_code='LocID', period='Period', year='Year', trajectory='Trajectory', tfr='TF')
	header <- get.traj.ascii.header(meta)
	nyears <- dim(trajectories)[1]
	pred.years <- get.prediction.years(meta, nyears, present.index)
	pred.period <- get.prediction.periods(meta, nyears, present.index)
	results <- NULL
	for (traj in 1:length(index)) {
		results <- rbind(results, cbind(country_code=rep(country.code, nyears), 
								period=pred.period, year=pred.years, 
								trajectory=rep(index[traj], nyears), 
								tfr=round(trajectories[,index[traj]], 5)))
	}
	#match column names and header
	colnames(results)[colnames(results)==names(header)] <- header
	write.table(results, file=file.path(output.dir, 'ascii_trajectories.csv'), sep=',', 
					quote=FALSE, row.names=FALSE, col.names=!append, append=append)
	return(results)
}

get.predORest.year.index <- function(pred, year) {
	projection.index <- get.prediction.year.index(pred, year)
	projection <- TRUE
	if(is.null(projection.index)) {
		projection <- FALSE
		projection.index <- get.estimation.year.index(pred$mcmc.set$meta, year)
	}
	return(c(index=projection.index, is.projection=projection))
}

get.prediction.year.index <- function(pred, year) {
	years <- get.all.prediction.years(pred)
	lyears <- length(years)
	breaks <- c(years-3, years[lyears]+2)
	h <- try(hist(year, breaks=breaks, plot=FALSE)$count, silent=TRUE)
	return(if(inherits(h, "try-error")) NULL else which(h > 0)[1])
}

get.all.prediction.years <- function(pred) {
	return(get.prediction.years(pred$mcmc.set$meta, pred$nr.projections+1, pred$present.year.index))
}

get.prediction.years <- function(meta, n, present.year.index=NULL) {
	if(is.null(present.year.index)) present.year.index <- nrow(get.data.matrix(meta))
	present.year <-  as.numeric(rownames(get.data.matrix(meta))[present.year.index])
	return (seq(present.year, length=n, by=5))
}

get.prediction.periods <- function(meta, n, ...) {
	mid.years <- get.prediction.years(meta, n, ...)
	return (paste(mid.years-3, mid.years+2, sep='-'))
}

get.estimation.years <- function(meta)
	return(as.numeric(rownames(get.data.matrix(meta))))
	
get.estimation.year.index <- function(meta, year) {
	years <- get.estimation.years(meta)
	lyears <- length(years)
	breaks <- c(years-3, years[lyears]+2)
	h <- try(hist(year, breaks=breaks, plot=FALSE)$count, silent=TRUE)
	return(if(inherits(h, "try-error")) NULL else which(h > 0)[1])
}

get.tfr.periods <- function(meta) {
	mid.years <- get.estimation.years(meta)
	return (paste(mid.years-3, mid.years+2, sep='-'))
}

do.convert.trajectories <- function(pred, n, output.dir, countries=NULL, verbose=FALSE) {
	# Converts all trajectory rda files into UN ascii, selecting n trajectories by equal spacing.
	if(n==0) return(NULL)
	nr.simu <- pred$nr.traj
	has.na <- rep(FALSE, nr.simu)
	has.na[pred$na.index] <- TRUE
	if (n=='all') n <- nr.simu
	selected.traj <- get.thinning.index(n, nr.simu)
	is.selected <- rep(FALSE, nr.simu)
	is.selected[selected.traj$index] <- TRUE
	is.sel.has.na <- is.selected & has.na
	for(NAidx in (1:nr.simu)[is.sel.has.na]) { #for selected NA-spots, find the closest neighbours that are not NA
		is.selected[NAidx] <- FALSE
		if(n==nr.simu) next
		i <- NAidx-1
		dist <- 1
		while(TRUE) {
			if(i>0) {
				if (!is.selected[i] & !has.na[i]) { # looking at lower index
					is.selected[i] <- TRUE
					break
					}
				}
			i <- NAidx + dist
			if (i > nr.simu) break 
			if (is.selected[i]) break
			if (!has.na[i]) { # looking at higher index
				is.selected[i] <- TRUE
				break
			}
			dist <- dist + 1
			i <- NAidx - dist
		}
	}
	index <- (1:nr.simu)[is.selected]
	country.codes <- country.names <- c()
	result.wide <- c()
	header <- get.traj.ascii.header(pred$mcmc.set$meta)
	convert.countries <- if(is.null(countries)) pred$mcmc.set$meta$regions$country_code else countries
	for (country in convert.countries) {
		country.obj <- get.country.object(country, pred$mcmc.set$meta)
		if(verbose) cat('Converting trajectories for', country.obj$name, '(code', country.obj$code, ')\n')
		trajectories <- get.trajectories(pred, country.obj$code)$trajectories
		if (is.null(trajectories)) {
			warning('No trajectories for ', country.obj$name, ' (code ', country.obj$code, ')')
		} else {
			append <- length(country.codes) > 0
			country.codes <- c(country.codes, country.obj$code)
			country.names <- c(country.names, country.obj$name)			
			result <- store.traj.ascii(trajectories, n, output.dir, country.obj$code, 
							pred$mcmc.set$meta, index=index, append=append, present.index=pred$present.year.index)
			if(!append) {
				result.wide <- result[,2:5]
			} else {
				result.wide <- cbind(result.wide, result[,header$tfr])
			}
		}
	}
	# order result.wide by country name
	o <- order(country.names)
	result.wide[,4:ncol(result.wide)] <- result.wide[,3+o]
	# write transposed version
	file.wide <- file.path(output.dir, 'ascii_trajectories_wide.csv')
	colnames(result.wide) <- c('Period', 'Year', 'Trajectory', country.names[o])
	write.table(rbind(c(' ', ' ', 'LocID', country.codes[o]), colnames(result.wide)), 
					file=file.wide, sep=',', 
					quote=TRUE, row.names=FALSE, col.names=FALSE)
	write.table(result.wide, file=file.wide, sep=',', 
					quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

	if(verbose) cat('Number of trajectories stored for each country:', length(index), '\n')
	cat('Converted trajectories stored into', file.path(output.dir, 'ascii_trajectories(_wide).csv'), '\n')
}


convert.tfr.trajectories <- function(dir=file.path(getwd(), 'bayesTFR.output'), 
								 n=1000, output.dir=NULL, 
								 verbose=FALSE) {
	# Converts all trajectory rda files into UN ascii, selecting n trajectories by equal spacing.
	if(n <= 0) return()
	pred <- get.tfr.prediction(sim.dir=dir)
	if (is.null(output.dir)) output.dir <- pred$output.directory
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	cat('Converting trajectories from', dir, '\n')
	if (is.null(pred$na.index)) {
		if(verbose) cat('Finding NA values in each country ...\n')
		for (country in 1:pred$mcmc.set$meta$nr_countries) {
			country.obj <- get.country.object(country, pred$mcmc.set$meta, index=TRUE)
			trajectories <- get.trajectories(pred, country.obj$code)$trajectories
			if (country==1) hasNAs <- rep(FALSE, dim(trajectories)[2])
			this.hasNAs <- apply(is.na(trajectories), 2, any)
			hasNAs[this.hasNAs] <- TRUE
		}
		pred$na.index <- (1:pred$nr.traj)[hasNAs]
	}
	do.convert.trajectories(pred=pred, n=n, output.dir=output.dir, verbose=verbose)
}

write.projection.summary <- function(dir=file.path(getwd(), 'bayesTFR.output'), 
									 output.dir=NULL, revision=NULL, adjusted=FALSE) {
# Writes three prediction summary files, one in a user-friendly format, one in a UN-format,
# and one parameter file.
	pred <- get.tfr.prediction(sim.dir=dir)
	if (is.null(output.dir)) output.dir <- pred$output.directory
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	tfr.write.projection.summary.and.parameters(pred, output.dir, revision=revision, adjusted=adjusted)
}

tfr.write.projection.summary.and.parameters <- function(pred, output.dir, revision=NULL, adjusted=FALSE) {
	# two summary files
	do.write.projection.summary(pred, output.dir, revision=revision, adjusted=adjusted)
	# third file about MCMC parameters
	do.write.parameters.summary(pred, output.dir, adjusted=adjusted)
}

do.write.parameters.summary <- function(pred, output.dir, adjusted=FALSE) {
	meta <- pred$mcmc.set$meta
	tfr.years <- get.tfr.periods(meta)
	tfr <- get.data.imputed(pred)
	if(!is.null(pred$present.year.index)) {
		tfr.years <- tfr.years[1:pred$present.year.index]
		tfr <- tfr[1:pred$present.year.index,]
	}
	all.years <- c(tfr.years, get.prediction.periods(meta, pred$nr.projections+1, present.year.index=pred$present.year.index)[-1])
	nr.data <- pred$nr.projections+dim(tfr)[1]
	# write parameters file
	par.header <- list(country.name="country_name", country.code="country_code", 
					tau.c="TF_time_start_decline", Uc="TF_max", dc="TF_max_decrement",  
					Triangle.c4="TF_end_level", Triangle.c4.low="TF_end_level_low", 
					Triangle.c4.high="TF_end_level_high", Tend="TF_time_end_decline")
	result <- NULL
	precision<-4
	con <- textConnection("sout", "w", local=TRUE) # redirect output (to get rid of coda's messages)
	for (country in get.countries.index(meta)) {
		country.obj <- get.country.object(country, meta, index=TRUE)
		tfr.and.pred.median <- c(tfr[,country], 
								get.median.from.prediction(pred, country.obj$index, 
												country.obj$code, adjusted=adjusted)[-1])
		sink(con, type='message')
		s <- summary(coda.list.mcmc(pred$mcmc.set, country=country.obj$code, 
					par.names=NULL, par.names.cs=tfr.parameter.names.cs(trans=FALSE, back.trans=FALSE), 
					thin=1, burnin=0))
		sink(type='message')
		lambda_c <- find.lambda.for.one.country(tfr.and.pred.median, nr.data)
		result <- rbind(result, c(country.obj$name, country.obj$code, 
			if(meta$tau_c[country.obj$index] > 0) tfr.years[meta$tau_c[country.obj$index]] else -1, #tau_c
			round(s$statistics[paste('U_c',country.obj$code, sep=''),1],precision), # TFR at tau_c
			round(s$statistics[paste('d_c',country.obj$code, sep=''),1],precision),
			round(s$statistics[paste('Triangle_c4_c',country.obj$code, sep=''),1],precision),
			round(s$quantiles[paste('Triangle_c4_c',country.obj$code, sep=''),'2.5%'],precision),
			round(s$quantiles[paste('Triangle_c4_c',country.obj$code, sep=''),'97.5%'],precision),
			all.years[lambda_c]
			))
	}
	close(con)
	colnames(result) <- par.header
	file.suffix <- if(adjusted) '_adjusted' else ''
	file.name <- file.path(output.dir, paste('projection_summary_parameters', file.suffix, '.csv', sep=''))
	write.table(result, file=file.name, sep=',', row.names=FALSE, col.names=TRUE)
	cat('Parameter summary stored into: \n\t\t', file.name, '\n')
}

"get.projection.summary.header" <- function(pred, ...) UseMethod("get.projection.summary.header")
get.projection.summary.header.bayesTFR.prediction <- function(pred, ...) 
	return (list(revision='RevID', variant='VarID', country='LocID', year='TimeID', indicator='IndicatorID', sex='SexID', tfr='Value'))

"get.UN.variant.names" <- function(pred, ...) UseMethod("get.UN.variant.names")
get.UN.variant.names.bayesTFR.prediction <- function(pred, ...) 
	return(c('BHM median', 'BHM80 lower',  'BHM80 upper', 'BHM95 lower',  'BHM95 upper', 'Low', 
					'High', 'Constant fertility'))
					
"get.friendly.variant.names" <- function(pred, ...) UseMethod("get.friendly.variant.names")
get.friendly.variant.names.bayesTFR.prediction <- function(pred, ...)
	return(c('median', 'lower 80', 'upper 80', 'lower 95', 'upper 95', '-0.5child', '+0.5child', 'constant'))

get.wpp.revision.number <- function(pred) {
	wpps <- seq(2008, pred$mcmc.set$meta$wpp.year, by=2)
	lwpps <- length(wpps)
	return(seq(13, length=lwpps)[lwpps])
}

do.write.projection.summary <- function(pred, output.dir, revision=NULL, indicator.id=19, sex.id=3, adjusted=FALSE) {
	cat('Creating summary files ...\n')
	e <- new.env()
	# R check does not like the two lines below; not sure why
	#data('UN_time', envir=e)
	#data('UN_variants', envir=e)
	do.call("data", list('UN_time', envir=e))
	do.call("data", list('UN_variants', envir=e))
	nr.proj <- pred$nr.projections+1
	tfr <- get.data.imputed(pred)
	tfr.years <- get.tfr.periods(pred$mcmc.set$meta)
	if(!is.null(pred$present.year.index)) {
		tfr.years <- tfr.years[1:pred$present.year.index]
		tfr <- tfr[1:pred$present.year.index,]
	}
	ltfr <- dim(tfr)[1] - 1
	nr.proj.all <- nr.proj + ltfr
	pred.period <- get.prediction.periods(pred$mcmc.set$meta, nr.proj, present.year.index=pred$present.year.index)
	header1 <- list(country.name='country_name',  country.code='country_code', variant='variant')
	un.time.idx <- c()
	un.time.label <- as.character(e$UN_time$TLabel)
	l.un.time.label <- length(un.time.label)
	for (i in 1:ltfr) 
		un.time.idx <- c(un.time.idx, which(un.time.label==tfr.years[i])[1])
	for (i in 1:nr.proj) {
		header1[[paste('year', i, sep='')]] <- pred.period[i]
		un.time.idx <- c(un.time.idx, which(un.time.label==pred.period[i]))
	}
	if(is.null(revision)) revision <- get.wpp.revision.number(pred)
	header2 <- get.projection.summary.header(pred)
	UN.variant.names <- get.UN.variant.names(pred)
	friendly.variant.names <- get.friendly.variant.names(pred)
	nr.var <- length(UN.variant.names)
	result1 <- result2 <- NULL
	for (country in 1:get.nr.countries(pred$mcmc.set$meta)) {
		country.obj <- get.country.object(country, pred$mcmc.set$meta, index=TRUE)
		this.tfr <- tfr[,country.obj$index][1:ltfr]
		this.result1 <- cbind(
				country.name=rep(country.obj$name, nr.var), 
				country.code=rep(country.obj$code, nr.var),
				variant=friendly.variant.names)
		median <- get.median.from.prediction(pred, country.obj$index, country.obj$code, adjusted=adjusted)
		proj.result <- rbind(median, 
							   get.traj.quantiles(pred, country.obj$index, country.obj$code, pi=80, adjusted=adjusted),
							   get.traj.quantiles(pred, country.obj$index, country.obj$code, pi=95, adjusted=adjusted))
		if(any(friendly.variant.names == '-0.5child'))
			proj.result <- rbind(proj.result,
					   get.half.child.variant(median))
		proj.result <- round(rbind(proj.result,
							   		rep(median[1], nr.proj)), 4)
		colnames(proj.result) <- grep('year', names(header1), value=TRUE)
		this.result1 <- cbind(this.result1, proj.result)
		result1 <- rbind(result1, this.result1)
		for(ivar in 1:nr.var) {
			result2 <- rbind(result2, cbind(revision=rep(revision, nr.proj.all), 
								   variant=rep(e$UN_variants[e$UN_variants$Vshort==UN.variant.names[ivar],'VarID'], nr.proj.all),
								   country=rep(country.obj$code, nr.proj.all),
								   year=e$UN_time[un.time.idx,'TimeID'],
								   indicator=rep(indicator.id, nr.proj.all),
  									sex=rep(sex.id, nr.proj.all),
								   tfr=c(this.tfr, proj.result[ivar,])))
		}
	}
	colnames(result1)[colnames(result1)==names(header1)] <- header1
	colnames(result2)[colnames(result2)==names(header2)] <- header2
	file.suffix <- if(adjusted) '_adjusted' else ''
	file1 <- paste('projection_summary_user_friendly', file.suffix, '.csv', sep='')
	file2 <- paste('projection_summary', file.suffix, '.csv', sep='')
	write.table(result1, file=file.path(output.dir, file1), sep=',', 
				row.names=FALSE, col.names=TRUE)
	write.table(result2, file=file.path(output.dir, file2), sep=',', 
				quote=FALSE, row.names=FALSE, col.names=TRUE)
	cat('Projection summaries stored into: \n\t\t', 
			file.path(output.dir, file1), '\n\t\t',
			file.path(output.dir, file2), '\n')
}

get.tfr.reconstructed <- function(tfr, meta) {
	tfr_matrix_reconstructed <- tfr
	if(is.null(tfr_matrix_reconstructed)) tfr_matrix_reconstructed <- meta$tfr_matrix_all
	return(tfr_matrix_reconstructed)
}

get.tfr.shift.all <- function(pred, projection.index) {
	# Return shift for all countries in one vector
	meta <- pred$mcmc.set$meta
	nr.countries <- get.nr.countries(meta)
	shift <- rep(0, nr.countries)
	if(is.null(pred$median.shift)) return(shift)
	codes <- meta$regions$country_code
	for(code in names(pred$median.shift)) {
		idx <- which(code == codes)
		shift[idx] <- pred$median.shift[[code]][projection.index]
	}
	return(shift)
}

get.tfr.shift <- function(country.code, pred) {
	if(is.null(pred$median.shift)) return(NULL)
	return(pred$median.shift[[as.character(country.code)]])
}

.bdem.median.shift <- function(pred, type, country, reset=FALSE, shift=0, from=NULL, to=NULL) {
	meta <- pred$mcmc.set$meta
	country.obj <- get.country.object(country, meta=meta)
	if(is.null(country.obj$name)) stop('Country not found.')
	bdem.shift <- do.call(paste('get.', type, '.shift', sep=''), list(country.obj$code, pred))
	pred.years <- as.numeric(dimnames(pred$quantiles)[[3]])
	nr.proj <- pred$nr.projections+1 
	if(is.null(from)) from <- pred.years[2]
	if(is.null(to)) to <- pred.years[nr.proj]
	which.years <- (pred.years >= from) & (pred.years <= to)
	all.years <- FALSE
	if(reset) { # reset to 0
		if (sum(which.years) >= nr.proj-1) {
			bdem.shift <- NULL # reset everything
			all.years <- TRUE
		} else bdem.shift[which.years] <- 0 
		action <- 'reset to BHM values'
	} else { # shift by given amount
		if(is.null(bdem.shift)) bdem.shift <- rep(0, nr.proj)
		bdem.shift[which.years] <- bdem.shift[which.years] + shift
		action <- 'modified'
	}
	if(sum(bdem.shift) == 0) bdem.shift <- NULL
	pred$median.shift[[as.character(country.obj$code)]] <- bdem.shift
	cat('\nMedian of', country.obj$name, action, 
		if(all.years) 'for all years' else c('for years', pred.years[which.years]), '.\n')
	return(pred)
}

tfr.median.reset <- function(sim.dir, countries) {
	for(country in countries) pred <- tfr.median.shift(sim.dir, country, reset=TRUE)
	invisible(pred)
}

tfr.median.shift <- function(sim.dir, country, reset=FALSE, shift=0, from=NULL, to=NULL) {
	pred <- get.tfr.prediction(sim.dir)
	pred <- .bdem.median.shift(pred, type='tfr', country=country, reset=reset, 
				shift=shift, from=from, to=to)
	store.bayesTFR.prediction(pred)
	invisible(pred)
}

.bdem.median.set <- function(pred, type, country, values, years=NULL) {
	meta <- pred$mcmc.set$meta
	country.obj <- get.country.object(country, meta=meta)
	if(is.null(country.obj$name)) stop('Country not found.')
	bdem.shift <- do.call(paste('get.', type, '.shift', sep=''), list(country.obj$code, pred))
	pred.years <- as.numeric(dimnames(pred$quantiles)[[3]])
	nr.proj <- pred$nr.projections+1 
	if(is.null(years)) years <- pred.years[2:nr.proj]
	mid.years <- cut(years, labels=pred.years, 
					breaks=seq(from=pred.years[1]-3, to=pred.years[nr.proj]+2, by=5))
	which.years <- is.element(pred.years, mid.years)
	lvalues <- length(values)
	if(lvalues > sum(which.years)) {
		start <- which.max(cumsum(which.years))+1
		end <- min(start + lvalues - sum(which.years), nr.proj)
		if(start > end) stop ('Mismatch in length of values and years.')
		which.years[start:end] <- TRUE
	}
	if(lvalues < sum(which.years)) { 
		start <- which(cumsum(which.years)==lvalues)[1]+1
		which.years[start:nr.proj] <- FALSE
	}
	if(is.null(bdem.shift)) bdem.shift <- rep(0, nr.proj)
	medians <- pred$quantiles[country.obj$index, '0.5',]
	bdem.shift[which.years] <- values - medians[which.years]
	if(sum(bdem.shift) == 0) bdem.shift <- NULL
	pred$median.shift[[as.character(country.obj$code)]] <- bdem.shift
	cat('\nMedian of', country.obj$name, 'modified for years', pred.years[which.years], '\n')
	return(pred)
}

tfr.median.set <- function(sim.dir, country, values, years=NULL) {
	pred <- get.tfr.prediction(sim.dir)
	pred <- .bdem.median.set(pred, 'tfr', country=country, values=values, years=years)
	store.bayesTFR.prediction(pred)
	invisible(pred)
}

tfr.median.adjust <- function(sim.dir, countries, factor1=2/3, factor2=1/3, forceAR1=FALSE) {
	pred <- get.tfr.prediction(sim.dir)
	if (is.null(pred)) stop('No valid prediction in ', sim.dir)
	mcmc.set <- pred$mcmc.set
	if(is.null(countries)) {
		cat('\nNo countries given. Nothing to be done.\n')
		return(invisible(pred))
	}
	codes <- c()
	for(country in countries) codes <- c(codes, get.country.object(country, mcmc.set$meta)$code)
	countries.idx <- which(is.element(mcmc.set$meta$regions$country_code, codes))
	if(length(countries.idx) == 0) {
		cat('\nNo valid countries given. Nothing to be done.\n')
		return(invisible(pred))	
	}
	m3.set <- if(pred$use.tfr3) get.tfr3.mcmc(sim.dir) else NULL
	new.pred <- make.tfr.prediction(mcmc.set, start.year=pred$start.year, end.year=pred$end.year, replace.output=FALSE,
									nr.traj=NULL, burnin=0, 
									use.tfr3=pred$use.tfr3, mcmc3.set=m3.set, burnin3=pred$burnin3,
									mu=pred$mu, rho=pred$rho, sigmaAR1=pred$sigmaAR1, 
									countries=countries.idx, adj.factor1=factor1, adj.factor2=factor2,
									forceAR1=forceAR1, save.as.ascii=0, output.dir=NULL,
									write.summary.files=FALSE, is.mcmc.set.thinned=TRUE, 
									write.trajectories=FALSE, verbose=FALSE)
	new.means <- new.pred$traj.mean.sd[,1,2:dim(new.pred$traj.mean.sd)[3]]
	for(icountry in 1:length(countries)) {
		tfr.median.set(sim.dir, countries[icountry], new.means[get.country.object(countries[icountry], mcmc.set$meta)$index,])
	}
	# reload adjusted prediction
	invisible(get.tfr.prediction(sim.dir))
}

tfr.correlation <- function(meta, cor.pred=NULL, low.coeffs=c(0.11, 0.26, 0.05, 0.09),
								high.coeffs=c(0.05, 0.06, 0.00, 0.02), kappa=5) {
	# kappa is the TFR cut-off for low and high coefs 
	nr_countries <- get.nr.countries(meta)
	low.eps.cor <- matrix(NA,nrow=nr_countries,ncol=nr_countries)
	high.eps.cor <-  matrix(NA,nrow=nr_countries,ncol=nr_countries)
	country.codes <- meta$regions$country_code
	if(is.null(cor.pred)) {
		e <- new.env()
		# the following code causes a NOTE in R check; not sure why
		#data("correlation_predictors", envir=e)
		do.call("data", list("correlation_predictors", envir=e))
		cor.pred <- e$correlation_predictors
	}
	for(i in 1:(nr_countries-1)) {
		i_is_in <- (cor.pred[,1] == country.codes[i]) | (cor.pred[,2] == country.codes[i])
		if(sum(i_is_in)==0) {
			warning('No records found for country ', country.codes[i])
			next
		}
  		for(j in (i+1):nr_countries) {
        	pred.row <- which(i_is_in & ((cor.pred[,1] == country.codes[j]) | (cor.pred[,2] == country.codes[j])))
        	if(length(pred.row) <= 0) {
        		warning('No records found for pair ', paste(country.codes[c(i,j)], collapse=', '))
        		next
        	}
        	if(length(pred.row)>1) pred.row <- pred.row[1]
        	reg.values <- c(1,as.numeric(cor.pred[pred.row,-c(1:2)]))
        	low.eps.cor[i,j] <- low.eps.cor[j,i] <- sum(reg.values*low.coeffs)
        	high.eps.cor[i,j] <- high.eps.cor[j,i] <- sum(reg.values*high.coeffs)        
        	if(is.na(low.eps.cor[i,j]) || is.na(high.eps.cor[i,j]))
        		warning('Correlation resulted in NA for pair ', paste(country.codes[c(i,j)], collapse=', '), immediate.=TRUE)  
  		}
	}
	diag(low.eps.cor) <- 1
	diag(high.eps.cor) <- 1
	return(list(low=low.eps.cor, high=high.eps.cor, kappa=kappa))
}

tfr.correlation.subnat <- function(mcmc.set, 
									method=c('meth12', 'meth11', 'meth10', 'meth9', 'meth8', 'meth7', 'meth6', 'meth5', 'bayes.mean', 'bayes.mode', 'median', 'mean'), 
									burnin=0, thin=1, burnin3=0, 
									use.external.phase3=FALSE, use.phase3.from = NULL, 
									ar1pars=NULL, kappa=5, cor.start.year=NULL, verbose=FALSE, ...) {
	meta <- mcmc.set$meta
	method <- match.arg(method)
	errs <- is.low <- matrix(NA, nrow=meta$T_end, ncol=meta$nr_countries, 
						dimnames=list(dimnames(meta$tfr_matrix)[[1]], meta$regions$country_code))
	m3 <- NULL
	has.phase3 <- has.tfr3.mcmc(mcmc.set$meta$output.dir)
	if(has.phase3) {
		m3 <-  get.tfr3.mcmc(sim.dir=mcmc.set$meta$output.dir)
		if(length(m3$mcmc.list)==0) {
			m3 <- NULL
			has.phase3 <- FALSE
		}
	}
	if(!has.phase3 && use.external.phase3) {
		m3 <- get.tfr3.mcmc(sim.dir=use.phase3.from)
		if (has.tfr.prediction(sim.dir=use.phase3.from)) 
			burnin3 <- get.tfr.prediction(use.phase3.from)$burnin3
	}	
	if(!is.null(m3)) { 
		iter3 <- get.total.iterations(m3$mcmc.list, burnin3)
		thin3 <- max(floor(iter3/2000),1)			
	}
	.get.err <- function(i, mu, rho, sigma) return((tfr.all[i]-mu-rho*(tfr.all[i-1]-mu))/sigma)
	.get.err.distr <- function(i, traces) 
						return(apply(traces, 1, function(x) return(.get.err(i, x[1], x[2], x[3]))))
	# get prediction errors
	for (country in 1:meta$nr_countries) {
		country.obj <- get.country.object(country, meta, index=TRUE)
		tfr.all <- bayesTFR:::get.observed.tfr(country, meta, 'tfr_matrix_observed')
		is.low[,country] <- tfr.all < kappa
		dl.obs.idx <- if(max(meta$tau_c[country],1) >= meta$lambda_c[country]) c() 
						else seq(max(meta$tau_c[country],1), meta$lambda_c[country])
		start.idx <- if(is.null(cor.start.year)) 1 else which(as.integer(names(tfr.all)) >= cor.start.year)[1]
		dl.obs.idx <- dl.obs.idx[dl.obs.idx>=start.idx]
		if(length(dl.obs.idx) > 0) { # Phase II
			tfr <- tfr.all[dl.obs.idx]
			not.na.idx <- which(!is.na(tfr))
			if(length(not.na.idx) > 1) { # needs at least two points for decrements
				tfr <- tfr[not.na.idx]
				decr <- -diff(tfr)
				dlc.sigma <- tfr.get.dlcurves(tfr, mcmc.set$mcmc.list, country.obj$code, country, 
									burnin=burnin, nr.curves=0, return.sigma=TRUE)
				terr <- apply(-dlc.sigma$dl[,-1, drop=FALSE], 1, "+", decr)
				if(length(decr) > 1) terr <- t(terr)
				dl.err <- terr/dlc.sigma$sigma[,-1, drop=FALSE]
				res.idx <- dl.obs.idx[not.na.idx]
				errs[res.idx[-1],country] <- apply(dl.err, 2, mean)
				#if(country.obj$code==549) stop('')
			}
		}
		if(meta$lambda_c[country] < meta$T_end_c[country])  { # Phase III
			p3.idx <- seq(meta$lambda_c[country]+1, meta$T_end_c[country])
			p3.idx <- p3.idx[p3.idx >= start.idx]		
			if(!is.null(m3)) {
				if(!use.external.phase3 && is.element(country, m3$meta$id_phase3))
					tr <- get.tfr3.parameter.traces.cs(m3$mcmc.list, country.obj, c('mu.c', 'rho.c'), 
								burnin=burnin3, thin=thin3)
				else 
					tr <- get.tfr3.parameter.traces(m3$mcmc.list, c('mu', 'rho'), 
								burnin=burnin3, thin=thin3)
				tr <- cbind(tr, get.tfr3.parameter.traces(m3$mcmc.list, c('sigma.eps'), 
										burnin=burnin3, thin=thin3))
				eps.k <- NULL 
				for(j in 1:length(p3.idx)) 
					eps.k <- cbind(eps.k, .get.err.distr(p3.idx[j], tr))
				errs[p3.idx, country] <- apply(eps.k, 2, mean)
			} else  # use fixed parameters
				errs[p3.idx, country] <- sapply(p3.idx, .get.err, mu=ar1pars$mu, rho=ar1pars$rho, sigma=ar1pars$sigmaAR1)
		}
	}
	ind <- apply(errs, 1, function(x) all(is.na(x)))
	errs <- errs[!ind,]
	is.low <- is.low[!ind,]
	cormat <- switch(method,
				bayes.mean = cor.bayes(errs, is.low, method='mean', verbose=verbose, ...),
				bayes.mode = cor.bayes(errs, is.low, method='mode', verbose=verbose, ...),
				meth5 = cor.bayes(errs, is.low, method='mean', scale.errors=TRUE, verbose=verbose, ...),
				meth6 = cor.bayes(errs, is.low, method='mode', scale.errors=TRUE, verbose=verbose, ...),
				meth7 = cor.modified(errs, is.low, verbose=verbose, ...),
				meth8 = cor.moments2(errs, is.low, verbose=verbose, ...),
				meth9 = cor.method9(errs, is.low, verbose=verbose, ...),
				meth10 = cor.bayes.meth10(errs, is.low, arcsin.prior=FALSE, verbose=verbose, ...),
				meth11 = cor.bayes.meth10(errs, is.low, arcsin.prior=TRUE, verbose=verbose, ...),
				meth12 = cor.constant(errs, verbose=verbose, ...),
				cor.moments(errs, is.low, method=method, verbose=verbose, ...)
				)
	return(c(cormat, list(kappa=kappa)))
}

tfr.correlation.subnat2 <- function(mcmc.set, mcmc3.set, meta, country,
									method=c('meth12', 'meth11', 'meth10', 'meth9', 'meth8'), 
									burnin=0, thin=1, burnin3=0, 
									ar1pars=NULL, kappa=5, cor.start.year=NULL, verbose=FALSE, ...) {
	method <- match.arg(method)
	errs <- is.low <- matrix(NA, nrow=meta$T_end, ncol=meta$nr_countries, 
						dimnames=list(dimnames(meta$tfr_matrix)[[1]], meta$regions$country_code))
	country.obj <- get.country.object(country, mcmc.set$meta) # from world meta
	has.phase3 <- !is.null(mcmc3.set)
	if(has.phase3) {
		if(length(mcmc3.set$mcmc.list)==0) {
			mcmc3.set <- NULL
			has.phase3 <- FALSE
		} else {
			iter3 <- get.total.iterations(mcmc3.set$mcmc.list, burnin3)
			thin3 <- max(floor(iter3/2000),1)
		}
		if(is.element(country.obj$index, mcmc3.set$meta$id_phase3))
			tr3 <- get.tfr3.parameter.traces.cs(mcmc3.set$mcmc.list, country.obj, c('mu.c', 'rho.c'), 
								burnin=burnin3, thin=thin3)
		else 
			tr3 <- get.tfr3.parameter.traces(mcmc3.set$mcmc.list, c('mu', 'rho'), 
								burnin=burnin3, thin=thin3)
		tr3 <- cbind(tr3, get.tfr3.parameter.traces(mcmc3.set$mcmc.list, c('sigma.eps'), 
										burnin=burnin3, thin=thin3))
	}
											
	.get.err <- function(i, mu, rho, sigma) return((tfr.all[i]-mu-rho*(tfr.all[i-1]-mu))/sigma)
	.get.err.distr <- function(i, traces) 
						return(apply(traces, 1, function(x) return(.get.err(i, x[1], x[2], x[3]))))
	# get prediction errors
	for (reg in 1:meta$nr_countries) {
		reg.obj <- get.country.object(reg, meta, index=TRUE)
		tfr.all <- bayesTFR:::get.observed.tfr(reg, meta, 'tfr_matrix_observed')
		is.low[,reg] <- tfr.all < kappa
		dl.obs.idx <- if(max(meta$tau_c[reg],1) >= meta$lambda_c[reg]) c() 
						else seq(max(meta$tau_c[reg],1), meta$lambda_c[reg])
		start.idx <- if(is.null(cor.start.year)) 1 else which(as.integer(names(tfr.all)) >= cor.start.year)[1]
		dl.obs.idx <- dl.obs.idx[dl.obs.idx >= start.idx]
		if(length(dl.obs.idx) > 0) { # Phase II
			tfr <- tfr.all[dl.obs.idx]
			not.na.idx <- which(!is.na(tfr))
			if(length(not.na.idx) > 1) { # needs at least two points for decrements
				tfr <- tfr[not.na.idx]
				decr <- -diff(tfr)
				dlc.sigma <- tfr.get.dlcurves(tfr, mcmc.set$mcmc.list, country.obj$code, country.obj$index, 
									burnin=burnin, nr.curves=0, return.sigma=TRUE)
				terr <- apply(-dlc.sigma$dl[,-1, drop=FALSE], 1, "+", decr)
				if(length(decr) > 1) terr <- t(terr)
				dl.err <- terr/dlc.sigma$sigma[,-1, drop=FALSE]
				res.idx <- dl.obs.idx[not.na.idx]
				errs[res.idx[-1],reg] <- apply(dl.err, 2, mean)
				#if(country.obj$code==549) stop('')
			}
		}
		if(meta$lambda_c[reg] < meta$T_end_c[reg])  { # Phase III
			p3.idx <- seq(meta$lambda_c[reg]+1, meta$T_end_c[reg])
			p3.idx <- p3.idx[p3.idx >= start.idx]		
			if(has.phase3) {			
				eps.k <- NULL 
				for(j in 1:length(p3.idx)) 
					eps.k <- cbind(eps.k, .get.err.distr(p3.idx[j], tr3))
				errs[p3.idx, reg] <- apply(eps.k, 2, mean)
			} else  # use fixed parameters
				errs[p3.idx, reg] <- sapply(p3.idx, .get.err, mu=ar1pars$mu, rho=ar1pars$rho, sigma=ar1pars$sigmaAR1)
		}
	}
	ind <- apply(errs, 1, function(x) all(is.na(x)))
	errs <- errs[!ind,]
	is.low <- is.low[!ind,]
	cormat <- switch(method,
				bayes.mean = cor.bayes(errs, is.low, method='mean', verbose=verbose, ...),
				bayes.mode = cor.bayes(errs, is.low, method='mode', verbose=verbose, ...),
				meth5 = cor.bayes(errs, is.low, method='mean', scale.errors=TRUE, verbose=verbose, ...),
				meth6 = cor.bayes(errs, is.low, method='mode', scale.errors=TRUE, verbose=verbose, ...),
				meth7 = cor.modified(errs, is.low, verbose=verbose, ...),
				meth8 = cor.moments2(errs, is.low, verbose=verbose, ...),
				meth9 = cor.method9(errs, is.low, verbose=verbose, ...),
				meth10 = cor.bayes.meth10(errs, is.low, arcsin.prior=FALSE, verbose=verbose, ...),
				meth11 = cor.bayes.meth10(errs, is.low, arcsin.prior=TRUE, verbose=verbose, ...),
				meth12 = cor.constant(errs, verbose=verbose, ...),
				cor.moments(errs, is.low, method=method, verbose=verbose, ...)
				)
	return(c(cormat, list(kappa=kappa)))
}


"get.data.imputed" <- function(pred, ...) UseMethod("get.data.imputed")

get.data.imputed.bayesTFR.prediction <- function(pred, ...)
	return(get.tfr.reconstructed(pred$tfr_matrix_reconstructed, pred$mcmc.set$meta))
	
"get.data.imputed.for.country" <- function(pred, country.index, ...) UseMethod("get.data.imputed.for.country")

get.data.imputed.for.country.bayesTFR.prediction <- function(pred, country.index, ...)
	return(get.observed.with.supplemental(country.index, pred$tfr_matrix_reconstructed, pred$mcmc.set$meta$suppl.data, 'tfr_matrix_all'))
	