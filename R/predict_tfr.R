tfr.predict <- function(mcmc.set=NULL, end.year=2100,
						sim.dir=file.path(getwd(), 'bayesTFR.output'),
						replace.output=FALSE,
						nr.traj = NULL, thin = NULL, burnin=2000, use.diagnostics=FALSE,
						save.as.ascii=1000, output.dir = NULL,
						low.memory=TRUE,
						# parameters for the AR1 stuff
						mu=2.1, rho=0.8859, sigmaAR1=0.1016,
						seed=NULL, verbose=TRUE) {
	if(!is.null(mcmc.set)) {
		if (class(mcmc.set) != 'bayesTFR.mcmc.set') {
			stop('Wrong type of mcmc.set. Must be of type bayesTFR.mcmc.set.')
			}
	} else {		
		mcmc.set <- get.tfr.mcmc(sim.dir, low.memory=low.memory, verbose=verbose)
	}
	if (is.null(rho) || is.na(rho) || is.null(sigmaAR1) || is.na(sigmaAR1)) {
		res <- get.ar1.parameters(mu = mu, mcmc.set$meta)
		if (is.null(rho) || is.na(rho)) 
			rho <- res$rho
		if(is.null(sigmaAR1) || is.na(sigmaAR1))
			sigmaAR1 <- res$sigmaAR1
	}
	if(verbose)
		cat('\nAR(1) parameters used: mu=', mu, ', rho=', rho, ', sigma=', sigmaAR1, '\n')
	if(!is.null(seed)) set.seed(seed)
	
	# Get argument settings from existing convergence diagnostics
	if(use.diagnostics) {
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
		nr.traj.idx <- if(sum(larger2T)>0) (1:ldiag)[larger2T][which.min(use.nr.traj[larger2T])] else (1:ldiag)[traj.is.notna][which.max(use.nr.traj[traj.is.notna])]
		nr.traj <- use.nr.traj[nr.traj.idx]
		burnin <- use.burnin[nr.traj.idx]
		if(verbose)
			cat('\nUsing convergence settings: nr.traj=', nr.traj, ', burnin=', burnin, '\n')
	}
	invisible(make.tfr.prediction(mcmc.set, end.year=end.year, replace.output=replace.output,  
					nr.traj=nr.traj, burnin=burnin, thin=thin,
					mu=mu, rho=rho,  sigmaAR1 = sigmaAR1,
					save.as.ascii=save.as.ascii,
					output.dir=output.dir, verbose=verbose))			
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
	new.pred <- make.tfr.prediction(mcmc.set, end.year=pred$end.year, replace.output=FALSE,
									nr.traj=pred$nr.traj, burnin=pred$burnin,
									mu=pred$mu, rho=pred$rho, sigmaAR1=pred$sigmaAR1, 
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


make.tfr.prediction <- function(mcmc.set, end.year=2100, replace.output=FALSE,
								nr.traj = NULL, burnin=0, thin = NULL, 
								mu=2.1, rho=0.9057, sigmaAR1 = 0.0922, countries = NULL,
								adj.factor1=NA, adj.factor2=0, forceAR1=FALSE,
							    save.as.ascii=1000, output.dir = NULL, write.summary.files=TRUE, 
							    is.mcmc.set.thinned=FALSE, force.creating.thinned.mcmc=FALSE,
							    write.trajectories=TRUE,
							    verbose=verbose){
	# if 'countries' is given, it is an index
	# sigmaAR1 can be a vector. The last element will be repeated up to nr.projections

	nr_project <- ceiling((end.year - mcmc.set$meta$present.year)/5)
	suppl.T <- if(!is.null(mcmc.set$meta$suppl.data$regions)) mcmc.set$meta$suppl.data$T_end else 0
#	if (verbose)
		cat('\nPrediction from', mcmc.set$meta$present.year, 
			'(excl.) until', end.year, '(i.e.', nr_project, 'projections)\n\n')
	l.sigmaAR1 <- length(sigmaAR1)
	sigma.end <- rep(sigmaAR1[l.sigmaAR1], nr_project + mcmc.set$meta$T_end-l.sigmaAR1)
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
	if (is.null(output.dir)) output.dir <- mcmc.set$meta$output.dir
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
											
	prediction.countries <- if(is.null(countries)) 1:mcmc.set$meta$nr_countries else countries
	nr_countries <- mcmc.set$meta$nr_countries
	tfr_matrix_reconstructed <- get.tfr.reconstructed(mcmc.set$meta$tfr_matrix_observed, mcmc.set$meta)
	ltfr_matrix <- dim(tfr_matrix_reconstructed)[1]
	ltfr_matrix.all <- ltfr_matrix + suppl.T
	
	#quantiles.to.keep <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
	#keep these defaults for checking the out-of-sample projections
    quantiles.to.keep <- c(0,0.025,0.05,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,0.95,0.975,1)
	PIs_cqp <- array(NA, c(nr_countries, length(quantiles.to.keep), nr_project+1))
	dimnames(PIs_cqp)[[2]] <- quantiles.to.keep
	proj.middleyears <- seq(max(as.numeric(dimnames(tfr_matrix_reconstructed)[[1]])), by=5, length=nr_project+1)
	dimnames(PIs_cqp)[[3]] <- proj.middleyears
	mean_sd <- array(NA, c(nr_countries, 2, nr_project+1))
	hasNAs <- rep(FALSE, nr_simu)
	adjust.true <- if(is.na(adj.factor1)) FALSE else TRUE
	if (verbose) cat('Load parameters mean_eps_tau and sd_eps_tau.\n')
    tau.par.names <- c('mean_eps_tau', 'sd_eps_tau')
    tau.par.values <- get.tfr.parameter.traces(load.mcmc.set$mcmc.list, tau.par.names, burnin=0)

	if (verbose) cat('Load hierarchical parameters.\n')
	alpha.vars <- paste('alpha_',1:3, sep='')
	delta.vars <- paste('delta_',1:3, sep='')
	other.vars <- c('chi', 'psi', 'Triangle4', 'delta4')
	cs.par.values_hier <- get.tfr.parameter.traces(load.mcmc.set$mcmc.list, 
										c(alpha.vars, delta.vars, other.vars), burnin=0)
										
	const_sd_dummie_Tc <- matrix(0, mcmc.set$meta$T_end+suppl.T, nr_countries)
	mid.years <- as.integer(c(if(suppl.T > 0) rownames(mcmc.set$meta$suppl.data$tfr_matrix) else c(), rownames(tfr_matrix_reconstructed)))
	const_sd_dummie_Tc[mid.years < 1975,] <- 1
	
	country.counter <- 0
	status.for.gui <- paste('out of', length(prediction.countries), 'countries.')
	gui.options <- list()
	#########################################
	for (country in prediction.countries){
	#########################################
		if(getOption('bDem.TFRpred', default=FALSE)) {
			# This is to unblock the GUI, if the run is invoked from bayesDem
			# and pass info about its status
			# In such a case the gtk libraries are already loaded
			country.counter <- country.counter + 1
			gui.options$bDem.TFRpred.status <- paste('finished', country.counter, status.for.gui)
			unblock.gtk('bDem.TFRpred', gui.options)
		}

		country.obj <- get.country.object(country, mcmc.set$meta, index=TRUE)
		if (verbose) {			
 			cat('TFR projection for country', country, country.obj$name, 
 						'(code', country.obj$code, ')\n')
 		}
	cs.par.values <- c()
	if (is.element(country,mcmc.set$meta$id_DL)){
		U.var <- paste('U_c', country.obj$code, sep='')
		d.var <- paste('d_c', country.obj$code, sep='')
		Triangle_c4.var <- paste("Triangle_c4_c", country.obj$code, sep = "")
		gamma.vars <- paste('gamma_',1:3,'_c', country.obj$code, sep='')
		if(country <= mcmc.set$meta$nr_countries_estimation) {
			cs.par.values <- get.tfr.parameter.traces.cs(load.mcmc.set$mcmc.list, country.obj, 
								tfr.parameter.names.cs(trans=FALSE), burnin=0)
		} else {
			# if it's country that was not included in the estimation, determine the posterior index
			# again, since the MCMC might be of different length
			if (is.element(country.obj$code, load.mcmc.set$meta$regions$country_code)) {
				cs.par.values <- get.tfr.parameter.traces.cs(load.mcmc.set$mcmc.list, country.obj, 
								tfr.parameter.names.cs(trans=FALSE), burnin=0)
			} else { # there are no thinned traces for this country, use the full traces 
				cs.par.values <- get.tfr.parameter.traces.cs(mcmc.set$mcmc.list, country.obj, 
								tfr.parameter.names.cs(trans=FALSE), burnin=burnin)
					selected.simu <- get.thinning.index(nr_simu, dim(cs.par.values)[1])
					if (length(selected.simu$index) < nr_simu)
						selected.simu$index <- sample(selected.simu$index, nr_simu, replace=TRUE)
					cs.par.values <- cs.par.values[selected.simu$index,]
			}
		}
		pc_si <- matrix(NA, nr_simu, 3)
		for (i in 1:3){
 			pc_si[,i] <- exp(cs.par.values[,gamma.vars[i]])/
 								apply(exp(cs.par.values[,gamma.vars]), 1,sum)
		}
		theta_si <- cbind(
				(cs.par.values[,U.var] - cs.par.values[, Triangle_c4.var])*pc_si[,1],
				(cs.par.values[,U.var] - cs.par.values[, Triangle_c4.var])*pc_si[,2],
				(cs.par.values[,U.var] - cs.par.values[, Triangle_c4.var])*pc_si[,3],
				cs.par.values[, Triangle_c4.var], 
				cs.par.values[,d.var]) 
	} else { #Tistau countries
             # sample decline parameters from the hier distributions                 
		Triangle4_tr_s = rnorm(1, mean = cs.par.values_hier[,'Triangle4'], 
									sd = cs.par.values_hier[, 'delta4'])
		Triangle_c4_s <- ( mcmc.set$meta$Triangle_c4.up*exp(Triangle4_tr_s) + mcmc.set$meta$Triangle_c4.low)/(1+exp(Triangle4_tr_s))

		# need U and Triangle_c4 in cs... later in loop for start of phase III and prior on f_t
		cs.par.values = rep(get.observed.tfr(country, mcmc.set$meta, 'tfr_matrix_all')[mcmc.set$meta$tau_c[country]], nr_simu)
		Triangle_c4.var <- 'Triangle_c4'
		U.var <- 'U'
		cs.par.values = cbind(cs.par.values, Triangle_c4_s)
		colnames(cs.par.values) = c(U.var, Triangle_c4.var)

		d_tr_s = rnorm(1, mean = cs.par.values_hier[,'chi'], 
							sd = cs.par.values_hier[,'psi'])
		d_s =  (mcmc.set$meta$d.up*(exp(d_tr_s) + mcmc.set$meta$d.low)/(1+exp(d_tr_s)))
		#print(d_s) # check!

		gamma_si = matrix(NA, nr_simu, 3)
		for (i in 1:3){
 			gamma_si[,i] <- rnorm(1, mean = cs.par.values_hier[,alpha.vars[i]], 
 									sd = cs.par.values_hier[,delta.vars[i]])
		}

		pc_si = matrix(NA, nr_simu, 3)
		for (i in 1:3) pc_si[,i] <- exp(gamma_si[,i])/apply(exp(gamma_si), 1,sum)
		theta_si <- cbind((cs.par.values[,U.var] - Triangle_c4_s)*pc_si[,1],
                          (cs.par.values[,U.var] - Triangle_c4_s)*pc_si[,2],
                          (cs.par.values[,U.var] - Triangle_c4_s)*pc_si[,3],
                          Triangle_c4_s, d_s) 
	}
		all.tfr <- get.observed.tfr(country, mcmc.set$meta, 'tfr_matrix_all')
		this.T_end <- mcmc.set$meta$T_end_c[country]
		nmissing <- ltfr_matrix.all - this.T_end
		this.nr_project <- nr_project + nmissing
			
		#### projections might start before 1975 ...
		dummie_c1975 = rep(0, this.nr_project)
		# if this.T_end is before 1975, need to add const_sd to sd of distortion terms
		if (sum(const_sd_dummie_Tc[this.T_end:(mcmc.set$meta$T_end+suppl.T),country]) !=0){
 			# gives number of add. periods with const
 			dummie_c1975[1:sum(const_sd_dummie_Tc[this.T_end:(mcmc.set$meta$T_end+suppl.T),country])] <- 1
		}
                
		nr_obs_in_phaseIII = (this.T_end - mcmc.set$meta$lambda_c[country])
		sigma_t  = sigmas_all[(nr_obs_in_phaseIII+1):(nr_obs_in_phaseIII+this.nr_project)]
		if (verbose & nmissing > 0) 
			cat('\t', nmissing, 'data points reconstructed.\n')
		# the projections
		f_ps <- matrix(NA, this.nr_project+1,nr_simu)
		# fs also includes last estimate!
		f_ps[1,] <- all.tfr[this.T_end]
		if(adjust.true) {
			D11 <- (all.tfr[this.T_end-1] - all.tfr[this.T_end])
 		  	E11 <- (all.tfr[this.T_end-1] - (mu + rho*(all.tfr[this.T_end-1]-mu)))
 		  	S11pIII <- D11 - E11
 		}
 		S11 <- 0

		for (s in 1:nr_simu){
			year <- 2
			if(forceAR1) repl_stop <- TRUE
			else
  				repl_stop <- ifelse(
  					(mcmc.set$meta$lambda_c[country] != this.T_end) || 
                	((min(all.tfr, na.rm=TRUE) <= cs.par.values[s, Triangle_c4.var]) && 
                 	(all.tfr[this.T_end] > all.tfr[this.T_end-1])), TRUE, FALSE)

			if(adjust.true) {
			 	if(!repl_stop) {
					# country in Phase II
            		d11 <- DLcurve(theta_si[s,], all.tfr[this.T_end-1], mcmc.set$meta$dl.p1, mcmc.set$meta$dl.p2)
		  			S11 <- D11 - d11
		  		} else {
					# country already in Phase III in last obs. period
		  			S11 <- S11pIII
		  		}
		  	}
	       	# vector W with the weight for the first two periods:
            W <- c(adj.factor1, adj.factor2, rep(0, this.nr_project-2))
            W[is.na(W)] <- 0
            # first projection year
                                
			if (is.element(country, mcmc.set$id_Tistau)){
				### need to start with sampling eps_tau!!
       			 while(TRUE) {
                 	eps_tmin1 <- rnorm(1, tau.par.values[s, 'mean_eps_tau'], tau.par.values[s, 'sd_eps_tau'])
                    fps <- (f_ps[year -1,s]- DLcurve(theta_si[s,], f_ps[year -1,s], 
                                                        mcmc.set$meta$dl.p1, mcmc.set$meta$dl.p2)+ eps_tmin1 
                                                        - W[year-1]*S11)

                    #if( fps > 0 && fps <= cs.par.values[s,U.var] ) break
                    if( fps > 0.5 && fps <= cs.par.values[s,U.var] ) break
                 }
        		f_ps[year,s] <- fps
        		year = year +1 # no repl stop possible if decline has just started!
			}
 			while (!repl_stop && year <(this.nr_project+2)){ 
 				# year goes upto nr_project+1
  				# sample the distortion
  				sigma_eps <- var.par.values[s,'sigma0'] + 
  					(f_ps[year -1,s] - var.par.values[s,'S_sd'])*
  						ifelse(f_ps[year -1,s] > var.par.values[s,'S_sd'], 
  						-var.par.values[s,'a_sd'], var.par.values[s,'b_sd'])
  				sigma_eps <- max(sigma_eps, mcmc.set$meta$sigma0.min)
  				while(TRUE) {
  					eps_tmin1 <- rnorm(1, 0, sigma_eps)
  					fps <- (f_ps[year -1,s]- DLcurve(theta_si[s,], f_ps[year -1,s], 
  								mcmc.set$meta$dl.p1, mcmc.set$meta$dl.p2)+ eps_tmin1 - W[year-1]*S11)
  					#if( fps > 0 && fps <= cs.par.values[s,U.var] ) break
  					if( fps > 0.5 && fps <= cs.par.values[s,U.var] ) break
  				}
  				f_ps[year,s] <- fps
  				if ((min(f_ps[1:year, s]) <= cs.par.values[s, Triangle_c4.var]) && 
                                                           (f_ps[year, s] > f_ps[year-1, s]))
                	repl_stop <- TRUE 
				if(adjust.true && year == 2){
					# if a country with adjustment enters Phase III in second proj. step (corresponds to year ==3)
					# then the adjustment needs to be changed, based on observed diff in last proj step and AR(1) decrement
		  			S11 <- ( f_ps[year-1,s] - f_ps[year,s]) - ( f_ps[year-1,s] - (mu + rho*( f_ps[year-1,s]-mu)))
			 	}
			 	year = year + 1
 			} # end while loop, continue if repl_stop ==TRUE and year<nr_project+2
 			
 			if (repl_stop && year <(this.nr_project+2)){
 				# use sigma_t
                index_sigma <- 0
 				for (period in year:(this.nr_project+1)){
 					index_sigma <- index_sigma + 1
 					while (TRUE){
						f_ps[period,s] <- (mu + rho*
								(f_ps[period-1,s] - mu) + 
								rnorm(1, 0, sigma_t[index_sigma]) - W[period-1]*S11)
						if (f_ps[period, s] > 0.5 )   break
					}
				}
 			}# end AR(1) loop
		} # end simu loop
		#isnotNA <- apply(1-(f_ps<0), 2, prod) 
		isnotNA <- apply(1-(f_ps<0.5), 2, prod) 
		isnotNA <- ifelse(is.na(isnotNA),0,isnotNA)
		f_ps[,isnotNA==0] <- NA
		if (nmissing > 0) {
			f_ps_future <- f_ps[(nmissing+1):nrow(f_ps),]
			f_ps_future[1,] <- quantile(f_ps_future[1,], 0.5, na.rm = TRUE)
			tfr_matrix_reconstructed[(this.T_end-suppl.T+1):ltfr_matrix,country] <- apply(matrix(f_ps[2:(nmissing+1),], nrow=nmissing), 
												1, quantile, 0.5, na.rm = TRUE)
		} else {
			f_ps_future <- f_ps
		}
		this.hasNAs <- apply(is.na(f_ps_future), 2, any)
		hasNAs[this.hasNAs] <- TRUE
		if(write.trajectories) {
			trajectories <- f_ps_future # save only trajectories simulated for the future time
  			save(trajectories, file = file.path(outdir, paste('traj_country', country.obj$code, '.rda', sep='')))
  		}
 		PIs_cqp[country,,] = apply(f_ps_future, 1, quantile, quantiles.to.keep, na.rm = TRUE)
 		mean_sd[country,1,] <- apply(f_ps_future, 1, mean, na.rm = TRUE)
 		mean_sd[country,2,] = apply(f_ps_future, 1, sd, na.rm = TRUE)
 		
	} # end countries loop
	
	##############
	
	mcmc.set <- remove.tfr.traces(mcmc.set)
	bayesTFR.prediction <- structure(list(
				quantiles = PIs_cqp,
				traj.mean.sd = mean_sd,
				nr.traj=nr_simu,
				tfr_matrix_reconstructed = tfr_matrix_reconstructed,
				output.directory=outdir,
				na.index=(1:nr_simu)[hasNAs],
				mcmc.set=load.mcmc.set,
				nr.projections=nr_project,
				burnin=burnin,
				end.year=end.year,
				mu=mu, rho=rho,  sigma_t = sigma_t, sigmaAR1 = sigmaAR1),
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
		
store.traj.ascii <- function(trajectories, n, output.dir, country.code, meta, index, append=FALSE) {
	# Store trajectories into ASCII files of a specific UN format 
	#header <- list(country_code='LocID', period='Period', year='Year', trajectory='Trajectory', tfr='TF')
	header <- get.traj.ascii.header(meta)
	nyears <- dim(trajectories)[1]
	pred.years <- get.prediction.years(meta, nyears)
	pred.period <- get.prediction.periods(meta, nyears)
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
	return(get.prediction.years(pred$mcmc.set$meta, pred$nr.projections+1))
}

get.prediction.years <- function(meta, n) {
	present.year <- as.numeric(rownames(get.data.matrix(meta))[nrow(get.data.matrix(meta))])
	return (seq(present.year, length=n, by=5))
}

get.prediction.periods <- function(meta, n) {
	mid.years <- get.prediction.years(meta, n)
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
							pred$mcmc.set$meta, index=index, append=append)
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
									 output.dir=NULL, revision=14, adjusted=FALSE) {
# Writes three prediction summary files, one in a user-friendly format, one in a UN-format,
# and one parameter file.
	pred <- get.tfr.prediction(sim.dir=dir)
	if (is.null(output.dir)) output.dir <- pred$output.directory
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	tfr.write.projection.summary.and.parameters(pred, output.dir, revision=revision, adjusted=adjusted)
}

tfr.write.projection.summary.and.parameters <- function(pred, output.dir, revision=14, adjusted=FALSE) {
	# two summary files
	do.write.projection.summary(pred, output.dir, revision=revision, adjusted=adjusted)
	# third file about MCMC parameters
	do.write.parameters.summary(pred, output.dir, revision=revision, adjusted=adjusted)
}

do.write.parameters.summary <- function(pred, output.dir, revision=14, adjusted=FALSE) {
	meta <- pred$mcmc.set$meta
	tfr.years <- get.tfr.periods(meta)
	all.years <- c(tfr.years, get.prediction.periods(meta, pred$nr.projections+1)[-1])
	tfr <- get.data.imputed(pred)
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
					par.names=NULL, par.names.cs=tfr.parameter.names.cs(trans=FALSE), 
					thin=1, burnin=0))
		sink(type='message')
		lambda_c <- find.lambda.for.one.country(tfr.and.pred.median, nr.data)
		result <- rbind(result, c(country.obj$name, country.obj$code, 
			if(meta$tau_c[country] > 0) tfr.years[meta$tau_c[country]] else -1, #tau_c
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
	return (list(revision='RevID', variant='VarID', country='LocID', year='TimeID', tfr='TFR'))

"get.UN.variant.names" <- function(pred, ...) UseMethod("get.UN.variant.names")
get.UN.variant.names.bayesTFR.prediction <- function(pred, ...) 
	return(c('BHM median', 'BHM80 lower',  'BHM80 upper', 'BHM95 lower',  'BHM95 upper', 'Low', 
					'High', 'Constant fertility'))
					
"get.friendly.variant.names" <- function(pred, ...) UseMethod("get.friendly.variant.names")
get.friendly.variant.names.bayesTFR.prediction <- function(pred, ...)
	return(c('median', 'lower 80', 'upper 80', 'lower 95', 'upper 95', '-0.5child', '+0.5child', 'constant'))

do.write.projection.summary <- function(pred, output.dir, revision=14, adjusted=FALSE) {
	cat('Creating summary files ...\n')
	data(UN_time)
	data(UN_variants)
	nr.proj <- pred$nr.projections+1
	tfr <- get.data.imputed(pred)
	ltfr <- dim(tfr)[1] - 1
	nr.proj.all <- nr.proj + ltfr
	tfr.years <- get.tfr.periods(pred$mcmc.set$meta)
	pred.period <- get.prediction.periods(pred$mcmc.set$meta, nr.proj)
	header1 <- list(country.name='country_name',  country.code='country_code', variant='variant')
	un.time.idx <- c()
	un.time.label <- as.character(UN_time[,'TLabel'])
	l.un.time.label <- length(un.time.label)
	for (i in 1:ltfr) 
		un.time.idx <- c(un.time.idx, which(un.time.label==tfr.years[i])[1])
	for (i in 1:nr.proj) {
		header1[[paste('year', i, sep='')]] <- pred.period[i]
		un.time.idx <- c(un.time.idx, which(un.time.label==pred.period[i]))
	}
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
								   variant=rep(UN_variants[UN_variants[,'Vshort']==UN.variant.names[ivar],'VarID'], nr.proj.all),
								   country=rep(country.obj$code, nr.proj.all),
								   year=UN_time[un.time.idx,'TimeID'],
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

.bdem.median.shift <- function(type, sim.dir, country, reset=FALSE, shift=0, from=NULL, to=NULL) {
	pred <- do.call(paste('get.', type, '.prediction', sep=''), list(sim.dir))
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
	do.call(paste('store.', class(pred), sep=''), list(pred))
	cat('\nMedian of', country.obj$name, action, 
		if(all.years) 'for all years' else c('for years', pred.years[which.years]), '.\n')
	return(pred)
}

tfr.median.shift <- function(sim.dir, country, reset=FALSE, shift=0, from=NULL, to=NULL) {
	invisible(.bdem.median.shift(type='tfr', sim.dir=sim.dir, country=country, reset=reset, 
				shift=shift, from=from, to=to))
}

.bdem.median.set <- function(type, sim.dir, country, values, years=NULL) {
	pred <- do.call(paste('get.', type, '.prediction', sep=''), list(sim.dir))
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
	do.call(paste('store.', class(pred), sep=''), list(pred))
	cat('\nMedian of', country.obj$name, 'modified for years', pred.years[which.years], '\n')
	return(pred)
}

tfr.median.set <- function(sim.dir, country, values, years=NULL) {
	invisible(.bdem.median.set(type='tfr', sim.dir=sim.dir, country=country, values=values, years=years))
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
	new.pred <- make.tfr.prediction(mcmc.set, end.year=pred$end.year, replace.output=FALSE,
									nr.traj=NULL, burnin=0,
									mu=pred$mu, rho=pred$rho, sigmaAR1=pred$sigmaAR1, 
									countries=countries.idx, adj.factor1=factor1, adj.factor2=factor2,
									forceAR1=forceAR1, save.as.ascii=0, output.dir=NULL,
									write.summary.files=FALSE, is.mcmc.set.thinned=TRUE, 
									write.trajectories=FALSE, verbose=FALSE)
	new.means <- new.pred$traj.mean.sd[,1,2:dim(new.pred$traj.mean.sd)[3]]
	for(icountry in 1:length(countries)) {
		tfr.median.set(sim.dir, countries[icountry], new.means[countries.idx[icountry],])
	}
	# reload adjusted prediction
	invisible(get.tfr.prediction(sim.dir))
}

"get.data.imputed" <- function(pred, ...) UseMethod("get.data.imputed")

get.data.imputed.bayesTFR.prediction <- function(pred, ...)
	return(get.tfr.reconstructed(pred$tfr_matrix_reconstructed, pred$mcmc.set$meta))
	
"get.data.imputed.for.country" <- function(pred, country.index, ...) UseMethod("get.data.imputed.for.country")

get.data.imputed.for.country.bayesTFR.prediction <- function(pred, country.index, ...)
	return(get.observed.with.supplemental(country.index, pred$tfr_matrix_reconstructed, pred$mcmc.set$meta$suppl.data, 'tfr_matrix_all'))
	