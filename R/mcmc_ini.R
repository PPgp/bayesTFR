
DLcurve <- function(DLpar, tfr, p1, p2){ 
 # gives the DL-decrement
 # DLpar is c(Delta1, Delta2, Delta3, Delta4, d_c)
 # tfr is a vector for which the decrements for this curve need to  	be calculated
 	dlvalue <- rep(0.0, length(tfr))
	res <- .C("doDLcurve", as.numeric(DLpar), as.numeric(tfr), p1, p2, length(tfr), dl_values=dlvalue)
	return(res$dl_values)
#    t_mid1 <- DLpar[4] + DLpar[3] + DLpar[2] + 0.5 * DLpar[1]
#    t_mid3 <- DLpar[4] + 0.5 * DLpar[3]
#    DLcurve <- 5 * DLpar[5] * (-1/(1 + exp(-log(p1^2)/DLpar[1] * 
#        (tfr - t_mid1))) + 1/(1 + exp(-log(p2^2)/DLpar[3] * (tfr - 
#        t_mid3))))
#    return(ifelse((DLcurve < 0)|(tfr <= 1), 0, DLcurve))

}

##################################################################
# function to get the distortion for country for a given set of DLparameters
# note: this function gives only the eps's in (tau, lambda-1)
# (because rest is NA!!)
get_eps_T = function (DLpar, country, tfr_matrix, start, lambda, p1, p2) 
{
#    eps_T <- NULL
    tfr <- tfr_matrix[start:lambda, country]
    ldl <- length(tfr)-1
    dl <- DLcurve(DLpar, tfr[1:ldl], p1, p2)
#    for (t in start:(lambda - 1)) {
#        eps_T <- c(eps_T, tfr_matrix[t + 1, country] - tfr_matrix[t, 
#            country] + DLcurve(DLpar, tfr_matrix[t, country], 
#            p1, p2))
#    }
    eps_T <- tfr[2:(ldl+1)] - tfr[1:ldl] + dl
    return(eps_T)
}

get_eps_T_all <- function (mcmc) {
	eps_Tc <- matrix(NA, mcmc$meta$T_end-1, mcmc$meta$nr_countries)
    for (country in mcmc$meta$id_DL){
    	theta <- c((mcmc$U_c[country]-mcmc$Triangle_c4[country])*exp(mcmc$gamma_ci[country,])/                                     
                                sum(exp(mcmc$gamma_ci[country,])), mcmc$Triangle_c4[country], mcmc$d_c[country])
        eps_Tc[mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country]-1), country] <- 
                                        get_eps_T(theta, country, mcmc$meta$tfr_matrix, 
                                                          mcmc$meta$start_c[country], 
                                                          mcmc$meta$lambda_c[country], mcmc$meta$dl.p1, mcmc$meta$dl.p2)
        }
    return(eps_Tc)
}

find.lambda.for.one.country <- function(tfr, T_end) {
	lambda <- T_end
	if ( sum(tfr<2, na.rm=TRUE)>2 ){
		period = 3
		while (period<=T_end){
			if (( (tfr[period] - tfr[period-1]) >0 )& 
                    ( (tfr[period-1] - tfr[period-2]) >0) & 
                    (prod(tfr[(period-2):period]<2)==1) 
             	) {
				lambda = period-1
				period = T_end+1
			} else { 
				period = period +1
            }
         }
	}
	return(lambda)
}

get.lambda <- function(tfr_matrix, T_end_c){
	# find lambda_c's based on definition
	# tfr increased twice, below 2
    nr_countries = length(tfr_matrix[1,])
    lambda_c <- T_end_c
    for (country in 1:nr_countries)
    	lambda_c[country] <- find.lambda.for.one.country(tfr_matrix[, country], T_end_c[country])
	return(lambda_c)
}


find.tau.and.DLcountries = function(tfr_matrix, min.TFRlevel.for.start.after.1950 = 5.5, 
												max.diff.local.and.global.max.for.start.at.loc = 0.5) 
# gets tau_c and puts NAs before tau_c
# gets ids of DL (where decline has been observed), excl countries with min <1 (macao and HK)
# and divides those into early and not early
{
    T_end <- length(tfr_matrix[, 1])
    nr_countries <- dim(tfr_matrix)[2]
    T_end_c <- rep(T_end, nr_countries)
    decr_matrix = -tfr_matrix[-1, ] + tfr_matrix[-T_end, ]
    tau_c = rep(NA, nr_countries)
    for (country in 1:nr_countries) {
    	T_end_c[country] = sum(!is.na(tfr_matrix[,country]))
    	local_max_indices = rep(NA, T_end_c[country])
        does_tfr_decrease = ifelse(diff(tfr_matrix[1:T_end_c[country],country]) < 0, 1, 0)
        local_max_indices[1] = does_tfr_decrease[1]
   		# in middle only a local max if increase is followed by decrease
        local_max_indices[-c(1, T_end_c[country])] = diff(does_tfr_decrease)
   		# at end local max if preceded by an increase 
        local_max_indices[T_end_c[country]] = 1 - does_tfr_decrease[T_end_c[country] - 1]
 		value_global_max = max(tfr_matrix[, country], na.rm = TRUE)
 		max_index = max(seq(1, T_end_c[country]) * (local_max_indices >
            0) * (ifelse(tfr_matrix[1:T_end_c[country], country] >
            value_global_max - max.diff.local.and.global.max.for.start.at.loc, 1, 0)))
        tau_c[country] = max_index
        if ((tfr_matrix[tau_c[country], country] < min.TFRlevel.for.start.after.1950)) {
            tau_c[country] = -1
        }
        if (tau_c[country] > 1) {
            tfr_matrix[1:(tau_c[country] - 1), country] = NA
        }
    }

    id_Tistau = seq(1, nr_countries)[tau_c == T_end_c]
        # not needed in fit
    id_DL = seq(1, nr_countries)[(tau_c != T_end_c)
    #& (apply(tfr_matrix, 2, min, na.rm = TRUE) >1) # excludes Macao and Hong Kong
                ]
        # for par in BHM
   id_early = seq(1, nr_countries)[(tau_c == -1)
    #& (apply(tfr_matrix, 2, min, na.rm = TRUE) >1)
    			]
    # needed for which U_c's to update, this is updated in mcmc
   id_notearly = seq(1, nr_countries)[(tau_c != -1)  & (tau_c != T_end_c)
    #& (apply(tfr_matrix, 2, min, na.rm = TRUE) >1)
    			]
    return(list(tau_c = tau_c, id_Tistau = id_Tistau, id_DL = id_DL, id_early = id_early,
        		id_notearly = id_notearly, tfr_matrix = tfr_matrix, T_end_c=T_end_c
         ))
}

mcmc.meta.ini <- function(...,
						U.c.low,
						start.year=1950, present.year=2005, 
						wpp.year=2008, my.tfr.file = NULL, 
						proposal_cov_gammas = NULL, # should be a list with elements 'values' and 'country_codes'
						verbose=FALSE
					 ) {
	# Initialize meta parameters - those that are common to all chains.
	args <- list(...)
	mcmc.input <- list()
	for (arg in names(args)) mcmc.input[[arg]] <- args[[arg]]
	mcmc.input$U.c.low.base <- U.c.low
	mcmc.input$start.year <- start.year
	mcmc.input$present.year <- present.year
	mcmc.input$wpp.year <- wpp.year
	tfr.with.regions <- set_wpp_regions(start.year=start.year, present.year=present.year, 						wpp.year=wpp.year, my.tfr.file = my.tfr.file, verbose=verbose)

	meta <- do.meta.ini(mcmc.input, tfr.with.regions, my.tfr.file=my.tfr.file, 
						proposal_cov_gammas=proposal_cov_gammas, verbose=verbose)
	return(structure(c(mcmc.input, meta), class='bayesTFR.mcmc.meta'))
}
	
	
do.meta.ini <- function(meta, tfr.with.regions, my.tfr.file=NULL, proposal_cov_gammas = NULL, 
						use.average.gammas.cov=FALSE, burnin=200, verbose=FALSE) {
	results_tau <- find.tau.and.DLcountries(tfr.with.regions$tfr_matrix)
	tfr_matrix_all <- tfr.with.regions$tfr_matrix_all
	tfr_matrix_observed <- tfr.with.regions$tfr_matrix
	updated.tfr.matrix <- results_tau$tfr_matrix
    
    lambda_c = get.lambda(tfr_matrix_observed, results_tau$T_end_c)
    nr_countries = length(lambda_c)
    nr_countries_estimation <- tfr.with.regions$nr_countries_estimation
    
    for (country in 1:nr_countries) {
        if (lambda_c[country] < results_tau$T_end_c[country]) {
            updated.tfr.matrix[(lambda_c[country] + 1):results_tau$T_end_c[country],country] = NA
        }
    }
                                               
    # uniform prior for U_c, make lower bound country specific
    tfr_min_c <- apply(updated.tfr.matrix, 2, min, na.rm = TRUE)
    lower_U_c <- ifelse(tfr_min_c > meta$U.c.low.base, tfr_min_c, meta$U.c.low.base)
    # start_c is the max of tau and 1
    start_c = results_tau$tau_c
    start_c[results_tau$id_early] = 1
	
	T_end <- dim(tfr.with.regions$tfr_matrix)[1]
	const_sd_dummie_Tc <- matrix(0, T_end-1, nr_countries)
	const_sd_dummie_Tc[1:5,] <- 1
	
	prop_cov_gammas <- array(NA, c(nr_countries,3,3))
	if(use.average.gammas.cov) {
		cov.to.average <- get.cov.gammas(sim.dir=meta$output.dir, burnin=burnin)$values
		if (all(is.na(cov.to.average))) {
			warning('Covariance of gamma is NA for all countries. Average from default covariance will be used.', 
						immediate.=TRUE)
			data(proposal_cov_gammas_cii)
			cov.to.average <- proposal_cov_gammas_cii$values
		}
	} else {
		# get default proposal_cov_gammas_cii and match with country codes of this run
    	data(proposal_cov_gammas_cii)
    	current.country.codes <- tfr.with.regions$regions$country_code
    	matched.index <- match(proposal_cov_gammas_cii$country_codes, current.country.codes)
    	is.notNA <- !is.na(matched.index)
    	prop_cov_gammas[matched.index[is.notNA],,] <- proposal_cov_gammas_cii$values[is.notNA,,]
    
		if (!is.null(proposal_cov_gammas)) { #user-specified, overwrites defaults for given countries
			matched.index <- match(proposal_cov_gammas$country_codes, current.country.codes)
			is.notNA <- !is.na(matched.index)
			prop_cov_gammas[matched.index[is.notNA],,] <- proposal_cov_gammas$values[is.notNA,,]
		}		
		cov.to.average <- prop_cov_gammas
	}
	# where NAs, put averages
	isNA <- apply(is.na(prop_cov_gammas), 1, any)
	if (any(isNA)) {
		avg <- matrix(NA, 3, 3)
		for(i in 1:3)
			avg[,i] <- apply(cov.to.average[,,i], 2, mean, na.rm=TRUE)
		for(is.na.country in 1:sum(isNA)) 
			prop_cov_gammas[(1:nr_countries)[isNA][is.na.country],,] <- avg
	}

	return(list(
			tfr_matrix=updated.tfr.matrix, 
			tfr_matrix_all=tfr.with.regions$tfr_matrix_all,
			tfr_matrix_observed=tfr_matrix_observed,
            tau_c = results_tau$tau_c, lambda_c = lambda_c,
            const_sd_dummie_Tc=const_sd_dummie_Tc,
            proposal_cov_gammas_cii = prop_cov_gammas,
            start_c = start_c,
            id_Tistau = results_tau$id_Tistau, 
            id_DL = results_tau$id_DL, 
            id_early = results_tau$id_early,
            id_notearly = results_tau$id_notearly,
            id_notearly_estimation = results_tau$id_notearly[results_tau$id_notearly <= nr_countries_estimation],
			U.c.low=lower_U_c,
            nr_countries=nr_countries,
            nr_countries_estimation=nr_countries_estimation,
            T_end=T_end, T_end_c=results_tau$T_end_c,
            regions=tfr.with.regions$regions))

}

# ini MCMC for UN estimates

mcmc.ini <- function(chain.id, mcmc.meta, iter=100,
					 S.ini=5, 
					 a.ini=0, 
					 b.ini=a.ini, 
					 sigma0.ini=0.1, 
					 const.ini=1, 				 
					 gamma.ini=1, Triangle_c4.ini = 1.85,
					 d.ini=0.17,
					 save.all.parameters=FALSE,
					 verbose=FALSE
					 ) {
				 		 	
	nr_countries <- mcmc.meta$nr_countries

    ############################################
    # U_c, the starting levels of the decline
    U_c <- runif(nr_countries, mcmc.meta$U.c.low, mcmc.meta$U.up)
    for (country in mcmc.meta$id_notearly){
		U_c[country] = mcmc.meta$tfr_matrix[mcmc.meta$tau_c[country],country ]
	}
	
	##############################################
	# NON-CONST SD
    ##############################################
    # for non-constant sd: sd decreases linearly from TFR of f_sd (S in report) to both sides
    ### notation in report:
    ### a = a_sd, b = b_sd, S = f_sd, sigma0 = sigma0, c = const_sd
    S_sd <- S.ini # max SD
    a_sd <- a.ini
    b_sd <- b.ini
    sigma0 <- sigma0.ini
    const_sd <- const.ini

    sd_eps_tau <- mcmc.meta$sd.eps.tau0
    mean_eps_tau <- mcmc.meta$mean.eps.tau0


    ######################################################
    # for tranformed d: dt ~ N(chi, psi^2)
    ######################################################
    # initiate psi and chi
    chi <- mcmc.meta$chi0
    psi <- mcmc.meta$psi0

	# initiate d:
	d_c= rep(d.ini, nr_countries)

	##################################################################
	# alpha and deltas
    alpha <- mcmc.meta$alpha0.p
    delta <- rep(mcmc.meta$delta0, 3)

    ##################################################################
    # gammas
    gamma_ci <- matrix(gamma.ini, nrow=nr_countries, ncol=3)
                
	### for Triangle4's:
    Triangle4 <- mcmc.meta$Triangle4.0
    delta4 <- mcmc.meta$delta4.0
    T_end = mcmc.meta$T_end
   	Triangle_c4 <- rep(Triangle_c4.ini, nr_countries)
    for (country in mcmc.meta$id_DL) {
    	minf <- min(mcmc.meta$tfr_matrix[, country], na.rm = TRUE)
        if (minf < mcmc.meta$Triangle_c4.up) {
        	Triangle_c4[country] = max(mcmc.meta$Triangle_c4.low+0.0001, minf)
          }
    }
   	dontsave.pars <- c('add_to_sd_Tc', 'meta')
    if (!save.all.parameters) dontsave.pars <- c(dontsave.pars, 'eps_Tc')
    if (!exists(".Random.seed")) runif(1)	    	
	mcmc <- structure(list(
						meta=mcmc.meta,
                        U_c=U_c, d_c=d_c, gamma_ci=gamma_ci, 
                        Triangle_c4 = Triangle_c4,
                        delta4 = delta4, Triangle4 = Triangle4,
                        alpha=alpha, delta=delta, psi=psi, chi=chi, 
                        a_sd=a_sd, b_sd=b_sd, const_sd=const_sd,
                        S_sd=S_sd, sigma0=sigma0,sd_eps_tau = sd_eps_tau, 
                        mean_eps_tau = mean_eps_tau,
                        d.ini=d.ini, gamma.ini=gamma.ini,Triangle_c4.ini=Triangle_c4.ini,
                        iter=iter, finished.iter=1, length = 1,
                        id=chain.id,
                        output.dir=paste('mc', chain.id, sep=''),
                        traces=0, traces.burnin=0, 
                        rng.state = .Random.seed,
                        dontsave=dontsave.pars
                        ),
                   class='bayesTFR.mcmc')
                   
	##################################################################
	# distortions
	##################################################################
	# note: the eps will always be NA outside (tau_c, lambda-1)!!
	# ini the epsilons
	mcmc$eps_Tc <- get_eps_T_all(mcmc)
	return(mcmc)
}

mcmc.meta.ini.extra <- function(mcmc.set, countries=NULL, my.tfr.file = NULL, 
									burnin = 200, verbose=FALSE) {
	meta <- mcmc.set$meta
	#create tfr matrix only for the extra countries
	tfr.with.regions <- set.wpp.extra(meta, countries=countries, 
									  my.tfr.file = my.tfr.file, verbose=verbose)
	if(is.null(tfr.with.regions)) return(list(meta=meta, index=c()))
	Emeta <- do.meta.ini(meta, tfr.with.regions=tfr.with.regions, 
								my.tfr.file=my.tfr.file, 
								use.average.gammas.cov=TRUE, burnin=burnin,
						 		verbose=verbose)
						 		
	# join the new meta with the existing one
	is.old <- tfr.with.regions$is_processed
	is.new <- !tfr.with.regions$is_processed
	nold <- sum(is.old)
	nr_countries.all <- meta$nr_countries + Emeta$nr_countries - nold
	if (nold > 0) {
		codes.replace <- tfr.with.regions$regions$country_code[is.old]
		id.replace <- unlist(sapply(codes.replace, get.country.object, meta=meta)['index',])
	} else {id.replace <- c()}
	proposal_cov_gammas_cii <- array(NA, c(nr_countries.all, 3, 3))
	for (i in 1:3) {
		meta$proposal_cov_gammas_cii[id.replace,i,] <- matrix(
								Emeta$proposal_cov_gammas_cii[is.old,i,], ncol=3)
		proposal_cov_gammas_cii[,i,] <- rbind(meta$proposal_cov_gammas_cii[,i,], 
							matrix(Emeta$proposal_cov_gammas_cii[is.new,i,], ncol=3))
	}
	new.meta <- list(proposal_cov_gammas_cii = proposal_cov_gammas_cii,
					 nr_countries=nr_countries.all
					)
					
	for (name in c('tfr_matrix', 'tfr_matrix_all', 'tfr_matrix_observed', 'const_sd_dummie_Tc')) {
		meta[[name]][,id.replace] <- Emeta[[name]][,is.old]
		new.meta[[name]] <- cbind(meta[[name]], Emeta[[name]][,is.new])
	}
	for (name in c('tau_c', 'lambda_c', 'start_c', 'U.c.low', 'T_end_c')) {
		meta[[name]][id.replace] <- Emeta[[name]][is.old]
		new.meta[[name]] <- c(meta[[name]], Emeta[[name]][is.new])
	}
	idx.old <- (1:length(is.old))[is.old]
	idx.new.wo.old <- cumsum(is.new)
	for (name in c('id_Tistau', 'id_DL', 'id_early', 'id_notearly')) {
		is.inold <- is.element(Emeta[[name]], idx.old)
		if(any(is.inold)) {
			codes <- unlist(sapply(Emeta[[name]][is.inold], 
								get.country.object, meta=Emeta, index=TRUE)['code',])
			idx <- unlist(sapply(codes, get.country.object, meta=meta)['index',])
			not.incl <- (1:length(idx))[!is.element(meta$regions$country_code[idx], codes)]
			if (length(not.incl) > 0)
				meta[[name]] <- meta[[name]][-not.incl]
		}
		new.meta[[name]] <- meta[[name]]
		idx2 <- !is.inold
		if(any(idx2))
			new.meta[[name]] <- c(new.meta[[name]], idx.new.wo.old[Emeta[[name]][idx2]] + meta$nr_countries)
	}
	new.meta[['regions']] <- list()
	for (name in c('code', 'area_code', 'country_code')) {
		meta$regions[[name]][id.replace] <- Emeta$regions[[name]][is.old]
		new.meta$regions[[name]] <- c(meta$regions[[name]], Emeta$regions[[name]][is.new])
	}
	for (name in c('name', 'area_name', 'country_name')) {
		meta$regions[[name]][id.replace] <- as.character(Emeta$regions[[name]])[is.old]
		new.meta$regions[[name]] <- c(as.character(meta$regions[[name]]), 
									  as.character(Emeta$regions[[name]])[is.new])
	}
	index <- id.replace
	if (new.meta$nr_countries > meta$nr_countries) 
		index <- c(index, seq(meta$nr_countries+1, new.meta$nr_countries))
	for (item in names(new.meta)) {
		meta[[item]] <- new.meta[[item]]
	}

	return(list(meta=meta, index=index, index.replace=id.replace, 
				index_DL=index[is.element(index, new.meta$id_DL)]))
}

mcmc.ini.extra <- function(mcmc, countries, index.replace=NULL) {
	nr.countries.extra <- length(countries)
	nreplace <- length(index.replace)
	Uc <- runif(nr.countries.extra, mcmc$meta$U.c.low, mcmc$meta$U.up)
	if(nreplace > 0) {
		mcmc$U_c[index.replace] <- Uc[1:nreplace]
		mcmc$d_c[index.replace] <- mcmc$d.ini
		mcmc$Triangle_c4[index.replace] <- mcmc$Triangle_c4.ini
		mcmc$gamma_ci[index.replace,] <- matrix(mcmc$gamma.ini, nrow=nreplace, ncol=3)
	}
	U_c <- mcmc$U_c
	if(nr.countries.extra > nreplace)
		U_c <- c(U_c, Uc[(nreplace+1):nr.countries.extra])
    for (country in mcmc$meta$id_notearly[is.element(mcmc$meta$id_notearly, countries)]){
		U_c[country] = mcmc$meta$tfr_matrix[mcmc$meta$tau_c[country],country ]
	}
	mcmc.update <- list(U_c=U_c, d_c=c(mcmc$d_c, rep(mcmc$d.ini, nr.countries.extra-nreplace)), 
						gamma_ci=rbind(mcmc$gamma_ci, 
									   matrix(mcmc$gamma.ini, nrow=nr.countries.extra-nreplace, ncol=3)), 
                        Triangle_c4 = c(mcmc$Triangle_c4, rep(mcmc$Triangle_c4.ini, nr.countries.extra-nreplace))
                        )
	for (item in names(mcmc.update)) {
		mcmc[[item]] <- mcmc.update[[item]]
	}
	return(mcmc)
}