
DLcurve <- function(DLpar, tfr, p1, p2, annual = FALSE){ 
 # gives the DL-decrement
 # DLpar is c(Delta1, Delta2, Delta3, Delta4, d_c)
 # tfr is a vector for which the decrements for this curve need to be calculated
 	dlvalue <- rep(0.0, length(tfr))
 	perfct <- if(annual) 1 else 5
	res <- .C("doDLcurve", as.numeric(DLpar), as.numeric(tfr), p1, p2, length(tfr), dl_values=dlvalue, period_multiplier=as.integer(perfct),
	          PACKAGE = "bayesTFR")
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

get.eps.T <- function (DLpar, country, meta, ...) 
{
    tfr <- get.observed.tfr(country, meta, ...)[meta$start_c[country]:meta$lambda_c[country]]
    ldl <- length(tfr)-1
    dl <- DLcurve(DLpar, tfr[1:ldl], meta$dl.p1, meta$dl.p2, annual = meta$annual.simulation)
    eps <- tfr[2:(ldl+1)] - tfr[1:ldl] + dl
    if (!is.null(meta$ar.phase2) && meta$ar.phase2) 
    {
      args <- list(...)
      if ('rho.phase2' %in% names(args) && length(eps) > 1) eps <- c(eps[1], eps[2:ldl]-args[['rho.phase2']] * eps[1:(ldl-1)])
    }

    # Put NAs on eps indexed by meta$indices.outliers[[country]]
    if (as.character(country) %in% names(meta$indices.outliers))
    {
      outlier_indices <- meta$indices.outliers[[as.character(country)]]
      if (!is.null(meta$ar.phase2) && meta$ar.phase2)
        outlier_indices <- sort(unique(c(outlier_indices, outlier_indices+1)))
      
      eps <- eps[-outlier_indices]
    }
      
    return (eps)
}

get_eps_T_all <- function (mcmc, ...) {
  suppl.T <- if(!is.null(mcmc$meta$suppl.data$regions)) mcmc$meta$suppl.data$T_end else 0
	eps_Tc <- matrix(NA, mcmc$meta$T_end-1 + suppl.T, mcmc$meta$nr_countries)
  for (country in mcmc$meta$id_DL){
    theta <- c((mcmc$U_c[country]-mcmc$Triangle_c4[country])*exp(mcmc$gamma_ci[country,])/                                     
      sum(exp(mcmc$gamma_ci[country,])), mcmc$Triangle_c4[country], mcmc$d_c[country])
    idx <- mcmc$meta$start_c[country]:(mcmc$meta$lambda_c[country]-1)
    raw.outliers <- mcmc$meta$indices.outliers[[as.character(country)]]
    if (!is.null(mcmc$meta$ar.phase2) && mcmc$meta$ar.phase2) 
      raw.outliers <- sort(unique(c(raw.outliers, raw.outliers+1)))
    idx <- setdiff(idx, raw.outliers)
    eps_Tc[idx, country] <- get.eps.T(theta, country, mcmc$meta, ...)
  }
	
  return(eps_Tc)
}

find.lambda.for.one.country <- function(tfr, T_end, annual = FALSE) {
    # find the end of Phase II
	if(annual) { # compute 5-year average
	    tfrorig <- tfr
	    Tendorig <- T_end
	    years <- as.integer(names(tfr))
	    ranges <- range(years[years %% 5 == 0])
	    mid.points <- c(0, seq(ranges[1]-2, ranges[2]+3, by = 5))
	    brks <- seq(ranges[1]-5, ranges[2] + 5, by = 5)
	    year.bin <- findInterval(years, brks, left.open = TRUE)
	    tfr <- aggregate(tfr, by = list(year.bin), FUN = mean, na.rm = TRUE)[,"x"]
	    T_end <- year.bin[T_end]
	}
    lambda <- T_end
	if ( sum(tfr<2, na.rm=TRUE)>2 ){
		period <- .get.T.start.end(tfr)[1]+2
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
	if(annual)  # convert lambda from 5-year scale to annual scale
	{
	  lambda <- min(which(year.bin == lambda) + 2, Tendorig)
	  if (length(year.bin) - lambda < 5) { # if in the last time period, set it to the end of the period
	      lambda <- length(year.bin) 
	      while(is.na(tfrorig[lambda])) lambda <- lambda - 1 # move it before the last NA if any
	  }
	}

	return(lambda)
}


.get.T.start.end <- function(tfr) {
	# Return first index after NAs at the beginning of the time series and the last index before any NA's at the end
	isna.tfr <- is.na(tfr)
	start <- if(sum(isna.tfr) > 0 && isna.tfr[1]) which(diff(cumsum(isna.tfr))==0)[1]+1 else 1
	return(c(start, sum(!isna.tfr) + start - 1))
} 

get.observed.tfr <- function(country.index, meta, matrix.name='tfr_matrix', matrix.name.suppl=matrix.name, ...)
	return(get.observed.with.supplemental(country.index, meta[[matrix.name]], meta$suppl.data, matrix.name.suppl))

get.observed.with.supplemental <- function(country.index, matrix, suppl.data, matrix.name='tfr_matrix') {
	data <- matrix[,country.index]
	if(is.null(names(data))) names(data) <- rownames(matrix) # the names attribute is stripped if there is just one row
	if(!is.null(suppl.data[[matrix.name]])) {
    	supp.c.idx <- suppl.data$index.from.all.countries[country.index]
    	if(is.na(supp.c.idx)) {sdata <- rep(NA, nrow(suppl.data[[matrix.name]])); names(sdata) <- rownames(suppl.data[[matrix.name]])}
    	else sdata <- suppl.data[[matrix.name]][,supp.c.idx]
    	data <- c(sdata, data)
    }
	return(data)
}

find.tau.lambda.and.DLcountries <- function(tfr_matrix, min.TFRlevel.for.start.after.1950 = 5, #5.5, 
												max.diff.local.and.global.max.for.start.at.loc = 0.53, #0.5,
												delta.for.local.max = 0.001, annual = FALSE,
												suppl.data=NULL) {
    # gets tau_c and puts NAs before tau_c (start of Phase II)
    # gets ids of DL (where decline has been observed)
    # and divides those into early (tau_c has been not observed) and not early (tau_c observed)
    # find lambda_c's based on definition tfr increased twice, below 2 (start of Phase III)

    post.v2 <- TRUE
	if(!is.null(getOption("TFRphase2.pre.v2", NULL)) && getOption("TFRphase2.pre.v2")==TRUE) {
        # This is for backward-compatibility to make some publications reproducible.
		min.TFRlevel.for.start.after.1950 = 5.5
		max.diff.local.and.global.max.for.start.at.loc = 0.5
		delta.for.local.max = 0
		post.v2 <- FALSE
		warning("Phase II is searched for using pre-v2.0 methodology. To switch to the current method set 'options(TFRphase2.pre.v2=FALSE)'.")
	}
    T_end <- dim(tfr_matrix)[1]
    nr_countries <- dim(tfr_matrix)[2]
    T_end_c <- lambda_c <-rep(T_end, nr_countries)
    T.suppl <- if(is.null(suppl.data$regions)) 0 else dim(suppl.data$tfr_matrix)[1]
    tau_c <- start_c <- rep(NA, nr_countries)
    for (country in 1:nr_countries) {
        data <- get.observed.with.supplemental(country, tfr_matrix, suppl.data)
    	has.suppl <- length(data) > T_end && !is.na(suppl.data$index.from.all.countries[country])
    	# Find the beginning of Phase II (index tau_c)
    	# ignoring NAs at the beginning and at the end
    	T.start.end <- .get.T.start.end(data)
    	T.start <- T.start.end[1]
    	T_end_c[country] <- T.start.end[2]
    	lT <- T_end_c[country] - T.start + 1
    	local_max_indices <- rep(NA, lT)
    	d <- diff(data[T.start:T_end_c[country]])
        does_tfr_decrease <- ifelse(d < delta.for.local.max, 1, 0)
        local_max_indices[1] <- does_tfr_decrease[1]
   		# in the middle only a local max if increase is followed by decrease
        local_max_indices[-c(1, lT)] <- diff(does_tfr_decrease)
   		# at the end local max if preceded by an increase 
        local_max_indices[lT] <- 1 - does_tfr_decrease[lT - 1]
 		value_global_max <- max(data, na.rm = TRUE)
 		max_index <- max(seq(T.start, T_end_c[country]) * (local_max_indices > 0) * 
 						ifelse(data[T.start:T_end_c[country]] >
            				value_global_max - max.diff.local.and.global.max.for.start.at.loc, 1, 0))
        if(post.v2) {
            # move the point to the right if there are more recent points with the same values
        	is.same <- c(0, ifelse(abs(d) < delta.for.local.max, 1, 0), 0)
        	cs.same <- cumsum(is.same[(max_index-T.start+1):(lT+1)])
        	max_index <- max_index + cs.same[min(which(diff(cs.same)==0))]
        }
        tau_c[country] <- max_index
        start_c[country] <- tau_c[country] # start index of phase II

		if(data[tau_c[country]] < min.TFRlevel.for.start.after.1950) {
        	if ((post.v2 && as.integer(names(data)[T.start]) > 1855) || !post.v2) {
        		tau_c[country] <- -1
        		start_c[country] <- which(!is.na(data))[1] # first data point that is not NA
        	}
        }
        if (tau_c[country] > 1 && T.suppl > 0 && has.suppl) 
        	suppl.data$tfr_matrix[1:min(tau_c[country] - 1, T.suppl), suppl.data$index.from.all.countries[country]] <- NA
        if(tau_c[country] > T.suppl + 1)
            tfr_matrix[1:(tau_c[country] - 1 - T.suppl), country] <- NA

        # end index of Phase II
        lambda_c[country] <- do.call(getOption("TFRphase3findfct", "find.lambda.for.one.country"), 
                                     list(data, T_end_c[country], annual = annual))
        
        if (lambda_c[country] < T_end_c[country]) { # set NA all values between lambda_c and T_c_end
         	if(lambda_c[country] < T.suppl) {
         		suppl.data$tfr_matrix[(lambda_c[country] + 1):min(T.suppl, T_end_c[country]),
         										suppl.data$index.from.all.countries[country]] <- NA
         		if(T_end_c[country] > T.suppl) tfr_matrix[1:(T_end_c[country]-T.suppl),country] <- NA
         	} else tfr_matrix[(lambda_c[country] - T.suppl + 1):(T_end_c[country]-T.suppl),country] <- NA
        }
    }

    id_Tistau <- seq(1, nr_countries)[tau_c == T_end_c] # have not started Phase II yet
    id_DL <- seq(1, nr_countries)[(tau_c != T_end_c) & (lambda_c != 1)] # have started Phase II
    # for par in BHM
    id_early <- seq(1, nr_countries)[(tau_c == -1) & (lambda_c != 1)] # start of Phase II not observed (happened prior to observed data)
    # needed for which U_c's to update, this is updated in mcmc
    id_notearly <- seq(1, nr_countries)[(tau_c != -1)  & (tau_c != T_end_c) & (lambda_c != 1)] # start of Phase II observed

    return(list(tau_c = tau_c,  id_Tistau = id_Tistau, id_DL = id_DL, id_early = id_early,
        		id_notearly = id_notearly, tfr_matrix = tfr_matrix, T_end_c=T_end_c, 
        		lambda_c=lambda_c, start_c=start_c, suppl.matrix=suppl.data$tfr_matrix
         ))
}


find.raw.data.outliers <- function(raw.tfr, iso.unbiased, max.drop=1, max.increase=1, 
                                   source.col.name="source"){
  if(is.null(iso.unbiased)) return(NULL)
  MAX_TFR <-  10
  outliers  <- list()
  for (iso.code in iso.unbiased){
    vr.ts <- raw.tfr[(raw.tfr[,source.col.name]=="VR") & (raw.tfr$country_code==iso.code),]
    if (nrow(vr.ts) == 0) next
    year.ts <- round(vr.ts$year - 0.01)
    tfr.ts <- vr.ts$tfr
    annual.tfr.diff <- diff(tfr.ts) / diff(year.ts)
    outlier <- (annual.tfr.diff > max.increase & abs(annual.tfr.diff) < MAX_TFR) |
      (annual.tfr.diff < -max.drop & abs(annual.tfr.diff) < MAX_TFR)
    outlier[is.na(outlier)] <- FALSE
    year.outlier <- year.ts[outlier]
    # for all year in year.outlier, TFR of year->year+1 is an outlier
    if (length(year.outlier) > 0) outliers[[as.character(iso.code)]] <- year.outlier
  }
  return (outliers)
}

mcmc.meta.ini <- function(...,
						U.c.low,
						start.year=1950, present.year=2020, 
						wpp.year=2019, my.tfr.file = NULL, my.locations.file = NULL,
						proposal_cov_gammas = NULL, # should be a list with elements 'values' and 'country_codes'
						verbose=FALSE, uncertainty=FALSE,
						my.tfr.raw.file=NULL, 
						ar.phase2=FALSE, iso.unbiased=NULL, source.col.name = "source",
						use.wpp.data = TRUE) {
	# Initialize meta parameters - those that are common to all chains.
    
  args <- list(...)
	mcmc.input <- list()
	for (arg in names(args)) mcmc.input[[arg]] <- args[[arg]]
	mcmc.input$U.c.low.base <- U.c.low
	mcmc.input$start.year <- start.year
	mcmc.input$present.year <- present.year
	mcmc.input$wpp.year <- wpp.year
	if(present.year-3 > wpp.year) warning("present.year is much larger then wpp.year. Make sure WPP data for present.year are available.")
	if(!use.wpp.data && is.null(my.tfr.file)) {
	    warning("If use.wpp.data is set to FALSE, my.tfr.file should be given. The simulation will use default WPP data.")
	    use.wpp.data <- TRUE
	}
	tfr.with.regions <- set_wpp_regions(start.year=start.year, present.year=present.year, wpp.year=wpp.year, 
										my.tfr.file = my.tfr.file, my.locations.file=my.locations.file, 
										annual = mcmc.input$annual.simulation, ignore.last.observed = uncertainty,
										use.wpp.data = use.wpp.data, verbose=verbose)
	meta <- do.meta.ini(mcmc.input, tfr.with.regions,  
						proposal_cov_gammas=proposal_cov_gammas, verbose=verbose, 
						uncertainty=uncertainty, my.tfr.raw.file=my.tfr.raw.file, ar.phase2=ar.phase2, 
						iso.unbiased=iso.unbiased, source.col.name = source.col.name)
	return(structure(c(mcmc.input, meta), class='bayesTFR.mcmc.meta'))
}
	
	
do.meta.ini <- function(meta, tfr.with.regions, proposal_cov_gammas = NULL, 
						use.average.gammas.cov=FALSE, burnin=200, verbose=FALSE, uncertainty=FALSE, 
						my.tfr.raw.file=NULL, 
						ar.phase2=FALSE, iso.unbiased=NULL, source.col.name = "source") {
  results_tau <- find.tau.lambda.and.DLcountries(tfr.with.regions$tfr_matrix, annual = meta$annual.simulation,
	                                               suppl.data=tfr.with.regions$suppl.data)
	tfr_matrix_all <- tfr.with.regions$tfr_matrix_all
	tfr_matrix_observed <- tfr.with.regions$tfr_matrix
	updated.tfr.matrix <- results_tau$tfr_matrix
	suppl.data <- tfr.with.regions$suppl.data
	if(!is.null(suppl.data$regions)) {
		suppl.data$tfr_matrix_all <- tfr.with.regions$suppl.data$tfr_matrix
		suppl.data$tfr_matrix <- results_tau$suppl.matrix
		suppl.data$T_end <- dim(suppl.data$tfr_matrix)[1]
		suppl.data$nr_countries <- dim(suppl.data$tfr_matrix)[2]
	}
    lambda_c = results_tau$lambda_c
    nr_countries = length(lambda_c)
    nr_countries_estimation <- tfr.with.regions$nr_countries_estimation
                                                   
    # uniform prior for U_c, make lower bound country specific
    if(any(apply(updated.tfr.matrix, 2, function(x) all(is.na(x))))) { # some countries start Phase III before 1950
    	tfr_min_c <- c()
 		# loop over countries to find minimum
 		for (country in 1:nr_countries){
    		data <- get.observed.with.supplemental(country, updated.tfr.matrix, suppl.data)
    		tfr_min_c <- c(tfr_min_c, min(data, na.rm=TRUE))
    	}
    } else
    	tfr_min_c <- apply(updated.tfr.matrix, 2, min, na.rm = TRUE)
    lower_U_c <- ifelse(tfr_min_c > meta$U.c.low.base, tfr_min_c, meta$U.c.low.base)
	prop_cov_gammas <- array(NA, c(nr_countries,3,3))
	if(use.average.gammas.cov) {
			cov.to.average <- get.cov.gammas(sim.dir=meta$output.dir, burnin=burnin)$values
			if (all(is.na(cov.to.average))) {
				warning('Covariance of gamma is NA for all countries. Average from default covariance will be used.', 
						immediate.=TRUE)
			}
			e <- new.env(parent = emptyenv())
			data('proposal_cov_gammas_cii', envir=e)
			cov.to.average <- e$proposal_cov_gammas_cii$values
	} else {
		# get default proposal_cov_gammas_cii and match with country codes of this run
		e <- new.env(parent = emptyenv())
    	data('proposal_cov_gammas_cii', envir=e)
    	current.country.codes <- tfr.with.regions$regions$country_code
    	matched.index <- match(e$proposal_cov_gammas_cii$country_codes, current.country.codes)
    	is.notNA <- !is.na(matched.index)
    	prop_cov_gammas[matched.index[is.notNA],,] <- e$proposal_cov_gammas_cii$values[is.notNA,,]
		if (!is.null(proposal_cov_gammas)) { #user-specified, overwrites defaults for given countries
			matched.index <- match(proposal_cov_gammas$country_codes, current.country.codes)
			is.notNA <- !is.na(matched.index)
			prop_cov_gammas[matched.index[is.notNA],,] <- proposal_cov_gammas$values[is.notNA,,]
		}		
		cov.to.average <- prop_cov_gammas
		if(all(is.na(prop_cov_gammas))) {
		    # if no country match, average everything
		    cov.to.average <- e$proposal_cov_gammas_cii$values
		}
	}
	prop_cov_gammas <- .impute.prop.cov.gamma(prop_cov_gammas, cov.to.average, nr_countries)
  
	output <- list(
	  tfr_matrix=updated.tfr.matrix, 
	  tfr_matrix_all=tfr.with.regions$tfr_matrix_all,
	  tfr_matrix_observed=tfr_matrix_observed,
	  tau_c = results_tau$tau_c, lambda_c = lambda_c,
	  proposal_cov_gammas_cii = prop_cov_gammas,
	  start_c = results_tau$start_c, 
	  id_Tistau = results_tau$id_Tistau, 
	  id_DL = results_tau$id_DL, 
	  id_early = results_tau$id_early,
	  id_notearly = results_tau$id_notearly,
	  id_notearly_estimation = results_tau$id_notearly[results_tau$id_notearly <= nr_countries_estimation],
	  U.c.low=lower_U_c,
	  nr_countries=nr_countries,
	  nr_countries_estimation=nr_countries_estimation,
	  T_end=dim(tfr.with.regions$tfr_matrix)[1], T_end_c=results_tau$T_end_c, 
	  regions=tfr.with.regions$regions,
	  suppl.data=suppl.data
	)
	if(uncertainty)
	{
	  if (is.null(my.tfr.raw.file)) 
	  {
	    ertfr <- new.env(parent = emptyenv())
	    data("rawTFR", envir = ertfr)
	    rawTFR <- ertfr$rawTFR
	  }
	  else
	  {
	    file.type <- substr(my.tfr.raw.file, nchar(my.tfr.raw.file)-2, nchar(my.tfr.raw.file))
	    if (file.type == 'txt')
	      rawTFR <- read.delim(my.tfr.raw.file, sep='\t')
	    else if (file.type == 'csv')
	      rawTFR <- read.delim(my.tfr.raw.file, sep=',')
	    else 
	    {
	      stop("File type not detectible. Please use txt or csv files.")
	    }
	  }
	  rawTFR <- rawTFR[rawTFR$country_code %in% as.numeric(colnames(output$tfr_matrix_all)),]
	  rawTFR <- rawTFR[rawTFR$year < meta$present.year + 0.48,] # to guarantee that 0.5 is included and rounded downward even when 0.01 is added
	  rawTFR <- rawTFR[rawTFR$year > meta$start.year,]
	  output$raw_data.original <- rawTFR
	  output$raw_data.original <- merge(output$raw_data.original, 
	                                    data.frame(country_code=tfr.with.regions$regions$country_code, country_index = 1:nr_countries), 
	                                    by = "country_code")
	  index <- 1:nrow(output$raw_data.original)
	  country.ind.by.year <- list()
	  ind.by.year <- list()
	  
	  # Defining a list containing indices of outliers by countries 
	  # (only countries with outliers will have an entry in the list)
	  # outliers are those where the annual delta is outside of the interval (meta$raw.outliers[1], meta$raw.outliers[2])
	  
	  year.outliers <- find.raw.data.outliers(rawTFR, iso.unbiased, source.col.name = source.col.name,
	                                          max.drop = -meta$raw.outliers[1], max.increase = meta$raw.outliers[2])
	  indices.outliers <- list()
	  yearly.outliers <- list()
	  for (iso.code in names(year.outliers))
	  {
	    country.index <-  get.country.object(country = as.numeric(iso.code), meta = output)$index
	    indices <- year.outliers[[iso.code]] - meta$start.year + 1
	    # We may remove outliers in Phase I?
	    indices <- indices[indices > results_tau$start_c[country.index]]
	    indices <- indices[indices <= results_tau$lambda_c[country.index]]
	    if (length(indices) > 0) {
	      indices.outliers[[as.character(country.index)]] <- indices
	      for (ind in indices) {
	        if (as.character(ind) %in% year.outliers) yearly.outliers[[as.character(ind)]] <- c(yearly.outliers[[as.character(ind)]], country.index)
	        else yearly.outliers[[as.character(ind)]] <- country.index
	      }
	        
	    }
	  }
	  output$indices.outliers <- indices.outliers
	  output$yearly.outliers <- yearly.outliers
	  if (meta$annual.simulation)
	  {
	    for (year in meta$start.year:meta$present.year)
	    {
	      country.ind.by.year[[year-meta$start.year+1]] <- output$raw_data.original$country_index[which(round(output$raw_data.original$year-0.01) == year)]
	      ind.by.year[[year-meta$start.year+1]] <- index[which(round(output$raw_data.original$year-0.01) == year)]
	    }
	    output$raw_data.original$year <- round(output$raw_data.original$year - 0.01)
	    output$country.ind.by.year <- country.ind.by.year
	    output$ind.by.year <- ind.by.year
	  }
	  else
	  {
	    left.distance <- list()
	    years <- as.numeric(rownames(output$tfr_matrix))
	    count <- 1
	    for (year in c(years, rev(years)[1]+5))
	    {
	      inds <- which((output$raw_data.original$year <= year) & (output$raw_data.original$year > year - 5))
	      country.ind.by.year[[count]] <- output$raw_data.original$country_index[inds]
	      ind.by.year[[count]] <- index[inds]
	      left.distance[[count]] <- output$raw_data.original$year[inds] - year + 5
	      count <- count + 1
	    }
	    output$country.ind.by.year <- country.ind.by.year
	    output$ind.by.year <- ind.by.year
	    output$left.distance <- left.distance
	  }
	  output$tfr_all <- output$tfr_matrix_all
	  output$tfr_mu_all <- output$tfr_all
	  output$tfr_sd_all <- matrix(0, nrow = nrow(output$tfr_all), ncol=ncol(output$tfr_all))
	  output$id_phase1_by_year <- list()
	  output$id_phase2_by_year <- list()
	  output$id_phase3_by_year <- list()
	  for (year in 1:output$T_end)
	  {
	    output$id_phase1_by_year[[year]] <- which(output$start_c > year)
	    output$id_phase3_by_year[[year]] <- which(output$lambda_c <= year)
	    output$id_phase2_by_year[[year]] <- setdiff(1:nr_countries, c(output$id_phase1_by_year[[year]], output$id_phase3_by_year[[year]]))
	  }
	}
	
	if (meta$annual.simulation) output$ar.phase2 <- ar.phase2
	
	return(output)

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
					 verbose=FALSE, uncertainty=FALSE, iso.unbiased=NULL,
					 covariates=c('source', 'method'), cont_covariates=NULL, source.col.name = "source") {
				 		 	
	nr_countries <- mcmc.meta$nr_countries

    ############################################
    # U_c, the starting levels of the decline
    U_c <- runif(nr_countries, mcmc.meta$U.c.low, mcmc.meta$U.up)
    for (country in mcmc.meta$id_notearly){
		U_c[country] = get.observed.tfr(country, mcmc.meta)[mcmc.meta$tau_c[country]]
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
    	data <- get.observed.tfr(country, mcmc.meta)
    	minf <- min(data, na.rm = TRUE)
        if (minf < mcmc.meta$Triangle_c4.up) {
        	Triangle_c4[country] = max(mcmc.meta$Triangle_c4.low+0.0001, minf)
          }
    }
   	dontsave.pars <- c('add_to_sd_Tc', 'const_sd_dummie_Tc', 'meta')
   	if(uncertainty) dontsave.pars <- c(dontsave.pars, 'meta3')
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
                        output.dir=paste0('mc', chain.id),
                        traces=0, traces.burnin=0, 
                        rng.state = .Random.seed,
                        compression.type=mcmc.meta$compression.type,
                        dontsave=dontsave.pars, uncertainty=uncertainty
                        ),
                   class='bayesTFR.mcmc')
  ##################################################################
	# distortions
	##################################################################
	# note: the eps will always be NA outside (tau_c, lambda-1)!!
	# ini the epsilons
	if (!is.null(mcmc.meta$ar.phase2) && mcmc.meta$ar.phase2) mcmc$rho.phase2 <- 0.7
	mcmc$eps_Tc <- get_eps_T_all(mcmc, rho.phase2=mcmc$rho.phase2)
	if (uncertainty)
	{
	  # mcmc <- get.obs.estimate.diff(mcmc)
	  # mcmc <- estimate.bias.sd.raw(mcmc)
	  mcmc <- get.obs.estimate.diff.original(mcmc)
	  mcmc <- estimate.bias.sd.original(mcmc, iso.unbiased, covariates, cont_covariates, source.col.name)
	}
	
	# mcmc <- as.list(mcmc)
	return(mcmc)
}

.impute.prop.cov.gamma <- function(prop_cov_gammas, cov.to.average, nr_countries) {
    # where NAs, put averages
    isNA <- apply(is.na(prop_cov_gammas), 1, any)
    if (!any(isNA)) return(prop_cov_gammas)
    avg <- apply(cov.to.average, c(2,3), mean, na.rm=TRUE)
    for(is.na.country in 1:sum(isNA)) 
        prop_cov_gammas[(1:nr_countries)[isNA][is.na.country],,] <- avg
    return(prop_cov_gammas)
}

mcmc.meta.ini.extra <- function(mcmc.set, countries=NULL, my.tfr.file = NULL, 
									my.locations.file=NULL, burnin = 200, verbose=FALSE, uncertainty=FALSE, 
									my.tfr.raw.file=ifelse(uncertainty, file.path(find.package("bayesTFR"), "data", "rawTFR.csv"), NULL),
									average.gammas.cov = TRUE, iso.unbiased=NULL, use.wpp.data = TRUE) {
	update.regions <- function(reg, ereg, id.replace, is.new, is.old) {
		nreg <- list()
		for (name in c('code', 'area_code', 'country_code')) {
			reg[[name]][id.replace] <- ereg[[name]][is.old]
			nreg[[name]] <- c(reg[[name]], ereg[[name]][is.new])
		}
		for (name in c('name', 'area_name', 'country_name')) {
			reg[[name]][id.replace] <- as.character(ereg[[name]])[is.old]
			nreg[[name]] <- c(as.character(reg[[name]]), 
									  as.character(ereg[[name]])[is.new])
		}
		return(nreg)
	}
	meta <- mcmc.set$meta
	#create tfr matrix only for the extra countries
	tfr.with.regions <- set.wpp.extra(meta, countries=countries, annual = meta$annual.simulation,
									  my.tfr.file = my.tfr.file, my.locations.file=my.locations.file, 
									  verbose=verbose, uncertainty=uncertainty, use.wpp.data = use.wpp.data)
	if(is.null(tfr.with.regions)) return(list(meta=meta, index=c()))
	has.mock.suppl <- FALSE
	if(is.null(tfr.with.regions$suppl.data$regions) && !is.null(meta$suppl.data$regions)) {
		# create mock suppl.data in order to get the right data indices
		nrc <- length(tfr.with.regions$regions$code)
		mock.suppl <- list(regions=tfr.with.regions$regions, 
							tfr_matrix=matrix(NA, nrow=nrow(meta$suppl.data$tfr_matrix), ncol=nrc,
												dimnames=list(rownames(meta$suppl.data$tfr_matrix), NULL)),
							index.to.all.countries=1:nrc,
							index.from.all.countries=1:nrc)
		tfr.with.regions$suppl.data <- mock.suppl
		has.mock.suppl <- TRUE
	}
	Emeta <- do.meta.ini(meta, tfr.with.regions=tfr.with.regions, 
								use.average.gammas.cov=average.gammas.cov, burnin=burnin,
						 		verbose=verbose, uncertainty=uncertainty, my.tfr.raw.file = my.tfr.raw.file)
			 		
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
	proposal_cov_gammas_cii <- .impute.prop.cov.gamma(proposal_cov_gammas_cii, proposal_cov_gammas_cii, nr_countries.all)

	new.meta <- list(proposal_cov_gammas_cii = proposal_cov_gammas_cii,
					 nr_countries=nr_countries.all
					)
					
	for (name in c('tfr_matrix', 'tfr_matrix_all', 'tfr_matrix_observed')) {
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
		remove.from.old <- !is.element(idx.old, Emeta[[name]])
		if(any(remove.from.old)) {
			idx <- which(is.element(meta[[name]], id.replace[remove.from.old]))
			if (length(idx) > 0) meta[[name]] <- meta[[name]][-idx]
		}
		new.meta[[name]] <- meta[[name]]
		idx2 <- !is.inold
		if(any(idx2))
			new.meta[[name]] <- c(new.meta[[name]], idx.new.wo.old[Emeta[[name]][idx2]] + meta$nr_countries)
	}
	new.meta[['regions']] <- update.regions(meta$regions, Emeta$regions, id.replace, is.new, is.old)
	if(!is.null(Emeta$suppl.data$regions) && !has.mock.suppl) {
		suppl.id.replace <- meta$suppl.data$index.from.all.countries[id.replace]
		suppl.id.replace <- suppl.id.replace[!is.na(suppl.id.replace)]
		suppl.is.old <- which(is.old)[which(is.element(meta$suppl.data$index.from.all.countries[id.replace], suppl.id.replace))]
		suppl.old <- Emeta$suppl.data$index.from.all.countries[suppl.is.old]
		suppl.is.new <- which(is.new & !is.na(Emeta$suppl.data$index.from.all.countries))
		suppl.new <- Emeta$suppl.data$index.from.all.countries[suppl.is.new]
		for (name in c('tfr_matrix', 'tfr_matrix_all')) {
			meta$suppl.data[[name]][,suppl.id.replace] <- Emeta$suppl.data[[name]][,suppl.old]
			new.meta$suppl.data[[name]] <- cbind(meta$suppl.data[[name]], Emeta$suppl.data[[name]][,suppl.new])
		}
		suppl.is.old.tmp <- rep(FALSE, Emeta$suppl.data$nr_countries)
		suppl.is.old.tmp[suppl.is.old] <- TRUE
		new.meta$suppl.data$regions <- update.regions(meta$suppl.data$regions, Emeta$suppl.data$regions, 
												suppl.id.replace, suppl.new, suppl.old)
		n.new <- ncol(new.meta$suppl.data$tfr_matrix) - ncol(meta$suppl.data$tfr_matrix)
		new.meta$suppl.data$nr_countries <- ncol(new.meta$suppl.data$tfr_matrix)
		new.meta$suppl.data$T_end <- nrow(new.meta$suppl.data$tfr_matrix)
		new.meta$suppl.data$index.from.all.countries <- meta$suppl.data$index.from.all.countries
		new.meta$suppl.data$index.to.all.countries <- meta$suppl.data$index.to.all.countries
		if (n.new > 0) {
			new.meta$suppl.data$index.from.all.countries <- c(new.meta$suppl.data$index.from.all.countries, rep(NA, sum(is.new)))
			new.meta$suppl.data$index.from.all.countries[meta$nr_countries + suppl.is.new] <- seq(meta$suppl.data$nr_countries + 1, 
												length=n.new)
			new.meta$suppl.data$index.to.all.countries <- c(new.meta$suppl.data$index.to.all.countries, 
											seq(meta$nr_countries+1, new.meta$nr_countries)[suppl.is.new])
		} 
	}
	index <- id.replace
	if (new.meta$nr_countries > meta$nr_countries) 
		index <- c(index, seq(meta$nr_countries+1, new.meta$nr_countries))
	for (item in names(new.meta)) {
		meta[[item]] <- new.meta[[item]]
	}
	if (uncertainty)
	{
	  meta$tfr_all <- meta[['tfr_matrix_all']]
	  meta$tfr_mu_all <- meta[['tfr_matrix_all']]
	  meta$tfr_sd_all <- matrix(0, nrow = nrow(meta$tfr_all), ncol=ncol(meta$tfr_all))
	  meta$raw_data.original <- Emeta$raw_data.original
	  meta$country.ind.by.year <- Emeta$country.ind.by.year
	  meta$ind.by.year <- Emeta$ind.by.year
	  if (!meta$annual.simulation) meta$left.distance <- Emeta$left.distance
	  for (i in 1:length(meta$country.ind.by.year))
	  {
	    meta$country.ind.by.year[[i]] <- index[meta$country.ind.by.year[[i]]]
	  }
	  for (year in 1:meta$T_end)
	  {
	      meta$id_phase1_by_year[[year]] <- which(meta$start_c > year)
	      meta$id_phase3_by_year[[year]] <- which(meta$lambda_c <= year)
	      meta$id_phase2_by_year[[year]] <- setdiff(1:meta$nr_countries, c(meta$id_phase1_by_year[[year]], meta$id_phase3_by_year[[year]]))
	  }
	  meta[['id_phase3']] <- which(meta$lambda_c < meta$T_end_c)
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
		U_c[country] = get.observed.tfr(country, mcmc$meta)[mcmc$meta$tau_c[country]]
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

mcmc.meta.ini.subnat <- function(meta, country,
                                 start.year=1950, present.year=2010, annual = NULL,
                                 my.tfr.file = NULL, buffer.size=1000,
                                 verbose=FALSE
                                ) {
  # Initialize meta parameters - those that are common to all chains.
  meta$start.year <- start.year
  meta$present.year <- present.year
  meta$buffer.size <- buffer.size
  if(is.null(annual)) annual <- meta$annual.simulation
  meta$annual.simulation <- annual
  
  tfr.with.regions <- set.wpp.subnat(country=country, start.year=start.year, present.year=present.year,  
                                     my.tfr.file = my.tfr.file, annual = annual, verbose=verbose)
  this.meta <- do.meta.ini.subnat(meta, tfr.with.regions)
  for (item in names(meta))
    if(!(item %in% names(this.meta))) this.meta[[item]] <- meta[[item]]
  return(structure(this.meta, class='bayesTFR.mcmc.meta'))
}

do.meta.ini.subnat <- function(meta, tfr.with.regions) {
    tfr_matrix_all <- tfr.with.regions$tfr_matrix_all
    tfr_matrix_observed <- tfr.with.regions$tfr_matrix
    nr_countries = ncol(tfr_matrix_observed)
    nr_countries_estimation <- tfr.with.regions$nr_countries_estimation

    #tfr_min_c <- apply(tfr_matrix_observed, 2, min, na.rm = TRUE)

    output <- list(
        tfr_matrix=tfr_matrix_observed, 
        tfr_matrix_all=tfr_matrix_all,
        tfr_matrix_observed=tfr_matrix_observed,
        #tau_c = results_tau$tau_c, lambda_c = lambda_c,
        #proposal_cov_gammas_cii = prop_cov_gammas,
        #start_c = results_tau$start_c, 
        #id_Tistau = results_tau$id_Tistau, 
        #id_DL = results_tau$id_DL, 
        #id_early = results_tau$id_early,
        #id_notearly = results_tau$id_notearly,
        #id_notearly_estimation = results_tau$id_notearly[results_tau$id_notearly <= nr_countries_estimation],
        #U.c.low=lower_U_c,
        nr_countries=nr_countries,
        nr_countries_estimation=nr_countries_estimation,
        T_end=dim(tfr.with.regions$tfr_matrix)[1], #T_end_c=results_tau$T_end_c, 
        regions=tfr.with.regions$regions#,
        #suppl.data=suppl.data
    )

    return(output)
    
}

