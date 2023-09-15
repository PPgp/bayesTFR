# Read in the UN estimates

set_wpp_regions <- function(start.year=1950, present.year=2010, wpp.year=2012, my.tfr.file=NULL, 
							my.locations.file=NULL, annual = FALSE, ignore.last.observed = FALSE, 
							use.wpp.data = TRUE, verbose=FALSE) {
# outputs:
# tfr_matrix_all, with each column one countries UN estimates
# tfr_matrix with NAs after last observed data point
# reg_name, reg_code and area_name, area_code: vectors with UN regions and continents

	########################################
	# set data and match with areas
	########################################
	un.object <- read.UNtfr(wpp.year=wpp.year, my.tfr.file=my.tfr.file, 
								present.year=present.year, annual = annual, 
								use.wpp.data = use.wpp.data, verbose=verbose)
	tfr_data <- un.object$data.object$data
	# not just countries, includes areas etc as well
	# get region and area data
	locations <- read.UNlocations(tfr_data, wpp.year=wpp.year, my.locations.file=my.locations.file, verbose=verbose)
	loc_data <- locations$loc_data
	include <- locations$include
	prediction.only <- locations$prediction.only

	tfr_data_countries <- tfr_data[include,]
	nr_countries_estimation <- length(tfr_data_countries[,1])
	if(any(!is.na(prediction.only))) { # move prediction countries at the end of tfr_data
		tfr_data_prediction <- tfr_data[prediction.only,]
		tfr_data_countries <- rbind(tfr_data_countries, tfr_data_prediction)
	}
	
	TFRmatrix.regions <- get.TFRmatrix.and.regions(tfr_data_countries, loc_data, 
												start.year=start.year, 
												present.year=present.year,
												annual = annual, verbose=verbose, 
												interpolate = wpp.year < 2022 && annual && is.null(my.tfr.file),
												ignore.last.observed = ignore.last.observed)
	if(!annual) {
	    TFRmatrixsuppl.regions <- .get.suppl.matrix.and.regions(un.object, TFRmatrix.regions, loc_data, 
									start.year, present.year)
	    if(!is.null(un.object$suppl.data.object) && verbose) 
		    cat('Dimension of the supplemental TFR matrix:', dim(TFRmatrixsuppl.regions$obs_matrix), '\n')
	} else TFRmatrixsuppl.regions <- NULL		
	return(list(tfr_matrix=TFRmatrix.regions$tfr_matrix, 
				tfr_matrix_all=TFRmatrix.regions$tfr_matrix_all, 
				regions=TFRmatrix.regions$regions, 
				nr_countries_estimation=nr_countries_estimation,
				suppl.data=.get.suppl.data.list(TFRmatrixsuppl.regions)
				))
}

.get.suppl.data.list <- function(matrix.regions, matrix.name='tfr_matrix'){
	suppl.data <- list(regions=NULL, index.to.all.countries=NULL, index.from.all.countries=NULL)
	suppl.data[[matrix.name]] <- NULL 
	if(!is.null(matrix.regions)) {
		suppl.data[[matrix.name]] <- matrix.regions$obs_matrix
		suppl.data$regions <- matrix.regions$regions
		suppl.data$index.to.all.countries <- matrix.regions$all.countries.index
		suppl.data$index.from.all.countries <- matrix.regions$index.from.all.countries
	}
	return(suppl.data)
}

load.from.wpp <- function(dataset, wpp.year, annual = FALSE, ...){
    if(wpp.year >= 2022) {
        dsname <- paste0(dataset, if(annual) 1 else 5)
        ds <- load.bdem.dataset(dsname, wpp.year, ..., check.if.exists = TRUE)
        if(!is.null(ds)) return(ds)
        # otherwise check the original name
    }
    load.bdem.dataset(dataset, wpp.year, ...)
}

load.bdem.dataset <- function(dataset, wpp.year, envir=NULL, verbose=FALSE, check.if.exists = FALSE) {
	pkg <- paste0('wpp', wpp.year)
	do.call('require', list(pkg))
	if(check.if.exists && !dataset %in% data(package = pkg)$results[, "Item"]) return(NULL)
	if(verbose) cat('Loading ', dataset, ' from ', pkg, '.\n')
	if(is.null(envir)) envir <- new.env()
	do.call('data', list(dataset, package=pkg, envir=envir))
	return(envir[[dataset]])
}

read.tfr.file <- function(file) return(read.delim(file=file, comment.char='#', check.names=FALSE))

do.read.un.file <- function(un.file.name, wpp.year, my.file=NULL, present.year=2012, 
                            annual = FALSE, countries = NULL, use.wpp.data = TRUE, verbose=FALSE) {
	tfr_data <- load.from.wpp(un.file.name, wpp.year, annual = annual, verbose=verbose)
	my.tfr.file <- my.file
	if(!is.null(tfr_data) || !is.null(my.tfr.file)) {
		if(!is.element('last.observed', colnames(tfr_data)))
			tfr_data <- cbind(tfr_data, last.observed=present.year)
		if(!is.element('include_code', colnames(tfr_data)))
			tfr_data <- cbind(tfr_data, include_code=rep(-1, nrow(tfr_data)))
		colnames(tfr_data)[colnames(tfr_data)=='name'] <- 'country'
	}
	replaced <- c()
	added <- c()
	# read user-defined TFRs and replace the UN data with it
	if(!is.null(my.tfr.file)) {
		if(verbose) cat('Reading file ', my.tfr.file, '.\n')
		my.tfr_data <- read.tfr.file(file=my.tfr.file)
		colnames(my.tfr_data)[colnames(my.tfr_data)=='name'] <- 'country'
		if(use.wpp.data)
		    cols.to.use <- colnames(my.tfr_data)[is.element(colnames(my.tfr_data), colnames(tfr_data))]
		else {
		    # if default wpp data is not used, trust the user that the columns are legitimate
		    # and completely overwrite the default dataset
		    cols.to.use <- colnames(my.tfr_data)
		    tfr_data <- my.tfr_data
		    if(!is.element('last.observed', colnames(tfr_data)))
		        tfr_data <- cbind(tfr_data, last.observed=present.year)
		    if(!is.element('include_code', colnames(tfr_data)))
		        tfr_data <- cbind(tfr_data, include_code=rep(-1, nrow(tfr_data)))
		}
		# don't overwrite country_name
		cols.wo.name <- setdiff(cols.to.use, 'country')
		if (!is.element('country_code', cols.to.use))
			stop('Column "country_code" must be included in ', my.tfr.file)
		codes <- my.tfr_data$country_code
		if (length(codes) > 0) {
			for (icode in 1:length(codes)) {
			    if(!is.null(countries) && !codes[icode] %in% countries) next
				idx <- tfr_data$country_code == codes[icode]
				if(!any(idx)){ # country not found, add a row
					row <- tfr_data[1,]
					row[] <- NA
					if(is.element('country', cols.to.use)) row['country'] <- my.tfr_data[icode,'country']
					tfr_data <- rbind(tfr_data, row)
					idx <- nrow(tfr_data)
					# set defaults
					tfr_data[idx,'last.observed'] <- present.year
					tfr_data[idx,'include_code'] <- -1
					added <- c(added, codes[icode])
				} else replaced <- c(replaced, codes[icode])
				# overwrite the UN data
				colsnotna <- which(!is.na(my.tfr_data[icode,cols.wo.name]))
				tfr_data[idx,cols.wo.name[colsnotna]] <- my.tfr_data[icode,cols.wo.name[colsnotna]]
			}
		}
		if(length(replaced) > 0) 
			if(verbose) cat('\tValues in', paste(cols.to.use, collapse=', '), '(possibly partially) replaced for countries', 
						paste(replaced, collapse=', '), '.\n') 
		if(length(added) > 0) 
		    if(verbose) cat('\tCountries', paste(added, collapse=', '), 'added to the UN data.\n')
		if(length(replaced) == 0 && length(added) == 0)
		    cat('No country in', my.tfr.file, 'matched the UN file.\n')
	} 
	return(list(data=tfr_data, replaced=replaced, added=added))
}

read.UNtfr <- function(wpp.year, my.tfr.file=NULL, annual = FALSE, ...) {
	data <- do.read.un.file('tfr', wpp.year, my.file=my.tfr.file, annual = annual, ...)
	suppl.data <- NULL
	if(!annual) {
	    suppl.data <- do.read.un.file('tfr_supplemental', wpp.year, my.file = my.tfr.file, annual = annual, ...)
	    if(is.null(suppl.data$data)) suppl.data <- NULL
	}
	return(list(data.object=data, suppl.data.object=suppl.data))
}

read.UNlocations <- function(data, wpp.year, package="bayesTFR", my.locations.file=NULL, verbose=FALSE) {
	loc_data <- load.bdem.dataset('UNlocations', wpp.year, verbose=verbose)
	if(!is.null(my.locations.file)) {
		my.locations <- read.tfr.file(file=my.locations.file)
		#if(!all(colnames(loc_data) %in% colnames(my.locations)))
		required.columns <- c("country_code", "name", "location_type")
		if(!all(required.columns %in% colnames(my.locations)))
			stop('my.locations.file must contain columns: ', paste(required.columns, collapse=', '))
		col.not.present <- colnames(loc_data)[!colnames(loc_data)%in%colnames(my.locations)]
		if(length(col.not.present) > 0) {
			my.locations <- cbind(my.locations, as.data.frame(matrix(NA, nrow=nrow(my.locations), 
							ncol=length(col.not.present), dimnames=list(NULL, col.not.present))))
		}
		my.locations <- my.locations[,colnames(loc_data)]
		overlap <- my.locations$country_code %in% loc_data$country_code
		if(any(overlap)) { # overwrite UN locations
			woverlap <- which(overlap)
			loc_data[sapply(my.locations$country_code[woverlap], function(x) which(loc_data$country_code == x)),] <- my.locations[woverlap,]
		}
		wnoverlap <- which(!overlap)
		if(length(wnoverlap)>0) { # add to UN locations
			loc_data <- rbind(loc_data, my.locations[wnoverlap,])
		}
	}
	
	icd <- data(list = paste0('include_', wpp.year), envir = environment(), package = package)
	include_codes <- get(icd)
    loc_data <- merge(loc_data, include_codes, by='country_code', all.x=TRUE)
	# this include some areas that are not in the tfr file
	loc_data[is.na(loc_data$include_code),"include_code"] <- 0
	# get the include_code from this file to get the countries from the tfr file
	nr_tfr_outcomes <- length(data[,1])
	include <- rep(NA, nr_tfr_outcomes)
	prediction.only <- rep(NA, nr_tfr_outcomes)
	n_loc_data <- length(loc_data[,1])

	for (i in 1:nr_tfr_outcomes){
 		loc_index <- (1:n_loc_data)[loc_data$country_code == data$country_code[i]]
 		if (length(loc_index)<=0)
 			stop('Country ', data$country_code[i], ' not found in the location dataset of wpp', wpp.year, '.')
 		if (length(loc_index) > 1) {
 		    warning("Duplicates found in the location file for country ", data$country_code[i])
 		    loc_index <- loc_index[1]
 		}
 		incl.code <- if(data$include_code[i] >= 0) data$include_code[i] else loc_data$include_code[loc_index]
 		include[i] <- incl.code == 2
 		prediction.only[i] <- incl.code == 1	
	}
	return(list(loc_data=loc_data, include=include, prediction.only=prediction.only#, include.code=incl.code)
			))
}


get.observed.time.matrix.and.regions <- function(data, loc_data, start.year=1950, present.year=2010, annual = FALSE,
												 datacolnames=c(country.code='country_code', country.name='country', reg.name='reg_name',
															reg.code='reg_code', area.name='area_name', area.code='area_code'), 
												 interpolate = FALSE, ignore.last.observed = FALSE) {
	tfr_data <- data
	nr_countries <- length(tfr_data[,1])
	if (annual && interpolate) # interpolate 5-year data
	{
	    names.tfr.data <- names(tfr_data)
	    num.columns <- grep('^[0-9]{4}.[0-9]{4}$', names.tfr.data)
	    ncol.tfr <- length(num.columns)
	    cols.starty <- as.integer(substr(names.tfr.data[num.columns], 1,4))
	    cols.endy <- as.integer(substr(names.tfr.data[num.columns], 6,9))
	    years.to.interp.to <- (cols.starty[1]):cols.endy[length(cols.endy)]
	    years.to.interp.from <- seq(cols.starty[1]+3, cols.endy[length(cols.endy)]-2, by = 5)
	    tfr_data_new <- matrix(NA, nrow = nrow(tfr_data), ncol = length(years.to.interp.to))
	    for(row in 1:nrow(tfr_data_new)) 
	        tfr_data_new[row, ] <- approx(years.to.interp.from, tfr_data[row, num.columns], 
	                                      xout = years.to.interp.to, rule=2)$y
	    tfr_data_new <- cbind(tfr_data[, -num.columns], tfr_data_new) 
	    colnames(tfr_data_new) <- c(colnames(tfr_data)[-num.columns], years.to.interp.to)
	  tfr_data <- tfr_data_new
	}
	names.tfr.data <- names(tfr_data)
	if(!annual) { # index of year-columns 
	    num.columns <- grep('^[0-9]{4}.[0-9]{4}$', names.tfr.data) 
	} else 
	    num.columns <- grep('^[0-9]{4}$', names.tfr.data)
	ncol.tfr <- length(num.columns)
	cols.starty <- as.integer(substr(names.tfr.data[num.columns], 1,4))
	if(!annual) {
	    cols.endy <- as.integer(substr(names.tfr.data[num.columns], 6,9))
	} else
	    cols.endy <- cols.starty+0.5
	    #start.index <- (1:ncol.tfr)[cols.starty == start.year]
	start.index <- (1:ncol.tfr)[(cols.starty <= start.year) & (cols.endy > start.year)]
	if(length(start.index) <= 0) {
		if(cols.starty[1] > start.year)	start.index <- 1
		else if(cols.endy[length(cols.endy)] <= start.year) return(NULL)
	}
	start.col <- names.tfr.data[num.columns][start.index[1]]
	present.index <- (1:ncol.tfr)[(cols.endy >= present.year) & (cols.starty <= present.year)]
	if(length(present.index) <= 0) {
		if(cols.endy[length(cols.endy)] < present.year) present.index <- length(cols.endy)
		else if(cols.starty[1] > present.year) return(NULL)
	}
	present.col <- names.tfr.data[num.columns][present.index[1]]

	tfr_matrix <- t(tfr_data[,which.max(names(tfr_data)==start.col):which.max(names(tfr_data)==present.col), drop = FALSE])
	start.years <- as.integer(substr(rownames(tfr_matrix), 1,4))
	if (!annual) {
	    end.years <- as.integer(substr(rownames(tfr_matrix), 6,9))
	    mid.years <- start.years + ceiling((end.years- start.years)/2)
	} else mid.years <- start.years
	
	rownames(tfr_matrix) <- mid.years				
	tfr_matrix_all <- tfr_matrix

	reg_name <- reg_code <- area_name <- area_code <- NULL
	has.reg.name <- is.element(datacolnames['reg.name'], colnames(loc_data))
	has.reg.code <- is.element(datacolnames['reg.code'], colnames(loc_data))
	has.area.name <- is.element(datacolnames['area.name'], colnames(loc_data))
	has.area.code <- is.element(datacolnames['area.code'], colnames(loc_data))
	all.na <- na.array <- rep(NA, nr_countries)
	for (i in 1:nr_countries){
		loc_index <- which.max(loc_data[[datacolnames['country.code']]] == tfr_data[[datacolnames['country.code']]][i])
 		if(has.reg.name) reg_name <- c(reg_name, paste(loc_data$reg_name[loc_index]))
 		if(has.reg.code) reg_code <- c(reg_code, loc_data$reg_code[loc_index])
 		if(has.area.name) area_name <- c(area_name, paste(loc_data$area_name[loc_index]))
 		if(has.area.code) area_code <- c(area_code, loc_data$area_code[loc_index])
		# set NAs for entries that are not observed data (after last.observed) 
		if(!ignore.last.observed) tfr_matrix[(tfr_data[i,'last.observed'] < mid.years), i] <- NA
		all.na[i] <- all(is.na(tfr_matrix[,i]))
	}
	country.codes <- tfr_data[[datacolnames['country.code']]]
	colnames(tfr_matrix) <- colnames(tfr_matrix_all) <- country.codes
	regions <- list(name=if(has.reg.name) reg_name else na.array, 
					code=if(has.reg.code) reg_code else na.array, 
					area_name=if(has.area.name) area_name else na.array, 
					area_code=if(has.area.code) area_code else na.array,
				country_name=tfr_data[[datacolnames['country.name']]], 
				country_code=country.codes)
	return(list(obs_matrix=tfr_matrix, obs_matrix_all=tfr_matrix_all, regions=regions, all.na=all.na))
}

get.TFRmatrix.and.regions <- function(tfr_data, ..., verbose=FALSE){
	result <- get.observed.time.matrix.and.regions(data=tfr_data, ...)
	if (verbose) 
		cat('Dimension of the TFR matrix:', dim(result$obs_matrix), '\n')
	return(list(tfr_matrix=result$obs_matrix, tfr_matrix_all=result$obs_matrix_all, regions=result$regions))
}


.get.suppl.matrix.and.regions <- function(un.object, matrix.regions, loc_data, start.year, present.year, ...) {
	matrixsuppl.regions <- NULL
	if(is.null(un.object$suppl.data.object)) return(NULL)
	suppl.data <- un.object$suppl.data.object$data
	include <- which(is.element(suppl.data[,'country_code'], matrix.regions$regions$country_code))
	if(length(include) > 0)
		matrixsuppl.regions <- get.observed.time.matrix.and.regions(
							suppl.data[include,], loc_data, 
							start.year=start.year, 
							present.year=present.year, ...)
	if(!is.null(matrixsuppl.regions)) {
		# remove rows where all tfr entries are NA
		matrixsuppl.regions$obs_matrix <- matrixsuppl.regions$obs_matrix[,!matrixsuppl.regions$all.na, drop=FALSE]
		matrixsuppl.regions$obs_matrix_all <- matrixsuppl.regions$obs_matrix_all[,!matrixsuppl.regions$all.na, drop=FALSE]
		for(name in c(names(matrixsuppl.regions$regions))) 
			matrixsuppl.regions$regions[[name]] <-  matrixsuppl.regions$regions[[name]][!matrixsuppl.regions$all.na]
		if(ncol(matrixsuppl.regions$obs_matrix) > 0) {
			matrixsuppl.regions$all.countries.index <- c()
			index.from.all <- rep(NA, length(matrix.regions$regions$country_code))
			for(i in 1:length(matrixsuppl.regions$regions$country_code)) {
				incl.idx <- which(is.element(matrix.regions$regions$country_code,
												matrixsuppl.regions$regions$country_code[i]))
				matrixsuppl.regions$all.countries.index <- c(matrixsuppl.regions$all.countries.index, incl.idx)
				index.from.all[incl.idx] <- i
			}
			matrixsuppl.regions$index.from.all.countries <- index.from.all
		} else {
			matrixsuppl.regions$regions <- NULL
			matrixsuppl.regions$obs_matrix <- NULL
			matrixsuppl.regions$obs_matrix_all <- NULL
		}
	}			
	return(matrixsuppl.regions)
}


.extra.matrix.regions <- function(data, countries, meta, package="bayesTFR", my.locations.file=NULL, 
                                  verbose=FALSE, uncertainty=FALSE, ...) {
    tfrs <- data
	country.codes.processed <- meta$regions$country_code
	ncountries <- get.nr.countries(meta)
	ncountries.est <- get.nrest.countries(meta)
	replaced.processed <- intersect(tfrs$replaced, country.codes.processed)
	idx.replaced.processed <- is.element(country.codes.processed, replaced.processed)
	replaced.processed.est <- country.codes.processed[1:ncountries.est][idx.replaced.processed[1:ncountries.est]]
	
	added.processed <- intersect(tfrs$added, country.codes.processed)
	idx.added.processed <- is.element(country.codes.processed, added.processed)
	added.processed.est <- country.codes.processed[1:ncountries.est][idx.added.processed[1:ncountries.est]]
	
	countries.processed <- intersect(countries, country.codes.processed)
	idx.countries.processed <- is.element(country.codes.processed, countries.processed)
	countries.processed.est <- country.codes.processed[1:ncountries.est][idx.countries.processed[1:ncountries.est]]

	if((length(replaced.processed.est)+length(added.processed.est)+length(countries.processed.est)> 0) && !(uncertainty)) {
		cat('\nCountries', paste(c(replaced.processed.est, added.processed.est, 
									countries.processed.est), collapse=', '), 
					'used for estimation and will be excluded from processing.\n')
	}
	
	include.codes <- c(setdiff(tfrs$replaced, replaced.processed.est), 
					   setdiff(tfrs$added, added.processed.est))
	if (uncertainty)
	{
	  include.codes <- c(include.codes, countries)
	}
	else
	{
	  include.codes <- c(include.codes, setdiff(countries, countries.processed.est))
	}
					   
	if(length(include.codes) > 0) {
		include <- is.element(tfrs$data$country_code, include.codes)
		locations <- read.UNlocations(tfrs$data, wpp.year=meta$wpp.year, package=package, my.locations.file=my.locations.file, verbose=verbose)
		TFRmatrix.regions <- get.TFRmatrix.and.regions(tfrs$data[include,], locations$loc_data, 
												start.year=meta$start.year, 
												present.year=meta$present.year,
												ignore.last.observed = uncertainty, verbose=verbose, ...)
		processed.include.codes <- intersect(include.codes, country.codes.processed)
		return(list(tfr_matrix=TFRmatrix.regions$tfr_matrix,
				tfr_matrix_all=TFRmatrix.regions$tfr_matrix_all, 
				regions=TFRmatrix.regions$regions,
				nr_countries_estimation=0,
				is_processed=is.element(TFRmatrix.regions$regions$country_code, processed.include.codes)))
	}
	return(NULL)
}

set.wpp.extra <- function(meta, countries=NULL, my.tfr.file=NULL, my.locations.file=NULL, 
                          annual = FALSE, verbose=FALSE, uncertainty=FALSE, use.wpp.data = TRUE) {
	#'countries' is a vector of country or region codes 
	un.object <- read.UNtfr(wpp.year=meta$wpp.year, my.tfr.file=my.tfr.file, 
							present.year=meta$present.year, annual = annual, countries = countries,
							use.wpp.data = use.wpp.data, verbose=verbose)
	data <- un.object$data.object
	extra.wpp <- .extra.matrix.regions(data=data, countries=countries, meta=meta, my.locations.file=my.locations.file, 
	                                   verbose=verbose, annual=annual, uncertainty=uncertainty, 
	                                   interpolate = meta$wpp.year < 2022 && annual && is.null(my.tfr.file)
	                                   )
	if(!is.null(extra.wpp) && !annual) {
		locations <- read.UNlocations(data$data, wpp.year=meta$wpp.year, my.locations.file=my.locations.file, verbose=verbose)
		suppl.wpp <- .get.suppl.matrix.and.regions(un.object, extra.wpp, locations$loc_data, 
									meta$start.year, meta$present.year)
		extra.wpp$suppl.data <- .get.suppl.data.list(suppl.wpp)
	} 
	return(extra.wpp)
}

create.sublocation.dataset <- function(data) {
	loc.idx <- which(is.element(colnames(data), c('name',  'country_code', 'include_code', 'reg_code')))
	loc_data <- data[, loc.idx]
	if(! 'include_code' %in% colnames(data)) data$include_code <- 2
	return(list(loc_data=loc_data, include=data$include_code == 2, 
				prediction.only=data$include_code == 1, include.code=data$include_code))
}

do.read.subnat.file <- function(file.name, present.year=2012) {
	tfr_data <- read.tfr.file(file=file.name)
	if(! 'last.observed' %in% colnames(tfr_data) && 'last_observed' %in% colnames(tfr_data))
	    colnames(tfr_data)[colnames(tfr_data) == 'last_observed'] <- 'last.observed' # guessing a typo, so rename
	if(!is.element('last.observed', colnames(tfr_data)))
		tfr_data <- cbind(tfr_data, last.observed=present.year)
	if(!is.element('include_code', colnames(tfr_data)))
		tfr_data <- cbind(tfr_data, include_code=rep(2, nrow(tfr_data)))
	return(tfr_data)
}

set.wpp.subnat <- function(country, start.year=1950, present.year=2010, annual = FALSE, my.tfr.file=NULL, verbose=FALSE) {
	tfr_data <- do.read.subnat.file(my.tfr.file, present.year = present.year)
	tfr_data <- tfr_data[tfr_data$country_code == country,]
	locations <- create.sublocation.dataset(tfr_data)
	loc_data <- locations$loc_data
	include <- locations$include
	prediction.only <- locations$prediction.only

	tfr_data_countries <- tfr_data[include | locations$prediction.only,]
	nr_countries_estimation <- sum(include)
	TFRmatrix.regions <- get.TFRmatrix.and.regions(tfr_data_countries, loc_data, 
												start.year=start.year, 
												present.year=present.year,
												annual = annual, 
												datacolnames=c(country.code='reg_code', country.name='name', reg.name='reg_name',
															reg.code='NA', area.name='country', area.code='country_code'),
												verbose=verbose)

	return(list(tfr_matrix=TFRmatrix.regions$tfr_matrix, 
				tfr_matrix_all=TFRmatrix.regions$tfr_matrix_all, 
				regions=TFRmatrix.regions$regions, 
				nr_countries_estimation=nr_countries_estimation,
				suppl.data=.get.suppl.data.list(NULL)
				))

}
