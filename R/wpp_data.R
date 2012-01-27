# Read in the UN estimates

set_wpp_regions <- function(start.year=1950, present.year=2010, wpp.year=2010, my.tfr.file=NULL, 
							verbose=FALSE) {
# outputs:
# tfr_matrix_all, with each column one countries UN estimates
# tfr_matrix with NAs after last observed data point
# reg_name, reg_code and area_name, area_code: vectors with UN regions and continents

	########################################
	# set data and match with areas
	########################################
	tfr_data <- read.UNtfr(wpp.year=wpp.year, my.tfr.file=my.tfr.file, 
								present.year=present.year, verbose=verbose)$data
	# not just countries, includes areas etc as well
	# get region and area data
	locations <- read.UNlocations(tfr_data, wpp.year=wpp.year, verbose=verbose)
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
												verbose=verbose)
												
	return(list(tfr_matrix=TFRmatrix.regions$tfr_matrix, 
				tfr_matrix_all=TFRmatrix.regions$tfr_matrix_all, 
				regions=TFRmatrix.regions$regions, 
				nr_countries_estimation=nr_countries_estimation))
}

read.tfr.file <- function(file) return(read.delim(file=file, comment.char='#', check.names=FALSE))

do.read.un.file <- function(un.file.name, wpp.year, my.file=NULL, present.year=2010, verbose=FALSE) {
	cat('Reading file ', un.file.name, '.\n')
	tfr_data <- read.tfr.file(file=un.file.name)
	if(!is.element('last.observed', colnames(tfr_data)))
		tfr_data <- cbind(tfr_data, last.observed=present.year)
	if(!is.element('include_code', colnames(tfr_data)))
		tfr_data <- cbind(tfr_data, include_code=rep(-1, nrow(tfr_data)))
	
	replaced <- c()
	added <- c()
	my.tfr.file <- my.file
	# read user-defined TFRs and replaced the UN data with it
	if(!is.null(my.tfr.file)) {
		cat('Reading file ', my.tfr.file, '.\n')
		my.tfr_data <- read.tfr.file(file=my.tfr.file)
		cols.to.use <- colnames(my.tfr_data)[is.element(colnames(my.tfr_data), colnames(tfr_data))]
		# don't overwrite country_name
		cols.wo.name <- setdiff(cols.to.use, 'country')
		if (!is.element('country_code', cols.to.use))
			stop('Column "country_code" must be included in ', my.tfr.file)
		codes <- my.tfr_data$country_code
		if (length(codes) > 0) {
			for (icode in 1:length(codes)) {
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
				tfr_data[idx,cols.wo.name] <- my.tfr_data[icode,cols.wo.name]
			}
		}
		if(length(replaced) > 0) { 
			cat('\tValues in', paste(cols.to.use, collapse=', '), 'replaced for countries', 
						paste(replaced, collapse=', '), '.\n') 
		} else { 
			cat('\tNo country matched the UN TFR file.\n') 
		}
	}
	return(list(data=tfr_data, replaced=replaced, added=added))
}

read.UNtfr <- function(wpp.year, my.tfr.file=NULL, ...) {
	un.file.name <- file.path(.find.package("bayesTFR"), "data", paste('UN', wpp.year, '.txt', sep=''))
	result <- do.read.un.file(un.file.name, wpp.year, my.file=my.tfr.file, ...)
	return(list(data=result$data, replaced=result$replaced, added=result$added))
}

read.UNlocations <- function(data, wpp.year, package="bayesTFR", verbose=FALSE) {
	loc.file.name <- file.path(.find.package(package), "data", 
								paste('WPP', wpp.year, '_LOCATIONS', '.txt', sep=''))
	cat('Reading file ', loc.file.name, '.\n')
	loc_data <- read.tfr.file(file=loc.file.name)
	
	# this include some areas that are not in the tfr file
	# get the include_code from this file to get the countries from the tfr file
	nr_tfr_outcomes <- length(data[,1])
	include <- rep(NA, nr_tfr_outcomes)
	prediction.only <- rep(NA, nr_tfr_outcomes)
	n_loc_data <- length(loc_data[,1])

	for (i in 1:nr_tfr_outcomes){
 		loc_index <- (1:n_loc_data)[loc_data$country_code == data$country_code[i]]
 		if (length(loc_index)<=0)
 			stop('Country ', data$country_code[i], ' not found in the file ', loc.file.name)
 		incl.code <- if(data$include_code[i] >= 0) data$include_code[i] else loc_data$include_code[loc_index]
 		include[i] <- incl.code == 2
 		prediction.only[i] <- incl.code == 1
	}
	return(list(loc_data=loc_data, include=include, prediction.only=prediction.only, include.code=incl.code))
}


get.observed.time.matrix.and.regions <- function(data, loc_data, start.year=1950, 
										present.year=2010) {
	tfr_data <- data
	nr_countries <- length(tfr_data[,1])
	names.tfr.data <- names(tfr_data)
	#num.columns <- grep('^X[0-9]{4}.[0-9]{4}$', names.tfr.data) # index of year-columns
	num.columns <- grep('^[0-9]{4}.[0-9]{4}$', names.tfr.data) # index of year-columns 
	ncol.tfr <- length(num.columns)
	#cols.starty <- as.integer(substr(names.tfr.data[num.columns], 2,5))
	cols.starty <- as.integer(substr(names.tfr.data[num.columns], 1,4))
	#cols.endy <- as.integer(substr(names.tfr.data[num.columns], 7,10))
	cols.endy <- as.integer(substr(names.tfr.data[num.columns], 6,9))
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
	proj.start.col <- names.tfr.data[num.columns][present.index[1]+1]

	tfr_matrix <- t(tfr_data[,which.max(names(tfr_data)==start.col):which.max(names(tfr_data)==present.col)])
	#start.years <- as.integer(substr(rownames(tfr_matrix), 2,5))
	start.years <- as.integer(substr(rownames(tfr_matrix), 1,4))
	#end.years <- as.integer(substr(rownames(tfr_matrix), 7,10))
	end.years <- as.integer(substr(rownames(tfr_matrix), 6,9))
	mid.years <- start.years + ceiling((end.years- start.years)/2)
	rownames(tfr_matrix) <- mid.years				
	tfr_matrix_all <- tfr_matrix

	reg_name <- reg_code <- area_name <- area_code <- NULL
	for (i in 1:nr_countries){
		loc_index <- which.max(loc_data$country_code == tfr_data$country_code[i])
 		reg_name <- c(reg_name, paste(loc_data$reg_name[loc_index]))
 		reg_code <- c(reg_code, loc_data$reg_code[loc_index])
 		area_name <- c(area_name, paste(loc_data$area_name[loc_index]))
 		area_code <- c(area_code, loc_data$area_code[loc_index])
		# set NAs for entries that are not observed data (after last.observed) 
		tfr_matrix[(tfr_data[i,'last.observed'] < mid.years), i] <- NA
	}
	regions <- list(name=reg_name, code=reg_code, area_name=area_name, area_code=area_code,
				country_name=tfr_data$country, country_code=tfr_data$country_code)
	return(list(obs_matrix=tfr_matrix, obs_matrix_all=tfr_matrix_all, regions=regions))
}

get.TFRmatrix.and.regions <- function(tfr_data, ..., verbose=FALSE){
	result <- get.observed.time.matrix.and.regions(data=tfr_data, ...)
	if (verbose) 
		cat('Dimension of the TFR matrix:', dim(result$obs_matrix), '\n')
	return(list(tfr_matrix=result$obs_matrix, tfr_matrix_all=result$obs_matrix_all, regions=result$regions))
}

.extra.matrix.regions <- function(data, countries, meta, package="bayesTFR", verbose=FALSE) {
	tfrs <- data
	country.codes.processed <- meta$regions$country_code
	ncountries <- get.nr.countries(meta)
	ncountries.est <- get.nr.countries.est(meta)
	replaced.processed <- intersect(tfrs$replaced, country.codes.processed)
	idx.replaced.processed <- is.element(country.codes.processed, replaced.processed)
	replaced.processed.est <- country.codes.processed[idx.replaced.processed[1:ncountries.est]]
	
	added.processed <- intersect(tfrs$added, country.codes.processed)
	idx.added.processed <- is.element(country.codes.processed, added.processed)
	added.processed.est <- country.codes.processed[idx.added.processed[1:ncountries.est]]
	
	countries.processed <- intersect(countries, country.codes.processed)
	idx.countries.processed <- is.element(country.codes.processed, countries.processed)
	countries.processed.est <- country.codes.processed[idx.countries.processed[1:ncountries.est]]

	if(length(replaced.processed.est)+length(added.processed.est)+length(countries.processed.est)> 0) {
		cat('\nCountries', paste(c(replaced.processed.est, added.processed.est, 
									countries.processed.est), collapse=', '), 
					'used for estimation and will be excluded from processing.\n')
	}
	include.codes <- c(setdiff(tfrs$replaced, replaced.processed.est), 
					   setdiff(tfrs$added, added.processed.est),
					   setdiff(countries, countries.processed.est))
	if(length(include.codes) > 0) {
		include <- is.element(tfrs$data$country_code, include.codes)
		locations <- read.UNlocations(tfrs$data, wpp.year=meta$wpp.year, package=package, verbose=verbose)
		TFRmatrix.regions <- get.TFRmatrix.and.regions(tfrs$data[include,], locations$loc_data, 
												start.year=meta$start.year, 
												present.year=meta$present.year,
												verbose=verbose)
		processed.include.codes <- intersect(include.codes, country.codes.processed)
		return(list(tfr_matrix=TFRmatrix.regions$tfr_matrix,
				tfr_matrix_all=TFRmatrix.regions$tfr_matrix_all, 
				regions=TFRmatrix.regions$regions,
				nr_countries_estimation=0,
				is_processed=is.element(TFRmatrix.regions$regions$country_code, processed.include.codes)))
	}
	return(NULL)
}

set.wpp.extra <- function(meta, countries=NULL, my.tfr.file=NULL, verbose=FALSE) {
	#'countries' is a vector of country or region codes 
	tfrs <- read.UNtfr(wpp.year=meta$wpp.year, my.tfr.file=my.tfr.file, 
							present.year=meta$present.year, verbose=verbose)
	return(.extra.matrix.regions(tfrs, countries, meta, verbose=verbose))
}
