if(getRversion() >= "2.15.1") utils::globalVariables(c("counter", "counter3"))

store.mcmc <- local({
	# Writes parameter values into ascii files - one file per parameter and country (if country-specific)
	##########################
	par.names <- c(tfr.parameter.names(trans=FALSE))
	par.cs.names <- c(tfr.parameter.names.cs(trans=FALSE, back.trans=FALSE), 'eps_T')
	var.names <- list(gamma='gamma_ci', d='d_c', Triangle_c4='Triangle_c4', eps_T='eps_Tc', U='U_c')
	
	default.buffer.size <- 10
	buffer <- buffer.cs <- NULL
	
	# Country-specific parameters are buffered in one matrix per parameter, 
	# with one row per iteration and a block of columns per country 
	# (countries given by cs.index; each block has cs.width[[par]] columns).
	cs.index <- cs.index.tfr <- cs.width <- NULL
	
	get.cs.values <- function(mcmc, par, country.index) {
		# values of a country-specific parameter for the given countries, concatenated by country
		return(switch(par,
			gamma = as.vector(t(mcmc$gamma_ci[country.index,,drop=FALSE])),
			eps_T = as.vector(mcmc$eps_Tc[,country.index,drop=FALSE]),
			tfr = as.vector(mcmc$meta[['tfr_all']][,country.index,drop=FALSE]),
			mcmc[[var.names[[par]]]][country.index]))
	}
	
	get.cs.width <- function(mcmc, par) {
		# number of values per country
		return(switch(par, 
			gamma = ncol(mcmc$gamma_ci),
			eps_T = nrow(mcmc$eps_Tc),
			tfr = nrow(mcmc$meta[['tfr_all']]),
			1))
	}
	
	has.tfr <- function(mcmc) return(!is.null(mcmc$uncertainty) && mcmc$uncertainty)
	
	buffers.insert <- function(mcmc, countries=NULL) {
	  counter <<- counter + 1
		if (is.null(countries)) {
			for (par in par.names) {
				if (is.element(par, mcmc$dontsave)) next
				buffer[[par]][counter,] <<- mcmc[[par]]
			}
		  if (!is.null(mcmc$meta$ar.phase2) && mcmc$meta$ar.phase2) buffer[['rho.phase2']][counter,] <<- mcmc[["rho.phase2"]]
		}
		for (par in names(buffer.cs))
			buffer.cs[[par]][counter,] <<- get.cs.values(mcmc, par, if(par == 'tfr') cs.index.tfr else cs.index)
	}
		
	buffers.ini <- function(mcmc, size, countries=NULL) {
		buffer <<- list()
		if (is.null(countries)) {
			for (par in par.names) {
				if (is.element(par, mcmc$dontsave)) next
				buffer[[par]] <<- matrix(NA, ncol=length(mcmc[[par]]), nrow=size)
			}
		    if (!is.null(mcmc$meta$ar.phase2) && mcmc$meta$ar.phase2) 
		        buffer[['rho.phase2']] <<- matrix(NA, ncol=1, nrow=size)
			cs.index <<- mcmc$meta$id_DL
			cs.index.tfr <<- 1:mcmc$meta$nr_countries
		} else {
			cs.index <<- cs.index.tfr <<- countries
		}
		buffer.cs <<- list()
		cs.width <<- list()
		cs.pars <- par.cs.names[!sapply(par.cs.names, function(par) is.element(var.names[[par]], mcmc$dontsave))]
		if (has.tfr(mcmc)) cs.pars <- c(cs.pars, 'tfr')
		for (par in cs.pars) {
			cs.width[[par]] <<- get.cs.width(mcmc, par)
			ncountries <- length(if(par == 'tfr') cs.index.tfr else cs.index)
			buffer.cs[[par]] <<- matrix(NA, ncol=cs.width[[par]] * ncountries, nrow=size)
		}
		counter <<- 0
	}
	
	
	do.flush.buffers <- function(mcmc, append=FALSE, countries=NULL, verbose=FALSE) {
		if (verbose)
			cat("Flushing results into disk.\n")
		output.dir <- file.path(mcmc$meta$output.dir, mcmc$output.dir)
		if(!file.exists(output.dir)) 
			dir.create(output.dir)
		open <- if(append) 'a' else 'w'
		if (is.null(countries)) {
		  for(par in par.names) { # write country-independent parameters
				if (is.null(buffer[[par]])) next
				if (counter == 1) {
					values <- t(buffer[[par]][1:counter,])
				} else {
					values <- buffer[[par]][1:counter,]
				}
				write.values.into.file.cindep(par, values, output.dir, mode=open, 
												compression.type=mcmc$compression.type)
		  }
		  if (!is.null(mcmc$meta$ar.phase2) && mcmc$meta$ar.phase2)
		  {
		    if (counter == 1) {
		      values <- t(buffer[['rho.phase2']][1:counter,])
		    } else {
		      values <- buffer[['rho.phase2']][1:counter,]
		    }
		    write.values.into.file.cindep('rho_phase2', values, output.dir, mode=open, 
		                                  compression.type=mcmc$compression.type)
		  }
		}
		for (par in names(buffer.cs)) { # write country-specific parameters (and tfr if uncertainty)
			w <- cs.width[[par]]
			country.index <- if(par == 'tfr') cs.index.tfr else cs.index
			for (i in seq_along(country.index)){
				values <- buffer.cs[[par]][1:counter, (i-1)*w + (1:w), drop=FALSE]
				write.values.into.file.cdep(par, values, output.dir, 
						get.country.object(country.index[i], meta=mcmc$meta, index=TRUE)$code, mode=open, 
											compression.type=mcmc$compression.type)
			}
		}
		
		resmc <- as.list(mcmc)
		class(resmc) <- 'bayesTFR.mcmc'
		store.bayesTFR.object(resmc, output.dir)
	}
	
	store <- function(mcmc, append=FALSE, flush.buffer=FALSE, countries=NULL, verbose=FALSE) {
		# If countries is not NULL, only country-specific parameters 
		# for those countries (given as index) are stored
		buffer.size <- mcmc$meta$buffer.size
		if (is.null(buffer.size)) buffer.size <- default.buffer.size
		if (is.null(buffer)) buffers.ini(mcmc, buffer.size, countries=countries)
		buffers.insert(mcmc, countries=countries)
		flushed <- FALSE
		if (flush.buffer || (counter >= buffer.size)) {
			do.flush.buffers(mcmc, append=append, countries=countries, verbose=verbose)
			buffer <<- buffer.cs <<- NULL
			flushed <- TRUE
		}
		return(flushed)
	}

})

store.mcmc3 <- local({
	# Writes parameter values into ascii files - one file per parameter and country (if country-specific)
	##########################
	par.names <- tfr3.parameter.names()
	par.cs.names <- tfr3.parameter.names.cs()
	
	default.buffer.size <- 10
	buffer3 <- buffer3.cs <- NULL
		
	# Country-specific parameters are buffered in one matrix per parameter, 
	# with one row per iteration and a block of columns per country 
	# (countries given by cs3.index; each block has cs3.width[[par]] columns).
	cs3.index <- cs3.width <- NULL
	
	get.cs.values <- function(mcmc, par, country.index) {
		# values of a country-specific parameter for the given countries, concatenated by country
		if(is.null(dim(mcmc[[par]]))) return(mcmc[[par]][country.index])
		return(as.vector(mcmc[[par]][,country.index,drop=FALSE]))
	}
	
	buffers.insert <- function(mcmc, countries=NULL) {
	  counter3 <<- counter3 + 1
		if (is.null(countries)) {
			for (par in par.names) buffer3[[par]][counter3,] <<- mcmc[[par]]
		}
		for (par in par.cs.names)
			buffer3.cs[[par]][counter3,] <<- get.cs.values(mcmc, par, cs3.index)
	}
		
	buffers.ini <- function(mcmc, size, countries=NULL) {
		buffer3 <<- list()
		if (is.null(countries)) {
			for (par in par.names) 
				buffer3[[par]] <<- matrix(NA, ncol=length(mcmc[[par]]), nrow=size)
			cs3.index <<- 1:mcmc$meta$nr.countries
		} else cs3.index <<- countries
		buffer3.cs <<-list()
		cs3.width <<- list()
		for (par in par.cs.names) {
			cs3.width[[par]] <<- if(is.null(dim(mcmc[[par]]))) 1 else nrow(mcmc[[par]])
			buffer3.cs[[par]] <<- matrix(NA, ncol=cs3.width[[par]] * length(cs3.index), nrow=size)
		}
		counter3 <<- 0
	}
	
	do.flush.buffers <- function(mcmc, append=FALSE, countries=NULL, verbose=FALSE) {
	  if (verbose)
			cat("Flushing results into disk.\n")
	  if (!is.null(mcmc$uncertainty) && (mcmc$uncertainty) )
	  {
	    output.dir <- file.path(mcmc$meta$output.dir, "phaseIII", mcmc$output.dir)
	  }
	  else
	  {
	    output.dir <- file.path(mcmc$meta$output.dir, mcmc$output.dir)
	  }
	  if(!file.exists(output.dir)) 
			dir.create(output.dir)
		open <- if(append) 'a' else 'w'
		
		if (is.null(countries)) {
			for(par in par.names) { # write country-independent parameters
				if (is.null(buffer3[[par]])) next
				values <- if (counter3 == 1) t(buffer3[[par]][1:counter3,])
				 			else buffer3[[par]][1:counter3,]
				write.values.into.file.cindep(par, values, output.dir, mode=open, 
												compression.type=mcmc$compression.type)
			}
		}

		for (par in par.cs.names) { # write country-specific parameters
			if (is.null(buffer3.cs[[par]])) next
			w <- cs3.width[[par]]
			for (i in seq_along(cs3.index)){
				values <- buffer3.cs[[par]][1:counter3, (i-1)*w + (1:w), drop=FALSE]
				write.values.into.file.cdep(par, values, output.dir, 
						get.country.object(mcmc$meta$id_phase3[cs3.index[i]], meta=mcmc$meta$parent, index=TRUE)$code, mode=open, 
											compression.type=mcmc$compression.type)
			}
		}
		resmc <- as.list(mcmc)
		class(resmc) <- 'bayesTFR.mcmc'
		store.bayesTFR.object(resmc, output.dir)
	}
	
	store <- function(mcmc, append=FALSE, flush.buffer=FALSE, countries=NULL, verbose=FALSE) {
	  # If countries is not NULL, only country-specific parameters 
		# for those countries (given as index) are stored
	  buffer.size <- mcmc$meta$buffer.size
		if (is.null(buffer.size)) buffer.size <- default.buffer.size
		if (is.null(buffer3)) buffers.ini(mcmc, buffer.size, countries=countries)
		buffers.insert(mcmc, countries=countries)
		flushed <- FALSE
		if (flush.buffer || (counter3 >= buffer.size)) {
			do.flush.buffers(mcmc, append=append, countries=countries, verbose=verbose)
			buffer3 <<- buffer3.cs <<- NULL
			flushed <- TRUE
		}
		return(flushed)
	}

})

.get.compression.settings.obsolete <- function(compression.type='None') {
    if(is.null(compression.type)) compression.type <- 'None'
    return(switch(compression.type,
                  None=c('file', '', ''),
                  xz = c('xzfile', '.xz', 'b'),
                  bz = c('bzfile', '.bz2','b'),
                  gz = c('gzfile', '.gz', 'b')))
}

.get.compression.settings <- function(compression.type='None') {
	if(is.null(compression.type)) compression.type <- 'None'
	return(switch(compression.type,
							None=c('none', '', ''),
							xz = c('xzfile', '.xz', 'b'),
							bz = c('bzfile', '.bz2','b'),
							gz = c('gzip', '.gz', '')))
}

do.write.values.into.file <- function(filename, data, mode, compression.type='None') {
	cmd.suffix.mode <- .get.compression.settings(compression.type)
	# con <- do.call(cmd.suffix.mode[1], list(paste(filename, cmd.suffix.mode[2], sep=''), 
	#                                         open=paste(mode, cmd.suffix.mode[3], sep='')))
	# write.table(data, file=con, row.names=FALSE, col.names = FALSE, sep=" ")
	# close(con)
	# # return()
	if(cmd.suffix.mode[1] %in% c("xzfile", "bzfile")) {
	    #con <- bzfile(filename, open=mode)
	    con <- do.call(cmd.suffix.mode[1], list(paste(filename, cmd.suffix.mode[2], sep=''), 
	    			open=paste(mode, cmd.suffix.mode[3], sep='')))
	    write.table(data, file=con, row.names=FALSE, col.names = FALSE, sep=" ")
	    close(con)
	} else
    data.table::fwrite(data.table::data.table(data), file = paste0(filename, cmd.suffix.mode[2]), sep = "\t",
                       showProgress = FALSE, compress = cmd.suffix.mode[1], append = mode == "a", col.names = FALSE)
}

write.values.into.file.cindep <- function(par, data, output.dir, mode='w', compression.type='None') {
	do.write.values.into.file(file.path(output.dir, paste(par,'txt', sep='.')), data, mode=mode, 
									compression.type=compression.type)
}

write.table.into.file.cindep <- function(data, ...) {
	for (par in colnames(data))
		write.values.into.file.cindep(par, data[,par], mode='w', ...)
}

write.list.into.file.cindep <- function(mcmc, data.list) {
	output.dir <- file.path(mcmc$meta$output.dir, mcmc$output.dir)
	if(!file.exists(output.dir)) 
		dir.create(output.dir)
	for (par in names(data.list))
		write.values.into.file.cindep(par, data.list[[par]], output.dir, mode='w', compression.type=mcmc$meta$compression.type)
}

write.values.into.file.cdep <- function(par, data, output.dir, country.code, mode='w', compression.type='None') {
	do.write.values.into.file(file.path(output.dir, paste(par,"_country", country.code, ".txt",sep = "")), 
									data, mode=mode, compression.type=compression.type)
}

write.table.into.file.cdep <- function(data, ...) {
	for (par in colnames(data))
		write.values.into.file.cdep(par, data[,par], mode='w', ...)
}

store.bayesTFR.object <- function(mcmc, output.dir) {
	bayesTFR.mcmc <- mcmc
	for (item in bayesTFR.mcmc$dontsave)  # don't save meta and some other data
		bayesTFR.mcmc[[item]] <- NULL
	bayesTFR.mcmc$meta <- NULL
	bayesTFR.mcmc$eps.T.index <- NULL # caches used during sampling
	bayesTFR.mcmc$proposal.gamma.factor <- NULL
	save(bayesTFR.mcmc, file=file.path(output.dir, 'bayesTFR.mcmc.rda'))
}

store.bayesTFR.meta.object <- function(meta, output.dir) {
	bayesTFR.mcmc.meta <- meta
	save(bayesTFR.mcmc.meta, file=file.path(output.dir, 'bayesTFR.mcmc.meta.rda'))
}

store.bayesTFR.prediction <- function(pred, output.dir=NULL) {
	bayesTFR.prediction <- pred
	if (is.null(output.dir)) output.dir <- pred$output.directory
	save(bayesTFR.prediction, file=file.path(output.dir, 'prediction.rda'))
}

store.bayesTFR.convergence <- function(diag, thin, burnin, output.dir){
	save.file <- file.path(output.dir, paste('bayesTFR.convergence_', thin, '_', burnin, '.rda', sep=''))
	bayesTFR.convergence <- diag
	save(bayesTFR.convergence, file=save.file)
	return(save.file)
}