tfr.predict.subnat <- function(countries, my.tfr.file, sim.dir=file.path(getwd(), 'bayesTFR.output'),
                               end.year=2100, start.year=NULL, output.dir = NULL, annual = NULL,
                               nr.traj=NULL, seed = NULL, min.tfr = 0.5, ar.pars = NULL, 
                               save.as.ascii = 0, verbose = TRUE) {
  # Run subnational projections, using the Scale AR(1) model applied to a national bayesTFR simulation 
  # sim.dir is the world-national simulation. Set output.dir to store results somewhere else.  

  wpred <- get.tfr.prediction(sim.dir) # contains national projections
  wdata <- wpred$tfr_matrix_reconstructed
  wmeta <- wpred$mcmc.set$meta
  if(!is.null(seed)) set.seed(seed)
  if (is.null(output.dir)) output.dir <- wmeta$output.dir
  quantiles.to.keep <- as.numeric(dimnames(wpred$quantiles)[[2]])
  wannual <- wmeta$annual.simulation
  if(is.null(annual)) annual <- wannual
  
  ar.pars.default <- c(mu = 1, rho = 0.92464, sigma = 0.04522)
  if(annual) { # Raftery's conversion from 5-year AR(1) parameters to 1-year parameters
    ar.pars.default["sigma"] <- ar.pars.default["sigma"] * sqrt((1-ar.pars.default["rho"]^(2/5))/(1-ar.pars.default["rho"]^2))
    ar.pars.default["rho"] <- ar.pars.default["rho"]^(1/5)
  }
  
  result <- list()
  orig.nr.traj <- nr.traj
  for (country in countries) {
    country.obj <- get.country.object(country, wmeta)
    if(is.null(country.obj$code)) {
      warning("Country ", country, " not found in the national projections. Use numerical codes or exact names.")
      next
    }
    if(verbose) 
      cat('\nPredicting TFR for ', country.obj$name, '\n')
    year.step <- if(annual) 1 else 5
    wyear.step <- if(wannual) 1 else 5
    
    wtfr <- wdata[,country.obj$index]
    wtfrobsy <- as.integer(names(wtfr))[-length(wtfr)]
    
    # set start year and present year of the national and subnational simulations
    wstarty <- if(is.null(start.year)) as.integer(dimnames(wpred$tfr_matrix_reconstructed)[[1]][wpred$present.year.index])+wyear.step else start.year
    wpresenty <- wstarty - wyear.step
    
    starty <- wstarty # set subnational start year to the national start year (might get overwritten below)
    presenty <- starty - year.step # subnational present year (might get overwritten below)
    
    if(!wannual){ # national prediction is 5-year
      # snap start year and present year to be in the middle of the corresponding 5-years interval
      seqy <- seq(min(wtfrobsy)-3, wpred$end.year, by=5)
      wstarty <- seqy[cut(wstarty, seqy, labels=FALSE)] + 3
      wpresenty <- seqy[cut(wstarty, seqy, labels=FALSE)-1] + 3
      if(!annual){ # subnational prediction is 5-year
        starty <- wstarty
        presenty <- wpresenty
      }
    } else { # national prediction is annual
        if(!annual){ # subnational prediction is 5-year
          # adjust start year and present year to be in the middle of the corresponding 5-years interval
          seqy <- seq(1700, wpred$end.year, by=5)
          starty <- seqy[cut(wstarty, seqy, labels=FALSE)] + 3
          presenty <- starty - 5
        }
    }
    
    # get various meta values for the subnational simulation
    meta <- mcmc.meta.ini.subnat(wmeta, country = country.obj$code, my.tfr.file = my.tfr.file, 
                                 start.year = 1750, present.year = presenty, annual = annual, verbose = verbose)
    
    if(is.null(meta$tfr_matrix_all))
      stop("No data available. It might be a mismatch between the available observed years and start.year. ",
           "The data file should contain year ", if(annual) presenty else paste(presenty - 3, presenty + 2, sep = "-"),
           ". Or change the argument start.year.")
    
    this.output.dir <- file.path(output.dir, 'subnat', paste0('c', country.obj$code))
    outdir <- file.path(this.output.dir, 'predictions')
    meta$output.dir <- this.output.dir
    
    # extract national trajectories and thin them if needed
    wtrajs <- get.tfr.trajectories(wpred, country.obj$code) 
    nr.traj <- orig.nr.traj
    if(is.null(nr.traj)) nr.traj <- ncol(wtrajs)
    thinning.index <- round(seq(1, ncol(wtrajs), length=nr.traj))
    wtrajs <- wtrajs[as.integer(rownames(wtrajs)) <= end.year, thinning.index]
    nr.traj <- ncol(wtrajs)
    
    # attach observed data to the trajectories, so that imputation can be done using all data 
    addtfr <- matrix(wtfr, ncol = nr.traj, nrow = length(wtfr), 
                        dimnames = list(names(wtfr), colnames(wtrajs)))
    wtrajs.all <- rbind(addtfr[! rownames(addtfr) %in% rownames(wtrajs),], wtrajs)

    if(annual && !wannual) { # interpolate national trajectories
      wtrajs.all <- interpolate.trajectories(wtrajs.all)
    } else {
      if(!annual && wannual) # average annual national trajectories
        wtrajs.all <- average.trajectories(wtrajs.all)
    }
    # remove everything from the national trajectories before present year
    yrs <- as.integer(rownames(wtrajs.all))  
    wtrajs <- wtrajs.all[-which(yrs < presenty),]
  
    if(!presenty %in% rownames(meta$tfr_matrix_all)) { # add NAs to tfr matrix
      addyears <- sort(seq(presenty, max(as.integer(rownames(meta$tfr_matrix_all))+year.step), by=-year.step))
      adddata <- matrix(NA, nrow=length(addyears), ncol=ncol(meta$tfr_matrix_all),
                        dimnames=list(addyears, colnames(meta$tfr_matrix_all)))
      for(tfrname in c("tfr_matrix_all", "tfr_matrix_observed", "tfr_matrix"))
        meta[[tfrname]] <- rbind(meta[[tfrname]], adddata)
    }
    # save meta to disk
    if(file.exists(this.output.dir)) unlink(this.output.dir, recursive=TRUE)
    dir.create(outdir, recursive=TRUE)
    bayesTFR.mcmc.meta <- meta
    store.bayesTFR.meta.object(bayesTFR.mcmc.meta, this.output.dir)
    
    # prepare for subnational simulation
    this.nr.project <- nrow(wtrajs)-1
    nr.reg <- get.nr.countries(meta)
    PIs_cqp <- array(NA, c(nr.reg, length(quantiles.to.keep), nrow(wtrajs)),
                     dimnames=list(meta$regions$country_code, dimnames(wpred$quantiles)[[2]], dimnames(wtrajs)[[1]]))
    mean_sd <- array(NA, c(nr.reg, 2, nrow(wtrajs)))
    meta$T_end <- nrow(meta$tfr_matrix_all)
    meta$T_end_c <- rep(meta$T_end, nr.reg)
    country.char <- as.character(country.obj$code)
    
    # Make sure all AR(1) parameters are available. If not, take the default values.
    if(is.null(ar.pars)) ar.pars <- ar.pars.default
    for(par in names(ar.pars.default)) if(! par %in% names(ar.pars)) ar.pars[par] <- ar.pars.default[par]
      
    # Get observed variance (SD) from a time point that has at least one half of the regions non-NA
    for(i in meta$T_end:1) if(sum(!is.na(meta$tfr_matrix_observed)[i,]) >= max(2, nr.reg/2)) break # check NA's at the end
    widx <- which(rownames(wtrajs.all) == rownames(meta$tfr_matrix_observed)[i])
    arsd <- min(ar.pars['sigma'], sqrt((1-ar.pars['rho']^2)*var(meta$tfr_matrix_observed[i,]/mean(wtrajs.all[widx,]), na.rm = TRUE)))
    tfr_reconstructed <- meta$tfr_matrix_observed
    
    # iterate over subnational units
    for(region in 1:nr.reg) {
        reg.obj <- get.country.object(region, meta, index=TRUE)
        regcode.char <- as.character(reg.obj$code)			
        regtfr <- get.observed.tfr(region, meta, 'tfr_matrix_observed')
        regtfr.last <- rep(regtfr[length(regtfr)], nr.traj)
        c.first <- regtfr.last/wtrajs.all[rownames(wtrajs.all) %in% names(regtfr.last)] # set of initial scales
        if(is.na(regtfr.last[1])) { # impute where NA's at the end of the time series
            for(i in length(regtfr):1) if(!is.na(regtfr[i])) break # find the last non-NA value
            widx <- which(rownames(wtrajs.all) %in% names(regtfr[i])) # match the year with the national TFR
            c.first <- regtfr[i]/wtrajs.all[widx,] # set of initial scales
            meta$T_end_c[region] <- i
            tmptraj <- matrix(NA, nrow = length(regtfr) - i, ncol = nr.traj) # trajectory matrix for imputation
            for(tr in 1:nr.traj) { # iterate over trajectories
                scale.prev <- c.first[tr]
                wi <- 1 # counter for moving along the time of world trajectories
                for(j in (i+1):length(regtfr)) { # iterate over imputation time
                    for(k in 1:100)	{ # loop for getting tfr above the required minimum (0.5)
                        scale <- ar.pars['rho'] * (scale.prev - ar.pars['mu']) + ar.pars['mu'] + rnorm(1, 0, sd=arsd)
                        tmptraj[j-i, tr] <- wtrajs.all[widx+wi,tr]*scale
                        if(tmptraj[j-i, tr] > min.tfr) break
                    }
                    tmptraj[j-i, tr] <- max(min.tfr, tmptraj[j-i, tr])
                    scale.prev <- scale
                    wi <- wi + 1
                }
                # set the starting scale of the projection to the end scale of the imputation (for this trajectory)
                c.first[tr] <- scale 
            }
            # impute using median
            regtfr[(i+1):length(regtfr)] <- tfr_reconstructed[(i+1):length(regtfr),region] <- apply(tmptraj, 1, median)
            regtfr.last <- tmptraj[nrow(tmptraj),]
        } # end of imputation
        
        # initiate matrix for holding the projections
        tfr.pred <- matrix(NA, nrow=this.nr.project+1, ncol=nr.traj)
        tfr.pred[1,] <- regtfr.last
      
      for(s in 1:ncol(tfr.pred)) { # iterate over trajectories
        scale.prev <- c.first[s]
        for(year in 2:(this.nr.project+1)) { # iterate over time
          for(i in 1:100)	{ # keep simulating if tfr below the limit
            scale <- ar.pars['rho'] * (scale.prev - ar.pars['mu']) + ar.pars['mu'] + rnorm(1, 0, sd=arsd) # AR1
            tfr.pred[year, s] <- wtrajs[year,s] * scale 
            if(tfr.pred[year, s] > min.tfr) break # lower limit for tfr is 0.5
          }
          tfr.pred[year, s] <- max(min.tfr, tfr.pred[year, s])
          scale.prev <- scale
        }
      }
      # save resulting trajectories
      trajectories <- tfr.pred
      rownames(trajectories) <- rownames(wtrajs)
      save(trajectories, file = file.path(outdir, paste0('traj_country', reg.obj$code, '.rda')))
      
      # compute quantiles
      PIs_cqp[region,,] <- apply(trajectories, 1, quantile, quantiles.to.keep, na.rm = TRUE)
      mean_sd[region,1,] <- apply(trajectories, 1, mean, na.rm = TRUE)
      mean_sd[region,2,] <- apply(trajectories, 1, sd, na.rm = TRUE) 	
    }
    present.year.index <- which(rownames(meta$tfr_matrix_all) == rownames(wtrajs)[1])

    bayesTFR.prediction <- structure(list(
      quantiles = PIs_cqp,
      traj.mean.sd = mean_sd,
      nr.traj=nr.traj,
      tfr_matrix_reconstructed = tfr_reconstructed,
      output.directory = normalizePath(outdir),
      na.index=wpred$na.index,
      mcmc.set=list(meta=meta, mcmc.list=list()),
      nr.projections=this.nr.project,
      burnin=NA, thin=NA,
      end.year = max(as.integer(rownames(wtrajs))), use.tfr3=FALSE, burnin3=NA, thin3=NA,
      mu=ar.pars['mu'], rho=ar.pars['rho'], sigmaAR1 = arsd, mu.c=NA, rho.c=NA, min.tfr=min.tfr,
      use.correlation=NA, start.year=starty,
      present.year.index=present.year.index,
      present.year.index.all=present.year.index),
      class='bayesTFR.prediction')
    store.bayesTFR.prediction(bayesTFR.prediction, outdir)
    do.convert.trajectories(pred=bayesTFR.prediction, n=save.as.ascii, output.dir=outdir, verbose=verbose)
    result[[as.character(country.obj$code)]] <- bayesTFR.prediction
  }
  cat('\nPrediction stored into', output.dir, '\n')
  invisible(result)
}

interpolate.trajectories <- function(trajs) {
  years <- as.integer(rownames(trajs))
  years.to.interp.to <- (min(years) - 3):max(years) + 2
  traj.new <- matrix(NA, ncol = ncol(trajs), nrow = length(years.to.interp.to))
  rownames(traj.new) <- years.to.interp.to
  for(col in 1:ncol(traj.new)) 
    traj.new[, col] <- approx(years, trajs[, col], xout = years.to.interp.to, rule=2)$y
  return(traj.new)
}

average.trajectories <- function(trajs) {
  # convert annual TFR into 5-year series
  years <- as.integer(rownames(trajs))
  ranges <- range(years[years %% 5 == 0])
  mid.points <- c(0, seq(ranges[1]-2, ranges[2]+3, by = 5))
  brks <- seq(ranges[1]-5, ranges[2] + 5, by = 5)
  year.bin <- findInterval(years, brks, left.open = TRUE)
  traj.new <- apply(trajs, 2, function(tfr) aggregate(tfr, by = list(year.bin), FUN = mean, na.rm = TRUE)[,"x"])
  rownames(traj.new) <- mid.points[-c(1, length(mid.points))]
  return(traj.new)
}
