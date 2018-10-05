tfr.predict.subnat <- function(countries, my.tfr.file, sim.dir=file.path(getwd(), 'bayesTFR.output'),
                               end.year=2100, start.year=NULL, output.dir = NULL, nr.traj=NULL, seed = NULL, 
                               min.tfr=0.5, ar.pars=c(mu = 1, rho = 0.92464, sigma = 0.04522), 
                               save.as.ascii=0, verbose=TRUE) {
  # Run subnational projections, using the Scale AR(1) model applied to a national bayesTFR simulation 
  # sim.dir is the world-national simulation. Set output.dir to store results somewhere else.  

  wpred <- get.tfr.prediction(sim.dir) # contains national projections
  wdata <- wpred$tfr_matrix_reconstructed
  wmeta <- wpred$mcmc.set$meta
  if(!is.null(seed)) set.seed(seed)
  if (is.null(output.dir)) output.dir <- wmeta$output.dir
  quantiles.to.keep <- as.numeric(dimnames(wpred$quantiles)[[2]])
  e <- new.env()
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
    starty <- if(is.null(start.year)) as.integer(dimnames(wpred$tfr_matrix_reconstructed)[[1]][wpred$present.year.index])+5
    else start.year
    meta <- mcmc.meta.ini.subnat(wmeta, country=country.obj$code, my.tfr.file=my.tfr.file, 
                                 start.year=1750, present.year=starty-5, verbose=verbose)
    this.output.dir <- file.path(output.dir, 'subnat', paste0('c', country.obj$code))
    outdir <- file.path(this.output.dir, 'predictions')
    meta$output.dir <- this.output.dir
    
    wtrajs <- get.tfr.trajectories(wpred, country.obj$code)
    nr.traj <- orig.nr.traj
    if(is.null(nr.traj)) nr.traj <- ncol(wtrajs)
    thinning.index <- round(seq(1, ncol(wtrajs), length=nr.traj))
    wtrajs <- wtrajs[as.integer(rownames(wtrajs)) <= end.year, thinning.index]
    nr.traj <- ncol(wtrajs)
    wyears <- as.integer(rownames(wtrajs))
    wend.year <- max(wyears)
    wtfr <- wdata[,country.obj$index]
    wtfrobsy <- as.integer(names(wtfr))[-length(wtfr)]
    seqy <- seq(min(wtfrobsy)-3, max(wyears)+2, by=5)
    midy <- seqy+3
    presenty <- midy[cut(starty, seqy, labels=FALSE)-1]
    if(any(wyears < presenty)) { # remove time periods from national trajectories before present year
      wtrajs <- wtrajs[-which(wyears < presenty),]
    }
    if(presenty < min(wyears)) { # add observed data to national trajectories if present year is not there
      addtfr <- wtfr[-length(wtfr)][wtfrobsy >= presenty]
      adddata <- matrix(addtfr, ncol=nr.traj, nrow=length(addtfr), 
                        dimnames=list(names(addtfr), colnames(wtrajs)))
      wtrajs <- rbind(adddata, wtrajs)
    }
    if(!presenty %in% rownames(meta$tfr_matrix_all)) { # add NAs to tfr matrix
      addyears <- sort(seq(presenty, max(as.integer(rownames(meta$tfr_matrix_all))+5), by=-5))
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
    
    this.nr.project <- nrow(wtrajs)-1
    nr.reg <- get.nr.countries(meta)
    PIs_cqp <- array(NA, c(nr.reg, length(quantiles.to.keep), nrow(wtrajs)),
                     dimnames=list(meta$regions$country_code, dimnames(wpred$quantiles)[[2]], dimnames(wtrajs)[[1]]))
    mean_sd <- array(NA, c(nr.reg, 2, nrow(wtrajs)))
    meta$T_end <- nrow(meta$tfr_matrix_all)
    meta$T_end_c <- rep(meta$T_end, nr.reg)
    country.char <- as.character(country.obj$code)
    # Get observed variance (SD) from a time point that has at least one half of the regions non-NA
    for(i in meta$T_end:1) if(sum(!is.na(meta$tfr_matrix_observed)[i,]) >= max(2, nr.reg/2)) break # check NA's at the end
    widx <- which(names(wtfr) == rownames(meta$tfr_matrix_observed)[i])
    arsd <- min(ar.pars['sigma'], sqrt((1-ar.pars['rho']^2)*var(meta$tfr_matrix_observed[i,]/wtfr[widx], na.rm = TRUE)))
    tfr_reconstructed <- meta$tfr_matrix_observed
    for(region in 1:nr.reg) {
        reg.obj <- get.country.object(region, meta, index=TRUE)
        regcode.char <- as.character(reg.obj$code)			
        regtfr <- get.observed.tfr(region, meta, 'tfr_matrix_observed')
        regtfr.last <- rep(regtfr[length(regtfr)], nr.traj)
        c.first <- regtfr.last/wtfr[names(wtfr) %in% names(regtfr.last)] # set of initial scales
        if(is.na(regtfr.last[1])) { # impute where NA's at the end
            for(i in length(regtfr):1) if(!is.na(regtfr[i])) break
            widx <- which(names(wtfr) %in% names(regtfr[i]))
            c.first <- rep(regtfr[i]/wtfr[widx], nr.traj) # set of initial scales
            meta$T_end_c[region] <- i
            tmptraj <- matrix(NA, nrow = length(regtfr) - i, ncol = nr.traj) # trajectory matrix for imputation
            for(tr in 1:nr.traj) { # iterate over trajectories
                scale.prev <- c.first[tr]
                wi <- 1 # counter for moving along the time of world trajectories
                for(j in (i+1):length(regtfr)) { # iterate over imputation time
                    for(k in 1:100)	{ # loop for getting tfr above the required minimum (0.5)
                        scale <- ar.pars['rho'] * (scale.prev - ar.pars['mu']) + ar.pars['mu'] + rnorm(1, 0, sd=arsd)
                        tmptraj[j-i, tr] <- wtfr[widx+wi]*scale
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
            regtfr.last <- tmptraj[,length(regtfr)]
        } # end of imputation
        
        tfr.pred <- matrix(NA, nrow=this.nr.project+1, ncol=nr.traj)
        tfr.pred[1,] <- regtfr.last
      
      for(s in 1:ncol(tfr.pred)) { # iterate over trajectories
        #is.in.phase3 <- is.reg.phase3
        scale.prev <- c.first[s]
        for(year in 2:(this.nr.project+1)) {		
          for(i in 1:100)	{
            scale <- ar.pars['rho'] * (scale.prev - ar.pars['mu']) + ar.pars['mu'] + rnorm(1, 0, sd=arsd) # AR1
            tfr.pred[year, s] <- wtrajs[year,s] * scale 
            if(tfr.pred[year, s] > min.tfr) break # lower limit for tfr is 0.5
          }
          tfr.pred[year, s] <- max(min.tfr, tfr.pred[year, s])
          scale.prev <- scale
        }
      }
      trajectories <- tfr.pred
      save(trajectories, file = file.path(outdir, paste0('traj_country', reg.obj$code, '.rda')))
      # compute quantiles
      PIs_cqp[region,,] <- apply(trajectories, 1, quantile, quantiles.to.keep, na.rm = TRUE)
      mean_sd[region,1,] <- apply(trajectories, 1, mean, na.rm = TRUE)
      mean_sd[region,2,] <- apply(trajectories, 1, sd, na.rm = TRUE) 	
    }
    present.year.index <- which(rownames(meta$tfr_matrix_all) == rownames(wtrajs)[1])
    #stop('')
    bayesTFR.prediction <- structure(list(
      quantiles = PIs_cqp,
      traj.mean.sd = mean_sd,
      nr.traj=nr.traj,
      tfr_matrix_reconstructed = tfr_reconstructed,
      output.directory=outdir,
      na.index=wpred$na.index,
      mcmc.set=list(meta=meta, mcmc.list=list()),
      nr.projections=this.nr.project,
      burnin=NA, thin=NA,
      end.year=wend.year, use.tfr3=FALSE, burnin3=NA, thin3=NA,
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

