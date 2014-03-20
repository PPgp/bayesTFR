# 12/13/13: Posterior mean of the ICC for subnational projections. See C8686.
# Here e is a vector of region- and year-specific centered prediction
#   errors divided by their predictive standard deviations,
#   and year is a vector of years

##################################################################
rho.knownvar.bayes <- function (e,year,nrho=100,plot.posterior=FALSE) {
# 12/16/2013 (AER): Function to compute posterior expectation of rho
#   for the intraclass correlation model from data with known means 0 
#   and variances 1.

# Arguments:
#  e=	 vector of region- and year-specific centered prediction
#        errors divided by their predictive standard deviations
#  year= a corresponding vector of years
#  nrho= the length of a vector of values of rho used for computing the
#        posterior mean and mode. Default 100.
#  plot.posterior=TRUE if the posterior is to be plotted.

# Value:
#  postmean = posterior mean of rho
#  postmode = posterior mode of rho
#  Pearson = Pearson correlation coefficient
#  Pearson.known = Pearson correlation coefficient with known mean and
#    variance plugged in instead of empirical ones.

#print(sqrt(mean(e^2)))
#e <- e/sqrt(mean(e^2))

# Compute the likelihood scaling constant 
loglikscale.out <- loglikscalefun (e=e,year=year)
loglikscale <- loglikscale.out$loglikscale

# Pearson correlation coefficient
rhoPearson <- loglikscale.out$rhohat

# Pearson correlation coefficient with plugged-in means and variances
rhoPearson.known <- loglikscale.out$rhohat.known

# Compute posterior mean
rhovec <- seq(from=0,to=1,length.out=nrho+1)[-(nrho+1)]
likvec <- rep(NA,nrho)
for (i in 1:nrho) {
rhoval <- rhovec[i]
likvec[i] <- likrho (rhoval,e=e,year=year,loglikscale=loglikscale)
}
likvec <- likvec - max(likvec)
likvec <- exp(likvec)

if (plot.posterior) {
	dev.new()
	plot (rhovec,likvec,type="l",xlab="rho",ylab="Posterior")
}


# Bayes posterior mean
rho.postmean <- sum (rhovec*likvec)/sum(likvec)

# Bayes posterior mode
rho.postmode <- rhovec [likvec==max(likvec)]
#stop('')
list (postmean=rho.postmean, postmode=rho.postmode, Pearson=rhoPearson,
 Pearson.known=rhoPearson.known)
}

##################################################################
likrho <- function (rho, e=e, year=year, loglikscale=0) {
# likrho: Function to evaluate the likelihood function for the ICC 
#   with given variance.
# Initialization
tmp1 <- 0
years <- sort (unique(year))

# Compute loglikelihood - loglikscale
for (tt in years) {
et <- e[year==tt]
Rt <- length (et)
I <- diag (rep(1,Rt))
J <- matrix (rep(1,Rt^2),ncol=Rt)
Vinv <- (I - (rho/(1+(Rt-1)*rho)) * J) / (1-rho)
logdetV <- (Rt-1)*log(1-rho) + log (1+(Rt-1)*rho)
tmp1 <- tmp1 - (Rt/2)*log(2*pi) - 0.5*logdetV - 0.5 * (t(et) %*% Vinv %*% et)
}
#tmp1 <- tmp1 - loglikscale
#exp (tmp1)
tmp1
}

##################################################################
loglikscalefun <- function (e=e, year=year) {
#likescale: Function to compute the scale for the loglikelihood function.
# This is taken to be the loglikelihood function at a standard MLE of rho.

# First, compute rhohat.
tmp3 <- 0
tmp4 <- 0
tmp5 <- 0
tmp6 <- 0
years <- sort (unique(year))
for (tt in years) {
et <- e[year==tt]
Rt <- length (et)
I <- diag (rep(1,Rt))
J <- matrix (rep(1,Rt^2),ncol=Rt)
tmp3 <- tmp3 + (t(et) %*% J %*% et) - sum(et^2)
tmp4 <- tmp4 + (Rt^2 - Rt)
tmp5 <- tmp5 + sum (et^2)
tmp6 <- tmp6 + Rt}
rhohat <- as.numeric ((tmp3/tmp4) / (tmp5/tmp6) )
rhohat.known <- as.numeric (tmp3/tmp4)

# Second, compute loglik(rhohat)
# Compute loglikelihood - loglikscale
rhohat.bounded <- max (min (rhohat,0.99), .01)
tmp8 <- 0
for (tt in years) {
et <- e[year==tt]
Rt <- length (et)
I <- diag (rep(1,Rt))
J <- matrix (rep(1,Rt^2),ncol=Rt) 

Vinv <- (I - (rhohat.bounded/(1+(Rt-1)*rhohat.bounded))*J) / (1-rhohat.bounded)
logdetV <- (Rt-1)*log(1-rhohat.bounded) + log (1+(Rt-1)*rhohat.bounded)
tmp8 <- tmp8 - (Rt/2)*log(2*pi) - 0.5*logdetV - 0.5 * (t(et) %*% Vinv %*% et)
}
list (loglikscale=tmp8, rhohat=rhohat, rhohat.known=rhohat.known)
}
##################################################################

cor.bayes <- function(errs, is.low, method=c('mean', 'mode'), scale.errors=FALSE, verbose=FALSE, ...) {
	years <- rep(rownames(errs), each=ncol(errs))
	if (scale.errors) errs <- errs/sqrt(mean(errs^2, na.rm=TRUE)) 
	errs.high <- errs.low <- errs
	errs.high[is.low] <- NA
	errs.low[!is.low] <- NA
	errs.list <- list(high=errs.high, low=errs.low)
	errs.all <- as.vector(t(errs))
	keep.all <- which(!is.na(errs.all))
	cor.res <- obs <- cormat <- list()
	use.all <- list(high=FALSE, low=FALSE)
	meth <- paste0('post', method)
	for(what in c('high', 'low')) {
		terrs <- as.vector(t(errs.list[[what]]))
		obs[[what]] <- which(!is.na(terrs))
		cor.res[[what]] <- NULL
		if(length(obs[[what]]) > 0) 
			cor.res[[what]] <- rho.knownvar.bayes(terrs[obs[[what]]], years[obs[[what]]], ...)
		if(is.null(cor.res[[what]]) || is.na(cor.res[[what]][[meth]])) {
			cor.res[[what]] <- rho.knownvar.bayes(errs.all[keep.all], years[keep.all], ...)
			use.all[[what]] <- TRUE
		}
		cormat[[what]] <- matrix(cor.res[[what]][[meth]], nrow=ncol(errs), ncol=ncol(errs))
		diag(cormat[[what]]) <- 1
	}
	if(verbose) {
		cat('\nCor bayes:    high =', round(cor.res[['high']][[meth]],2), ', low =', round(cor.res[['low']][[meth]],2))
		cat('\nUse all pts.: high -', use.all[['high']], ', low -', use.all[['low']])
		cat('\nObservations: high =', length(obs[['high']]), ', low =', length(obs[['low']]), '\n')
	}
	#stop('')
	if(is.na(cor.res[['high']][[meth]])) cormat[['high']] <- cormat[['low']]
	if(is.na(cor.res[['low']][[meth]])) cormat[['low']] <- cormat[['high']]
	return(list(low=cormat[['low']], high=cormat[['high']], corcoef=cor.res$low[[meth]]))
}

cor.bayes.meth10 <- function(errs, is.low, arcsin.prior=FALSE, standardize=TRUE, verbose=FALSE, ...) {
	years <- rep(rownames(errs), each=ncol(errs))
	errs.st <- if(standardize) errs/sqrt(mean(errs^2, na.rm=TRUE))  else errs
	errs.high <- errs.low <- errs.st
	errs.high[is.low] <- NA
	errs.low[!is.low] <- NA
	errs.list <- list(high=errs.high, low=errs.low)
	cor.res <- obs <- cormat <- list()
	ncountries <- ncol(errs)
	ss <- rep(NA, 3)
	for(what in c('high', 'low')) {
		cor.res[[what]] <- matrix(NA, ncol=ncountries, nrow=ncountries)
		err <- errs.list[[what]]
		for(r in 1:(ncountries-1)) {
			if(all(is.na(err[,r]))) next
			for(s in (r+1):ncountries) {
				if(all(is.na(err[,s]))) next
				idx <- !is.na(err[,r]) & !is.na(err[,s])
				ss[1] <- sum(err[idx,r]^2)
				ss[2] <- sum(err[idx,s]^2)
				ss[3] <- sum(err[idx,r] * err[idx,s])
				n <- sum(idx)
				if(arcsin.prior) n <- n+1
				integ <- try(integrate(rhofcn.nom, 0, 1, ss=ss, n=n)$value / integrate(rhofcn.denom, 0, 1, ss=ss, n=n)$value)
				cor.res[[what]][r,s] <- if(!inherits(integ, "try-error")) integ else NA
				cor.res[[what]][s,r] <- cor.res[[what]][r,s] 
			}
		}
		diag(cor.res[[what]]) <- NA
		med <- NA
		if(all(is.na(cor.res[[what]]))) cor.res[[what]] <- NULL
		else {
			m <- mean(cor.res[[what]], na.rm=TRUE)
			if(what == 'low')
				med <- median(cor.res$low, na.rm=TRUE)
			diag(cor.res[[what]]) <- 1
			cor.res[[what]][is.na(cor.res[[what]])] <- m	
		}
	}
	if(verbose) 
		cat('\nCor 10: ', round(med,2))
	if(is.null(cor.res$low) && is.null(cor.res$high)) {
		warning('Correlation method failed. Switching to method 9.')
		return(cor.method9(errs, is.low, verbose=verbose))
	}
	return(list(low=if(is.null(cor.res$low)) cor.res$high else cor.res$low, 
				high=if(is.null(cor.res$high)) cor.res$low else cor.res$high, corcoef=med))
}

rhofcn.nom <- function(rho, ss, n) {
	rho * rhofcn.denom(rho, ss, n)
} 

rhofcn.denom <- function(rho, ss, n) {
	(1/sqrt(1-rho^2))^n * exp(-1/(2*(1-rho^2)) * (ss[1] - 2*rho*ss[3] + ss[2]))
}

cor.moments <- function(errs, is.low, method=c('median', 'mean'), verbose=FALSE) {
##################################################################
# Create two correlation matrices consisting of median or mean of the 
# correlations of supplied errors. Matrix is.low determines which of the 
# errors belong into the 'low' and 'high' correlation matrix, respecitvely.
	method <- match.arg(method)
	corm <- cor(errs, use="pairwise.complete.obs")
	corm[corm<0] <- 0 # truncate
	diag(corm) <- NA
	errs.low <- errs
	errs.low[!is.low] <- NA
	corm.low <- cor(errs.low, use="pairwise.complete.obs")
	corm.low[corm.low<0] <- 0 # truncate
	diag(corm.low) <- NA
	errs.high <- errs
	# This way of doing the high correlation is wrong but ok, because it affects only a few subregions.
	# Right way is to filter out pairs that are both low. 
	errs.high[is.low] <- NA
	corm.high <- cor(errs.high, use="pairwise.complete.obs")
	corm.high[corm.high<0] <- 0 # truncate
	diag(corm.high) <- NA
	avg.cor <- do.call(method, list(corm, na.rm=TRUE))
	avg.cor.high <- do.call(method, list(corm.high, na.rm=TRUE))
	avg.cor.low <- do.call(method, list(corm.low, na.rm=TRUE))
	if(verbose) {
		cat('\nCor ', method, ': high =', round(avg.cor.high,2), ', low =', round(avg.cor.low,2), ', all = ', round(avg.cor,2))
		cat('\nObservations: high =', sum(!is.na(corm.high)), ', low =', sum(!is.na(corm.low)), '\n')
	}
	corm.high[] <- if(is.na(avg.cor.high)) avg.cor else avg.cor.high
	corm.low[] <- if(is.na(avg.cor.low)) avg.cor else avg.cor.low
	diag(corm.high) <- 1
	diag(corm.low) <- 1
	return(list(low=corm.low, high=corm.high, corcoef=if(is.na(avg.cor.low)) avg.cor else avg.cor.low))
}

cor.modified <- function(errs, is.low, verbose=FALSE, ...) {
	errs.low <- errs.high <- errs
	errs.low[!is.low] <- NA
	errs.high[is.low] <- NA
	errs.list <- list(high=list(errs.high, errs), low=list(errs.low, errs.low))
	cor.res <- obs <- cormat <- list()
	nyears <- nrow(errs)
	ncountries <- ncol(errs)
	for(what in c('high', 'low')) {
		counterA <- counterB <- 0
		sumprod <- sumsq <- 0
		for(iyear in 1:nyears) {
			e1 <- errs.list[[what]][[1]][iyear,1:(ncountries-1)]
			if(all(is.na(e1))) next
			for(r in 1:length(e1)) {			
				e2 <- errs.list[[what]][[2]][iyear,(r+1):ncountries]
				if(all(is.na(e2))) next
				p <- e1[r] * e2
				counterA <- counterA + sum(!is.na(p))
				sumprod <- sumprod + sum(p, na.rm=TRUE)
				#stop('')
			}
			err.sq <- errs.list[[what]][[1]][iyear,]^2
			counterB <- counterB + sum(!is.na(err.sq))
			sumsq <- sumsq + sum(err.sq, na.rm=TRUE)
		}
		A <- sumprod/counterA
		B <- sumsq/counterB
		cor.res[[what]] <- A/B
		cormat[[what]] <- matrix(cor.res[[what]], nrow=ncountries, ncol=ncountries)
		diag(cormat[[what]]) <- 1
		obs[[what]] <- c(counterA, counterB)
	}
	if(verbose) {
		cat('\nCor meth7: high =', round(cor.res[['high']],2), ', low =', round(cor.res[['low']],2))
		cat('\nObservations: high =', obs[['high']], ', low =', obs[['low']], '\n')
	}
	if(is.na(cor.res[['high']])) cormat[['high']] <- cormat[['low']]
	if(is.na(cor.res[['low']])) cormat[['low']] <- cormat[['high']]

	return(list(low=cormat[['low']], high=cormat[['high']], corcoef=cor.res$low))
}

cor.moments2 <- function(errs, is.low, verbose=FALSE) {
##################################################################
# Method 8
	errs.low <- errs.high <- errs
	errs.low[!is.low] <- NA
	errs.high[is.low] <- NA
	errs.list <- list(high=errs.high, low=errs.low)
	corm <- list()
	for(what in c('high', 'low')) {
		cormat <- cor(errs.list[[what]], use="pairwise.complete.obs")
		corm[[what]] <- .get.cor.meth8(cormat, errs.list[[what]], verbose=verbose && what=='low')
	}
	return(list(low=if(is.null(corm$low$corm)) corm$high$corm else corm$low$corm, 
				high=if(is.null(corm$high$corm)) corm$low$corm else corm$high$corm, corcoef=corm$low$corcoef))
}

cor.constant <- function(errs, constant=0.99, verbose=FALSE) {
##################################################################
# Method 12
	ncountries <- ncol(errs)
	cormat <- list()
	cormat[['high']] <- cormat[['low']] <- matrix(constant, nrow=ncountries, ncol=ncountries)
	diag(cormat[['high']]) <- 1
	diag(cormat[['low']]) <- 1
	return(list(low=cormat[['low']], high=cormat[['high']], corcoef=constant))
}

.get.cor.meth8 <- function(corm, errs, verbose=FALSE) {
	corm[corm<0] <- 0 # truncate
	diag(corm) <- 1
	Ti <- mean(apply(errs, 2, function(x) sum(!is.na(x))))
	cormat <- ((Ti-1)/Ti)*corm + 1/(2*Ti)
	diag(cormat) <- 1
	cormat.na <- which(apply(is.na(errs), 2, sum) > dim(errs)[1]-2)	
	cormat.no.na <- if(length(cormat.na)==0) cormat else cormat[-cormat.na, -cormat.na]
	if(length(cormat.no.na)<=1) return(list(corm=NULL))
	if(any(is.na(cormat.no.na))) { # replace remaining NAs with averages
		diag(cormat.no.na) <- NA
		avg <- mean(cormat.no.na, na.rm=TRUE)
		cormat.no.na[is.na(cormat.no.na)] <- avg
		diag(cormat.no.na) <- 1
	}
	if(!is.cor.positive.definite(cormat.no.na)) 
        cormat.no.na <- zero.neg.evals(cormat.no.na)   
    if(length(cormat.na)==0) cormat<-cormat.no.na else cormat[-cormat.na, -cormat.na] <- cormat.no.na
    diag(cormat) <- NA
	cor.mom <- list()
	for(method in c('mean', 'median')) {
		cor.mom[[method]] <- do.call(method, list(cormat, na.rm=TRUE))
		if(verbose) 
			cat('\nCor 8,9 ', method, ': ', round(cor.mom[[method]],2))
	}
	if(verbose) cat('\n')
	diag(cormat) <- 1
	cormat[is.na(cormat)] <- cor.mom$mean
	return(list(corm=cormat, corcoef=cor.mom$median))
}

cor.method9 <- function(errs, is.low, verbose=FALSE) {
# Method 9
	errs.low <- errs.high <- errs
	errs.low[!is.low] <- NA
	errs.high[is.low] <- NA
	errs.list <- list(high=errs.high, low=errs.low)
	nyears <- nrow(errs)
	ncountries <- ncol(errs)
	corm <- list()
	for(what in c('high', 'low')) {
		err <- errs.list[[what]]
		cormat <- matrix(NA, ncol=ncountries, nrow=ncountries)
		nom <- denom1 <- denom2 <- counter <- matrix(0, ncol=ncountries, nrow=ncountries)
		for(iyear in 1:nyears) {
			e1 <- err[iyear,1:(ncountries-1)]
			if(all(is.na(e1))) next
			for(r in 1:length(e1)) {
				if(is.na(e1[r])) next			
				e2 <- err[iyear,(r+1):ncountries]
				if(all(is.na(e2))) next
				not.na <- which(!is.na(e2))
				idx <- ((r+1):ncountries)[not.na]
				nom[idx,r] <- nom[r, idx] <- nom[idx,r] + e1[r] * e2[not.na]
				denom1[idx,r] <- denom1[r, idx] <- denom1[idx,r] + rep(e1[r]^2, length(not.na))
				denom2[idx,r] <- denom2[r, idx] <- denom2[idx,r] + e2[not.na]^2
				counter[idx,r] <- counter[r, idx] <- counter[idx,r] + 1
			}
		}
		has.data <- counter > 0
		cormat[has.data] <- (nom[has.data]/counter[has.data])/(sqrt(denom1[has.data]/counter[has.data])*sqrt(denom2[has.data]/counter[has.data]))
		corm[[what]] <- .get.cor.meth8(cormat, err, verbose=verbose && what=='low')
	}
	return(list(low=if(is.null(corm$low$corm)) corm$high$corm else corm$low$corm, 
				high=if(is.null(corm$high$corm)) corm$low$corm else corm$high$corm, corcoef=corm$low$corcoef))
}