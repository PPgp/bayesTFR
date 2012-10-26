start.test <- function(name) cat('\n<=== Starting test of', name,'====\n')
test.ok <- function(name) cat('\n==== Test of', name, 'OK.===>\n')

test.load.UNtfr <- function(wpp.year=2008) {
	# read the UN TFR input file
	test.name <- 'loading UN TFR file'
	start.test(test.name)
	tfr <- bayesTFR:::read.UNtfr(wpp.year)
	stopifnot(length(dim(tfr$data.object$data))==2)
	stopifnot(dim(tfr$data.object$data)[1] > 150)
	stopifnot(is.element('last.observed', colnames(tfr$data.object$data)))
	stopifnot(length(tfr$data.object$replaced) == 0)
	stopifnot(length(tfr$data.object$added) == 0)
	test.ok(test.name)
}

test.load.UNlocations <- function(wpp.year=2008) {
	test.name <- 'loading WPP location file'
	start.test(test.name)
	tfr <- bayesTFR:::read.UNtfr(wpp.year)
	locs <- bayesTFR:::read.UNlocations(tfr$data.object$data, wpp.year)
	stopifnot(length(dim(locs$loc_data)) == 2)
	stopifnot(all(is.element(c('country_code', 'include_code'), colnames(locs$loc_data))))
	stopifnot(dim(locs$loc_data)[1] > 200)
	stopifnot(all(is.element(intersect(c(0,1,2), locs$loc_data[,'include_code']), c(0,1,2))))
	test.ok(test.name)
}

test.create.tfr.matrix <- function(wpp.year=2008) {
	test.name <- 'creating TFR matrix'
	start.test(test.name)
	tfr <- bayesTFR:::read.UNtfr(wpp.year)
	locs <- bayesTFR:::read.UNlocations(tfr$data.object$data, wpp.year)
	tfr.and.regions <- bayesTFR:::get.TFRmatrix.and.regions(tfr$data.object$data, locs$loc_data, 
												present.year=2009)
	tfr.matrix <- tfr.and.regions$tfr_matrix
	stopifnot(dim(tfr.matrix)[1] == 12)
	stopifnot(rownames(tfr.matrix)[12] == '2008')
	test.ok(test.name)
}

test.run.mcmc.simulation <- function() {
	sim.dir <- tempfile()

	# run MCMC
	test.name <- 'running MCMC'
	start.test(test.name)
	m <- run.tfr.mcmc(iter=5, nr.chains=1, output.dir=sim.dir)
	stopifnot(m$mcmc.list[[1]]$finished.iter == 5)
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 5)
	stopifnot(bayesTFR:::tfr.set.identical(m, get.tfr.mcmc(sim.dir)))
	test.ok(test.name)

	# continue MCMC
	test.name <- 'continuing MCMC'
	start.test(test.name)
	m <- continue.tfr.mcmc(iter=5, output.dir=sim.dir)
	stopifnot(m$mcmc.list[[1]]$finished.iter == 10)
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 10)
	stopifnot(!is.element(900, m$meta$regions$country_code)) # 'World' should not be included
	test.ok(test.name)

	# run MCMC for an aggregation
	test.name <- 'running MCMC for extra areas'
	start.test(test.name)
	data.dir <- file.path(.find.package("bayesTFR"), 'data')
	m <- run.tfr.mcmc.extra(sim.dir=sim.dir, 
			my.tfr.file=file.path(data.dir, 'my_tfr_template.txt'), burnin=0)
	stopifnot(is.element(900, m$meta$regions$country_code)) # 'World' should be included
	test.ok(test.name)

	# run prediction
	test.name <- 'running projections'
	start.test(test.name)
	pred <- tfr.predict(m, burnin=0, rho=NULL, sigmaAR1=NULL, verbose=FALSE)
	spred <- summary(pred)
	stopifnot(spred$nr.traj == 10)
	stopifnot(!is.element(903, pred$mcmc.set$regions$country_code))
	test.ok(test.name)
	
	npred <- dim(pred$tfr_matrix_reconstructed)[2]
	rho <- pred$rho
	sigma <- pred$sigmaAR1
	mu <- pred$mu
	

	# run MCMC for another aggregation
	test.name <- 'running projections on extra areas'
	start.test(test.name)
	m <- run.tfr.mcmc.extra(sim.dir=sim.dir, countries=903, burnin=0)
	# run prediction only for the area 903
	pred <- tfr.predict.extra(sim.dir=sim.dir, verbose=FALSE)
	stopifnot(is.element(903, pred$mcmc.set$meta$regions$country_code))
	stopifnot(dim(pred$tfr_matrix_reconstructed)[2] == npred+1)
	stopifnot(!is.null(bayesTFR:::get.trajectories(pred, 903)$trajectories))
	test.ok(test.name)

	test.name <- 'estimating AR(1) parameters'
	start.test(test.name)
	ar1pars <- get.ar1.parameters(mu, m$meta)
	stopifnot(rho==ar1pars$rho)
	stopifnot(sigma==ar1pars$sigmaAR1)
	test.ok(test.name)

	test.name <- 'shifting the median'
	start.test(test.name)
	projs <- summary(pred, country='Uganda')$projections
	tfr.median.shift(sim.dir, country='Uganda', shift=1.5, from=2051, to=2080)
	shifted.pred <- get.tfr.prediction(sim.dir)
	shifted.projs <- summary(shifted.pred, country='Uganda')$projections
	stopifnot(all(projs[10:15,c(1,3:dim(projs)[2])]+1.5 == shifted.projs[10:15,c(1,3:dim(projs)[2])]))
	stopifnot(all(projs[c(1:9, 16:19),c(1,3:dim(projs)[2])] == shifted.projs[c(1:9, 16:19),c(1,3:dim(projs)[2])]))
	test.ok(test.name)

	test.name <- 'resetting the median'
	start.test(test.name)
	shifted.pred <- tfr.median.shift(sim.dir, country='Uganda', reset = TRUE)
	shifted.projs <- summary(shifted.pred, country='Uganda')$projections
	stopifnot(all(projs[,c(1,3:dim(projs)[2])] == shifted.projs[,c(1,3:dim(projs)[2])]))
	test.ok(test.name)

	test.name <- 'setting the median'
	start.test(test.name)
	expert.values <- c(2.3, 2.4, 2.4)
	shift <- expert.values - pred$quantiles[15, '0.5',4:6] # Uganda has index 15
	mod.pred <- tfr.median.set(sim.dir, country='Uganda', values=expert.values, years=2024)
	mod.projs <- summary(mod.pred, country='Uganda')$projections
	stopifnot(all(mod.projs[4:6, c(1,3:dim(projs)[2])]==projs[4:6, c(1,3:dim(projs)[2])]+shift))
	stopifnot(all(mod.projs[c(1:3,7:19), c(1,3:dim(projs)[2])]==projs[c(1:3,7:19), c(1,3:dim(projs)[2])]))
	test.ok(test.name)
	
	unlink(sim.dir, recursive=TRUE)
}

test.thinned.simulation <- function() {
	sim.dir <- tempfile()
	# run MCMC
	test.name <- 'running thinned MCMC'
	start.test(test.name)
	m <- run.tfr.mcmc(iter=10, nr.chains=2, output.dir=sim.dir, thin=2)
	stopifnot(m$mcmc.list[[1]]$finished.iter == 10)
	stopifnot(m$mcmc.list[[1]]$length == 6)
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 20)
	stopifnot(get.stored.mcmc.length(m$mcmc.list, burnin=4) == 6)
	test.ok(test.name)

	# continue MCMC
	test.name <- 'continuing thinned MCMC'
	start.test(test.name)
	m <- continue.tfr.mcmc(iter=10, output.dir=sim.dir)
	stopifnot(m$mcmc.list[[1]]$finished.iter == 20)
	stopifnot(m$mcmc.list[[1]]$length == 11)
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 40)
	stopifnot(get.stored.mcmc.length(m$mcmc.list, burnin=4) == 16)
	test.ok(test.name)

	# run MCMC for an aggregation
	test.name <- 'running thinned MCMC for extra areas'
	start.test(test.name)
	data.dir <- file.path(.find.package("bayesTFR"), 'data')
	m <- run.tfr.mcmc.extra(sim.dir=sim.dir, my.tfr.file=file.path(data.dir, 'my_tfr_template.txt'), burnin=0)
	stopifnot(is.element(900, m$meta$regions$country_code)) # 'World' should be included
	test.ok(test.name)

	# run prediction
	test.name <- 'running thinned projections'
	start.test(test.name)
	pred <- tfr.predict(m, burnin=5, verbose=FALSE)
	spred <- summary(pred)
	stopifnot(spred$nr.traj == 16) # 2x8
	stopifnot(pred$mcmc.set$mcmc.list[[1]]$finished.iter == 16)
	stopifnot(length(pred$mcmc.set$mcmc.list) == 1)
	npred <- dim(pred$tfr_matrix_reconstructed)[2]
	test.ok(test.name)

	test.name <- 'running thinned projections on extra areas'
	start.test(test.name)
	m <- run.tfr.mcmc.extra(sim.dir=sim.dir, countries=903, burnin=5)
	# run prediction only for the area 903
	pred <- tfr.predict.extra(sim.dir=sim.dir, verbose=FALSE)
	stopifnot(dim(pred$tfr_matrix_reconstructed)[2] == npred+1)
	test.ok(test.name)

	test.name <- 'plotting DL curves with thinned MCMC'
	start.test(test.name)
	filename <- tempfile()
	png(filename=filename)
	DLcurve.plot(m, 903, burnin=5, nr.curves=5)
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)

	test.name <- 'plotting TFR trajectories with thinned MCMC'
	start.test(test.name)
	filename <- tempfile()
	png(filename=filename)
	tfr.trajectories.plot(pred, 900, nr.traj=5, pi=c(70,65))
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)

	test.name <- 'getting parameter traces with thinned MCMC'
	start.test(test.name)
	traces <- get.tfr.parameter.traces(m$mcmc.list, burnin=5, 
					thinning.index=c(1, 3, 10, 15))
	stopifnot(nrow(traces)==4)
	test.ok(test.name)

	test.name <- 'plotting parameter density'
	start.test(test.name)
	filename <- tempfile()
	png(filename=filename)
	tfr.pardensity.plot(m, burnin=10)
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)

	unlink(sim.dir, recursive=TRUE)
}

test.run.mcmc.simulation.auto <- function() {
	sim.dir <- tempfile()
	# run MCMC
	test.name <- 'running auto MCMC'
	start.test(test.name)
	m <- run.tfr.mcmc(iter='auto', output.dir=sim.dir,
					auto.conf=list(iter=10, iter.incr=5, max.loops=3, nr.chains=2, thin=1, burnin=5))
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 40)
	test.ok(test.name)

	test.name <- 'continuing auto MCMC'
	start.test(test.name)
	m <- continue.tfr.mcmc(iter='auto', output.dir=sim.dir, auto.conf=list(max.loops=2))
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 60)
	test.ok(test.name)

	unlink(sim.dir, recursive=TRUE)
}

test.imputation <- function() {
	sim.dir <- tempfile()

	# run MCMC
	test.name <- 'running MCMC with missing values'
	start.test(test.name)
	my.tfr.file <- file.path(.find.package('bayesTFR'), 'data', 'UN2010_with_last_obs.txt')
	m <- run.tfr.mcmc(iter=5, nr.chains=1, output.dir=sim.dir, my.tfr.file=my.tfr.file)
	stopifnot(m$mcmc.list[[1]]$finished.iter == 5)
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 5)
	stopifnot(bayesTFR:::tfr.set.identical(m, get.tfr.mcmc(sim.dir)))
	# some countries are not DL because of the missing data
	# This was the case in UN2008 but not in UN2010
	# stopifnot(length(get.countries.index(m$meta)) != get.nr.countries(m$meta))
	test.ok(test.name)
	
	test.name <- 'running projections with imputation'
	start.test(test.name)
	pred <- tfr.predict(m, burnin=0, verbose=FALSE)
	spred <- summary(pred)
	stopifnot(spred$nr.traj == 5)
	test.ok(test.name)
	
	test.name <- 'plotting imputed TFR trajectories'
	start.test(test.name)
	filename <- tempfile()
	png(filename=filename)
	tfr.trajectories.plot(pred, 'Eritrea', pi=c(90, 54))
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)

	test.name <- 'creating parameter maps with imputation'
	filename <- tempfile()
	tfr.map(pred, device='png', par.name='Triangle_c4', device.args=list(filename=filename))
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
	
	unlink(sim.dir, recursive=TRUE)
}

test.existing.simulation <- function() {
	test.name <- 'retrieving MCMC results'
	start.test(test.name)
	sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
	m <- get.tfr.mcmc(sim.dir, low.memory=FALSE, burnin=25, chain.ids=c(2))
	stopifnot(length(m$mcmc.list)==1)
	stopifnot(dim(m$mcmc.list[[1]]$traces)[1]==5)
	test.ok(test.name)
}

test.DLcurve <- function() {
	test.name <- 'plotting DL curves'
	start.test(test.name)
	sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
	m <- get.tfr.mcmc(sim.dir)
	filename <- tempfile()
	png(filename=filename)
	DLcurve.plot(m, 'Burkina Faso')
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
}

test.TFRtrajectories <- function() {
	test.name <- 'plotting TFR trajectories'
	start.test(test.name)
	sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
	pred <- get.tfr.prediction(sim.dir=sim.dir)
	filename <- tempfile()
	png(filename=filename)
	tfr.trajectories.plot(pred, 'Australia')
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)

	test.name <- 'tabulating TFR trajectories'
	start.test(test.name)
	t <- tfr.trajectories.table(pred, 'Australia', pi=c(90, 80, 70))
	stopifnot(all(dim(t) == c(30, 9)))
	test.ok(test.name)
}

test.plot.all <- function() {
	test.name <- 'plotting TFR trajectories and DL curves for all countries'
	start.test(test.name)
	sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
	pred <- get.tfr.prediction(sim.dir=sim.dir)
	mc <- get.tfr.mcmc(sim.dir)
	dir <- tempdir()
	tfr.trajectories.plot.all(pred, output.dir=dir, main='XXX trajs')
	trajf <- length(list.files(dir, pattern='png$', full.names=FALSE))
	DLcurve.plot.all(mc, output.dir=dir, main='DL XXX', output.type='jpeg')
	dlf <- length(list.files(dir, pattern='jpeg$', full.names=FALSE))
	unlink(dir, recursive=TRUE)
	stopifnot(trajf == get.nr.countries(mc$meta))
	stopifnot(dlf == length(mc$meta$id_DL))
	test.ok(test.name)
}


test.plot.density <- function() {
	test.name <- 'plotting parameter density'
	start.test(test.name)
	sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
	pred <- get.tfr.prediction(sim.dir=sim.dir)
	filename <- tempfile()
	png(filename=filename)
	tfr.pardensity.cs.plot('Ireland', pred)
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
}

test.plot.map <- function() {
	test.name <- 'creating TFR maps'
	start.test(test.name)
	sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
	pred <- get.tfr.prediction(sim.dir=sim.dir)
	filename <- tempfile()
	tfr.map(pred, year=2043, device='png', device.args=list(filename=filename))
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
	
	test.name <- 'creating parameter maps'
	filename <- tempfile()
	tfr.map(pred, year=1985, device='png', par.name='gamma_2', device.args=list(filename=filename))
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
	
	test.name <- 'creating lambda maps'
	filename <- tempfile()
	tfr.map(pred, device='png', par.name='lambda', quantile=0.7, 
			device.args=list(filename=filename), numCats=20, catMethod="pretty")
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)
}
	
test.get.parameter.traces <- function() {
	test.name <- 'getting parameter traces'
	start.test(test.name)
	sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
	m <- get.tfr.mcmc(sim.dir, low.memory=TRUE)
	traces <- get.tfr.parameter.traces(m$mcmc.list, burnin=15, 
					thinning.index=c(4, 25, 29))
	stopifnot(nrow(traces)==3)
	m.check <- get.tfr.mcmc(sim.dir, low.memory=FALSE, burnin=15)
	stopifnot(traces[1,'chi']==m.check$mcmc.list[[1]]$traces[4,'chi'])
	stopifnot(all(traces[c(2,3),'chi']==m.check$mcmc.list[[2]]$traces[c(10,14),'chi']))
	
	traces <- get.tfr.parameter.traces(m$mcmc.list, burnin=15, thin=8)
	stopifnot(nrow(traces)==4)
	stopifnot(traces[2,'psi']==m.check$mcmc.list[[1]]$traces[9,'psi'])
	stopifnot(traces[4,'psi']==m.check$mcmc.list[[2]]$traces[10,'psi'])
	test.ok(test.name)
}

test.median.adjust <- function() {
	test.name <- 'adjusting median'
	start.test(test.name)
	sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
	pred.dir <- tempfile()
	m <- get.tfr.mcmc(sim.dir)
	pred.dir <- tempfile()
	pred <- tfr.predict(m, burnin=0, output.dir=pred.dir, nr.traj=10, save.as.ascii=0, verbose=FALSE)
	country.obj1 <- get.country.object('Kenya', m$meta)
	country.obj2 <- get.country.object('Mongolia', m$meta)
	new.pred <- tfr.median.adjust(pred.dir, countries=c('Kenya', 'Mongolia'))
	orig.median1 <- get.median.from.prediction(new.pred, country.obj1$index, country.obj1$code, adjusted=FALSE)
	new.median1 <- get.median.from.prediction(new.pred, country.obj1$index, country.obj1$code, adjusted=TRUE)
	orig.median2 <- get.median.from.prediction(new.pred, country.obj2$index, country.obj2$code, adjusted=FALSE)
	new.median2 <- get.median.from.prediction(new.pred, country.obj2$index, country.obj2$code, adjusted=TRUE)
	stopifnot(sum(orig.median1==new.median1)==1) # only the present year should be equal
	stopifnot(sum(orig.median2==new.median2)==1)
	test.ok(test.name)

	test.name <- 'plotting adjusted TFR trajectories'
	start.test(test.name)
	filename <- tempfile()
	png(filename=filename)
	tfr.trajectories.plot(new.pred, 'Kenya', adjusted.only=FALSE)
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	test.ok(test.name)

	test.name <- 'writing summary files'
	start.test(test.name)
	write.projection.summary(sim.dir, adjusted=FALSE)
	write.projection.summary(sim.dir, adjusted=TRUE)
	test.ok(test.name)

	unlink(pred.dir)
}

test.estimate.mcmc.with.suppl.data <- function() {
	sim.dir <- tempfile()
   	# run MCMC
    test.name <- 'estimating MCMC using supplemental data'
    start.test(test.name)
    m <- run.tfr.mcmc(nr.chains=1, iter=30, thin=1, output.dir=sim.dir, start.year=1750, seed=1)
    stopifnot(length(m$meta$suppl.data$regions$country_code) == 103)
	stopifnot(all(dim(m$meta$suppl.data$tfr_matrix) == c(40, 103)))
    test.ok(test.name)
        
	# continue MCMC
	test.name <- 'continuing MCMC with supplemental data'
	start.test(test.name)
	m <- continue.tfr.mcmc(iter=10, output.dir=sim.dir)
	stopifnot(m$mcmc.list[[1]]$finished.iter == 40)
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 40)
	stopifnot(!is.element(900, m$meta$regions$country_code)) # 'World' should not be included
	test.ok(test.name)
        
	# run MCMC for an aggregation
	test.name <- 'estimating MCMC for extra areas with supplemental data'
	start.test(test.name)
	data.dir <- file.path(.find.package("bayesTFR"), 'data')
	m <- run.tfr.mcmc.extra(sim.dir=sim.dir, 
                            my.tfr.file=file.path(data.dir, 'my_tfr_template.txt'), burnin=0, verbose=TRUE)
	stopifnot(is.element(900, m$meta$regions$country_code)) # 'World' should be included
	stopifnot(is.element(900, m$meta$suppl.data$regions$country_code))
	test.ok(test.name)
        
	# run prediction
	test.name <- 'running projections for simulation with supplemental data'
	start.test(test.name)
	pred <- tfr.predict(m, burnin=10, verbose=FALSE, save.as.ascii=0, rho=NULL, sigmaAR1=NULL)
	spred <- summary(pred)
	stopifnot(spred$nr.traj == 30)
	stopifnot(!is.element(903, pred$mcmc.set$regions$country_code))
	npred <- dim(pred$tfr_matrix_reconstructed)[2]
	test.ok(test.name)
	unlink(sim.dir, recursive=TRUE)
}
	
