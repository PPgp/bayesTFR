library(bayesTFR)
source('test_functions.R')

cran <- FALSE

wpp <- 2019
test.load.UNtfr(wpp)
test.load.UNlocations(wpp)
test.create.tfr.matrix(wpp)
test.load.UNtfr.and.my.tfr.file()
test.existing.simulation()

## disable the following tests when submitting to CRAN
## to speed-up the checking procedure
if(!cran) {
	test.DLcurve()
	test.TFRtrajectories()
	test.plot.density()
	test.get.parameter.traces()
	test.median.adjust()
	test.plot.map()
	test.estimate.mcmc.with.suppl.data()
	test.run.mcmc.simulation(compression='xz')
	test.imputation()
	test.thinned.simulation(compression='bz')
	test.run.mcmc.simulation.auto()
	test.run.mcmc.simulation.auto.parallel()
	test.plot.all()
	test.reproduce.simulation()
	test.subnational.predictions()
	for (wpp in c(2019, 2022)){
	    test.run.mcmc.simulation(wpp.year = wpp)
	    test.run.mcmc.simulation.with.uncertainty(wpp.year = wpp)
	    test.thinned.simulation(wpp.year = wpp)
	    test.run.annual.simulation(wpp.year = wpp)
	}
	for(wpp in rev(c(2010, 2012, 2015, 2017, 2022))) { # these are either suggested packages or not on CRAN
	    test.load.UNtfr(wpp)
	    test.load.UNlocations(wpp)
	    test.create.tfr.matrix(wpp)
	}
}