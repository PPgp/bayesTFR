library(bayesTFR)
source('test_functions.R')

wpp <- 2010
test.load.UNtfr(wpp)
test.load.UNtfr.and.my.tfr.file()
test.load.UNlocations(wpp)
test.create.tfr.matrix(wpp)
test.existing.simulation()

## disable the following tests when submitting to CRAN
## to speed-up the checking procedure
test.DLcurve()
test.TFRtrajectories()
test.plot.density()
test.get.parameter.traces()
test.median.adjust()
test.plot.map()
test.estimate.mcmc.with.suppl.data()
test.run.mcmc.simulation()
test.run.mcmc.simulation(compression='xz')
test.imputation()
test.thinned.simulation()
test.thinned.simulation(compression='bz')
test.run.mcmc.simulation.auto()
test.plot.all()
