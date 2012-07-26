library(bayesTFR)
source('test_functions.R')

test.load.UNtfr()
test.load.UNlocations()
test.create.tfr.matrix()
test.run.mcmc.simulation()
test.existing.simulation()
test.DLcurve()
test.TFRtrajectories()
test.plot.density()
test.get.parameter.traces()
test.median.adjust()
test.plot.map()

# disable the following tests when submitting to CRAN
# to speed-up the checking procedure
test.imputation()
test.thinned.simulation()
test.run.mcmc.simulation.auto()
test.plot.all()
