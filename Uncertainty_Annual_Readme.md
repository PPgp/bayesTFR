# bayesTFR-Annual Version with Uncertainty

This is the beta version for running bayesTFR with uncertainty into consideration and the annual version.

## Install:

Users must install **bayesTFR** with the version **modulize** to test the function of this part. To install this, run the following part of the code:
```R
library(devtools)
install_github("PPgp/bayesTFR", ref = "modulize")
```

Then if bayesTFR has been loaded previously, users could switch to the new version either by re-starting the r session or
```R
detach(package:bayesTFR, unload=TRUE); library(bayesTFR)
```

## Usage
### Original 5-year version
We should keep the same as:
```R
mcmc.list <- run.tfr.mcmc(output.dir="bayesTFR.output", iter=10, nr.chains=1)
```

If we want to run the phase III part, then:
```R
mcmc3.list <- run.tfr3.mcmc(sim.dir="bayesTFR.output", iter=10, nr.chains=1, thin=1)
```

For prediction:
```R
pred.list <- tfr.predict(sim.dir="bayesTFR.output", burnin=0, burnin3=0)
```

The process is the same with the original **bayesTFR** package.

### Annual Version
There is an **annual** option implemented. To be specific, if we want to run annual version of bayesTFR, then 
```R
mcmc.list <- run.tfr.mcmc(output.dir="bayesTFR.output", iter=10, nr.chains=1, annual=TRUE)
```
This will take UN estimates, and interpolate TFR estimates with linear interpolation.

If users are interested in using their own annual TFR files, then they could call:
```R
mcmc.list <- run.tfr.mcmc(output.dir="bayesTFR.output", iter=10, nr.chains=1, annual=TRUE, my.tfr.file="sample_tfr.txt")
```

My.tfr.file should be look like:

|country_code|name|1950|1951|...|2020|
|---------------|--------|----|---------|------------------|-----------|
|4|Afghanistan|7.45|7.45|...|4.176|


Then, for phase III MCMC and prediction, we don't need to specify annual again and could call with the same setting as 5-year version.

### Estimation with Uncertainty
An uncertainty option is provided, and we could call them by:
```R
run.tfr.mcmc(output.dir="bayesTFR.output", iter=10, replace.output=TRUE, uncertainty=TRUE, my.tfr.raw.file="TFR_cleaned_2019.csv")
```

When **uncertainty=TRUE**, the function will by default take the "TFR_cleaned_2019.csv" in the modulize branch. Users are happy to provide their own data with the same format: (and provide the name in my.tfr.raw.file)

|Country.or.area|ISO.code|Year|DataValue|Estimating.Methods|DataProcess|
|---------------|--------|----|---------|------------------|-----------|
|Afghanistan|4|1965|7.97|Indirect|Census|

If specify **annual=TRUE**, the function will estimate the past annual TFR for countries with uncertainty. (** Currently, running with the full MCMC with annual=TRUE and uncertainty=TRUE will take a long time.**) Otherwise it will estimate 5-year average.

Note that when uncertainty is taken into consideration, phase III needs to be run simultaneously, and this part has already be included so users will **NOT** need run **run.tfr3.mcmc** again.

For prediction, and option for uncertainty is provided. If someone wants to use the prediction with past TFR estimation as trajectories, then call:
```R
pred.list <- tfr.predict(sim.dir="bayesTFR.output", burnin=0, burnin3=0, uncertainty=TRUE)
```

Otherwise, it will take the past TFR as the input in my.tfr.file, and if that is not provided it will take the linearly interpolated UN estimates as input (but still use the estimated parameters with uncertainty considered).

## More to come
- Faster version for Annual estimation
- Run with simple countries
- More user-friendly interface
- Email me for for requirements


