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

### Obtain estimation and plotting
To obtain the estimated uncertainty, users could call:
```R
get.tfr.estimation(mcmc.list=NULL, country.code=840, ISO.code=NULL, sim.dir="bayesTFR.output", burnin=0, thin = 1, probs = c(0.05, 0.5, 0.95))
```

Here, if **mcmc.list** is **NULL**, then the user should provide the simulation directory in **sim.dir** specified for MCMC runs. **country.code** is the UN code (for example USA is 840), and if the user wants to use the ISO.code (say "USA"), then put it in **ISO.code** and leave the **country.code** as **NULL**. If **probs** is **NULL**, only TFR estimated trajectories are selected, otherwise the required quantiles will be obtained.

To make plot easily, users could call:
```R
tfr.estimation.plot(mcmc.list=NULL, country.code=840, ISO.code=NULL, sim.dir="bayesTFR.output", burnin=0, thin = 1, 
  pis = c(80,95), plot.raw=TRUE, grouping="DataProcess", save.image=TRUE, plot.dir = 'Estimation.plot')
```

Here, all parameters mentioned in previous function is as the same. **pis** represents probability interval, where it could be length 1 or 2, and should be values between 50 and 100. **plot.raw** is to control whether the raw data (used for estimation) should be plot, and grouping is for how the colors should be assigned based on groups (could be **"DataProcess"** or **"Estimating.Methods"**). **save.image** will control whether the image is saved, and if yes, users could provide their own plot directory or using the default.

### Run with simple countries
Since running the complete MCMC process takes long, a simpler version for estimation TFR is available. To do that, take annual version as an example:
```R
output.dir <- 'bayesTFR.output'
nr.chains <- 3
total.iter <- 62000
annual <- TRUE
burnin <- 0
thin <- 1
### Specify with the same setting
run.tfr.mcmc(output.dir = output.dir, nr.chains = nr.chains, iter = total.iter, annual = annual, burnin = burnin, thin = thin)
run.tfr3.mcmc(sim.dir = output.dir, nr.chains = nr.chains, iter = total.iter, thin = thin, burnin = burnin)

## Finish running MCMC for both phase II and phase III
## Specify countries of interest, for example Nigeria and the United States
countries <- c(566, 840)
run.tfr.mcmc.extra(sim.dir = output.dir, countries = countries, iter= 50000, burnin=0, uncertainty=TRUE)
```

With this method, users could get access to TFR estimations easier without updating the hyper parameter estimation (or can use with the previous estimates). Users could provide raw data with the same format as described before. 

## More to come
- Faster version for Annual estimation
- More user-friendly interface
- Email me for for requirements


