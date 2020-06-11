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
run.tfr.mcmc(output.dir="bayesTFR.output", iter=10, replace.output=TRUE, uncertainty=TRUE, my.tfr.raw.file=NULL)
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

### Change of plotting forecast with uncertainty
After users make forecast of TFR with uncertainty considered in the past, they could making the plot showing uncertainty for both estimation and projections. This could be done simply by calling:
```R
tfr.pred <- get.tfr.prediction(sim.dir="bayesTFR.output")
tfr.trajectories.plot(tfr.pred=tfr.pred, country =566, nr.traj = 20, uncertainty=TRUE)
```

### Run with single or a few countries
Since running the complete MCMC process takes long, a simpler version for estimation TFR is available. To do that, take annual version as an example:
```R
output.dir <- 'bayesTFR.extra'
nr.chains <- 3
total.iter <- 100
annual <- TRUE
burnin <- 0
thin <- 1
### Specify with the same setting
run.tfr.mcmc(output.dir = output.dir, nr.chains = nr.chains, iter = total.iter, annual = annual, 
  burnin = burnin, thin = thin, uncertainty = TRUE)

## Finish running MCMC for both phase II and phase III
## Specify countries of interest, for example Nigeria and the United States
countries <- c(566, 840)
run.tfr.mcmc.extra(sim.dir = output.dir, countries = countries, iter= 100, burnin=0, uncertainty=TRUE)
```

With this method, users could get access to TFR estimations easier without updating the hyper parameter estimation (or can use with the previous estimates). Users could provide raw data with the same format as described before. 

It is worth mentioning that if someone want a faster version, the user could run the first step with setting **uncertainty=FALSE**. In the next step, the function **run.tfr.mcmc.extra** could still work. To make it clear, the user could:
```R
run.tfr.mcmc(output.dir = output.dir, nr.chains = nr.chains, iter = total.iter, annual = annual, burnin = burnin, thin = thin)
run.tfr3.mcmc(sim.dir = output.dir, nr.chains = nr.chains, iter = total.iter, thin = thin, burnin = burnin)

countries <- c(566, 840)
run.tfr.mcmc.extra(sim.dir = output.dir, countries = countries, iter= 100, burnin=0, uncertainty=TRUE)

```

To make predictions on these countries, one need to **make predictions for all countries first** and have a stored version, then call:
```R
tfr.predict(sim.dir = output.dir, nr.traj = total.iter, burnin = 0, burnin3 = 0)

tfr.predict.extra(sim.dir = output.dir, countries = countries, uncertainty = TRUE)
```

Here, similar to **tfr.predict**, the uncertainty parameter is used to control the TFR used in prediction, but not for the MCMC parameters.

### Special Treatment of VR data
With our observation, some countries, especially for some OECD countries such as the United States, the VR records are very accurate. Since the current UN estimates are 5-year based, the major source of bias estimated for these VR records is the gap between calendar year records and five-year average, which will over-smooth the phase III estimates. Thus, we add another option **iso.unbiased** for this purpose. Users could call:
```R
output.dir <- 'bayesTFR.extra'
nr.chains <- 3
total.iter <- 100
annual <- TRUE
burnin <- 0
thin <- 1
### Specify with the same setting
run.tfr.mcmc(output.dir = output.dir, nr.chains = nr.chains, iter = total.iter, annual = annual, 
  burnin = burnin, thin = thin, uncertainty = TRUE, iso.unbiased = c(36,40,56,124,203,208,246,250,276,300,352,372,380,392,410,428, 
                                         442,528,554,578,620,724,752,756,792,826,840))
```

The example showed above includes all OECD countries, which may be further analyzed. For countries listed in **iso.unbiased**, the estimated bias of VR records and official estimates is 0, and the standard deviation of VR records and official estimates is 0.0161.

*With the same setting as above, the low uncertainty in phase III prediction is partly resolved.*

### User-input covariates
As Patrick requires, we may include raw data sets with different covariates, including continuous covariates. Thus, we add two options, **covariates** for categorical variables and **cont_covariates** for continuous variables. 

|Country.or.area|ISO.code|Year|DataValue|Estimating.Methods|DataProcessShort|TimeLag|
|---------------|--------|----|---------|------------------|-----------|-------|
|Afghanistan|4|1965|7.97|Indirect|Census|3.5|

users could call:
```R
mcmc.list <- run.tfr.mcmc(output.dir = output.dir, nr.chains = nr.chains, iter = total.iter, annual = annual, 
  burnin = burnin, thin = thin, uncertainty = TRUE, iso.unbiased = c(36,40,56,124,203,208,246,250,276,300,352,372,380,392,410,428, 
                                         442,528,554,578,620,724,752,756,792,826,840), my.tfr.raw.file = 'tfr_raw_v1.csv'),
                                         covariates = c('Estimating.Methods', 'DataProcessShort'), cont_covariates=c('TimeLag'))
```
Note that if users specified their own covariates, then when using **tfr.estimation.plot** with **plot.raw=TRUE**, users need to specify grouping as one of the covariates. For example:
```R
tfr.estimation.plot(mcmc.list = mcmc.list, country.code=566, grouping='DataProcessShort')
```

### Auto-correlation for phase II
Currently, the model for phase II is
$$ f_{c,t+1} = f_{c,t} - g(f_{c,t}|\theta_c) + \varepsilon_{c,t} $$
where $$g_{c,t} = g(f_{c,t}|\theta_c)$$ is the double-logistic decline.

With five-year version of the data, this is OK. However, with annual version of the data, the auto-correlation of the $\varepsilon_{c,t}$ for phase II is severe. Thus, we introduce a different model:
$$ f_{c,t+1} = f_{c,t} - d_{c,t} \\ d_{c,t+1} - g_{c,t+1} = \phi(d_{c,t} - g_{c,t}) + \varepsilon_{c,t}$$

We introduce this as an option **ar.phase2** for **run.tfr.mcmc**. Users who are interested in using this could call:
```R
mcmc.list.ar <- run.tfr.mcmc(output.dir = output.dir, nr.chains = nr.chains, iter = total.iter, annual = annual, 
  burnin = burnin, thin = thin, uncertainty = TRUE, iso.unbiased = c(36,40,56,124,203,208,246,250,276,300,352,372,380,392,410,428, 
                                         442,528,554,578,620,724,752,756,792,826,840), ar.phase2=TRUE)

```
Note that if users set **ar.phase2=TRUE**, they will keep it the same when calling **run.tfr.mcmc.extra** at the same directory.

### Discussion of Uncertainty and Smoothness


### Faster Version for bayesTFR with uncertainty
With the new version of the code, users could run **run.tfr.mcmc** at a fast speed. For 3 chains with 62000 iterations each without parallel, it will take at most 1.5 days to finish running. For example, with the following data input with name **tfr_raw_v1.csv**,


## More to come
- Email me for for requirements


