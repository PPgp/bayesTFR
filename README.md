# bayesTFR

[![Travis-CI Build Status](https://travis-ci.org/PPgp/bayesTFR.svg?branch=dev)](https://travis-ci.org/PPgp/bayesTFR)


R package for projecting probabilistic total fertility rate. 

# Annual Version with Uncertainty

The dev branch of the bayesTFR repository implements an annual version of the model with taking uncertainty about the estimation into account. 


## Install:

To install this version, run the following:

```R
library(devtools)
install_github("PPgp/bayesTFR", ref = "dev")
```

Then if bayesTFR has been loaded previously, users could switch to the new version either by re-starting the r session or

```R
detach(package:bayesTFR, unload=TRUE); library(bayesTFR)
```

## Usage

There are two ways to estimate the model:

1. **Two-step estimation**: 
	1. Estimate Phase II using `run.tfr.mcmc`.
	2. Estimate Phase III using `run.tfr3.mcmc`.
2. **One-step estimation**:
	* Phase II and Phase III are estimated jointly via `run.tfr.mcmc`.

Moreover, additional autocorrelation in Phase II can be taken into account using an AR(1) model. This will be denoted as "Phase II AR(1)".

The following table shows the availability of the various options depending on whether an annual version is used and whether uncertainty is taken into account:


|Annual |Uncertainty  ||
|------|:----:|:----:|
| |TRUE|FALSE|
|TRUE|**A**<br/>one-step estimation<br/> Phase-II AR(1) allowed|**B**<br/>two-step estimation<br/> Phase-II AR(1) allowed|
|FALSE|**C**<br/>one-step estimation | **D**<br/>two-step estimation |

Cell **D** represents the original **bayesTFR** (available on CRAN and the master branch), while the other cells represent new functionalities available on this branch.

###  5-year version (Cell D)

Here is a toy example of estimating Phase II:

```R
mcmc.list <- run.tfr.mcmc(output.dir = "bayesTFR.output", iter = 10, nr.chains = 1)
```

Then to estimate Phase III, one can do:

```R
mcmc3.list <- run.tfr3.mcmc(sim.dir = "bayesTFR.output", iter = 10, nr.chains = 1, thin = 1)
```

To generate TFR predictions for this example, do:

```R
pred <- tfr.predict(sim.dir = "bayesTFR.output", burnin = 0, burnin3 = 0)
```

See [Sevcikova et al (2011)](https://www.jstatsoft.org/article/view/v043i01) for more details.

### Annual Version (Cell A and B)

For the rest of this tutorial, we will use a toy example with the following settings:

```R
output.dir <- 'bayesTFR.tutorial'
nr.chains <- 3
total.iter <- 100
annual <- TRUE
burnin <- 0
thin <- 1
```

The **annual** option is implemented in the logical argument `annual` which is passed to the `run.tfr.mcmc` function, regardless if two-step or one-step estimation is used. To be specific, if we want to estimate using annual data, then 

```R
mcmc.list <- run.tfr.mcmc(output.dir = output.dir, iter = total.iter, 
	nr.chains = nr.chains, annual = annual)
```

By default the WPP UN 5-year TFR estimates are linearly interpolated and used for the estimation. If users are interested in using their own annual TFR file, the name of the file is passed into the `my.tfr.file` argument:

```R
mcmc.list <- run.tfr.mcmc(output.dir = output.dir, iter = total.iter, 
	nr.chains = nr.chains, annual = annual, my.tfr.file = "sample_tfr.txt")
```

The `my.tfr.file` should be a tab-separated file of the following format:

|country_code|name|1950|1951|...|2020|
|------|----|----|---------|------------------|-----------|
|4|Afghanistan|7.45|7.45|...|4.176|


For Phase III estimation and for prediction, the `annual` argument is not needed. Thus, the same settings as for the 5-year version can be used.

#### Phase II AR(1) option

Cells A and B can account for additional auto-correlation in Phase II. 
Currently, the model for phase II is
$$ f_{c,t+1} = f_{c,t} - g(f_{c,t}|\theta_c) + \varepsilon_{c,t} $$
where $$g_{c,t} = g(f_{c,t}|\theta_c)$$ is the double-logistic decline.

With five-year version of the data, this is OK. However, with annual version of the data, the auto-correlation of the $\varepsilon_{c,t}$ for phase II is severe. Thus, we introduce a different model:
$$ f_{c,t+1} = f_{c,t} - d_{c,t} \\ d_{c,t+1} - g_{c,t+1} = \phi(d_{c,t} - g_{c,t}) + \varepsilon_{c,t}$$

We introduce this as a logical argument `ar.phase2` of `run.tfr.mcmc`. Users who are interested in using this could call:

```R
mcmc.list.ar <- run.tfr.mcmc(output.dir = output.dir, iter = total.iter, 
	nr.chains = nr.chains, annual = annual, ar.phase2 = TRUE)
```

Note that this setting is inherited by the `run.tfr.mcmc.extra` function, see below. 

### Estimation with Uncertainty (Cell A and C)

To account for uncertainty in the input data, the logical argument `uncertainty` can be used: 

```R
run.tfr.mcmc(output.dir = "bayesTFR.output", iter = 10, 
	uncertainty = TRUE, my.tfr.raw.file = NULL)
```

The package contains a dataset `rawTFR` which is used by default to assess the uncertainty. Users can provide their own data in a file with the same format, see `data(rawTFR)`: 

|country_code|year|tfr|method|source|
|--------|----|---------|------------------|-----------|
|4|1965|7.97|Indirect|Census|

The name of such file should be passed to the `my.tfr.raw.file` argument.

Note that columns **country_code, year, tfr** should always be provided.

If specify **annual=TRUE**, the function will estimate the past annual TFR for countries with uncertainty. Otherwise it will estimate 5-year average.

Note that when uncertainty is taken into consideration, phase III needs to be run simultaneously, and this part has already be included so users will **NOT** need run **run.tfr3.mcmc** again. This is the reason why we call estimation with uncertainty **one step estimation**.

For prediction, and option for uncertainty is provided. If someone wants to use the prediction with past TFR estimation as trajectories, then call:

```R
pred.list <- tfr.predict(sim.dir="bayesTFR.output", burnin=0, burnin3=0, uncertainty=TRUE)
```

Otherwise, it will take the past TFR as the input in my.tfr.file, and if that is not provided it will take the linearly interpolated UN estimates as input (but still use the estimated parameters with uncertainty considered).

#### Special Treatment of VR data
With our observation, some countries, especially for some OECD countries such as the United States, the VR records are very accurate. Since the current UN estimates are 5-year based, the major source of bias estimated for these VR records is the gap between calendar year records and five-year average, which will over-smooth the phase III estimates. Thus, we add another option **iso.unbiased** for this purpose. Users could call:
```R
### Specify with the same setting
run.tfr.mcmc(output.dir = output.dir, nr.chains = nr.chains, iter = total.iter, annual = annual, 
  burnin = burnin, thin = thin, uncertainty = TRUE, iso.unbiased = c(36,40,56,124,203,208,246,250,276,300,352,372,380,392,410,428, 
                                         442,528,554,578,620,724,752,756,792,826,840))
```

The example showed above includes all OECD countries, which may be further analyzed. For countries listed in **iso.unbiased**, the estimated bias of VR records and official estimates is 0, and the standard deviation of VR records and official estimates is 0.0161.

*With the same setting as above, the low uncertainty in phase III prediction is partly resolved.*

#### User-input covariates
We may include raw data sets with different covariates, including continuous covariates. Thus, we add two options, **covariates** for categorical variables and **cont_covariates** for continuous variables. 

|country_code|year|tfr|method|source|lag|
|--------|----|---------|------------------|-----------|------|
|4|1965|7.97|Indirect|Census|3.5|

users could call:
```R
mcmc.list <- run.tfr.mcmc(output.dir = output.dir, nr.chains = nr.chains, iter = total.iter, annual = annual, 
  burnin = burnin, thin = thin, uncertainty = TRUE, my.tfr.raw.file = 'tfr_raw_v1.csv'),
                                         covariates = c('source', 'method'), cont_covariates=c('lag'))
```
Note that if users specified their own covariates, then when using **tfr.estimation.plot** with **plot.raw=TRUE**, users need to specify grouping as one of the covariates. For example:
```R
tfr.estimation.plot(mcmc.list = mcmc.list, country.code=566, grouping='source')
```

**Note: If users want to use the iso.unbiased and covariates together, they need to have one covariate with the name "source" including the source of the data (where we could find VR in categories). Otherwise the iso.unbiased will not work (but the process could continue).**

#### Run with single or a few countries
Since running the complete MCMC process takes long, a simpler version for estimation TFR is available. To do that, take annual version as an example:
```R
output.dir.extra <- "bayesTFR.tutorial.extra"
### Specify with the same setting
run.tfr.mcmc(output.dir = output.dir.extra, nr.chains = nr.chains, iter = total.iter, annual = annual, 
  burnin = burnin, thin = thin, uncertainty = TRUE)

## Finish running MCMC for both phase II and phase III
## Specify countries of interest, for example Nigeria and the United States
countries <- c(566, 840)
run.tfr.mcmc.extra(sim.dir = output.dir.extra, countries = countries, iter= 100, burnin=0, uncertainty=TRUE)
```

With this method, users could get access to TFR estimations easier without updating the hyper parameter estimation (or can use with the previous estimates). Users could provide raw data with the same format as described before. 

It is worth mentioning that if someone want a faster version, the user could run the first step with setting **uncertainty=FALSE**. In the next step, the function **run.tfr.mcmc.extra** could still work. To make it clear, the user could:
```R
run.tfr.mcmc(output.dir = output.dir.extra, nr.chains = nr.chains, iter = total.iter, annual = annual, burnin = burnin, thin = thin, replace.output = TRUE)
run.tfr3.mcmc(sim.dir = output.dir, nr.chains = nr.chains, iter = total.iter, thin = thin, burnin = burnin, replace.output = TRUE)

countries <- c(566, 840)
run.tfr.mcmc.extra(sim.dir = output.dir, countries = countries, iter= 100, burnin=0, uncertainty=TRUE)
```

To make predictions on these countries, one need to **make predictions for all countries first** and have a stored version, then call:
```R
tfr.predict(sim.dir = output.dir, nr.traj = total.iter, burnin = 0, burnin3 = 0)

tfr.predict.extra(sim.dir = output.dir, countries = countries, uncertainty = TRUE)
```

**Note**, if users want to have different raw data files, or different covariates, it is allowed to run:
```R
run.tfr.mcmc.extra(sim.dir = output.dir, countries = countries, iter= 100, burnin=0, uncertainty=TRUE, 
    my.tfr.raw.file='tfr_raw_v1.csv', covariates=c('cov1', 'cov2'), cont_covariates=c('cov3'))

```

Here, similar to **tfr.predict**, the uncertainty parameter is used to control the TFR used in prediction, but not for the MCMC parameters.

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

Moreover, **tfr.estimation.table** is updated for both forecast and estimation uncertainties.

#### Change of plotting forecast with uncertainty
After users make forecast of TFR with uncertainty considered in the past, they could making the plot showing uncertainty for both estimation and projections. This could be done simply by calling:
```R
tfr.pred <- get.tfr.prediction(sim.dir="bayesTFR.output")
tfr.trajectories.plot(tfr.pred=tfr.pred, country =566, nr.traj = 20, uncertainty=TRUE)
```

#### Obtain model estimation for bias and measurement error variance
After running MCMC steps with uncertainty, users could obtain the bias and measurement error standard deviation used in the estimation process with **get.bias.model** and **get.std.model** as follows:
```R
return_list <- get.bias.model(sim.dir = "bayesTFR.output", country.code=566)
# return_list <- get.std.model(sim.dir = "bayesTFR.output", country.code=566)
summary(return_list$model)

head(return_list$table)
```
The function will return a list with two elements. **model** is the estimated linear model for bias or standard deviation. **table** is the bias (or standard deviation) used in the model estimation. If we only used categorical variables, then we only kept those distinct combinations in the table.

Note that for standard deviation, the final estimation is different from the model by a factor of $\sqrt{\frac{2}{\pi}}$, and we also made some adjustments and thus the estimation could be different between the model and the table for some cases. We adjusted those cases with large bias but little standard deviations estimated, which usually because we have so little data in that group. Moreover, if some countries are included in the **iso.unbiased**, then for the VR records, the bias and standard deviation estimated in the model are not used, but the value in the table is used.

### Change in Diagnosis
Use of **tfr.diagnose** is the same as before. However, we have updated for doing diagnosis for TFR estimation. Instead of requiring all parameters to converge, now we only require 95% TFR estimations to converge as it is almost impossible to achieve convergence for all of the parameters together.

### Discussion of Uncertainty and Smoothness
Here we discuss briefly for the users to understand the bias, standard deviations and uncertainty in estimation in our experiments. Since this is not a formal paper, I will discuss about the ideas of how things work for the whole bayesTFR MCMC process.

The model used in BayesTFR with uncertainty is a Bayesian Hierarchical Model. Roughly speaking, we are assuming that the underlying TFR is unknown, and all data are used to provide evidences for obtaining those values. On the other hand, based on the knowledge as well as the historical observation of the total fertility rates, demographers have concluded three phase model which could be used to describe (and predict) the movement of the total fertility rates. This means that in general for most countries where the population is not extremely small, they should follow this pattern, but with different speed, shape and other features.

Therefore, the idea of the model is the balance between data and the model. For each countries, we assume that their past TFR should follow this model, and the data could help us identify the detailed shape and variability (of their TFR). For countries with better data quality as well as more data, the data will tell us more than the model, and vice versa. For example for the United States, the vital registration system is nearly perfect, and thus, we don't need any model to estimate their past TFR. However, for example for Nigeria, all historical data sources were not perfect with possible large bias and measurement error variance, and thus the model will contribute more on the estimation.

Thus, we could conclude: **With better data with little variance, the estimated median will become less smooth, and the estimated model variance would become larger.** This is observed when we set VR records unbiased for most OECD countries. 

It could be inferred that the opposite should also hold, that is **When the estimated model variance is smaller, the estimated median will become smoother, and the data will contribute less in the estimation process.** This is observed when we introduce **ar.phase2=TRUE** in the model. Since we have bring another source of freedom in the model, which based on our analysis could explain about at least 50% uncertainty in phase II transition, the estimated standard deviation of the model goes low, which is universal for all countries (since this parameter is not country-specific, and also don't need). This problem become extreme when TFR is high. Note that the MCMC estimation process is the balance of likelihoods, when the uncertainty of the model goes extremely low (less than 0.01 for standard deviation), the likelihood change for a small change in TFR will be huge for model part, but less for data parts. This will make the estimation trust the model instead of data, unless the data is extremely convincing even more than VR. This is unstable, and it will also bring problems, since we won't expect that the uncertainty of next year's TFR for countries at the starting part of fertility transition, is much smaller than countries which has gone into phase III for a long time. Therefore, **I suggest that when running the MCMC process, we set the sigma0.min=0.04, especially when auto-correlation of the noise in phase II is considered.** For example,
```R
run.tfr.mcmc(output.dir = output.dir, nr.chains = nr.chains, iter = total.iter, annual = annual, 
  burnin = burnin, thin = thin, uncertainty = TRUE, sigma.min=0.04, ar.phase2=TRUE)
```

In summary, the way to control the smoothness of the estimated median, is to balance the standard deviation of the model and the data. If we want a smoother median, we should avoid extreme low standard deviation estimates in the data (*in the current code, the measurement error standard deviation is not allowed to be smaller than half of the absolute value of bias*), and if we want to avoid over-confidence in model, we should set the minimum of the model standard deviation.

Anyone interested in this part could email me for details.


## More to come
- Email prliu@uw.edu for for requirements



