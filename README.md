# bayesTFR

[![R build status](https://github.com/PPgp/bayesTFR/workflows/R-CMD-check/badge.svg?branch=cran)](https://github.com/PPgp/bayesTFR/actions?workflow=R-CMD-check)


R package for projecting probabilistic total fertility rate. See [Documentation on CRAN](https://cran.r-project.org/web/packages/bayesTFR/bayesTFR.pdf). 

For basic usage see [Sevcikova et al (2011)](https://www.jstatsoft.org/article/view/v043i01). Below we describe the package workflow after a major update in version 7.0, which adds the option of taking into account uncertainty about the observed data, as well as estimation and projections for annual data.


# Annual Version with Uncertainty


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

Cell **D** represents the original **bayesTFR** (available on CRAN and on the master branch), while the other cells represent new functionalities available on this branch.

###  5-year version (Cell D)

Here is a toy example of estimating Phase II:

```R
mc2 <- run.tfr.mcmc(output.dir = "bayesTFR.output", iter = 10, nr.chains = 2)
```

Then to estimate Phase III, one can do:

```R
mc3 <- run.tfr3.mcmc(sim.dir = "bayesTFR.output", iter = 10, nr.chains = 2, thin = 1)
```

To generate TFR predictions for this example, do:

```R
pred <- tfr.predict(sim.dir = "bayesTFR.output", burnin = 5, burnin3 = 5)
```

See [Sevcikova et al (2011)](https://www.jstatsoft.org/article/view/v043i01) for more details.

To obtain a previously generated estimation and prediction object, do

```R
mc2 <- get.tfr.mcmc(sim.dir = "bayesTFR.output") # Phase II
mc3 <- get.tfr3.mcmc(sim.dir = "bayesTFR.output") # Phase III
pred <- get.tfr.prediction(sim.dir = "bayesTFR.output") # prediction
```
These three commands can be used for any of the four cells. 

### Annual Version (Cell A and B)

For the rest of this tutorial, we will use a toy example with the following settings:

```R
output.dir <- 'bayesTFR.tutorial'
nr.chains <- 3
total.iter <- 100
thin <- 1
```

The **annual** option is implemented in the logical argument `annual` which is passed to the `run.tfr.mcmc` function, regardless if two-step or one-step estimation is used. To be specific, if we want to estimate using annual data, then 

```R
mc <- run.tfr.mcmc(output.dir = output.dir, iter = total.iter, 
	nr.chains = nr.chains, annual = TRUE)
```

By default the WPP UN 5-year TFR estimates are linearly interpolated and used for the estimation. If users are interested in using their own annual TFR file, the name of the file is passed into the `my.tfr.file` argument:

```R
mc <- run.tfr.mcmc(output.dir = output.dir, iter = total.iter, 
	nr.chains = nr.chains, annual = TRUE, my.tfr.file = "sample_tfr.txt")
```

The `my.tfr.file` should be a tab-separated file of the following format:

|country_code|name|1950|1951|...|2020|
|------|----|----|---------|------------------|-----------|
|4|Afghanistan|7.45|7.45|...|4.176|


For Phase III estimation and for prediction, the `annual` argument is not needed as it is inherited from Phase II. Thus, the same settings as for the 5-year version can be used.

#### Phase II AR(1) option

<!--Currently, the model for phase II is
$$ f_{c,t+1} = f_{c,t} - g(f_{c,t}|\theta_c) + \varepsilon_{c,t} $$
where $$g_{c,t} = g(f_{c,t}|\theta_c)$$ is the double-logistic decline.

With five-year version of the data, this is OK. However, with annual version of the data, the auto-correlation of the $\varepsilon_{c,t}$ for phase II is severe. Thus, we introduce a different model:
$$ f_{c,t+1} = f_{c,t} - d_{c,t} \\ d_{c,t+1} - g_{c,t+1} = \phi(d_{c,t} - g_{c,t}) + \varepsilon_{c,t}$$-->

Using annual data revealed that there is an additional auto-correlation in Phase II that the model does not acount for, since the model was designed for 5-year series. Thus, cells A and B can account for this auto-correlation using a logical argument `ar.phase2` in the `run.tfr.mcmc` function. Users who are interested in using this can call:

```R
mc <- run.tfr.mcmc(output.dir = output.dir, iter = total.iter, 
	nr.chains = nr.chains, annual = TRUE, ar.phase2 = TRUE)
```


### Estimation and Projection with Uncertainty (Cell A and C)

#### One-step Estimation 

To account for uncertainty in the input data (see [Liu & Rafery (2019)](https://projecteuclid.org/euclid.aoas/1593449321)  for details on the method), the logical argument `uncertainty` can be used. In this case, both Phase II and Phase III are estimated jointly and thus, there is no need to invoke `run.tfr3.mcmc`.

```R
mc <- run.tfr.mcmc(output.dir = output.dir, iter = 10, 
	uncertainty = TRUE, my.tfr.raw.file = NULL)
```

The package contains a dataset `rawTFR` which is used by default to assess the uncertainty. It contains the World Fertility Data 2019 provided by the UN Population Division. Users can provide their own data in a file with the same format, see `data(rawTFR)`: 

|country_code|year|tfr|method|source|
|--------|----|---------|------------------|-----------|
|4|1965|7.97|Indirect|Census|

The name of such file should be passed to the `my.tfr.raw.file` argument. Note that columns **country_code, year, tfr** should always be provided.

If it is an annual simulation (i.e. `annual=TRUE`), the function will estimate the past annual TFR, otherwise it will estimate 5-year average.

#### View Estimation Results 
To obtain results with the estimated uncertainty for example for the USA, users could call:

```R
estUS <- get.tfr.estimation(mcmc.list = NULL, country = 840, 
	sim.dir = output.dir, burnin = 0, thin = 1, probs = c(0.05, 0.5, 0.95))
```

Here, if `mcmc.list` is `NULL`, the user should provide the simulation directory in `sim.dir`. `country` is the UN numerical code (USA is 840). It could also be specified as ISO-2 or ISO-3 code ("USA" or "US").

The function returns two elements, `tfr_table` and `tfr_quantile`. The former is a set of estimated TFR trajectories, while the latter is a table of quantiles defined by the argument `probs`.  If `probs` is `NULL`, only trajectories are returned.

In addition, function `tfr.estimation.table` now contains uncertainty for both, the forecast and estimation.

Estimation results can be graphically viewed using e.g.

```R
pl <- tfr.estimation.plot(country = 566, sim.dir = output.dir,
	burnin = 0, thin = thin, pis = c(80, 95), plot.raw = TRUE)
print(pl)
```

The `plot.raw = TRUE` option causes the raw data used in the estimation to be plotted.

Here is different variant of the plotting function which saves the plot into a file:

```R
tfr.estimation.plot(mcmc.list = NULL, country = 840, 
	sim.dir = output.dir, burnin = 0, thin = 1, 
	pis = c(80,95), plot.raw = TRUE, grouping = "DataProcess", 
	save.image = TRUE, plot.dir = 'Estimation.plot')
```

Here, all parameters mentioned in previous functions are as the same. `pis` represents probability interval, where it could be length 1 or 2, and should be values between 50 and 100. `grouping` is for how the colors should be assigned based on groups (could be "DataProcess" or "Estimating.Methods"). `save.image` controls whether the image is saved, and if yes, users could provide a directory where to save it or use the default.


#### Prediction

For prediction, an argument `uncertainty` is provided. If it is set to `TRUE` (which is recommended if the model was also estimated with uncertainty), the projection starts from each estimated trajectory:


```R
pred <- tfr.predict(sim.dir = output.dir, burnin = 0, burnin3 = 0, 
	uncertainty = TRUE)
```

If `uncertainty` is set to `FALSE`, the projection starts from a single point. By default it is the last UN WPP estimate (possibly linearly interpolated if the annual option is used), or data provided in the `my.tfr.file` argument. These are datasets that the estimation function uses as starting values. However, even with uncertainty switched off in the prediction, if the estimation did consider uncertainty, those parameters are still used by the prediction function.

See results using e.g.

```R
tfr.trajectories.plot(pred, country = "Nigeria", uncertainty = TRUE, nr.traj = 20)
```

#### Special Treatment of VR data

There are countries, especially OECD countries such as the United States, for which the Vital Registration (VR) records are very accurate. Since the current UN estimates are 5-year based, the major source of bias estimated for these VR records is the gap between calendar year records and five-year average, which will over-smooth the Phase III estimates. Thus, we added an option `iso.unbiased` which allows to specify which countries have unbiased VR records. Users could call:

```R
mc <- run.tfr.mcmc(output.dir = output.dir, nr.chains = nr.chains, 
	iter = total.iter, annual = TRUE, thin = thin, uncertainty = TRUE, 
	iso.unbiased = c(36,40,56,124,203,208,246,250,276,300,352,372,
		380,392,410,428,442,528,554,578,620,724,752,756,792,826,840))
```

The example includes all OECD countries. For countries listed in `iso.unbiased`, the estimated bias and standard deviation of VR records is 0 and near 0, respectively (in practice the standard deviation is set to 0.0161). VR records are identified as those for which the "source" column in the `my.tfr.raw.file` or `rawTFR` is set to "VR".


#### User-defined input covariates

Users can include different covariates in the raw dataset, including continuous covariates. For example, we can add a lag variable to the dataset:

|country_code|year|tfr|method|source|lag|
|--------|----|---------|------------------|-----------|------|
|4|1965|7.97|Indirect|Census|3.5|


The `run.tfr.mcmc` function has two arguments for this purpose, namely
`covariates` for categorical variables and `cont_covariates` for continuous variables. In this case, if the above dataset is stored in "tfr\_raw\_v1.csv", we can use:


```R
mc <- run.tfr.mcmc(output.dir = output.dir, nr.chains = nr.chains, 
	iter = total.iter, annual = TRUE, thin = thin, uncertainty = TRUE,
	my.tfr.raw.file = 'tfr_raw_v1.csv',
	covariates = c('source', 'method'), cont_covariates = c('lag'))
```

Note that if users specifies their own covariates, then when using `tfr.estimation.plot` with `plot.raw = TRUE`, a grouping of the covariates needs to be specified. For example:

```R
tfr.estimation.plot(mcmc.list = mc, country = 566, grouping = 'source')
```

Note again that if the `iso.unbiased` option is used, the raw dataset should contain the column "source" containing values "VR" for the unbiased records. Alternatively, pass the column name which contains the "VR" entries into the argument `source.col.name` in the `run.tfr.mcmc()` call.

#### Exploring bias and measurement error variance

After running MCMC steps with uncertainty, users could obtain the bias and measurement error standard deviation used in the estimation process using the function `tfr.bias.sd` as follows:

```R
bias_sd <- tfr.bias.sd(sim.dir = output.dir, country = 566)
summary(bias_sd$model_bias)
head(bias_sd$table)
```

The `tfr.bias.sd` function will return a list with elements `model_bias`, `model_sd` and `table`. The `model_bias` and `model_sd` objects are of class `lm` and contain the linear models used to estimate the bias and standard deviation, respectively, while the `table` object includes the observed data points, data quality covariates, and the actual estimates for the specified country, here for Nigeria. If we only use categorical variables, then only distinct combinations are kept in the table.


<!--Note that for standard deviation, the final estimation is different from the model by a factor of $\sqrt{\frac{2}{\pi}}$, and we also made some adjustments and thus the estimation could be different between the model and the table for some cases. We adjusted those cases with large bias but little standard deviations estimated, which usually because we have so little data in that group. Moreover, if some countries are included in the `iso.unbiased`, then for the VR records, the bias and standard deviation estimated in the model are not used, but the value in the table is used.-->



### Run with a single or a few countries

Since running the complete MCMC process takes a long time, an option for experimenting with TFR estimation of one or a few countries is available. However, a full run for all countries has to be available. We strongly recommend to copy the directory with the full run, as the experimenting option overwrites the directory with the new results. 

Let's assume we copied a run for all countries into a directory called "bayesTFR.tutorial.extra". Now we want to experiment with estimation for Nigeria and the USA. We could do:

```R
output.dir.extra <- "bayesTFR.tutorial.extra"
countries <- c(566, 840) # Nigeria and the United States
run.tfr.mcmc.extra(sim.dir = output.dir.extra, countries = countries, 
	iter = total.iter, burnin = 0, uncertainty = TRUE)
```

The function updates only the country-specific parameters of the specified countries, while re-using the existing distribution of the global parameters. Users are allowed to change various estimation settings in this call, e.g. supply raw data (argument `my.tfr.raw.file`) or change covariates (arguments `covariates` and `cont_covariates`). However, several global settings are not subject to change, such as switching between annual and five-year estimation, or changing the `ar.phase2` argument.

It is worth mentioning that one can also change the `uncertainty` setting. For example, for time efficiency reasons one could use a full simulation obtained using `uncertainty = FALSE` while running `run.tfr.mcmc.extra` with `uncertainty = TRUE`. 

To make predictions for these countries, make sure that the main simulation directory contains predictions for all countries (i.e. `tfr.predict` has been processed). Then call:

```R
tfr.predict.extra(sim.dir = output.dir.extra, countries = countries, 
	uncertainty = TRUE)
```

Here, similar to `tfr.predict`, the `uncertainty` argument is used to control if the observed TFR used in prediction are multiple trajectories (`TRUE`) or point estimates (`FALSE`). It does not influence the estimated parameters. 


### Checking Convergence
If estimation is performed with `uncerainty = TRUE`, TFR for each country and time period becomes a parameter to be estimated. This adds over 14200 parameters to the already existing large set of parameters. In such a situation it is almost impossible to achieve convergence for all of the parameters together. Therefore instead of requiring all parameters to converge, the `tfr.diagnose` function now requires that 95% of the TFR parameters to converge. The usage of the `tfr.diagnose` function is the same as in the previous version. If an estimation was obtained via a one-step estimation, the `tfr.diagnose` function checks parameters from Phase II as well as Phase III. 


### Discussion of Uncertainty and Smoothness
Here we discuss briefly for the users to understand the bias, standard deviations and uncertainty in estimation in our experiments. Since this is not a formal paper, I will discuss about the ideas of how things work for the whole bayesTFR MCMC process.

The model used in BayesTFR with uncertainty is a Bayesian Hierarchical Model. Roughly speaking, we are assuming that the underlying TFR is unknown, and all data are used to provide evidences for obtaining those values. On the other hand, based on the knowledge as well as the historical observation of the total fertility rates, demographers have concluded three phase model which could be used to describe (and predict) the movement of the total fertility rates. This means that in general for most countries where the population is not extremely small, they should follow this pattern, but with different speed, shape and other features.

Therefore, the idea of the model is the balance between data and the model. For each country, we assume that its past TFR should follow this model, and the data could help us identify the detailed shape and variability (of its TFR). For countries with better data quality as well as more data, the data will tell us more than the model, and vice versa. For example for the United States, the vital registration system is nearly perfect, and thus, we don't need any model to estimate their past TFR. However, for example for Nigeria, all historical data sources were not perfect with possible large bias and measurement error variance, and thus the model will contribute more to the estimation.

Thus, we could conclude: **With better data with little variance, the estimated median will become less smooth, and the estimated model variance would become larger.** This is observed when we set VR records unbiased for most OECD countries. 

It could be inferred that the opposite should also hold, that is, **when the estimated model variance is smaller, the estimated median will become smoother, and the data will contribute less in the estimation process.** This is observed when we use `ar.phase2=TRUE` in the model. Since we have brought another source of freedom into the model, which based on our analysis could explain about at least 50% uncertainty in the Phase II transition, the estimated standard deviation of the model goes low, which is universal for all countries (since this parameter is not country-specific). This problem becomes extreme when TFR is high. Note that the MCMC estimation process is the balance of likelihoods, when the uncertainty of the model goes extremely low (less than 0.01 for standard deviation), the likelihood change for a small change in TFR will be huge for model part, but less for data parts. This will make the estimation trust the model instead of data, unless the data is extremely convincing even more than VR. This is unstable, and it will also bring problems, since we won't expect that the uncertainty of next year's TFR for countries at the starting part of fertility transition, is much smaller than countries which has gone into Phase III for a long time. Therefore, **I suggest that when running the MCMC process, we set the sigma0.min=0.04, especially when auto-correlation of the noise in Phase II is considered.** For example,

```R
run.tfr.mcmc(output.dir = output.dir, nr.chains = nr.chains, iter = total.iter, 
	annual = TRUE, thin = thin, uncertainty = TRUE, sigma.min = 0.04, ar.phase2 = TRUE)
```

In summary, the way to control the smoothness of the estimated median, is to balance the standard deviation of the model and the data. If we want a smoother median, we should avoid extreme low standard deviation estimates in the data (*in the current code, the measurement error standard deviation is not allowed to be smaller than half of the absolute value of bias*), and if we want to avoid over-confidence in model, we should set the minimum of the model standard deviation.

Anyone interested in this part could email prliu@uw.edu for details.





