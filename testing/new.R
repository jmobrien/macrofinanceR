## NEW VERSION:

# Load packages ############
require(macrofinanceR)
require(tidyverse)


# CALCULATE_GAUSSIAN_MODE (new code analogous to this file) --------------------

## Load data ###################
rm(list = ls()) ## clear workspace

new_df <-
  read_csv("~/Freelancing/2023_02_replication_analyses/use_case_data.csv")

## Set parameters #################
nvar = 4 ## number of variables
lc_A0 = matrix(TRUE, nvar, nvar) ## no restrictions on A0

# Getting rough even break points for periods:
# nrow(new_df) %/% 4
# new_df$time[c(60, 120, 180)]


# 1962Q1-1979Q3 & Oil crisis and stagflation
# 1979Q4-1982Q4 & Volcker disinflation
# 1983Q1-1989Q4 & Major S\&L crisis defaults
# 1990Q1-2007Q4 & Great Moderation
# 2008Q1-2010Q4 & Financial crisis
# 2011Q1-2019Q4 & Zero Lower Bound, Recovery from crisis
# 2020Q1-2022Q3 & Covid-19 pandemic and war in Ukraine

period_breaks <-
  c("1979q3", "1982q4",  "1989q4", "2007q4",  "2010q4", "2019q4")

# Produces a few warnings about NaN's
# appears to be fine (iterations that will be discarded),
# though this needs better handling

optout <-
  fit_tvv(dat = new_df,
          date_var = time,
          vars = c(fedfunds, lgdp, lpce, vfci), # not required, as this is all vars in data
          ## object loaded from data file
          n_lags = 3,
          n_cores = 4,
          lc_A0 = lc_A0, # no restrictions
          lmdblock = NULL, # default option: all variances change in all periods
          period_breaks = period_breaks,
          v_prior = 0, # default options, weight 0.01 on each row
          # Startdate and enddate not specified, using all data
          verbose = TRUE,
          maxit = 500 # Seems to converge w/in 500 iterations
  )


impulseplots(optout$ir[,,,1],
             filename = './testing/new_output/IR_gaussian_mode_pack',
             nsteps = 60,
             varnames = c("fedfunds",  "lgdp",  "lpce",    "vfci"),
             height = 10, width = 6)


## Save output
save(optout, file = './testing/new_output/gmode.rda')

# Reload for workshopping:
# testdat <- readRDS(file = './testing/new_data/test_gmode.rds')


# CALCULATE_POSTERIOR_DRAWS ----------------------------------------------------

## This script is for running MCMC, starting from
## the Gaussian posterior mode in "output/gmode.Rdata"

## The last part of the script does a full 10x10 plot of impulse responses
## The script plot_impulse_responses.R generates some other interesting plots,
## including the credit data responses with monetary policy maintaining zero
## interest rates

## Karthik Sastry
## January 2018

## This script is structured to recreate the Gaussian, t, and 3 normal mixture
## models with the exact parameters in the paper. With minor modifications,
## it could handle other cases (e.g., t with different degrees of freedom)

model_choices = c('gaussian','t','mix')
my_choice = model_choices[2] ## by default, draw from the t model


## Choose Tuning Parameters  #################

ncore = 4 ## number of cores to use in parallelization. Default is 1.
## Note that Windows machines don't work with the parallel package with ncore > 1

hessian_scaling = 0.05 ## scaling of inverse hessian as covariance matrix
disperse = TRUE ## if TRUE, disperse the starting point using inverse hessian as covariance matrix
## Note that recommended practice would be to run multiple MCMC chains and assess convergence
## by comparing their posterior draws. One could do that by running this script multiple times
## in parallel (bearing in mind whether R is using a separate random seed each time)

## How long to run MCMC Chain
## This is for "proof of concept": Of course increase this for a proper sample!
ndraw = 1000 ## total draws to record
nsep = 1 ## number of draws per recorded draws
## To clarify futher: you have done nsep * ndraw MCMC draws, but only recorded ndraw of them.

savespots = seq(min(100,ndraw),ndraw,length.out=10) ## draws at which to save output
## NOTE: 100 draws is just for testing the code! You want many thousands to get a reliable
## posterior sample. See the Appendix 3 about algorithm convergence for suggestions (and
## convergence diagnostics)

filename = 'testing/new_output/mcmc_out' ## filename in which to save output


## Set MCMC Starting point #################


# If needed, reload data from CALCULATE_POSTERIOR_DRAWS above:
load('testing/new_output/gmode.rda')

if (disperse){
  optout$x = optout$x + mvrnorm(n = 1, mu = rep(0,length(optout$x)), Sigma = optout$opt$H)
}


## Do Gibbs Sampling  ###################


if (my_choice == 'gaussian'){
  mcmc_output = gdraw(optout,
                      ndraw = ndraw,
                      nburn = 0,
                      nsep = nsep,
                      filename = paste0(filename, ".rda"),
                      savespots = savespots,
                      hsnfactor = hessian_scaling,
                      alpha = 1,
                      k = 1,
                      drawa = TRUE, ## if TRUE, save the A_{1:nlag}
                      dscore = TRUE, ## if TRUE, save the log inverse gamma densities, which is useful for MDD calculation
                      drawe = TRUE  ## If TRUE, save the structural residuals
  )
} else if (my_choice == 't'){

  # JMO - Not sure where code for the calculation of this df is,
  #   though for description see p.1852 footnote 11, as well as  appendix
  df = 5.703233415698;                # our preferred (calibrated) choice of df
  scale = 1;                          # normalization; variances are free to adjust

  mcmc_output = gdraw(optout,
                      ndraw = ndraw,
                      nburn = 0,
                      nsep = nsep,
                      filename = paste0(filename, ".rda"),
                      savespots = savespots,
                      hsnfactor = hessian_scaling,
                      tparam  = df / 2,               # This is the right parameter to feed in for t(df)
                      tscale = scale, ## sacle of the t distribution. So the draw is from IG(tparam,tscale^2 * tparam)
                      drawa = TRUE, ## if TRUE, save the A_{1:nlag}
                      dscore = TRUE, ## if TRUE, save the log inverse gamma densities, which is useful for MDD calculation
                      drawe = TRUE  ## If TRUE, save the structural residuals
  )
} else if (my_choice == 'mix'){

  ## Our calibrated choices for the parameters
  params =  c(0.585819818107, 0.675631367346, 0.392593120275, 1.168238423432, 2.550145348080)
  alpha = c(params[c(1,3)],1-sum(params[c(1,3)])) ## probabilities, which sum to 1
  k = params[c(2,4,5)] ## variance scale for each

  mcmc_output = gdraw(optout,
                      ndraw = ndraw,
                      nburn = 0,
                      nsep = nsep,
                      filename = paste0(filename, ".rda"),
                      savespots = savespots,
                      hsnfactor = hessian_scaling,
                      alpha = alpha,
                      k = k,
                      tscale = 1, ## sacle of the t distribution. So the draw is from IG(tparam,tscale^2 * tparam)
                      drawa = TRUE, ## if TRUE, save the A_{1:nlag}
                      dscore = TRUE, ## if TRUE, save the log inverse gamma densities, which is useful for MDD calculation
                      drawe = TRUE  ## If TRUE, save the structural residuals
  )
}


## Checking the MCMC chain ######################


# number of transitions in /recorded/ draws
ntrans = sum(mcmc_output$tout)
# a "sample estimate" of the transition probability (JMO - may be high?)
accept = ntrans/ndraw
# .2 to .3 is a reasonable range to shoot for
## See appendix 3 for more on convergence diagnostics. the CODA package is also very useful
## (see its documentation online).


## Table of relative variances ####################


lambda <-
  array(0,
        c(4, # of variables
          7, # of regimes
          ndraw # of draws
        ))


# JMO - mcmc_output$xout is vector whose elements combining these estimates:
  # n_var * n_var elements (here, 4 * 4 = 16)
  # n_var * (n_regimes - 1) (here, 4 * 6 = 24)
  # total here of 40

# Write gibbs sampling results for regimes into respective spots in the lambda array
lambda[,1:6,] = mcmc_output$xout[17:40,]
# Get normalized variance for final regime column set in array:
lambda[,7,] = 7 - apply(lambda[,1:6,], c(1,3), sum) # retrieving the normalized variance

lambda_median = apply(lambda, c(1,2), median); # taking the median across draws
# can do other posterior stats, credible sets, etc.


## Table of largest shocks (posterior median) #####################


e = mcmc_output$eout # Residuals up to scaling by Lambda

## Deriving right scaling

#JMO - convenience additions:
n_vars <- ncol(optout$dat_ts)
n_obs <- nrow(optout$dat_ts)
n_lags <- optout$n_lags
regime_breaks <- optout$breaks_pos
datevar <- new_df$time

lmdndx = rep(1:7, times = diff(c(regime_breaks, n_obs)))
# JMO - don't fully understand this yet, appears to be scaling all but the first vars
lambda_scaling = sqrt(lambda[,lmdndx[(n_lags + 1):n_obs],])
lambda_scaling = aperm(lambda_scaling,c(2,1,3))
e = e * lambda_scaling

## Median
e_median = apply(e, c(1,2), median)

## Attaching date sequences (JMO - remove 3 lags from beginning of date var)
dseq = datevar[4:length(datevar)]

topN = 4
nShocks = dim(e)[2]

big_shocks = list()
## Making table
for (i in 1:nShocks){
  ## Rank shocks
  myOrder = order(-abs(e_median[,i]))[1:topN]
  ## Put information in table
  shock_df = data.frame(dseq[myOrder],e_median[myOrder,i])
  colnames(shock_df) = c('Date','Shock Value')
  big_shocks[[i]] = shock_df
}



## Get Impulse Responses    ##########################

thin = seq(1,dim(mcmc_output$xout)[2],by=1) ## if you would prefer to do IR for a subsample

## to scale ir to an "average" period, we set the variances in the first period to 1
mcmc_output$xout[17:20,] = 1

ir = McmcIr(t(mcmc_output$xout)[thin,],
            optout, ## this is the posterior mode file. as written, the function extracts the names of the variables and a few other things from here (not actually used in the calculation)
            lrange = 1, ## which variance regimes ("lambdas") for which to report IR. this was useful when I was graphing the IR separately for each regime --- if the goal is to look at an "average" regime with the new rescaling, there's no reason not to just 1
            cores = ncore, ## this will parallelize the sampling of the reduced form coefficients with n cores. If you're already parallelizing the main job, it's probably smarter to use one core here (calculating the IR takes relatively little time either way)
            ##oweights = mcmcout$dout ## these are the extra weights (the variance shock draws) that you need to take account of to get the correct reduced form residuals
            aplus = mcmc_output$aout[,,thin]
)

## warning: this takes up a lot of space on your hard drive!
save(ir,file = paste0(filename,'_ir.rda'))


## Plot all IR  ##########################

vnames = c("FED", "LGDP", "LPCE", "VFCI")

impulseplots(
  ir = apply(ir$ir[,,,1,],c(1:3),median),
  irdraws = ir$ir[,,,1,],
  conf = c(0.68,0.90), ## posterior uncertainty bands
  shocknames = as.character(1:10),
  filename = paste0(filename,'_irplot'),
  varnames = vnames,
  width = 6, height = 8)



## Get MDD #############################

draw_seq = 1:200                        # what draws to use
mout <-
  get_mdd(
    x = t(mcmc_output$x[,draw_seq]),  # draws
    lh = mcmc_output$lh[draw_seq],                # un-normalized posterior
    efac = mcmc_output$dsout[draw_seq],      # density component from non-normal adj parameters
    trunc = 0.95                        # where to truncate Gaussian for GD method
)
## OUTPUT: with and without correction for this GD truncation


# CALCULATE_FORECASTS ----------------------------------------------------------

# JMO - NOT RUN, NEED MORE INFO:

## This script estimates models (posterior modes thereof) using data available up to each date
## t, and then calculates 48-month-out forecasts with each model. It is provided to allow the user
## to tweak model assumptions for the forecasting exercise of the paper. If you just want to re-create
## the plots, we included an .Rdata file with the calculated forecasts that we use in the paper
## and a script "forecasting_plot.R" to generate the relevant plots.

## As written, this script proceeds from the beginning of the sample to the end
## sequentially, using the date t-1 model to ease computation of the date t model.
## You can run separate "chunks" of the date range simultaneously to make the
## calculations a good deal quicker

## Karthik Sastry
## August 2019






## Set tuning parameters #######################


filename = 'testing/new_output/forecast_out'

nvar = 4
lc_A0 = matrix(TRUE, nvar, nvar)

## dates at which to calculate POOS
test_dates = seq(
  from = as.Date('1979/10/01'),
  to = as.Date('2015/05/01'),
  by = '1 month')

## Locations of variance breaks
# 1962Q1-1979Q3 & Oil crisis and stagflation
# 1979Q4-1982Q4 & Volcker disinflation
# 1983Q1-1989Q4 & Major S\&L crisis defaults
# 1990Q1-2007Q4 & Great Moderation
# 2008Q1-2010Q4 & Financial crisis
# 2011Q1-2019Q4 & Zero Lower Bound, Recovery from crisis
# 2020Q1-2022Q3 & Covid-19 pandemic and war in Ukraine

period_breaks <-
  c("1979q3", "1982q4",  "1989q4", "2007q4",  "2010q4", "2019q4")

ndates = length(test_dates)



## Run the estimation and forecasts ######################


## Warning: This takes a very long time!

models_out = list() ## list of the A0,Lambda of each model
fc_out = list() ## list of forecasts
retcode = list() ## return code of the optimization, so you can check if a posterior mode was maybe not found

## seedx = NULL ## could specify starting points for x and H, based on other versions of the model
## seedh = NULL

seed_model = NULL ## starting point
ntsbold = 0 ## this will be used by the loop to determine when to add new variance regime

for (imodel in 1:ndates){

  idate = test_dates[imodel]
  # JMO - (this won't work w/o actual dates)
  ntsb = sum(period_breaks < idate) ## number of variance regimes

  seedx = seed_model$x
  seedh = seed_model$opt$H

  if ((ntsb > ntsbold) & (imodel > 1)){
    ## new variance
    seed_model$lambda = cbind(seed_model$lambda,rep(1,10))
    seedx = c(seedx[1:100],c(seed_model$lambda[,-(dim(seed_model$lambda)[2])]))
    seedh = NULL
  }

  ntsbold = ntsb

  if (idate <= as.Date('1982/12/01')){
    nLags = 2 ## not enough data for 10 lags
  } else {
    nLags = 3
  }

  ## get posterior mode
  seed_model = fit_tvv(
    dat = df_out,
    date_var = time,
    vars = c(fedfunds, lgdp, lpce, vfci), # not required, as this is all vars in data
    n_lags = n_lags,
    lc_A0 = lc_A0, # no restrictions
    lmdblock = NULL, # default option: all variances change in all periods
    period_breaks = period_breaks,
    v_prior = 0, # default options, weight 0.01 on each row
    startdate = as.Date('1973/1/1'),
    enddate = idate,
    verbose = TRUE,
    seedx = seedx,
    seedH = seedh,
    cosprior = NULL,
    critval = 1e-6
  )

  models_out[[imodel]] = seed_model$x ## save only x
  fc_out[[imodel]] = get_forecast(seed_model,48)
  retcode[[imodel]] = seed_model$opt$retcodeh

  ## Save as we go, just in case of an interruption
  save(models_out, file = paste(filename,'_models.Rdata',sep=''))
  save(fc_out, file = paste(filename,'_fc_list.Rdata',sep=''))
  save(retcode, file = paste(filename,'_ret.Rdata',sep=''))

}

fc_array = array(unlist(fc_out),
                 c(dim(fc_out[[1]]),ndates))
save(fc_array,file = paste(filename,'_fc_array.Rdata',sep=''))


# PLOT_IMPULSE_RESPONSE (Figs 1-8) ---------------------------------------------------


## Plots a few different permutations of impulse responses
## All output is saved as pdf files in the plots/ subdiectory

## Karthik Sastry
## January 2018


## Load in an impulse response file ########################


filename = 'testing/new_output/mcmc_out_ir.rda' # if you ran get_posterior_draws.R
# filename = 'testing/new_output/replication_irf_draws.Rdata' # if you want to exactly replicate plots in the draft

irname = 'testing/new_output/model'               # stem for all the titles

load(filename)

if (!is.null(names(ir))){
  ir = ir$ir[,,,1,]                       # easier to read code
}

## Names of all the variables
vnames = c("FED", "LGDP", "LPCE", "VFCI")




## Plot all IR  #############################


## Duplicates the same action at the end of get_posterior_draws.R

impulseplots(
  ir = apply(ir,c(1:3),median),
  irdraws = ir,
  conf = c(0.68,0.90), ## posterior uncertainty bands
  shocknames = as.character(1:10),
  filename = paste(irname,'_irplot',sep=''),
  varnames = vnames,
  width = 6, height = 8)



## 5x5 Blocks (Figures 1-4)  ##########################

# only one 4x4 block here, making it the same as the one just above

blocks = list()
blocks[[1]] = list(x = 1:4,y=1:4)
# blocks[[2]] = list(x = 1:5,y=6:10)
# blocks[[3]] = list(x = 6:10,y=1:5)
# blocks[[4]] = list(x = 6:10,y=6:10)

impulseplots(
  ir = apply(ir,c(1:3),median),
  irdraws = ir,
  shocknames = as.character(1:10),
  blocks = blocks,
  filename = paste(irname,'_block',sep=''),
  varnames = vnames,
  width = 6, height = 6)


## Monetary Policy Shock #########################


blocks = list(list(x=c(1:5),y=6))
blocks[[2]] = list(x=c(6:10),y=6)

impulseplots(
  ir = apply(ir,c(1:3),median),
  irdraws = ir,
  shocknames = NULL,
  blocks = blocks,
  filename = paste(irname,'_mopo',sep=''),
  varnames = vnames,
  width = 3, height = 6)


## Credit Shocks  ###########################


blocks = list(list(x=c(1:5),y=c(3:4)))
blocks[[2]] = list(x=c(6:10),y=c(3:4))

impulseplots(
  ir = apply(ir,c(1:3),median),
  irdraws = ir,
  shocknames = as.character(1:10),
  blocks = blocks,
  filename = paste(irname,'_credit',sep=''),
  varnames = vnames,
  width = 3, height = 6)


## Spread Shocks ###########################


blocks = list(list(x=c(1:5),y=c(9:10)))
blocks[[2]] = list(x=c(6:10),y=c(9:10))
impulseplots(
  ir = apply(ir,c(1:3),median),
  irdraws = ir,
  shocknames = as.character(1:10),
  blocks = blocks,
  filename = paste(irname,'_spread',sep=''),
  varnames = vnames,
  width = 3, height = 6)


## Credit + Monetary Shocks ###########################

ir_small = ir[,c(3,4,6),,]
ndraw = dim(ir_small)[4]

## allocating space
irmix = array(0,c(10,5,60,ndraw)) ## the mixed impulse response

## Constructing the IR mix
nperiod = 60

for (iperiod in 1:nperiod){
  raw_response = ir_small[,1:2,iperiod,] + irmix[,4:5,iperiod,] ## without compensating policy

  ## for shock 3
  policy_weight = -raw_response[6,1,] / ir_small[6,3,1,] ## amount of monetary policy shock to add
  irmix[,4,iperiod,] = raw_response[,1,]
  irmix[,4,iperiod:nperiod,] = irmix[,4,iperiod:nperiod,] + rep(policy_weight,each=10 * (nperiod-iperiod+1)) * ir_small[,3,1:(nperiod-iperiod+1),]

  ## for shock 4
  policy_weight = -raw_response[6,2,] / ir_small[6,3,1,] ## amount of monetary policy shock to add
  irmix[,5,iperiod,] = raw_response[,2,]
  irmix[,5,iperiod:nperiod,] = irmix[,5,iperiod:nperiod,] + rep(policy_weight,each=10 * (nperiod-iperiod+1)) * ir_small[,3,1:(nperiod-iperiod+1),]
}

irmix[,1:3,,] = ir_small

## Plotting this

sn = c('Shock 3 (HHC)','Shock 4 (BC)','MP','3 + MP','4 + MP')

blocks = list()
blocks[[1]] = list(x = c(1:4,6),y=c(1,2,4,5))

impulseplots(
  ir = apply(irmix,c(1:3),median),
  irdraws = irmix,
  shocknames = sn,
  blocks = blocks,
  filename = paste(irname,'_creditMP',sep=''),
  varnames = vnames,
  width = 6, height = 6)



## Variance decomposition ##########################


vd = ir;

for (iv in 1:dim(ir)[4]){
  vd[,,,iv] = vdecomp(ir[,,,iv])
}

blocks = list(list(x=c(1:4),y=c(1:10)))

## Portrait orientation
## impulseplots(
##     ir = apply(vd,c(1:3),median),
##     irdraws = vd,
##     shocknames = as.character(1:10),
##     blocks = blocks,
##     filename = paste(irname,'_vd_portrait',sep=''),
##     varnames = vnames,
##     width = 6, height = 3)

impulseplots(
  ir = apply(vd,c(1:3),median),
  irdraws = vd,
  shocknames = as.character(1:10),
  blocks = blocks,
  filename = paste(irname,'_vd',sep=''),
  varnames = vnames,
  width = 8, height = 5)

# PLOT_FORECASTS (Figs 9-13)----------------------------------------------------------

## Generates two plots in the paper: a comparison of data and rolling forecasts during the crisis,
## and a plot of root mean squared error of different models through the whole sample.

## As written, it loads forecasts from the baseline, no credit, and no spreads models from an included
## file "fc_combined.Rdata." To combine this script with a different set of forecasts created after
## running "forecasting_run.R", change the data loading block at the beginning.



## Loading data and forecasts  #########################


## Example posterior mode, from which I grab "cleaned" (log transformed, correct
## date range) data
load('testing/old_output_replic/postmode.Rdata')

## Forecasts
load('testing/new_output/fc_combined.Rdata')

forecast_full = fc$main                 # all variables model
forecast_nos = fc$nos                   # no spreads
forecast_noc = fc$noc                   # no credit variables


## Plot forecasts in the financial crisis #######################


ydata = optout$listData$Y[-(1:81),] # 1 = October 1979

## What sequence of dates at which to plot forecasts
ds = seq(as.Date('2007/1/1'),by='3 month',length = 16) # sequence of dates to plot
fulldates = seq(from = as.Date('1979/10/01'),          # date sfor ydata
                to = as.Date('2015/05/01'),
                by = '1 month')

ids = rep(0,length(ds))                 # indexes of the forecast dates
for (i in 1:length(ds)){
  ids[i] = which(fulldates == ds[i])
}

padding = c(3,15)                       # how many extra months of data to put on either side
ydata = ydata[(ids[1] - padding[1]):(ids[length(ids)] + padding[2]),]
fulldates = fulldates[(ids[1] - padding[1]):(ids[length(ids)] + padding[2])]


vnames = c('IP','P','HHC','BC','M1','R','PCM','TS','GZ','ES') # short names



## NOTE: This takes an optional argument ymat that specifies axis limits
## I calibrated these "by hand" to include the data plus the forecasts plus some white space.
## If you change the time period for this plots, remember to change the limits too!

plotfc(forecast_full[1:12,,ids],      # full variables model
       ydata, ds,
       vnames, fulldates,
       filename = 'testing/new_output/crisis_forecast_full')

plotfc(forecast_nos[1:12,,ids],       # no spread model
       ydata[,1:7], ds,
       vnames[1:7], fulldates,
       filename = 'testing/new_output/crisis_forecast_nos',
       ymat = c(4.35,4.70,
                4.55,4.65,
                8.2, 8.5,
                7.00,7.40,
                7.20,7.60,
                -0.02,0.06,
                5.6,6.4))

plotfc(forecast_noc[1:12,,ids],       # no credit model
       ydata[,c(1:2,5:10)], ds,
       vnames[c(1:2,5:10)], fulldates,
       filename = 'testing/new_output/crisis_forecast_noc',
       ymat = c(4.35,4.70,
                4.55,4.65,
                7.20,7.60,
                -0.02,0.06,
                5.6,6.4,
                -0.01,0.05,
                0.00,0.12,
                -0.01,0.05))




## Get and plot RMSE #########################


yd = optout$listData$Y[-(1:81),]

myH = c(1,6,12,24,48)                   # what horizons for RMSE
mySeries = c(1:2,5:7)      # what variables to use; this is overlap of
# all three models

rmse_full = getrmse(forecast_full[,c(1:2,5:7),],yd[,mySeries],h=myH)
rmse_nos = getrmse(forecast_nos[,c(1:2,5:7),],yd[,mySeries],h=myH)
rmse_noc = getrmse(forecast_noc[,c(1:5),],yd[,mySeries],h=myH)

plotrmse(c(1979,10),
         rmse_full,
         rmse_nos,
         rmse_noc,
         'testing/new_output/rmse_6mo',
         vname = vnames[mySeries],
         ih = 2)                          # 6-month horizon (2nd one in the list)

plotrmse(c(1979,10),
         rmse_full,
         rmse_nos,
         rmse_noc,
         'testing/new_output/rmse_24mo',
         vname = vnames[mySeries],
         ih = 4)                          # 24-month horizon (2nd one in the list)



# MAKE_TABLE_3 ------------------------------------------------------------



## Plots a few different permutations of impulse responses
## All output is saved as pdf files in the plots/ subdiectory

## Karthik Sastry
## January 2018


rm(list = ls()) ## clear workspace

## setwd('path/to/files')

source('load_functions.R')
load_functions() ## also loads packages



## Load in an impulse response file ############################


## filename = 'testing/new_output/mcmc_out_ir.Rdata' # if you ran get_posterior_draws.R
filename = 'testing/new_output/replication_irf_draws.Rdata' # if you want to exactly replicate plots in the draft

irname = 'testing/new_output/model'               # stem for all the titles

load(filename)
if (!is.null(names(ir))){
  ir = ir$ir[,,,1,]                       # easier to read code
}

## Names of all the variables
vnames = c('IP','P','HHC','BC','M1','R','PCM','TS','GZ','ES')




## Plot all IR  ###########################


## Duplicates the same action at the end of get_posterior_draws.R

impulseplots(
  ir = apply(ir,c(1:3),median),
  irdraws = ir,
  conf = c(0.68,0.90), ## posterior uncertainty bands
  shocknames = as.character(1:10),
  filename = paste(irname,'_irplot',sep=''),
  varnames = vnames,
  width = 6, height = 8)



## 5x5 Blocks (Figures 1-4)  ###########################


blocks = list()
blocks[[1]] = list(x = 1:5,y=1:5)
blocks[[2]] = list(x = 1:5,y=6:10)
blocks[[3]] = list(x = 6:10,y=1:5)
blocks[[4]] = list(x = 6:10,y=6:10)

impulseplots(
  ir = apply(ir,c(1:3),median),
  irdraws = ir,
  shocknames = as.character(1:10),
  blocks = blocks,
  filename = paste(irname,'_block',sep=''),
  varnames = vnames,
  width = 6, height = 6)


## Monetary Policy Shock ##########################


blocks = list(list(x=c(1:5),y=6))
blocks[[2]] = list(x=c(6:10),y=6)

impulseplots(
  ir = apply(ir,c(1:3),median),
  irdraws = ir,
  shocknames = NULL,
  blocks = blocks,
  filename = paste(irname,'_mopo',sep=''),
  varnames = vnames,
  width = 3, height = 6)


## Credit Shocks  ##########################


blocks = list(list(x=c(1:5),y=c(3:4)))
blocks[[2]] = list(x=c(6:10),y=c(3:4))

impulseplots(
  ir = apply(ir,c(1:3),median),
  irdraws = ir,
  shocknames = as.character(1:10),
  blocks = blocks,
  filename = paste(irname,'_credit',sep=''),
  varnames = vnames,
  width = 3, height = 6)


## Spread Shocks ##########################


blocks = list(list(x=c(1:5),y=c(9:10)))
blocks[[2]] = list(x=c(6:10),y=c(9:10))
impulseplots(
  ir = apply(ir,c(1:3),median),
  irdraws = ir,
  shocknames = as.character(1:10),
  blocks = blocks,
  filename = paste(irname,'_spread',sep=''),
  varnames = vnames,
  width = 3, height = 6)


## Credit + Monetary Shocks ############################

ir_small = ir[,c(3,4,6),,]
ndraw = dim(ir_small)[4]

## allocating space
irmix = array(0,c(10,5,60,ndraw)) ## the mixed impulse response

## Constructing the IR mix
nperiod = 60

for (iperiod in 1:nperiod){
  raw_response = ir_small[,1:2,iperiod,] + irmix[,4:5,iperiod,] ## without compensating policy

  ## for shock 3
  policy_weight = -raw_response[6,1,] / ir_small[6,3,1,] ## amount of monetary policy shock to add
  irmix[,4,iperiod,] = raw_response[,1,]
  irmix[,4,iperiod:nperiod,] = irmix[,4,iperiod:nperiod,] + rep(policy_weight,each=10 * (nperiod-iperiod+1)) * ir_small[,3,1:(nperiod-iperiod+1),]

  ## for shock 4
  policy_weight = -raw_response[6,2,] / ir_small[6,3,1,] ## amount of monetary policy shock to add
  irmix[,5,iperiod,] = raw_response[,2,]
  irmix[,5,iperiod:nperiod,] = irmix[,5,iperiod:nperiod,] + rep(policy_weight,each=10 * (nperiod-iperiod+1)) * ir_small[,3,1:(nperiod-iperiod+1),]
}

irmix[,1:3,,] = ir_small

## Plotting this

sn = c('Shock 3 (HHC)','Shock 4 (BC)','MP','3 + MP','4 + MP')

blocks = list()
blocks[[1]] = list(x = c(1:4,6),y=c(1,2,4,5))

impulseplots(
  ir = apply(irmix,c(1:3),median),
  irdraws = irmix,
  shocknames = sn,
  blocks = blocks,
  filename = paste(irname,'_creditMP',sep=''),
  varnames = vnames,
  width = 6, height = 6)



## Variance decomposition #################################################


vd = ir;

for (iv in 1:dim(ir)[4]){
  vd[,,,iv] = vdecomp(ir[,,,iv])
}

blocks = list(list(x=c(1:4),y=c(1:10)))

## Portrait orientation
## impulseplots(
##     ir = apply(vd,c(1:3),median),
##     irdraws = vd,
##     shocknames = as.character(1:10),
##     blocks = blocks,
##     filename = paste(irname,'_vd_portrait',sep=''),
##     varnames = vnames,
##     width = 6, height = 3)

impulseplots(
  ir = apply(vd,c(1:3),median),
  irdraws = vd,
  shocknames = as.character(1:10),
  blocks = blocks,
  filename = paste(irname,'_vd',sep=''),
  varnames = vnames,
  width = 8, height = 5)


# MAKE_TABLE_4 ------------------------------------------------------------

## This script generates Table 4 in the main paper
## using the exact same set of draws from the posterior distribution

## Karthik Sastry
## October 2020


## Constructs epsilon_it and ranks by that

rm(list = ls())                         # Clear workspace

## setwd('path/to/files')

source('load_functions.R')
load_functions() ## also loads packages

load('testing/new_output/replication_draws.Rdata')
ndraw = dim(output$xout)[2]

lambda = array(0,c(10,6,ndraw))
lambda[,1:5,] = output$xout[101:150,]
lambda[,6,] = 6 - apply(lambda[,1:5,], c(1,3), sum) # retrieving the normalized variance

lambda_median = apply(lambda, c(1,2), median); # taking the median across draws
# can do other posterior stats, credible sets, etc.

print(lambda_median)


# MAKE_TABLE_5 ------------------------------------------------------------


## This script generates Table 5 in the main paper
## using the exact same set of draws from the posterior distribution

## NOTE: Because this file does not save the structural residuals,
## (which greatly increase the file size), there is an extra step to
## calculate the residuals. This step is skipped in the corresponding calculation
## in get_posterior_draws.R

## Karthik Sastry
## October 2020


## Constructs epsilon_it and ranks by that

rm(list = ls())                         # Clear workspace

## setwd('path/to/files')

source('load_functions.R')
load_functions() ## also loads packages

load('testing/new_output/replication_draws.Rdata')

## Define function that generates residuals from draws
get_epsilon = function(Aplus,A0,y,nlag = 10){

  ## Forms e_t = A(L)y_t

  ## Aplus is (nvar * nlag + 1 constant) x (neq) matrix

  ## Construct the data with lags
  nv = dim(y)[2] ## number of variables
  T = dim(y)[1]
  Tnl = T - nlag ## time periods, adjusting for lags

  ydata = y[nlag + 1:Tnl,] ## data used in calculations
  ylag = matrix(0,Tnl,nv*nlag + 1)
  for (ilag in 1:nlag){
    ylag[,(ilag-1) * nv + 1:nv] =
      y[nlag - ilag + 1:Tnl,] ## lagged data
  }
  ylag[, nv*nlag + 1] = 1 ## constant row

  ## -A+(L) y_t
  rhs = t(t(Aplus) %*% t(ylag))

  ## A_0 y_t
  lhs = t(A0 %*% t(ydata))

  ## putting it together
  eps_mat = lhs - rhs

  return(eps_mat)
}

## Get the cleaned version of the data
load('testing/new_output/postmode.Rdata')
yout = optout$listData$Y



## Calculate residuals #################################################

xout = output$xout
aout = output$aout

ndraw = dim(xout)[2]
T = 500
eps_out = array(0,c(500,10,ndraw))
eps_scaled = eps_out ## divided by the regime variance
dimpact = matrix(0,10,ndraw)

Tsigbrk = diff(c(0,optout$Tsigbrk[2:6] - 10,500))
Tlist = rep(1:6,times=Tsigbrk)


for (idraw in 1:ndraw){
  eps_out[,,idraw] = get_epsilon(
    Aplus = aout[,,idraw],
    A0 = matrix(xout[1:100,idraw],10,10),
    y = yout,
    nlag = 10)

  dimpact[,idraw] = diag(solve(matrix(xout[1:100,idraw],10,10)))

  lambda = matrix(c(xout[101:150,idraw],rep(0,10)),10,6)
  lambda[,6] = 6 - apply(lambda[,1:5],1,sum)

  eps_scaled[,,idraw] = eps_out[,,idraw] / sqrt(t(lambda[,Tlist]))
}

## For keeping track of dates
year = c(rep(1973,2),
         rep(1974:2014,each=12),
         rep(2015,6))
month = c(11:12,
          rep(1:12,times=41),
          1:6)

## Matrix of eps_out
dmat = matrix(eps_out, 500 * 10, ndraw)
dmed = apply(dmat,1,median)

dmat_scaled = matrix(eps_scaled, 500 * 10, ndraw)
dmed_scaled = apply(dmat_scaled,1,median)

## Identifying the shock and time period
nshock = length(dmed)
shock = rep(1:10,each=500)
period = (1:nshock) - ((shock - 1)*500)
sdate = cbind(month[period],year[period])

## diagonal impacts
median_impact = rep(0,nshock)
for (itop in 1:nshock){
  median_impact[itop] = median(
    dimpact[shock[itop],] *
      dmat[itop,])
}



#### Table with all shocks together


ntop = 10

biggest = order(abs(dmed),decreasing=TRUE)[1:ntop]
bigvals = dmed[biggest]

## table for unscaled values
table = cbind(
  sdate[biggest,],
  shock[biggest],
  round(bigvals,3),
  round(median_impact[biggest],3))

biggest = order(abs(dmed_scaled),decreasing=TRUE)[1:ntop]
bigvals = dmed_scaled[biggest]

table_scaled = cbind(
  sdate[biggest,],
  shock[biggest],
  round(bigvals,3),
  round(2*(1-pt(abs(bigvals),df=5.7))*5000,3)
)

table_combined = cbind(table,table_scaled)



#### Table separated by shocks ( what is printed in the paper )


ntop = 4
tab_list = list()

for (ishock in 1:10){
  biggest = order(abs(dmed * (shock == ishock)),decreasing=TRUE)[1:ntop]
  bigvals = dmed[biggest]

  ## table for unscaled values
  table = cbind(
    sdate[biggest,],
    ## shock[biggest],
    round(bigvals,3),
    round(median_impact[biggest],3))
  tab_list[[ishock]] = table
}


mytab = matrix(0,ntop * 5, 8)

for (ish in 1:5){
  rows = (ish-1)*ntop + 1:ntop
  mytab[rows,1:4] = tab_list[[ish]]
  mytab[rows,5:8] = tab_list[[ish + 5]]
}

print(mytab)


# MAKE_TABLE_6 ------------------------------------------------------------


## A Monte Carlo experiment to see how reproducible are results from smaller regressions
## Karthik Sastry, Jan 2020

## Reproduces results in Section V
## Note that this requires that you have already made the posterior draws

rm(list = ls()) ## clear workspace

## setwd('path/to/files')

source('load_functions.R')
load_functions() ## also loads packages


## Set parameters #################################################


ntrial = 6000 ## number of individual trials

useSeed = TRUE                 # If True, same random seed as main analysis
loadMonteCarlo = TRUE          # If True, don't actually do a new calculation, load saved calculation that replicates the main table

ncore = 1 ## number of cores to use in parallelization. Default is 1.
## Note that Windows machines don't work with the parallel package with ncore > 1

## draws_filename = 'testing/new_output/mcmc_out.Rdata' # Note: use get_posterior_draws.R first to replicate completely
draws_filename = 'testing/new_output/replication_draws.Rdata' # Short-cut for exact replication: same draws we used for main analysis
mode_filename = 'testing/new_output/postmode.Rdata' # Code uses the posterior mode output to grab the data
# (for initial conditions)

load(draws_filename) ## load in a list called 'output' with elements xout, aout
load(mode_filename)

nss = dim(output$xout)[2] ## number of draws total


mysample = seq(1,nss,length.out=ntrial) # the draws we use


## X and A draws that are needed for the simulation
myx = output$xout[,mysample]
mya = output$aout[,,mysample]

## Initial conditions are real data
ic = optout$listData$Y[1:10,]

## Periods for variance regimes
lmdp =  c(diff(optout$Tsigbrk),
          dim(optout$listData$Y)[1] - optout$Tsigbrk[length(optout$Tsigbrk)])
lmdp[1] = lmdp[1] - 10

## Set random seed for exact replication
if (useSeed){
  set.seed(100)
}


if (loadMonteCarlo){
  mc_filename = 'testing/new_output/replication_montecarlo.Rdata'
  load(mc_filename)
  fullr = saveout$fullr
  simout = saveout$simout
} else {
  ## Run simluation
  simout = mclapply(1:ntrial,
                    FUN = \(ii) {

                      myout = tryCatch({

                        ## Get correct posterior draw
                        hpd_model = list()
                        hpd_model$A = matrix(myx[1:100,ii],10,10)

                        hpd_model$vout$var$By = matrix(mya[1:100,,ii],
                                                       10,100,byrow = TRUE)
                        hpd_model$vout$var$By = array(hpd_model$vout$var$By,c(10,10,10))
                        hpd_model$vout$var$Bx = mya[101,,ii]

                        ## lambda = matrix(c(myx[101:140,ii],rep(0,10)),10,5)
                        ## lambda[,5] = 5 - apply(lambda[,1:4],1,sum)

                        lambda = matrix(c(myx[101:150,ii],rep(0,10)),10,6)
                        lambda[,6] = 6 - apply(lambda[,1:5],1,sum)


                        ## generate data
                        simout = sim_model(hpd_model,
                                           lambda = lambda,
                                           lmdp = lmdp,
                                           burn = 0,
                                           y0 = ic,
                                           delta = 5.703233415698/2)

                        mydata = t(simout$sim_out$xout[c(1,2,3,4),])

                        ## quarterly data
                        nq = floor(dim(mydata)[1] / 3)
                        mydata = array(t(mydata),c(4,3,nq))
                        mydata = log(t(apply(exp(mydata),c(1,3),mean)))

                        ## Annual data
                        ny = floor(dim(mydata)[1] / 4)
                        mydata = array(t(mydata),c(4,4,ny))
                        mydata = log(t(apply(exp(mydata),c(1,3),mean)))

                        T = dim(mydata)[1]

                        ## Normalize by lag GDP
                        ngdp = mydata[1:(T-1),1] + mydata[1:(T-1),2]
                        mydata = mydata[2:T,]

                        ## 3Y average regression
                        lm1 = msvreg(mydata[,c(1,3:4)], hi = 3, hd = 3, nlag = 0,
                                     ngdp = ngdp, ratio = TRUE)
                        lm2 = msvreg(mydata[,c(1,3:4)], hi = 3, hd = 3, nlag = 3,
                                     ngdp = ngdp, ratio = TRUE)
                        c1 = (coef(summary(lm1))[2:3,1])
                        c2 = (coef(summary(lm2))[2:3,1])

                        coefs = c(c1,c2)

                        coefs2 = coefs # Normalized coefficients
                        coefs2[c(1,3)] = coefs2[c(1,3)] * lm1$sd[1]
                        coefs2[c(2,4)] = coefs2[c(2,4)] * lm1$sd[2]

                        return(c(coefs,coefs2))
                      }, error = function(error_condition){
                        return(rep(NA,8))
                      })
                    },
                    mc.cores = 1)

  fullr = matrix(unlist(simout),8,ntrial) ## all the coefficients stacked into one matrix

  ## Save outpput
  saveout = list(simout=simout,fullr = fullr)
  save(saveout, file = 'testing/new_output/montecarlo_output.Rdata')
}


## printed results


## NA results suggest an explosive impulse response. We always count these as
## "failures" in our tests
fullr[is.na(fullr)] = 1e99

print('Probability of Negative Coefficient')
print(matrix(apply(fullr[1:4,] < rep(0,4),1,sum) / ntrial,2,2,byrow=TRUE))
print('Probability of Economically Significant Negative Coefficient')
print(matrix(apply(fullr[5:8,] < rep(-0.02,4),1,sum) / ntrial,2,2,byrow=TRUE))

