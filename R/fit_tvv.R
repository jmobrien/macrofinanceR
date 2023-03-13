#' Finding maximum posterior density model w/single structural matrix
#' (inverse(A0)) and time-varying variances for structural innovations
#' (exp(-lmd[,timePeriod]))
#'
#' @param dat Matrix or data frame. Input data
#' @param date_var Variable in data representing date term.
#' @param vars Vector of variables from *dat* to use in
#'   analysis.  Defaults to all variables.
#' @param vars_log Vector of variables to log transform in
#'   analysis. Defaults to all variables untransformed.
#' @param startdate,enddate
#' @param timedummy formatted like start/endDate, but currently not used, could
#'   input a dummy for one/more period,
#' @param period_breaks Date or character vector indicating when periods end.
#'   Should have length of ([# of periods to be specified] - 1).
#' @param n_lags Integer (or integer-like). Number of lags.
#' @param lc_A0 An *n* x *n* matrix, where *n* is the number of variables used.
#'   Each element indicates TRUE for a fitted value and FALSE for a non-fitted
#'   value. If NULL (default) uses lower triangular with TRUE on diagonal, and
#'   zero to other FALSE fields.
#' @param lmdblock n_vars x nPeriods x n_blocks array that identifies, with TRUE,
#'   blocks of parameters that are constant (presumably across periods). The
#'   most obvious application is keeping certain parameters constant over
#'   periods, but it could also be used to have one value in 1/2 the periods,
#'   another in the other half, etc. If left null, all variances change in every
#'   period.specifies exactly how to restrict some variances not to be
#'   time-varying. Not used in the most recent calculations.
#' @param verbose Boolean. if TRUE (default), gives additional output during run
#' @param seed_model a list that looks like the output of this function from
#'   which a first iteration x vector can be extracted intended usage: if you
#'   tweak the prior a bit, the mode will be close to old mode downside:
#'   estiamte of the inverse hessian from optimization will be poor
#' @param n_step
#' @param extracheck
#' @param n_cores Length-1 Integer (or integer-like numeric). Number of cores used for parallelization
#' @param seed_x
#' @param lmd_Prior
#' @param mn_start
#' @param mn_tight
#' @param mn_decay
#' @param v_prior
#' @param cos_prior
#' @param nlt
#' @param hparam_nl
#' @param hparam
#' @param critval
#' @param seed_H
#' @param tvA
#' @param noLmd
#'
#' @details Automatically 'cleans' but does not load the data. Need to run
#'   'datacat' or something similar to grab data from different
#'   sources/spreadsheets and format appropriately

#'
#' @return If `verbose = TRUE`: list output of the optimization routine (in
#'   $opt), plus information about the VAR estimated at the posterior mode. The
#'   whole structure is an input for the programs that run MCMC
#' @export
#'
fit_tvv <-
  function(
    dat,
    date_var,
    vars = NULL,
    vars_log = NULL,
    startdate = NULL, enddate = NULL, timedummy = NULL,
    period_breaks = NULL,
    n_lags = 10,
    lc_A0 = NULL,
    lmdblock = NULL,
    verbose = FALSE,
    seed_model = NULL,
    n_step = 60,
    extracheck = TRUE,
    n_cores = 1,
    seed_x = NULL, lmdPrior = FALSE,
    mn_start = 1, mn_tight = 3, mn_decay = 0.5,
    v_prior = 0, cos_prior = NULL,
    nlt = NULL,
    hparam_nl = rep(0,3), hparam = rep(0,7),
    critval = 1e-10,
    seed_H = NULL,
    tvA = FALSE,
    noLmd = FALSE
  )

  {


    # Data prep (previously dataprep() function) ----

    ## Get date variable ----
    date_q <- rlang::enquo(date_var)
    date_pos <- tidyselect::eval_select(date_q, dat)
    datevar_orig <- dat[[date_pos]]

    ## Make data subset ----
    if(!missing(vars)){
      vars_q <- rlang::enquo(vars)
      vars_pos <- tidyselect::eval_select(vars_q, dat)
      dat_touse <- dat[vars_pos]
    } else {
      # if no selection just remove the date variable:
      dat_touse <- dat[-date_pos]
    }

    ## Log transforms ----
    if(!missing(vars_log)){
      log_q <- rlang::enquo(vars_log)
      # NB This works on the already subsetted data, so positions match:
      log_pos <- tidyselect::eval_select(log_q, dat_touse)
      dat_touse[log_pos] <- log(dat_touse[log_pos]) # needs checks for introduction of NA's

      # Names for output later:
      log_vars <- names(log_pos)
    } else {
      log_vars <- NULL
    }

    ## Start/end dates----

    if(!missing(startdate)){
      # Needs checks - no duplicates, matches, compatible types, etc.
      begin_pos <- which(datevar_orig == startdate)
    } else {
      begin_pos <- purrr::detect(dat_touse, complete.cases)
    }

    if(!missing(enddate)){
      # Needs checks - no duplicates, matches, compatible types, etc.
      end_pos <- which(datevar_orig == enddate)
    } else {
      end_pos <- length(datevar_orig)
    }

    # Slice down to just the time-range specified:
    dat_touse <- dat_touse[begin_pos:end_pos,]
    datevar <- datevar_orig[begin_pos:end_pos]
    ## Time dummy codes (inactive) ----

    ## JMO - taking this out for now until rework is understood:

    # if (!is.null(timedummy)){
    #
    #   n_dummies <- length(timedummy)
    #   dummyMatrix <- matrix(0, nrow = length(dates), ncol = n_dummies)
    #
    #
    #   for (i_dum in (1:n_dummies)) {
    #     if (grepl(':',timedummy[i_dum])){
    #       ##dummy for time range
    #       time_range <- strsplit(timedummy[i_dum],':')
    #       start_point <- which(dates == time_range[1], ind = TRUE)
    #       end_point <- which(dates == time_range[2], ind = TRUE)
    #       dummyMatrix[startPoint:endPoint, i_dum] <- 1
    #     } else {
    #       ##dummy for single period
    #       dummyMatrix[dates == endpoints, i_dum] <- 1
    #     }
    #   }
    #   colnames(dummyMatrix) <- timedummy
    #   X <- ts(dummyMatrix)
    #
    # ## Empty X if time dummies not provided:
    # } else {
    #   X <- NULL
    # }

    ## Breaks for periods ----

    ## Index of breaks position based on breaks date input (formerly breaks_pos):
    ## JMO - needs checks:
    breaks_pos <- which(datevar %in% period_breaks)


    ## MAIN DATA - construct timeseries ----
    dat_ts <- ts(dat_touse)


    # Initial optimization params ----

    # # of variables, variable names, & # of time points (total - lags):
    n_vars <- length(dat_touse)
    var_names <- names(dat_touse)
    n_obs <- nrow(dat_touse) - n_lags


    ##Current method: Estimate a full-sample model and use linear projections
    ##to get individual lmd

    ##Restrictions on A0, lmd by brute force (seed model has no restrictions)
    # Default A0 matrix if needed:
    if (missing(lc_A0)) {
      ##lower triangular default structure:
      lc_A0 <- matrix(FALSE, nrow = n_vars, ncol = n_vars)
      lc_A0[lower.tri(lc_A0)] <- TRUE
    }


    # Setting up prior parameters (formerly pparams()) ----

    # JMO - previously unaddressed arguments in pparams
    a_diag = 1
    ur_lambda <- 5
    ur_mu <- 1

    # JMO - Needs better checks for both mn_decay and mn_tight
    if (!is.null(mn_tight)){
      #### if mn prior is specified
      mn_prior <- list(tight = mn_tight, decay = mn_decay)
    } else {
      mn_prior <- NULL #### no mn_prior
    }

    # JMO - need to review this as a check:
    ## cos_prior just passes through
    ## specify one or the other, probably not both

    ## JMO - not clear what this is all about, review:
    ##these options are grandfathered in, from a time before we thought it was necessary
    ##to tweak v_prior. Else v_prior is just specified as a vector
    ##Could raise an error if length v_prior != n_varss
    if (length(v_prior) == 1 && v_prior == 0){
      v_prior_sig <- rep(.01, n_vars)
    } else if (length(v_prior) == 1 && v_prior == 1){
      ##quick fix: v_prior needs to recognize that some variables are percentage points
      v_prior_sig <- c(.01,.01,.01,1,.01,1,1,1)
      ## three log, one pct pt rate, one log, 3 pct pt rate
    } else {
      ##if full v_prior is providedL:
      v_prior_sig <- v_prior
    }

    # Add variable names:
    names(v_prior_sig) <- var_names

    # Build v_prior
    v_prior <- list(sig = v_prior_sig, w = 0)

    # asig, asd, & urprior
    a_sig <- 2
    a_sd <- outer(v_prior$sig, 1/v_prior$sig)
    ur_prior <- list(lambda = ur_lambda, mu = ur_mu)

    # JMO - these were in pprior, but weren't returned by that function?
    # sigfix <- diag(v_prior$sig^2)
    # dimnames(sigfix) <- list(var_names, var_names)

    # Construct list of params (might drop later unless pparam reintroduced):
    prior_params <-
      list(
        mn_prior = mn_prior, v_prior = v_prior, cos_prior = cos_prior,
        a_sig = a_sig,
        a_sd = a_sd,
        a_diag = a_diag,
        ur_prior = ur_prior,
        mn_start = mn_start,
        cos_prior = cos_prior
      )

    # Build model seed ----

    ## Use previous output of routine
    if (!is.null(seed_model)){

      ##program will re-apply zero restrictions as appropriate, but number of variables
      ##must be correct
      seed_lmd <- seed_model$lambda ## not lmd for this model
      seed_A <- seed_model$A

      n_regimes <- dim(seed_lmd)[2]

      ### New seed not provided, generate
    } else if (is.null(seed_x)){

      ##project a reduced form with fixed variances onto the format we want
      if (any(hparam_nl > 0)){
        ## apply the nonlinear transformation via helper function:
        dat_seed <- transform_nl(dat_ts,nlt)$data
      } else{
        dat_seed <- dat_ts
      }

      ## JMO - needs review?
      ## Hack to deal with dimensionality problem
      if ((n_lags * n_vars) > (n_obs - n_lags)){
        n_lags_rf <- 4 ## crude
      } else {
        n_lags_rf <- n_lags
      }

      ## Call rfvar3 to build reduced-form model:
      seed_rfmodel <-
        rfvar3(
          ydata = dat_seed,
          lags = n_lags_rf,
          lambda = prior_params$ur_prior$lambda,
          mu = prior_params$ur_prior$mu
        )


      temp_A_inv <-
        (t(chol(crossprod(seed_rfmodel$u) / dim(seed_rfmodel$u)[1]))) ##unscaled

      u <- seed_rfmodel$u
      regimes <- c(0, breaks_pos - n_lags_rf, n_obs)
      n_regimes <- length(regimes) - 1
      seed_lmd <- matrix(1, n_vars, n_regimes)

      seed_A <- solve(temp_A_inv)

      seed_A[!lc_A0] <- 0 ##restrictions
      seed_A[upper.tri(seed_A)] <- 0 ##avoiding floating point junk; unclear if important

      seed_A_inv <- solve(seed_A)
      seed_A_inv[upper.tri(seed_A_inv)] <- 0

      ## if (n_regimes == 1) { ##one regime, no need for below
      ## 	seed_lmd <- fullLmd
      ## } else{
      ## 	for (iRegime in (1:n_regimes)) {
      ## 		iRange <- (regimes[iRegime] + 1) : (regimes[iRegime + 1])
      ## 		##breaks_pos indicates new regime starts at next time
      ## 		iU <- u[iRange,]
      ## 		iOmega <- crossprod(iU) / dim(iU)[1] ##reduced form variances
      ## 		iLambda <- seed_Ainv %*% iOmega %*% t(seed_Ainv) ##structural variances, in a sense
      ## 		##Just taking diagonal elements
      ## 		seed_lmd[,iRegime] <- -log(diag(iLambda))
      ## 	}
      ## }

      ## fullLmd <- kronecker(matrix(fullLmd,nrow = length(fullLmd)),t(rep(1,n_regimes)))
    }

    ##like lc_A0, tells program what values to fill from x
    n_regimes <- length(period_breaks) + 1
    lc_lmd <- matrix(TRUE, nrow = n_vars, ncol = n_regimes)
    lc_lmd[,n_regimes] <- FALSE

    if (!is.null(lmdblock)){
      ##generate matrix that tells optimizer what to fill in lmd
      n_blocks <- dim(lmdblock)[3]
      for (iBlock in 1:n_blocks){
        ##could not get logical indexing to work here at all, so I'm doing an awkward loop
        reps <- which(lmdblock[,,iBlock])
        nReps <- length(reps)

        for (iRep in 2:nReps) {
          lc_lmd[reps[iRep]] <- FALSE
        }
        ##(lc_lmd[which(lmdblock[,,iBlock])[-1]]) =
        ##rep(FALSE, length(which(lmdblock[,,iBlock])) - 1) ##fill in first of the block

        if (is.null(seed_model)){
          ## seed_lmd[lmdblock[,,iBlock]] <- fullLmd[lmdblock[,,iBlock]] ##seed with full sample lmd estimates
          seed_lmd[lmdblock[,,iBlock]] <- 1 ##seed with full sample lmd estimates
        }
      }
    }

    ## If number of lambda parameters = number of variables, don't include any lambda parameters
    ## ie variances are constant
    if (sum(lc_lmd) <= n_vars){
      lc_lmd[,] <- FALSE
    }

    ##construct an x vector, if one not already provided
    if (is.null(seed_x)){
      if (tvA == FALSE){              #A0 is constant
        seed_x <- c(c(seed_A[lc_A0]), c(seed_lmd[lc_lmd]))
      } else {
        if (noLmd == TRUE){
          ## Estimate without a prior on lmd --- set constant to 1
          seed_x <- rep(c(seed_A[lc_A0]),times=n_regimes)
        } else {
          ## Estimate with the prior on lmd
          seed_x <- c(rep(seed_A[lc_A0],times=n_regimes), seed_lmd[lc_lmd])
        }

      }
    }


    if (any(hparam_nl > 0)){
      ## add params for nl transformation if desired
      if (length(seed_x) == (sum(lc_A0) + sum(lc_lmd))){
        ## add new params
        sval <- c(.3,.4,1,3,2)
        ## some guesses for a,c,beta,iv,k
        seed_x <- c(seed_x, sval[hparam_nl > 0])

        ## n_aparams <- sum(lc_A0)
        ## n_lmd_params <- sum(lc_lmd)
        ## seed_H <- diag(c(rep(1, n_aparams),
        ##                .1 * c(seed_lmd[lc_lmd]),
        ##                rep(.001,sum(hparam_nl > 0)),
        ##                rep(1e-2, sum(hparam > 0))))
      }
    }

    if (any(hparam > 0)){
      ## add regular hyperparameters if theyre not there
      if (length(seed_x) == (sum(lc_A0) + sum(lc_lmd) + sum(hparam_nl > 0))){


        ## JMO - stop assuming, fix this with warning:
        ## assume we have one of the mn or cosine priors
        ## this is a quick fix
        if (!is.null(cos_prior)){
          pv <-
            c(0,
              0,
              prior_params$ur_prior$lambda,
              prior_params$ur_prior$mu,
              prior_params$cos_prior$tight,
              prior_params$cos_prior$smooth,
              prior_params$cos_prior$damp)
        } else { ## asume mn prior

          pv <-
            c(prior_params$mn_prior$tight,
              prior_params$mn_prior$decay,
              prior_params$ur_prior$lambda,
              prior_params$ur_prior$mu,
              0,
              0,
              0)
        }

        seed_x <- c(seed_x, pv[hparam > 0])
      }
    }


    ## } else {
    ##     ## no nl, or hypers

    ##     ##initial inverse hessian; want log stuff to move less in dx = -H0 * g
    ##     seed_H <- diag(c(rep(1, n_aparams),
    ##                    .1 * c(seed_lmd[lc_lmd])))
    ##     ##seed_H <- diag(length(seed_x))
    ## }

    n_aparams <- sum(lc_A0)

    if (tvA == TRUE){
      ## More A0 parameters
      n_aparams <- n_aparams * n_regimes
    }

    n_lmd_params <- sum(lc_lmd)

    ## seed_H <- diag(c(rep(1, n_aparams),
    ##                .1 * c(seed_lmd[lc_lmd]),
    ##                rep(.001,sum(hparam_nl > 0)),
    ##                rep(1e-2, sum(hparam > 0))))


    if (is.null(seed_H)){
      seed_H <- diag(c(rep(1, n_aparams),
                       1e-5 * rep(1, n_lmd_params),
                       rep(.001,sum(hparam_nl > 0)),
                       rep(1e-2, sum(hparam > 0))))
    }


    # Running optimization ----

    ##Adding zero at beginning of breaks_pos
    breaks_pos <- c(0, breaks_pos)


    ##optoutput <- csminwelNew(bvarwrap5, seed_x, seed_H, nit = 200,
    ##listData <- listData, lc_A0 = lc_A0, breaks_pos = breaks_pos,
    ##prior_params <- prior_params, n_lags = n_lags, grad = seedG)

    # Select function to use:
    if (tvA == FALSE){
      lhfcn <- bvarwrap5
    } else {
      lhfcn <- bvarwrap_tvA
    }


    optoutput <-
      csminwelNew(
        fcn = lhfcn,
        x0 = seed_x,
        H0 = seed_H,
        nit = maxit,
        dat = dat_ts,
        lc_A0 = lc_A0, lc_lmd = lc_lmd,
        lmdblock = lmdblock,
        breaks_pos = breaks_pos,
        prior_params = prior_params,
        nLags = n_lags,
        crit = critval,
        Verbose = FALSE,
        n_cores = n_cores,
        nlt = nlt,
        hparam_nl = hparam_nl,
        hparam = hparam
      )

    ##cleaner pointer names
    lh <- optoutput$fh
    x <- optoutput$xh

    ## checking reason for termination above

    term1 <- optoutput$retcodeh

    # Checking optimization, if appropriate ----

    ## currently, this is triggered for ANYTHING except critical value
    ## termination, or if specified in original TvvVar() call. Does not
    ## have different strategies for different errors (e.g., change gradient
    ## if problem was bad gradient, etc.)

    ## alternatively could make this recursive (keep trying until we reach max iterations
    ## or get a benign termination)

    if (term1 != 0) extracheck <- TRUE

    if (extracheck){
      ## checks if we had problems related to bad Hessian by trying with
      ## 'agnostic' Hessian we started with
      optoutput2 <-
        csminwelNew(fcn = lhfcn,
                    x0 = x,
                    H0 = seed_H,
                    nit = maxit,
                    n_cores = n_cores,
                    dat = dat_ts,
                    lc_A0 = lc_A0,
                    lc_lmd = lc_lmd,
                    lmdblock = lmdblock,
                    breaks_pos = breaks_pos,
                    prior_params = prior_params,
                    nlt = nlt,
                    hparam_nl = hparam_nl,
                    n_lags = n_lags,
                    hparam = hparam,
                    crit = critval,
                    Verbose = FALSE,
                    n_varscores = n_varscores)

      term2 <- optoutput2$retcodeh
      different <- any((x - optoutput2$xh) < (1e4 * .Machine$double.eps))

      # If different, replace lh & x with second run's output (JMO - review this strategy):
      if (different) {
        lh <- optoutput$fh
        x <- optoutput2$xh
      }

    } else {
      # If no extracheck, just do NULL for these as the output:
      different <- NULL
      optoutput2 <- NULL
      term2 <- NULL
    }


    ## Retrieving model and parameters ----

    ##sigpar <- list(A = A, lmd = lmd, breaks_pos = strTsigbrk)
    mlmodel <-
      lhfcn(
        x, verbose = TRUE,
        dat = dat_ts,
        n_lags = n_lags,
        lc_A0 = lc_A0, lmdblock = lmdblock, lc_lmd = lc_lmd, breaks_pos = breaks_pos,
        prior_params = prior_params, n_varscores = n_varscores, nlt = nlt,
        hparam_nl = hparam_nl, hparam = hparam
      )

    ##cleaning up pointer names, etc.
    A <- mlmodel$A
    lmd <- mlmodel$llmd ##ACTUAL lmd
    lambda <- mlmodel$lambda ##exponentiated
    relLambda <- lambda / lambda[,1] ##variances relative to the first period
    vout <- mlmodel$vout
    lh <- mlmodel$lh

    ## Calculating an impulse response ----

    ## This is a somewhat clumsy way of getting the reduced form parameters,
    ## but this step is very quick anyways

    nLambdas <- dim(lambda)[2]
    ir <- array(0, dim = c(n_vars, n_vars, n_step, nLambdas))

    for (iLambda in (1:nLambdas)){

      if (tvA == TRUE){
        iA <- A[,,iLambda]
      } else{
        iA <- A
      }


      smat <- solve(iA) %*% diag(sqrt(lambda[,iLambda]))
      By <- vout$var$By
      n_lags <- dim(By)[3]
      for (iLag in (1:n_lags)) {
        By[,,iLag] <- solve(iA) %*% By[,,iLag]
      }

      tempvar <- vout$var
      tempvar$By <- By
      # JMO - is this nstep different from the nStep (now n_step) used elsewhere?
      ir[,,,iLambda] <- impulsdtrf(tempvar, smat = smat, nstep = 60)
      dimnames(ir[,,,iLambda])[[1]] <- vars
    }

    output <-
      list(A = A, lmd = lmd, lambda = lambda, relLambda = relLambda,
           vout = vout, ir = ir, x = x, breaks_pos = breaks_pos,
           startdate = startdate, enddate = enddate, lc_lmd = lc_lmd,
           lc_A0 = lc_A0, log_vars = log_vars, prior_params = prior_params,
           dat_ts = dat_ts, lmdblock = lmdblock,
           different = different,
           extracheck = extracheck, term1 = term1,
           term2 = term2, n_lags = n_lags, lh = lh)

    if (!verbose){
      return(output)
    } else {
      ##return full optimization output
      ##return(c(output, list(opt = optoutput, opt2 = optoutput2)))
      output$opt <- optoutput
      output$opt2 <- optoutput2
      return(output)
    }
  }


