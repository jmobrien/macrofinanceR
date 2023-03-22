
McmcIr <-
  function(xout, model, rootcheck = FALSE, lrange = 1:6, cores = 8,
           oweights = NULL,
           aplus = NULL, ## don't need both oweights and aplus
           nStep = 60, hparam = rep(0,7))
  {

    ## Main wrapper function for drawing impulse responses

    ## --------------------INPUTS, for all uses--------------------
    ## xout : A0 and Lambda draws in matrix form
    ## rootcheck : if TRUE, will check if roots are stable
    ## lrange : which variance regimes you want IR for. The IR are all the same up to scale
    ## cores: number of cores to use for mclapply. This is easily parallelized, so use many!
    ## oweights: draws of delta_it
    ## aplus : draws of A+. Code will run a lot faster if you specify this
    ## nStep : number of steps of impulse responses



    ## --------------------OUTPUTS--------------------
    ## list with element ir, the impulse responses [nvar x nshock x nperiod x nregime x ndraws]
    ## other stuff which is less useful (mainly related to orderings/permutations)

    ## Karthik Sastry
    ## R 3.1.2, 64 Bit
    ## January 2017

    ## END PREAMBLE

    vars <- colnames(model$A)

    nTrials <- dim(xout)[1]

    nX <- dim(xout)[1]
    if (is.null(nX)) {
      nX <- 1
      dim(xout) <- c(1,length(xout))
    }

    if (!is.null(oweights)){ ## add these to x
      ## xout <- cbind(xout, matrix(0,dim(xout)[1],
      ##                           dim(oweights)[1] * dim(oweights)[2]))
      ## xout[,-c(1:dim(xout)[2])] <- matrix(oweights,dim(xout)[1],
      ##                                    dim(oweights)[1] * dim(oweights)[2],
      ##                                    byrow = TRUE)

      owmat <- matrix(oweights,dim(xout)[1],
                      dim(oweights)[1] * dim(oweights)[2],
                      byrow = TRUE)
      xout <- cbind(xout, owmat)

      owflag <- TRUE
    } else {
      owflag <- FALSE
    }

    if (!is.null(aplus)){
      na <- dim(aplus)[1] ## need to remove last row
      amat <- matrix(aplus[-na,,],
                     dim(aplus)[3],
                     (na-1) * dim(aplus)[2],
                     byrow = TRUE)
      xout <- cbind(xout, amat)
      aflag <- TRUE
    } else{
      aflag <- FALSE
    }


    ## if (!is.null(aplus)){ ## add these to x
    ##     xout <- cbind(xout, matrix(0,dim(xout)[1],
    ##                               dim(oweights)[1] * dim(oweights)[2]))
    ##     xout[,-c(1:dim(xout)[2])] <- matrix(oweights,dim(xout)[1],
    ##                                        dim(oweights)[1] * dim(oweights)[2],
    ##                                        byrow = TRUE)
    ## }


    listXout <- lapply(1:nX, function(iRow){xout[iRow,]})

    listIrtrials <- parallel::mclapply(listXout, IrRun2, model$A,
                                       model$lmd, model$dat_ts, model$lcA,
                                       model$lc_lmd, nLags = model$nLags,
                                       owflag = owflag, aflag = aflag,
                                       lmdblock = model$lmdblock,
                                       hparam = hparam,
                                       breaks_pos = model$breaks_pos, prior_params = model$prior_params,
                                       lmdPrior = model$lmdPrior, rootcheck = rootcheck,
                                       lrange = lrange, nStep = nStep, mc.cores = cores)

    ##save(listIrtrials, file = 'rawIR.Rdata')

    ##irdim <- dim(listIrtrials[[1]]$ir)

    ##arrayIrtrials = array(sapply(listIrtrials, '[[', 1), dim = c(irdim, nTrials))

    ##percentiles <- apply(arrayIrtrials, 1:4, FUN = quantile, probs = c(.05,.5,.95))

    elements <- c('ir','ch','orderings','badshift','trimprove','x','lh')
    if (rootcheck) elements <- c(elements, 'isgood','maxmod', 'maxRev', 'iscomplex','modmod')
    output <- list()
    for (iElement in elements){
      iObject <- sapply(listIrtrials,'[[',iElement)
      if (is.null(dim(iObject))) iObject <- unlist(iObject)
      output[[iElement]] <- iObject
    }

    output$ir <- array(output$ir, dim = c(length(vars), length(vars), 60, length(lrange), nTrials))

    if (inherits(xout, 'mcmc')){
      output$x <- coda::mcmc(t(output$x), thin = thin(xout))
    }
    output$ess <- coda::effectiveSize(output$x)
    if (rootcheck) output$iscomplex <- t(output$iscomplex)

    return(output)
  }





IrRun2 <-
  function(x, modeA, modeLmd, dat_ts, lc_A0, lc_lmd, nLags = nLags,
           owflag = FALSE, aflag = FALSE, ..., nStep = 60, rootcheck = FALSE, lrange = 1:6){
    ## Main function for drawing impulse responses
    ## It's written in a funny way, with a single draw input x, to facilitate use with
    ## lapply and its multicore version, mclapply

    ## --------------------INPUTS, for all uses--------------------
    ## x : in the base case, this includes only A0 and Lambda
    ## modeA,modeLmd : for unused parts of code that checked ordering. Set to NULL
    ## lc_A0 and lc_lmd : logical arrays showing restrictions on A0 and Lmd. Easiest to let these
    ## be populated automatically by the parent function McmcIr
    ## nStep : number of steps of impulse responses
    ## lrange : which variance regimes you want IR for. The IR are all the same up to scale

    ## --------------------INPUTS, if you've drawn A+--------------------
    ## x : vector of (A0,Lambda,A+)
    ## aflag : set to TRUE
    ## rest of stuff: doesn't really matter

    ## --------------------INPUTS, if you haven't drawn A+--------------------
    ## dat_ts : data in correct format
    ## owflag : if TRUE, delta_it are at the end of x
    ## rootcheck : if TRUE, will check if roots are stable
    ## ... : extra args passed to LH evalulation function

    ## --------------------OUTPUTS--------------------
    ## list with element ir, the impulse responses [nvar x nshock x nperiod x nregime]
    ## other stuff which is not so useful


    ## Karthik Sastry
    ## R 3.1.2, 64 Bit
    ## January 2017

    ## END PREAMBLE


    ## if (length(x) > (sum(lc_A0) + sum(lc_lmd))){ ## oweights in the mix
    if (owflag){
      oweights <- matrix(x[-(1:(sum(lc_A0) + sum(lc_lmd)))],
                         dim(dat_ts)[1] - nLags,dim(lc_A0)[1])
    } else {
      oweights <- NULL
    }

    if (aflag){ ## aplus is already included
      dimAplus <- c(nrow(lc_A0),ncol(lc_A0),nLags)
      Aplus <- matrix(x[-(1:(sum(lc_A0) + sum(lc_lmd)))],(nrow(lc_A0))^2,nrow(lc_A0))

      A <- matrix(0,nrow(lc_A0),ncol(lc_A0))
      A[lc_A0] <- x[1:(sum(lc_A0))]


      lambda <- matrix(0, dim(lc_lmd)[1],dim(lc_lmd)[2])
      lambda[lc_lmd] <- x[sum(lc_A0) + (1 : sum(lc_lmd))]
      lambda[,dim(lambda)[2]] <- dim(lambda)[2] -
        apply(lambda[,1:(dim(lambda)[2]-1),drop=FALSE],1,sum)

      nVar <- dim(A)[1]
      vars <- colnames(lc_A0)

      lh <- lh


    } else { ## needs to be calculated

      mlmodel <- bvarwrap5(x, dat = dat_ts, lc_A0 = lc_A0, lc_lmd = lc_lmd, oweights = oweights, verbose = TRUE, nLags = nLags, ...)

      lh <- mlmodel$lh

      ##add gaussian noise to A+
      ##covariance matrix
      Aplus <- mlmodel$vout$var$By
      u <- mlmodel$vout$var$u
      xx <- mlmodel$vout$var$xx
      nEq <- dim(xx)[3]
      nLags <- dim(Aplus)[3]
      xxi <- list(nEq)
      nX <- nEq * nLags * nEq + nEq
      ##arraySigma <- array(0, dim = c(nX, nX, nEq))
      ##xxi = solve(xx[,,1])
      for (iEq in (1:dim(xx)[3])){
        xxi[[iEq]] <- solve(xx[,,iEq])
      }
      ##getting correct Sigma matrix
      Sigma <- bdiag(xxi)
      ##coercing to matrix class....
      Sigma <- matrix(Sigma, nrow = nX)

      dimAplus <- dim(Aplus)
      ##for (iEq in 1:nEq) {
      ##	xxi[,,iEq] = solve(xx[,,iEq])
      ##}
      covu <- cov(u)
      ##Sigma <- kronecker(covu, xxi)
      ##Sigma <- kronecker(diag(dim(covu)[1]), xxi)

      ##calculating the noise
      nDraws <- dim(Sigma)[1]
      ##sndraws <- rnorm(nDraws) ##standard normal draws
      ##draws <- t(chol(Sigma)) %*% sndraws
      nAplus <- length(Aplus)
      ##mlmodel$vout$var$By <- array(c(mlmodel$vout$var$By) +  draws[1:nBy],
      ##  dim <- dim(mlmodel$vout$var$By))
      ##should not need to permute dimensions, change to matrix, etc.

      ##just to make sure, this is the postdraw() style implementation
      sndraws <- matrix(rnorm(nDraws), c(nDraws, 1))
      matAplus <- t(matrix(Aplus, nEq, nEq * nLags))
      #### Aplusdraws <- t(chol(Sigma)) %*% sndraws
      Aplusdraws <- MASS::mvrnorm(1,rep(0,dim(Sigma)[1]),Sigma) #### new method?
      Aplusnoise <- matrix(Aplusdraws, c(nEq * nLags + 1, nEq)) ##i.e., the matrix of noise that I am adding to the posterior mean/mode
      Aplus  <- matAplus + Aplusnoise[1:(nEq*nLags),] ##ignoring the constant, which is last row

      ##cleaning pointers, etc.
      A <- mlmodel$A
      vout <- mlmodel$vout
      lambda <- mlmodel$lambda
      nVar <- dim(A)[1]
      vars <- colnames(A)

    }


    Aplus  <- array(t(Aplus), dim = dimAplus)

    ##Aplus <- mlmodel$vout$var$By ##cancelling variation in Aplus
    ##above I only unocmmented to see if X'X singularites were creating problems


    ##nLambdas <- dim(lambda)[2]
    ##Now a lambda range is set in the function call

    ir <- array(0, dim = c(nVar, nVar, nStep, length(lrange)))

    ##        isgood <- rep(0, nLambdas)
    ##        maxmod <- rep(0, nLambdas) ##potentially are same across variance regimes, but I want to check
    ##        maxRev <- rep(0, nLambdas)
    ##        iscomplex <- rep(0, nLambdas)

    for (iLambda in lrange){
      By <- array(0, dim = dimAplus)
      Ainv <- solve(A)
      smat <- Ainv %*% diag(sqrt(lambda[,iLambda]))
      for (iLag in (1:nLags)) {
        By[,,iLag] <- Ainv %*% Aplus[,,iLag]
      }

      tempvar <- list()
      tempvar$By <- By
      ir[,,,which(lrange == iLambda)] <- impulsdtrf(tempvar, smat = smat, nstep = nStep)
      dimnames(ir[,,,which(lrange == iLambda)])[[1]] <- vars
    }

    newordering <- NULL ## didn't check
    badshift <- NULL
    if (rootcheck){
      rc <- CheckRoots(By)
      isgood <- rc$isgood
      maxmod <- rc$maxmod
      iscomplex <- rc$iscomplex
      maxRev <- rc$maxRev
      modmod <- rc$modmod
      return(list(ir=ir, noloop = newordering$noloop, badshift = badshift,
                  orderings = newordering$ordrng, ch = newordering$changes,
                  trimprove = newordering$trimprove, x = x, isgood = isgood,
                  maxmod = maxmod, iscomplex = iscomplex,
                  maxRev = maxRev, modmod = modmod, lh = lh))
    } else{
      return(list(ir = ir,noloop = newordering$noloop, badshift = badshift,
                  orderings = newordering$ordrng, ch = newordering$changes,
                  trimprove = newordering$trimprove, x = x, lh = lh))
    }
  }


