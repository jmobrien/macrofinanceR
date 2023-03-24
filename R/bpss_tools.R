transform_nl <-
  function(dat_ts, nl){

    ## applies a simple nonlinear transformation (see paper for details), and calculates
    ## appropriate jacobian term to put in likelihood (for correct model comparison, or
    ## inference over the nl parameters).

    att <- nl[1] ## lower treshold
    ct <- nl[2] ## upper
    beta <- nl[3] ## slope increase
    iv <- nl[4] ## which variable to worry about
    ## alpha <- nl[5] ## slope of the original log transformation (x-p) - ky

    iseries <- dat_ts[,iv]

    jterm <- rep(0,length(iseries)) ## jacobian term

    ## quadratic range
    a1 <- beta / (2 * (ct - att))
    a2 <- 1 - 2 * att * a1
    a3 <- att - a2 *att - a1 * att^2

    iseries[(iseries > att) & (iseries < ct)] <-
      a3 + a2 * (iseries[(iseries > att) & (iseries < ct)]) +
      a1 * (iseries[(iseries > att) & (iseries < ct)])^2

    jterm <-
      sum(log(a2 + 2 * a1 * (dat_ts[(dat_ts[,iv] > att) & (dat_ts[,iv] < ct),iv])))


    ## upper linear range
    ct2 <- a3 + a2*ct + a1 * ct^2 ## f(ct)

    iseries[(dat_ts[,iv] >= ct)] <-
      (1+beta) * (iseries[(dat_ts[,iv] >= ct)] - ct) + ct2
    jterm <- jterm + sum(dat_ts[,iv] > ct) * log(1+beta)

    ## grDevices::pdf('~/check.pdf')
    ## ## irange <- (mydata[,iv] > att) & (mydata[,iv] < ct)
    ## ## plot(mydata[irange,iv],iseries[irange])
    ## plot(mydata[,iv],iseries,type <- 'l')
    ## dev.off()


    dat_ts[,iv] <- iseries

    return(list(dat_ts = dat_ts, jterm = jterm))
  }




#' gets a forecast nperiods into the future
#'
#' @param tout full output of optimization
#' @param nperiods number of periods into future to forecast
#' @param Tobs Total number of observations
#'
#' @return
#' @export
#'
get_forecast <-
  function(tout,
           nperiods = 8,
           Tobs = dim(ydata)[1]){

    ## gets a forecast nperiods into the future
    ## starts at cap Tobs

    ## input to this function is full output of optimization


    ydata <- as.matrix(tout$dat_ts)

    ## get the reduced form coefficients
    A <- tout$A


    Ai <- solve(A)

    vout <- tout$vout

    By <- vout$var$By
    nv <- dim(By)[1]
    nlags <- dim(By)[3]

    Bx <- Ai %*% matrix(vout$var$Bx,nv,1)

    for (ilag in 1:nlags){
      By[,,ilag] <- Ai %*% By[,,ilag]
    }

    ## get the system matrix
    sys <- sysmat(By = By, Bx = Bx)
    ## get the data vector

    ## Tobs <- dim(ydata)[1]
    y <- matrix(0,nv,nlags)

    for (ilag in 1:nlags){
      y[,ilag] <- ydata[Tobs - ilag + 1,]
    }

    y <- c(y,1) # constant at the end

    ## calculate the predictions
    pmat <- matrix(0,nperiods,nv)

    for (iperiod in 1:nperiods){
      y <- sys %*% y
      pmat[iperiod,] <- y[1:nv]
    }

    return(pmat)
  }


#' Title Small wrapper function for calculated mdd with MHM methods
#'
#' @param x draws in a matrix
#' @param lh vector of -log posterior densities (conditional on certain params)
#' @param trunc extra piece numerator density for MHM
#' @param delta delta_it in an array, if appropriate
#' @param alpha,beta parameters for the distribution of the delta (maybe you should change this)
#' @param efac extra piece numerator density for MHM
#' @param covsub whether to calc covariance matrix with all xout, or 5000 evenly spaced draws
#'
#' @return
#' @export
#'
get_mdd <-
  function(x, lh, trunc = .95, delta = NULL,
           alpha =3, beta = 3, efac = NULL,
           covsub = FALSE){


    ## Small wrapper function for calculated mdd with MHM methods
    ## automatically puts gaussian distribution on x truncated for trunc density
    ## can put extra density in efac (e.g., on the delta_it)

    ## x is the draws in a matrix
    ## lh is a vector of -log posterior densities (conditional on certain params)
    ## efac is extra piece numerator density for MHM
    ## covsub determines whether to calc covariance matrix with all xout, or 5000 evenly spaced draws
    ## of it (mainly to save time)

    ## trunc is the truncation of the Gaussian (what percentage is kept) for A0 and Lambda
    ## delta is the delta_it in an array, if appropriate
    ## alpha and beta are for the distribution of the delta (maybe you should change this)

    ## END PREAMBLE




    ## gets estimates of the MDD for the model

    ## first, get a Gaussian approx
    xmean <- apply(x,2,mean)
    Tobs <- dim(x)[1]
    nv <- dim(x)[2]
    ## covmat <- crossprod(x) / Tobs

    if (covsub){
      ## use subsample for the covariance matrix
      xsub <- seq(from=1,to=Tobs,length.out = 5000)
      covmat <- cov(x[xsub,])
    } else {
      covmat <- cov(x)
    }

    ## If covariance matrix is near singular, add noise
    if (min(eigen(covmat)$values) < 1e-8){
      covmat <- covmat + 1e-3 * diag(dim(covmat)[1])
    }


    ## get det
    cdet <- sum(log(eigen(covmat)$values))

    ## ## evaluate gaussian approx
    ci <- (chol(solve(covmat)))
    xdm <- t(x) - xmean
    cixdm <- ci %*% xdm
    qscore <- apply(cixdm^2,2,sum) ## quadratic part

    ## qscore <- diag(t(xdm) %*% solve(covmat) %*% (xdm))

    ## qscore <- rep(0,dim(xdm)[1])
    ## for (i in 1:dim(xdm)[1]){
    ##     qscore[i] <- (xdm[i,]) %*% covmat %*% t(xdm[1,])
    ## }


    dout <- -.5 * qscore - .5 * (nv * log(2 *pi) + cdet) ## density

    kval <- qscore < qchisq(trunc,nv) ## truncate the tails of the gaussian

    if (!is.null(delta)){ ## put some distributio on the delta

      ## The conditional posterior is inverse gamma for the t case,
      ## and a known multinomial for the normal mix case
      ## so we should probably just use those distributions independent for each
      ## not what's implemented here


      dmat <- matrix(delta, dim(delta)[1] * dim(delta)[2], dim(delta)[3])

      ## ## trial 1: ivnerse gamma
      ## alpha <- 3
      ## beta <- 3 ## inverse gamma parameters
      ## dweight <- -(alpha-1) * (2 * dmat) - beta / exp(2 * dmat)
      ## dweight <- dweight + alpha * log(beta) - lgamma(alpha)

      ## ## trial 2: regular gamma
      dweight <- log(dgamma(exp(2 * dmat), shape = alpha, rate = 1/alpha))
      dweight <- apply(dweight,2,sum)
      dout <- dout + dweight
    }

    gdout <- dout + lh - log(trunc) ## before truncation

    gdmean <- mean(gdout)

    gdfix <- gdout
    gdfix[!kval] <- -1e10 ## zero values

    ##return(-(gdmean + log(mean(exp(gdfix-gdmean)))))
    if (is.null(efac)){
      return(-(gdmean + log(mean(exp(gdfix-gdmean)))))
    } else {
      m1 <- -(gdmean + log(mean(exp(gdfix-gdmean))))

      ## with the additonal correction
      gdout <- dout + efac + lh - log(trunc)


      gdmean <- mean(gdout)
      gdfix <- gdout
      gdfix[!kval] <- -1e10 ## zero values
      m2 <- -(gdmean + log(mean(exp(gdfix-gdmean))))

      outlist <- list(mout = c(m1,m2), cv = covmat)
      return(outlist)
    }
  }



#' Title get root mean squared error
#'
#' @param fcast fcast array
#' @param ydata model data for Y
#' @param h
#'
#' @return
#' @export
#'
getrmse <-
  function(fcast,ydata,h = c(1,6,12,24,48)){

    ## get root mean squared error
    ## input is "fcast" array

    nv <- dim(fcast)[2] ## number of variables
    nfc <- dim(fcast)[3] ## number of forecasts
    nh <- length(h) ## number of horizons to consider
    Tobs <- dim(ydata)[1]

    rmse <- array(0, c(nh,nv,nfc))

    for (it in 1:nfc){
      if (sum(abs(fcast[,,it])) > 1e-10){ ## make sure fc is not blank
        for (ih in 1:nh){
          hzn <- h[ih] ## what horizon
          if ((hzn + it) <= Tobs){ ## make sure there is space
            df <- fcast[1:hzn,,it] - ydata[it + 1:hzn,]
            if (hzn > 1){
              rmse[ih,,it] <- sqrt(apply(df^2,2,sum) / hzn)
            } else {
              rmse[ih,,it] <- abs(df)
            }

          }
        }
      }
    }
    return(rmse)
  }















