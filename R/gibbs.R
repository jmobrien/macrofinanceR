#' gibbs draws
#'
#' @param tvv_mod
#' @param ndraw
#' @param lhfcn
#' @param nburn
#' @param nsep
#' @param filename
#' @param hsnfactor
#' @param alpha
#' @param k
#' @param tparam
#' @param tscale
#' @param savespots
#' @param dscore
#' @param drawdout
#' @param drawe
#' @param drawa
#' @param hparam_nl
#' @param nlt
#'
#' @return
#' @export
#'
#' @examples
gdraw <-
  function(tvv_mod, ndraw,
           lhfcn = bvarwrap5,
           nburn = 1e3,
           nsep = 1e2,
           filename = 'gibbsout.Rdata',
           hsnfactor = 0.2,
           alpha = .003, k = 4,
           tparam = NULL,
           tscale = 1,
           savespots = NULL,
           dscore = FALSE,
           drawdout = FALSE, ## if TRUE, keep the delta_it (the xi_it)
           drawe = FALSE, ## if TRUE, keep the e_it
           drawa = TRUE,  ## if TRUE, keep the A+ draws
           hparam_nl = rep(0,5),        # indicates if any hyperparams for nonlinear transformation are used
           nlt = rep(0,5)               # constant parameters for nlt
  ){

    ## Main wrapper function for Gibbs Sampling from posterior of normal mixture models
    ## Note that this can also be used to sample from the Normal errors model
    ## (for instance, by setting alpha = 1 and k = 1)

    ## The real sampling happens in gstep.R

    ## --------------------INPUTS--------------------
    ## tvv_mod : Output of optimization for normal errors model
    ##            <<IMPORTANT>> that this structure includes the data and some other model parameters!
    ##            So if you want to customize everything, please see TvvDir.R
    ## ndraw : Number of draws that will be //recorded//
    ## nburn : Number of initial draws discarded as burn-in. As discussed in the Appendix, etc.,
    ##         it is probably smarter to save everything and monitor convergence "manually"
    ## nsep : The separation between draws. The code is written to make nsep draws and then save
    ##        one draw, so the actual "separation" of draws is nsep-1
    ## filename : remember the extension
    ## hsnfactor : by default, the code looks at tvv_mod$opt$H for a covariance matrix (inv. Hessian)
    ##             and scales it by hsnfactor (in //variance// not standard deviation units).
    ## alpha, k : specify these for the n-normal mixture model. alpha is a vector of probabilities.
    ##            k is a vector of standard deviations for each normal
    ## tparam, tscale : specify these for t errors model. tparam is df/2 (which is the first param
    ##                  for the inverse gamma distribution). tscale is the scale of the t
    ##                  which is not identified by the model and easiest to leave as unit
    ## savespots : a list of draws (between 1 asnd ndraw) at which to save. Code will generate some
    ##             automatically if left null. Note the code create a connection for output that
    ##             can be flushed at savespots to reduce memory burden
    ## dscore : save the inverse gamma density evaluated at the variance shocks, useful for MDD calculations
    ## drawdout : if TRUE, keep (save) the delta_it (called xi_it in the paper).
    ## drawe : if TRUE, keep the structural residuals
    ## drawa : if TRUE, keep the A+ (necessary for impulse responses, etc!)

    ## --------------------OUTPUT--------------------
    ## a list with elements:
    ## xout : the A0, Lambda draws
    ## lhout : the negative log posterior from metropolis step
    ## tout : 0 if there was no metropolis move, 1 if metropolis move
    ## dout : array with the delta_it (xi_it)
    ## dsout : matrix of log ig density evaluated at dout
    ## aout : A+ draws
    ## eout : structural residual draws



    ## Karthik Sastry
    ## R 3.1.2, 64 Bit
    ## August 2016

    ## END PREAMBLE


    hsn <- tvv_mod$opt$H
    x0 <- tvv_mod$x

    lc_A0 <- tvv_mod$lc_A0
    lc_lmd <- tvv_mod$lc_lmd
    ## ndraw <- nsep* ndraw + nburn
    nvar <- dim(lc_A0)[1]

    breaks_pos <- tvv_mod$breaks_pos
    prior_params <- tvv_mod$prior_params
    dat_ts <- tvv_mod$dat_ts
    lmdblock <- tvv_mod$lmdblock
    n_lags <- tvv_mod$n_lags


    Sigma <- hsnfactor * hsn

    ## decide where to save
    if (is.null(savespots)){
      nsaves <- 5
      inc <- ndraw %/% nsaves
      if (inc == 0){
        savespots <- ndraw ## no need to save until end
      } else {
        savespots <- c(seq(1, ndraw, inc), ndraw)
      }
    }


    ## allocate space for output
    xout <- matrix(0,length(x0),ndraw)
    lhout <- rep(0,ndraw)
    tout <- rep(0,ndraw)

    ## always save delta, a, e in savespots
    if (drawdout){ ## delta_it
      dout <- array(0,
                    c(dim(dat_ts)[1] - n_lags,
                      dim(dat_ts)[2],
                      ndraw))
    }  else {
      dout <- array(0,
                    c(dim(dat_ts)[1] - n_lags,
                      dim(dat_ts)[2],
                      length(savespots)))
    }

    if (drawa){ ## A+
      aout <- array(0,
                    c(n_lags * nvar + 1,
                      nvar,
                      ndraw))
    }  else {
      aout <- array(0,
                    c(n_lags * nvar + 1,
                      nvar,
                      ndraw))
    }

    if (drawe){ ## epsilon_it
      eout <- array(0,
                    c(dim(dat_ts)[1] - n_lags,
                      dim(dat_ts)[2],
                      ndraw))
    }  else {
      eout <- array(0,
                    c(dim(dat_ts)[1] - n_lags,
                      dim(dat_ts)[2],
                      length(savespots)))
    }

    isave <- 1 ## first save

    if (dscore){
      dsout <- rep(0,ndraw)
    } else {
      dsout <- NULL
    }

    output <- NULL

    if (!is.null(tvv_mod$dout)){
      dout[,,1] <- tvv_mod$dout
    }


    gout <- list(xout = x0, deltaout = dout[,,1])


    ## Run the Gibbs Sampler
    if (nburn > 0){
      for (iburn in 1:nburn){
        gout <- gstep(gout, lhfcn, Sigma, lc_A0, lc_lmd, alpha, k,
                      lmdblock = lmdblock, breaks_pos = breaks_pos,
                      prior_params = prior_params, dat = dat_ts, n_lags = n_lags,
                      tparam = tparam, tscale = tscale, drawbe = TRUE,
                      hparam_nl = hparam_nl, nlt = nlt)
      }
    }
    for (idraw in 1:ndraw){
      for (isep in 1:nsep){
        gout <- gstep(gout, lhfcn, Sigma, lc_A0, lc_lmd, alpha, k,
                      lmdblock = lmdblock, breaks_pos = breaks_pos,
                      prior_params = prior_params, dat = dat_ts, n_lags = n_lags,
                      tparam = tparam, dscore = dscore, tscale = tscale, drawbe =TRUE,
                      hparam_nl = hparam_nl, nlt = nlt)
      }

      ## update
      xout[,idraw] <- gout$xout
      lhout[idraw] <- gout$lhout
      tout[idraw]  <- gout$trans
      if(drawdout) {
        dout[,,idraw] <- gout$deltaout
      }

      if(drawa) {
        aout[,,idraw] <- gout$aout
      }

      if(drawe) {
        eout[,,idraw] <- gout$eout
      }

      if (dscore){
        dsout[idraw] <- gout$dsout
      }
      if(idraw %% 100 ==1) print(paste(filename, "draw number", idraw))
      if (idraw %in% savespots){
        dout[,,isave] <- gout$deltaout
        isave <- isave + 1
        output <- list(xout = xout, lhout = lhout,
                       tout = tout, dout = dout,
                       dsout = dsout, aout = aout, eout = eout)
        outcon <- file(filename, open="wb")
        save(output, file = outcon)
        flush(outcon)
        close(outcon)
      }

    }

    if (is.null(output)){
      output <- list(xout = xout, lhout = lhout,
                     tout = tout, dout = dout,
                     dsout = dsout, aout = aout, eout = eout)
    }
    outcon <- file(filename, open="wb")
    save(output, file = outcon)
    flush(outcon)
    close(outcon)
    return(output)
  }


#' Title
#'
#' @param ydata
#' @param xdata
#' @param olsfun
#' @param ndraw
#' @param x0
#' @param nburn
#' @param nsep
#' @param filename
#' @param sigfactor
#' @param tparam
#' @param tscale
#' @param savespots
#' @param dscore
#' @param drawdout
#'
#' @return
#' @export
#'
#' @examples
gdraw_1v <-
  function(ydata, xdata, olsfun, ndraw,
           x0,
           nburn = 1e3, nsep = 1e2,
           filename = 'gdraw_1v.Rdata',
           sigfactor = 0.2,
           tparam = NULL,
           tscale = 1,
           savespots = NULL,
           dscore = FALSE,
           drawdout = FALSE){

    ## Gibbs sampling for normal mixture single equation models
    ## see gdraw.R for more details

    ## Karthik Sastry
    ## R 3.1.2, 64 Bit
    ## August 2016

    ## END PREAMBLE

    ## decide where to save
    if (is.null(savespots)){
      nsaves <- 5
      inc <- ndraw %/% nsaves
      if (inc == 0){
        savespots <- ndraw ## no need to save until end
      } else {
        savespots <- c(seq(1, ndraw, inc), ndraw)
      }
    }


    ## allocate space for output
    xout <- matrix(0,length(x0),ndraw)
    lhout <- rep(0,ndraw)
    tout <- rep(0,ndraw)

    if (drawdout){
      dout <- matrix(0,length(ydata),ndraw)
    }  else {
      dout <- matrix(0,length(ydata),length(savespots))
    }

    isave <- 1 ## first save

    if (dscore){
      dsout <- rep(0,ndraw)
    } else {
      dsout <- NULL
    }

    output <- NULL

    gout <- list(xout = x0, deltaout = dout[,,1])


    ## Run the Gibbs Sampler
    if (nburn > 0){
      for (iburn in 1:nburn){
        gout <- gstep_1v(gout, olsfun, ydata, xdata,
                         tparam, tscale, dscore)
      }
    }
    for (idraw in 1:ndraw){
      for (isep in 1:nsep){
        gout <- gstep_1v(gout, olsfun, ydata, xdata,
                         tparam, tscale, dscore)
      }

      ## update
      xout[,idraw] <- gout$xout
      lhout[idraw] <- gout$lhout
      tout[idraw]  <- gout$trans
      if(drawdout) {
        dout[,idraw] <- gout$deltaout
      }

      if (dscore){
        dsout[idraw] <- gout$dsout
      }

      if (idraw %in% savespots){
        dout[,,isave] <- gout$deltaout
        isave <- isave + 1
        output <- list(xout = xout, lhout = lhout,
                       tout = tout, dout = dout,
                       dsout = dsout)
        save(output, file = filename)
      }

    }


    if (is.null(output)){
      output <- list(xout = xout, lhout = lhout,
                     tout = tout, dout = dout,
                     dsout = dsout)
    }

    save(output, file = filename)
    return(output)
  }


gstep_1v <-
  function(gout, olsfun, ydata, xdata,
           tparam, tscale, dscore, lmdparam = NULL,Sigma = NULL){

    ## Gibbs sampling step for normal mixture single equation models
    ## see gdraw_1.R for more details

    ## Karthik Sastry
    ## R 3.1.2, 64 Bit
    ## August 2016

    ## END PREAMBLE

    if (is.null(lmdparam)){
      ## OLS step
      olsout <- olsfun(ydata,xdata,gout$deltaout)
    } else{ ## mcmc step
      lh0 <- olsfun(gout$xout,ydata,xdata,gout$deltaout,lmdparam)
      olsout <- GsnMove(olsfun,gout$xout,lh0,Sigma,
                        gout$deltaout,verbose = TRUE)
    }

    uout <- olsout$u * exp(gout$deltaout) ## multiply by last time's deltaout

    if (is.null(tparam)){
      deltaout <- drawdelta(u = uout,alpha = alpha, k = k)
    } else { ## draw inverse gamma
      deltaout <- drawt(u = uout, alpha = tparam, beta = (tscale^2) * tparam)
    }

    if (dscore){ ## reduce delta to a f(delta) for mdd calculation
      dsout <- fd(deltaout, tparam, tscale = tscale)
    } else {
      dsout <- NULL
    }

    return(list(xout = olsout$x, deltaout = deltaout,
                tout = tout, lhout = olsout$lh, dsout = dsout))

  }


gstep <-
  function(gout, lhfcn, Sigma, lc_A0, lc_lmd,
           alpha, k,
           tparam = NULL,
           tscale = 1,
           dscore = FALSE, ## if TRUE, calculate the IG density of delta_it
           ...){

    ## The main iterative step of our Gibbs sampler
    ## Note that most model-related parameters are passed in the ellipses! So see
    ## gdraw.R for more info

    ## --------------------INPUTS--------------------
    ## gout : a list which contains
    ##        xout, draw of A0 and Lambda
    ##        deltaout, draw of delta_it
    ##        trans, 1 iff metropolis step moved
    ##        lhout, negative log posterior
    ##        dsout, ig density evaluated at delta
    ##        eout, structural resids
    ##        aout, A+
    ## tparam, tscale : specify these for t errors model. tparam is df/2 (which is the first param
    ##                  for the inverse gamma distribution). tscale is the scale of the t
    ##                  which is not identified by the model and easiest to leave as unit
    ## dscore : save the inverse gamma density evaluated at the variance shocks, useful for MDD calculations

    ## --------------------OUTPUT--------------------
    ## see gout input


    ## Karthik Sastry
    ## R 3.1.2, 64 Bit
    ## August 2016

    ## END PREAMBLE


    ## markov step
    x <- gout$xout

    model <- lhfcn(x,lc_A0,lc_lmd, oweights = gout$deltaout,...,verbose = TRUE)
    almdout <- GsnMove(lhfcn, x, model$lh, Sigma,
                       modeA = NULL,
                       modeLmd = NULL,
                       lc_A0,
                       lc_lmd,
                       model = model,
                       oweights = gout$deltaout,
                       verbose = TRUE, ...)


    ## drawing the new variance weightso
    if (almdout$lh > 1e4){ ## bad draw
      ## this option should never come into play --- it would only happen if both the initial
      ## AND proposed draws were bad for some reason. It is probably possible if the variance weight
      ## from the inverse gamma is large enough to give the linear regression numerical problems, in
      ## which case you will get a singularity in the LH calculation. So this step imposes some kind of
      ## truncation on the prior for those (not completely clear what).
      return(gout)
    } else {

      ## next line does not account for the fact that we drew A+
      ## uout <- almdout$u[1:dim(gout$deltaout)[1],] * exp(gout$deltaout) ## multiply by last time's deltaout

      aout <- almdout$model$vout$var$Bdraw ## A+ draws
      eout <- almdout$model$vout$var$udraw[1:dim(gout$deltaout)[1],] * exp(gout$deltaout) ## multiply by last time's deltaout to get "real" structural resids

      if (is.null(tparam)){
        deltaout <- drawdelta(u = eout,alpha = alpha, k = k)
      } else { ## draw inverse gamma
        deltaout <- drawt(u = eout, alpha = tparam, beta = (tscale^2) * tparam)
      }

      if (dscore){ ## reduce delta to a f(delta) for mdd calculation
        dsout <- fd(deltaout, tparam)
      } else {
        dsout <- NULL
      }

      return(list(xout = almdout$x, deltaout = deltaout,
                  trans = almdout$trans, lhout = almdout$lh, dsout = dsout,
                  eout = eout, aout = aout))
    }
  }




drawdelta <-
  function(uout, alpha = .01, k = 4){

    ## draws from mixture of normals
    ## alpha is the probability OF an outlier

    ## Karthik Sastry
    ## R 3.1.2, 64 Bit
    ## August 2016

    ## END PREAMBLE



    u <- c(uout) ## treat everything symmetrically

    if (length(k) == 1){
      nprob <- dnorm(u) ## if standard
      oprob <- dnorm(u/k) ## if outlier

      post <- (alpha * oprob) / (alpha * oprob + k * (1-alpha) * nprob) ## posterior proability

      delta <- log(k) * (runif(length(u)) < post) ## take delta with this probability
    } else { ## n generalization

      scalemat <- array(rep(k, times = rep(length(u),length(k))),
                        dim = c(dim(uout),length(k)))
      alphamat <- array(rep(alpha, times = rep(length(u),length(k))),
                        dim = c(dim(uout),length(alpha)))

      ubig <- array(uout, dim = c(dim(uout),length(k)))

      ratio <- alphamat * dnorm(ubig / scalemat) / scalemat

      prob <- ratio / c(apply(ratio,c(1,2),sum)) ## normalized to probability
      cprob <- aperm(apply(prob,c(1,2),cumsum),c(2,3,1))


      ## assume alpha are ordered small to large
      delta <- matrix(0,dim(uout)[1],dim(uout)[2])
      dk <- c(log(k[1]),diff(log(k)))
      rnum <- matrix(runif(length(u)),dim(uout)[1],dim(uout)[2])

      delta  <- delta + log(scalemat[,,1]) * (rnum < cprob[,,1])

      ## dkmat <- array(rep(dk, times = rep(length(u),length(k))),
      ##                  dim = c(dim(uout),length(k)))

      for (iv in 2:length(k)){
        delta <- delta + log(scalemat[,,iv]) * (rnum > cprob[,,iv-1]) *
          (rnum < cprob[,,iv])
      }


      ## for (iv in 1:length(k)){
      ##     delta <- delta + dk[iv] * (rnum < prob[,,iv])
      ## }

      ## dkmat <- array(rep(dk, times = rep(length(u),length(k))),
      ##                  dim = c(dim(uout),length(k)))


    }

    return(matrix(delta,dim(uout)[1],dim(uout)[2]))
  }





drawt <-
  function(uout, alpha = 4, beta = alpha){

    ## draws igamma weights, for model with T errors
    ## beta is a rate, not a scale

    ## Karthik Sastry
    ## R 3.1.2, 64 Bit
    ## August 2016

    ## END PREAMBLE


    u <- c(uout) ## treat everything symmetrically

    ## delta <- rgamma(rep(alpha + 1/2,length(u)), 1/(1/beta + .5 * u^2)) ## this gives 1/sigma^2w


    palpha <- rep(alpha + 1/2,length(u))
    pbeta <-  beta + .5 * u^2

    delta <- rgamma(n = length(palpha), shape = palpha, rate = pbeta) ##1/variance units

    delta <- -log(delta)/2 ## log of standard deviation units
    return(matrix(delta,dim(uout)[1],dim(uout)[2]))
  }




fd <-
  function(delta,tparam,tscale = 1){ ## igamma tparam/2, 2/tparam

    ## returns inverse gamma density evaluated at values of matrix delta

    alpha <- tparam/2
    beta <- tparam/2 * tscale^2 ## scale, not rate

    cterm <- alpha * log(beta) - lgamma(alpha)
    xterm <- -(alpha + 1) * (2 * delta) - beta / exp(2 *  delta)

    return(sum(xterm + cterm))
  }



GsnMove <-
  function(lhfcn, x0, lh0, Sigma, modeA, modeLmd, lc_A0, lc_lmd,
           verbose = FALSE, model = NULL, ordercheck = FALSE,...){

    ## One iteration of an MCMC run, using a Gaussian transition density
    ## Based off code in Chris Sims SVN directory

    ## Contains additional options to check the permutation of shocks relative
    ## to some baseline. We don't use these features in our final analysis in the paper

    ## END PREAMBLE



    x1 <- x0 + MASS::mvrnorm(n = 1, mu = rep(0, length(x0)), Sigma = Sigma)

    newmodel <- lhfcn(x1, lc_A0 = lc_A0, lc_lmd = lc_lmd, ..., verbose = TRUE)
    lh1 <- newmodel$lh



    ## trans <- ((lh0 - lh1 > log(runif(1))) && !attr(lh1, "bad"))
    trans <- (lh0 - lh1 > log(runif(1)))

    ## this logic could be better, but trying to avoid checking for reorderings if transition not accepted
    if (trans) {
      if (ordercheck) {
        mxOutput <- GetALmd(x1)
        perm <- normAlmd(mxOutput$A, mxOutput$lmd, modeA, modeLmd)
        x1 <- GetX(perm$A, perm$lmd, lcA, lc_lmd)

        reordering <- perm$noloop
      } else {
        reordering <- FALSE
      }
      lh0 <- lh1
      x0 <- x1
      model <- newmodel
    }

    if (!verbose){
      return(list(x = x0, lh = lh0, trans = trans ))
    } else {
      return(list(x = x0, lh = lh0, trans = trans, aplusmode = model$vout$var$By,
                  xxi = model$vout$var$xxi, u = model$vout$var$u,model=model))
    }

  }


normAlmd <- function(Aml, lmdml, A, lmd) {
  if(is.null(dim(lmd))) lmd <- matrix(lmd, length(lmd), 1)
  if(is.null(dim(lmdml))) lmdml <- matrix(lmdml, length(lmdml), 1)
  nsig <- dim(lmd)[2]
  nv <- dim(lmd)[1]
  ## normalize diagonal of A, just in case
  sf <- diag(A)
  A <- (1/sf) * A
  lmd <- lmd - 2 * c(log(abs(sf)))        #vector of log sf's gets reused, col by col
  Alml <- array(0, c(nv, nv, nsig))
  Al <- Alml
  for (il in 1:nsig) {
    Alml[ , , il] <- exp(-.5 * lmdml[ , il]) * Aml
    Al[ , , il] <- exp(-.5 * lmd[ , il]) * A
  }
  Alml <- matrix(Alml, nv)
  Al <- matrix(Al, nv)
  xp <- abs(Al %*% t(Alml))
  xp <- log(xp) #better approach to avoiding zeros on diagonal, orthogonal rows, etc.
  xpo <- xp #to check later if trace actually increased
  ## xp <- abs(cor(t(Al), t(Alml)))
  ## Algorithm tries reordering up to nv times to find an invariant ordering,
  ## then gives up and returns nv'th reordering and noloop=FALSE
  ordrng <- 1:nv
  crit <- vector("numeric", nv)
  noloop <- 0
  for (ntrial in 1:nv) {
    thisOrdrng <- 1:nv
    ## Make any switch with 1 that increases trace(xp), then any with 2, etc.
    for (iv in 1:nv) {
      for (iv2 in iv:nv) {
        crit[iv2] <- xp[iv2,iv] - xp[iv,iv] + xp[iv,iv2] - xp[iv2,iv2]
      }
      idtr <- which.max(crit[iv:nv])
      newiv <- thisOrdrng[iv:nv][idtr]
      thisOrdrng[iv:nv][idtr] <- thisOrdrng[iv]
      thisOrdrng[iv] <- newiv
      Al <- Al[thisOrdrng, ]
      xp <- xp[thisOrdrng, ]
    }
    ordrng <- ordrng[thisOrdrng]
    if (all(thisOrdrng == 1:nv)) {
      noloop <- ntrial
      break
    }
  }

  trimprove <- sum(diag(xp)) >= sum(diag(xpo))


  A <- A[ordrng, ]
  sf <- diag(A)

  badshift <- any(sf < .Machine$double.eps)

  if (badshift){
    A <- Aml
    lmd <- lmdml
  } else {
    A <- (1/sf) * A
    lmd <- lmd[ordrng, ] - 2 *c(log(abs(sf)))
  }

  changes <- any(ordrng != 1:nv)
  return(list(Anormed=A , lmdnormed=lmd, ordrng=ordrng, noloop=noloop, badshift = badshift, changes = changes, trimprove = trimprove))
}
