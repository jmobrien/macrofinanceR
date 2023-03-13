impulsdtrf <-
  function(vout=NULL, smat=NULL, nstep=40, order=NULL)
    ### vout:           output structure from rfvar3.
    ###                 To use this with output from postdraw, create a dummy vout with vout$By=pout$By[ , , ,id] and provide smat=pout$smat[ , ,id]
    ### smat:           if !is.null(vout) and order and smat are NULL, the impulse responses will be for a cholesky decomp with variables
    ###                 ordered as in the input to rfvar3.  More generally, can be any set of initial
    ###                 values for the shocks.  To keep the shocks orthogonal and scaled properly,
    ###                 smat should be such that smat %*% t(smat) == crossprod(vout$u)/dim(u)[1].
    ###                 However, the routine works with singular smat or with smat's column dimension
    ###                 less than its row dimension.
    ### order:          To get a cholesky decomp with a different ordering, set order to an integer
    ###                 vector giving the desired ordering.
    ### response:       nvar x nshocks x nstep array of impulse responses.
###
###                 with vout from rfvarKF, smat argument is required, since there is no vout$u.
###
### Code written by Christopher Sims, based on 6/03 matlab code.  This version 3/27/04.
### Added dimension labeling, 8/02/04.  Allow non-square smat, integrate with rfvar3 output, 4.7.10.
  {
    B <- vout$By
    neq <- dim(B)[1]
    nvar <- dim(B)[2]
    lags <- dim(B)[3]
    dimnB <- dimnames(B)
    if (is.null(smat)) {
      if (is.null(order) ) {
        order <- 1:neq
      }
      smat <- t(pchol(crossprod(vout$u)/dim(vout$u)[1], order)) # makes first shock affect all variables
    }
    nshock <- dim(smat)[2]
    if(dim(smat)[1] != dim(B)[1]) stop("B and smat conflict on # of equations") #
    response <- array(0,dim=c(neq,nshock,nstep+lags-1));
    response[ , , lags] <- smat
    response <- aperm(response, c(1,3,2))
    irhs <- 1:(lags*nvar)
    ilhs <- lags * nvar + (1:nvar)
    response <- matrix(response, ncol=nshock)
    B <- B[, , seq(from=lags, to=1, by=-1)] #reverse time index to allow matrix mult instead of loop
    B <- matrix(B,nrow=nvar)
    for (it in 1:(nstep-1)) {
      response[ilhs, ] <- B %*% response[irhs, ]
      irhs <- irhs + nvar
      ilhs <- ilhs + nvar
    }
    ## for (it in 2:nstep)
    ##       {
    ##         for (ilag in 1:min(lags,it-1))
    ##           response[,,it] <- response[,,it]+B[,,ilag] %*% response[,,it-ilag]
    ##       }
    dim(response) <- c(nvar, nstep + lags - 1, nshock)
    response <- aperm(response[ , -(1:(lags-1)), ,drop=FALSE], c(1, 3, 2)) #drop the zero initial conditions; array in usual format
    dimnames(response) <- list(dimnB[[1]], dimnames(smat)[[2]], NULL)
    ## dimnames(response)[2] <- dimnames(smat)[1]
    ## dimnames(response)[1] <- dimnames(B)[2]
    return(response)
  }


kfVC <-
  function(y, X, shat, sig, M) {
    ## s is the state, and the plant equation is s(t)=G %*% s(t-1)+t(M) %*% e, where e is
    ## N(0,I).  The observation equation is y(t)=Hs(t).  The prior distribution for
    ## s is N(shat,sig).  To handle the standard Kalman Filter setup where the observation
    ## equation is y(t)=Hs(t)+Nu(t), expand s to become [s;v] (where v=Nu), expand H to [H I], replace
    ## G by [G 0;0 0], and replace M with [M 0;0 N].  The resulting model has no "error" in the
    ## observation equation but is equivalent to the original model with error in the state equation.
    ## The posterior on the new state is N(shatnew,signew) and lh is a two-dimensional vector containing
    ## the increments to the two component terms of the log likelihood function.  They are added
    ## to form the log likelihood, but are used separately in constructing a concentrated or marginal
    ## likelihood. fcsterr is the error in the forecast of y based on the prior.
    ## ....
    ## This version specializes to the case of a constant-coefficient VAR, where H= cbind(kronecker(I, X[it, ], I)
    ## and G is the identity (for the constant coefficients) with zeros (for the shocks in the state vector)
    ## appended in the lower right.  Also, kf2's M is all zeros except for the lower right nvar x nvar.  By
    ## bypassing lots of multiplication by zero, this version is faster than generic kf2by a factor of four for
    ## a 7-variable, 13-lag VAR.  Here M is the
    ## transposed cholesky factor of the VAR shock covariances, not the full state equation M (whichis full of
    ## zeros).  With M constant, plain rfvar3 is more efficient.
    ## ....
    SMALLSV <- 1e-7
    nv <- length(y)
    nx <- length(X)
    nXX <- nv*nx
    nstate <- length(shat)
    omega <- matrix(0, nstate, nstate)
    omega[1:nXX, 1:nXX] <- sig[1:nXX,1:nXX]
    sigObs <- crossprod(M)
    omega[nXX + (1:nv), nXX + (1:nv)] <- sigObs
    nstate <- length(shat)
    ## stopifnot (nstate >= nobs)
    ##------------ Don't need separate treatment of H == 0.  H %*% G %*% t(H) = 0 covers it.
    ##   if (isTRUE(all.equal(H, 0))) { # No observation case.  Just propagate the state.
    ##     lh <- c(0,0)
    ##     shatnew <- G %*% shat
    ##     signew <- omega
    ##     fcsterr <- y                        # y had better be 0
    ##     if (!all.equal(y,0) ) warning("zero H but non-zero y")
    ##   } else {
    ## ho <- H %*% omega
    ho <- matrix(0, nv, nstate)
    for (iv in 1:nv)
      ho[iv, ] <- c(X %*% omega[ ((iv - 1) * nx) + 1:nx, 1:nXX],   sigObs[iv, ])
    hoh <- matrix(0, nv, nv)
    for (jv in 1:iv)
      hoh[, jv] <- ho[ , nx * (jv-1) + (1:nx)] %*% X + sigObs[ , jv]
    svdhoh <- svd( hoh )
    shatmat <- matrix(shat[1:(nv * nx)], nv, nx, byrow=TRUE)
    if (all(svdhoh$d < SMALLSV)) { # Observation is uninformative. Propagate state.
      lh <- c(0,0)
      shatnew <- shat
      signew <- omega
      fcsterr <- y - shatmat %*% X     # had better be 0
      if (!all(abs(fcsterr) < 1e-7)) warning("Uninformative H but non-zero fcsterr")
    } else {
      first0 <- match(TRUE, svdhoh$d < SMALLSV)
      if (is.na(first0)) first0 <- nv + 1
      u <- svdhoh$u[ , 1:(first0-1), drop=FALSE]
      v <- svdhoh$v[ , 1:(first0-1), drop=FALSE]
      d <- svdhoh$d[1:(first0-1), drop=FALSE]
      fcsterr <- y - shatmat %*% X
      hohifac <- (1/sqrt(d)) * t(u)       #diag(1/sqrt(d)) %*% t(u)
      ferr <- hohifac %*% fcsterr
      lh <- c(0,0)
      lh[1] <- -.5 * crossprod(ferr)
      lh[2] <- -.5 * sum( log(d) )
      hohoifac <-hohifac %*% ho
      shatnew <- crossprod(hohoifac, ferr) + c(shat[1:nXX], rep(0,nv))
      signew <- omega - crossprod(hohoifac)
    }
    return(list(shat=shatnew, sig=signew, lh=lh, fcsterr=fcsterr))
  }


kfVCx <-
  function(y, X, shat, sig, M) {
    ## This is for the case of the state  constant, with the observation
    ## equation y = X %*% s + e, sqrt(Var(e))=M.
    SMALLSV <- 1e-30
    if (!is.null(dim(y))) {
      nv <- dim(y)[2]
      nx <- dim(X)[2]
      Tobs <- dim(y)[1]
    } else {
      nv <- length(y)
      nx <- length(x)
      Tobs <- 1
      y <- matrix(Tobs, nv)
      X <- matrix(Tobs, nx)
    }
    nXX <- nv*nx
    nstate <- length(shat)
    omega <- array(sig, c(nx, nv, nx, nv))
    sigObs <- crossprod(M)
    ## ho: kron(I,X) %*% sighat
    ## ho <- array(0, c(Tobs, nv, nx, nv))
    ## sv <- array(0, c(Tobs, nv, Tobs, nv))
    ##cat("\nstart double X prod loop", proc.time())
    ## for (iv in 1:nv) {
    ##     ho[ , iv, , ] <- tensor::tensor(X, omega[ , iv, ,  ], 2, 1)
    ##     for (jv in 1:nv)
    ##         sv[ , iv, , jv] <- tensor::tensor(ho[ , iv, , jv ], X, 2, 2) + sigObs[iv, jv] * diag(Tobs)
    ## }
    ho <- tensor::tensor(X, omega, 2, 1)  #Tobs, nv, nx, nv
    sv <- tensor::tensor(ho, X, 3, 2)     #Tobs, nv, nv, Tobs
    sv <- aperm(sv, c(1,2,4,3))
    for (irow in 1:nv)
      for (icol in 1:nv)
        sv[ , irow, , icol] <- sv[ , irow, , icol] + sigObs[irow, icol] * diag(Tobs)
    ##svdsv <- svd( matrix(sv, nv * Tobs, nv * Tobs ))
    ##cat("\nstart qr", proc.time)
    qrsv <- qr( matrix(sv, nv * Tobs, nv * Tobs ), LAPACK=TRUE)
    rankr <- match(TRUE, abs(diag(qrsv$qr)) < SMALLSV)
    if (is.na(rankr)) {
      rankr <- nv * Tobs
    } else {
      rankr <- rankr - 1
    }
    shatmat <- matrix(shat,nx, nv )
    fcsterr <- y - X %*% shatmat
    lh <- c(0,0)
    ## we're relying on qr( ,LAPACK=TRUE) always putting zeros in diagonal of R on bottom.
    ## Note that we have to account for possible pivot
    ## in QR
    if (rankr == 0) { # Observation is uninformative. Propagate state.
      shatnew <- shat
      signew <- omega
      if (!all(abs(fcsterr) < 1e-7)) warning("Uninformative H but non-zero fcsterr")
    } else if (rankr < nv *Tobs) {
      R <- qr.R(qrsv)[1:rankr, ]
      d <- diag(R)
      R <- R[ , sort.list(qrsv$pivot)]
      Q <- qr.Q(qrsv)[ , 1:rankr]
      RQ <- R %*% Q
      oie <- Q %*% solve(RQ, crossprod(Q, c(fcsterr)))
      lh[1] <- -.5 * c(fcsterr) %*% oie
      lh[2] <- -.5 * sum(log(abs(d)))
      ho <- matrix(ho, Tobs*nv, nx*nv)
      shatnew <- shat + oie %*% ho   # works because oie is a vector, not a matrix
      qho <- crossprod(Q, ho)
      signew <- sig - crossprod(qho , solve(RQ, qho))
    } else {
      ## first0 <- match(TRUE, svdsv$d < SMALLSV)
      ## if (is.na(first0)) first0 <- nv * Tobs + 1
      ## u <- svdsv$u[ , 1:(first0-1), drop=FALSE]
      ## v <- svdsv$v[ , 1:(first0-1), drop=FALSE]
      ## d <- svdsv$d[1:(first0-1), drop=FALSE]
      ## svifac <- (1/sqrt(d)) * t(u)   #diag(1/sqrt(d)) %*% t(u)
      ## ferr <- svifac %*% c(fcsterr)
      ##lh[1] <- -.5 * crossprod(ferr)
      oie <- solve(qrsv, c(fcsterr))
      lh[1] <- -.5 * c(fcsterr) %*% oie
      lh[2] <- -.5 * sum( log(abs(diag(qr.R(qrsv)))) )
      ## svoifac <-svifac %*% matrix(ho, nv * Tobs, nx * nv)
      ## shatnew <- crossprod(svoifac, ferr) + shat
      ## signew <- sig - crossprod(svoifac)
      ho <- matrix(ho, Tobs*nv, nx*nv)
      shatnew <- shat + oie %*% ho
      signew <- sig - crossprod(ho, solve(qrsv, ho))
    }
    return(list(shat=shatnew, sig=signew, lh=lh, fcsterr=fcsterr))
  }


matrictint <-
  function(S,XXi,Tobs,cx=NULL,cs=NULL)
    ###  S:  usually sample cross product matrix of LS residuals
    ###  cs: instead of S, can provide UT cs s.t. cs'cs==S, as from chol(S) or qr.R(qr(S))
    ### XXi:  inv(X'X) matrix for rhs variables
    ### cx:  instead of XXi, can provide UT cx
    ###  Tobs:  number of observations
    ###  w:  log of integrated posterior for SUR or RF VAR with det(Sigma)^(-(m+1)/2) Jeffreys-like prior
    ###  To get the log of the integral of the likelihood for a VAR with Tobs observations,
    ###   k rhs variables in each equation, and m equations, set Tobs=Tobs-m-1 and subtract .5*m*(m+1)*log(2*pi).
    ### We are integrating the exponential of -.5*Tobs*m*log(2*pi)-.5*(Tobs+m+1)*log(det(Sigma))-.5*trace(Sigma\S(beta)).
  { if(is.null(cx))
  {
    ## browser()
    k<-dim(XXi)[1]
    ##cx <- chol(XXi)
    cx <- try(chol(XXi))
    if(inherits(cx,"try-error")) stop("XXI not p.d.")
  }
    else
    {
      k <- dim(cx)[1]
    }
    if(is.null(cs))
    {
      m<-dim(S)[1]
      ##cs <- chol(S)
      cs<-try(chol(S));
      if(inherits(cs,"try-error")) stop("S not p.d.")
    }
    else
    {
      m <- dim(cs)[1]
    }

    w<-(-Tobs+k+(m-1)/2)*m*.5*log(pi)-(Tobs-k)*sum(log(diag(cs)))+m*sum(log(diag(cx)))+ggammaln(m,(Tobs-k)/2)
    return(w)
  }

ggammaln <-
  function (m,ndf)
    ###function gg<-ggamma(m,ndf)
    ### From 8.2.22 on p.427 of Box and Tiao, this is the log of generalized
    ### gamma divided by gamma(.5)^(.5*m*(m-1))
  {
    if( ndf<=(m-1)/2)
      stop('too few df in ggammaln')
    else
      ##lgg=.5*m*(m-1)*gammaln(.5); % normalizing factor not used in Wishart integral
      garg<-ndf+.5*seq(0,1-m,by=-1)
    lgg<-sum(lgamma(garg))
    return(lgg)
  }


mgnldnsty <-
  function(ydata,lags,xdata=NULL, const=TRUE, breaks=NULL,lambda=5,mu=1,mnprior=list(tight=3,decay=.5),
           vprior=list(sig=NULL,w=1),train=0,flat=FALSE,nonorm=FALSE,ic=NULL)
    ### ydata:        endogenous variable data matrix, including initial condition dates.
    ### xdata:        exogenous variable data matrix, including initial condition dates.
    ### const:        Constant term is added automatically if const=TRUE.
    ### breaks:       breaks in the data.  The first lags data points after a break are used
    ###               as new initial conditions, not data points for the fit.
    ### lambda:       weight on the co-persistence prior dummy observation.  (5 is reasonable)
    ###               lambda>0 => x variables included; lambda<0 => x variables excluded;
    ### mu:           weight on variable-by-variable sum of coeffs dummy obs. (1 is reasonable)
    ### mnprior$tight:weight on the Minnesota prior dummies.  Prior std dev on first lag is
    ###               1/mnprior$tight
    ### mnprior$decay:prior std dev on own lag j is 1/j^decay
### vprior$sig:   vector of nv prior std dev''s of equation shocks.  vprior$sig is needed
###               to scale other components of the prior, even if vprior$w=0. Not needed for a pure training
###               sample prior.
### vprior$w:     weight on vcv dummies.  (1 is reasonable; higher values tighten up.)
### train:        If non-zero, this is the point in the sample at which the
###               "training sample" ends.  Prior x likelihood to this point is weighted to
###               integrate to 1, and therefore is treated as if it were itself the prior.
###               To do a pure training sample prior, set lambda=mu=0, mnprior=NULL, vprior$w=0,
###               train>lags.
### flat:         Even with lambda=mu=vprior$w=0, mnprior=NULL, det(Sigma)^(-(nv+1)/2) is used
###               as a "prior", unless flat=TRUE. flat=TRUE is likely not to work unless train is reasonably large.
### nonorm:       If true, use dummy observations but do not normalize posterior to make them a
###               proper prior.  Useful to duplicate results obtained by others, to use
###               dummy observations that do not imply a proper prior, or to save computing time in case only the
###               posterior on this model's parameters, not the weight on the model, is needed.
### ic:           Initial conditions matrix for use in forming the sums of coefficients dummy observations.
###               If ic=NULL, the means of the first lags observations in ydata are used.  If !is.null(ic),
###               ic should be a single "observation" on the y's and x's that will be used as the persistent
###               values entering the sums of coefficients dummies.
###
###               Note that to enter a prior directly as dummy observations, one can treat the
###               Dummy observations as a training sample.
###
  {
    if (is.null(dim(ydata)))  ydata <- matrix(ydata, ncol=1)
    Tobs <- dim(ydata)[1]
    nv <- dim(ydata)[2]
    if (const) {
      xdata <- cbind(xdata, matrix(1,Tobs,1))
    }
    ## looks likely that const=FALSE, xdata=NULL case crashes.  (2012.9.24)
    if (!is.null(xdata) ) stopifnot( dim(xdata)[1] == Tobs)
    Tx <- dim(xdata)[1]
    nx <- dim(xdata)[2]
    ## 2013.8 fix:  added urprior here, set lambda and mu to NULL in rfvar3 call, so
    ## prior dummies treated correctly in normalizing the prior.
    if (is.null(ic)) {
      ybar <- apply(ydata[1:lags, , drop=FALSE], 2, mean)
    } else {
      ybar <- ic
    }
    vp <- varprior(nv,nx,lags,mnprior,vprior, urprior=list(lambda=lambda, mu=mu), ybar=ybar)
    ## vp$: ydum,xdum,pbreaks
    var <- rfvar3(ydata=rbind(ydata, vp$ydum), lags=lags, xdata=rbind(xdata,vp$xdum), breaks=matrix(c(breaks, Tobs, Tobs + vp$pbreaks), ncol=1),
                  const=FALSE, lambda=NULL, mu=NULL, ic=ic) # const is FALSE in this call because ones alread put into xdata
    Tu <- dim(var$u)[1]
    if ( var$snglty > 0 ) error( var$snglty, " redundant columns in rhs matrix")
    w <- matrictint(crossprod(var$u),var$xxi,Tu-flat*(nv+1))-flat*.5*nv*(nv+1)*log(2*pi);
    if(train!=0) {
      if(train <= lags) {
        cat("end of training sample <= # of lags\n")  #
        return
      }
      Tp <- train
      tbreaks <- c(breaks[breaks<train],Tp)
    } else {
      Tp <- lags
      ## because need initial conditions to form lambda/mu prior dummy obs
      tbreaks <- Tp
    }
    ytrain <- ydata[1:Tp,,drop=FALSE]
    xtrain <- xdata[1:Tp,,drop=FALSE]
    if (!nonorm) {
      ## fixed 2013.8.14:  Looks as if dummy obs from urprior are missed here.  Should include
      ## non-null lambda, mu in call to varprior, not in rfvar3 call.
      ## varp <- rfvar3(ydata=rbind(ytrain, vp$ydum), lags=lags, xdata=rbind(xtrain, vp$xdum),
      ##               breaks=c(tbreaks, Tp+vp$pbreaks), lambda=lambda, mu=mu, const=FALSE, ic=ic)
      varp <- rfvar3(ydata=rbind(ytrain, vp$ydum), lags=lags, xdata=rbind(xtrain, vp$xdum),
                     breaks=c(tbreaks, Tp+vp$pbreaks), lambda=NULL, mu=NULL, const=FALSE, ic=ic)
      ## const is FALSE here because xdata already has a column of ones.
      if (varp$snglty > 0) {
        warning("Prior improper, short ", varp$snglty, " df.  Results likely nonsense.")
      } else {
        Tup <- dim(varp$u)[1]
        wp <- matrictint(crossprod(varp$u),varp$xxi,Tup-flat*(nv+1)/2)-flat*.5*nv*(nv+1)*log(2*pi)
        w=w-wp
      }
    } else {
      varp <- NULL
    }
    return(list(w=w,var=var,varp=varp,prior=list(lambda=lambda,mu=mu,vprior=vprior,mnprior=mnprior), call=match.call()))
  }


pchol <-
  function(sig, porder) {
    invporder <- match(1:length(porder), porder)
    return(chol(sig[porder, porder])[invporder, invporder])
  }


postdraw <-
  function(vout,n,nosigprior=FALSE){
    ## Provides n draws from the joint posterior on residual covariance and coefficients of a
    ## reduced form VAR.
    ## vout:	return value from rfvar3()
    ## n:		Number of draws
    ## nosigprior:	If true, don't use the Jeffreys-motivated improper prior.
    ## smat:        t(chol(sigma draw)).  Use directly as smat in impulsdtrf.
    ## By,Bx:       coefficient draws
    ## 11/25/09 Bugnote:  If is.null(vout$Bx) (no constant, no exog vbles), code below
    ## doesn't work.  Need to fix.
    xxi <- chol(vout$xxi)
    ncf <- dim(xxi)[1]
    S <- crossprod(vout$u)
    df <- dim(vout$u)[1]-dim(xxi)[1]      # This effectively uses a |Sigma|^{-(nvar+1)/2} prior
    neq <- dim(vout$u)[2]
    lags <- (ncf-dim(vout$Bx)[2])/neq
    if(nosigprior){df <- df-dim(S)[1]-1}	# This undoes the |Sigma|^{-(nvar+1)/2} prior
    wmat <- rwwish(df,solve(S),n)
    for (it in 1:n){wmat[,,it] <- chol(solve(wmat[,,it]))}
    nmat <- array(rnorm(n*neq*ncf),c(ncf,neq,n))
    cfmat <- t(cbind(matrix(vout$By,neq,neq*lags),vout$Bx))
    for (ir in 1:n){
      nmat[,,ir] <- crossprod(xxi,nmat[,,ir]) %*% wmat[,,ir]+cfmat
    }
    Byx <- aperm(nmat, c(2,1,3))
    By <- Byx[ , 1:(neq*lags), ]
    dim(By) <- c(neq,neq,lags,n)
    ## Bx <- as.vector(vout$Bx)+aperm(nmat,c(2,1,3))[,(neq*lags+1):ncf,]  # Bx added in both here and in cfmat. Bug caught by A.Zawadwoski
    Bx <- Byx[,(neq*lags+1):ncf,]
    ## Note that if reordered as below, wmat[,,i] would be  ready for use as input to impulsdtrf.
    ## but as is, needs transposition.
    ## wmat <- aperm(wmat, c(2,1,3))
    return(list(By=By,Bx=Bx,smat=wmat))
  }



rfvar3 <-
  function(ydata=NA,n_lags=6,xdata=NULL,const=TRUE,breaks=NULL,
           lambda=5,mu=2,ic=NULL, sigpar=NULL, # cores = 1,
           oweights = NULL, drawbe = FALSE) {
    #### This algorithm goes for accuracy without worrying about memory requirements.
    ####
    #### The standard prior it implements is NOT APPROPRIATE for seasonally unadjused data, even
    #### if seasonal dummies are included in xdata.  The prior shrinks toward simple persistence, so it
    #### will tend to prevent the dummies from picking up all the seasonality.
    ####
    #### ydata:   T x nvar dependent variable data matrix.
    #### xdata:   T x nx exogenous variable data matrix.
    ####          Note that if either ydata or xdata has only one column, it must still have a dim vector.  In
    ####          other words it must be a Tx1 array, not a vector of length T.
    ####
    #### const:   If TRUE, a column of ones is added to (or becomes, if xdata is NULL) the xdata matrix.
    #### n_lags:    number of lags
    #### breaks:  rows in ydata and xdata after which there is a break.  This allows for
    ####          discontinuities in the data (e.g. war years) and for the possibility of
    ####          adding dummy observations to implement a prior.  This must be a column vector.
    ####          Note that a single dummy observation becomes lags+1 rows of the data matrix,
    ####          with a break separating it from the rest of the data.  The function treats the
    ####          first lags observations at the top and after each "break" in ydata and xdata as
    ####          initial conditions.
    #### lambda:  weight on "co-persistence" prior dummy observations.  This expresses
    ####          belief that when all variables are at a fixed initial level, they tend to
    ####          stay there.  This is consistent with stationarity and with nonstationarity with or
    ####          without cointegration.  With lambda < 0 , the
    ####          constant term is not included in the dummy observation, so that stationary models
    ####          with means equal to initial ybar do not fit the prior mean.  With lambda>0, the prior
    ####          implies that large constants are unlikely if unit roots are present.  To omit this type of
    ####          dummy observation, use lambda=NULL.
    #### mu:      weight on "own persistence" prior dummy observation.  Expresses belief
    ####          that when y_i has been stable at its initial level, it will tend to persist
    ####          at that level, regardless of the values of other variables.  There is
    ####          one of these for each variable.  A reasonable first guess is mu=2.
    ####          To omit this type of dummy observation, use mu=NULL
    #### ic:      for direct input of the initial conditions mean that is used in the persistence dummy observations,
    ####          as ic$ybar and ic$xbar.
    ####          If is.null(ic), the mean of the first lags observations in ydata, xdata are used.
    ####      The program assumes that the first lags rows of ydata and xdata are real data, not dummies.
    ####      Dummy observations should go at the end, if any.  If pre-sample x's are not available,
    ####      repeating the initial xdata(lags+1,:) row or copying xdata(lags+1:2*lags,:) into
    ####      xdata(1:lags,:) are reasonable subsititutes.  These values are used in forming the
    ####      persistence priors.
    #### sigpar: list(A0, lmd, breaks_pos) Allow SVAR with time varying shock variances.  See below.
    #### returns:
    #### By:      nvar x nvar x lags matrix of coefficients on lagged y's.  1st dimension is "equation number"
    #### Bx:      nvar x nx matrix of coefficients on x's
    #### u:       (T-6+ (number of dummy obs)) x nvar matrix of residuals.  If ydata is a ts object, u will be also, and will
    ####          be correctly dated.  u observations dated after end(ydata) are dummy observations.
    #### xxi:     X'X inverse, same for all equations.  kronecker(cov(u),xxi) is the full covariance matrix of the regression coefficients.
    #### snglty:  Usually 0.  If the rhs variable matrix is not full column rank, this is the gap between the number of columns and the
    ####          number of non-zero singular values.
    #### Code written by Christopher Sims.  This version 8/13/04.
    #### 12/18/05:  added ts properties for u, better comments.
    ####
    #### Modified 2013.8.12 to allow use of A0, lmd, breaks_pos.  With non-null A0, By is A+ from
    #### A0 %*% y(t) = A+(L) %*% y(t) + exp(.5 lmd(t)) * eps(t) .  This works even with
    #### lmd constant, but in that case running a single rf estimate (A0=I), then iterating
    #### on (A0, lmd) alone makes more sense. With lmd varying, rf estimates change with lmd.
    ####
    if (is.null(dim(ydata))) {
      dim(ydata) <- c(length(ydata),1)
    }

    Tobs <- dim(ydata)[1]
    #### Note that if rfvar3() has been called with dummy obs's already in place, this Tobs
    #### includes the dummies.
    nvar <- dim(ydata)[2]
    ####nox <- isempty(xdata)
    if (const) {
      xdata <- cbind(xdata,matrix(1,Tobs,1))
    }
    nox <- identical(xdata,NULL)
    if(!nox){
      Tobs2 <- dim(xdata)[1]
      nx <- dim(xdata)[2]
    } else {
      Tobs2 <- Tobs; nx <- 0; xdata <- matrix(0,Tobs2,0)
    }
    #### note that x must be same length as y, even though first part of x will not be used.
    #### This is so that the lags parameter can be changed without reshaping the xdata matrix.
    ####
    if (!identical(Tobs2,Tobs)) {
      print('Mismatch of x and y data lengths')
      return()
    }
    if (identical(breaks,NULL))
      nbreaks <- 0
    else {
      #### if (is.ts(ydata)) {                ## Can use Yr, month-or-quarter pairs, or real number dates.
      ####   if (is.matrix(breaks) ) {
      ####     breaks <- breaks[ , 1] + (breaks[ ,2] - 1) / frequency(ydata)
      ####   } else {
      ####     if (any(abs(breaks - round(breaks))) > 1e-8) {
      ####       breaks <- match(breaks, time(ydata))
      ####     }
      ####   }                               ##if not real numbers, not yr-month pairs, it's just obs number
      #### }
      #### Any use of tsp(ydata) has to be in external processing functions.
      nbreaks <- length(breaks)
    }

    breaks <- c(0,breaks,Tobs)
    if(any(breaks[2:length(breaks)] < breaks[1:(length(breaks)-1)]))
      stop("list of breaks must be in increasing order\n")
    #### initialize smpl as null if initial observations are only there for lambda/mu prior.
    #### matlab code uses the fact that in matlab a:b is null if b<a, which is not true for R.
    #### if(breaks[2]>lags)
    ####   smpl <- (lags+1):breaks[2]
    #### else
    ####   smpl <- NULL
    #### if(nbreaks>0){
    ####   for (nb in 2:(nbreaks+1))
    ####     smpl <- c(smpl,(breaks[nb]+lags+1):breaks[nb+1])
    #### }
    smpl <- NULL
    for (nb in 2:(nbreaks + 2)) {
      if ( breaks[nb] > breaks[nb-1] + n_lags )
        smpl <- c(smpl, (breaks[nb-1] + n_lags + 1):breaks[nb])
    }
    #### With logic above, one can use an mts-type ydata and omit sections of it by including sequences of breaks separated by
    #### less than lags+1.  E.g. with lags=6, monthly data, breaks=rbind(c(1979,8), c(1980,2), c(1980,8), c(1980,12)) omits
    #### Sep 1979 through Dec 1981, plus 6 months after that, which are initial conditions for the next sample segment.
    Tsmpl <- length(smpl)
    X <- array(0,dim=c(Tsmpl,nvar,n_lags))
    for(ix in seq(along=smpl))
      X[ix,,] <- t(ydata[smpl[ix]-(1:n_lags),,drop=FALSE])
    dim(X) <- c(Tsmpl,nvar*n_lags)
    X <- cbind(X, xdata[smpl,,drop=FALSE])
    y <- ydata[smpl,,drop=FALSE]
    #### Everything now set up with input data for y=Xb+e
    #### ------------------Form persistence dummies-------------------
    if (! (is.null(lambda) & is.null(mu) ) ) {
      if(is.null(ic)) {
        ybar <- apply(as.array(ydata[1:n_lags,,drop=FALSE]),2,mean)
        dim(ybar) <- c(1,dim(ydata)[2])
        if (!nox) {
          xbar <- apply(array(xdata[1:n_lags,,drop=FALSE],dim=c(n_lags,dim(xdata)[2])),2,mean)
          dim(xbar) <- c(1,dim(xdata)[2])
        } else {
          xbar <- NULL
        }
      } else {
        ybar <- ic$ybar
        xbar <- ic$xbar
      }
      if (!is.null(lambda)){
        if (lambda<0){
          lambda <- -lambda
          xbar <- array(0,c(1,dim(xdata)[2]))
        }
        xdum <- lambda * cbind(array(rep(ybar,n_lags),dim=c(1,n_lags*length(ybar))), xbar)
        ydum <- array(0,c(1,nvar))
        ydum[1,] <- lambda*ybar
        y <- rbind(y,ydum)
        X <- rbind(X,xdum)
      }
      if (!is.null(mu)) {
        xdum <- cbind(
          array(rep(diag(as.vector(ybar),nrow=length(ybar)),n_lags),
                dim=c(dim(ybar)[2],dim(ybar)[2]*n_lags)),
          array(0,dim=c(nvar,dim(xdata)[2])))*mu
        ydum <- mu*diag(as.vector(ybar),nrow=length(ybar))
        X <- rbind(X,xdum)
        y <- rbind(y,ydum)
      }
    }
    if (!is.null(sigpar)) {
      breaks_pos <- sigpar$breaks_pos
      lmd <- sigpar$lmd
      A0 <- sigpar$A0
      if (!is.null(breaks_pos)) {
        #### breaks_pos <- invtime(breaks_pos, ydata) ##so breaks_pos given as dates
        nsig <- length(breaks_pos)
      } else {
        nsig <- 1
      }
      breaks_pos <- c(breaks_pos, Tobs)
      lmdndx <- rep(1:nsig, times=diff(breaks_pos))
      lmdseries <- lmd[ , lmdndx]
      if ( Tsmpl < dim(y)[1] ) {      ##dummy obs formed in rfvar3
        #### Should not be combining this branch with dummy obs's from varprior()
        #### already included in ydata.
        lmdp <- apply(lmdseries[ ,smpl], 1, mean)
        lmdseries <- cbind(lmdseries[ , smpl], matrix(lmdp, nvar, dim(y)[1] - Tsmpl))
      } else {
        lmdseries <- lmdseries[ , smpl]
      }

      if (!is.null(oweights)){
        ## add outlier weights
        ## scaled as log standard deviation
        lmdseries[,1:dim(oweights)[1]]  <-
          lmdseries[,1:dim(oweights)[1]] - 2*(t(oweights))
      }

      #### i.e., use mean of lmdseries for dummy observation weights.  Note that
      #### since lmd is logged, this is geometric mean, maybe not best.

      nX <- dim(X)[2]

      ## Apply A0
      if (length(dim(A0)) == 2){      # constant A0
        ya0 <- y %*% t(A0)
      } else {
        ## time-varying A0
        A0list <- A0[,,lmdndx]
        A0list[,,smpl]

        ny <- dim(y)[1]              # number of observations
        ya0 <- array(0, dim(y))
        for (iy in 1:ny){           # apply A0 to each observation
          ya0[iy,] <- y[iy,] %*% t(A0list[,,iy])
        }
      }

      ## Save frequency of each regime
      freqs <- lmdndx[smpl]




      ##B <- matrix(0,  nX, nvar)
      ##u <- matrix(0, Tsmpl, nvar)
      ##uraw <- u
      uraw <- NULL ##why does it matter?
      ##xxi <- array(0, c(nX, nX, nvar))
      ##logdetxxi <- vector("numeric", nvar)
      ##snglty <- vector("numeric", nvar)

      ## Parallel implementation
      ##linreg in a separate file
      ##listOutput <- parallel::mclapply(1:nvar, linreg, lmdseries, X, ya0, mc.cores = cores)
      listOutput <- lapply(1:nvar, linreg, lmdseries, X, ya0, drawbe = drawbe)
      errorflags <- matrix(unlist(sapply(listOutput, '[[', 1)))
      if (any(errorflags)) return(NULL) ##tells program to blow up the likelihood
      B <- matrix(unlist(sapply(listOutput, '[[', 2)), ncol = nvar)
      u <- matrix(unlist(sapply(listOutput, '[[', 3)), ncol = nvar)
      logdetxxi <- matrix(unlist(sapply(listOutput, '[[', 4)))
      snglty <- matrix(unlist(sapply(listOutput, '[[', 5)))
      ## xx, not xxi!
      xx <- array(unlist(sapply(listOutput, '[[', 6)), dim = c(nX, nX, nvar))

      ## draws of coefs and epsilon_it
      Bdraw <- matrix(unlist(sapply(listOutput, '[[', 7)), ncol = nvar)
      udraw <- matrix(unlist(sapply(listOutput, '[[', 8)), ncol = nvar)

      ## B <- matrix(0,nX,nvar)
      ## u <- matrix(0,dim(X)[1],nvar)
      ## logdetxxi <- matrix(0,nvar,1)
      ## snglty <- matrix(0,nvar,1)
      ## xx <- array(0,c(nX,nX,nvar))
      ## for (iv in 1:nvar){
      ##     linout <- linreg(iv,lmdseries,X,ya0)
      ##     B[,iv] <- linout[[2]]
      ##     u[,iv] <- linout[[3]]
      ##     logdetxxi[iv] <- linout[[4]]
      ##     snglty <- linout[[5]]
      ##     xx[,,iv] <- linout[[6]]
      ## }

    } else {
      #### Instead of svd below, could invoke lsfit.  Faster?
      vldvr <- svd(X)
      dfx <- sum(vldvr$d > 100*.Machine$double.eps)
      di <- 1./vldvr$d[1:dfx]
      vldvr$u <- vldvr$u[, 1:dfx]
      vldvr$v <- vldvr$v[, 1:dfx]
      snglty <- dim(X)[2] - dfx
      logdetxxi <- 2 * sum(log(abs(di)))
      ####B <- vldvr$v %*% diag(di,nrow=length(di)) %*% t(vldvr$u) %*% y (line below is just more efficient)
      B <- vldvr$v %*% (di * (t(vldvr$u) %*% y))
      u <-  y-X %*% B
      xxi <-  di * t(vldvr$v)
      xxi <-  crossprod(xxi)
      uraw <- NULL       ## so it won't be missing in the list of outputs
    }
    if (!is.null(tsp(ydata))) u <- ts(u, start=start(ydata)+c(0,n_lags),freq=frequency(ydata))
    #### dates at end of sample are for dummy obs, meaningless.  If there are other
    #### nontrivial breaks, the dates for u are also meaningless.
    #### dim(B) <-  c(nvar*lags+nx,nvar) ## rhs variables, equations (this was redundant)
    By <-  B[1:(nvar*n_lags),]
    dim(By) <-  c(nvar,n_lags,nvar)       ## variables, lags, equations
    By <-  aperm(By,c(3,1,2)) ##equations, variables, lags to match impulsdt.m
    #### label all the output, if the data matrices had labels
    if(!is.null(dimnames(ydata)[2]))
    {
      ynames <- dimnames(ydata)[[2]]
    }else
    {
      ynames <- rep("",times=nvar)
    }
    if(!nox)
    {
      if(!is.null(dimnames(xdata)[[2]]))
      {
        xnames <- dimnames(xdata)[[2]]
      } else {
        xnames <- rep(" ",times=nx)
      }
    }
    dimnames(By) <- list(ynames,ynames,as.character(1:n_lags))
    ##xxinames <- c(paste(rep(ynames,lags),rep(1:lags, each=length(ynames)),sep=""),xnames)
    ##dimnames(xxi) <- list(xxinames,xxinames)
    if (nox)
      Bx <-  NULL
    else
    {
      Bx <-  matrix(B[nvar*n_lags+(1:nx),],dim(B)[2],nx)
      dimnames(Bx) <- list(ynames,xnames)
    }
    ###### logintlh <-  matrictint(u'*u,xxi,size(X,1)-nvar-1)-.5*nvar*(nvar+1)*log(2*pi);
    ###### Might want to create a version without the dimnames if using this in a program.
    ####
    #### returns some things that are not available with sigpar=NULL.  Either split to
    #### separate programs, or create alternate return lists.
    if(!is.null(sigpar)) {
      return(list(By=By, Bx=Bx, u=u, uraw=uraw, xx = xx, snglty=snglty, logdetxxi=logdetxxi,
                  lmdseries=lmdseries,
                  Bdraw = Bdraw, udraw = udraw,
                  freqs = freqs,
                  call=match.call())) ##var.logintlh =  logintlh
    } else {
      return(list(By=By, Bx=Bx, u=u, xxi= xxi, snglty=snglty, logdetxxi=logdetxxi, call=match.call()))
    }
  }


SVARlh <-
  function(A0,sigma,Tobs) {
    ## Calculates log likelihood (or posterior density, if conjugate prior has been used) for an
    ## overidentified structural VAR, assuming restrictions on the contemporaneous coefficient
    ## matrix A0 only.
    ##
    ## Note that determinant() returns the log of the determinant, as a list.
    lh <- -.5 * Tobs * log(2*pi) + Tobs * with(determinant(A0), modulus) - .5 * Tobs * sum(crossprod(A0) * sigma)
    return(lh)
  }


SVARlh0 <- function(pvec, idmat, sigma, Tobs) {
  ## idmat is a logical matrix, TRUE where A0 has a non-zero coefficient
  ## pvec  is the vector of parameter values that fill A0[idmat].
  ## sigma is the covariance matrix of estimated residuals from the reduced
  ##       form VAR
  ## Tobs     is the sample size.
  ## This function returns minus the likelihood, so it can be used directly in csminwel
  n <- dim(idmat)[1] # better be same as dim(idmat[2])
  A0 <- matrix(0,n,n)
  A0[idmat] <- pvec
  lh <- SVARlh(A0, sigma, Tobs)
  Tdum <- 60
  lh <- lh + 2*Tdum*log(abs(A0[1,1])) - Tdum*((.005*A0[1,1] + .003*A0[1,2])^2/2 + .01 * A0[1,5]^2) #prior isn't normalized
  return(-lh)
}


sysmat <- function(By, Bx=NULL) {
  ## Constructs the lags*nv x lags*nv system matrix for the "stacked" first-order
  ## version, from the nv x nv x lags array of coefficients By returned by rfvar3.
  ## If there's a constant (constant terms in Bx from rfvar3), the matrix is expanded
  ## to include the trivial constant dynamics.
  n <- dim(By)
  ny <- n[1]
  nnl <- n[2] * n[3]
  dim(By) <- c(ny, nnl)
  By <- rbind(By,diag(1, nrow = nnl - ny, ncol=nnl))
  if( !is.null(Bx)) {
    By <- cbind(By, c(Bx, rep(0, nnl -ny)))
    By <- rbind(By, c(rep(0, nnl), 1))
  }
  return(By)
}

varprior <-
  function(nv=1,nx=0,lags=1,mn_prior=list(tight=5,decay=.5),
           v_prior=list(sig=1,w=1),
           ur_prior=list(lambda=NULL, mu=NULL), xsig=NULL, ybar=NULL, xbar=1, nstat=rep(TRUE,nv),
           mn_start = 1,
           cos_prior = NULL
  )
    ### ydum, xdum:   dummy observation data that implement the prior
    ### breaks:       vector of points in the dummy data after which new dummy obs start
    ###                   Set breaks=Tobs+matrix(c(0,breaks),ncol=1), ydata=rbind(ydata,ydum), xdum=rbind(xdata,xdum), where
    ###                   actual data matrix has Tobs rows, in preparing input for rfvar3
    ### nv,nx,lags: VAR dimensions
    ### mn_prior$tight:Overall tightness of Minnesota prior. 1/tight ~ own lag std dev
    ### mn_prior$decay:Standard deviations of lags shrink as lag^(-decay)
    ### v_prior$sig:   Vector of prior modes for square roots of diagonal elements of r.f. covariance matrix
    ###                  Names of this vector name columns of output ydum.
    ### v_prior$w:     Weight on prior on vcv.  1 corresponds to "one dummy observation" weight
    ###                   v_prior$sig is needed
###                   to scale the Minnesota prior, even if the prior on sigma is not used itself.
###                   Set v_prior$w=0 to achieve this.
###                   mn_prior and v_prior.w can each be set to NULL, thereby eliminating the corresponding
###                   dummy observations.
### xsig:          rough scale of x variances.  names of this vector name output xdum
### ur_prior:       Parameters of the "unit roots" and "co-persistence" priors that are
###                   implemented directly in rfvar3.  lambda and mu should be NULL here if
###                   the dummy observations generated here are used with rfvar3 and lanbda and mu
###                   are not NULL in rfvar3.   lambda < 0 means x'st not included.  Note that constant
###                   is assumed to be last element of x.  If you want lambda < 0 to be the only source
###                   of a prior on the constant, but xsig is not null, set the last element of xsig
###                   to zero.
### ybar,xbar:        estimates of data means, used in constructing ur_prior component, but not otherwise.
###                   The default xbar=1 is correct when the constant is the only x.
### nstat:         Set components corresponding to non-persistent variables to FALSE.
### Note:          The original Minnesota prior treats own lags asymmetrically, and therefore
###                   cannot be implemented entirely with simple dummy observations.  It is also usually
###                   taken to include the sum-of-coefficients and co-persistence components
###                   that are implemented directly in rfvar3.R.  The diagonal prior on v, combined
###                   with sum-of-coefficients and co-persistence components and with the unit own-first-lag
###                   prior mean generates larger prior variances for own than for cross-effects even in
###                   this formulation, but here there is no way to shrink toward a set of unconstrained
###                   univariate ARs.
###
###
  {
    ## nx=0 case messes up, at least at the end (2012.9.23)
    if (!is.null(mn_prior)) ## implement an MN prior
    { ## single-coefficient prior dummy obs.
      ## each vbl and each lag has a dummy observation, and each dummy obs has values for current and lagged
      ## y's  and current x's. we separate the y's and the x's into two arrays.  The last two indexes, lag
      ## and rhsy, index the dummy observations.
      xdum <- if(nx > 0) {
        array(0, dim=c(lags + 1, nx, lags, nv), dimnames=list(obsno=1:(lags + 1), xvbl=1:nx, lag=1:lags, rhsy=1:nv))
      } else {
        NULL
      }
      ydum <- array(0,dim=c(lags+1,nv,lags,nv),dimnames=list(obsno=1:(lags+1),rhsy=1:nv,lag=1:lags, rhsy=1:nv))
      for (il in 1:lags) {


        #ydum[il + 1,,il,] <- il^mn_prior$decay * diag(v_prior$sig,nv,nv)

        # KS alteration, July 2014: allowing the mn prior to start the decay at the mn_start + 1 lag (default 1, so decay starts with 2)
        decayfactor <- (il <= mn_start) * 1 + ((il > mn_start) * (il - mn_start + 1))^(mn_prior$decay)
        # il <= mn_start gets decayfactor = 1; il > mn_start gets decay normally associated with (il - mn_start + 1)
        # e.g., start prior at ilag = 2; no decay on lag 1, 2 variances; lag 3 variance is treated like lag 2 variance in previous set-up

        ydum[il+1,,il,] <- decayfactor*diag(v_prior$sig,nv,nv)
      }
      ## If we have non-trivial x's, need dobs's for them, also.
      if(!is.null(xsig)) {
        ydumx <-  array(0, dim=c(lags + 1, nv, nx), dimnames=list(obsno=1:(lags + 1), rhsy=1:nv, dx=1:nx))
        xdumx <-  array(0, dim=c(lags + 1, nx, nx), dimnames=list(obsno=1:(lags + 1), xvbl=nx, dx=1:nx))
        xdumx[1, , ] <- diag(xsig, nx, nx)
        ## note that xvalues for obsno 2:(lags+1) don't matter.  This is one dummy obseervation,
        ## so only the "current" x is used.
      }
      ydum[1,,1,] <- diag(v_prior$sig * nstat, nv, nv) # so own lag has mean zero if nstat FALSE
      ydum <- mn_prior$tight * ydum
      dim(ydum) <- c(lags+1,nv,lags*nv)
      ydum <- ydum[seq(lags+1,1,by=-1),,]
      xdum <- mn_prior$tight*xdum
      dim(xdum) <- c(lags+1,nx,lags*nv)
      xdum <- xdum[seq(lags+1,1,by=-1),,]
    } else {
      ydum <- NULL;
      xdum <- NULL;
      breaks <- NULL;
      lbreak <- 0;
    }

    if (!is.null(cos_prior)){ ## Implement the cosine transformation prior

      ## same size arrays as with the MN prior
      xdumc <- if(nx > 0) {
        array(0, dim=c(lags + 1, nx, lags, nv), dimnames=list(obsno=1:(lags + 1), xvbl=1:nx, lag=1:lags, rhsy=1:nv))
      } else {
        NULL
      }


      ydumc <- array(0,dim=c(lags+1,nv,lags,nv),dimnames=list(obsno=1:(lags+1),rhsy=1:nv,lag=1:lags, rhsy=1:nv))

      cosmat <- ctmat(lags) ## cosine transformation matrix, common for all coefficients
      dmat <- cos_prior$smooth^ (0:(lags-1)) * cosmat %*% diag(cos_prior$damp ^ (0:(lags-1))) ## matrix of dummies, except for variance scaling for rhs

      for (iv in 1:nv){ ## iterate over lhs variables
        ydumc[1 + 1:lags, iv, 1:lags, iv] <- t(dmat) * v_prior$sig[iv]
      }

      for (iv in 1:nv){
        ydumc[1,iv,,iv] <- ydumc[2,iv,,iv]
      }

      ydumc <- cos_prior$tight * ydumc
      dim(ydumc) <- c(lags+1,nv,lags*nv)
      ydumc <- ydumc[seq(lags+1,1,by=-1),,]

      ## ydumc <- array(0,c(lags + 1, nv, lags*nv))

      ## for (irhs in 1:nv){ ## each rhs variable
      ##     idmat <- dmat * v_prior$sig[irhs] ## scaled dummies
      ##     xseq <- seq(irhs,irhs + (lags - 1)*nv,nv)
      ##     for (ilhs in 1:nv){
      ##         ydumc[1 + 1:lags,ilhs, xseq] <- idmat
      ##     }

      ##     ydumc[1,irhs,1:nv] <- ydumc[1 + 1:lags, irhs, xseq]
      ## }


      dim(xdumc) <- c(lags+1,nx,lags*nv)
      xdumc <- cos_prior$tight*xdumc
      xdumc <- xdumc[seq(lags+1,1,by=-1),,,drop = FALSE]

    } else {
      ydumc <- NULL
      xdumc <- NULL

    }

    if (!is.null(ur_prior$lambda) ) {
      ## lambda obs.  just one
      ydumur <- matrix(ybar, nrow=lags+1, ncol=nv, byrow=TRUE) * abs(ur_prior$lambda)
      if(ur_prior$lambda > 0) {
        xdumur <- matrix(xbar, lags + 1, nx, byrow=TRUE) * ur_prior$lambda # (all but first row redundant)
      } else {
        xdumur <- matrix(0, lags + 1, nx)
      }
    } else {
      ydumur <- NULL
      xdumur <- NULL
    }
    ## mu obs. sum(nstat) of them
    if (!is.null(ur_prior$mu)) {
      ydumuri <-array(0, c(lags+1, nv, nv))
      for (iv in which(nstat)) {
        ydumuri[ , iv, iv] <- ybar[iv]
      }
      ydumur <- abind::abind(ydumur, ur_prior$mu *ydumuri, along=3)
      xdumur <- abind::abind(xdumur, array(0, c(lags+1, nx, nv)), along=3)
    }
    if (!is.null(v_prior) && v_prior$w > 0)
    {
      ydum2 <- array(0,dim=c(lags+1,nv,nv))
      xdum2 <- array(0,dim=c(lags+1,nx,nv))
      ydum2[lags+1,,] <- diag(v_prior$sig,nv,nv)*v_prior$w #The v_prior$w factor was missing until 11/29/06
      # Original idea, not implemented, was probably that w be an integer
      # repetition count for variance dobs.
      # Now it's just a scale factor for sig. in variance prior.
    } else {
      ydum2 <- NULL
      xdum2 <- NULL
    }
    ## stack everything up.
    if (!is.null(ydum)){
      dim(ydum) <- c(lags + 1, nv, lags * nv) # merge all the individual mn dobs
      dim(xdum) <- c(lags + 1, nx, lags * nv)
    }

    ydum <- abind::abind(ydum, ydumc, ydumur, ydum2, along=3)
    xdum <- abind::abind(xdum, xdumc, xdumur, xdum2, along=3)
    breaks <- (lags+1) * (1:(dim(ydum)[3] -1)) # end of sample is not a "break".
    ydum <- aperm(ydum, c(1, 3, 2))
    ydum <- matrix(ydum, ncol=dim(ydum)[3])
    xdum <- aperm(xdum, c(1,3,2))
    xdum <- matrix(xdum, ncol=dim(xdum)[3])
    ##   dim(ydum2) <- c((lags+1)*nv,nv)
    ##   dim(ydum) <- c((lags+1)*nv,lags*nv)
    ##   ydum <- cbind(ydum,ydum2)
    ##   dim(xdum2) <- c((lags+1)*nx,nv)
    ##   dim(xdum) <- c((lags +1)*nx,lags*nv)
    ##   xdum <- cbind(xdum,xdum2)
    ##   dim(ydum) <- c(lags+1,nv,dim(ydum)[2])
    ##   ydum <- aperm(ydum,c(1,3,2))
    ##   dim(ydum) <- c(dim(ydum)[1]*dim(ydum)[2],nv)
    ##   dim(xdum) <- c(lags+1,nx,dim(xdum)[2])
    ##   xdum <- aperm(xdum,c(1,3,2))
    ##   dim(xdum) <- c(dim(xdum)[1]*dim(xdum)[2],nx)
    ##   if(nv>1){
    ##     breaks <- c(breaks, (lags+1)*(0:(nv-1))+lbreak)
    ##   }
    ## } else {
    ##   if (!is.null(ydum)) { # case with mn_prior non-null, but v_prior null
    ##     ydum <- aperm(ydum, c(1, 3, 2))
    ##     dim(ydum) <- c(prod(dim(ydum)[1:2]), dim(ydum)[3])
    ##     xdum <- aperm(xdum, c(1,3,2))
    ##     dim(xdum) <- c(prod(dim(xdum)[1:2]), dim(xdum)[3])
    ##   }
    ## }
    dimnames(ydum) <- list(NULL, names(v_prior$sig))
    dimnames(xdum) <- list(NULL, names(xsig))
    return(list(ydum=ydum,xdum=xdum,pbreaks=breaks))
    ## data here in the form of Tobs by nv y, and Tobs x nx x.  Lagged y's not put in to a rhs
    ## regression matrix, so a "breaks" vector is needed.
    ## rfvar3 adds persistence and sum of coeffs dummy observations at end of  data in lhs and rhs
    ## regression matrix form.  So to combine this with rfvar3, set lambda and mu to NULL in one or the
    ## other program.
  }


varpriorN <-
  function(nv=1,nx=1,lags=1,mnprior=list(tight=5,decay=.5),vprior=list(sig=1,w=1),
           urprior=list(lambda=NULL, mu=NULL), xsig=NULL, ybar, xbar=1, nstat=rep(TRUE,nv), sigfix=NULL) {
    ## varprior produces dommy observations, interpreted as scaled relative to residual variances.
    ## This function uses the var prior output to produce a mean and covariance matrix for a normal prior.
    if( !is.null(sigfix)) {
      vprior$w <- 0
      vprior$sig <- sqrt(diag(sigfix))
    }
    prior <- varprior(nv, nx, lags, mnprior,vprior, urprior, xsig, ybar, xbar, nstat)
    ##
    ## contruct prior mean and variance from varprior output
    vpout <- rfvar3(prior$ydum, lags=lags, xdata=prior$xdum, const=FALSE,
                    breaks=prior$pbreaks, lambda=NULL, mu=NULL)
    shat <- with(vpout, c(rbind(matrix(aperm(By, c(2,3,1)), prod(dim(By)[2:3]), dim(By)[1]), t(Bx))))
    ## lags,vbls for y, then x, then eq
    ## Indexing in shat: ((vbl x lag), const) x eqn
    if(is.null(sigfix)) {
      sighat <- with(vpout, kronecker(crossprod(u), xxi))
    } else {
      sighat <- kronecker(sigfix, vpout$xxi)
    }
    ##
    ## crossprod(u) rather than var(u) above because in the prior,only nv dummy observations contribute
    ## to the variance prior.  The others "fit perfectly" and should not contribute to the sigma prior.
    ##
    ## scale by vprior$sig, since this is an absolute normal prior, not dummy observations that would be
    ## implicitly scaled by equation variances.
    ##
    ## This rescaling looks wrong to me today  (2012.9.24)
    ## wtvec <- rep(vprior$sig, each=nv * lags + nx)
    ## sighat <- wtvec * t(wtvec * sighat)
    return(list(shat=shat, sighat=sighat, call=match.call()))
  }


vdecomp <-
  function(resp) {
    vdc <- apply(resp^2, c(1,2), cumsum)
    for (it in 1:dim(resp)[3]) {
      vdc[it , , ] <- diag(1/apply(vdc[it, , ], 1, sum)) %*% vdc[it , , ]
    }
    vdc <- aperm(vdc, c(2,3,1))
    dimnames(vdc) <- dimnames(resp)
    return(vdc)
  }

