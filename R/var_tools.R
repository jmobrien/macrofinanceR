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

