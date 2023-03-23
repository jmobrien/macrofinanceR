# These functions do not appear to have a role in package processes


# From bpss_tools -----


logmean <-
  function(x){

    ## Little function to calculate harmonic means of really big or really small log quantities

    ## Karthik Sastry
    ## R 3.1.2, 64 Bit
    ## August 2016

    ## END PREAMBLE

    lm <- rep(0,4)

    ## 1. simple
    lm[1] <- log(mean(exp(x)))

    ## 2. demeaned
    mn <- mean(x)
    xdm <- x - mn
    lm[2] <- (mn) + log(mean(exp(xdm)))

    ## 3. trimmed
    lm[3] <- log(mean(exp(x), trim = 0.2))

    ## 4. trimmed dm
    lm[4] <- (mn) + log(mean(exp(xdm), trim = 0.2))

    return(lm)
  }

msvreg <- function(dmat,
                   hd = 3,
                   hi = 3,
                   rel = TRUE,
                   ratio = FALSE,
                   nlag = 0,
                   verbose = FALSE,
                   ngdp = NULL){

  ## Runs single equation OLS inspired by Mian, Sufi, and Verner [and other
  ## papers in the "projection regression" literature]
  ## y_{t + hd} - y_{t} = lhs variable
  ## (c/y)_{t-1} - (c/y)_{t-1-hi} = main rhs variable
  ## for each ilag, add y_{t - ilag} - y_{t - ilag - hd}

  ## -------------------- INPUTS --------------------
  ## dmat : matrix of dependent variables. It's assumed that
  ##        column 1 is REAL output
  ##        columns 2 to N are the credit variables
  ## ngdp: specify this as a Tobs x 1 vector of nominal gdp, or log real + log price level
  ## hd and hi: size of difference for dependent and independent variables respectively
  ## rel : if TRUE, take independent variable relative to GDP
  ## ratio: if TRUE, take ind. variable as ratio to GDP
  ## nlag: number of lags of real GDP to add to right hand side
  ## verbose : if true, return (fake) data


  Tobs <- dim(dmat)[1]

  ## preparing the dep. var
  depvar <- matrix(NA,Tobs,1)
  depvar[1:(Tobs-hd)] <- dmat[(hd+1):Tobs,1] - dmat[1:(Tobs-hd),1]



  ## preparing the independent variables
  ncredit <- dim(dmat)[2] - 1
  indvar <- matrix(NA,Tobs,ncredit)

  creditvar <- dmat[,-1,drop=FALSE] ## this works no matter how many there are

  if (is.null(ngdp)){
    ## scale with real
    scalegdp <- dmat[,1]
  } else {
    scalegdp <- ngdp
  }

  if (rel){
    ## relative to GDP
    ## log ratio is one option
    creditvar <- creditvar - scalegdp

    if (ratio){
      ## but can also use the absolute ratio
      creditvar <- exp(creditvar)
    }

  }

  indvar[(hi+2):Tobs,] =
    as.matrix(creditvar[(hi+1):(Tobs-1),] - ## t-1
                creditvar[1:(Tobs-1-hi),]) ## t-1-hi
  sd <- c(sd(indvar[,1],na.rm=TRUE),sd(indvar[,2],na.rm=TRUE))

  ## adding lags if desired
  if (nlag > 0){
    lagmat <- matrix(NA,Tobs,nlag)
    diffY <- c(0,diff(c(dmat[,1])))
    for (ilag in 1:nlag){
      ## lagmat[(hd +ilag):Tobs,ilag]  <- depvar[1:(Tobs-hd-ilag)]
      lagmat[(1 +ilag):Tobs,ilag]  <- diffY[1:(Tobs-ilag)]
    }
    indvar <- cbind(indvar,lagmat)
  }

  ## running the regression
  regout <- lm(depvar ~ indvar)
  regout$sd <- sd


  if (!verbose){
    return(regout)
  } else {
    return(list(regout = regout, depvar = depvar, indvar = indvar))
  }
}



restrictVAR <-
  function(vout, type=c("3", "KF","SVhtskd"), rmat=NULL, yzrone=NULL, xzrone=NULL,
           const=NULL, cyzr=NULL, cxzr=NULL) {
    ## restrictions can be specified as rows of rmat, with coefficients applied to elements of By and Bx
    ## stacked as they are in xxi (and then repeated across the equation index), or they can be specified
    ## in yzrone, xzrone.  Each zero element of yzrone or xzrone generates a restriction that sets the corresponding
    ## coefficient in By or Bx to zero (or to a constant, if !is.null(const)).  Both kinds of restrictions
    ## can be non-trivial in the same call.

    ## type:     vout as from rfvar3 ("3") or as from rfvarKF ("KF")
    ## const:    the right hand side of rmat %*% coeff = const, not the constant in the var.
    ## cyzr, cxzr:  If using yzrone, xzrone with non-trivial constants, leave const=NULL and specify
    ##           constants with cyzr and cxzr
    ## sc:       The Schwarz criterion rejects the restriction if the chisq value plus the sc value
    ##           is positive.  This version of the sc is scale-sensitive.  Variables with higher
    ##           variance are penalized more strongly, as with a prior that expects higher-variance
    ##           variables to explain more variance.
    ##
    ## Note 2013-3-4:  Try eliminating scale effects by converting X'X to correlation matrix
    if (length(type) > 1) type <- type[1]
    if (type == "SVhtskd") {
      bvw <- vout
      vout <- bvw$vout$var
    }
    ncf <- dim(vout$By)[2] * dim(vout$By)[3] + dim(vout$Bx)[2]
    neq <- dim(vout$By)[1]
    ny <- dim(vout$By)[2]
    lags <- dim(vout$By)[3]
    nx <- dim(vout$Bx)[2]
    if (is.null(rmat)) {
      rmat <- matrix(0, 0, ncf *neq)
    }
    if (!is.null(yzrone)) {
      byz <- which(yzrone == 0, arr.ind=TRUE)
      nrstr <- dim(byz)[1]
      if (is.null( cyzr)) cyzr <- array(0, dim(yzrone))
      for (ir in 1:nrstr ) {
        newrow <- rep(0, neq * ncf)
        newrow[(byz[ir,1] - 1) * ncf + (byz[ir, 3] -1) * ny + byz[ir, 2]] <- 1
        rmat <- rbind(rmat,newrow)
      }
      const <- c(const, cyzr[byz])
    }
    if (!is.null(xzrone)) {
      bxz <- which(xzrone == 0, arr.ind=TRUE )
      nrstr <- dim(bxz)[1]
      if (is.null(cxzr)) cxzr <- matrix(0, neq, nx)
      for (ir in 1:nrstr)  {
        newrow <- rep(0,ncf * neq)
        newrow[(bxz[ir,1] - 1) * ncf + ny * lags + bxz[ir, 2]] <- 1
        rmat <- rbind(rmat, newrow)
      }
      const <- c(const, cxzr[bxz])
    }
    svdr <- svd(rmat)
    if (max(abs(svdr$d)) > 1e10 * min(abs(svdr$d))){
      error("restrictions not full rank")
    }
    ## Note that t(rv) spans the same space as rmat, so the restrictiosn are crossprod(v,coeffs)=gamma
    ## rv <- svdr$v    #2013.5.9
    if (length(type) > 1) type <- type[1]
    Tobs <- if (type == "3" || type == "SVhtskd") dim(vout$u)[1] else dim(vout$fcsterr)[1]
    if (type == "3") {
      sig <- cov(vout$u)
      svdsig <- svd(sig)
      singsig <- (max(svdsig$d) > 1e10 * min(svdsig$d))
      if(singsig) warning("Near singular sig matrix in restrictVAR")
      svdxxi <- svd(vout$xxi)
      singxxi <- (max(svdxxi$d) > 1e10 * min(svdxxi$d))
      ## if(singxxi) warning("Near singular xxi matrix in restrictVAR")
      ## schwarz <- rmat %*% kronecker(svdsig$u %*% diag(1/sqrt(svdsig$d)), svdxxi$u %*% diag(1/sqrt(svdxxi$d)))
      ##schwarz <- kronecker((1/sqrt(svdsig$d)) * t(svdsig$u), (1/sqrt(svdxxi$d)) * t(svdxxi$u)) %*% rv  #2013.5.9
      ## sqrtVb <- kronecker(sqrt(svdsig$d) * t(svdsig$u), 1/sqrt(svdxxi$d)) * t(svdxxi$u)
      ## line above seems to be a mistake, since xxi is already x'x-inverse
      sqrtVb <- kronecker(sqrt(svdsig$d) * t(svdsig$u), sqrt(svdxxi$d) * t(svdxxi$u))
      dgVb <- apply(sqrtVb^2, 2, sum)
      rmatC <- rmat %*% diag(sqrt(Tobs * dgVb))
      sqrtVbC <- sqrtVb %*% diag(1/sqrt(Tobs * dgVb))
      lndetVb <- sum(log(svdsig$d)) * dim(vout$xxi)[1] + sum(log(svdxxi$d)) * dim(sig)[1]
      lndetVbC <- lndetVb - sum(log(dgVb * Tobs))
    } else if (type == "KF") {          #type=="KF"
      svdVb <- svd(vout$Vb)
      sqrtVb <- sqrt(diag(svdVb$d)) %*% t(svdVb$u)
      dgVb <- diag(vout$Vb)
      rmatC <- rmat %*% diag(sqrt(Tobs * dgVb))
      sqrtVbC <- sqrtVb %*% diag(1/sqrt(Tobs * dgVb))
      lndetVb <- sum(log(svdVb$d))
      lndetVbC <- lndetVb - sum(log(dgVb * Tobs))
      ## schwarz <- rmat %*% svdVb$u %*% diag(1/sqrt(svdVb$d)) #below is more efficient version for large Vb
      ## schwarz <- (1/sqrt(svdVb$d)) * (t(svdVb$u) %*% rv)
    } else {                            #type="SVhtskd"
      nv <- dim(vout$By)[1]
      nX <- dim(vout$xxi)[1]
      if (is.null(nX)) nX = dim(vout$xx)[1]
      Vb <- matrix(0, nX * nv, nX * nv)


      #check if xxi is missing
      if (is.null(vout$xxi)){
        #instead we have xx, which can be inverted to get xxi
        vout$xxi = vout$xx
        for (iv in 1:nv){
          vout$xxi[,,iv] = solve(vout$xx[,,iv])
        }
      }

      for (iq in 1:nv) {
        Vb[nX * (iq-1) + 1:nX, nX * (iq-1) + 1:nX] <- vout$xxi[ , , iq]
      }

      A0i <- solve(bvw$A)
      Vb <- kronecker(A0i, diag(nX)) %*% Vb %*% kronecker(t(A0i), diag(nX))
      svdVb <- svd(Vb)

      sqrtVb <- sqrt(diag(svdVb$d)) %*% t(svdVb$u)
      # KS note: is this correct?
      sqrtVb2 = svdVb$u %*% sqrt(diag(svdVb$d)) %*% t(svdVb$u)
      dgVb <- diag(Vb)
      rmatC <- rmat %*% diag(sqrt(Tobs * dgVb))
      sqrtVbC <- sqrtVb %*% diag(1/sqrt(Tobs *dgVb))

      lndetVb <- sum(log(svdVb$d))
      lndetVbC <- lndetVb - sum(log(dgVb * Tobs))
    }

    vr2 <- (sqrtVb2 %*% t(rmat))
    cvr2 <- crossprod(vr2)
    svdvr2 <- svd(vr2)

    svdvr <- svd(sqrtVb %*% t(rmat))
    svdvrC <- svd(sqrtVbC %*% t(rmatC)) #result == line above?
    vdim1 <- dim(svdvr$u)[1]
    svdvrp <- svd(diag(vdim1) - svdvr$u %*% t(svdvr$u), nu=vdim1 - dim(rmat)[1])
    svdvrpC <- svd(diag(vdim1) - svdvrC$u %*% t(svdvrC$u), nu=vdim1 - dim(rmat)[1])
    svdvrpuv <- svd(crossprod(svdvrp$u, t(sqrtVb)))
    svdvrpuvC <- svd(crossprod(svdvrpC$u, t(sqrtVbC)))
    lndetUR <- sum(log(svdvrpuv$d))
    lndetURC <- sum(log(svdvrpuvC$d))
    df <- dim(rmat)[1]
    ## schwarz <- -2 * sum(log(diag(chol(crossprod(schwarz)))))   +  df * log(2 * pi)
    schwarz <- lndetVb - 2 * lndetUR + df * log(2 * pi)
    schwarzC <- lndetVbC - 2 * lndetURC + df * log(2 * pi)
    if(is.null(const)) const <- rep(0, dim(rmat)[1])


    if(type == "SVhtskd") {
      vout$By <- tensor::tensor(A0i, vout$By, 2, 1)
      vout$Bx <- A0i %*% vout$Bx
    }
    stackedcf <- c(t(cbind(matrix(vout$By, nrow=neq), vout$Bx)))
    gap <- rmat %*% stackedcf - const
    ##svdv <- svd(rmat %*% vout$Vb %*% t(rmat))
    chstat <- (1/svdvr$d) * t(svdvr$v) %*%  gap
    chstat <- crossprod(chstat)

    #trying with other sqrt matrix -- doesn't seem to make a difference?
    chstat4 = (1/svdvr2$d) * t(svdvr2$v) %*% gap
    chstat4 = crossprod(chstat4)

    ## # alternate method?
    ## # seems to work better
    iR <- apply(rmat == 1, 2, sum)
    rcoefs <- stackedcf
    lcR <- (iR == 1)
    rcoefs[!lcR] <- 0
    ## # rcoefs is zero for everything unrestricted, equal to coefficient for everything restricted
    chstat2 <- t(rcoefs) %*% solve(Vb) %*% rcoefs

    ## ## another alternate method : looking at implied difference in structural form coefficients
    ##    newrfc <- stackedcf
    ##    newrfc[lcR] <- 0
    ##    stackedmat <- t(cbind(matrix(vout$By, nrow = neq), vout$Bx))
    ##    newstmat <- array(newrfc, dim = dim(stackedmat))
    ##    By <- array(c(t(newstmat)), dim  = dim(vout$By))

    ##    for (i in dim(By)[3]){
    ##        By[,,i] <- bvw$A %*% By[,,i]
    ##    }

    ##    listXxi <- lapply(1:dim(vout$xxi)[3], function(i) vout$xxi[,,i])
    ##    covmat <- matrix(bdiag(listXxi), nrow = length(listXxi) * dim(listXxi[[1]])[1])

    ##    oldstack <- t(cbind(matrix(oldBy, nrow=neq), oldBx))
    ##    newstack <- t(cbind(matrix(By, nrow=neq), vout$Bx))

    ##    difference <- oldstack - newstack

    ##    chstat3 <- t(c(difference)) %*% covmat %*% c(difference)
    chstat3 <- NULL

    return(list(chiSquared=chstat, chstat2 = chstat2, chstat3 = chstat3,
                chstat4= chstat4, df=df, sc=schwarz, pval=pchisq(chstat,df),
                sc2 = schwarz - (ncf*neq-df)*log(1 - df/(neq*ncf)), scC=schwarzC))
  }


rfrun <-
  function(ydata,
           xdata = NULL,
           const = TRUE,
           lags = 8,
           nstep = 16,
           mnprior = list(tight = 3, decay = .5),
           urprior = list(lambda = 5, mu = 1),
           vprior = list(sigma = rep(.01,nv), w = 1)){

    ## runs reduced form regressions
    ## this specifies prior "externally" instead of in rfvar3

    ## -------------------- INPUTS --------------------
    ## ydata = main data for VAR
    ## xdata = any extra data for x side
    ## const = if TRUE, add constant to xdata
    ## lags = number of lags in VAR
    ## nstep = number of steps of IR to calculate
    ## mnprior, urprior, vprior: parameters for var prior


    ## so R doesnt complain later
    vnames <- colnames(ydata)
    ydata <- as.matrix(ydata)
    nv <- dim(ydata)[2]

    ## for the unit root prior
    ybar <- apply(ydata[1:lags, ], 2, mean)

    Tobs <- dim(ydata)[1]
    nv <- dim(ydata)[2]
    if (const) {
      xdata <- cbind(xdata, matrix(1,Tobs,1))
    }
    if (!is.null(xdata) ) stopifnot( dim(xdata)[1] == Tobs)
    Tx <- dim(xdata)[1]
    nx <- dim(xdata)[2]


    ## prior
    vp <- varprior(nv,nx,lags,mnprior,vprior, urprior=urprior, ybar=ybar) # vp$: ydum,xdum,pbreaks
    varp <- rfvar3(ydata = vp$ydum,
                   lags = lags,
                   xdata = vp$xdum,
                   lambda = NULL, mu = NULL, ic = NULL, const = FALSE)

    ## posterior mode of standard model
    var <- rfvar3(ydata = rbind(ydata, vp$ydum),
                  lags = lags,
                  xdata = rbind(xdata, vp$xdum),
                  lambda = NULL, mu = NULL, ic = NULL, const = FALSE,
                  breaks = matrix(c(dim(ydata)[1], dim(ydata)[1]+ vp$pbreaks),ncol=1))

    ## p. mode of model with "flat prior"
    vnp <- rfvar3(ydata = rbind(ydata),
                  lags = lags,
                  xdata = rbind(xdata),
                  lambda = NULL, mu = NULL, ic = NULL, const = FALSE)

    ## p. mode of model with only unit root prior
    vp2 <- rfvar3(ydata = rbind(ydata),
                  lags = lags,
                  xdata = rbind(xdata),
                  lambda = 5, mu = 1, ic = NULL, const = FALSE)


    ## impulse responses
    irmn <- impulsdtrf(vout = var, nstep = nstep) ## standard model
    irout <- impulsdtrf(vout = vp2, nstep = nstep) ## with just unit root prior
    irnp <- impulsdtrf(vout = vnp, nstep = nstep) ## no prior IR


    ## marginal likelihood (code is copied from a different program)
    Tu <- dim(var$u)[1]
    Tup <- dim(varp$u)[1]
    flat <- FALSE
    w <- matrictint(crossprod(var$u),var$xxi,Tu-flat*(nv+1))-flat*.5*nv*(nv+1)*log(2*pi);
    wp <- matrictint(crossprod(varp$u),varp$xxi,Tup-flat*(nv+1)/2)-flat*.5*nv*(nv+1)*log(2*pi)
    w=w-wp

    return(list(w = w, vp2 = vp2, var = var, varp = varp, vnp = vnp, irmn = irmn,
                irout = irnp, irp = irout, vnames = vnames))
  }

scaleplots <-
  function(ir1, ir2,
           ird1 = NULL,
           ird2 = NULL,
           shockscale = rep(0,10),
           filename = 'impulse_plot', format = 'pdf',
           shockvars = NULL, responsevars = NULL,
           varnames = rep('',dim(ir)[1]),
           color1 = c(0, 0, 1), gr = 0.7,
           color2 = c(1,0,0),
           width = 5, height = 8,
           conf = .68, nSteps = 60,
           ymat = NULL,
           logseries = rep(0,10),
           shocknames = NULL,
           addtitle = FALSE)
  {

    ## Plots all shocks/responses for an SVAR

    ## Code is specific to one (not very general) organization of input -- key things to preserve in a more
    ## general edit would preserve how the plots look and change how the data is processed

    ## Karthik Sastry
    ## R 3.0.2, 64 Bit
    ## First edit, June 2014


    ## END Preamble


    #### filename <- paste('/home/karthik/R/summer/macrofin/plotting','/plots/', filename,sep='')

    filename <- paste('plots/', filename,sep='')
    if (format == 'pdf'){
      filename <- paste(filename, '.pdf', sep = '')
      grDevices::pdf(filename, width = width, height = height)
      color1 <- grDevices::rgb(color1[1], color1[2], color1[3], 1)
      color2 <- grDevices::rgb(color2[1], color2[2], color2[3], 1)
    } else if (format == 'eps'){
      filename <- paste(filename, '.eps', sep = '')
      lattice::trellis.device(device="postscript", color = TRUE)
      ##setEPS()
      grDevices::postscript(filename, width = width, height = height)

      ## Cannot do transparent colors, so here is a crude workaround
      alphafy <- function(col,alpha=1) {
        rr <- 1-alpha*(1-c(col/255))
        return(grDevices::rgb(rr[1],rr[2],rr[3]))
      }
      color <- alphafy(color, alpha)

    } ##else raise some kind of error?


    ## Determining format of output ----

    ## Defaults to all structural shocks, all responses
    nVar <- dim(ir)[1]

    if (is.null(shockvars)) shockvars <- 1:nVar
    if (is.null(responsevars)) responsevars <- 1:nVar

    nShockvars <- length(shockvars)
    nRespvars <- length(responsevars)

    if (is.null(varnames)){
      varnames <- rownames(ir)
    }

    gRange <- 1:(nSteps) ##graph range
    ##In a model sense, you are only going to nSteps - 1 periods (index 1 is really period 0)

    ## Converting probability "range" into quantiles
    probs <- c(rev((1 - conf) / 2), .5 + conf/2)


    ## Scale every shock correctly
    for (ishock in shockvars){
      sval <- ir1[shockscale[ishock],
                  ishock,1]

      ir2[,ishock,] <- ir2[,ishock,] * sval / ir2[shockscale[ishock],ishock,1]
      ird1[,ishock,,] <- ird1[,ishock,,] * sval / rep(
        ird1[shockscale[ishock],ishock,1,], each = length(ird1[,ishock,,1]))

      ird2[,ishock,,] <- ird2[,ishock,,] * sval / rep(
        ird2[shockscale[ishock],ishock,1,], each = length(ird2[,ishock,,1]))
    }



    ## Plot arrangement ----

    arr <- c(nRespvars, nShockvars)

    ## par(mfrow = arr,col.lab="black",col.main="black",
    ##     oma=c(1,3,1,2), mar=c(.5,.25,.5,.25), tcl=-0.1, mgp=c(0,0,0))

    par(mfrow = arr,
        col.lab="black",col.main="black",
        oma=c(1,5,1,2), mar=c(.5,.25,.5,.25), tcl=-0.1, mgp=c(3,1,0))

    for (irsvar in responsevars) {

      irsseries_1 <- ir1[irsvar,,]
      irstrials_1 <- ird1[irsvar,,,]

      irsseries_2 <- ir2[irsvar,,]
      irstrials_2 <- ird2[irsvar,,,]


      ## Getting error bands, if appropriate
      if (!is.null(irstrials_1) & length(irstrials_1> 0)){
        irsbounds_1 <- array(NA, c(nSteps, length(probs), 10))
        for (iShock in shockvars){
          iShocktrials <- irstrials_1[iShock,,]
          iShockbounds <- t(apply(iShocktrials,1,quantile,probs = probs))
          irsbounds_1[,,iShock] <- iShockbounds[1:(nSteps),]
        }

        irsbounds_2 <- array(NA, c(nSteps, length(probs), 10))
        for (iShock in shockvars){
          iShocktrials <- irstrials_2[iShock,,]
          iShockbounds <- t(apply(iShocktrials,1,quantile,probs = probs))
          irsbounds_2[,,iShock] <- iShockbounds[1:(nSteps),]
        }

      } else {
        irsbounds <- NULL
      }

      ##Determining plot scaling/ axis size
      ##A few different reasonable approaches...

      ##yMax <- max(irsseries, irsbounds, 0) ##+ .1e-5 * abs(max(irsseries, irsbounds))
      ##yMin <- min(irsseries, irsbounds, 0) ##- 1e-5 * abs(min(irsseries, irsbounds))

      if (is.null(ymat)){
        yMax <- max(irsseries_1[shockvars,], irsbounds_1[1:nSteps,,shockvars], irsseries_2[shockvars,], irsbounds_2[1:nSteps,,shockvars], 0) ##+ .1e-5 * abs(max(irsseries, irsbounds))
        yMin <- min(irsseries_1[shockvars,], irsbounds_1[1:nSteps,,shockvars], irsseries_2[shockvars,], irsbounds_2[1:nSteps,,shockvars], 0) ##+ .1e-5 * abs(max(irsseries, irsbounds))
      } else {
        yMin <- ymat[1,irsvar]
        yMax <- ymat[2,irsvar]
      }


      for (iShockvar in shockvars){

        ##Plotting each series

        plottitle <- paste(varnames[irsvar], ' (log)'[logseries[irsvar]],
                           sep = '')

        ## if (length(probs) > 2){ #### 2 sets
        ##     upper <- irsbounds[,4,iShockvar]
        ##     lower <- irsbounds[,1,iShockvar]

        ##     p84 <- irsbounds[,3,iShockvar]
        ##     p16 <- irsbounds[,2,iShockvar]

        ## } else {
        upper1 <- irsbounds_1[,3,iShockvar]
        lower1 <- irsbounds_1[,2,iShockvar]

        upper2 <- irsbounds_2[,3,iShockvar]
        lower2 <- irsbounds_2[,2,iShockvar]

        ## p16 <- NULL
        ## p84 <- NULL
        ## }


        plot(irsseries_1[iShockvar, 1:nSteps], ylim = c(yMin, yMax), type = 'l', lwd = 1,
             xlab = '', ylab = '', yaxt = 'n', xaxt = 'n',
             ## fg = gray(gr),
             xlim = c(1,nSteps),
             xaxs = 'i',
             col = color1)
        lines(irsseries_2[iShockvar, 1:nSteps], lwd = 1, col = color2)

        abline(a = 0, b = 0, lwd = 0.75)

        lines(upper1, lwd = 1, col = color1, lty = 2)
        lines(lower1, lwd = 1, col = color1, lty = 2)

        lines(upper2, lwd = 1, col = color2, lty = 2)
        lines(lower2, lwd = 1, col = color2, lty = 2)


        ##Adding variable name and axis on leftmost plot
        if (which(shockvars == iShockvar, arr.ind = TRUE) == 1) {
          ## mtext(plottitle, side = 2, line = 2, cex = .5)
          mtext(plottitle, side = 2, line = 5, cex = 0.5, las = 1, adj = 0)
          axis(side = 2, cex.axis = .75, las = 1)
        }
        ##Right side axis labels on right most plot
        ## if (which(shockvars == iShockvar, arr.ind = TRUE) == nShockvars) {
        ##     axis(side = 4, cex.axis = .5, fg = gray(gr))
        ## }

        ##Shock name if appropriate
        if (!is.null(shocknames) && (which(responsevars == irsvar, arr.ind = TRUE) == 1)) {
          mtext(shocknames[iShockvar], side = 3, line = 0, cex = .5)
        }
      }
    }

    ## Titles ----

    if (addtitle){
      bigtitle <- paste(type, 'over', as.character(nSteps), 'periods', sep = ' ')
      title(bigtitle, outer = TRUE, cex = 1.2)
    }

    dev.off()
    ##dev.copy2pdf(file = filename)
  }

#' Priors for time series regression
#'
#'     Creates dummy observations for a prior symmetric in variables,
#'     favoring persistence.
#'
#'     The regression right-hand-side is assumed to contain lagged values, possibly
#'     of both dependent and independent variables.  By default it centers on a
#'     random walk with no effects of independent variables. The results can be used
#'     directly in \code{\link{blm}}
#'
#'    @param  vlist Character vector of names of the rhs variables
#'    @param  dvname Name of dependent variable
#'    @param  lags How many lags of each variable.  Can be a single integer,
#'            if all variables appear with the same number of lags, or a
#'            vector of integers, one for each variable. Note that this
#'            program does not care whether, when \code{lags[iv]=3}, this means that
#'            lags 1 to 3, 0 to 2, or 10 to 12 are included.
#'    @param  scale A vector, of the length of \code{vlist} plus 1, of
#'            reasonable guesses as to the standard deviations of
#'            the variables.  The last element scales the levels dummy that relates the constant to the
#'            other parameters, which should be 1.0, or smaller if guesses of variable means in vmeans are
#'            quite uncertain.
#'   @param   bbar mean of the regression coefficient vector.  Default of 1
#'            followed by zeros is reasonable when the first variable is a lagged dependent variable and this
#'            is a regression with a persistent dependent variable.  Last element
#'            is prior mean of constant (usually 0).
#'   @param   smooth The rate at which the scale of dummy observations drops as we go from low to
#'            high frequencies.  Should usually be in (0,1).  If close to zero, prior on the
#'            variable's effect is stronger at low than at high frequencies.  With smooth=1
#'            and damp=1, the prior just shrinks all coefficients toward their prior means
#'            symmetrically and independently.
#'   @param   damp The rate at which the dummy observations increase with increasing lag length. Should
#'            equal or exceed one if distant lags are more likely to be large than near-in lags.
#'   @param   erratio ratio of error variance from parameter uncertainty to error variance from the
#'            residual for the first sample observation. Keeping this constant across models of
#'            different size avoids a prior that implies big models should fit much better.
#'   @param   vmeans A priori means for the variables.  They are used in forming the last dummy observation,
#'            connecting the constant to the other coefficients.  Could be set as sample means of initial conditions.
#'            Using full sample means is problematic: the contamination of prior by data may not be a big
#'            problem for stationary data, but could be substantially distorting for non-stationary data.
#'   @return  A list with components:
#'   \itemize{
#'     \item y. Dependent variable value for dummy observations.
#'     \item X. Dummy observations implementing the prior.
#'     \item scalefac. Amount by which the prior has been scaled to match \code{erratio}.
#'   }
#'   @details Note that the last column of \code{X} is the
#'            "constant", which will *not* be a column of ones.  Generating the sample
#'            \code{X} matrix using \code{lagts()} will not by itself create a constant term in the last column.
#'            \code{lagts()} does put the variables and lags in the same order as this function.  But you have
#'            to tack on the constant vector with \code{cbind()}.  \code{X} is a matrix with column names
#'            constructed from \code{vnames} (repeated for lags).  Note that \code{lagts()} allows non-sequential lists
#'            of lags.  This program allows different lag lengths, but always \code{1:lag[iv]}.
#'  @seealso \code{\link{lagts}} for preparation of a data matrix in a form that is compatible with this
#'           function.  \code{\link{blm}} to combine the prior with data to obtain estimates.
tsregPrior <- function(vlist, dvname="y", lags=rep(0, length(vlist)), ldv=1, scale,
                       bbar=NULL, smooth, damp, vmeans, erratio) {
  nv <- length(vlist)
  nx <- if (length(lags) > 1) {
    sum(lags)
  } else {
    nv * lags
  }
  if (length(lags) == 1) {
    lags <- rep(lags, nv)
  }
  X <- matrix(0, nx + 1, nx + 1)
  nv <- length(vlist)
  js <- 0
  for (iv in 1:nv) {
    X[js + (1:lags[iv]), js + (1:lags[iv])] <- smooth^(0:(lags[iv]-1)) * ctmat(lags[iv]) %*% diag(damp^(0:(lags[iv]-1))) * scale[iv]
    js <- js + lags[iv]
  }
  X[ ,  nx+1] <- 0
  vmeans <- c(rep(vmeans, times=lags), 1)
  X[nx+1, ] <- vmeans * scale[nv + 1]
  erratio0 <- solve(t(X), vmeans)
  erratio0 <- sum(erratio0)^2
  X <- sqrt(erratio0/erratio) * X
  w <- -.5 * determinant(crossprod(X))$modulus
  if (is.null(bbar)) bbar <- c(1, rep(0,nx))
  y <- X %*% bbar
  y <- matrix(y, ncol=1, dimnames=list(NULL, dvname))
  dimnames(X)[[2]] <- c(rep(vlist, times=lags), "const")           # no indication of lags in names
  return(list(y=y, X=X, scalefac=sqrt(erratio0/erratio), call=match.call()))
}


# From opt_tools: ----

numHess <-
  function(fcn, x, ...) {
    f1 <- fcn
    n <- length(x)
    h <- matrix(0, n, n)
    f2 <- function(z, ...) { numgrad(fcn=f1, z, ...)$g}
    h <- numgrad(fcn=f2, x=x, ...)$g
    return(h)
  }


# From var_tools: ----

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
