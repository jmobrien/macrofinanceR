bfgsi <-
  function(H0,dg,dx) {
    ### dg is previous change in gradient; dx is previous change in x;
    ### 6/8/93 version that updates inverse hessian instead of hessian
    ### itself.
    ### Copyright by Christopher Sims 1996.  This material may be freely
    ### reproduced and modified.
    n <- length(dg)
    dim(dg) <- c(n,1)
    dim(dx) <- c(n,1)
    Hdg <- H0 %*% dg
    dgdx <- as.numeric(crossprod(dg,dx))
    dxHdg <- drop(dx %o% Hdg) # drops are needed to get rid of redundant unit-dimensions
    x1 <- as.numeric(1+crossprod(dg,Hdg)/dgdx)*drop(dx %o% dx)
    x2 <- dxHdg+t(dxHdg)
    ## code below causes problems when x1 and x2 have matching zeros, and I can't now (2005-12-22)
    ## figure out why I thought it was a good idea
    ##   if ( max(abs(x1-x2)/(abs(x1)+abs(x2))) <= 1e-12 ) {
    ##     cat("bfgs update failed.\n")
    ##     cat("x1, x2 = ",x1,x2,"\n")
    ##     H=H0
    ##     return(H)
    ##   }
    if (abs(dgdx) <= 1e-12)   {
      cat("bfgs update failed.\n")
      cat("|dg| =", sqrt(sum(dg^2)), "|dx| = ", sqrt(sum(dx^2)),"\n")
      cat("crossprod(dg,dx) =", dgdx,"\n")
      cat("|H*dg| =", crossprod(Hdg),"\n")
      H=H0
      return(H)
    }
    ## otherwise,
    H <- H0 + x1/dgdx - x2/dgdx
    save(file="H.dat", H)
    return(H)
  }


csminit <-
  function(fcn,x0,f0,g0,badg,H0,...){
    ### retcodes: 0, normal step.  5, largest step still improves too fast.
    ### 4,2 back and forth adjustment of stepsize didn't finish.  3, smallest
    ### stepsize still improves too slow.  6, no improvement found.  1, zero
    ### gradient.
    ###---------------------X
    ### Fixed 7/17/93 to use inverse-hessian instead of hessian itself in bfgs
    ### update.
    ###
    ### Fixed 7/19/93 to flip eigenvalues of H to get better performance when
    ### it's not psd.
    ###
    ## ANGLE <- .0005
    ANGLE <- 1e-7
    THETA <- .01 #(0<THETA<.5) THETA near .5 makes long line searches, possibly fewer iterations.
    FCHANGE <- 1000
    MINLAMB <- 1e-9
    ### fixed 7/15/94
    ### MINDX <- .0001;
    ### MINDX <- 1e-6;
    MINDFAC <- .01
    fcount<-0
    lambda<-1
    xhat<-x0
    f<-f0
    fhat<-f0
    g <- g0
    gnorm <- sqrt(sum(g^2))
    ###
    if (!badg && (gnorm < 1.e-12)) {      # put !badg 8/4/94
      retcode <- 1
      dxnorm <- 0
      ## gradient convergence
    } else {
      ## with badg true, we don't try to match rate of improvement to directional
      ## derivative.  We're satisfied just to get some improvement in f.
      ##
      ##if(badg)
      ##   dx = -g*FCHANGE/(gnorm*gnorm);
      ##  dxnorm = norm(dx)
      ##  if dxnorm > 1e12
      ##     disp('Bad, small gradient problem.')
      ##     dx <- dx*FCHANGE/dxnorm;
      ##   end
      ##else
      ## Gauss-Newton step;
      ##---------- Start of 7/19/93 mod ---------------X
      ##[v d] = eig(H0);
      ##toc
      ##d <- max(1e-10,abs(diag(d)));
      ##d <- abs(diag(d));
      ##dx <- -(v.*(ones(size(v,1),1)*d'))*(v'*g);
      dx <- -H0 %*% g
      dxnorm <- sqrt(sum(dx^2))
      if (dxnorm > 1e12){
        cat("Near-singular H problem.\n")
        dx <- dx*FCHANGE/dxnorm
      }
      dfhat <- crossprod(dx,g0)
      if(!badg){
        ## test for alignment of dx with gradient and fix if necessary
        a <- -dfhat/(gnorm*dxnorm)
        if(a<ANGLE){
          if (a < 0) {
            dx <- -g / gnorm^2
            dfhat <- -1
            cat("H unused\n")
            ## Don't bother with H.  It's not psd. Step length here appropriate for log LH,
            ## where 1.0 is a reasonable scale for changes.
          } else {
            dx <- dx - as.numeric(ANGLE*dxnorm/gnorm+dfhat/(gnorm*gnorm))*g
            dx <- dx * dxnorm / sqrt(sum(dx^2)) # This line added 2/18/2004
            dfhat <- crossprod(dx,g)
            ## dxnorm <- sqrt(sum(dx^2)) # No longer needed with 2/18/2004 change
            cat("Correct for low angle:" ,a,"\n")
          }
        }
      }
      cat(sprintf("Predicted improvement            = %18.9f\n",-dfhat/2))
      ## cat("Predicted improvement:", sprintf("%18.9f",-dfhat/2),"\n")
      ##
      ## Have OK dx, now adjust length of step (lambda) until min and
      ## max improvement rate criteria are met.
      done <- 0
      factor <- 3
      shrink <- 1
      lambdaMin <- 0
      lambdaMax <- Inf
      lambdaPeak <- 0
      fPeak <- f0
      lambdahat <- 0
      while(!done){
        ## argument of fcn retains its dim (or lack thereof), but g
        ## always emerges as n x 1
        ddx <- dx*lambda
        dim(ddx) <- dim(x0)
        dxtest <- x0+ddx
        f <- fcn(dxtest,...)
        cat(sprintf("lambda = %10.5g; f = %20.7f",lambda,f ),"\n")

        #test code
        if (is.nan(f)) f <- Inf #this helps with explosive regressions
        #end test code

        if(f < fhat){
          fhat <- f
          xhat <- dxtest
          lambdahat <- lambda
        }
        fcount <- fcount+1
        shrinkSignal <- (!badg && (f0-f < max(-THETA*dfhat*lambda,0))) || (badg && ((f0-f) < 0))
        growSignal <- (!badg && ( (lambda > 0)  &&  (f0-f > -(1-THETA)*dfhat*lambda) ))
        if(  shrinkSignal  &&   ( (lambda>lambdaPeak) || (lambda<0) )){
          if ((lambda>0) && ((!shrink || (lambda/factor <= lambdaPeak)))){
            shrink <- 1
            factor <- factor^.6
            while(lambda/factor <= lambdaPeak){
              factor <- factor^.6
            }
            if( abs(factor-1)<MINDFAC){
              if( abs(lambda) < 4){
                retcode <- 2
              }else{
                retcode <- 7
              }
              done <- 1
            }
          }
          if(  (lambda<lambdaMax) && (lambda>lambdaPeak) ){
            lambdaMax <- lambda
          }
          lambda <- lambda/factor
          if( abs(lambda) < MINLAMB ){
            if( (lambda > 0) & (f0 <= fhat) ){
              ## try going against gradient, which may be inaccurate
              lambda <- -lambda*factor^6
            }else{
              if( lambda < 0 ){
                retcode <- 6
              }else{
                retcode <- 3
              }
              done <- 1
            }
          }
        }else{
          if(  (growSignal && lambda>0) ||  (shrinkSignal && ((lambda <= lambdaPeak) && (lambda>0)))  ) {
            if( shrink ){
              shrink <- 0
              factor <- factor^.6
              if( abs(factor-1)<MINDFAC ) {
                if( abs(lambda)<4 ) {
                  retcode <- 4
                }else{
                  retcode <- 7
                }
                done <- 1
              }
            }
            if( ( f<fPeak ) && (lambda>0) ) {
              fPeak <- f
              lambdaPeak <- lambda
              if( lambdaMax<=lambdaPeak ) {
                lambdaMax <- lambdaPeak*factor*factor
              }
            }
            lambda <- lambda*factor
            if( abs(lambda) > 1e20 ) {
              retcode <- 5
              done <-1
            }
          } else {
            done <- 1
            if( factor < 1.2 ){
              retcode <- 7
            } else {
              retcode <- 0
            }
          }
        }
      }
    }
    cat(sprintf("Norm of dx %10.5g\n", dxnorm))
    return(list(fhat=fhat,xhat=xhat,fcount=fcount,retcode=retcode))
  }


csminwelNew <-
  function(fcn,x0,H0,...,grad=NULL,crit=1e-7,nit,Verbose=TRUE,Long=FALSE,nCores = 1) {
    ### fcn:   the objective function to be minimized.  If it has a "gradient" attribute (like output of deriv), that is used
    ###        as analytic gradient.  If it has a "hessian" attribute, that is used as the hessian.
    ### fcn0:  Lean version of fcn, without grad or hessian attributes.  May save time to provide this. (Removed this option for now.)
    ### x0:    initial value of the parameter vector
    ### H0:    initial value for the inverse Hessian.  Must be positive definite, if used.  (Not used if attr(fcn,"hessian") exists.)
    ### grad:  If this is a numerical vector and attr(fcn,"gradient") is not present, then grad is used as an initial gradient vector.
    ###        Useful for restarting if numerical gradient calculation is slow.
    ### crit:  Convergence criterion.  Iteration will cease when it proves impossible to improve the
    ###        function value by more than crit.
    ### nit:   Maximum number of iterations.
    ### ...:   A list of optional length of additional parameters that get handed off to fcn each
    ###        time it is called.
    ###        Note that if the program ends abnormally, it is possible to retrieve the current x,
    ###        f, and H from the files g1.mat and H.mat that are written at each iteration and at each
    ###        hessian update, respectively.  (When the routine hits certain kinds of difficulty, it
    ###        writes g2.mat and g3.mat as well.  If all were written at about the same time, any of them
    ###        may be a decent starting point.  One can also start from the one with best function value.)

    ### nCores: specifies number of cores to use in parallel computation of gradient (can be useful
    ### if there are lots of variables). Added by Karthik Sastry, July 2014.

    dots <- list(...) # (need this to save these arguments in case of cliffs)
    nx <- length(x0)
    done <- 0
    itct <- 0
    fcount <- 0
    snit <- 100
    badg1 <- badg2 <- badg3 <- FALSE
    f0 <- fcn(x0,...)
    NumGrad <- is.null(attr(f0,"gradient"))
    NumHess <- is.null(attr(f0,"hessian"))
    badg <- FALSE
    ## browser()
    if( f0 > 1e50 ) {
      stop(paste("Bad initial parameter. f0 =", f0))
    }
    if( NumGrad ) {
      if( !is.numeric(grad) ) {
        gbadg <- numgrad(fcn,x0,nCores = nCores,...)
        g <- gbadg$g
        badg <- gbadg$badg
      } else {
        badg <- FALSE
        ## This is dangerous if you use a saved g file and it
        ## turns out to have been "bad".  We used to set badg to TRUE if
        ## grad contained any zeros.
        g <- grad
      }
    } else {
      g <- attr(f0,"gradient")
      badg <- FALSE
      gbadg <- list(g=g, badg=badg)
    }
    retcode3 <- 101
    x <- x0
    f <- f0
    if (is.null(attr(f0,"hessian"))) {
      H <- H0
    }else{
      H <- attr(f0,"hessian")
    }
    cliff <- 0
    while( !done ) {
      g1 <- NULL; g2 <- NULL; g3 <- NULL

      if (Verbose) {
        ##addition fj. 7/6/94 for control
        cat('-----------------\n')
        cat('-----------------\n')
        cat(sprintf('f at the beginning of new iteration, %20.10f',f),"\n")
      }
      if (!Long && Verbose) { # set Long=TRUE if parameter vector printouts too long
        cat("x =\n")
        print(x,digits=12)
      }
      ##-------------------------X
      itct <- itct+1

      #	if (itct == 9) browser()

      itout <- csminit(fcn,x0=x, f0=f, g0=g,badg=badg, H0=H,...)
      f1 <- itout$fhat
      x1 <- itout$xhat
      fc <- itout$fcount
      retcode1 <- itout$retcode
      fcount <- fcount+fc
      if( retcode1 != 1 ) {         # Not gradient convergence
        ## if( retcode1==2 || retcode1==4) {
        ##   wall1 <- TRUE; badg1 <- TRUE
        ## } else {                          # not a back-forth wall, so check gradient
        if( NumGrad ) {
          gbadg <- numgrad(fcn, x1, nCores = nCores,...)
        } else {
          gbadg <- list(g=attr(f1,"gradient"),badg=FALSE)
        }
        g1 <- gbadg$g
        badg1 <- gbadg$badg
        wall1 <- (badg1 || retcode1==2 || retcode1 == 4) # A wall is back-forth line search close or bad gradient
        ## g1
        save(file="g1", g1, x1, f1, dots)
        ## }
        if( wall1 && dim(H)[1] > 1) {
          ## Bad gradient or back and forth on step length.  Possibly at
          ## cliff edge.  Try perturbing search direction, if problem is not unidimensional
          ##
          Hcliff <- H+diag(diag(H) * rnorm(nx))
          cat("Cliff.  Perturbing search direction. \n")
          itout <- csminit(fcn,x0=x,f0=f,g0=g, badg=badg, H0=Hcliff,...)
          f2 <- itout$fhat
          x2 <- itout$xhat
          fc <- itout$fcount
          retcode2 <- itout$retcode
          fcount <- fcount+fc   # put by Jinill
          ## if(  f2 < f ) {
          ##   if( retcode2 == 2 || retcode2==4 ){
          ##     wall2 <- 1
          ##     badg2 <- 1
          ##   } else {
          if( NumGrad ){
            gbadg <- numgrad(fcn, x2, nCores = nCores,...)
          }else{
            gbadg <- list(g=attr(f2,"gradient"),badg=FALSE)
          }
          g2 <- gbadg$g
          badg2 <- gbadg$badg
          wall2 <- (badg2 || retcode2==2 || retcode2==4)
          ## g2
          ## badg2
          save(file="g2", g2, x2, f2, dots)
          if( wall2 ){
            cat('Cliff again.  Try traversing\n')
            normdx <- sqrt(sum(x2-x1)^2)
            if( normdx < 1e-13 ) { # two trial x's too close, can't traverse
              f3 <- f; x3 <- x; badg3 <- NumGrad;retcode3 <- 101
            }else{
              ## as.numeric below for robustness against f's being 1x1 arrays
              gcliff <- (x2 - x1) * (as.numeric(f2 - f1)/(normdx^2))
              dim(gcliff) <- c(nx,1)
              itout <- csminit(fcn,x0=x,f0=f, g0=gcliff, badg=FALSE,H0=diag(nx),...)
              f3 <- itout$fhat
              x3 <- itout$xhat
              fc <- itout$fc
              retcode3 <- itout$retcode
              fcount <- fcount+fc  # put by Jinill
              ## if( retcode3==2 || retcode3==4 ) {
              ##   wall3 <- 1
              ##   badg3 <- 1
              ## } else {
              if( NumGrad ) {
                gbadg <- numgrad(fcn, x3, nCores = nCores,...)
              }else{
                gbadg <- list(g=attr(f3,"gradient"),badg=FALSE)
              }
              g3 <- gbadg$g
              badg3 <- gbadg$badg
              wall3 <- (badg3 || retcode3==2 || retcode3==4)
              ## g3
              ## badg3
              save(file="g3", g3, x3, f3, dots)
            }
          } else { # wall1, but not wall2, so pack f3, etc with initial values
            f3 <- f
            x3 <- x
            badg3 <- NumGrad              #i.e., use the gradient if it's analytic
            g3 <- g
            retcode3 <- 101
          }
        } else {     # no walls, or one-dimensional, so no perturbations
          f3 <- f
          x3 <- x
          badg3 <- NumGrad
          g3 <- g
          retcode3 <- 101
          f2 <- f
          x2 <- x
          badg2 <- NumGrad
          g2 <- g
          retcode2 <- 101
        }
      } else { # gradient convergence --- csminit just looked at gradient and quit
        f1 <-  f2 <-  f3 <- f; g1 <- g <- FALSE; badg2 <-  badg3 <- TRUE; retcode1 <- 1; retcode2 <- 101; retcode3 <- 101
      }
      ## how to pick gh and xh:
      ## pick first fj that improves by crit and has badg==FALSE, starting with j=3
      ## may pick other than min fj to stay away from cliff
      if( f3 < f-crit & badg3==0 ) {
        ih <- 3
        fh <- f3;xh <- x3;gh <- g3;badgh <- badg3;retcodeh <- retcode3
      } else {
        if( f2 < f-crit && badg2==0 ) {
          ih <- 2
          fh <- f2;xh <- x2;gh <- g2;badgh <- badg2;retcodeh <- retcode2
        } else {
          if( f1 < f-crit && badg1==0 ) {
            ih <- 1
            fh <- f1;xh <- x1;gh <- g1;badgh <- badg1;retcodeh <- retcode1
          } else {
            ## none qualify, so pick the one that improves most.
            fh <- min(f1,f2,f3)
            ih <- which.min(c(f1,f2,f3))
            cat("ih =",ih,"\n")
            xh <- switch(ih,x1,x2,x3)
            retcodeh <- switch(ih,retcode1,retcode2,retcode3)
            nogh <- (!exists("gh",inherits=FALSE) || is.null(gh))
            if( nogh ) {
              if( NumGrad ) {
                gbadg <- numgrad(fcn,xh, nCores = nCores,...)
              } else {
                gbadg <- list(g=switch(ih,attr(f1,"gradient"),attr(f2,"gradient"),attr(f3,"gradient")),badg=FALSE)
              }
              gh <- gbadg$g
              badgh <- gbadg$badg
            }
            badgh <- NumGrad
          }
        }
      }
      ## end of picking
      ##ih
      ##fh
      ##xh
      ##gh
      ##badgh
      stuck <- (abs(fh-f) < crit)
      if ( (!badg) && (!badgh) && (!stuck) ) {
        if(NumHess){
          H <- bfgsi(H,gh-g,xh-x)
        } else {
          H <- attr(fh,"hessian")
        }
      } else {
        cat("skipped bfgsi\n")
      }
      if( Verbose ) {
        cat("----\n")
        cat(sprintf('Improvement on iteration %8.0f = %18.9f\n',itct,f-fh))
      }
      if( itct > nit ) {
        cat("iteration count termination\n")
        done <- 1
      } else {
        if( stuck ) {
          cat("improvement < crit termination\n")
          done <- 1
        }
        rc <- retcodeh
        switch( rc ,
                cat("zero gradient\n"),    #1
                cat("back and forth on step length never finished\n"), #2
                cat("smallest step still improving too slow\n"),       #3
                cat("back and forth on step length never finished\n"), #4
                cat("largest step still improving too fast\n"),        #5
                cat("smallest step still improving too slow, reversed gradient\n"), #6
                cat("warning: possible inaccuracy in H matrix\n"), #7
        )
      }
      f <- fh
      x <- xh
      g <- gh
      badg <- badgh
    }                                     # while !done
    return(list(fh=fh,xh=xh,gh=gh,H=H,itct=itct,fcount=fcount,retcodeh=retcodeh,...))
  }


numgrad <-
  function(fcn, x, nCores = 1,...) {
    ## fcn can return a vector, in which case numgrad returns a matrix.
    # delta <- 1e-5
    ## delta <- 1e-8
    nVar <- length(x)
    ## we tolerate x's that may be n x 1, 1 x n, or R vectors (with no dim),
    ## but note that g comes out as n x k matrix regardless.
    #tvec <- delta*diag(n)
    f0 <- fcn(x,...)
    # g <- matrix(0,n,k)
    # badg <- FALSE


    # for (i in 1:n){
    # scale <- 1
    # tvecv <- tvec[,i]
    # if(is.null(dim(x))){
    # tvecv <- as.vector(tvecv)
    # }else{
    # dim(tvecv) <- dim(x)
    # }
    # g0 <- (fcn(x+scale*tvecv,...) - f0)/(scale*delta)
    # if (max(abs(g0))< 1e50){
    # g[i, ] <- as.vector(g0)
    # }else{
    # cat("bad gradient ------------------------\n")
    # badg <- TRUE
    # }
    # }

    # Parallel implementation:
    vars <- 1:nVar
    listOutput <- mclapply(vars, pderv, fcn, x, f0, ..., mc.cores = nCores)

    g <- t(matrix(unlist(listOutput), nrow = length(f0)))
    badg <- any(is.nan(g))


    return(list(g=g,badg=badg))
  }


numHess <-
  function(fcn, x, ...) {
    f1 <- fcn
    n <- length(x)
    h <- matrix(0, n, n)
    f2 <- function(z, ...) { numgrad(fcn=f1, z, ...)$g}
    h <- numgrad(fcn=f2, x=x, ...)$g
    return(h)
  }



pderv <-
  function(iVar, fcn, x, f0, ...){

    # Calculates partial derivative
    # Pulled out of numgrad to allow an easier multi-core implementation

    delta <- 1e-5 #can be adjusted
    x[iVar] <- x[iVar] + delta

    g0 <- (fcn(x, ...) - f0) / (delta)

    #checking for very large gradient value
    if (max(abs(g0) < 1e50)) {
      derv <- as.vector(g0)
    } else {
      c(derv) <- NaN #all values set to NaN
      #	dervbad <- TRUE
      cat("bad gradient ------------------------\n")
    }

    return(derv)
  }

