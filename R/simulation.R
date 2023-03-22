sim_model <-
  function(tvvout, ## output of optimization
           lambda = matrix(1,nv,6), ## - log variances
           lmdp = rep(10,6), ## number of periods in each regime
           burn = 0, ## number of periods in the first regime to run as "burn in"
           y0 = c(rep(0,nv*nlag),1), ## initial state
           delta = NULL ## variance modifiers
  ) {

    ## Generates fake data from our SVAR model

    ## NOTE: if delta is a scalar, the model assumes you wanted t distribution with
    ## 2*delta df, and simulates delta_it as appropriate inverse gamma



    ## Get coefficients ----

    A <- tvvout$A ## A0 matrix
    Ai <- solve(A)

    vout <- tvvout$vout

    Ap <- vout$var$By ## structural form coefficients A+
    nv <- dim(Ap)[1]
    nlag <- dim(Ap)[3]
    C <- matrix(vout$var$Bx,nv,1) ## constant

    ##
    ## Generate a system matrix
    ##

    By <- Ap ## reduced form coefs
    for (ilag in 1:nlag){
      By[,,ilag] <- Ai %*% Ap[,,ilag]
    }

    rsys <- sysmat(By,Ai %*% C) ## system for the reduced form
    ##
    ## Define a (sub) function that moves the system n steps forward
    ##

    moven <- function(x0,sysmat,n,nv,nlag,lambda,Ai){
      ## noise <- exp(-lambda) * matrix(rnorm(n * nv),nv,n)
      noise <- Ai %*% (sqrt(lambda) * matrix(rnorm(n * nv),nv,n))
      xout <- matrix(0,nv,n) ## put draws in here
      for (ii in 1:n){
        x0 <- sysmat %*% x0 + c(noise[,ii],rep(0,1 + (nlag - 1)*nv))
        xout[,ii] <- x0[1:nv]
      }
      return(list(xout = xout, xf = x0, shocks = noise))
    }

    ## translate the data into a vector, if necessary
    if (length(dim(y0)) > 1){
      ## assume the format is like
      ##       v1  v2 ....
      ## t <- 1
      ## t <- 2
      ## ...

      y0 <- c(t(y0[nlag:1,]),1)
    }

    ##
    ## Burn in, if desired
    ##

    if (burn > 0){
      lambda_burn <- matrix(lambda[,1],dim(lamda)[1],burn) ## repeat period 1 variances
      burn_in <- moven(y0,rsys,burn,nv,nlag,lambda_burn,Ai)
      y0 <- burn_in$xf
    } else {
      burn_in <- NULL
    }
    ##
    ## Run the main simulation
    ##

    replambda <- lambda[,rep(1:dim(lambda)[2],times = lmdp)]
    ## edit: drop the first 10 obs
    replambda <- replambda

    if (!is.null(delta)){ ## t case
      if (length(delta) > 1) { ## multiply by t adjusters
        replambda <- replambda * t(exp(delta)) ## log variancs
      } else { ## unconditional simulation
        gd <- 1 / rgamma(n = length(replambda), shape = delta, rate = delta)
        replambda <- replambda * gd
      }
    }

    sim_out <- moven(y0, rsys, dim(replambda)[2],
                     nv, nlag,
                     replambda,Ai)

    ##
    ## get reduced form
    ##
    ## rfsim <- t(Ai %*% sim_out$xout)

    return(list(sim_out = sim_out,
                sshocks = sim_out$shocks,
                burn_in = burn_in))
  }
