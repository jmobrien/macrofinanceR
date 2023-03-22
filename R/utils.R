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
