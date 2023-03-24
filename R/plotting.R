#' Impluse plots
#'
#' @param ir
#' @param irdraws
#' @param conf
#' @param blocks
#' @param filename
#' @param format
#' @param addtitle
#' @param nsteps
#' @param shocknames
#' @param varnames
#' @param ptype
#' @param color
#' @param alpha
#' @param gr
#' @param width
#' @param height
#' @param response_ylim
#' @param xtick
#'
#' @return
#' @export
#'
#' @examples
impulseplots <-
  function(ir, ## array of impulse response objects, nvar x nshock x nperiod
           irdraws = NULL, ## nvar x nshock x nperiod x ndraw
           conf = c(.68,.90), ## what confidence bands to put on the plot
           blocks = NULL, ## possible block structure, see below for example
           filename = 'impulse_plot',
           format = 'pdf',
           addtitle = FALSE, ## put a big title on the top?
           nsteps = 60, ## number of steps to use
           shocknames = NULL, ## names of the shocks to print as titles
           varnames = rep('',dim(ir)[1]), ## names of the response variables
           ptype = 'Impulse response', ## part of the title in some cases
           color = c(0, 0, 1), ## base color (default blue)
           alpha = c(0.5,0.3), ## how much alpha for each oonf band color
           gr = 0.7, ## what shade of grey
           width = 5, ## inches
           height = 8,
           response_ylim = NULL, ## can specify the size of each response row
           xtick = NULL ## where to put xticks; default, every 12 months
  ) {

    ## Plots impulse response draws
    ## more current than other code in the same folder
    ## Handles blocking of the impulse responses into several minigarphs




    ## Some formatting issues
    ## filename = paste('/home/karthik/R/svar_kit/plotting/plots/', filename,sep='')


    if (format == 'pdf'){ ## different methods for the "light color" for pdf and eps
      shades <- grDevices::rgb(color[1], color[2], color[3], alpha)
    } else if (format == 'eps'){
      lattice::trellis.device(device="postscript", color = TRUE)
      #setEPS()
      ## grDevices::postscript(filename, width = width, height = height)

      # Cannot do transparent colors, so here is a crude workaround
      alphafy <-
        function(col,alpha=1) {
          rr <-
            1-alpha*(1-c(col/255))
          return(grDevices::rgb(rr[1],rr[2],rr[3]))
        }
      color <- alphafy(color, alpha)

    } #else raise some kind of error?

    if (is.null(varnames)){ ## get the names of each shock
      varnames <- rownames(ir)
    }

    if (is.null(xtick)){
      xtick <- seq(1,nsteps,by=12)
    }

    ##
    ## Calculating the IR quantiles and
    ## y limits for plots
    ##

    nv <- dim(ir)[1] ## variables
    ns <- dim(ir)[2] ## shocks

    ir  <- ir[,,1:nsteps] ## trimming unnecessary steps

    if (!is.null(irdraws)){
      ## get the correct quantiles
      irq <- apply(irdraws[,,1:nsteps,],1:3,quantile,
                   ##probs <- c(rev((1 - conf) / 2), .5 + conf/2))
                   probs <- c((1-conf)/2,0.5 + conf/2))
      irq <- aperm(irq,perm=c(2,3,4,1))

      nconf <- length(conf)
    } else{
      irq <- array(0,c(dim(ir),1))
      nconf <- NULL
    }

    ## determine the ylimits for each variable, if necessary
    if (is.null(response_ylim)){
      response_ylim <- matrix(0,nv,2) ## idea is that this should be identical for all blocks
      for (iv in 1:nv){
        response_ylim[iv,1] <- min(ir[iv,,],irq[iv,,,],0)
        response_ylim[iv,2] <- max(ir[iv,,],irq[iv,,,],0)
      }
    }

    ##
    ## Blocking prep
    ##

    if (is.null(blocks)){
      blocks <- list()
      blocks[[1]] <- list(x = 1:nv, y = 1:nv)
    }

    nb <- length(blocks)

    for (ib in 1:nb){

      ## open the file
      if (format == 'pdf'){
        fname <- paste(filename,'_',ib,'.pdf', sep = '')
        grDevices::pdf(fname, width = width, height = height)
      } else if (format == 'eps'){
        fname <- paste(filename,'_',ib,'.pdf', sep = '')
        lattice::trellis.device(device="postscript", color = TRUE)
        grDevices::postscript(fname,width=width,height=height)
      }


      par(mfrow = c(length(blocks[[ib]]$x),length(blocks[[ib]]$y)),
          col.lab="black",col.main="black",
          oma=c(1,5,1,2), mar=c(.5,.25,.5,.25), tcl=-0.1, mgp=c(3,1,0))


      for (rv in blocks[[ib]]$x){ ## responses

        ptitle <- varnames[rv] ## name of data series

        for (sv in blocks[[ib]]$y){ ## shocks
          plot(ir[rv,sv,],
               ylim = response_ylim[rv,],
               type = 'l',
               lwd = 0.5,
               xlab = '',
               ylab = '',
               yaxt  = 'n',
               xaxt = 'n',
               ## fg = gray(gr),
               xlim = c(1,nsteps), xaxs = 'i')

          ytick <- pretty(response_ylim[rv,],4)

          abline(h = ytick, lty = 'dotted')## ,col=gray(gr))
          abline(v = xtick, lty = 'dotted') ## ,col=gray(gr))

          abline(a=0,b=0,lwd=0.75)

          if (!is.null(nconf)){ ## plot confidence bands
            for (ic in 1:nconf){
              polygon(c(1:nsteps, nsteps:1),
                      c(irq[rv,sv,,ic],rev(irq[rv,sv,,ic+nconf])),
                      col = shades[ic],
                      border = NA)
            }
          }

          ##Adding variable name and axis on leftmost plot
          if (sv == blocks[[ib]]$y[1]){
            mtext(ptitle, side = 2, line = 5, cex = 0.5, las = 1, adj = 0)
            axis(side = 2, cex.axis = .75, las = 1,at=ytick)
          }

          ## ##Right side axis labels on right most plot
          ## if (sv == blocks[[ib]]$y[length(blocks[[ib]]$y)]){
          ##     axis(side = 4, cex.axis = .75, las = 1)
          ## }

          ##Shock name if appropriate
          if (!is.null(shocknames) &&
              rv == blocks[[ib]]$x[1]) {
            mtext(shocknames[sv], side = 3, line = 0, cex = .5)
          }
        }
      }

      if (addtitle){
        bigtitle <- paste(type, 'over', as.character(nSteps), 'periods', sep = ' ')
        title(bigtitle, outer = TRUE, cex = 1.2)
      }

      dev.off()

    }
  }


#' Title forecast plotting
#'
#' @param fcout forecasting results
#' @param ydata Original variable dataframe
#' @param dateseq Vector of points from which to draw forecasting lines
#' @param vnames Names of variables to use for each plot
#' @param fulldates full range of dates associated with time series
#' @param frequency how often time points are sampled per year (e.g. 12 for
#'   monthly data, 4 for quarterly data)
#' @param ymat
#' @param filename
#' @param cushion
#'
#' @return
#' @export
plotfc <- function(fcout, ydata, dateseq,
                   vnames, fulldates,
                   frequency = 12,
                   ymat = c(4.35,4.70,
                            4.55,4.65,
                            8.2,8.5,
                            7.00,7.40,
                            7.20,7.60,
                            -0.02,0.06,
                            5.6,6.4,
                            -0.01,0.05,
                            0.00,0.12,
                            -0.01,0.05),
                   start_date = NULL,
                   filename = 'fcout',
                   cushion = 0){

  ## Come up with a nice way of plotting all the forecasts
  ## assume ydata starts at datseq[1]

  nfwd <- dim(fcout)[1] ## how far fwd does each fc go
  nx <- dim(fcout)[2]
  nd <- dim(fcout)[3] ## number of dates

  ## fulldates <- seq(from = datseq[1], to = datseq[nd], by = 'month') ## dates for yvec

  ymat <- matrix(ymat, nrow = 2)
  fclist <- list()

  ## make everything a time series
  for (idate in 1:nd){
    ## sdate <- seq(from = datseq[idate], length = 2, by = 'month')
    ## sdate <- sdate[2] ## month after

    sdate <- dateseq[idate]
    iy <- which(fulldates == dateseq[idate])
    yt <- ydata[iy,] # y data

    ## JMO this is a terrible hack but I'm putting it in to get this to work right now with quarters
    if(frequency == 12){
      startpoint <-
        c(as.numeric(format(sdate,format = '%Y')),
          as.numeric(format(sdate,format = '%m')))
    } else if (frequency == 4) {
      startpoint <-
        sdate |>
        strsplit(" ?[qQ] ?") |>
        unlist() |>
        as.numeric()
    }


    fclist[[idate]] <-
      unlist(yt) |>
      rbind(fcout[,,idate]) |>
      ts(start = startpoint,
         frequency = frequency)

    ## ylist[[idate]] = ts(yt,
    ##                     start =
    ##                         c(as.numeric(format(sdate,format = '%Y')),
    ##                           as.numeric(format(sdate,format = '%m'))),
    ##                     frequency = 12)
  }

  ## ystart = which(fulldates == dateseq[1])
  ## yend = which(fulldates == dateseq[nd]) + 1 + nfwd

  ## ytrim = ts(ydata[ystart:yend,],
  ##            start =
  ##                c(as.numeric(format(dateseq[1],format = '%Y')),
  ##                  as.numeric(format(dateseq[1],format = '%m'))),
  ##            frequency = 12)

  ## JMO this is a terrible hack but I'm putting it in to get this to work right now with quarters
  if(frequency == 12){
    trim_startpoint <-
      c(as.numeric(format(fulldates[1],format = '%Y')),
        as.numeric(format(fulldates[1],format = '%m')))
  } else if (frequency == 4) {
    trim_startpoint <-
      (fulldates[1]) |>
      strsplit(" ?[qQ] ?") |>
      unlist() |>
      as.numeric()
  }

  ytrim <-
    ydata |>
    ts(
      start = trim_startpoint,
      frequency = frequency
    )






  grDevices::pdf(paste(filename,'.pdf',sep = ''),width = 6.5, height = 8)

  ## plotting, could be adjusted for nv
  par(mfrow = c(5,2),col.lab="black",col.main="black",
      oma=c(4,4,4,4), mar=c(2,2,2,2), tcl=-0.1, mgp=c(0,0,0))

  for (iv in 1:nx){
    plot(
      ## ts(tsdata[fullrange,iv],start = c(2006,9),frequency = 12),
      ytrim[,iv],
      type = 'l', lwd = 1.5,
      ylab = '',xlab = '' #,
      # ylim = ymat[,iv]
    )

    title(main = vnames[iv])

    grid(nx = NULL, ny = nx, col = "lightgray", lty = "dotted",
         lwd = 1, equilogs = TRUE)

    for (idate in 1:nd){
      lines(fclist[[idate]][,iv], lwd = .75, col = 'red')
    }
  }


  dev.off()
}


#' generate rmse plots
#'
#' @param ... rmse's to include, from getrmse()
#' @param sdate start date for time series, as 2 element numeric: c([time unit], [time subunit in specified frequency])
#' @param frequency how often time points are sampled per year (e.g. 12 for
#'   monthly data, 4 for quarterly data)
#' @param filename file name to output pdf to
#' @param vname names of variables to be plotted
#' @param ih which horizon to plot?  Indicates position in getrmse()'s array output, not horizon length itself)
#' @param height plot height
#' @param width plot width
#' @param diff if TRUE, calculate as differences
#' @param rshade if TRUE, shade recession periods
#' @param rlist numeric position (in years) of start/end points of recession periods
#' @param colors what colors to use for each of the rmse's added?  Will generate from rainbow() if unspecified
#'
#'
#' @return
#' @export
plotrmse <- function(...,
                     sdate,
                     frequency = 12,
                     filename,
                     vname = rep(NULL,nv),
                     ih = 1,
                     height = 8, width = 6.5,
                     diff = FALSE,
                     rshade = TRUE,
                     rlist = c(1980, 1980 + 6/12,
                               1981 + 6/12,  1982 + 10/12,
                               1990 + 6/12, 1991 + 2/12,
                               2001 + 2/12, 2001 + 10/12,
                               2007 + 11/12, 2009 + 5/12),
                     colors = NULL
                     ){

  # Pull out the rmse's
  rmses <- list(...)

  names(rmses) <- paste0("rmse", seq_along(rmses))

  if(missing(colors)){
    line_colors <- rainbow(length(rmses))
  } else {
    line_colors <- colors
  }

  # list2env(rmses)

  nv <- dim(rmses[[1]])[2]

  grDevices::pdf(paste(filename,'.pdf',sep = ''),height = 8, width = 6.5)

  par(mfrow = c(nv,1),col.lab="black",col.main="black",
      oma=c(1,3,1,2), mar=c(2,1,1,1), tcl=-0.1, mgp=c(0,0,0))

  nrmse <- c(1,6,12,24,48) * (frequency/12) ## number of months (JMO - adjusted)
  nrmse <- nrmse[ih]
  ## snip these off the end
  Tobs <- dim(rmses[[1]])[3]

  irange <- 1:(Tobs-nrmse)

  ## rmse1 <- rmse1[,,1:(Tobs-nrmse)]
  ## rmse2 <- rmse1[,,1:(Tobs-nrmse)]
  ## rmse3 <- rmse1[,,1:(Tobs-nrmse)]

  for (iv in 1:nv){
    ## plot each variable

    if (diff){
      ## difference mode

      drmse <- rmses[[1]][ih,iv,] - rmse[[2]][ih,iv,]

      plot(ts(drmse,sdate,frequency = frequency), type = 'l', lwd = 1,
           ylab = NULL, xlab = NULL)
      abline(a = 0, b = 0, lwd = 0.75) ## axis
      grid(lwd = .5)

    } else {
      ##ymin <- min(rmse1[ih,iv,],rmse2[ih,iv,])
      ymin <- 0
      # ymax <- max(rmses[[1]][ih,iv,],rmse2[ih,iv,]) * 1.1 ## breathing room
      ymax <- (sapply(rmses, \(x){x[ih,iv,]}) |> max()) * 1.1 ## breathing room

      ## plot(ts(rmse1[ih,iv,],sdate,frequency=12), type = 'l', lwd = .75, col = 'blue',
      ##      ylab = NULL, xlab = NULL, bty = 'l',ylim = c(ymin,ymax))

      plot(ts(rmses[[1]][ih,iv,],sdate,frequency=frequency), type = 'n', lwd = .75, col = line_colors[1],
           ylab = NULL, xlab = NULL, bty = 'l',ylim = c(ymin,ymax), panel.first = {
             grid(lwd = 1)
           })

      if (rshade){
          rlist <- matrix(rlist,2)
        for (ir in 1:5){

          polygon(c(rlist[,ir],rev(rlist[,ir])),
                  c(-1e10,-1e10,1e10,1e10),
                  col = 'grey',border = NA)
        }
      }

      for(i in seq_along(rmses)){

        linedat <-
          ts(rmses[[i]][ih,iv,irange],sdate,frequency=frequency)

      lines(linedat, lwd = .75, col = line_colors[i])

      }
      ## abline(a = 0, b = 0, lwd = 0.75) ## axis


    }


    title(main = vname[iv])


  }

  dev.off()
}



