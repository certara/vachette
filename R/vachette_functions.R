
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import prospectr
#' @importFrom Hmisc approxExtrap
#' @importFrom purrr map_dfr
#' @importFrom magrittr %>%
#' @importFrom minpack.lm nls.lm nls.lm.control
#' @importFrom utils askYesNo capture.output head
#' @rawNamespace import(stats, except = c(filter, lag))
#'
NULL

#' Multi approx
#'
#' Find monotonic increasing and decreasing segments
#' Used for both derivative and asymptote identification
#'
#' @param x x values
#' @param y y values
#' @param yout y value to approximate
#' @param tol tolerance for max and min identification
#'
#' @noRd
#'
multi.approx <- function(x,y,yout,tol=1e-9) {

  seg <- NULL
  yrange <- max(y,na.rm=T)-min(y,na.rm=T)
  tolval <- tol*yrange

  # Checks
  if(length(x) != length(y)) return(NA)
  if(length(x) <= 1)         return(NA)

  # Based on previous two and next two data points (two to prevent impact noise)
  # "True" max's or min's if both 2 preceding and 2 proceeding data points
  # point to that max or min. With sufficient difference between all 5 points (tolval)

  # initialize
  df      <- data.frame(x=x,y=y)
  df$max    <- 0
  df$min    <- 0

  # First data point always part of segment
  for(i in c(3:(dim(df)[1]-2)))
  {
    if(df$y[i] > df$y[i+1] & df$y[i+1] > df$y[i+2] &
       df$y[i] > df$y[i-1] & df$y[i-1] > df$y[i-2] &
       (abs(df$y[i-2] - df$y[i-1]) > tolval |
        abs(df$y[i-1] - df$y[i])   > tolval |
        abs(df$y[i]   - df$y[i+1]) > tolval |
        abs(df$y[i+1] - df$y[i+2]) > tolval)) df$max[i] <- 1

    if(df$y[i] < df$y[i+1] & df$y[i+1] < df$y[i+2] &
       df$y[i] < df$y[i-1] & df$y[i-1] < df$y[i-2] &
       (abs(df$y[i-2] - df$y[i-1]) > tolval |
        abs(df$y[i-1] - df$y[i])   > tolval |
        abs(df$y[i]   - df$y[i+1]) > tolval |
        abs(df$y[i+1] - df$y[i+2]) > tolval)) df$min[i] <- 1
  }

  # All potential segments
  iseg      <- 1
  df$seg[1] <- iseg
  df.add    <- NULL
  for(i in c(2:dim(df)[1]))
  {
    # If a max or min was identified
    if(df$min[i] == 1 | df$max[i] == 1)
    {
      # Save last rec to add (avoid discontinuities)
      df.add <- rbind(df.add,df[i,] %>% mutate(seg=iseg))
      # Start new segment
      iseg <- iseg+1
    }
    df$seg[i] <- iseg
  }
  df <- rbind(df,df.add) %>%
    arrange(x,seg)

  # Check for each segment if it crosses yout, if not, set to zero.
  # Test if all points above or below yout, accounting for tolval
  for(iseg in c(1:length(unique(df$seg))))
  {
    # Either all points +tolval are above yout, or all points -tolval are below yout
    if((sum(df$y[df$seg==iseg] > (yout+tolval))   == 0) |
       (sum(df$y[df$seg==iseg] < (yout-tolval))   == 0))   df$seg[df$seg==iseg] <- 0
  }

  # Now find approx's (only if 2 or more segments points which should always be the case)
  segs <- unique(df$seg[df$seg!=0])
  xcollect <- NULL
  if(length(segs)>=1)
  {
    for(iseg in segs)
    {
      xhit     <- approx(df$y[df$seg==iseg],df$x[df$seg==iseg],xout=yout)$y
      xcollect <- c(xcollect,xhit)
    }
  }

  # Remove NAs
  xcollect <- xcollect[!is.na(xcollect)]

  # Return NA if no hits at all
  if(length(xcollect)==0) xcollect <- c(NA)

  return(xcollect)
}


# Get min/max or inflec for curve via polynomial fit
extremes <- function(x,y,type='minmax',polyorder=6, tol.poly=0.001)
{
  warning_out <- get("warning_out", envir = vachette_env)
  mydata <- data.frame(x=x,y=y)

  # Polynomial fit
  fit=lm(y~poly(x,polyorder,raw=T),data=mydata)
  f0.coeff <- as.numeric(fit$coefficients)

  # Set NA terms to 0
  f0.coeff[is.na(f0.coeff)] <- 0

  mydata$fit <- predict(fit)

  # y'(x)  (f1)
  f1.coeff <- NULL
  for(icoeff in c(2:length(f0.coeff))) # Skip intercept --> 0
  {
    add      <- f0.coeff[icoeff] * (icoeff-1)
    f1.coeff <- c(f1.coeff,add)
  }

  # inflection points. solve y''(x) = 0 (f2 = 0)
  f2.coeff <- NULL
  for(icoeff in c(2:length(f1.coeff))) # Skip intercept --> 0
  {
    add      <- f1.coeff[icoeff] * (icoeff-1)
    f2.coeff <- c(f2.coeff,add)
  }

  # ---

  if (type=='minmax') solutions = polyroot(f1.coeff)
  if (type=='inflec') solutions = polyroot(f2.coeff)

  # AL 04182024
  # solve the poly
  fx.zero <- Re(solutions)[!(round(Im(solutions), 5))]
  if (type=='minmax') {
    poly.eval=array(0,length(fx.zero))
    for (j in 1:length(fx.zero)) {
      y.eval = f1.coeff[1]
      for (i in 1:(length(f1.coeff)-1)) {
        y.eval = y.eval + f1.coeff[i + 1]*fx.zero[j]^i
      }
      poly.eval[j] = y.eval
    }
  }
  if (type=='inflec') {
    poly.eval=array(0,length(fx.zero))
    for (j in 1:length(fx.zero)) {
      y.eval = f2.coeff[1]
      for (i in 1:(length(f2.coeff)-1)) {
        y.eval = y.eval + f2.coeff[i + 1]*fx.zero[j]^i
      }
      poly.eval[j] = y.eval
    }
  }

  # If there is a solution
  if (!is.na(fx.zero[1]))
  {
    out   <- fx.zero <  min(mydata$x) | fx.zero >  max(mydata$x)
    keep  <- fx.zero >= min(mydata$x) & fx.zero <= max(mydata$x)

    # Check domain
    # if(sum(out)>0)
    # {
    #   wmsg <- paste0("WARNING: ignoring inflection point solutions using polynomial fit at:\n",
    #                  "x = ", paste(fx.zero[out]), "\n", collapse= " ")
    #   log_output(wmsg)
    #   warning_out <- c(warning_out, wmsg)
    #   assign("warning_out", warning_out, envir = vachette_env)
    # }

    if(sum(keep)==0) fx.0 <- NA  # No solution
    if(sum(keep)==1) {
      fx.0 <- fx.zero[keep]
      poly.eval <- poly.eval[keep]
    }
    # if(sum(keep)>1) stop('Cannot handle: multiple inflection points detected using polynomial fit')
    if(sum(keep)>1)  {
      fx.0 <- c(fx.zero[keep])
      poly.eval <- c(poly.eval[keep])
    }


    # only for inflection point
    if(sum(keep)!=0 & type=='inflec')
    {
      # take the first solution to the poly
      fx.0 = min(fx.0)
      poly.eval = poly.eval[which(fx.0 == min(fx.0))]

      # if above tolerance is not a true zero
      if (min(abs(poly.eval))>tol.poly) {
        wmsg <- paste0("WARNING: ignoring inflection solution of poly as > tol.poly:\n",
                       "x = ", paste(fx.0[which(abs(poly.eval)==min(abs(poly.eval)))]),collapse=' ')
        log_output(wmsg)
        warning_out <- c(warning_out, wmsg)
        fx.0 <- NA
      }
    } else if (sum(keep)!=0 & type=='minmax') {
      if (length(which(abs(poly.eval)<tol.poly))!=0) {
        fx.0 <- fx.0[which(abs(poly.eval)<tol.poly)]
      } else {
        wmsg <- paste0("WARNING: ignoring minmax solution of poly as > tol.poly:\n",
                       "x = ", paste(fx.0[which(abs(poly.eval)>=tol.poly)]),collapse=' ')
        log_output(wmsg)
        warning_out <- c(warning_out, wmsg)
        fx.0 <- NA
      }
    }

    log_output(paste0("Searched for: ",type))
    log_output(paste0("Returned: ",fx.0))
    assign("warning_out", warning_out, envir = vachette_env)


      return(fx.0)


  }

  # No solution
  return(NA)

}

#' Get x multi landmarks
#'
#' Find initial landmark position using f1 and f2 derivatives and Savitzky Golay smoothing
#' Default window = 17 gridpoints
#'
#' @param x x values
#' @param y y values
#' @param w window size
#' @param tol tolerance for max and min identification
#'
#' @noRd
#'
get.x.multi.landmarks <- function(x,y,w=17,tol=1e-9) {
  # Using splinefun to determine derivatives
  # tol to be passed on to multi.approx function

  # JL230909. Two step process
  # 1. split by extremes
  # 2. between pair of extremes search for inflection points
  warning_out <- get("warning_out", envir = vachette_env)


  # Make numeric
  x <- as.numeric(x)
  y <- as.numeric(y)

  # z <- data.frame(x=x,y=y)
  # write.csv(z,'full-oral.csv')

  # Make sure no NA's
  z <- data.frame(x=x,y=y) %>% filter(!is.na(x),!is.na(y))
  x <- z$x
  y <- z$y

  # Checks
  if(length(x) != length(y)) return(NA)
  if(length(x) <= 1)         return(NA)

  # Keep first and last datapoint (before noise removal)
  xstart <- x[1]
  xend   <- x[length(x)]

  # JL 240313
  # Quick check if there should be a min or max
  # If not detected, too few grid points to detect and shortest piece is to be removed
  # E.g. short infusion
  SETMAX <- FALSE
  SETMIN <- FALSE
  if(max(y) > y[1] & max(y) > y[length(y)])
  {
    log_output("Max expected")
    SETMAX  <- TRUE
    gridmax <- c(1:length(y))[grepl(max(y),y)]
  }
  if(min(y) < y[1] & min(y) < y[length(y)])
  {
    log_output("Min expected")
    SETMIN <- TRUE
    gridmin <- c(1:length(y))[grepl(min(y),y)]
  }

  # ------------------------------------------------
  # Stop detecting extremes and inflection points if curve does not change anymore within
  #      tolerance with respect to last datapoint provided
  # y.noise <- data.frame(noise = abs(y - y[length(y)]), flag = (abs(y - y[length(y)])) < tol) %>% mutate(n=row_number())
  # Last non-noise data point:
  # y.last.data <- max((y.noise %>% filter(flag==F))$n)

  # ------- Start collecting landmarks -----------

  # First data point
  lm    <- data.frame(x=xstart,type='start')

  # Last data point - may be "asymptotic"
  add  <- data.frame(x=xend,type='end')
  lm   <- rbind(lm,add)

  # ------------ MIN/MAX -----------

  nminmax <- 0
  # (Local) extremes - not tested in case of multiple min/max's
  f1.0  <- photobiology::get_peaks(x,y,span=w)$x
  npeak <- length(f1.0)

  if(npeak == 0 & SETMAX)
  {
    add    <- data.frame(x=x[gridmax],type='max')
    lm     <- rbind(lm,add)
    nminmax <- nminmax + 1
    log_output("Raw max landmark added")
  }
  if(npeak>=1)
  {
    for(ipeak in c(1:npeak))
    {
      newpeak <- f1.0[ipeak]
      # JL240216 - refine, use (same) window of grid points
      wsub <- floor(w/2)
      x.grid  <- which(x == newpeak)
      x.start <- max(1,(x.grid-wsub))
      x.end   <- min(length(x),(x.grid+wsub))
      xsub    <- x[x.start:x.end]
      ysub    <- y[x.start:x.end]

      newpeak.refined <- extremes(xsub, ysub, type='minmax', polyorder=6)

      log_output(paste0("Max detected at x = ",newpeak," and refined to x = ",newpeak.refined,"\n"))

      if(is.na(newpeak.refined))
      {
        log_output("    No max found via polynomial fit")
        log_output("    Keeping initial max x-position")
      }
      if(!is.na(newpeak.refined)) newpeak <- newpeak.refined

      # tmp
      # newpeak <- f1.0[ipeak]

      add    <- data.frame(x=newpeak,type='max')
      lm     <- rbind(lm,add)

      nminmax <- nminmax + 1
    }
  }

  f1.0 <- photobiology::get_valleys(x,y,span=w)$x
  nvalley <- length(f1.0)

  if(nvalley == 0 & SETMIN)
  {
    add    <- data.frame(x=x[gridmin],type='min')
    lm     <- rbind(lm,add)
    nminmax <- nminmax + 1
    log_output("Raw min landmark added")
  }
  if(nvalley>=1)
  {
    for(ivalley in c(1:nvalley))
    {
      newvalley <- f1.0[ivalley]
      # JL240216 - refine, use (same) window of grid points
      wsub <- floor(w/2)
      x.grid  <- which(x == newvalley)
      x.start <- max(1,(x.grid-wsub))
      x.end   <- min(length(x),(x.grid+wsub))
      xsub    <- x[x.start:x.end]
      ysub    <- y[x.start:x.end]

      newvalley.refined <- extremes(xsub, ysub, type='minmax', polyorder=6)

      log_output(paste0("Min detected at x = ",newvalley," and refined to x = ",newvalley.refined,"\n"))

      if(is.na(newvalley.refined))
      {
        log_output("    No min found via polynomial fit")
        log_output("    Keeping initial max x-position")
      }
      if(!is.na(newvalley.refined)) newvalley <- newvalley.refined

      # tmp
      # newvalley <- f1.0[ivalley]

      add    <- data.frame(x=newvalley,type='min')
      lm     <- rbind(lm,add)

      nminmax <- nminmax + 1
    }
  }
  if (nminmax == 0) log_output(paste0('No minimum or maximum detected'))

  # ------------ INFLECTION POINTS -----------

  # Loop between extremes to detect inflection points
  # First arrange by increasing x
  lm <- lm %>% arrange(x)

  # Each interval between landmarks, incl. first and last datapoint
  for(ilm in c(2:nrow(lm)))
  {
    # Inflection point detection via polynomial fit

    # First grid point after max/min:
    xdiff1 <- x-lm$x[ilm-1]
    xpos1  <- xdiff1 >= 0
    xdiff2 <- x-lm$x[ilm]
    xneg2  <- xdiff2 <= 0

    istart <- which(x==min(x[xpos1]))
    iend   <- which(x==max(x[xneg2]))
    xsub   <- x[istart:iend]
    ysub   <- y[istart:iend]

    # Only data points before last non-noise data point
    #      and first noise point which may be the last one: y[length(y)]
    # Also remove first data point (if sharp peak, notably for i.v.)

    # first_point <- ifelse((ysub[1] < ysub[2] & ysub[3] < ysub[2]),2,1)  # y[2] is peak
    # if(ysub[1] < ysub[2] & ysub[3] < ysub[2])
    #   warning("First grid point ignored as y[1]<y[2] & y[3]<y[2]. (FOR DETECTION MAX/MIN/INFLEC)")
    # xsub <- xsub[first_point:length(xsub)]
    # ysub <- ysub[first_point:length(ysub)]

    log_output("Points kept in y-range of 2.5% to 97.5% only")

    miny  <- min(ysub)
    maxy  <- max(ysub)
    diffy <- maxy - miny

    # JL 240313: truncate only when first point is "start" or when last point is "end"
    # REMOVE_HIGH <- FALSE
    # REMOVE_LOW  <- FALSE
    # if(lm$type[ilm-1]=='start' & lm$type[ilm]=='min') REMOVE_HIGH <- TRUE # Trunc. at start
    # if(lm$type[ilm-1]=='start' & lm$type[ilm]=='max') REMOVE_LOW  <- TRUE # Trunc. at start
    # if(lm$type[ilm-1]=='max'   & lm$type[ilm]=='end') REMOVE_LOW  <- TRUE # Trunc. at end
    # if(lm$type[ilm-1]=='min'   & lm$type[ilm]=='end') REMOVE_HIGH <- TRUE # Trunc. at end

    # JL 240314: truncate only when first point is "start" or when last point is "end"
    #            PROBABLY BETTER TO TRUNC AT BOTH SIDES (THERE'S ALSO GRADUAL DECREASE FROM E.G. A MAX)
    REMOVE_HIGH <- TRUE
    REMOVE_LOW  <- TRUE
    mydata <- data.frame(x=xsub,y=ysub)
    if(REMOVE_LOW)
    {
      mydata <- mydata %>% filter(y >= (miny+0.025*diffy))
    }
    if(REMOVE_HIGH)
    {
      mydata <- mydata %>% filter(y <= (miny+0.975*diffy))
    }

    # print(paste0("Old y-range: ",miny," to ",maxy))
    # print(paste0("New y-range: ",min(mydata$y)," to ",max(mydata$y)))

    ndata     <- nrow(mydata)

    if (ndata<7)
    {
      log_output(paste0("Only ",ndata," simulated grid points for segment-",ilm, " available"))
      log_output(paste0("No attempt made to find inflection point for this segment"))
    }
    if (ndata>=7)
    {
      polyorder <- floor(ndata/4)

      # Odd order only:
      if(abs(polyorder/4 - floor(polyorder/4)) == 0) polyorder <- polyorder-1

      # Use high order, minimum=13
      if(polyorder < 13) log_output("!! Few grid points for inflection point determination !!")

      # Highest order: 70
      polyorder <- min(70,polyorder)

      log_output(paste0("Polynomial order for open end curve fitting = ",polyorder))

      # f2.0 <- inflec(mydata$x, mydata$y, polyorder=polyorder)
      f2.0 <- extremes(mydata$x, mydata$y, type='inflec', polyorder=polyorder)

      if(!is.na(f2.0[1]) & length(f2.0) > 1)
      {
        wmsg <- "WARNING: multiple inflection points detected for a segment. First inflection point only used"
        warning_out <- c(warning_out, wmsg)
        assign("warning_out", warning_out, envir = vachette_env)
        log_output(wmsg)
        f2.0 <- f2.0[1]
      }

      if(!is.na(f2.0[1]))
      {
        log_output(paste0("Inflec detected at x = ",f2.0[1],"\n"))
        add  <- data.frame(x=f2.0[1],type='inflec')
        lm   <- rbind(lm,add)
      }

      if(is.na(f2.0[1]))
      {
        log_output("No inflection point found")
      }
    }
  }

  # arrange by increasing x
  lm <- lm %>% arrange(x)

  return(lm)
}

#' Determine best shape-matching last.x of query
#'
#' @param ref reference curve
#' @param query query curve
#' @param lm.ref landmarks of reference curve
#' @param lm.query landmarks of query curve
#' @param ngrid number of points in last segment
#' @param scaling scaling of x-axis
#'
#' @noRd
#'
get.query.x.open.end <- function(ref,query,lm.ref,lm.query,ngrid=100,
                                 scaling='linear',polyorder=9) {

  x <- y <- NULL
  warning_out <- get("warning_out", envir = vachette_env)

  n.segment        <- nrow(lm.ref)-1
  n.segment.query  <- nrow(lm.query)-1
  ref.region       <- ref$region[1]
  query.region     <- query$region[1]

  if(n.segment != n.segment.query) {
    emsg <- "\nERROR:\nunequal number of segments reference/query\n*** Vachette processing stopped ***"
    log_output(emsg)
    stop(emsg)
  }

  # define starting landmark point (updated JL 240228, holds for first region only)
  if(n.segment==1 & ref.region  ==1) ref.x.start     <- 0
  if(n.segment==1 & query.region==1) query.x.start   <- 0
  if(n.segment>1 | ref.region   > 1) ref.x.start     <- lm.ref$x[n.segment]
  if(n.segment>1 | query.region > 1) query.x.start   <- lm.query$x[n.segment]
  # if(n.segment==1) ref.x.start     <- 0
  # if(n.segment==1) query.x.start   <- 0
  # if(n.segment>1)  ref.x.start     <- lm.ref$x[n.segment]
  # if(n.segment>1)  query.x.start   <- lm.query$x[n.segment]

  # Get open end
  t0          <- ref %>% filter(x >= ref.x.start)
  q0          <- query %>% filter(x >= query.x.start)

  #  ------ JL 20240227 -------------
  # Remove 2.5% smallest-y or highest-y from open end tail
  # Assuming monotonic increasing or decreasing curves

  frac.off <- 0.025  # Shorten curve by 2.5% in y-range

  # JL 240313 - no change
  if(t0$y[1] > t0$y[nrow(t0)])
  {
    diffy     <- t0$y[1] - t0$y[nrow(t0)]
    tmaxykeep <- t0$y[1]
    tminykeep <- frac.off*diffy + t0$y[nrow(t0)]
  }
  if(t0$y[1] < t0$y[nrow(t0)])
  {
    diffy     <- t0$y[nrow(t0)] - t0$y[1]
    tminykeep <- t0$y[1]
    tmaxykeep <- t0$y[nrow(t0)] - frac.off*diffy
  }
  if(q0$y[1] > q0$y[nrow(q0)])
  {
    diffy     <- q0$y[1] - q0$y[nrow(q0)]
    qmaxykeep <- q0$y[1]
    qminykeep <- frac.off*diffy + q0$y[nrow(q0)]
  }
  if(q0$y[1] < q0$y[nrow(q0)])
  {
    diffy     <- q0$y[nrow(q0)] - q0$y[1]
    qminykeep <- q0$y[1]
    qmaxykeep <- q0$y[nrow(q0)] - frac.off*diffy
  }

  t0 <- t0 %>% filter(y>=tminykeep,y<=tmaxykeep)
  q0 <- q0 %>% filter(y>=qminykeep,y<=qmaxykeep)


  # Check lengths:

  if(nrow(t0)<(polyorder+1))
  {
    e1 <- paste0("ERROR:\nLess than ",(polyorder+1)," gridpoints for reference last segment. Currently: ",nrow(t0)," only\n")
    log_output(e1)
    e2 <- paste("---> Please adjust typical curves (finer grid)\n")
    log_output(e2)
    e3 <- "*** Vachette processing stopped ***"
    log_output(e3)
    error_out <- c(e1, e2, e3)
    stop(error_out)
  }
  if(nrow(q0)<(polyorder+1))
  {
    e1 <- paste0("ERROR:\nLess than ",(polyorder+1)," gridpoints for query last segment. Currently: ",nrow(q0)," only\n")
    log_output(e1)
    e2 <- paste("---> Please adjust typical curves (finer grid)\n")
    log_output(e2)
    e3 <- ("*** Vachette processing stopped ***")
    log_output(e3)
    error_out <- c(e1, e2, e3)
    stop(error_out)
  }

  #  -----------------------------

  t1.tmp      <- t0  %>% mutate(x =  x-ref.x.start)
  q1.tmp      <- q0  %>% mutate(x =  x-query.x.start)
  # JL 240318 - to prevent issues??
  # t1.tmp      <- t0  %>% mutate(x =  x-t0$x[1])
  # q1.tmp      <- q0  %>% mutate(x =  x-q0$x[1])

  # Carry out optimization
  # print(paste0("Length t1.tmp ",nrow(t1.tmp)))
  # print(paste0("Length q1.tmp ",nrow(q1.tmp)))
  # print(paste0("Polyorder applied ",polyorder))

  t1.fit=lm(y~poly(x,polyorder,raw=F),data=t1.tmp)
  q1.fit=lm(y~poly(x,polyorder,raw=F),data=q1.tmp)

  # If t1.tmp / q1.tmp already contain (simulated) x=0, then copy:
  if(t1.tmp$x[1]==0) t1 <- t1.tmp
  if(q1.tmp$x[1]==0) q1 <- q1.tmp

  # Add landmark grid point at x=0
  if(t1.tmp$x[1]!=0)
  {
    add.t1   <- t1.tmp[1,]
    add.t1$x <- 0
    add.t1$y <- predict(t1.fit,newdata=add.t1)
    t1       <- rbind(add.t1,t1.tmp)
  }
  if(q1.tmp$x[1]!=0)
  {
    add.q1   <- q1.tmp[1,]
    add.q1$x <- 0
    add.q1$y <- predict(q1.fit,newdata=add.q1)
    q1       <- rbind(add.q1,q1.tmp)
  }

  # Initial x.scaling from segment preceeding last segment
  # if no landmarks, there is no previous scaling:
  # if (n.segment == 1) x.scaling.init <- 1
  # if (n.segment > 1)  x.scaling.init <- (lm.ref$x[n.segment] - lm.ref$x[n.segment-1])/((lm.query$x[n.segment] - lm.query$x[n.segment-1]))

  # Best estimate of y.scaling.start based on polynomial fits
  y.scaling.init.start <- t1$y[1]/q1$y[1]

  # y.scaling.init.end   <- y.scaling.init.start

  # JL 240227: Better to start with x.scaling.init=1 and y.scaling.init.end=1 to avoid extremes
  # JL 240325: assure positive x.scaling only
  x.scaling.init     <- log(1)
  y.scaling.init.end <- log(1)

  result <- optim(par=c(x.scaling.init, y.scaling.init.end), fn=funk,
                  y.scaling.start = y.scaling.init.start, t1=t1, q1=q1,
                  t1.fit=t1.fit, q1.fit=q1.fit, ngrid=10)

  ofv = result$value

  x.scaling.optim   <- exp(result$par[1])
  y.scale.optim.end <- exp(result$par[2])

  # message(paste("x.scaling = ",x.scaling.optim))
  if(x.scaling.optim<=0) {
    wmsg <- "WARNING:\nzero or negative open end x.scaling factor"
    warning_out <- c(warning_out, wmsg)
    assign("warning_out", warning_out, envir = vachette_env)
    log_output(wmsg)
  }

  log_output(paste0("Optimized last segment x.scaling ", signif(x.scaling.optim,4)))

  # Check if last.x fitted curve <= last.x template curve
  # (original translated curves)
  # If not, reverse the fitting and take inverse x.scaling
  # JL 240227 - try reverse also if x.scaling.optim out of range
  if (max(q1$x) * x.scaling.optim > max(t1$x) | x.scaling.optim <= 0)
  {
    log_output("  *** Reverse the last segment curve fitting ***")

    # Both x and y init scaling factors inverted
    y.scaling.init.start <- 1/y.scaling.init.start
    # JL 240325: assure positive x.scaling only
    x.scaling.init     <- log(1)
    y.scaling.init.end <- log(1)

    result <- optim(par=c(x.scaling.init, y.scaling.init.end), fn=funk,
                    y.scaling.start = y.scaling.init.start, t1=q1, q1=t1,
                    t1.fit=q1.fit, q1.fit=t1.fit, ngrid=10)

    ofv = result$value

    # Invert scaling result:
    x.scaling.optim      <- 1/exp(result$par[1])

    log_output(paste("x.scaling = ",x.scaling.optim))

    if(x.scaling.optim<=0) {
      emsg <- "\nERROR:\n(Still) zero or negative open end x.scaling factor\n*** Vachette processing stopped ***"
      log_output(emsg)
      stop(emsg)
    }

    log_output(paste0("Updated (reversed) last segment x.scaling ",signif(x.scaling.optim,4)))

  }

  # Just return updated query last.x (e.g. if x.scaling < 1, then "longer" query curve)
  ref.seg.length          <- lm.ref$x[n.segment+1] - lm.ref$x[n.segment]
  lm.query$x[n.segment+1] <- lm.query$x[n.segment] + (ref.seg.length/x.scaling.optim)

  return(lm.query)

}


#' 2-parameter optimization
#' With constant or linearly changing scaling factor
#'
#' @param x last.x of query
#' @param query.x.start first.x of query
#' @param query query curve
#' @param refGrid reference curve
#' @param queryMid.x middle ("ignored") landmark of query segment pair
#' @param refMid.x middle ("ignored") landmark of reference segement pair
#' @param ngrid number of points in last segment
#' @param ngrid number of points in last segment
#' @param scaling scaling of x-axis
#'
#' @noRd
#'
funk <- function(param,y.scaling.start,t1,q1,t1.fit,q1.fit,ngrid=10){

  # ---- Inits ------
  x <- NULL

  x.scaling     <- exp(param[1])
  y.scaling.end <- exp(param[2])

  # JL240318 Do not allow for negative x.scaling and y.scaling.end ...
  # x.scaling     <- abs(x.scaling)
  # y.scaling.end <- abs(y.scaling.end)

  # ----

  # Linear y-scaling equation
  icept.sc <- y.scaling.start  # Start is by default @ x.scaled=0
  slp.sc   <- (y.scaling.end - y.scaling.start)/(x.scaling*max(q1$x)-x.scaling*min(q1$x))

  # ----- Scaling the query curve ----

  # Copy template
  t2 <- t1

  # Scale query
  q2 <- q1
  q2$x.scaled <- q2$x * x.scaling
  q2$y.scaled <- q2$y * (icept.sc + slp.sc*q2$x.scaled)

  # ----- Create evenly space grid to calc y differences for OFV ----

  # Creat x-grid - use query domain
  Grid <- data.frame(x=seq(from=q2$x.scaled[1],to=max(q2$x.scaled),
                           length.out=ngrid),step=c(1:ngrid))

  # Template/Query - Use polynomal fits to get best y value approximation
  t3Grid   <- Grid %>% select(x)
  t3Grid$y <- predict(t1.fit,newdata=t3Grid)

  q3Grid   <- Grid %>% mutate(x.scaled = x, x = x/x.scaling)
  q3Grid$y <- predict(q1.fit,newdata=q3Grid)

  # Linear y-scaling:
  # q3Grid$y.scaled <- q3Grid$y * (icept.sc + slp.sc*q3Grid$x) # Wrong
  q3Grid$y.scaled <- q3Grid$y * (icept.sc + slp.sc*q3Grid$x.scaled)

  # ------ OFV -----

  fit=sum((t3Grid$y - q3Grid$y.scaled)**2)

  # print(paste(x.scaling,y.scaling.end,fit))

  return(fit)
}


#' Print \code{vachette_data}
#'
#' Print generic used to return information about \code{vachette_data} object
#'
#' @param x \code{vachette_data} object.
#' @param ... additional args.
#' @return \code{x} invisibly.
#' @export
print.vachette_data <- function(x, ...) {
  stopifnot(inherits(x, "vachette_data"))

  cat(sprintf("Model Name:\t\t%s", x$model.name), "\n")
  cat(sprintf("Covariate Names:\t%s", paste0(names(x$covariates), collapse=", ")), "\n")
  cat(sprintf("Reference Values:\t%s", paste0(paste0(names(x$covariates), "=", x$covariates), collapse = " , ")))

  # output error type
  invisible(x)
}


`%notin%` <- Negate(`%in%`)

.validate_columns <- function(mappings, obs.data, typ.data, sim.data = NULL, iiv.correction) {

  obs_req_cols <- c("ID", "OBS", "x")
  typ_req_cols <- c("ID", "PRED", "x")
  sim_req_cols <- c("REP", obs_req_cols)

  if (isTRUE(iiv.correction)) {
    obs_req_cols <- c(obs_req_cols, "PRED", "IPRED")
    #typ_req_cols <- c(typ_req_cols, "IPRED") #this should not be a requirement
  }

  obs_cols <- colnames(obs.data)
  typ_cols <- colnames(typ.data)
  sim_cols <- colnames(sim.data)

  if (is.null(mappings)) {
    if (any(obs_req_cols %notin% obs_cols)) {
      missing_cols <- setdiff(obs_req_cols, obs_cols)
      stop(
        paste0(missing_cols, collapse = " "),
        " column(s) not found in obs.data. Use the 'mappings' argument to manually specify column mappings."
        , call. = FALSE)
    }
    if (!is.null(sim.data) && any(sim_req_cols %notin% sim_cols)) {
      missing_cols <- setdiff(sim_req_cols, sim_cols)
      stop(
        paste0(missing_cols, collapse = " "),
        " column(s) not found in sim.data. Use the 'mappings' argument to manually specify column mappings."
        , call. = FALSE)
    }

    if (any(typ_req_cols %notin%  typ_cols)) {
      missing_cols <- setdiff(typ_req_cols, typ_cols)
      stop(
        paste0(missing_cols, collapse = " "),
        " column(s) not found in typ.data. Use the 'mappings' argument to manually specify column mappings."
        , call. = FALSE)
    }
  } else {
    diff_req_mappings_obs <- setdiff(obs_req_cols, names(mappings))
    cols_not_mapped_obs <- setdiff(diff_req_mappings_obs, obs_cols)

    obs_mappings <- mappings[names(mappings) %notin% c("REP")]
    sim_mappings <- mappings
    typ_mappings <- mappings[names(mappings) %notin% c("IPRED", "OBS", "REP")]


    if (length(cols_not_mapped_obs) > 0) {
      stop(
        paste0(cols_not_mapped_obs, collapse = " "),
        " column(s) not found in observed data and not provided in 'mappings'. Use the 'mappings' argument to manually specify columns mappings."
        , call. = FALSE)
    }

    if (any(obs_mappings %notin% obs_cols)) {
      missing_cols <- setdiff(obs_mappings, obs_cols)
      stop(
        paste0(missing_cols, collapse = " "),
        " column(s) provided in mappings argument are not found in obs.data."
        , call. = FALSE)
    }

    if (!is.null(sim.data)) {
      diff_req_mappings_sim <- setdiff(sim_req_cols, names(sim_mappings))
      cols_not_mapped_sim <- setdiff(diff_req_mappings_sim, sim_cols)

      if (length(cols_not_mapped_sim) > 0) {
        stop(
          paste0(cols_not_mapped_sim, collapse = " "),
          " column(s) not found in sim.data and not provided in 'mappings'. Use the 'mappings' argument to manually specify columns mappings."
          , call. = FALSE)
      }

      if (any(sim_mappings %notin% sim_cols)) {
        missing_cols <- setdiff(sim_mappings, sim_cols)
        stop(
          paste0(missing_cols, collapse = " "),
          " column(s) provided in mappings argument are not found in sim.data."
          , call. = FALSE)
      }
    }

    typ_mappings <- mappings[names(mappings) %notin% c("IPRED", "OBS", "REP")]

    if (any(typ_mappings %notin% typ_cols)) {
      missing_cols <- setdiff(typ_mappings, typ_cols)
      stop(
        paste0(missing_cols, collapse = " "),
        " column(s) provided in mappings argument are not found in typ.data."
        , call. = FALSE)
    }

    obs.data <- dplyr::rename(obs.data, dplyr::all_of(obs_mappings))
    if (!is.null(sim.data)) {
      sim.data <- dplyr::rename(sim.data, dplyr::all_of(sim_mappings))
    }
    typ.data <- dplyr::rename(typ.data, dplyr::all_of(typ_mappings))
  }

  return(
    list(
      obs.data = obs.data,
      sim.data = sim.data,
      typ.data = typ.data
    )
  )
}

.process_covariates <- function(covariates, obs.data) {
  if (is.null(names(covariates))) {
    names(covariates) <- covariates
    cov_names <- covariates
  } else {
    cov_names <- names(covariates)
  }
  for (i in seq_along(cov_names)) {
    cov_name <- cov_names[i]
    cov_ref_val <- covariates[i]
    # If name not supplied, expect value as name, then assign name
    if (cov_name == "") {
      cov_name <- cov_ref_val
      names(covariates)[i] <- cov_name
    }
    stopifnot(cov_name %in% colnames(obs.data))
    # Automatically set ref value to median/mode given cont/cat covariate type, if ref value not given
    # Account for case median is used with even number of values in dataset
    # - compute median for each subject first, then compute median from there

    if (cov_name == cov_ref_val) {
      if (suppressWarnings(all(is.na(as.numeric(unlist(obs.data[, cov_name])))))) {
        cov_ref_val <- .mode(unlist(obs.data[, cov_name]))
      } else {
        if (nrow(obs.data) %% 2 == 0) {
          #obs.data <- dplyr::slice(obs.data, -1) #before midpoint, drop first value
          obs.data <- dplyr::slice(obs.data, -nrow(obs.data)) #after midpoint, drop last value
        }
        cov_ref_val <- median(unlist(obs.data[, cov_name]))
      }
      stopifnot(cov_ref_val %in% as.character(unlist(obs.data[, cov_name])))
      covariates[i] <- cov_ref_val
    } else {
      if (cov_ref_val %notin% as.character(unlist(obs.data[, cov_name]))) {
        if (suppressWarnings(!is.na(as.numeric(cov_ref_val)))) {
          closest_index <- which.min(abs(obs.data[, cov_name] - as.numeric(cov_ref_val)))
          closest_value <- obs.data[, cov_name][closest_index]
          warning("Reference value of ", cov_ref_val, " for covariate '", cov_name, "' not found in data, setting reference value to ", closest_value,
                  call. = FALSE)
          covariates[i] <- closest_value
        } else {
          stop("Reference value of ", cov_ref_val, " for covariate ", cov_name, " not found in data", call. = FALSE)
        }
      }
    }
  }
  return(covariates)
}


.calculate_dose_number <- function(data, data_type = c("obs.data", "typ.data", "sim.data")) {

  II <- REP <- x <- ID <- AMT <- ADDL <- EVID <- NULL

  data_type <- match.arg(data_type)

  # If dose number provided, and column in data, use that
  if ("dosenr" %in% colnames(data)) {
    log_output(paste0(
      "`dosenr` column found in ",
      data_type,
      ", using `dosenr` column in data for corresponding ref.dosenr value")
    )
  } else if (all(c("ADDL", "II") %in% colnames(data))) {
    log_output(paste0(
      "`ADDL`, `II`, columns found in ",
      data_type,
      ", creating `dosenr` column in data for corresponding ref.dosenr value")
    )

    if ("AMT" %notin% colnames(data)) {
      data <- data %>%
        mutate(AMT = ifelse(II != 0, 1, 0))
    }
    if (data_type == "sim.data") {
      dose_data_full <- data %>%
        select(REP, x, ID, AMT, ADDL, II) %>%
        filter(ADDL > 0)
      dose_data_rep_split <- split(data, f = data$REP)
      for (j in seq_along(dose_data_rep_split)) {
        dose_data <- dose_data_rep_split[[j]]
        dosing <- list()

        for (i in 1:nrow(dose_data)) {
          dose_data_row <- slice(dose_data, i)
          dosing[[i]] <- data.frame(
            x = seq(
              from = dose_data_row$x,
              by = as.numeric(dose_data_row$II),
              length.out = as.numeric(dose_data_row$ADDL) + 1
            ),
            ID = dose_data_row$ID,
            AMT = dose_data_row$AMT
          )
        }
        dose_data_expanded <- do.call(rbind, dosing)
        data <- filter(data, ADDL == 0) %>%
          select(-ADDL,-II)

        data_new <- bind_rows(data, dose_data_expanded) %>%
          arrange(ID, x)

        data <- data_new %>%
          # JL 20-JUN-2023 fix dosenr problem: let AMT rec be first
          arrange(ID,x,-AMT) %>%
          # -------------------------------------------------------------------
        group_by(ID) %>%
          mutate(dosenr = cumsum(AMT != 0)) %>%
          ungroup() %>%
          filter(is.na(AMT) | AMT == 0) %>%
          select(-AMT)
        dose_data_rep_split[[j]] <- data
      }
      data <- do.call(rbind, dose_data_rep_split)
    } else {
      dose_data <- data %>%
        select(x, ID, AMT, ADDL, II) %>%
        filter(ADDL > 0)
      dosing <- list()
      for (i in 1:nrow(dose_data)) {
        dose_data_row <- slice(dose_data, i)
        dosing[[i]] <- data.frame(
          x = seq(
            from = dose_data_row$x,
            by = as.numeric(dose_data_row$II),
            length.out = as.numeric(dose_data_row$ADDL) + 1
          ),
          ID = dose_data_row$ID,
          AMT = dose_data_row$AMT
        )
      }
      dose_data_expanded <- do.call(rbind, dosing)
      data <- filter(data, ADDL == 0) %>%
        select(-ADDL,-II)

      data_new <- bind_rows(data, dose_data_expanded) %>%
        arrange(ID, x)

      data <- data_new %>%
        # JL 20-JUN-2023 fix dosenr problem: let AMT rec be first
        arrange(ID,x,-AMT) %>%
        # -------------------------------------------------------------------
      group_by(ID) %>%
        mutate(dosenr = cumsum(AMT != 0)) %>%
        ungroup() %>%
        filter(is.na(AMT) | AMT == 0) %>%
        select(-AMT)
    }
  }  else if ("EVID" %in% colnames(data)) {
    log_output(paste0(
      "`EVID` column found in ",
      data_type,
      ", creating `dosenr` column in data for corresponding ref.dosenr value")
    )
    if (data_type == "sim.data") {
      data <- data %>%
        # JL 20-JUN-2023 prevent dosenr problem: let EVID=1 rec be before EVID=0
        arrange(REP, ID,x,-EVID) %>%
        # -------------------------------------------------------------------
      group_by(REP, ID) %>%
        mutate(dosenr = cumsum(EVID == 1)) %>%
        ungroup() %>%
        filter(EVID == 0) %>%
        select(-EVID)
    } else {
      data <- data %>%
        # JL 20-JUN-2023 prevent dosenr problem: let EVID=1 rec be before EVID=0
        arrange(ID,x,-EVID) %>%
        # -------------------------------------------------------------------
      group_by(ID) %>%
        mutate(dosenr = cumsum(EVID == 1)) %>%
        ungroup() %>%
        filter(EVID == 0) %>%
        select(-EVID)
    }
  } else if ("AMT" %in% colnames(data)) {
    if (data_type == "sim.data") {
      data <- data %>%
        # JL 20-JUN-2023 prevent dosenr problem: let AMT rec be first
        arrange(REP, ID,x,-AMT) %>%
        # -------------------------------------------------------------------
      group_by(REP, ID) %>%
        mutate(dosenr = cumsum(AMT != 0)) %>% #if values NA: cumsum(!is.na(amt))
        ungroup() %>%
        filter(is.na(AMT) | AMT == 0) %>%
        select(-AMT)
    } else {
      data <- data %>%
        # JL 20-JUN-2023 prevent dosenr problem: let AMT rec be first
        arrange(ID,x,-AMT) %>%
        # -------------------------------------------------------------------
      group_by(ID) %>%
        mutate(dosenr = cumsum(AMT != 0)) %>% #if values NA: cumsum(!is.na(amt))
        ungroup() %>%
        filter(is.na(AMT) | AMT == 0) %>%
        select(-AMT)
    }
  } else {
    # final control flow should be understood as no dose information in the data, no ref.dosnr, message such e.g., sigmoid model
    log_output(paste0("`dosenr` column is required and not found in ", data_type," `dosenr` can be automatically calculated using `EVID`, `AMT`, or `ADDL/II` columns in the data. Use the 'mappings' argument if these columns are named differently in your input data."))
    warning("Setting `dosenr` column to 1 in ", data_type, call. = FALSE)
    data <- mutate(data, dosenr = 1)
  }

  return(data)
}


.mode <- function(x){
  uniqv <- unique(x)
  uniqv[which.max(tabulate(match(x, uniqv)))]
}

###################################################################
#                                                                 #
#                          MAIN                                   #
#                                                                 #
###################################################################

.calculate_transformations <- function(vachette_data,
                                       tol.end = 0.001,
                                       tol.noise = 1e-8,
                                       step.x.factor = 1.5,
                                       ngrid.fit = 100,
                                       window = 17,     # peak detection window - initial landmarks
                                       asymptote_right = TRUE,
                                       asymptote_left = FALSE,
                                       zero_asymptote_right = TRUE,
                                       zero_asymptote_left = TRUE) {



  # Collect all Vachette query curves and original/transformed observations (incl reference)
  ref.extensions.all   <- NULL
  curves.all           <- NULL
  curves.scaled.all    <- NULL
  obs.all              <- NULL
  sim.all              <- NULL
  obs.excluded         <- NULL

  my.ref.lm.all    <- NULL
  my.query.lm.all  <- NULL
  lm.all           <- NULL
  ucov <- noise.flag <- exclude <- seg <- x.shift.query <- x.scaling <- x.trans <- NULL
  x.shift.ref <- REP <- ID <- OBS <- y.scaling <- y.scaled <- x.scaled <- NULL

  tab.ucov <- vachette_data$tab.ucov
  stopifnot(!is.null(tab.ucov))
  typ.orig <- vachette_data$typ.orig
  stopifnot(!is.null(typ.orig))
  obs.orig <- vachette_data$obs.orig
  stopifnot(!is.null(obs.orig))
  sim.orig <- vachette_data$sim.orig
  stopifnot(!is.null(sim.orig))

  covariates <- vachette_data$covariates
  ref.region <- vachette_data$ref.region
  ref.dosenr <- vachette_data$ref.dosenr

  n_iter <- nrow(tab.ucov)

  vachette_env$summary_out <- list(errors = list(),
                                   warnings = list(),
                                   values = list()
  )

  vachette_env$summary_out$values[["n_iter"]] <- n_iter
  vachette_env$summary_out$values[["start_time"]] <- Sys.time()

  pb <- progress::progress_bar$new(
    format = "  vachette transformations [:bar] :percent :elapsed elapsed / :eta remaining",
    total = n_iter, clear = FALSE)

  init_time <- numeric(n_iter)
  end_time <- numeric(n_iter)

  for(i.ucov in c(1:n_iter)) {
    warning_out <- c()
    assign("warning_out", warning_out, envir = vachette_env)
    pb$tick()
    init_time[i.ucov] <- Sys.time()
    # Initialize
    REMOVE_OBS         <- FALSE
    remove.obs.after.x <- NULL
    obs.query.excluded <- NULL

    log_output('\n----------------------------------------------------------\n')
    log_output(paste('---------------------','i.ucov = ', i.ucov,'----------------------\n'))
    log_output(utils::capture.output(print(tab.ucov[i.ucov,], row.names = FALSE)))
    log_output('----------------------------------------------------------\n')

    # ---------------------------------------------------------------
    # ----  Define reference and query curves and observations    ---
    # ---------------------------------------------------------------

    # Ref: may change (by extensions), so define (again) every new combination
    filter_ref <-
      paste0(
        "!is.na(region) & ",
        paste0(
          'as.character(',
          names(covariates),
          ')',
          "==",
          "'",
          covariates,
          "'",
          collapse = " & "
        ),
        " & region == ref.region"
      )

    ref <- typ.orig %>% filter(eval(parse(text = filter_ref)))

    ref.ucov <- unique(ref$ucov)

    if(length(ref.ucov) != 1){
      emsg <- "\nERROR:\nmore than one reference covariate declared\n*** Vachette processing stopped ***"
      vachette_env$summary_out$errors[[i.ucov]] <- emsg
      log_output(emsg)
      stop(emsg)     # A single ref only possible
    }

    ref.region.type <- tab.ucov$region.type[tab.ucov$ucov==ref.ucov]

    # --------------------------------------------------------------------
    # JL 14-SEP-2023 - check if ref "closed", if so, keep next grid point of next region too

    if(ref.region.type=='closed')
    {
      ref.next <- NULL

      ref.region.next <- ref.region+1
      filter_ref_next <-
        paste0(
          "!is.na(region) & ",
          paste0(
            'as.character(',
            names(covariates),
            ')',
            "==",
            "'",
            covariates,
            "'",
            collapse = " & "
          ),
          " & region == ref.region.next"
        )

      ref.next      <- typ.orig %>% filter(eval(parse(text = filter_ref_next)))

      # Raise error if empty
      if(is.null(ref.next)) {
        emsg <- "\nERROR:\nNo typical curve for next reference region found\n*** Vachette processing stopped ***"
        vachette_env$summary_out$errors[[i.ucov]] <- emsg
        log_output(emsg)
        stop(emsg)     # A single ref only possible
      }

      ref.next.ucov <- unique(ref.next$ucov)

      if(length(ref.next.ucov) != 1) {
        emsg <- "\nERROR:\nlength ref.next.ucov != 1\n*** Vachette processing stopped ***"
        vachette_env$summary_out$errors[[i.ucov]] <- emsg
        log_output(emsg)
        stop(emsg)     # A single ref only possible
      }

      ref.next.grid.point.x <- ref.next$x[1]
    }

    # ------ Similar for selected query i.ucov  -----

    # query = Typical query curve
    query.region       <- tab.ucov$region[tab.ucov$ucov==i.ucov]
    query.region.type  <- tab.ucov$region.type[tab.ucov$ucov==i.ucov]

    # Typical query curve
    query <- typ.orig %>%
      filter(ucov == i.ucov) %>%
      mutate(ref  = tab.ucov$ref[i.ucov])    # Flag if reference

    # JL 14-SEP-2023 - check if query "closed", if so, keep next grid point of next region too
    if(query.region.type=='closed')
    {
      query.next <- NULL

      cov.query      <- tab.ucov[tab.ucov$ucov==i.ucov,]
      query.region.next <- cov.query$region + 1
      cov.query      <- cov.query %>% select(2:(dim(tab.ucov)[2]-3))

      filter_query_next <-
        paste0(
          "!is.na(region) & ",
          paste0(
            'as.character(',
            names(cov.query),
            ')',
            "==",
            "'",
            cov.query,
            "'",
            collapse = " & "
          ),
          " & region == query.region.next"
        )

      # Next typical query region
      query.next      <- typ.orig %>% filter(eval(parse(text = filter_query_next)))

      # Raise error if empty
      if(is.null(query.next)) {
        emsg <- "\nERROR:\nNo typical curve for next query region found\n*** Vachette processing stopped ***"
        vachette_env$summary_out$errors[[i.ucov]] <- emsg
        log_output(emsg)
        stop(emsg)
      }

      query.next.ucov <- unique(query.next$ucov)

      if(length(query.next.ucov) != 1) {
        emsg <- "\nERROR:\nlength query.next.ucov != 1\n*** Vachette processing stopped ***"
        vachette_env$summary_out$errors[[i.ucov]] <- emsg
        log_output(emsg)
        stop(emsg)
      }

      query.next.grid.point.x <- query.next$x[1]

    }

    # ---------------------------------------------------------------
    # ---- Remove noisy ends of ref and query curves (if present) ---
    # ---------------------------------------------------------------

    log_output("Remove noisy ends of the reference and query curves, if applicable")

    # Set ref tolerance to tol.noise * y_range * grid_step
    y_range <- max(ref$y)-min(ref$y)
    x_step  <- (max(ref$x)-min(ref$x))/(length(ref$x)-1)  # Maybe approx if irregular grid
    tol_cut <- tol.noise * y_range * x_step

    y.noise <- data.frame(noise.flag = (abs(ref$y - ref$y[length(ref$y)])) < tol_cut) %>%
      mutate(n=row_number())
    # Last non-noise data point:
    n.last.data <- max((y.noise %>% filter(noise.flag==F))$n)+1    # +1, since last data point always "noise"

    ref.no.noise <- ref[1:n.last.data,]

    # Set query tolerance to tol.noise * y_range * grid_step
    y_range <- max(query$y)-min(query$y)
    x_step  <- (max(query$x)-min(query$x))/(length(query$x)-1)  # Maybe approx if irregular grid
    tol_cut <- tol.noise * y_range * x_step

    y.noise <- data.frame(noise.flag = (abs(query$y - query$y[length(query$y)])) < tol_cut) %>%
      mutate(n=row_number())
    # Last non-noise data point:
    n.last.data <- max((y.noise %>% filter(noise.flag==F))$n)+1    # +1, since last data point always "noise"

    query.no.noise <- query[1:n.last.data,]

    # Replace
    ref   <- ref.no.noise
    query <- query.no.noise

    # ---------------------------------------------------------------
    #               ---- Interregion Gap Extrapolation ---
    # ---------------------------------------------------------------

    # Extrapolate region last-x region by one gridstep size if region.type = 'closed'
    # Currently simple extra x value with same y value ("horizontal" extrapolation - LOCF)
    if(ref.region.type == 'closed')
    {
      log_output(paste0("Reference region-",ref.region," gap extrapolation"))

      ref.curve.next       <- ref[dim(ref)[1],]
      ref.curve.next$x     <- ref.next.grid.point.x

      # Last 6 datapoints
      last6 <- ref %>% slice((n()-5):n()) %>% select(x,y)

      # Exponential fit
      y  = last6$y
      x  = last6$x
      exp.model <-lm(log(y) ~ x)

      # Extrapolate y
      ref.curve.next$y     <- exp(predict(exp.model,list(x=ref.curve.next$x)))

      ref <- rbind(ref,ref.curve.next)

    }
    if(query.region.type == 'closed')
    {
      log_output(paste0("Query region-",query.region," gap extrapolation"))

      query.curve.next       <- query[dim(query)[1],]
      query.curve.next$x     <- query.next.grid.point.x

      # Last 6 datapoints
      last6 <- query %>% slice((n()-5):n()) %>% select(x,y)

      # Exponential fit
      y  = last6$y
      x  = last6$x
      exp.model <-lm(log(y) ~ x)

      # Extrapolate y
      query.curve.next$y     <- exp(predict(exp.model,list(x=query.curve.next$x)))

      query <- rbind(query,query.curve.next)
    }

    # Observations dataframe
      obs.query   <- obs.orig %>%
        filter(ucov == i.ucov) %>%
        mutate(ref = tab.ucov$ref[i.ucov]) %>%  # Flag for reference data
        mutate(exclude = 0)

    # Create PRED for obs.query data points by linear interpolation only if not available
    # May be needed for reference extrapolation, see obs.query$PRED
    if (sum(grepl("PRED",names(obs.query)))==0) obs.query$PRED <- approx(x=query$x,y=query$y,xout=obs.query$x)$y

    # ---------------------------------------------------------------
    # ---- Flag when obs cannot be associated with query curve     --
    # ---------------------------------------------------------------

      if(i.ucov != ref.ucov)
      {
        log_output("Some obs cannot be associated with query curve")
        obs.query.excluded <- rbind(obs.query.excluded,
                                    obs.query %>%
                                      filter(x <= min(query$x) | x > max(query$x)) %>% # AL, 09 May 2024
                                      mutate(exclude=1) %>%
                                      mutate(reason="Missing typical curve"))
        # Keep the rest
        obs.query          <- obs.query %>%
          filter(x >= min(query$x) & x <= max(query$x))
      }
      if(i.ucov == ref.ucov)
      {
        out <- obs.query %>% filter(x <= min(query$x) | x > max(query$x)) # AL, 09 May 2024
        if(nrow(out)>0)
        {
          log_output("Some obs cannot be associated with reference curve")
          obs.query.excluded <- rbind(obs.query.excluded,
                                      obs.query %>%
                                        filter(x <= min(query$x) | x > max(query$x)) %>% # AL, 09 May 2024
                                        mutate(exclude=1) %>%
                                        mutate(reason="Missing typical curve"))
          # Keep the rest
          obs.query          <- obs.query %>%
            filter(x <= max(query$x))
        }
      }


    # ---------------------------------------------------------------
    # ----               Determine landmark positions             ---
    # ---------------------------------------------------------------

    # Tolerance to apply for finding maximums, minimums and inflection points
    # Small stepsize -> small tolerance
    # Large stepsize -> large tol
    tolapply <- tol.noise*(typ.orig$x[2]-typ.orig$x[1])

    # get ref landmarks
    my.ref.lm.init    <- get.x.multi.landmarks(ref$x,ref$y,w=window,tol=tolapply) #contiguous inflec cause issue?
    my.ref.lm.init$y  <- approx(ref$x,ref$y, xout=my.ref.lm.init$x)$y #interpolation
    # get query landmarks
    my.query.lm.init    <- get.x.multi.landmarks(query$x,query$y,w=window,tol=tolapply)
    my.query.lm.init$y  <- approx(query$x,query$y, xout=my.query.lm.init$x)$y

    # Number of landmarks
    n.lm.ref   <- nrow(my.ref.lm.init)
    n.lm.query <- nrow(my.query.lm.init)

    # New (temporary) landmark dataframe
    my.ref.lm.transf   <- my.ref.lm.init
    my.query.lm.transf <- my.query.lm.init

    # Check order of Landmarks up to second last (why up to second last?):
    # n.lm.check <- min(n.lm.ref, n.lm.query) - 1
    # JL 240327
    n.lm.check <- min(n.lm.ref, n.lm.query)
    # Match
    if(!identical(my.ref.lm.transf$type[1:n.lm.check],
                  my.query.lm.transf$type[1:n.lm.check]))
    {
      log_output("Reference landmarks")
      log_output(capture.output(print(my.ref.lm.transf)))
      log_output("Query landmarks")
      log_output(capture.output(print(my.query.lm.transf)))
      # stop("No matching order of landmarks between reference and query")
      # JL 240327
      log_output("!! No matching order of landmarks between reference and query !!")
    }

    # ---------------------------------------------------------------
    # ----                Check number of landmarks               ---
    # ----      Remove observations which cannot be transformed   ---
    # ---------------------------------------------------------------

    # If missing reference segment, we cannot transform
    if(n.lm.query > n.lm.ref)
    {
      log_output("Too few reference landmarks (reference curve too short)")

      # List of un-transformable (if any, corresponding to query segments after last ref landmark number)
      obs.query.out <- obs.query %>% filter(x >= my.query.lm.transf$x[n.lm.ref])

      if(nrow(obs.query.out)>0) log_output("Not all observations can be transformed [n.lm.query > n.lm.ref]")
      if(nrow(obs.query.out)>0) log_output(capture.output(print(obs.query.out)))

      # Exclude from transformation:
      obs.query$exclude[obs.query$x >= my.query.lm.transf$x[n.lm.ref]] <- 1

      # Keep same number of landmarks reference as for query
      my.query.lm.transf <- my.query.lm.transf[1:n.lm.ref,]
    }
    if(n.lm.query < n.lm.ref)
    {
      log_output("More reference landmarks than query landmarks")
      # Keep same number of landmarks reference as for query
      my.ref.lm.transf <- my.ref.lm.transf[1:n.lm.query,]
    }

    # Update
    n.lm.ref   <- nrow(my.ref.lm.transf)
    n.lm.query <- nrow(my.query.lm.transf)

    # Do for last segment of each curve
    if(n.lm.ref <= 2 | n.lm.query <= 2)
    {
      log_output("First and last landmark available only")
    }

    my.ref.lm.refined   <- my.ref.lm.transf
    my.query.lm.refined <- my.query.lm.transf

    # # No need to exclude if observations belong to reference
    if(i.ucov == ref.ucov)
    {
      # Update obs.query
      log_output(paste0("REFERENCE! No need to remove ",nrow(obs.query %>% filter(exclude==1)),
                     " observations for i.ucov=",i.ucov))

      # obs.query.excluded <- NULL
      obs.query          <- obs.query %>% filter(exclude==0)
    }

    # Remove observations if no matching reference segment available
    if(i.ucov != ref.ucov)
    {
      # Update obs.query
      log_output(paste0("Removing ",nrow(obs.query %>% filter(exclude==1))," observations for i.ucov=",i.ucov))

      obs.query.excluded <- rbind(obs.query.excluded,
                                  obs.query %>% filter(exclude==1) %>%
                                    mutate(reason="Missing corresponding ref segment"))
      obs.query          <- obs.query %>% filter(exclude==0)
    }

    # Update
    n.ref.landmarks   <- nrow(my.ref.lm.refined)
    n.query.landmarks <- nrow(my.query.lm.refined)

    # Count
    # if(n.ref.landmarks != n.query.landmarks)
    # {
    #   log_output("No matching number landmarks between reference and query")
    #
    #   if(n.ref.landmarks < n.query.landmarks)
    #   {
    #     # Shorten query curve
    #     log_output("Shortened query curve")
    #     my.query.lm.refined <- my.query.lm.refined[1:(n.ref.landmarks+1),]
    #     my.query.lm.refined$type[nrow(my.query.lm.refined)] <- 'end'
    #     # Remove obs after (new) "end"
    #     REMOVE_OBS <- TRUE
    #     remove.obs.after.x <- my.query.lm.refined$x[nrow(my.query.lm.refined)]
    #
    #     obs.removed <- obs.query %>% filter(x >  remove.obs.after.x)
    #     obs.query   <- obs.query %>% filter(x <= remove.obs.after.x)
    #
    #     log_output("New Query")
    #     print(my.query.lm.refined)
    #
    #     log_output("Observation ignored (not transformed):")
    #     print(obs.removed)
    #
    #   }
    #   if(n.ref.landmarks > n.query.landmarks)
    #   {
    #     # Shorten ref curve, no further action needed
    #     log_output("Shortened ref curve")
    #     my.ref.lm.refined <- my.ref.lm.refined[1:(n.query.landmarks+1),]
    #     my.ref.lm.refined$type[nrow(my.ref.lm.refined)] <- 'end'
    #
    #     log_output("New Ref")
    #     print(my.ref.lm.refined)
    #
    #   }
    # }

    # ---------------------------------------------------------------
    # ----          Calculate open end x.scaling factors          ---
    # ---------------------------------------------------------------

    scaling <- 'linear' # other options?

    # ---- Reference and query last x -----

    # JL 27-Jan-2024
    # Warn user in case we do not have landmarks. x=0 will be used as "surrogate" landmark.
    if(nrow(my.ref.lm.refined) == 0) {
      emsg <- "\nERROR:\n error in simulated data (1)\n*** Vachette processing stopped ***"
      vachette_env$summary_out$errors[[i.ucov]] <- emsg
      log_output(emsg)
      stop(emsg)
    }
    if(nrow(my.ref.lm.refined) == 1) {
      emsg <- "\nERROR:\nerror in simulated data (2)*** Vachette processing stopped ***"
      vachette_env$summary_out$errors[[i.ucov]] <- emsg
      log_output(emsg)
      stop(emsg)
    }

    if(nrow(my.ref.lm.refined) == 2)
    {
      log_output("No landmarks detected, x=0 assumed first landmark")
    }

    # JL 04-APR-2024
    polyorder <- 9
    SIGMOID = FALSE
    if (nrow(my.ref.lm.refined) == 3 & my.ref.lm.refined$type[2] == 'inflec')
    {
      log_output("Sigmoid-shaped curve detected")
      # Vachette will "mirror" first part of curve to estimate x.scaling
      SIGMOID = TRUE

      # log_output("Setting inflec points x=0 and x=log(10)")
      # if(my.query.lm.refined$x[2] < 0.1) my.query.lm.refined$x[2] = 0
      # if(my.query.lm.refined$x[2] > 2.2) my.query.lm.refined$x[2] = log(10)
      # if(my.ref.lm.refined$x[2] < 0.1)   my.ref.lm.refined$x[2] = 0
      # if(my.ref.lm.refined$x[2] > 2.2)   my.ref.lm.refined$x[2] = log(10)
    }

    scaling   <- 'linear'
    ngrid.fit <- 10
    # ADD tryCatch here to catch error, add it to appropriate icuv, then stop, can likely remove error logging within get.query.x.open.end
    my.query.lm     <- tryCatch({
      suppressWarnings(
        get.query.x.open.end(
          ref,
          query,
          my.ref.lm.refined,
          my.query.lm.refined,
          ngrid = ngrid.fit,
          scaling = scaling,
          polyorder = polyorder
        )
      )
    }, error = function(e) {
      emsg <- e$message
      vachette_env$summary_out$errors[[i.ucov]] <- emsg
      stop(emsg)
    })

    # If Sigmoid, also carry out mirrored curve for x.scaling of first segment
    if(SIGMOID)
    {
      log_output(" ... Now mirrored sigmoid ...")
      ref.mirror     <- ref
      ref.mirror$x   <- 2*my.ref.lm.refined$x[2] - ref$x    # Mirror around x=inflection point
      query.mirror   <- query
      query.mirror$x <- 2*my.query.lm.refined$x[2] - query$x    # Mirror around x=inflection point
      my.ref.lm.refined.mirror     <- my.ref.lm.refined
      my.ref.lm.refined.mirror$x   <- 2*my.ref.lm.refined$x[2] - my.ref.lm.refined$x
      my.ref.lm.refined.mirror$type[1] <- 'end'
      my.ref.lm.refined.mirror$type[3] <- 'start'
      my.ref.lm.refined.mirror         <- my.ref.lm.refined.mirror %>% arrange(x)

      my.query.lm.refined.mirror   <- my.query.lm.refined
      my.query.lm.refined.mirror$x <- 2*my.query.lm.refined$x[2] - my.query.lm.refined$x
      my.query.lm.refined.mirror$type[1] <- 'end'
      my.query.lm.refined.mirror$type[3] <- 'start'
      my.query.lm.refined.mirror         <- my.query.lm.refined.mirror %>% arrange(x)

      my.query.lm.mirror     <- suppressWarnings(
        get.query.x.open.end(ref.mirror,query.mirror,my.ref.lm.refined.mirror,my.query.lm.refined.mirror,
                             ngrid=ngrid.fit,scaling=scaling,polyorder=polyorder)
      )
      # Transfer obtained x.scaling to seg=1 of my.query.lm:

      # Obtained x.scaling:
      ref.mirror.seg.length      <- my.ref.lm.refined.mirror$x[3] - my.ref.lm.refined.mirror$x[2]
      query.mirror.seg.length    <- my.query.lm.mirror$x[3]       - my.query.lm.mirror$x[2]
      x.scaling.mirror           <- ref.mirror.seg.length/query.mirror.seg.length

      # Transfer to non-mirrored lm.query (calc with reference to inflection points)
      ref.seg.length   <- my.ref.lm.refined$x[2] - my.ref.lm.refined$x[1]
      my.query.lm$x[1] <- my.query.lm$x[2] - (ref.seg.length/x.scaling.mirror)

    }

    # x.scaling is based on unchanged ref
    my.ref.lm     <- my.ref.lm.refined

    # Recalc landmark y's
    my.ref.lm$y    <- approx(ref$x,ref$y,     xout=my.ref.lm$x)$y
    my.query.lm$y  <- approx(query$x,query$y, xout=my.query.lm$x)$y
    # Collect all landmarks
    my.ref.lm.all   <- rbind(my.ref.lm.all,   my.ref.lm   %>% mutate(i.ucov = i.ucov))
    my.query.lm.all <- rbind(my.query.lm.all, my.query.lm %>% mutate(i.ucov = i.ucov))

    lm.all <- rbind(lm.all, my.query.lm %>% mutate(ucov = i.ucov))

    log_output("\nReference")
    log_output(capture.output(print(my.ref.lm)))
    log_output("\nQuery")
    log_output(capture.output(print(my.query.lm)))

    # ---------------------------------------------------------------
    # ----              Map segments                              ---
    # ---------------------------------------------------------------

    nseg     <- dim(my.ref.lm)[1]-1

    # Initial assignment segments to typical curves, collect all query curves
    ref$seg   <- NA
    query$seg <- NA
    for(iseg in c(1:nseg))
    {
      if(iseg==1) # includes first datapoint
      {
        ref <- ref %>%
          mutate(seg = ifelse(x>=my.ref.lm$x[iseg] & x<=my.ref.lm$x[iseg+1],iseg,seg))
        query <- query %>%
          mutate(seg = ifelse(x>=my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1],iseg,seg))
      }
      if(iseg>1)
      {
        ref <- ref %>%
          mutate(seg = ifelse(x>my.ref.lm$x[iseg] & x<=my.ref.lm$x[iseg+1],iseg,seg))
        query <- query %>%
          mutate(seg = ifelse(x>my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1],iseg,seg))
      }
    }

    # ---------------------------------------------------------------
    # ----         Scale segments                                 ---
    # ---------------------------------------------------------------


    # 230418 New query.scaled curve data frame by contracting/expanding to x-range ref
    query.scaled <- NULL
    for(iseg in c(1:nseg))
    {
      # Ref segment length
      ref.seg.length <- my.ref.lm$x[iseg+1] - my.ref.lm$x[iseg]
      # Query segment length
      query.seg.length <- my.query.lm$x[iseg+1] - my.query.lm$x[iseg]

      # Scaling factor
      my.query.add <- query %>%
        filter((nseg==1 & x>=my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1]) |
                 (nseg>1 & x>my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1])) %>%
        mutate(seg=iseg) %>%
        mutate(y.scaling     = 1, y.scaled =1) %>%                # dummies, To remove
        mutate(x.shift.ref   = 0 - my.ref.lm$x[iseg]) %>%
        mutate(x.shift.query = 0 - my.query.lm$x[iseg]) %>%
        mutate(x.trans       = x + x.shift.query) %>%             # May be re-defined to my.ref.lm$x[iseg]
        mutate(x.scaling     = ifelse(query.seg.length>0,
                                      ref.seg.length/query.seg.length,
                                      NA)) %>%                    # NA causes troubles?
        mutate(x.scaled      = x.scaling*x.trans - x.shift.ref)   # Back to corresponding ref position

      query.scaled <- rbind(query.scaled,my.query.add)
    }

    # ---------------------------------------------------------------
    # ----     Add curve extension (if required)                  ---
    # ---------------------------------------------------------------

    # Default no extension for extrapolation
    EXTENSION       <- FALSE
    EXTENSION_RIGHT <- FALSE
    EXTENSION_LEFT  <- FALSE

    my.query.lm.extension <- NULL

    # Check if there are observation at all
    OBSERV <- ifelse(dim(obs.query)[1]>0,T,F)

    # If all observations are not before query last x, then extrapolate and add to reference region curve
    if (OBSERV)
    {
      # investigate
      # dummy <-  0
      # tmp.shift.last   <- -my.query.lm$x[2] + my.ref.lm$x[2]

      # print(my.ref.lm)
      # print(my.query.lm)
      #
      # tmp.scaling.last <- (my.query.lm$x[3] - my.query.lm$x[2]) / (my.ref.lm$x[3] - my.ref.lm$x[2])
      # # LAST query obs will be moved to:
      # last.scaled.query.x <- tmp.scaling.last*(max(obs.query$x) + tmp.shift.last)
      # print(last.scaled.query.x)
      # # Scaling of first seg:
      # tmp.scaling.first <- (my.query.lm$x[2] - my.query.lm$x[1]) / (my.ref.lm$x[2] - my.ref.lm$x[1])
      # # FIRST query obs will be moved to:
      # first.scaled.query.x <- tmp.scaling.last*(min(obs.query$x) + tmp.shift.last)
      # print(first.scaled.query.x)

      # Only needed when there are query observations found outside the typical query fitted last.x
      # [May also be activated if the typical curve has been cut off a bit to soon]

      # James: please move finctions to better place
      # Fixed x0 exp model
      exp.model2 <- function(parS,x,x0)
      {
        y = parS$A + parS$B*exp(-parS$k*(x-x0))
        return(y)
      }
      exp.model2.residuals <- function(p, observed, x, x0)
      {
        res <- observed - exp.model2(p, x, x0=x0)
        return(res)
      }

      if(max(my.query.lm$x) < max(obs.query$x))
      {

        EXTENSION       <- TRUE
        EXTENSION_RIGHT <- TRUE

        log_output('\n*** RIGHT Extension reference curve by exponential extrapolation ***')

        # ---- QUERY -----

        # max x observation on query
        max.obs.x <- max(obs.query$x)
        # y observation on query at x max
        max.obs.y <- obs.query$y[obs.query$x==max.obs.x]
        # y typical on query at x max
        max.obs.y.typ <- obs.query$PRED[obs.query$x==max.obs.x]

        # scaled max x observation: (x+x.shift)*scaling
        lastpoint        <- dim(query.scaled)[1]
        # JL 230621. Changes to scaled + translation to ref position
        # max.obs.x.scaled <- (max.obs.x + query.scaled$x.shift.query[lastpoint])*query.scaled$x.scaling[lastpoint]
        # (x.scaled           = x.scaling*x.trans - x.shift.ref)
        # with (x.trans       = x + x.shift.query)
        max.obs.x.scaled <- query.scaled$x.scaling[lastpoint]*(max.obs.x + query.scaled$x.shift.query[lastpoint]) - query.scaled$x.shift.ref[lastpoint]

        # QUERY
        cur.query.seg <- query$seg[dim(query[!is.na(query$seg),])[1]]
        # Adjust my.query.lm (always last point of last segment)
        my.query.lm.extension           <- my.query.lm
        my.query.lm.extension$x[nseg+1] <- max.obs.x
        my.query.lm.extension$y[nseg+1] <- max.obs.y.typ

        # ------- REFERENCE --------------

        # Same as for interregion gap extrapolation - multiple points
        cur.ref.seg   <- ref$seg[dim(ref[!is.na(ref$seg),])[1]]

        ref.grid.step.size   <- ref$x[dim(ref)[1]] - ref$x[dim(ref)[1]-1]

        # Increase x by ref grid steps up to last point
        # number of steps: ceiling((max.obs.x - max(ref$x))/ref.grid.step.size)
        # JL 230607. Correction, n.steps now defined as total number of steps, isteps as step counts
        # Incorrect: (we have to look at extension of query obs needed after x.scaling!)
        # i.steps.x          <- c(1:ceiling((max.obs.x - max(ref$x))/ref.grid.step.size))
        # JL 230621. Correction. Also shift query to ref position
        i.steps.x          <- c(1:ceiling((max.obs.x.scaled - max(ref$x))/ref.grid.step.size))

        n.steps.x          <- length(i.steps.x)
        steps.x            <- max(ref$x) + i.steps.x*ref.grid.step.size

        # ref.curve.next     <- ref %>% slice(rep(n(),max(i.steps.x)))
        # JL 230621
        ref.curve.next     <- ref %>% slice(rep(n(),n.steps.x))

        ref.curve.next$x   <- steps.x
        # Last x exactly matching max query obs x
        ref.curve.next$x[n.steps.x] <- max.obs.x.scaled

        # Last 6 datapoints
        last6 <- ref %>% slice((n()-5):n()) %>% select(x,y)

        last6$x  = as.numeric(last6$x)
        last6$y  = as.numeric(last6$y)

        last6.mirror <- last6 %>%
          mutate(x = -x) %>%
          arrange(x)   # Important!

        # JL 230717 Extrapolation suitable for non-zero asymptote

        # Simple Exponential fit
        if(asymptote_right)
        {
          log_output("Ref extension - assuming asymptote right")

          x <- last6$x
          y <- last6$y

          if(zero_asymptote_right)
          {
            log_output("  >> Applying ZERO asymptote reference curve RIGHT-side extrapolation")

            exp.model <-lm(log(y) ~ x)
            ref.curve.next$y <- exp(predict(exp.model,list(x=ref.curve.next$x)))
          }

          # Exponential fit for non-zero asymptote:
          if(!zero_asymptote_right)
          {
            log_output("  >> Applying non-zero asymptote reference curve RIGHT-side extrapolation")

            # 2. New exponent fit based on R = A + B*exp(-k*(t-t0)) with fixed t0
            xa <- (x[1] + x[2])/2
            xb <- (x[5] + x[6])/2
            Sa <- (y[1] + y[2])/2
            Sb <- (y[5] + y[6])/2

            x0    = xa
            Ainit = 0
            Binit = Sa-Sb
            if(Sa>=Sb) kinit = -log(Sb/Sa)/(xb-xa)
            if(Sa<Sb)  kinit = -log(Sa/Sb)/(xb-xa)
            Ainit = Sb   # Better estimate for asymptote value

            # (Add option for user to set asymptote to zero)

            p.init <- list(A=Ainit, B=Binit, k=kinit)

            ## Carry out exp model fit
            nls.out <-
              minpack.lm::nls.lm(
                par = p.init,
                fn = exp.model2.residuals,
                observed = y,
                x = x,
                x0 = x0,
                control = minpack.lm::nls.lm.control(maxiter = 500)
              )

            # Extrapolate y
            ref.curve.next$y <- exp.model2(as.list(coef(nls.out)), ref.curve.next$x, x0=x0)
          }
        }

        # Exponential fit for no-asymptote (can be programmed more efficiently):
        if(!asymptote_right)
        {
          log_output("  >> Applying NO-asymptote reference curve RIGHT-side extrapolation")

          # Mirrored x
          x <- last6.mirror$x
          y <- last6.mirror$y

          # 2. New exponent fit based on R = A + B*exp(-k*(t-t0)) with fixed t0
          xa <- (x[1] + x[2])/2
          xb <- (x[5] + x[6])/2
          Sa <- (y[1] + y[2])/2
          Sb <- (y[5] + y[6])/2

          x0    = xa
          Ainit = 0
          Binit = Sa-Sb
          if(Sa>=Sb) kinit = -log(Sb/Sa)/(xb-xa)
          if(Sa<Sb)  kinit = -log(Sa/Sb)/(xb-xa)
          Ainit = Sb   # Better estimate for asymptote value

          # (Add option for user to set asymptote to zero)

          p.init <- list(A=Ainit, B=Binit, k=kinit)

          ## Carry out exp model fit
          nls.out <-
            minpack.lm::nls.lm(
              par = p.init,
              fn = exp.model2.residuals,
              observed = y,
              x = x,
              x0 = x0,
              control = minpack.lm::nls.lm.control(maxiter = 500)
            )

          # Extrapolate y using mirrored x:
          ref.curve.next$y <- exp.model2(as.list(coef(nls.out)), -ref.curve.next$x, x0=x0)
        }


        # JL 230607 Assign extrapolation reference part to seg=-99 (right) for plotting purposes
        ref.curve.next$seg <- -99

        ref <- rbind(ref,ref.curve.next)

        # JL 230607
        # Keep extensions with associated query curve ucov
        ref.curve.next$ucov <- i.ucov  # Included reference curve extrapolation
        ref.extensions.all  <- rbind(ref.extensions.all,ref.curve.next)  # Included reference curve extrapolation

        # Check:
        if(cur.ref.seg != cur.query.seg) {
          emsg <- "\nERROR:\n segment number assignment in reference extrapolation block\n*** Vachette processing stopped ***"
          vachette_env$summary_out$errors[[i.ucov]] <- emsg
          log_output(emsg)
          stop(emsg)
        }
      }

      # Only needed when there are query observations found outside the typical query fitted last.x
      # [May also be activated if the typical curve has been cut off a bit to soon]
      # JL 240214
      if(min(my.query.lm$x) > min(obs.query$x))
      {

        EXTENSION      <- TRUE
        EXTENSION_LEFT <- TRUE

        log_output('\n*** LEFT Extension reference curve by exponential extrapolation ***')

        # ---- QUERY -----

        # min x observation on query
        min.obs.x <- min(obs.query$x)
        # y observation on query at x min
        min.obs.y <- obs.query$y[obs.query$x==min.obs.x]
        # y typical on query at x min
        min.obs.y.typ <- obs.query$PRED[obs.query$x==min.obs.x]

        # scaled min x observation: (x+x.shift)*scaling
        firstpoint        <- 1
        # JL 240214
        min.obs.x.scaled <- query.scaled$x.scaling[firstpoint]*(min.obs.x + query.scaled$x.shift.query[firstpoint]) - query.scaled$x.shift.ref[firstpoint]

        # QUERY (left extension --> seg=1)
        cur.query.seg <- 1
        # Adjust my.query.lm
        if(!EXTENSION_RIGHT) my.query.lm.extension <- my.query.lm
        if(EXTENSION_RIGHT)  my.query.lm.extension <- my.query.lm.extension
        my.query.lm.extension$x[1] <- min.obs.x
        my.query.lm.extension$y[1] <- min.obs.y.typ

        # ------- REFERENCE --------------

        # Same as for interregion gap extrapolation - multiple points
        cur.ref.seg   <- 1

        ref.grid.step.size   <- ref$x[2] - ref$x[1]

        # Decrease x by ref grid steps up to last point
        i.steps.x          <- c(1:ceiling((-min.obs.x.scaled + min(ref$x))/ref.grid.step.size))

        n.steps.x          <- length(i.steps.x)
        steps.x            <- min(ref$x) - (1+n.steps.x - i.steps.x)*ref.grid.step.size

        ref.curve.prev     <- ref %>% slice(rep(1,n.steps.x))

        ref.curve.prev$x   <- steps.x
        # First x exactly matching min query obs x
        ref.curve.prev$x[1] <- min.obs.x.scaled

        # First 6 datapoints
        first6 <- ref %>% slice(1:6) %>% select(x,y)

        first6$x  = as.numeric(first6$x)
        first6$y  = as.numeric(first6$y)

        # -- Mirror in order to use same extrapolation methodology/equations --
        first6.mirror <- first6 %>%
          mutate(x=-x) %>%
          arrange(x)
        # Mirror ref.curve.prev too, for extrapolation
        ref.curve.prev.mirror <- ref.curve.prev %>%
          mutate(x=-x) %>%
          arrange(x)
        # Mirror ref too, for extrapolation check
        ref.mirror <- ref %>%
          mutate(x=-x) %>%
          arrange(x)


        # JL 230717 Extrapolation suitable for non-zero asymptote

        # 1. Simple Exponential fit
        # Simple Exponential fit
        if(asymptote_left)
        {
          log_output("Ref extension - assuming asymptote left")

          x  = first6.mirror$x
          y  = first6.mirror$y

          if(zero_asymptote_left)
          {
            log_output("  >> Applying ZERO asymptote reference curve LEFT-side extrapolation")

            exp.model <-lm(log(y) ~ x)
            ref.curve.prev$y <- exp(predict(exp.model,list(x=ref.curve.prev$x)))
          }

          # 1. Exponential fit for non-zero asymptote:
          if(!zero_asymptote_left)
          {
            log_output("  >> Applying non-zero asymptote reference curve LEFT-side extrapolation")
            # Fixed x0 exp model
            exp.model2 <- function(parS,x,x0)
            {
              y = parS$A + parS$B*exp(-parS$k*(x-x0))
              return(y)
            }
            exp.model2.residuals <- function(p, observed, x, x0)
            {
              res <- observed - exp.model2(p, x, x0=x0)
              return(res)
            }

            # 2. New exponent fit based on R = A + B*exp(-k*(t-t0)) with fixed t0
            # left-hand side, non-zero asymptote
            xa <- (x[1] + x[2])/2
            xb <- (x[5] + x[6])/2
            Sa <- (y[1] + y[2])/2
            Sb <- (y[5] + y[6])/2

            # if Sa > Sb, then A<Sb
            # if Sa < Sb, then A>Sb
            # x0    = xa
            # if(Sa>Sb) Ainit = Sb*0.9         # Last 2 pnt
            # if(Sa<Sb) Ainit = Sb*1.1         # Last 2 pnt
            # Binit = Sa - Ainit               # Diff
            # kinit = -log(1 + (Sb - Sa)/Binit)/(xb-xa)

            x0    = xa
            Ainit = 0
            Binit = Sa-Sb
            if(Sa>=Sb) kinit = -log(Sb/Sa)/(xb-xa)
            if(Sa<Sb)  kinit = -log(Sa/Sb)/(xb-xa)
            Ainit = Sb   # Better estimate for asymptote value

            # (Add option for user to set asymptote to zero)

            # parS$A + parS$B*exp(-parS$k*(x-x0))
            p.init <- list(A=Ainit, B=Binit, k=kinit)

            ## Carry out exp model fit
            nls.out <-
              minpack.lm::nls.lm(
                par = p.init,
                fn = exp.model2.residuals,
                observed = y,
                x = x,
                x0 = x0,
                control = minpack.lm::nls.lm.control(maxiter = 500)
              )

            # B
            y.mirror <- exp.model2(as.list(coef(nls.out)), x, x0=x0)

            # Extrapolate y
            ref.curve.prev.mirror$y <- exp.model2(as.list(coef(nls.out)), ref.curve.prev.mirror$x, x0=x0)

          }

          # -- Back-mirror --
          first6 <- first6.mirror %>%
            mutate(x=-x) %>%
            arrange(x)
          # Mirror ref.curve.prev too, for extrapolation
          ref.curve.prev <- ref.curve.prev.mirror %>%
            mutate(x=-x) %>%
            arrange(x)
        }

        # Exponential fit for no-asymptote (can be programmed more efficiently):
        if(!asymptote_left)
        {
          log_output("  >> Applying NO-asymptote reference curve LEFT-side extrapolation")

          # Double mirror, so back to normal x-axis:
          x <- first6$x
          y <- first6$y

          # 2. New exponent fit based on R = A + B*exp(-k*(t-t0)) with fixed t0
          xa <- (x[1] + x[2])/2
          xb <- (x[5] + x[6])/2
          Sa <- (y[1] + y[2])/2
          Sb <- (y[5] + y[6])/2

          x0    = xa
          Ainit = 0
          Binit = Sa-Sb
          if(Sa>=Sb) kinit = -log(Sb/Sa)/(xb-xa)
          if(Sa<Sb)  kinit = -log(Sa/Sb)/(xb-xa)
          Ainit = Sb   # Better estimate for asymptote value

          # (Add option for user to set asymptote to zero)

          p.init <- list(A=Ainit, B=Binit, k=kinit)

          ## Carry out exp model fit
          nls.out <-
            minpack.lm::nls.lm(
              par = p.init,
              fn = exp.model2.residuals,
              observed = y,
              x = x,
              x0 = x0,
              control = minpack.lm::nls.lm.control(maxiter = 500)
            )

          # Extrapolate y
          ref.curve.prev$y <- exp.model2(as.list(coef(nls.out)), ref.curve.prev$x, x0=x0)

          # Implementation can be carried out similar to right-side no-asymptote extrapolation
        }

        # JL 230607 Assign extrapolation reference part to seg=-88 (left) for plotting purposes
        ref.curve.prev$seg <- -88

        ref <- rbind(ref.curve.prev,ref)

        # JL 230607
        # Keep extensions with associated query curve ucov
        ref.curve.prev$ucov <- i.ucov  # Included reference curve extrapolation
        ref.extensions.all  <- rbind(ref.extensions.all,ref.curve.prev)  # Included reference curve extrapolation

        # Check:
        if(cur.ref.seg != cur.query.seg) {
          emsg <- "\nERROR:\nsegment number assignment in reference extrapolation block\n*** Vachette processing stopped ***"
          vachette_env$summary_out$errors[[i.ucov]] <- emsg
          log_output(emsg)
          stop(emsg)
        }
      }

    }


    # Collect all typical curves with scaling factors and scaled x,y values
    query.scaled$y.scaled <- approx(ref$x,ref$y,xout=query.scaled$x.scaled)$y

    # JL 230607 Assign part of curve mapped to extrapolated part of reference curve, if required
    # Only extension part....
    # AL 18042024 left-continuous for step function
    # if(EXTENSION) query.scaled$seg <- approx(ref$x,ref$seg,
    #                                          xout=query.scaled$x.scaled,
    #                                          method = 'constant', f=1)$y
    # AL 22042024, assign the segment using reference position
    if(EXTENSION) {
      for (w in 1:(dim(my.ref.lm)[1]-1)){
        for (j in 1:dim(query.scaled)[1]) {
          if (query.scaled$x.scaled[j]>=my.ref.lm$x[w] & query.scaled$x.scaled[j]<=my.ref.lm$x[w+1]){
            query.scaled$seg[j]=w
          }
        }
      }
    }
    curves.all            <- rbind(curves.all,query)              # Included reference curve extrapolation
    curves.scaled.all     <- rbind(curves.scaled.all,query.scaled)

    if(dim(obs.query)[1]>0)
    {
      # Steps
      # a. (Re-)assign segments
      # b. Determine obs new x-position
      # c. Determine obs y-scaling ADDITIVE OR PROPORTIONAL

      # a. (Re-)assign segments, x.shift.ref, x.shift.query and x.scaling to observations and query curve (ref already done)
      obs.query$seg <- NA
      for(iseg in c(1:nseg))
      {
        # extended.query.x <- ifelse(length(tab.extension$x[tab.extension$seg==iseg])>0,
        #                            tab.extension$x[tab.extension$seg==iseg],NA)
        first.query.x    <- my.query.lm$x[iseg]
        # last.query.x     <- ifelse(!is.na(extended.query.x),extended.query.x,my.query.lm$x[iseg+1])
        if(!EXTENSION) last.query.x     <- my.query.lm$x[iseg+1]
        if(EXTENSION)  last.query.x     <- my.query.lm.extension$x[iseg+1]

        obs.query <- obs.query %>% mutate(seg = ifelse((nseg==1 & x>=first.query.x & x<=last.query.x) |
                                                         (nseg>1 & x>first.query.x & x<=last.query.x),iseg, seg))
        query     <- query     %>% mutate(seg = ifelse((nseg==1 & x>=first.query.x & x<=last.query.x) |
                                                         (nseg>1 & x>first.query.x & x<=last.query.x),iseg, seg))
      }

      # ---------------------------------------------------------------
      # ----     Vachettte x-transformation of observations         ---
      # ---------------------------------------------------------------

      # b. Determine obs new x-position
      obs.query.orig <- obs.query
      obs.query      <- NULL
      for(iseg in c(1:nseg))
      {
        # Original query.lm needed for scaling factor
        # Ref segment length
        ref.seg.length <- my.ref.lm$x[iseg+1] - my.ref.lm$x[iseg]
        # Query segment length
        query.seg.length <- my.query.lm$x[iseg+1] - my.query.lm$x[iseg]
        # Scaling factor
        if(!EXTENSION) obs.query.tmp <- obs.query.orig %>%
          filter((nseg==1 & x>=my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1]) |
                   (nseg>1 & x>my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1]))
        if(EXTENSION) obs.query.tmp <- obs.query.orig %>%
          filter((nseg==1 & x>=my.query.lm.extension$x[iseg] & x<=my.query.lm.extension$x[iseg+1]) |
                   (nseg>1 & x>my.query.lm.extension$x[iseg] & x<=my.query.lm.extension$x[iseg+1]))
        obs.query.add <- obs.query.tmp %>%
          mutate(seg=iseg) %>%
          mutate(y.scaling     = 1, y.scaled =1) %>%               # dummies, To remove
          mutate(x.shift.ref   = 0 - my.ref.lm$x[iseg]) %>%
          mutate(x.shift.query = 0 - my.query.lm$x[iseg]) %>%
          mutate(x.trans       = x + x.shift.query) %>%            # May be redefine to my.ref.lm$x[iseg]
          mutate(x.scaling     = ref.seg.length/query.seg.length) %>%
          mutate(x.scaled      = x.scaling*x.trans - x.shift.ref)  # Back to corresponding ref position

        obs.query <- rbind(obs.query,obs.query.add)
      }

      # can start here with sim.data, since everything before was x.scaling

      if (!is.null(sim.orig)) {
        sim.query.orig <- dplyr::inner_join(sim.orig %>% dplyr::select(REP, ID, x, y),
                                            obs.query.orig %>% dplyr::select(-REP, -y, -OBS),
                                            by = c("ID", "x"))
        sim.query      <- NULL
        for (iseg in c(1:nseg))
        {
          # Original query.lm needed for scaling factor
          # Ref segment length
          ref.seg.length <- my.ref.lm$x[iseg + 1] - my.ref.lm$x[iseg]
          # Query segment length
          query.seg.length <- my.query.lm$x[iseg + 1] - my.query.lm$x[iseg]
          # Scaling factor
          if (!EXTENSION)
            sim.query.tmp <- sim.query.orig %>%
            filter((nseg == 1 &
                      x >= my.query.lm$x[iseg] & x <= my.query.lm$x[iseg + 1]) |
                     (nseg > 1 &
                        x > my.query.lm$x[iseg] & x <= my.query.lm$x[iseg + 1])
            )
          if (EXTENSION)
            sim.query.tmp <- sim.query.orig %>%
            filter((
              nseg == 1 &
                x >= my.query.lm.extension$x[iseg] &
                x <= my.query.lm.extension$x[iseg + 1]
            ) |
              (
                nseg > 1 &
                  x > my.query.lm.extension$x[iseg] &
                  x <= my.query.lm.extension$x[iseg + 1]
              )
            )
          sim.query.add <- sim.query.tmp %>%
            mutate(seg = iseg) %>%
            mutate(y.scaling     = 1, y.scaled = 1) %>%               # dummies, To remove
            mutate(x.shift.ref   = 0 - my.ref.lm$x[iseg]) %>%
            mutate(x.shift.query = 0 - my.query.lm$x[iseg]) %>%
            mutate(x.trans       = x + x.shift.query) %>%            # May be redefine to my.ref.lm$x[iseg]
            mutate(x.scaling     = ref.seg.length / query.seg.length) %>%
            mutate(x.scaled      = x.scaling * x.trans - x.shift.ref)  # Back to corresponding ref position

          sim.query <- rbind(sim.query, sim.query.add)

        }
      }

      # ---------------------------------------------------------------
      # ----            Vachettte y-scaling of observations         ---
      # ---------------------------------------------------------------

      # ADDITIVE OR PROPORTIONAL

      # In case of proportional error model
      PROP_TR <- vachette_data$PROP_TR
      if(PROP_TR)
      {
        dummy <- 0
        # PROPORTIONAL - log scale addition. Take full ref and curve curves ....
        ref            <- ref   %>% mutate(ylog = ifelse(y>0,log(y),NA))
        query          <- query %>% mutate(ylog = ifelse(y>0,log(y),NA))

        # Exclude if y<=0:
        if(sum(obs.query <= 0) & i.ucov != ref.ucov)
        {
          # print(names(obs.query.excluded))
          # print(names(obs.query))
          log_output("\nZero observations cannot be transformed with proportional error model")
          obs.query.excluded <- rbind(obs.query.excluded,
                        obs.query %>% filter(y<=0) %>%
                          mutate(reason = "Less or equal to zero") %>%
                          mutate(exclude=1) %>%
                          select(-seg,-y.scaling,-y.scaled,-x.shift.ref,
                                 -x.shift.query,-x.trans,-x.scaling,-x.scaled))
          obs.query <- obs.query %>% filter(y>0)
        }

        obs.query      <- obs.query %>% mutate(ylog = ifelse(y>0,log(y),NA))

        # Position on query curve
        obs.query$ylog.query       <- ifelse(!is.na(obs.query$ylog),
                                             approx(x=query$x, y=query$ylog, xout = obs.query$x)$y,
                                             NA)
        # log difference obs.query to query curve
        obs.query$ylog.diff        <- obs.query$ylog - obs.query$ylog.query

        # Position on ref curve (position at x.scaled!!)
        obs.query$ylog.ref         <- ifelse(!is.na(obs.query$ylog),
                                             approx(x=ref$x, y=ref$ylog, xout = obs.query$x.scaled)$y,
                                             NA)
        # New position
        obs.query$ylog.scaled      <- obs.query$ylog.diff + obs.query$ylog.ref
        obs.query$y.scaled         <- exp(obs.query$ylog.scaled)

        if (!is.null(sim.orig)) {
          sim.query <-  sim.query %>% mutate(ylog = ifelse(y>0,log(y),NA))
          sim.query$ylog.query       <- ifelse(!is.na(sim.query$ylog),
                                               approx(x=query$x, y=query$ylog, xout = sim.query$x)$y,
                                               NA)
          # log difference sim.query to query curve
          sim.query$ylog.diff        <- sim.query$ylog - sim.query$ylog.query

          # Position on ref curve (position at x.scaled!!)
          sim.query$ylog.ref         <- ifelse(!is.na(sim.query$ylog),
                                               approx(x=ref$x, y=ref$ylog, xout = sim.query$x.scaled)$y,
                                               NA)
          # New position
          sim.query$ylog.scaled      <- sim.query$ylog.diff + sim.query$ylog.ref
          sim.query$y.scaled         <- exp(sim.query$ylog.scaled)
        }
      }

      # In case of additive error model
      ADD_TR <- vachette_data$ADD_TR
      if(ADD_TR)
      {
        # ADDITIVE - normal scale addition

        # Position on query curve
        obs.query$y.query       <- ifelse(!is.na(obs.query$y),
                                          approx(x=query$x, y=query$y, xout = obs.query$x)$y,
                                          NA)
        # difference obs.query to query curve
        obs.query$y.diff        <- obs.query$y - obs.query$y.query

        # Position on ref curve (position at x.scaled!!)
        obs.query$y.ref         <- ifelse(!is.na(obs.query$y),
                                          approx(x=ref$x, y=ref$y, xout = obs.query$x.scaled)$y,
                                          NA)
        # New position
        obs.query$y.scaled      <- obs.query$y.diff + obs.query$y.ref

        if (!is.null(sim.orig)) {
          sim.query$y.query       <- ifelse(!is.na(sim.query$y),
                                            approx(x=query$x, y=query$y, xout = sim.query$x)$y,
                                            NA)
          # difference sim.query to query curve
          sim.query$y.diff        <- sim.query$y - sim.query$y.query

          # Position on ref curve (position at x.scaled!!)
          sim.query$y.ref         <- ifelse(!is.na(sim.query$y),
                                            approx(x=ref$x, y=ref$y, xout = sim.query$x.scaled)$y,
                                            NA)
          # New position
          sim.query$y.scaled      <- sim.query$y.diff + sim.query$y.ref

        }
      }

      obs.all      <- rbind(obs.all,obs.query)
      if (!is.null(sim.orig)) {
        sim.all      <- rbind(sim.all,sim.query)
      }

    } # If there are observations to transform

    # Collect observations excluded from Vachette transformation
    obs.excluded <- rbind(obs.excluded,obs.query.excluded)

    log_output(paste0(nrow(obs.excluded)," observations excluded for all i.ucov"))

    ### Message Time
    # end_time[i.ucov] <- Sys.time()
    vachette_env$summary_out$warnings[[i.ucov]] <- get("warning_out", envir = vachette_env)
  }

  # Check additive/prop transformations
  obs.all$dist.add.orig <- NA
  obs.all$dist.add.transformed <- NA
  obs.all$dist.prop.orig <- NA
  obs.all$dist.prop.transformed <- NA
  # All ucov's
  n.ucov <- vachette_data$n.ucov
  for(i.ucov in c(1:n.ucov))
  {
    typ.curve <- typ.orig %>% filter(ucov==i.ucov)
    ref.curve <- typ.orig %>% filter(ucov==ref.ucov)
    obs       <- obs.all  %>% filter(ucov==i.ucov)

    # Additive distances to original curves
    obs.all$dist.add.orig[obs.all$ucov==i.ucov] <- approx(typ.curve$x,typ.curve$y,xout=obs$x)$y - obs$y
    # Vachette transformed - distances to ref curve
    obs.all$dist.add.transformed[obs.all$ucov==i.ucov] <- approx(ref.curve$x,ref.curve$y,xout=obs$x.scaled)$y - obs$y.scaled

    # Proportional distances to original curves (assume all y > 0)
    # obs.all$dist.prop.orig[obs.all$ucov==i.ucov]        <- log(approx(typ.curve$x,typ.curve$y,xout=obs$x)$y) - log(obs$y)
    obs.all$dist.prop.orig[obs.all$ucov==i.ucov]        <- approx(typ.curve$x,log(typ.curve$y),xout=obs$x)$y - log(obs$y)
    # Vachette transformed - distances to ref curve
    # obs.all$dist.prop.transformed[obs.all$ucov==i.ucov] <- log(approx(ref.curve$x,ref.curve$y,xout=obs$x.scaled)$y) - log(obs$y.scaled)
    obs.all$dist.prop.transformed[obs.all$ucov==i.ucov] <- approx(ref.curve$x,log(ref.curve$y),xout=obs$x.scaled)$y - log(obs$y.scaled)
  }

  obs.all <- obs.all %>%
    arrange(REP, ID, x)

  if (!is.null(sim.orig)) {
    sim.all$dist.add.orig <- NA
    sim.all$dist.add.transformed <- NA
    sim.all$dist.prop.orig <- NA
    sim.all$dist.prop.transformed <- NA
    # All ucov's
    n.ucov <- vachette_data$n.ucov
    for(i.ucov in c(1:n.ucov))
    {
      typ.curve <- typ.orig %>% filter(ucov==i.ucov)
      ref.curve <- typ.orig %>% filter(ucov==ref.ucov)
      sim       <- sim.all  %>% filter(ucov==i.ucov)

      # Additive distances to original curves
      sim.all$dist.add.orig[sim.all$ucov==i.ucov] <- approx(typ.curve$x,typ.curve$y,xout=sim$x)$y - sim$y
      # Vachette transformed - distances to ref curve
      sim.all$dist.add.transformed[sim.all$ucov==i.ucov] <- approx(ref.curve$x,ref.curve$y,xout=sim$x.scaled)$y - sim$y.scaled

      # Proportional distances to original curves (assume all y > 0)
      # sim.all$dist.prop.orig[sim.all$ucov==i.ucov]        <- log(approx(typ.curve$x,typ.curve$y,xout=obs$x)$y) - log(obs$y)
      sim.all$dist.prop.orig[sim.all$ucov==i.ucov]        <- approx(typ.curve$x,log(typ.curve$y),xout=sim$x)$y - log(sim$y)
      # Vachette transformed - distances to ref curve
      # sim.all$dist.prop.transformed[sim.all$ucov==i.ucov] <- log(approx(ref.curve$x,ref.curve$y,xout=obs$x.scaled)$y) - log(obs$y.scaled)
      sim.all$dist.prop.transformed[sim.all$ucov==i.ucov] <- approx(ref.curve$x,log(ref.curve$y),xout=sim$x.scaled)$y - log(sim$y.scaled)
    }

    sim.all <- sim.all %>%
      arrange(REP, ID, x)
  } else {
    sim.all <- NULL
  }



  return(
    list(
      obs.all      = tidyr::as_tibble(obs.all),
      sim.all      = sim.all,
      obs.excluded = obs.excluded,
      curves.all   = curves.all,
      curves.scaled.all  = curves.scaled.all,
      ref.extensions.all = ref.extensions.all,
      lm.all = lm.all,
      nseg   = nseg
    )
  )

}

