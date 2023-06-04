
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import prospectr
#' @importFrom Hmisc approxExtrap
#' @rawNamespace import(stats, except = c(filter, lag))
#'
NULL

`%>%`<- magrittr::`%>%`

#' Multi approx
#'
#' Find monotonic increasing and decreasing segments
#'
#' @param x x values
#' @param y y values
#' @param yout y value to approximate
#' @param tol tolerance for max and min identification
#'
#' @export
#'
multi.approx <- function(x,y,yout,tol=1e-9) {

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
#' @export
#'
get.x.multi.landmarks <- function(x,y,w=17,tol=1e-9) {
  # Using splinefun to determine derivatives
  # tol to be passed on to multi.approx function

  # Checks
  if(length(x) != length(y)) return(NA)
  if(length(x) <= 1)         return(NA)

  # Collect first data point
  f0.0 <- x[1]
  lm    <- data.frame(x=f0.0,type='start')

  # All first derivatives f1=0
  f1.sg <- savitzkyGolay(X = y,  m = 1,  p = 1,  w = w)
  x1.sg  <- x[ (1+(w-1)/2) : (length(x)-(w-1)/2) ]

  f1.0     <- multi.approx(x1.sg,f1.sg,yout=0,tol=tol)  # NA if no maximum
  if(sum(!is.na(f1.0))>=1)
  {
    add    <- data.frame(x=f1.0,type='max')
    lm     <- rbind(lm,add)
  }

  # All second derivatives f2=0
  f2.sg <- savitzkyGolay(X = y,  m = 2,  p = 2,  w = w)
  x2.sg  <- x[ (1+(w-1)/2) : (length(x)-(w-1)/2) ]

  f2.0 <- multi.approx(x2.sg,f2.sg,yout=0,tol=tol)   # NA if no inflection point
  if(sum(!is.na(f2.0))>=1)
  {
    add  <- data.frame(x=f2.0,type='inflec')
    lm   <- rbind(lm,add)
  }

  # Last data point - may be "asymptotic"
  f9.9 <- x[length(x)]
  add  <- data.frame(x=f9.9,type='end')
  lm   <- rbind(lm,add)

  # arrange by increasing x
  lm <- lm %>% arrange(x)

  return(lm)
}

#' Refine x multi landmarks
#'
#' @param x x values
#' @param y y values
#' @param lm landmarks
#' @param tol tolerance for max and min identification
#' @param w1 window size for f1
#' @param w2 window size for f2
#'
#' @export
#'
refine.x.multi.landmarks <- function(x,y,lm,tol=1e-9,w1=5,w2=7) {
  # Using splinefun to determine derivatives
  # Checks
  if(length(x) != length(y)) return(NA)
  if(length(x) <= 1)         return(NA)

  # Refined landmarks: First data point
  f0sub.0 <- x[1]
  lmsub   <- data.frame(x=f0sub.0,type='start')
  lmsub$y <- approx(x,y, xout=lmsub$x)$y

  # Loop to check/improve each landmark, ignore start and end points
  for(ilm in c(2:(dim(lm)[1]-1)))
  {
    # row number nearest points
    z    <- data.frame(row=c(1:length(x)),diff=abs(x - lm$x[ilm]))
    hit  <- z$row[z$diff==min(z$diff)]

    # Take take if near start or end of the curve
    w <- 17
    wmax <- ifelse((hit+17)>length(x),length(x),hit+w)
    wmin <- ifelse((hit-17)>0,hit-w,1)
    xsub <- x[wmin:wmax]
    ysub <- y[wmin:wmax]


    # **Maximum** so polynome now p=1+1=2
    if(lm$type[ilm] == 'max' | lm$type[ilm] == 'min')
    {
      w  <- w1
      po <- 2
      f1sub.sg  <- savitzkyGolay(X = ysub,  m = 1,  p = po,  w = w)
      x1sub.sg  <- xsub[ (1+(w-1)/2) : (length(xsub)-(w-1)/2) ]
      # How does f1 look like?
      pf1sub <- data.frame(x=x1sub.sg,y=f1sub.sg) %>% ggplot(aes(x=x,y=y))+geom_line()+geom_point()+
        geom_vline(xintercept =lm$x[ilm],col='red')+ggtitle('F1 derivative')
      plot(pf1sub + theme_bw())

      # remove y is present
      fsub1.0     <- multi.approx(x1sub.sg,f1sub.sg,yout=0,tol=tol)  # NA if no maximum
      if(sum(!is.na(fsub1.0))>=1)
      {
        add    <- data.frame(x=fsub1.0,type='max')
        add$y  <- approx(x,y, xout=add$x)$y
        lmsub  <- rbind(lmsub,add)
      }
    }

    # **Inflection** so polynome now p=2+1=3
    if(lm$type[ilm] == 'inflec')
    {
      w  <- w2
      po <- 3
      f2sub.sg  <- savitzkyGolay(X = ysub,  m = 2,  p = po,  w = w)
      x2sub.sg  <- xsub[ (1+(w-1)/2) : (length(xsub)-(w-1)/2) ]

      #f2
      pf2sub <- data.frame(x=x2sub.sg,y=f2sub.sg) %>% ggplot(aes(x=x,y=y))+geom_line()+geom_point()+
        geom_vline(xintercept =lm$x[ilm],col='red')+ggtitle('F2 derivative')
      plot(pf2sub + theme_bw())

      # remove y is present
      fsub2.0     <- multi.approx(x2sub.sg,f2sub.sg,yout=0,tol=tol)  # NA if no maximum
      if(sum(!is.na(fsub2.0))>=1)
      {
        add    <- data.frame(x=fsub2.0,type='inflec')
        add$y  <- approx(x,y, xout=add$x)$y
        lmsub  <- rbind(lmsub,add)
      }
    }
  }

  # Last data point - may be "asymptote"
  fsub9.9 <- x[length(x)]
  add     <- data.frame(x=fsub9.9,type='end')
  add$y   <- approx(x,y, xout=add$x)$y
  lmsub   <- rbind(lmsub,add)

  # arrange by increasing x
  lmsub <- lmsub %>% arrange(x)

  return(lmsub)
}

#' Determine last x value for open end reference curve
#'
#' @param x x values
#' @param y y values
#' @param lm landmarks
#' @param step.x.factor factor to increase step.x
#' @param tol tolerance for max and min identification
#'
#' @export
#'
get.ref.x.open.end <- function(x,y,lm,step.x.factor=1.5,tol=0.01) {
  # Must contain type="end" as last entry"
  if(lm$type[length(lm$type)] != "end") return("No \"end\" in landmarks")

  # landmarks min/max (incl. first, excl. last data point)
  min.lm <- min(lm$y[lm$type!='end'])
  max.lm <- max(lm$y[lm$type!='end'])
  # If one landmark only, then range = absolute difference from zero
  if (dim(lm[lm$type!='end',])[1] == 0) stop('ERROR in reference landmark dataframe, no \"end\" defined')
  if (dim(lm[lm$type!='end',])[1] == 1) ytol <- tol*(max(y) - min(y))  # Start and end point only available in lm
  if (dim(lm[lm$type!='end',])[1] >  1) ytol <- tol*(max.lm - min.lm)

  # Last segment, all data points from last landmark point to end
  xlast      <- x[x > lm$x[length(lm$x)-1]]         # last.x first point from last landmark
  xlastNext  <- step.x.factor*xlast                 # next.x which is step.x.factor further out to test for ytol difference
  ylast      <- approx(x,y,xout=xlast)$y            # Calc y for last.x
  ylastNext  <- approx(x,y,xout=xlastNext,rule=1)$y # Calc y for next.x. NA's if outside domain
  notna      <- !is.na(ylastNext)

  # Remove NA's
  xlast      <- xlast[notna]
  xlastNext  <- xlastNext[notna]
  ylast      <- ylast[notna]
  ylastNext  <- ylastNext[notna]

  # Vector of differences
  ydiff      <- abs(ylast - ylastNext)

  # Find xNext position
  xsol <- multi.approx(xlast,ydiff,yout=ytol)   # tol.end
  # Note: For emax we have 2 solutions: one at start when curve is slowly increasing, one at end when curve is slowly flattening out
  # Always take open end = last solution
  xhit <- xsol[length(xsol)]

  # TURNED OFF, OK?
  # if(is.na(xhit)) return(paste0("No solution within tol ",tol," with step.x.factor ",step.x.factor," from get.ref.x.open.end()"))

  yhit <- approx(x,y,xout=xhit)$y

  # Replace "end"
  lm$x[lm$type=='end']  = xhit
  lm$y[lm$type=='end']  = yhit

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
#' @export
#'
get.query.x.open.end <- function(ref,query,lm.ref,lm.query,ngrid=100,scaling='linear') {
  # Ref last segment represented by ngrid points
  refGrid      <- data.frame(x=seq(from=lm.ref$x[length(lm.ref$x)-1],
                                   to=lm.ref$x[length(lm.ref$x)],length.out=ngrid),
                             step=c(1:ngrid))
  refGrid$y    <- approx(ref$x,ref$y,xout=refGrid$x)$y

  # Query last segment
  query.x.start <- lm.query$x[length(lm.query$x)-1]
  query.x.last  <- lm.query$x[length(lm.query$x)]

  # Optimizer
  result <- optimize(get.y.diff.ref.query,
                     interval = c(query.x.start,query.x.last),
                     query.x.start=query.x.start,
                     query=query,
                     refGrid=refGrid,
                     ngrid=ngrid,
                     scaling=scaling)

  # Check if at boundary
  if(abs(result$minimum-query.x.start)<1e-3) message("WARNING: Query minimum x last near landmark before end in get.query.x.open.end")
  if(abs(result$minimum-query.x.last)<1e-3)  message("WARNING: Query minimum x last near last x value provided/explored in get.query.x.open.end")
  if(abs(result$minimum-query.x.last)<1e-3)  message("         Increase the tolerance or increase length simulated profile or decrease x-step factor")

  lm.query$x[lm.query$type=='end'] <- result$minimum
  lm.query$y[lm.query$type=='end'] <- approx(query$x,query$y,xout=result$minimum)$y

  return(lm.query)
}

#' With constant or linearly changing scaling factor
#'
#' @param x last.x of query
#' @param query.x.start first.x of query
#' @param query query curve
#' @param refGrid reference curve
#' @param ngrid number of points in last segment
#' @param scaling scaling of x-axis
#'
#' @export
#'
get.y.diff.ref.query <- function(x, query.x.start, query, refGrid, ngrid=100, scaling='linear') {
  # Query segment represented by ngrid points
  queryGrid   <- data.frame(x=seq(from=query.x.start,to=x,length.out=ngrid),
                            step=c(1:ngrid))
  queryGrid$y <- approx(query$x,query$y,xout=queryGrid$x)$y

  # Scale query by linear changing scaling factor
  if(scaling == 'linear')
  {
    scale.df          <- data.frame(x=c(1,ngrid),y=c(refGrid$y[1]/queryGrid$y[1],refGrid$y[ngrid]/queryGrid$y[ngrid]))
    queryGrid$scaling  <- approx(scale.df$x, scale.df$y, xout=c(1:ngrid))$y
    queryGrid$y.scaled <- queryGrid$scaling*queryGrid$y
  }
  # Scale query by constant scaling factor
  if(scaling == 'constant')
  {
    queryGrid$y.scaled <- refGrid$y[1]/queryGrid$y[1]*queryGrid$y
  }
  # No scaling
  if(scaling == '1')
  {
    queryGrid$y.scaled <- queryGrid$y
  }

  # Scale query by linear changing scaling factor
  if(scaling == 'quadratic')
  {
    scale.df           <- data.frame(x=c(1,ngrid),y=c(refGrid$y[1]/queryGrid$y[1],refGrid$y[ngrid]/queryGrid$y[ngrid]))
    queryGrid$scaling  <- approx(scale.df$x, scale.df$y, xout=c(1:ngrid))$y
    queryGrid$y.scaled <- queryGrid$scaling*queryGrid$y
  }
  # Scale query by constant scaling factor
  if(scaling == 'constant')
  {
    queryGrid$y.scaled <- refGrid$y[1]/queryGrid$y[1]*queryGrid$y
  }
  # No scaling
  if(scaling == '1')
  {
    queryGrid$y.scaled <- queryGrid$y
  }
  # Ignore gridpoints for y-values equal to zero
  if(scaling == 'linear.not.0')
  {
    # Ignore zeroes for scaling
    zeroes <- (refGrid$y == 0 | queryGrid$y==0)
    df  <- data.frame(x=refGrid$x, yref=refGrid$y, yquery=queryGrid$y)
    df2 <- df[!zeroes,]

    df2$scale <- df2$yref/df2$yquery
    df$scale  <- approx(x=df2$x,y=df2$scale,xout=df$x,rule=2)$y  # rule=2: If outside interval, then closest value
    scale.df  <- data.frame(x=c(1,ngrid),
                            y=c(df$scale[1],df$scale[ngrid]))

    # Linear scaling between first and last point
    queryGrid$scaling  <- approx(scale.df$x, scale.df$y, xout=c(1:ngrid), rule=2)$y
    queryGrid$y.scaled <- queryGrid$scaling*queryGrid$y
  }

  # To minimize y-differences of reference with query-after-scaling
  diff <- sum((refGrid$y-queryGrid$y.scaled)**2)/ngrid

  return(diff)
}


#' Print \code{vachette_data}
#'
#' Print generic used to return information about \code{vachette_data} object
#'
#' @param x An \code{vachette_data}.
#' @param ... Additional args.
#' @return Returns \code{x} invisibly.
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
  sim_req_cols <- c("isim", obs_req_cols)

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

    obs_mappings <- mappings[names(mappings) %notin% c("isim")]
    sim_mappings <- mappings
    typ_mappings <- mappings[names(mappings) %notin% c("IPRED", "OBS", "isim")]


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

    typ_mappings <- mappings[names(mappings) %notin% c("IPRED", "OBS", "isim")]

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
        stop("Reference value of ", cov_ref_val, " for covariate ", cov_name, " not found in data", call. = FALSE)
      }
    }
  }
  return(covariates)
}


.calculate_dose_number <- function(data, ref.dosenr, data_type = c("obs.data", "typ.data", "sim.data")) {

  data_type <- match.arg(data_type)

  # If dose number provided, and column in data, use that
  if ("dosenr" %in% colnames(data)) {
    message(
      "`dosenr` column found in ",
      data_type,
      ", using `dosenr` column in data for corresponding ref.dosenr value"
    )
  }  else if ("EVID" %in% colnames(data)) {
    message(
      "`EVID` column found in ",
      data_type,
      ", creating `dosenr` column in data for corresponding ref.dosenr value"
    )
    data <- data %>%
      group_by(ID) %>%
      mutate(dosenr = cumsum(EVID == 1)) %>%
      ungroup() %>%
      filter(EVID == 0) %>%
      select(-EVID)
  } else if (all(c("ADDL", "II") %in% colnames(data))) {
    message(
      "`ADDL`, `II`, columns found in ",
      data_type,
      ", creating `dosenr` column in data for corresponding ref.dosenr value"
    )

    if ("AMT" %notin% colnames(data)) {
      data <- data %>%
        mutate(AMT = ifelse(II != 0, 1, 0))
    }
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
      group_by(ID) %>%
      mutate(dosenr = cumsum(AMT != 0)) %>%
      ungroup() %>%
      filter(is.na(AMT) | AMT == 0) %>%
      select(-AMT)

  } else if ("AMT" %in% colnames(data)) {
    data <- data %>%
      group_by(ID) %>%
      mutate(dosenr = cumsum(AMT != 0)) %>% #if values NA: cumsum(!is.na(amt))
      ungroup() %>%
      filter(is.na(AMT) | AMT == 0) %>%
      select(-AMT)
  } else {
    # final control flow should be understood as no dose information in the data, no ref.dosnr, message such e.g., sigmoid model
    message("`dosenr` column is required and not found in ", data_type," `dosenr` can be automatically calculated using `EVID`, `AMT`, or `ADDL/II` columns in the data. Use the 'mappings' argument if these columns are named differently in your input data.")
    warning("Setting `dosenr` column to 1 in ", data_type, call. = FALSE)
    data <- mutate(data, dosenr = 1)
  }
  # ensure value supplied to ref.dosenr exists inside dosenr column
  if (ref.dosenr %notin% unique(data[["dosenr"]])) {
    stop("ref.dosenr value of ", ref.dosenr, " not found in dosenr column in ", data_type, call. = FALSE)
  }

  return(data)
}


.mode <- function(x){
  uniqv <- unique(x)
  uniqv[which.max(tabulate(match(x, uniqv)))]
}

