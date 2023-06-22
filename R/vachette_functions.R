
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

  data_type <- match.arg(data_type)

  # If dose number provided, and column in data, use that
  if ("dosenr" %in% colnames(data)) {
    message(
      "`dosenr` column found in ",
      data_type,
      ", using `dosenr` column in data for corresponding ref.dosenr value"
    )
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
      # JL 20-JUN-2023 fix dosenr problem: let AMT rec be first
      arrange(ID,x,-AMT) %>%
      # -------------------------------------------------------------------
      group_by(ID) %>%
      mutate(dosenr = cumsum(AMT != 0)) %>%
      ungroup() %>%
      filter(is.na(AMT) | AMT == 0) %>%
      select(-AMT)

  }  else if ("EVID" %in% colnames(data)) {
    message(
      "`EVID` column found in ",
      data_type,
      ", creating `dosenr` column in data for corresponding ref.dosenr value"
    )
    data <- data %>%
      # JL 20-JUN-2023 prevent dosenr problem: let EVID=1 rec be before EVID=0
      arrange(ID,x,-EVID) %>%
      # -------------------------------------------------------------------
      group_by(ID) %>%
      mutate(dosenr = cumsum(EVID == 1)) %>%
      ungroup() %>%
      filter(EVID == 0) %>%
      select(-EVID)
  } else if ("AMT" %in% colnames(data)) {
    data <- data %>%
      # JL 20-JUN-2023 prevent dosenr problem: let AMT rec be first
      arrange(ID,x,-AMT) %>%
      # -------------------------------------------------------------------
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

  return(data)
}


.mode <- function(x){
  uniqv <- unique(x)
  uniqv[which.max(tabulate(match(x, uniqv)))]
}

.calculate_transformations <- function(vachette_data,
                                       lm.refine = FALSE,
                                       tol.end = 0.001,
                                       tol.noise = 1e-8,
                                       step.x.factor = 1.5,
                                       ngrid.fit = 100,
                                       window = 17,     # Savitzky Golay smoothing - initial landmarks
                                       window.d1.refine = 7,     # Savitzky Golay smoothing - refine landmarks - first derivative
                                       window.d2.refine = 5,
                                       run_sim = FALSE) {

  # Collect all Vachette query curves and original/transformed observations (incl reference)
  ref.extensions.all   <- NULL
  curves.all           <- NULL
  curves.scaled.all    <- NULL
  obs.all              <- NULL

  my.ref.lm.all   <- NULL
  my.query.lm.all <- NULL
  lm.all <- NULL

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
  #lines 399-760
  # # Loop for all combinations of covariates
  for(i.ucov in c(1:dim(tab.ucov)[1])) {
    # A. ----- Define reference and query typical curves and observations to transform  ----------

    # Ref: may change (by extensions), so define (again) every new combination
    filter_query <-
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

    ref <- typ.orig %>% filter(eval(parse(text = filter_query)))

    ref.ucov <- unique(ref$ucov)

    if(length(ref.ucov) != 1) stop("Error: length ref.ucov != 1")     # A single ref only possible

    ref.region.type <- tab.ucov$region.type[tab.ucov$ucov==ref.ucov]
    # ref = Typical reference curve
    query.region       <- tab.ucov$region[tab.ucov$ucov==i.ucov]
    query.region.type  <- tab.ucov$region.type[tab.ucov$ucov==i.ucov]

    # Typical query curve
    query <- typ.orig %>%
      filter(ucov == i.ucov) %>%
      mutate(ref = tab.ucov$ref[i.ucov])    # Flag for reference

    # -------- INTERREGION GAP EXTRAPOLATION --------------

    # 230418 - extrapolate region last-x region by one gridstep size if region.type = 'closed'
    # Currently simple extra x value with same y value ("horizontal" extrapolation - LOCF)
    if(ref.region.type == 'closed')
    {
      #message(paste0("Reference region-",ref.region," gap extrapolation"))

      ref.grid.step.size   <- ref$x[dim(ref)[1]] - ref$x[dim(ref)[1]-1]
      ref.curve.next       <- ref[dim(ref)[1],]
      # Increase x:
      ref.curve.next$x     <- ref.curve.next$x + ref.grid.step.size

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
      #message(paste0("Query region-",query.region," gap extrapolation"))

      query.grid.step.size   <- query$x[dim(query)[1]] - query$x[dim(query)[1]-1]
      query.curve.next       <- query[dim(query)[1],]
      # Increase x:
      query.curve.next$x     <- query.curve.next$x + query.grid.step.size

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


    # Should we perform transformations on obs or sim data
    if (!run_sim){
      obs.query   <- obs.orig %>%
        filter(ucov == i.ucov) %>%
        mutate(ref = tab.ucov$ref[i.ucov])  # Flag for reference data
    } else {
      obs.query   <- sim.orig %>%
        filter(ucov == i.ucov)  %>%
        mutate(ref = tab.ucov$ref[i.ucov])  # Flag for reference data
    }

    # JL 230606
    # Create PRED for obs.query data points by linear interpolation
    # May be needed for reference extrapolation, see obs.query$PRED
    obs.query$PRED <- approx(x=query$x,y=query$y,xout=obs.query$x)$y

    #  B. -------- Get landmarks approximate positions ----

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

    # stop at this point, perhaps, allow user to visually inspect this plot, to assess what increase in stepsize
    # do not error out, return data that is available, and provide to user for plotting, and provide suggestions for updating argument values

    #Validation check
    #If y contiguous y values in series are the same, provide error, increase tol.noise
    #option 1: if multiple inflec in sequence, try decreasing tol.noise
    # optionally refine landmarks if multiple inflection points
    # automatically could raise tol.noise and/or window size incrementally to try to solve issue
    if(!lm.refine)
    {
      my.ref.lm.refined      <- my.ref.lm.init
      my.query.lm.refined    <- my.query.lm.init
    }
    if(lm.refine)
    {
      my.ref.lm.refined    <- refine.x.multi.landmarks(x=ref$x,y=ref$y,lm=my.ref.lm.init,tol=0.01*tolapply, w1=window.d1.refine,w2=window.d2.refine)
      my.ref.lm.refined$y  <- approx(ref$x, ref$y, xout=my.ref.lm.refined$x)$y

      my.query.lm.refined    <- refine.x.multi.landmarks(query$x,query$y,lm=my.query.lm.init,tol=0.01*tolapply, w1=window.d1.refine,w2=window.d2.refine)
      my.query.lm.refined$y  <- approx(query$x,query$y, xout=my.query.lm.refined$x)$y
    }

    # Check if query and ref have same landmarks in same order
    if(nrow(my.ref.lm.refined) != nrow(my.query.lm.refined) ||
       !sum(my.ref.lm.refined$type == my.query.lm.refined$type)>0)
      stop("No matching landmarks between reference and query")

    # C. ---------- Get last.x for reference and query typical curves ------------------

    scaling <- 'linear' # other options?

    # ------ Landmark last.x determination depends on segment type ---------

    # 1. ref and query both open ends ("default")
    # 2. ref is open end, query is closed --> get.query.open.end() for last x **ref** only
    # 3. ref is closed, query is open end --> get.query.open.end() for last x query only
    # 4. ref and query both closed        --> decide either fit ref.last.x or query.last.x

    query.region.type <- unique(query$region.type)
    if(length(query.region.type) != 1) stop("Error: length query.region.type != 1")

    # ---- Reference and query last x -----

    # 1. "Default" Find ref last x, Fit query last x
    if(ref.region.type == 'open' & query.region.type == 'open')
    {
      my.ref.lm   <- get.ref.x.open.end(ref$x,ref$y,my.ref.lm.refined,step.x.factor=step.x.factor,tol=tol.end)
      my.query.lm <- get.query.x.open.end(ref,query,my.ref.lm,my.query.lm.refined,ngrid=ngrid.fit,scaling=scaling)
    }

    # 2. Fit ref last x, Fix query last x
    if(ref.region.type == 'open' & query.region.type == 'closed')
    {
      my.ref.lm   <- get.query.x.open.end(query,ref,my.query.lm.refined,my.ref.lm.refined,ngrid=ngrid.fit,scaling=scaling)
      my.query.lm <- my.query.lm.refined
    }

    # 3. Fix ref last x, Fit query last x
    if(ref.region.type == 'closed' & query.region.type == 'open')
    {
      my.ref.lm   <- my.ref.lm.refined
      my.query.lm <- get.query.x.open.end(ref,query,my.ref.lm,my.query.lm.refined,ngrid=ngrid.fit,scaling=scaling)
    }

    # 4. Fix ref last x, Fit query last x *** OR **** reverse!!
    if(ref.region.type == 'closed' & query.region.type == 'closed')
    {
      # 1. Try fit query.last.x
      my.ref.lm1   <- my.ref.lm.refined
      my.query.lm1 <- suppressWarnings(
        get.query.x.open.end(ref,query,my.ref.lm1,my.query.lm.refined,ngrid=ngrid.fit,scaling=scaling)
      )

      # 2. Try fit ref.last.x
      my.query.lm2   <- my.query.lm.refined
      my.ref.lm2     <- suppressWarnings(
        get.query.x.open.end(query,ref,my.query.lm2,my.ref.lm.refined,ngrid=ngrid.fit,scaling=scaling)
      )

      # difference query/ref last.x with original
      diff.query.last.x.1 <- abs(my.query.lm1$x[dim(my.query.lm1)[1]] - my.query.lm.init$x[dim(my.query.lm.init)[1]])
      diff.ref.last.x.2   <- abs(my.ref.lm2$x[dim(my.ref.lm2)[1]]     - my.ref.lm.init$x[dim(my.ref.lm.init)[1]])

      # print(paste0("Diff last x query fit: ",diff.query.last.x.1))
      # print(paste0("Diff last x ref fit:   ",diff.ref.last.x.2))

      # Needs 3rd tolerance? Dependent on grid step size?
      ref.grid.step.size   <- ref$x[dim(ref)[1]] - ref$x[dim(ref)[1]-1]
      tol.fit.last.x       <- ref.grid.step.size*0.1
      # tol.fit.last.x       <- ref.grid.step.size*1.1
      # query last x fitted:
      if(diff.query.last.x.1>tol.fit.last.x & diff.ref.last.x.2 <= tol.fit.last.x)
      {
        my.ref.lm   <- my.ref.lm1
        my.query.lm <- my.query.lm1
      }
      # ref last x fitted:
      if(diff.ref.last.x.2>tol.fit.last.x & diff.query.last.x.1 <= tol.fit.last.x)
      {
        my.ref.lm   <- my.ref.lm2
        my.query.lm <- my.query.lm2
      }
      # Same curve
      if(diff.ref.last.x.2<=tol.fit.last.x & diff.query.last.x.1 <= tol.fit.last.x)
      {
        my.ref.lm   <- my.ref.lm1
        my.query.lm <- my.query.lm1
      }
      # Error
      if(diff.ref.last.x.2>tol.fit.last.x & diff.query.last.x.1 > tol.fit.last.x)
      {
        print(paste('diff.query.last.x.1',diff.query.last.x.1))
        print(paste('diff.ref.last.x.2',diff.ref.last.x.2))
        print(paste('tol.fit.last.x',tol.fit.last.x))

        stop("Error finding last.x for two closed segments")
      }
    }

    # JL 230608 - REMOVE: DUPLICATED CODE:
    # @James recent developments
    # Overrides above (above needs reprogramming)
    # Always FIX ref and FIT query to find best fitting x ("open"). Then extrapolate if required
    # if(ref.region.type == 'closed')
    # {
    #   my.ref.lm   <- my.ref.lm.refined
    #   my.query.lm <- get.query.x.open.end(ref,query,my.ref.lm,my.query.lm.refined,ngrid=ngrid.fit,scaling=scaling)
    # }

    # Recalc landmark y's
    my.query.lm$y  <- approx(query$x,query$y, xout=my.query.lm$x)$y
    my.ref.lm$y    <- approx(ref$x,ref$y, xout=my.ref.lm$x)$y

    # Collect all landmarks
    my.ref.lm.all   <- rbind(my.ref.lm.all,   my.ref.lm   %>% mutate(i.ucov = i.ucov))
    my.query.lm.all <- rbind(my.query.lm.all, my.query.lm %>% mutate(i.ucov = i.ucov))

    lm.all <- rbind(lm.all, my.query.lm %>% mutate(ucov = i.ucov))

    # D. ----- Define and map segments -----------

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

    # E. ---- x-scaled query segments -----------

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

    # F. ---- Check for need of extended **reference** curve by extrapolation -----------

    # Default no extension for extrapolation
    # tab.extension <- NULL
    EXTENSION <- FALSE

    # Check if there are observation at all
    OBSERV <- ifelse(dim(obs.query)[1]>0,T,F)

    # If all observations are not before query last x, then extrapolate and add to reference region curve
    if (OBSERV)
    {
      # Only needed when there are observation outside the segment last landmark
      if(max(my.query.lm$x) < max(obs.query$x))
      {

        EXTENSION <- TRUE

        message('*** Extension reference curve by simple exponential extrapolation ***')

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
        my.query.lm.extension <- my.query.lm
        my.query.lm.extension$x[nseg+1] <- max.obs.x
        my.query.lm.extension$y[nseg+1] <- max.obs.y.typ

        # REFERENCE
        # Reference curve has to be extrapolated up to max.obs.x.scaled
        # Last segment reference with last two typical datapoints extrapolated to scaled max x for obs

        # lastpoint        <- dim(ref)[1]
        # add.x <- max.obs.x.scaled
        # add.y <- approxExtrap(c(ref$x[lastpoint-1],ref$x[lastpoint]),
        #                       c(ref$y[lastpoint-1],ref$y[lastpoint]),
        #                       xout = add.x)$y
        #
        # ref.add <- ref %>%
        #   slice(1) %>%
        #   mutate(x=add.x,
        #          y=add.y) %>%
        #   mutate(seg=cur.ref.seg)
        # ref <- rbind(ref,ref.add)

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

        # Exponential fit
        y  = last6$y
        x  = last6$x
        exp.model <-lm(log(y) ~ x)

        # Extrapolate y
        ref.curve.next$y     <- exp(predict(exp.model,list(x=ref.curve.next$x)))

        # JL 230607 Assign extrapolation reference part to seg=-99 for plotting purposes
        ref.curve.next$seg <- -99

        ref <- rbind(ref,ref.curve.next)

        # JL 230607
        # Keep extensions with associated query curve ucov
        ref.curve.next$ucov <- i.ucov  # Included reference curve extrapolation
        ref.extensions.all  <- rbind(ref.extensions.all,ref.curve.next)  # Included reference curve extrapolation

        # Check:
        if(cur.ref.seg != cur.query.seg) stop("Error segment number assignment in reference extrapolation block")
      }
    }

    # Collect all typical curves with scaling factors and scaled x,y values
    query.scaled$y.scaled <- approx(ref$x,ref$y,xout=query.scaled$x.scaled)$y

    # JL 230607 Assign part of curve mapped to extrapolated part of reference curve, if required
    # Only extension part....
    if(EXTENSION) query.scaled$seg <- approx(ref$x,ref$seg,xout=query.scaled$x.scaled,
                                             method = 'constant')$y

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

      # obs.query still OK

      # ----- Steps 3: Vachette x transformation -----

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

      # ----- Steps 4: Vachette y transformation -----

      # c. Determine obs y-scaling ADDITIVE OR PROPORTIONAL

      # In case of proportional error model
      PROP_TR <- vachette_data$PROP_TR
      if(PROP_TR)
      {
        # PROPORTIONAL - log scale addition. Take full ref and curve curves ....
        ref            <- ref   %>% mutate(ylog = ifelse(y>0,log(y),NA))
        query          <- query %>% mutate(ylog = ifelse(y>0,log(y),NA))
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
      }

      obs.all <- rbind(obs.all,obs.query)

    } # If there are observations to transform

  }

  n.ucov <- vachette_data$n.ucov
  for(i.ucov in c(1:n.ucov))
  {
    typ.curve <- typ.orig %>% filter(ucov==i.ucov)
    ref.curve <- typ.orig %>% filter(ucov==ref.ucov)
    obs       <- obs.all  %>% filter(ucov==i.ucov)

    #Init columns
    obs.all$dist.add.orig <- NA
    obs.all$dist.add.transformed <- NA
    obs.all$dist.prop.orig <- NA
    obs.all$dist.prop.transformed <- NA
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

  return(
    list(
      obs.all = obs.all,
      curves.all = curves.all,
      curves.scaled.all = curves.scaled.all,
      ref.extensions.all = ref.extensions.all,
      lm.all = lm.all,
      nseg = nseg
    )
  )

}
