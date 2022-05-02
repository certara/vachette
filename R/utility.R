# -------------------------------------------------------------------------
#  Sponsor           : Merck
#  Project Number    : MRK-FTE-VRSV-640
#  Project Root Path : Local
# -------------------------------------------------------------------------
#  Program : vachette-functions-v5.R
#  Author  : Jos Lommerse - Certara
#  Date    : 17 January 2022
#  Purpose : Vachette functions for Visualization of Analyses with Covariates
# -------------------------------------------------------------------------
#  Software : R version 4.1.2 (2021-11-01)
#  Platform : ThinkPad P15s Gen 2 Model 20W6008AUS
#  Environment : Windows-10, 11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz
# -------------------------------------------------------------------------

#' @import ggplot2
#'
NULL

`%>%`<- magrittr::`%>%`

#' multi.approx
#'
#' @param x
#' @param y
#' @param yout
#' @param tol
#'
#' @return Return multiple x hits of (x,y) series for y-value yout and returns NA if no hit
#' @export
#'
#' @examples
multi.approx <- function(x,y,yout,tol=1e-9)
{

  yrange <- max(y,na.rm=T)-min(y,na.rm=T)
  tolval <- tol*yrange

  # Checks
  if(length(x) != length(y)) return(NA)
  if(length(x) <= 1)         return(NA)

  # Find monotonic increasing and decreasing segments which cross yout
  # Better to base on previous two and next two data points (two to prevent impact noise)
  df      <- data.frame(x=x,y=y)

  # "True" max or min if both 2 preceeding and 2 proceeding data points
  #  indicate that max or min. With sufficient difference between all 5 point (tolval)
  # initialize
  df$max    <- 0
  df$min    <- 0
  df$min[1] <- 1   # Can be min or max, does not matter
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
  # Last data point not considered as min or max
  # df$min[dim(df)[1]] <- 1   # Can be min or max, does not matter

  # all potential segments
  iseg      <- 1
  df$seg[1] <- iseg
  for(i in c(2:dim(df)[1]))
  {
    if(df$min[i] == 1 | df$max[i] == 1) iseg <- iseg+1 # start new segm
    df$seg[i] <- iseg
  }
  df[df$min==1 | df$max==1,]
  df %>% ggplot(aes(x=x,y=y,col=factor(seg)))+geom_line()+geom_point()+
    geom_point(data=df %>% filter(df$max==1),col='black',pch=1,size=3)+
    geom_point(data=df %>% filter(df$min==1),col='black',pch=1,size=3)+
    coord_cartesian(xlim=c(0,72))
  # OK

  # Check for each segment if it crosses yout, if not, set to zero. Test if all point above or below yout
  # Need to account for tolval!!
  for(iseg in c(1:length(unique(df$seg))))
  {
    if((sum(df$y[df$seg==iseg] > (yout+tolval))   == 0) |
       (sum(df$y[df$seg==iseg] < (yout-tolval))   == 0))   df$seg[df$seg==iseg] <- 0
  }
  df[df$min==1 | df$max==1,]
  df %>% ggplot(aes(x=x,y=y,col=factor(seg)))+geom_line()+geom_point()+
    geom_point(data=df %>% filter(df$max==1),col='black',pch=1,size=3)+
    geom_point(data=df %>% filter(df$min==1),col='black',pch=1,size=3)
  df %>% ggplot(aes(x=x,y=y,col=factor(seg)))+geom_line()+geom_point()+
    geom_point(data=df %>% filter(df$max==1),col='black',pch=1,size=3)+
    geom_point(data=df %>% filter(df$min==1),col='black',pch=1,size=3)+
    coord_cartesian(xlim=c(0,72))
  # OK!

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


#' multi.peaks
#'
#' @param x
#' @param y
#'
#' @return Return multiple x hits of (x,y) series for y-value yout and returns NA if no hit
#' @export
#'
#' @examples
multi.peaks <- function(x,y)
{
  # Checks
  if(length(x) != length(y)) return(NA)
  if(length(x) <= 1)         return(NA)

  findpeaks(y)

  return(xcollect)
}

#' get.x.multi.landmarks
#'
#' @param x
#' @param y
#' @param w
#' @param tol
#'
#' @return
#' @export
#'
#' @examples
get.x.multi.landmarks <- function(x,y,w=17,tol=1e-9) {
  # Using splinefun to determine derivatives
  # tol to be passed on to multi.approx function only
  # Checks
  if(length(x) != length(y)) return(NA)
  if(length(x) <= 1)         return(NA)

  # First data point
  f0.0 <- x[1]
  lm    <- data.frame(x=f0.0,type='start')

  # First derivative=0
  f1.sg <- savitzkyGolay(X = y,  m = 1,  p = 1,  w = w)
  x1.sg  <- x[ (1+(w-1)/2) : (length(x)-(w-1)/2) ]

  f1.0     <- multi.approx(x1.sg,f1.sg,yout=0,tol=tol)  # NA if no maximum
  if(sum(!is.na(f1.0))>=1)
  {
    add    <- data.frame(x=f1.0,type='max')
    lm     <- rbind(lm,add)
  }

  f2.sg <- savitzkyGolay(X = y,  m = 2,  p = 2,  w = w)
  x2.sg  <- x[ (1+(w-1)/2) : (length(x)-(w-1)/2) ]

  f2.0 <- multi.approx(x2.sg,f2.sg,yout=0,tol=tol)   # NA if no inflection point
  if(sum(!is.na(f2.0))>=1)
  {
    add  <- data.frame(x=f2.0,type='inflec')
    lm   <- rbind(lm,add)
  }

  # Last data point - may be "asymptote"
  f9.9 <- x[length(x)]
  add  <- data.frame(x=f9.9,type='end')
  lm   <- rbind(lm,add)

  # arrange by increasing x
  lm <- lm %>% arrange(x)

  return(lm)
}

# x    <- ref$x
# y    <- ref$y
# lm   <- my.ref.lm.init
# tol  <- 1e-9

# x <- my.ref$x
# y <- my.ref$y
# lm <- my.ref.lm.init
# tol <- 0.01*tolapply

#' refine.x.multi.landmarks
#'
#' @param x
#' @param y
#' @param lm
#' @param tol
#' @param w1
#' @param w2
#'
#' @return
#' @export
#'
#' @examples
refine.x.multi.landmarks <- function(x,y,lm,tol=1e-9,w1=5,w2=7)
{
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

      # f2
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

#' get.ref.x.open.end
#'
#' Determine last x value for open end reference curve with some defaults
#'
#' @param x
#' @param y
#' @param lm
#' @param step.x.factor
#' @param tol
#'
#' @return
#' @export
#'
#' @examples
get.ref.x.open.end <- function(x,y,lm,step.x.factor=1.5,tol=0.01) {
  # tolerance used, can probably be the same as for get.x.multi.landmarks and refine.x.multi.landmarks
  # Must contain type="end" as last entry"
  if(lm$type[length(lm$type)] != "end") return("No \"end\" in landmarks")

  # landmarks min/max (incl. first, excl. last data point)
  min.lm <- min(lm$y[lm$type!='end'])
  max.lm <- max(lm$y[lm$type!='end'])
  # If one landmark only, take absolute difference from zero
  if (dim(lm[lm$type!='end',])[1] == 0) stop('ERROR in reference landmark dataframe, no \"end\" defined')
  if (dim(lm[lm$type!='end',])[1] == 1) ytol <- tol*(max(y) - min(y))  # Start and end point only available in lm
  if (dim(lm[lm$type!='end',])[1] >  1) ytol <- tol*(max.lm - min.lm)

  # Last segment, all data point from inflection point to end
  xlast      <- x[x > lm$x[length(lm$x)-1]]         # first points from inflection point onwards
  xlastNext  <- step.x.factor*xlast
  ylast      <- approx(x,y,xout=xlast)$y
  ylastNext  <- approx(x,y,xout=xlastNext,rule=1)$y # NA's if outside domain
  notna      <- !is.na(ylastNext)

  xlast      <- xlast[notna]
  xlastNext  <- xlastNext[notna]
  ylast      <- ylast[notna]
  ylastNext  <- ylastNext[notna]

  # Vector of differences
  ydiff      <- abs(ylast - ylastNext)

  ggplot(data.frame(x=xlast,y=ylast),aes(x=x,y=y))+geom_line()
  ggplot(data.frame(x=xlast,y=ydiff),aes(x=x,y=y))+geom_line()+geom_hline(yintercept = ytol)

  # Find xNext position
  xsol <- multi.approx(xlast,ydiff,yout=ytol)   # using default "tolnoise"
  # For emax we have 2 solutions: one at start when curve is slowly increasing, one at end when curve is slowly flattening out
  # Always take open end = last solution
  xhit <- xsol[length(xsol)]
  # if(length(xhit)!=1) return("Multiple solutions from get.ref.x.open.end()")
  if(is.na(xhit)) return(paste0("No solution within tol ",tol," with step.x.factor ",step.x.factor," from get.ref.x.open.end()"))

  yhit <- approx(x,y,xout=xhit)$y

  # Replace "end"
  lm$x[lm$type=='end']  = xhit
  lm$y[lm$type=='end']  = yhit

  return(lm)
}

#' get.query.x.open.end
#'
#' @param ref
#' @param query
#' @param lm.ref
#' @param lm.query
#' @param ngrid
#' @param scaling
#'
#' @return
#' @export
#'
#' @examples
get.query.x.open.end <- function(ref,query,lm.ref,lm.query,ngrid=100,scaling='linear') {
  # Ref end segment represented by ngrid points
  refGrid      <- data.frame(x=seq(from=lm.ref$x[length(lm.ref$x)-1],
                                   to=lm.ref$x[length(lm.ref$x)],length.out=ngrid),
                             step=c(1:ngrid))
  refGrid$y    <- approx(ref$x,ref$y,xout=refGrid$x)$y

  query.x.start <- lm.query$x[length(lm.query$x)-1]
  query.x.last  <- lm.query$x[length(lm.query$x)]

  refGrid %>% ggplot(aes(x=x,y=y))+geom_line()+
    geom_line(data=query,col='red')

  result <- optimize(get.y.diff.ref.query,
                     interval = c(query.x.start,query.x.last),
                     query.x.start=query.x.start,
                     query=query,
                     refGrid=refGrid,
                     ngrid=ngrid,
                     scaling=scaling)

  # Check if at boundary
  if(abs(result$minimum-query.x.start)<1e-3) print("WARNING: Query minimum x last near landmark before end in get.query.x.open.end")
  if(abs(result$minimum-query.x.last)<1e-3)  print("WARNING: Query minimum x last near last x value provided/explored in get.query.x.open.end")
  if(abs(result$minimum-query.x.last)<1e-3)  print("         Increase the tolerance or increase length simulated profile or decrease x-step factor")

  if(abs(result$minimum-query.x.start)<1e-3) stop("WARNING: Query minimum x last near landmark before end in get.query.x.open.end")
  if(abs(result$minimum-query.x.last)<1e-3)  stop("WARNING: Query minimum x last near last x value provided/explored in get.query.x.open.end")
  if(abs(result$minimum-query.x.last)<1e-3)  stop("         Increase the tolerance or increase length simulated profile or decrease x-step factor")

  lm.query$x[lm.query$type=='end'] <- result$minimum
  lm.query$y[lm.query$type=='end'] <- approx(query$x,query$y,xout=result$minimum)$y

  return(lm.query)
}

#' get.y.diff.ref.query
#'
#' With constant or linearly changing scaling factor, with some default values
#'
#' @param x
#' @param query.x.start
#' @param query
#' @param refGrid
#' @param ngrid
#' @param scaling
#'
#' @return
#' @export
#'
#' @examples
get.y.diff.ref.query <- function(x, query.x.start, query, refGrid, ngrid=100, scaling='linear')
{
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
  # Ignore gridpoints when a y-value is equal to zero
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

  # To minimize y-differences with reference
  diff <- sum((refGrid$y-queryGrid$y.scaled)**2)/ngrid

  return(diff)
}


#' map.segment
#'
#' @param ref
#' @param query
#' @param lim.ref
#' @param lim.query
#' @param scale_by
#'
#' @return
#' @export
#'
#' @examples
map.segment <- function(ref,query,lim.ref,lim.query,scale_by=c('y','y')) {
  # First and last data landmark points to be matched and transformed accordingly

  # 1. scale along y using linear or constant scaling
  if(scale_by[1]=='1')  scale.start <- 1
  if(scale_by[1]=='y')  scale.start <- lim.ref$y[1]/lim.query$y[1]
  if(scale_by[1]=='c')  scale.start <- lim.ref$y[2]/lim.query$y[2]
  if(scale_by[2]=='1')  scale.end   <- 1
  if(scale_by[2]=='y')  scale.end   <- lim.ref$y[2]/lim.query$y[2]
  if(scale_by[2]=='c')  scale.end   <- lim.ref$y[1]/lim.query$y[1]

  scale.df <- data.frame(x=c(lim.query$x[1],lim.query$x[2]),
                         y=c(scale.start,scale.end))
  query$y.scaling <- approx(scale.df$x, scale.df$y, xout=query$x)$y
  query$y.scaled  <- query$y.scaling*query$y

  # 2. scale along x using linear matching (if just outside ref$y domian, than max ref$y is used)
  query$x.scaled  <- approx(ref$y, ref$x, xout=query$y.scaled,rule=2)$y
  query <- query %>%
    mutate(x.scaling = ifelse(x==0,0,x.scaled/x))

  return(query)
}


#' explore.sg
#'
#' @param x.orig
#' @param y.orig
#' @param xmin
#' @param xmax
#' @param ymin
#' @param ymax
#' @param w0
#' @param w1
#' @param w2
#' @param p0
#' @param p1
#' @param p2
#' @param tag
#'
#' @return
#' @export
#'
#' @examples
explore.sg <- function(x.orig,y.orig,xmin,xmax,ymin,ymax,
                       w0=3,w1=5,w2=7,p0=0,p1=1,p2=2,tag='tag') {
  # Savitzky-Golay Smoothing
  # NOTE: evenly spaced vector assumed for this method!

  x <- x.orig
  y <- y.orig

  y.sg <- savitzkyGolay(X = y,  m = 0,  p = p0,  w = w0)
  x.sg <- x[ (1+(w0-1)/2) : (length(x)-(w0-1)/2) ]

  # First derivative=0
  f1.sg <- savitzkyGolay(X = y,  m = 1,  p = p1,  w = w1)
  x1.sg <- x[ (1+(w1-1)/2) : (length(x)-(w1-1)/2) ]

  # Second derivative=0
  f2.sg <- savitzkyGolay(X = y,  m = 2,  p = p2,  w = w2)
  x2.sg <- x[ (1+(w2-1)/2) : (length(x)-(w2-1)/2) ]

  df <- rbind(data.frame(x=x.sg,y=y.sg,type='f0'),
              data.frame(x=x1.sg,y=f1.sg,type='f1'),
              data.frame(x=x2.sg,y=f2.sg,type='f2'))

  # Create matrix to investigate further original function
  mysg0 <- NULL
  pp <- c(0,1,2,3)
  ww <- c(3,7,11,15,19)
  m  <- 0
  for(ip in c(1:length(pp)))
  {
    p <- pp[ip]
    for(iw in c(1:length(ww)))
    {
      w    <- ww[iw]
      if(p>=m & w>p)
      {
        f1.sg <- savitzkyGolay(X = y,  m = m,  p = p,  w = w)
        x.sg <- x[ (1+(w-1)/2) : (length(x)-(w-1)/2) ]
        add  <- data.frame(x=x.sg,y=f1.sg,p=p,w=ifelse(w<10,paste0("0",w),paste0(w)))
        mysg0 <- rbind(mysg0, add)
      }
    }
  }

  # Create matrix to investigate further first order derivative
  mysg1 <- NULL
  pp <- c(1,2,3,4)
  ww <- c(3,7,11,15,19)
  m  <- 1
  for(ip in c(1:length(pp)))
  {
    p <- pp[ip]
    for(iw in c(1:length(ww)))
    {
      w    <- ww[iw]
      if(p>=m & w>p)
      {
        f1.sg <- savitzkyGolay(X = y,  m = m,  p = p,  w = w)
        x.sg <- x[ (1+(w-1)/2) : (length(x)-(w-1)/2) ]
        add  <- data.frame(x=x.sg,y=f1.sg,p=p,w=ifelse(w<10,paste0("0",w),paste0(w)))
        mysg1 <- rbind(mysg1, add)
      }
    }
  }

  # Create matrix to investigate second order derivative
  mysg2 <- NULL
  pp <- c(2,3,4,5)
  ww <- c(3,7,11,15,19)
  m  <- 2
  for(ip in c(1:length(pp)))
  {
    p <- pp[ip]
    for(iw in c(1:length(ww)))
    {
      w    <- ww[iw]
      if(p>=m & w>p)
      {
        f2.sg <- savitzkyGolay(X = y,  m = m,  p = p,  w = w)
        x.sg <- x[ (1+(w-1)/2) : (length(x)-(w-1)/2) ]
        add  <- data.frame(x=x.sg,y=f2.sg,p=p,w=ifelse(w<10,paste0("0",w),paste0(w)))
        mysg2 <- rbind(mysg2, add)
      }
    }
  }
}


