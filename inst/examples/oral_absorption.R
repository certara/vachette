# oral-absorption example from mrgsolve

library(mrgsolve)
library(dplyr)
library(vachette)

tag      <- paste0("v",packageVersion("vachette"))
model  <- "oral-absorption"

script       <- paste0("vachette/examples/oral_absorption_",tag,".R")
setup        <- "A"     # choice of ref and query/queries. Choose B is ref/query should be swapped
my.covariate <- "none"  # To be changed for each model example

extra  <- "-delta01"      # extra info to add to output file names
mrgdelta <- 0.1

set.seed(125)

oral1cmt.ipred <- '

[ PARAM ] @annotated
TVCL : 1    : Clearance (L/hr)
TVV  : 35   : Volume of distribution (L)
TVKA : 0.5  : Absorption rate constant (1/hr)
TVF  : 1    : Bioavailability
WT   : 70   : Weight (kg)

[ CMT ] GUT CENT

[ SET ] delta=1

[ MAIN ]

double CL = TVCL*pow(WT/70,0.75)*exp(ECL);
double V  = TVV*exp(EV);
double KA = TVKA*pow(WT/70,0.75)*exp(EKA);

F_GUT  = TVF*pow(WT/70,0.9);

[ OMEGA ] @annotated
ECL : 0.1 : Eta-CL
EV  : 0.1 : Eta-V
EKA : 0.1 : Eta-KA

[ SIGMA ] @labels PROP ADD
0.00 0.0

[ ODE ]
dxdt_GUT = -KA*GUT;
dxdt_CENT = KA*GUT - (CL/V)*CENT;

[ TABLE ]
double IPRED = CENT/V;
double DV    = IPRED*(1+PROP)+ADD;

while(DV < 0) {
simeps();
DV = IPRED*(1+PROP)+ADD;
}

[ CAPTURE ] CP = CENT/V WT DV TVCL ECL;

'

oral1cmt.pred <- '

[ PARAM ] @annotated
TVCL : 1    : Clearance (L/hr)
TVV  : 35   : Volume of distribution (L)
TVKA : 0.5  : Absorption rate constant (1/hr)
TVF  : 1    : Bioavailability
WT   : 70   : Weight (kg)

[ CMT ] GUT CENT

[ SET ] delta=1

[ MAIN ]

double CL = TVCL*pow(WT/70,0.75)*exp(ECL);
double V  = TVV*exp(EV);
double KA = TVKA*pow(WT/70,0.75)*exp(EKA);

F_GUT  = TVF*pow(WT/70,0.9);

[ OMEGA ] @annotated
ECL : 0.0 : Eta-CL
EV  : 0.0 : Eta-V
EKA : 0.0 : Eta-KA

[ SIGMA ] @labels PROP ADD
0.000 0.00

[ ODE ]
dxdt_GUT = -KA*GUT;
dxdt_CENT = KA*GUT - (CL/V)*CENT;

[ TABLE ]
double IPRED = CENT/V;
double DV    = IPRED*(1+PROP)+ADD;

while(DV < 0) {
simeps();
DV = IPRED*(1+PROP)+ADD;
}

[ CAPTURE ] CP = CENT/V WT DV TVCL ECL;

'
my.ipred <- mcode("oral1cmt.ipred", oral1cmt.ipred)
my.pred  <- mcode("oral1cmt.pred", oral1cmt.pred)

# ---------------- SIMULATED TYPICAL/POPULATION CURVES -------------

# (Typical prediction)

# oral 2-dose
# xstop         <- 96   # Actual observations end at xstop = 48 hr
# xlast.user    <- 640  # User provided x value to max. simulate out - replaces expand.factor
# #expand.factor <- 10   # factor to expand simulated curve to find last X values for ref and query
# tolend        <- 0.001  # Max tolerance allowed minimizing difference open end ref and query
# tolnoise      <- 1e-8   # Max tolerance allowed between gridpoints to identify landmarks for stepsize=1
# step.x.factor <- 1.5  # Factor x-steps between which y-values are lower than the SCALED tolerance
# ngrid.open.end<- 100  # Number of grid points to characterize open end curves

xstop         <- 48   # Actual observations end at xstop = 48 hr
xlast.user    <- 500  # User provided x value to max. simulate out - replaces expand.factor
#expand.factor <- 10   # factor to expand simulated curve to find last X values for ref and query
tolend        <- 0.001  # Max tolerance allowed minimizing difference open end ref and query
tolnoise      <- 1e-8   # Max tolerance allowed between gridpoints to identify landmarks for stepsize=1
step.x.factor <- 1.1  # Factor x-steps between which y-values are lower than the SCALED tolerance
ngrid.open.end<- 100  # Number of grid points to characterize open end curves

# Original situation (query contraction to fit ref)
if(setup=="A") cov.ref        <- 70
if(setup=="A") cov.query      <- 30
# Alternative situation (query expansion to fit ref)
if(setup=="B") cov.ref        <- 30
if(setup=="B") cov.query      <- 70

my.covariate <- "weigth"

# # population data
data <- expand.idata(ID=1,WT=c(cov.query,cov.ref))
#
# simulate pred
long.out <- my.pred %>%
  ev(amt=100,addl=0,ii=24,cmt=1) %>%
  idata_set(data) %>%
  mrgsim(delta=mrgdelta,end=step.x.factor*xlast.user)

long  <- as.data.frame(long.out) %>%
  filter(!(time==0 & GUT==0)) %>%     # No Dose rec
  mutate(x=time, y=CP, COV=WT)         # Standardize

# all in standard parameters x,y
longref   <- long %>% filter(WT==cov.ref)
longquery <- long %>% filter(WT==cov.query)

# (Representing observations)

nsample <- 20
#
# individual data
idata <- expand.idata(ID = c(1:50),
                      WT = c(cov.query,cov.ref))

# simulate ipred
set.seed(20220131)
imy.out <- my.ipred %>%
  ev(amt=100,addl=0,ii=24,cmt=1) %>%
  idata_set(idata) %>%
  mrgsim(delta=mrgdelta,end=xstop)

indiv <- as.data.frame(imy.out) %>%
  filter(!(time==0 & GUT==0)) %>%     # No Dose rec
  mutate(x=time, y=CP, COV=WT)         # Standardize

# Sample nsample reference and nsample query data points:
indivsam.query <- indiv %>%
  filter(WT==cov.query) %>%
  slice(sample(n(),nsample,replace=F))
indivsam.ref <- indiv %>%
  filter(WT==cov.ref) %>%
  slice(sample(n(),nsample,replace=F))
indivsam <- rbind(indivsam.query,indivsam.ref) %>%
  arrange(time)


# Start vachette functions


p0 <- longref %>%
  ggplot(aes(x=x,y=y))+
  geom_line(col='red')+
  geom_line(data=longquery,col='blue') +
  labs(title=paste0(model," Extrapolated Query (blue) and Reference (red) ",setup),
       caption=script)


my.ref            <- long[long$COV==cov.ref,]

x <- x.orig <- my.ref$x
y <- y.orig <- my.ref$y

# Using splinefun to determine derivatives
# Checks
if(length(x) != length(y)) stop('Error, length x != lenght y')
if(length(x) <= 1)         stop('Error, no x coordinates')

# How does f0 look like?
data.frame(x=x,y=y) %>% ggplot(aes(x=x,y=y))+geom_line()

# ============== new 2-step landmark finding function =================

# ======== A. Landmarks approximate positions =========================

# Small stepzise -> small tolerance
# Large stepsize -> large tol
tolapply <- tolnoise*(x[2]-x[1])   # Grid stepsize


ref               <- long[long$COV==cov.ref,]
my.ref.lm.init    <- get.x.multi.landmarks(ref$x,ref$y,w=17,tol=tolapply)
my.ref.lm.init$y  <- approx(my.ref$x,my.ref$y, xout=my.ref.lm.init$x)$y

query               <- long[long$COV==cov.query,]
my.query.lm.init    <- get.x.multi.landmarks(query$x,query$y,w=17,tol=tolapply)
my.query.lm.init$y  <- approx(query$x,query$y, xout=my.query.lm.init$x)$y

# ======== B. Landmarks position cleaning and refinement ==============

my.ref.lm.refined    <- refine.x.multi.landmarks(my.ref$x,my.ref$y,lm=my.ref.lm.init,tol=0.01*tolapply)
my.ref.lm.refined$y  <- approx(ref$x, ref$y, xout=my.ref.lm.refined$x)$y

my.query.lm.refined    <- refine.x.multi.landmarks(query$x,query$y,lm=my.query.lm.init,tol=0.01*tolapply)
my.query.lm.refined$y  <- approx(query$x,query$y, xout=my.query.lm.refined$x)$y

# Show approx and refined landmarks
p0.lm.init <- p0 +
  coord_cartesian(xlim=c(0,96))+
  geom_point(data=my.ref.lm.init,aes(x=x,y=y),pch=21,fill='white',size=4)+
  geom_point(data=my.query.lm.init,aes(x=x,y=y),pch=21,fill='white',size=4)+
  geom_point(data=my.ref.lm.refined,aes(x=x,y=y),pch=21,fill='yellow',size=3)+
  geom_point(data=my.query.lm.refined,aes(x=x,y=y),pch=21,fill='yellow',size=3)

# Check if query and ref have same landmarks in same order
if(dim(my.ref.lm.refined)[1] != dim(my.query.lm.refined)[1] |
   !sum(my.ref.lm.refined$type == my.query.lm.refined$type)>0)
  stop("No matching landmarks between reference and query")

# ------------------ LAST X FOR REF ------------------

# Different step.x.factor brings different last x values,
# but probably it does not matter for transformation (to be tested)
# my.ref.lm.new.old <- get.ref.open.end.old(longref$x,longref$y,my.ref.lm,tol=0.01,step.x.factor=1.5)

my.ref.lm     <- get.ref.x.open.end(longref$x,longref$y,my.ref.lm.refined,
                                    step.x.factor=step.x.factor,tol=tolend)

# ------------------ LAST X FOR QUERY ------------------

# my.query.lm   <- get.query.x.open.end(longref,longquery,my.ref.lm,my.query.lm,ngrid=100,scaling='linear')
scaling <- 'linear'
if (model =='sigmoid') scaling <- 1        # If Emin=Emin and Emax=Emax for ref/query
if (model =='pembro')  scaling <- "linear.not.0"  # Ignore y-values of 0 (first time point)

my.query.lm <- get.query.x.open.end(longref,longquery,my.ref.lm,my.query.lm.refined,
                                    ngrid=ngrid.open.end,scaling=scaling)


# Show landmarks
p0.lm <- p0 +
  geom_point(data=my.ref.lm,aes(x=x,y=y),pch=21,fill='white',size=2)+
  geom_point(data=my.query.lm,aes(x=x,y=y),pch=21,fill='white',size=2)
if(model=="sigmoid"  ) plot(p0.lm + scale_x_log10())

# ----------------------- MAPPING SEGMENTS -----------------------

nseg     <- dim(my.ref.lm)[1]-1

# ----- Scaling segments by superimposing landmarks ------
# Needs to be automated, always 'c' for (0,0)
# For the time being, if query y-value is zero, then constant scaling
scale.first <- ifelse(longquery$y[1]==0,'c','y')
scale.use   <- c(scale.first,rep('y',nseg)) # c = constant (using the other landmark point)

# if(model=='pembro') stop("check this scaling")

if (model=='sigmoid') scale.use <- c("1","1","1") # Only for matching min and max response?

# Map segements and determine y- and x-scaling factor for query
my.query.scaled <- NULL
for(iseg in c(1:nseg))
{
  # Segment to map
  segm.ref <- longref %>%
    filter(x>=my.ref.lm$x[iseg] & x<=my.ref.lm$x[iseg+1])
  segm.query <- longquery %>%
    filter(x>=my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1])

  # First data point scaled by derivatives (l'Hopital) --> did not work, use constant scaling
  my.query.add <- map.segment(segm.ref,
                              segm.query,
                              my.ref.lm[c(iseg,(iseg+1)),],
                              my.query.lm[c(iseg:(iseg+1)),],
                              scale_by=c(scale.use[iseg],scale.use[iseg+1]))
  my.query.scaled <- rbind(my.query.scaled,my.query.add)
}

# Plot scaling factors of curves
p.scaling1 <- my.query.scaled %>%
  ggplot(aes(x=x,y=y.scaling))+
  geom_line()+
  theme_bw()
p.scaling2 <- my.query.scaled %>%
  ggplot(aes(x=x,y=y.scaling))+
  geom_line()+
  scale_y_continuous(limits = c(0,NA))+
  theme_bw()

# --------------- VACHETTE TRANSFORMATION OF OBSERVATIONS ------------

head(my.query.scaled)

# All information is available:
p1 <- my.query.scaled %>%
  ggplot(aes(x=x,y=y)) +
  geom_line() +
  geom_line(aes(x=x.scaled,y=y.scaled),col='red')+
  labs(title="Query original (black) and query transformed (red)",
       caption=script)
p1

# Some ypred points from query to ref curve
p2 <- p1 +
  geom_point(data=my.query.scaled[c(50,200,400),])+
  geom_point(data=my.query.scaled[c(50,200,400),],aes(x=x.scaled,y=y.scaled),col='red')+
  coord_cartesian(xlim=c(0,xstop)) +
  labs(title="Query original (black) and query transformed (red) with some data points",
       caption=script)

# Steps
# 1. Assign segments to each query sample data point and to reference curve
# 2. Determine obs to query curve factor (and obtain prop error and/or add error?)
# 3. Determine obs y-scaling
# 4. For each segment determine obs x-scaling

# 1. Assign segments to each query sample data point and to reference curve
indivsam.query$seg <- longref$seg <- longquery$seg <- NA
for(iseg in c(1:nseg))
{
  indivsam.query <- indivsam.query %>% mutate(seg = ifelse(x>=my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1], iseg, seg))
  longref        <- longref        %>% mutate(seg = ifelse(x>=my.ref.lm$x[iseg]   & x<=my.ref.lm$x[iseg+1], iseg, seg))
  longquery      <- longquery      %>% mutate(seg = ifelse(x>=my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1], iseg, seg))
}
#
if(sum(is.na(indivsam.query$seg)>0))                                                 stop("Error, not all samples could be assigned to a segment")
if(sum(is.na(longref$seg[longref$x<=my.ref.lm$x[my.ref.lm$type=='end']])>0))         stop("Error, not all ref data could be assigned to a segment")
if(sum(is.na(longquery$seg[longquery$x<=my.query.lm$x[my.query.lm$type=='end']])>0)) stop("Error, not all query data could be assigned to a segment")

# Check
p3 <- p1 +
  geom_line(data=longquery,aes(x=x,y=y,col=factor(seg)),lwd=1.5) +
  geom_line(data=longref,aes(x=x,y=y,col=factor(seg)),lwd=1.5) +
  coord_cartesian(xlim=c(0,xstop)) +
  labs(title="Segments",
       caption=script) +
  guides(color=guide_legend(title="Segment"))


# 2. Determine obs to query curve factor (and obtain prop/add error?)
indivsam.query$ypred      <- approx(longquery$x,longquery$y,xout=indivsam.query$x,rule=2)$y
# Do not understand warning message:
#   In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) :
#   collapsing to unique 'x' values

# I think we can use this as an alternative way to account for prop and additive error: to be investigated
indivsam.query$prop.error <- (indivsam.query$y - indivsam.query$ypred)/indivsam.query$ypred
indivsam.query$add.error  <-  indivsam.query$y - indivsam.query$ypred

# Check
p4 <- p1 + geom_point(data=indivsam.query,aes(x=x,y=ypred,col=factor(seg))) +
  coord_cartesian(xlim=c(0,xstop)) +
  labs(title="Preds of samples",
       caption=script) +
  guides(color=guide_legend(title="Segment"))


# 3. Determine obs y-scaling
indivsam.query$y.scaled     <- indivsam.query$y*approx(my.query.scaled$x,my.query.scaled$y.scaling,xout=indivsam.query$x,rule=2)$y
indivsam.query$ypred.scaled <- indivsam.query$ypred*approx(my.query.scaled$x,my.query.scaled$y.scaling,xout=indivsam.query$x,rule=2)$y

# Check
p5 <- p1 + geom_point(data=indivsam.query,aes(x=x,y=ypred.scaled,col=factor(seg))) +
  geom_line(data=longref,aes(x=x,y=y,col=factor(seg)),lwd=1) +
  coord_cartesian(xlim=c(0,xstop)) +
  labs(title="Y-scaled Preds of samples",
       caption=script) +
  guides(color=guide_legend(title="Segment"))



indivsam.query %>% group_by(seg) %>% dplyr::summarise(y.scaled.max=max(y.scaled),y.scaled.min=min(y.scaled))
indivsam.query %>% group_by(seg) %>% dplyr::summarise(ypred.max=max(ypred),y.scaled.min=min(ypred))
longref       %>% group_by(seg) %>% dplyr::summarise(y.max=max(y),y.min=min(y))

# 4. Determine obs x-scaling segment by segment
indivsam.query$xpred.scaled <- NA
for(iseg in c(1:nseg))
{
  indivsam.query$xpred.scaled[indivsam.query$seg==iseg] <-
    approx(longref$y[!is.na(longref$seg) & longref$seg==iseg],
           longref$x[!is.na(longref$seg) & longref$seg==iseg],
           xout=indivsam.query$ypred.scaled[indivsam.query$seg==iseg],
           rule=2)$y

  indivsam.query$x.scaled[indivsam.query$seg==iseg] <- indivsam.query$xpred.scaled[indivsam.query$seg==iseg]
}

# ------------ Done ------------------

# Plot up to last original or extended query data point
max.x <- max(indivsam.query$x.scaled, indivsam.query$x, xstop)

# Check
p6 <- p1 +
  geom_line(data=longref,aes(x=x,y=y),lwd=1.5) +
  geom_point(data=indivsam.query,aes(x=xpred.scaled,y=ypred.scaled,col=factor(seg))) +
  coord_cartesian(xlim=c(0,max.x)) +
  labs(title="Y-scaled and X-scaled Preds of samples",
       caption=script) +
  guides(color=guide_legend(title="Segment"))


# Check
p7 <- p1 +
  geom_line(data=longref,aes(x=x,y=y),lwd=1) +
  geom_point(data=indivsam.query,aes(x=x,y=y,col=factor(seg))) +
  geom_point(data=indivsam.query,aes(x=x.scaled,y=y.scaled,col=factor(seg))) +
  geom_segment(data=indivsam.query,aes(x=x,y=y,xend=x.scaled,yend=y.scaled),
               arrow = arrow(length = unit(0.2, "cm")),col='darkgrey') +
  geom_vline(xintercept = xstop, col='grey', lty=2) +
  coord_cartesian(xlim=c(0,max.x)) +
  labs(title="Vachette transformation normal scale",
       caption=script) +
  guides(color=guide_legend(title="Segment"))
# Looks good

lowest.y <- 0.1
if(model=='sigmoid') lowest.y <- min(long$y)
if(model=='pembro')  lowest.y <- 10
p8 <- p1 +
  geom_line(data=longref[longref$y>=lowest.y,],aes(x=x,y=y),lwd=1) +
  geom_point(data=indivsam.query[indivsam.query$y>=lowest.y,],aes(x=x,y=y,col=factor(seg))) +
  geom_point(data=indivsam.query[indivsam.query$y>=lowest.y,],aes(x=x.scaled,y=y.scaled,col=factor(seg))) +
  geom_segment(data=indivsam.query[indivsam.query$y>=lowest.y,],aes(x=x,y=y,xend=x.scaled,yend=y.scaled),
               arrow = arrow(length = unit(0.2, "cm")),col='darkgrey') +
  geom_vline(xintercept = xstop, col='grey', lty=2) +
  coord_cartesian(xlim=c(0,max.x),ylim=c(lowest.y,NA)) +
  scale_y_log10() +
  labs(title="Vachette transformation semi-log scale",
       caption=script) +
  guides(color=guide_legend(title="Segment"))

# Or (better)
p8 <- ggplot() +
  geom_line(data=longref[longref$y>=lowest.y,],aes(x=x,y=y),lwd=1) +
  geom_line(data=longquery[longref$y>=lowest.y,],aes(x=x,y=y)) +
  geom_point(data=indivsam.query[indivsam.query$y>=lowest.y,],aes(x=x,y=y,col=factor(seg))) +
  geom_point(data=indivsam.query[indivsam.query$y>=lowest.y,],aes(x=x.scaled,y=y.scaled,col=factor(seg))) +
  geom_segment(data=indivsam.query[indivsam.query$y>=lowest.y,],aes(x=x,y=y,xend=x.scaled,yend=y.scaled),
               arrow = arrow(length = unit(0.2, "cm")),col='darkgrey') +
  geom_vline(xintercept = xstop, col='grey', lty=2) +
  coord_cartesian(xlim=c(0,max.x),ylim=c(lowest.y,NA)) +
  scale_y_log10() +
  labs(title="Vachette transformation semi-log scale",
       caption=script) +
  guides(color=guide_legend(title="Segment"))

# p8 - demostrate distances on log scale have been preserved
indivsam.query$y.query.curve.original.x <- approx(longquery$x, longquery$y, xout=indivsam.query$x)$y
indivsam.query$y.ref.curve.scaled.x     <- approx(longref$x, longref$y, xout=indivsam.query$x.scaled)$y
# Log pred curve points
indivsam.query$log.y.query.curve.original.x <- log(indivsam.query$y.query.curve.original.x)
indivsam.query$log.y.ref.curve.scaled.x     <- log(indivsam.query$y.ref.curve.scaled.x)
# Log original y
indivsam.query$log.y        <- log(indivsam.query$y)
# Log original y.scaled
indivsam.query$log.y.scaled <- log(indivsam.query$y.scaled)
# Differences:
indivsam.query$log.y.diff        <- indivsam.query$log.y        - indivsam.query$log.y.query.curve.original.x
indivsam.query$log.y.diff.scaled <- indivsam.query$log.y.scaled - indivsam.query$log.y.ref.curve.scaled.x

p8.log.y.diff <- indivsam.query %>%
  ggplot()+
  geom_abline(slope=1,col='red')+
  geom_point(aes(x=log.y.diff,y=log.y.diff.scaled))+
  labs(x="Log diff original data to typical query",
       y="Log diff Vachette transformed data to typical reference",
       title="Log-scale difference original data versus typical query curve to\nlog-scale Vachette transformed data to typical reference curve",
       caption=script)
# Perfect

# Vachette plot, both transformed query observations together with reference observations
p9 <- longref %>%
  ggplot() +
  geom_line(data=longref,aes(x=x,y=y),lwd=1) +
  geom_point(data=indivsam.ref,aes(x=x,y=y),col='red') +
  geom_point(data=indivsam.query,aes(x=x.scaled,y=y.scaled),col='purple') +
  coord_cartesian(xlim=c(0,max.x)) +
  labs(title="Final Vachette plot on normal scale",
       caption=script)

# Vachette plot, both transformed query observations together with reference observations
p10 <- longref %>%
  ggplot() +
  geom_line(data=longref,aes(x=x,y=y),lwd=1) +
  geom_point(data=indivsam.ref,aes(x=x,y=y),col='red') +
  geom_point(data=indivsam.query,aes(x=x.scaled,y=y.scaled),col='purple') +
  coord_cartesian(xlim=c(0,max.x),ylim=c(lowest.y,NA)) +
  scale_y_log10() +
  labs(title="Final Vachette plot on semi-log scale",
       caption=script)



# ------------ Saving plots --------------
# Save plots:
plot_dir <- "./plots"
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

pdf(paste0("./plots/vachette-main-",model,"-",my.covariate,"-",tag,extra,"-",setup,".pdf"))

plot(p0+theme_vachette())
plot(p0.lm.init+theme_vachette())
plot(p0.lm+theme_vachette())
plot(p.scaling1+theme_vachette())
plot(p.scaling2+theme_vachette())
plot(p1+theme_vachette())
plot(p2+theme_vachette())
plot(p3+theme_vachette())
plot(p4+theme_vachette())
plot(p5+theme_vachette())
plot(p6+theme_vachette())
plot(p7+theme_vachette())
plot(p8+theme_vachette())
plot(p8.log.y.diff+theme_vachette())
plot(p9+theme_vachette())
plot(p10+theme_vachette())

if(model=='sigmoid')
{
  plot(p0+theme_vachette()+scale_x_log10())
  plot(p0+theme_vachette()+coord_cartesian(xlim=c(0,100)))
  plot(p7+theme_vachette()+coord_cartesian(xlim=c(0.0001,NA)) + scale_x_log10())
  plot(p8+theme_vachette()+coord_cartesian(xlim=c(0.0001,NA)) + scale_x_log10())
  plot(p9+theme_vachette()+coord_cartesian(xlim=c(0.0001,NA)) + scale_x_log10())
  plot(p10+theme_vachette()+coord_cartesian(xlim=c(0.0001,NA)) + scale_x_log10())
}
if(model=='indirect-response')
{
  plot(p0.lm+theme_vachette()+annotate("text",x=40,y=1,label=paste0("Delta=",mrgdelta))+coord_cartesian(xlim=c(0,48),ylim=c(0,11)))
  plot(p7+theme_vachette()+annotate("text",x=40,y=1,label=paste0("Delta=",mrgdelta))+guides(color="none"))
}
if(model=='pembro')
{
  plot(p0+scale_y_log10()+theme_vachette())
  plot(p0.lm.init+scale_y_log10()+theme_vachette())
  plot(p0.lm+scale_y_log10()+theme_vachette())
  plot(p1+scale_y_log10()+theme_vachette())
  plot(p2+scale_y_log10()+theme_vachette())
  plot(p3+scale_y_log10()+theme_vachette())
  plot(p4+scale_y_log10()+theme_vachette())
  plot(p5+scale_y_log10()+theme_vachette())
  plot(p6+scale_y_log10()+theme_vachette())
  plot(p8+coord_cartesian(xlim=c(0,xstop),ylim=c(30,NA))+theme_vachette())
  plot(p10+coord_cartesian(xlim=c(0,xstop),ylim=c(30,NA))+theme_vachette())
}
if(model=='sigmoid') # Add V2acher scaling
{
  plot(p7+theme_vachette()+coord_cartesian(xlim=c(0.0001,NA)) + scale_x_log10() +
         geom_point(data=indivsam[indivsam$COV==cov.query,],aes(x=x.tr,y=y), pch=1, size=4))
  plot(p8+theme_vachette()+coord_cartesian(xlim=c(0.0001,NA)) + scale_x_log10() +
         geom_point(data=indivsam[indivsam$COV==cov.query,],aes(x=x.tr,y=y), pch=1, size=4))
  plot(p9+theme_vachette()+coord_cartesian(xlim=c(0.0001,NA)) + scale_x_log10() +
         geom_point(data=indivsam[indivsam$COV==cov.query,],aes(x=x.tr,y=y), pch=1, size=4))
  plot(p10+theme_vachette()+coord_cartesian(xlim=c(0.0001,NA)) + scale_x_log10() +
         geom_point(data=indivsam[indivsam$COV==cov.query,],aes(x=x.tr,y=y), pch=1, size=4))
}

dev.off()
