# -------------------------------------------------------------------------
#  Sponsor           : Merck
#  Project Number    : MRK-FTE-VRSV-640
#  Project Root Path : Local
# -------------------------------------------------------------------------
#  Program : vachette-main-v7.R
#  Author  : Jos Lommerse - Certara
#  Date    : 02 February 2022
#  Purpose : Vachette method for Visualization of Analyses with Covariates
# -------------------------------------------------------------------------
#  Software : R version 4.1.2 (2021-11-01)
#  Platform : ThinkPad P15s Gen 2 Model 20W6008AUS
#  Environment : Windows-10, 11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz
# -------------------------------------------------------------------------

library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(dtw)
library(rootSolve)
library(pracma)
library(gtools)
library(Rfast)
library(prospectr)

# v1: Based on vachette-functions-oral-absorption-v2.R
# v2: indirect response example with updated functions
# v3: same as v2 (align version with vachette-functions-v3.R)
# v4: change choice of xlast ref
#     add i.v. example (no landmark)
# v5: New folder structure
#     +Simple Emax model
# v6 Using splines for landmark finding
# v7 Using Savitzky-Golay for landmark finding
# v8 Two-dose
# v9 Inspection 2nd order derivatives for inflection points
# v10/v11 2-step landmark finding (approx followed by refined)
# v12 a. Tolerance noise dependent on grid size "tolnoise"
#     b. Tolerance last data point (how far off asymptote) by user "tolend"
#        Really pass tolerances to functions
# v13 Test several models
# v14 Repair single dose "reverse"
#     Repair Emax
#     Repair indirect response
# v15 Pembro model
# v16 test indirect response with step mrgdelta 0.1 but smaller tolerance
#     test ALB as covariate for pembro
# v17 cleaner versions
# v17 pembro-only Pembro only
#     continuous covariate
#     incl. RUV and IIV correction
# v18 Updated RUV (email Nele, )

# v20 Additive error

# Clear memory
rm(list=ls())

tag      <- "v20-pembro-only"
script   <- "vachette-main-v20-pembro-only.R"
# model  <- "sigmoid"
# model  <- "iv"
# model  <- "oral-absorption"
# model  <- "indirect-response"
# model  <- "oral-two-dose"
model <- "pembro"

script       <- paste0("vachette-main-",tag,".R")
IIV_CORR     <- T         # Apply IIV corection
# Either additive transformation properties or proportional
PROP_TR  <- T
ADD_TR   <- ifelse(PROP_TR,F,T)
if(PROP_TR) extra2   <- "-prop"      # extra info to add to output file names
if(ADD_TR)  extra2   <- "-addi"      # extra info to add to output file names

setup        <- "alb-ref-highest"  # extra label
my.covariate <- "none"    # To be changed for each model example

if(IIV_CORR) extra   <- "-iiv-corr"      # extra info to add to output file names
if(!IIV_CORR) extra  <- "-no-iiv-corr"      # extra info to add to output file names
mrgdelta <- 0.1

# Location of this script
home   <- dirname(rstudioapi::getActiveDocumentContext()$path)
getwd()
setwd(home)

source('vachette-functions-v20.R')
source('models-v20.R')

set.seed(121)

# For plots;
render <- theme_bw()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        plot.title = element_text(size=14),
        plot.subtitle=element_text(size=12))


##########################################################
#                                                        #
#                    PEMBRO                              #
#              DEFINE REFERENCE AND QUERY                #
#                                                        #
##########################################################

# ----- Step 1: Model/simulations -----

# (Typical prediction)

xstop         <- 21     # Actual observations end at xstop = 21 days
xlast.user    <- 700    # User provided x value to max. simulate out - replaces expand.factor
#expand.factor <- 10    # factor to expand simulated curve to find last X values for ref and query
tolend        <- 0.001  # Max tolerance allowed minimizing difference open end ref and query
tolnoise      <- 1e-8   # Max tolerance allowed between gridpoints to identify landmarks for stepsize=1
step.x.factor <- 1.5    # Factor x-steps between which y-values are lower than the SCALED tolerance
ngrid.open.end<- 100    # Number of grid points to characterize open end curves


# Select pembro covariate to assess:
# my.covariate <- "sex"
my.covariate <- "albumin"

# Original situation (query contraction to fit ref)
# SEX covariate
if(setup=="A") cov.ref        <- 1
if(setup=="A") cov.query      <- 2
# Alternative situation (query expansion to fit ref)
if(setup=="B") cov.ref        <- 2
if(setup=="B") cov.query      <- 1

# ------

# mod.read <- mread("../literature-examples/pembrolizumab/mrgsolve/pembro")
# mod <- mcode("pembro", pembro)
mod <- mcode("pembro", pembro)
# mod <- mcode("pembro", pembro.ruv.add)

#simulation settings
dose <- 10                 #set pembrolizumab dose in mg/kg
dose.interval <- c(14,21)  #set dosing interval in days
dose.n <- 12
inf.duration = 0.5/24
maxtime = 25*7             #set end of simulation time interval

# ---- JL  start with single dose
dose.n        <- 1
# Single dose
dose.interval <- c(999)  # e.g. c(14) for Q2W

# Full (long) simulation
samples <- merge(c(0:(1*xlast.user*24))/24, 0) %>%    # 1 doses, 14 days
  mutate(times = x + y)
samples2 <- merge(c(0:(1*21*24))/24, 0) %>%   # 1 doses, 21 days
  mutate(times = x + y)

# "Actual" data to be sampled
samples.iiv <- merge(c(0:(1*xstop*24))/24, 0) %>%    # 1 doses, 14 days
  mutate(times = x + y)
samples2.iiv <- merge(c(0:(1*xstop*24))/24, 0) %>%   # 1 doses, 21 days
  mutate(times = x + y)
# ------

#also update omega matrix to zero to simulate typical curves:
mod.typ     <- omat(mod, dmat(0,0))

# Since it is difficult (not possible?) to modify dosedat and ind.dat dataframes (mrgsolve format)
# just do a loop for the simulations
n.alb <- 6   # For typical individual
n.sim <- 3   # For iiv - 5 individuals per covariate

#typical covariate values
tSEX = 1
tIPIP = 0
tALB = 39.6
tEGFR = 88.47
tBSLD = 89.60
tBECOG = 1
tADIAGN = 1

#individual covariate distributions (to match Table 2 in publication)
ind.sex <- sample(c(1,2), size = n.alb, replace = TRUE, prob = c(0.591, 0.409))
ind.diag <- sample(c(1,2), size = n.alb, replace = TRUE, prob = c(0.747, 0.253))
ind.ecog <- sample(c(0,1), size = n.alb, replace = TRUE, prob = c(0.576, 0.424))
ind.ipi <- sample(c(0,1), size = n.alb, replace = TRUE, prob = c(0.655, 0.345))
ind.egfr <- exp(rnorm(n.alb, mean = log(88.7), sd = 0.4))
ind.alb <- exp(rnorm(n.alb, mean = log(30), sd = 0.2))
ind.bsld <- exp(rnorm(n.alb, mean = log(86), sd = 0.5))

# create simulation dataset for n.alb subjects
dosedat <- ev(ID = 1:n.alb, amt = dose, ii = dose.interval, addl = 49) %>%
  mutate(WT = rnorm(n.alb, 76.8, 8),
         amt = dose * WT,
         rate = amt/inf.duration,
         dose = dose)

# ALB variations only - other covariates using typical values
ind.dat <- dosedat %>%
  dplyr::select(ID, WT) %>%
  mutate(SEX = tSEX,
         ADIAGN = tADIAGN,
         BECOG = tBECOG,
         IPIP = tIPIP,
         EGFR = tEGFR,
         ALB = ind.alb,
         BSLD = tBSLD)

# times <- c(0:(14*24))

output.typ      <- NULL
output.iiv.orig <- NULL
output.iiv      <- NULL

# times<-c(1:5)
times <- c(0:(10*(14*24))/10)


# Simulate w/o IIV (mod.typ)
sim <- mod.typ %>%
  data_set(dosedat) %>%
  idata_set(ind.dat) %>%
  carry_out(EVID) %>%
  Req(CONC, CONC_RUV, SEX, IPIP, ALB, EGFR, BSLD, BECOG, ADIAGN) %>%
  mrgsim(tad = TRUE, tgrid = times)

output.typ <- as.data.frame(sim) %>%
  filter(EVID == 0, time!=0) %>%   # time=0 causes problems at scaling
  dplyr::select(-EVID) %>%
  mutate(sim.n = 0)%>%
  mutate(ID.orig = ID)

for(isim in c(1:n.sim))
{
  # Simulate with IIV (mod)
  sim <- mod %>%
    data_set(dosedat) %>%
    idata_set(ind.dat) %>%
    carry_out(EVID) %>%
    Req(CONC, CONC_RUV, SEX, IPIP, ALB, EGFR, BSLD, BECOG, ADIAGN) %>%
    mrgsim(tad = TRUE, tgrid = times)

  output <- as.data.frame(sim) %>%
    filter(EVID == 0, time!=0) %>%   # time=0 causes problems at scaling
    dplyr::select(-EVID) %>%
    mutate(ID.orig = ID) %>%
    mutate(sim.n = isim) %>%
    mutate(ID = ID + (isim-1)*length(unique(sim$ID)))

  output.iiv <- rbind(output.iiv,output)
}


# Factorize ALB covariate and keep covariate number
output.typ$ALB.f <- factor(output.typ$ALB,levels=sort(unique(output.typ$ALB)))
output.typ$ALB.n <- as.numeric(output.typ$ALB.f)
output.iiv$ALB.f <- factor(output.iiv$ALB,levels=sort(unique(output.iiv$ALB)))
output.iiv$ALB.n <- as.numeric(output.iiv$ALB.f)
# Numbers must align! - .... build-in a check ....

# Check same IDs typical and iiv
unique(output.typ[,c('ID','ID.orig','ALB','ALB.n')]) %>% slice(1:10)
unique(output.iiv[,c('ID','ID.orig','ALB','ALB.n')]) %>% slice(1:10)

# Standardize naming variables PRED/IPRED (w/o RUV) and OBS (IPRED+RUV)
output.typ <- output.typ %>%
  rename(PRED = CONC)
output.iiv <- output.iiv %>%
  rename(IPRED = CONC) %>%
  rename(OBS   = CONC_RUV) %>%
  left_join(output.typ[,c('ID.orig','time','PRED')],by=c('ID.orig','time'))


# Plot
output.plot <- output.iiv %>%
  filter(ID %in% sample(length(unique(output.iiv$ID)),12)) %>% # Sample IDs
  filter(time %in% sample(c(1:100),12))                        # sample time points

output.plot  %>%
  ggplot(aes(x=time,group=ID,col=factor(ID)))+
  geom_point(aes(y=PRED))+
  geom_point(data=output.plot, aes(y=IPRED))+
  geom_point(data=output.plot, aes(y=OBS),pch=1)+
  geom_line(aes(y=PRED),lwd=0.8,col='black')+
  geom_line(data=output.plot, aes(y=IPRED))+
  geom_line(data=output.plot, aes(y=OBS),lty=2)+
  facet_wrap(~paste("Albumin",round(ALB,1)),ncol=3)+
  guides(color=guide_legend(title="ID"))+
  render

# ------- Select reference covariate effect and get samples/observations -----

# Reference ALB, somewhere in the middle
unique(output.typ$ALB.n)
# ALB.ref    <- round(length(unique(output.typ$ALB.n))/2,0)
# Lowest
# ALB.ref    <- 1
# Highest
ALB.ref    <- n.alb

# Samples query (take times of first subject, up to xstop)
nsample      <- 5
sample.times <- sample(output.iiv$time[output.iiv$time<=xstop & output.iiv$ID==unique(output.iiv$ID)[1]],nsample,replace=F)

# Currently 20 IDs per ALB (10 ALB values)
indivsam.query <- output.iiv %>%
  filter(ALB.n != ALB.ref) %>%         # NOT ref ALB
  filter(time %in% sample.times) %>%
  group_by(ID,ALB) %>%
  slice(sample(1:10,5)) %>%            # 3 samples per ID
  mutate(x=time, y=OBS, COV=ALB.n) %>% # Standardize
  arrange(ID,time)                     # Sort ID/time
# 180 ID's with each 3 data point

# Samples for reference
indivsam.ref <- output.iiv %>%
  filter(ALB.n == ALB.ref) %>%         # NOT ref ALB
  filter(time %in% sample.times) %>%
  group_by(ID,ALB) %>%
  slice(sample(1:10,5)) %>%            # 3 samples per ID
  mutate(x=time, y=OBS, COV=ALB.n) %>% # Standardize
  arrange(ID,time)                     # Sort ID/time
# 20 ID's with each 3 data point

# ----- Step 2: Correct for IIV (or omit) -----

# (x,y) coordinates are used for the Vachette transformation

# df with simulated curves related to all covariates
# df with observations - associated covariates
# User points to a reference covariate effect

# Vachette --> create ref curve and query curve


# iiv correction for both for query and reference
if(IIV_CORR) indivsam.query$y <- indivsam.query$OBS - indivsam.query$IPRED + indivsam.query$PRED
if(IIV_CORR) indivsam.ref$y   <- indivsam.ref$OBS   - indivsam.ref$IPRED   + indivsam.ref$PRED

# or omit iiv correction:
if(!IIV_CORR) indivsam.query$y <- indivsam.query$OBS
if(!IIV_CORR) indivsam.ref$y   <- indivsam.ref$OBS

# Plot - sample 12 indivs
output.plot <- indivsam.query %>%
  filter(ID %in% sample(length(unique(indivsam.query$ID)),12)) %>%
  arrange(time)
output.plot  %>%
  ggplot(aes(x=time,col=factor(ID)))+
  geom_line(aes(y=PRED),lwd=0.8)+
  geom_line(data=output.plot, aes(y=IPRED))+
  geom_line(data=output.plot, aes(y=OBS),lty=2)+
  # geom_point(aes(y=PRED))+
  # geom_point(data=output.plot, aes(y=IPRED))+
  geom_point(data=output.plot, aes(y=OBS),pch=1)+
  geom_point(data=output.plot, aes(y=y),pch=19)+
  # geom_segment(x=time,y=OBS,xend=time,yend=OBS.iiv)+
  facet_wrap(~paste(ID," Alb:",round(ALB,1)),ncol=4)+
  guides(color=guide_legend(title="ID"))+
  render
# OK
output.plot  %>%
  ggplot(aes(x=time,col=factor(ID)))+
  geom_line(aes(y=PRED),lwd=0.8,col='black')+
  geom_line(data=output.plot, aes(y=IPRED))+
  geom_line(data=output.plot, aes(y=OBS),lty=2)+
  # geom_point(aes(y=PRED))+
  # geom_point(data=output.plot, aes(y=IPRED))+
  geom_point(data=output.plot, aes(y=OBS),pch=1)+
  geom_point(data=output.plot, aes(y=y),pch=19)+
  # geom_segment(x=time,y=OBS,xend=time,yend=OBS.iiv)+
  facet_wrap(~paste("Albumin",round(ALB,1)),ncol=3)+
  guides(color=guide_legend(title="ID"))+
  render

output.plot <- indivsam.ref %>%
  # filter(ID %in% sample(length(unique(indivsam.ref$ID)),12)) %>%
  filter(ID ==9) %>%
  arrange(time)
output.plot  %>%
  ggplot(aes(x=time,col=factor(ID)))+
  geom_line(aes(y=PRED),lwd=0.8,col='black')+
  geom_line(aes(y=IPRED),lty=2,col='black')+
  # geom_line(data=output.plot, aes(y=IPRED))+
  geom_line(data=output.plot, aes(y=y),lty=1)+
  geom_line(data=output.plot, aes(y=OBS),lty=2)+
  # geom_point(aes(y=PRED))+
  # geom_point(data=output.plot, aes(y=IPRED))+
  geom_point(data=output.plot, aes(y=OBS),pch=1)+
  geom_point(data=output.plot, aes(y=y),pch=19)+
  # geom_segment(x=time,y=OBS,xend=time,yend=OBS.iiv)+
  facet_wrap(~paste("Albumin",round(ALB,1)),ncol=3)+
  guides(color=guide_legend(title="ID"))+
  render
# Correct, OBS is moved from distributed around the IPRED to distributed around PRED

# ----- Steps 3+4: Vachette -----

# ---- Vachette ---

# Collect output for plotting
query.all          <- NULL
ref.lm.all         <- NULL
query.lm.all       <- NULL
query.scaled.all   <- NULL
indivsam.query.all <- NULL

# Typical curves for ALB
long  <- as.data.frame(output.typ) %>%
  mutate(x=time, y=PRED, COV=ALB.n)       # Standardize

# Selection of covariates
ALB.query  <- sort(unique(long$ALB.n))
# Remove ref
ALB.query  <- ALB.query[ALB.query != ALB.ref]

# Collect output and transformations
landmarks            <- NULL
indivsam.query.saved <- indivsam.query

# Loop for each query covariate
iALB <- 1
# for(iALB in c(1:length(ALB.query)))
# {
# Transformation base on typical curves

# Covariate values
cov.ref   <- ALB.ref
cov.query <- ALB.query[iALB]

# Select reference and query samples
indivsam.query <- indivsam.query.saved %>% filter(ALB.n == cov.query)

# Select typical reference and query curves
ref   <- long %>% filter(ALB.n == cov.ref)
query <- long %>% filter(ALB.n == cov.query)

# Quick plot
ref %>%
  ggplot(aes(x=time,y=PRED,col=factor(round(ALB,1)))) +
  geom_line()+
  geom_point(data=indivsam.ref)+  # PRED points
  geom_line(data=query)+
  geom_point(data=indivsam.query) # PRED points

# ======== A. Landmarks approximate positions =========================

# Small stepsize -> small tolerance
# Large stepsize -> large tol
tolapply <- tolnoise*(long$x[2]-long$x[1])   # Grid stepsize

ref               <- long[round(long$COV,1)==cov.ref,]
my.ref.lm.init    <- get.x.multi.landmarks(ref$x,ref$y,w=17,tol=tolapply)
my.ref.lm.init$y  <- approx(ref$x,ref$y, xout=my.ref.lm.init$x)$y

query               <- long[round(long$COV,1)==cov.query,]
my.query.lm.init    <- get.x.multi.landmarks(query$x,query$y,w=17,tol=tolapply)
my.query.lm.init$y  <- approx(query$x,query$y, xout=my.query.lm.init$x)$y

# Check by plot
ref %>% ggplot(aes(x=x,y=y))+geom_line(col='red')+geom_line(data=query,col='blue')

# ======== B. Landmarks position cleaning and refinement ==============

my.ref.lm.refined    <- refine.x.multi.landmarks(x=ref$x,y=ref$y,lm=my.ref.lm.init,tol=0.01*tolapply)
my.ref.lm.refined$y  <- approx(ref$x, ref$y, xout=my.ref.lm.refined$x)$y
# OK
my.query.lm.refined    <- refine.x.multi.landmarks(query$x,query$y,lm=my.query.lm.init,tol=0.01*tolapply)
my.query.lm.refined$y  <- approx(query$x,query$y, xout=my.query.lm.refined$x)$y


# Check if query and ref have same landmarks in same order
if(dim(my.ref.lm.refined)[1] != dim(my.query.lm.refined)[1] |
   !sum(my.ref.lm.refined$type == my.query.lm.refined$type)>0)
  stop("No matching landmarks between reference and query")

# ------------------ LAST X FOR REF ------------------

# Different step.x.factor brings different last x values,
# but probably it does not matter for transformation (to be tested)
# my.ref.lm.new.old <- get.ref.open.end.old(ref$x,ref$y,my.ref.lm,tol=0.01,step.x.factor=1.5)

my.ref.lm     <- get.ref.x.open.end(ref$x,ref$y,my.ref.lm.refined,
                                    step.x.factor=step.x.factor,tol=tolend)

if (iALB == 1)  landmarks <- my.ref.lm %>% mutate(curve='reference', cov = cov.ref)

# ------------------ LAST X FOR QUERY ------------------

# my.query.lm   <- get.query.x.open.end(ref,query,my.ref.lm,my.query.lm,ngrid=100,scaling='linear')
scaling <- 'linear'
if (model =='sigmoid') scaling <- 1        # If Emin=Emin and Emax=Emax for ref/query
if (model =='pembro')  scaling <- "linear.not.0"  # Ignore y-values of 0 (first time point)

my.query.lm <- get.query.x.open.end(ref,query,my.ref.lm,my.query.lm.refined,
                                    ngrid=ngrid.open.end,scaling=scaling)

add.query <- my.query.lm %>%
  mutate(curve='query',
         cov = cov.query)

landmarks <- rbind(landmarks,add.query)
print(landmarks)

# ----------------------- MAPPING SEGMENTS -----------------------

nseg     <- dim(my.ref.lm)[1]-1

# ----- Scaling segments by superimposing landmarks ------
# Needs to be automated, always 'c' for (0,0)
# For the time being, if query y-value is zero, then constant scaling
scale.first <- ifelse(query$y[1]==0,'c','y')
scale.use   <- c(scale.first,rep('y',nseg)) # c = constant (using the other landmark point)

if (model=='sigmoid') scale.use <- c("1","1","1") # Only for matching min and max response?

# Map segments and determine y- and x-scaling factor for query
my.query.scaled <- NULL
for(iseg in c(1:nseg))
{
  # Segment to map
  segm.ref <- ref %>%
    filter(x>=my.ref.lm$x[iseg] & x<=my.ref.lm$x[iseg+1]) %>%
    mutate(seg=iseg)
  segm.query <- query %>%
    filter(x>=my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1]) %>%
    mutate(seg=iseg)

  # First data point scaled by derivatives (l'Hopital) --> did not work, use constant scaling
  my.query.add <- map.segment(segm.ref,
                              segm.query,
                              my.ref.lm[c(iseg,(iseg+1)),],
                              my.query.lm[c(iseg:(iseg+1)),],
                              scale_by=c(scale.use[iseg],scale.use[iseg+1])) %>%
    mutate(seg=iseg)
  my.query.scaled <- rbind(my.query.scaled,my.query.add)
}

segm.ref %>% ggplot(aes(x=x,y=y)) +
  geom_line(data=segm.ref,col='red') +
  geom_line(data=segm.query,col='blue') +
  geom_point(data=segm.ref[segm.ref$x==min(segm.ref$x),],col='red') +
  geom_point(data=segm.query[segm.query$x==min(segm.query$x),],col='blue') +
  geom_point(data=segm.ref[segm.ref$x==max(segm.ref$x),],col='red') +
  geom_point(data=segm.query[segm.query$x==max(segm.query$x),],col='blue') +
  # scale_y_log10() +
  theme_bw()

# Segment colroing, but one seg only!
my.query.scaled %>%
  ggplot(aes(x=x,y=y,col=factor(seg)))+
  geom_line()+   # one segment only
  geom_line(data=my.query.scaled,aes(x=x.scaled,y=y.scaled),lty=2) +
  geom_segment(data=my.query.scaled[100*c(1:14),],aes(x=x,y=y,xend=x.scaled,yend=y.scaled),
               arrow = arrow(length = unit(0.2, "cm")))

# --------------- TRANSFORMATION OF OBSERVATIONS ------------

# Steps
# 1. Assign segments to each query sample data point and to reference curve
# 2. Determine obs to query curve factor (and obtain prop error and/or add error?)
# 3. Determine obs y-scaling
# 4. For each segment determine obs x-scaling
# 5. Determine obs new x-position
# 6. Determine obs y-scaling ADDITIVE OR PROPORTIONAL

# Initialize
indivsam.query$x.scaled <- NA
indivsam.query$y.scaled <- NA

# 1. Assign segments to each query sample data point and to reference curve
indivsam.query$seg <- ref$seg <- query$seg <- NA
# if(nseg == 1)  # Why? if nseg=1 then my.query.lm has 2 rows!
# {
#   indivsam.query$seg <- ref$seg <- query$seg <- 1
# }
# if(nseg > 1)
# {
for(iseg in c(1:nseg))
{
  indivsam.query <- indivsam.query %>%
    mutate(seg = ifelse(x>=my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1], iseg, seg))
  # Full query and reference curves:
  ref        <- ref        %>%
    mutate(seg = ifelse(x>=my.ref.lm$x[iseg]   & x<=my.ref.lm$x[iseg+1], iseg, seg))
  query      <- query      %>%
    mutate(seg = ifelse(x>=my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1], iseg, seg))
}
# }
#
if(sum(is.na(indivsam.query$seg)>0))                                         stop("Error, not all samples could be assigned to a segment")
if(sum(is.na(ref$seg[ref$x<=my.ref.lm$x[my.ref.lm$type=='end']])>0))         stop("Error, not all ref data could be assigned to a segment")
if(sum(is.na(query$seg[query$x<=my.query.lm$x[my.query.lm$type=='end']])>0)) stop("Error, not all query data could be assigned to a segment")

ref %>% ggplot(aes(x=x,y=y)) +
  geom_line(data=ref[!is.na(ref$seg),],col='red') +
  geom_line(data=query[!is.na(query$seg),],col='blue') +
  geom_point(data=ref[ref$x==min(ref$x[!is.na(ref$seg)]),],col='red') +
  geom_point(data=query[query$x==min(query$x[!is.na(query$seg)]),],col='blue') +
  geom_point(data=ref[ref$x==max(ref$x[!is.na(ref$seg)]),],col='red') +
  geom_point(data=query[query$x==max(query$x[!is.na(query$seg)]),],col='blue') +
  # scale_y_log10() +
  theme_bw()

# 2. Determine y scaling
#      obs to query curve factor (and obtain prop/add error?)
# indivsam.query$ypred      <- approx(query$x,query$y,xout=indivsam.query$x,rule=2)$y
# Do not understand warning message:
#   In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) :
#   collapsing to unique 'x' values
# We already have PRED
indivsam.query$ypred <-  indivsam.query$PRED

unique(query$ALB)
query %>% ggplot(aes(x=x,y=PRED))+
  geom_line(col='blue')+
  geom_line(data=ref,col='red')

# ?? I think we can use this as an alternative way to account for prop and additive error: to be investigated
# indivsam.query$prop.error <- (indivsam.query$y - indivsam.query$ypred)/indivsam.query$ypred
# indivsam.query$add.error  <-  indivsam.query$y - indivsam.query$ypred

# 3. Determine obs new x-position
indivsam.query$x.scaled         <- approx(my.query.scaled$x, my.query.scaled$x.scaled, xout = indivsam.query$x)$y

# 4. Determine obs y-scaling ADDITIVE OR PROPORTIONAL
if(PROP_TR)
{
  # PROPORTIONAL - log scale addition
  ref   <- ref   %>% mutate(ylog = ifelse(y>0,log(y),NA))
  query <- query %>% mutate(ylog = ifelse(y>0,log(y),NA))
  indivsam.query <- indivsam.query %>% mutate(ylog = ifelse(y>0,log(y),NA))

  # Position on query curve
  indivsam.query$ylog.query       <- ifelse(!is.na(indivsam.query$ylog),
                                            approx(x=query$x, y=query$ylog, xout = indivsam.query$x)$y,
                                            NA)
  # log difference indivsam.query to query curve
  indivsam.query$ylog.diff        <- indivsam.query$ylog - indivsam.query$ylog.query

  # Position on ref curve (position at x.scaled!!)
  indivsam.query$ylog.ref         <- ifelse(!is.na(indivsam.query$ylog),
                                            approx(x=ref$x, y=ref$ylog, xout = indivsam.query$x.scaled)$y,
                                            NA)
  # New position
  indivsam.query$ylog.scaled      <- indivsam.query$ylog.diff + indivsam.query$ylog.ref
  indivsam.query$y.scaled         <- exp(indivsam.query$ylog.scaled)
}

if(ADD_TR)
{
  # ADDITIVE - normal scale addition

  # Position on query curve
  indivsam.query$y.query       <- ifelse(!is.na(indivsam.query$y),
                                         approx(x=query$x, y=query$y, xout = indivsam.query$x)$y,
                                         NA)
  # difference indivsam.query to query curve
  indivsam.query$y.diff        <- indivsam.query$y - indivsam.query$y.query

  # Position on ref curve (position at x.scaled!!)
  indivsam.query$y.ref         <- ifelse(!is.na(indivsam.query$ylog),
                                         approx(x=ref$x, y=ref$ylog, xout = indivsam.query$x.scaled)$y,
                                         NA)
  # New position
  indivsam.query$y.scaled      <- indivsam.query$y.diff + indivsam.query$y.ref
}

# Same differences on normal scale NOT preserved if ypred is not scaled as additive

# For plotting:
query.all <- rbind(query.all,query)

ref.lm.all   <- rbind(ref.lm.all,my.ref.lm %>% mutate(ALB.n = iALB))
query.lm.all <- rbind(query.lm.all,my.query.lm %>% mutate(ALB.n = iALB))

query.scaled.all <- rbind(query.scaled.all,my.query.scaled)

indivsam.query.all <- rbind(indivsam.query.all,indivsam.query)

# }  # for each covariate iALB


# ------------ Done ------------------




# --------------------------------------------------------

p0 <- ref %>%
  ggplot(aes(x=x,y=y,group=ALB.n,col=factor(round(ALB,1))))+
  geom_line(lwd=1)+
  geom_line(data=query.all) +
  labs(title=paste0(model,". Covariate impacted typical curves "),
       caption=script) +
  guides(color=guide_legend(title="Albumin"))

p0.obs <- ref %>%
  ggplot(aes(x=x,y=y,group=ALB.n,col=factor(round(ALB,1))))+
  geom_line(lwd=1)+
  geom_point(data=indivsam.query.all) +
  geom_line(data=query.all) +
  geom_point(data=indivsam.ref) +
  coord_cartesian(xlim=c(0,xstop)) +
  labs(title=paste0(model,". Covariate impacted typical curves with observations"),
       caption=script) +
  guides(color=guide_legend(title="Albumin"))

# Show landmarks
ref.lm.all <- ref.lm.all %>%
  left_join(ref[,c('ALB.n','ALB')], by='ALB.n')
query.lm.all <- query.lm.all %>%
  left_join(query[,c('ALB.n','ALB')], by='ALB.n')

p0.lm <- p0 +
  geom_point(data=ref.lm.all,aes(x=x,y=y),pch=21,fill='white',size=2)+
  geom_point(data=query.lm.all,aes(x=x,y=y),pch=21,fill='white',size=2) +
  labs(title=paste0(model,". Covariate impacted typical curves with landmarks"),
       caption=script) +
  guides(color=guide_legend(title="Albumin"))


if(model=="sigmoid"  ) plot(p0.lm + scale_x_log10())


# All information is available:
p1 <- query.scaled.all %>%
  ggplot(aes(x=x,y=y,group=ALB.n)) +
  geom_line() +
  geom_line(aes(x=x.scaled,y=y.scaled),col='red')+
  labs(title="Query original (black) and query transformed (red) curves",
       caption=script)
p1

# Some ypred points from query to ref curve
p2 <- p1 +
  geom_point(data=query.scaled.all[c(5,10,20),])+
  geom_point(data=query.scaled.all[c(5,10,20),],aes(x=x.scaled,y=y.scaled),col='red')+
  coord_cartesian(xlim=c(0,xstop)) +
  labs(title="Query original (black) and query transformed (red) curves with some data points",
       caption=script)

# Check of segments
p3 <- p1 +
  geom_line(data=query,aes(x=x,y=y,col=factor(seg)),lwd=1.5) +
  geom_line(data=ref,aes(x=x,y=y,col=factor(seg)),lwd=1.5) +
  coord_cartesian(xlim=c(0,xstop)) +
  labs(title="Segments",
       caption=script) +
  guides(color=guide_legend(title="Segment"))

# Check of PREDS
p4 <- p1 + geom_point(data=indivsam.query.all,aes(x=x,y=ypred,col=factor(seg))) +
  coord_cartesian(xlim=c(0,xstop)) +
  labs(title="Preds of samples",
       caption=script) +
  guides(color=guide_legend(title="Segment"))

# Check
p5 <- p1 + geom_point(data=indivsam.query.all,aes(x=x,y=ypred.scaled,col=factor(seg))) +
  geom_line(data=ref,aes(x=x,y=y,col=factor(seg)),lwd=1) +
  coord_cartesian(xlim=c(0,xstop)) +
  labs(title="Y-scaled Preds of samples",
       caption=script) +
  guides(color=guide_legend(title="Segment"))

# Plot up to last original or extended query data point
max.x    <- max(indivsam.query.all$x.scaled, indivsam.query.all$x, xstop)
all.y    <- c(indivsam.query.all$y.scaled, indivsam.query.all$y)
max.y    <- max(all.y)
lowest.y <- min(all.y[all.y!=0])

# Check Vachette transformation PREDS
p6 <- p1 +
  geom_line(data=ref,aes(x=x,y=y),lwd=1.5) +
  geom_point(data=indivsam.query.all,aes(x=xpred.scaled,y=ypred.scaled,col=factor(seg))) +
  coord_cartesian(xlim=c(0,max.x)) +
  labs(title="Y-scaled and X-scaled Preds of samples",
       caption=script) +
  guides(color=guide_legend(title="Segment"))


# Normal scale Vachette transformation
p7 <- query.scaled.all %>%
  ggplot(aes(x=x,y=y,group=ALB.n,col=factor(round(ALB,1)))) +
  geom_line() +
  geom_line(aes(x=x.scaled,y=y.scaled),col='red')+
  geom_line(data=ref,aes(x=x,y=y),lwd=1) +
  geom_point(data=indivsam.query.all) +
  geom_point(data=indivsam.query.all,aes(x=x.scaled,y=y.scaled)) +
  geom_segment(data=indivsam.query.all,aes(x=x,y=y,xend=x.scaled,yend=y.scaled),
               arrow = arrow(length = unit(0.2, "cm")),col='darkgrey') +
  geom_vline(xintercept = xstop, col='grey', lty=2) +
  coord_cartesian(xlim=c(0,max.x)) +
  labs(title="Vachette transformation",
       caption=script) +
  guides(color=guide_legend(title="Albumin"))
# Looks good

# Log scale Vachette transformation
if(model=='sigmoid') lowest.y <- min(long$y)
# if(model=='pembro')  lowest.y <- 10
p8 <- p7 +
  scale_y_log10(limits=c(lowest.y,NA))

# Show that distances are preserved in log scale
df <- NULL
iALB <- 1
# for(iALB in unique(query.scaled.all$ALB.n))
# {
idf <- indivsam.query.all %>% filter(ALB.n == iALB)

# Before transformation, distance to query curve, proportional on log scales:
if(ADD_TR)
{
  # y is IIV-corrected (or not!) data point
  idf$y.diff             <- idf$ypred - idf$y

  # After transformation, distance to reference curve, log scales
  idf$ypred.scaled       <- approx(ref$x, ref$y, xout=idf$x.scaled)$y
  idf$y.diff.scaled      <- idf$ypred.scaled - idf$y.scaled
}
if(PROP_TR)
{
  # y is IIV-corrected (or not!) data point
  idf$log.ypred[idf$ypred>0] <- log(idf$ypred[idf$ypred>0])
  idf$log.y[idf$y>0]         <- log(idf$y[idf$y>0])
  idf$log.y.diff             <- idf$log.ypred - idf$log.y

  # After transformation, distance to reference curve, log scales
  idf$ypred.scaled                         <- approx(ref$x, ref$y, xout=idf$x.scaled)$y
  idf$log.ypred.scaled[idf$ypred.scaled>0] <- log(idf$ypred.scaled[idf$ypred.scaled>0])
  idf$log.y.scaled[idf$y.scaled>0]         <- log(idf$y.scaled[idf$y.scaled>0])
  idf$log.y.diff.scaled                    <- idf$log.ypred.scaled - idf$log.y.scaled
}

# Collect
df <- rbind(df,idf)
# }
if(PROP_TR)
{
  p8.log.y.diff <- df %>%
    ggplot()+
    geom_abline(slope=1,col='black')+
    geom_point(aes(x=log.y.diff,y=log.y.diff.scaled,col=factor(round(ALB,1))))+
    labs(x="Log diff original data to typical query",
         y="Log diff Vachette transformed data to typical reference",
         title="Log-scale difference original data versus typical query curve to\nlog-scale Vachette transformed data to typical reference curve",
         caption=script) +
    guides(color=guide_legend(title="Albumin"))
}
if(ADD_TR)
{
  p8.y.diff <- df %>%
    ggplot()+
    geom_abline(slope=1,col='black')+
    geom_point(aes(x=y.diff,y=y.diff.scaled,col=factor(round(ALB,1))))+
    labs(x="Difference original data to typical query",
         y="Difference Vachette transformed data to typical reference",
         title="Difference original data versus typical query curve to\nVachette transformed data to typical reference curve",
         caption=script) +
    guides(color=guide_legend(title="Albumin"))
}

if(IIV_CORR)  write.csv(indivsam.ref,"indivsam-ref-iiv-corr.csv")
if(!IIV_CORR) write.csv(indivsam.ref,"indivsam-ref-no-iiv-corr.csv")
if(IIV_CORR)  write.csv(indivsam.query.all,"indivsam.query.all-iiv-corr.csv")
if(!IIV_CORR) write.csv(indivsam.query.all,"indivsam.query.all-no-iiv-corr.csv")
if(IIV_CORR)  write.csv(ref,"ref-iiv-corr.csv")
if(!IIV_CORR) write.csv(ref,"ref-no-iiv-corr.csv")

indivsam.ref.iiv.corr <- read.csv("indivsam-ref-iiv-corr.csv")
indivsam.ref.no.iiv.corr <- read.csv("indivsam-ref-no-iiv-corr.csv")
indivsam.query.all.iiv.corr <- read.csv("indivsam.query.all-iiv-corr.csv")
indivsam.query.all.no.iiv.corr <- read.csv("indivsam.query.all-no-iiv-corr.csv")
ref.iiv.corr <- read.csv("ref-iiv-corr.csv")
ref.no.iiv.corr <- read.csv("ref-no-iiv-corr.csv")

head(indivsam.ref.iiv.corr)
head(indivsam.ref.no.iiv.corr)

max.y    <- 300
lowest.y <- 10

# Vachette plot, both transformed query observations together with reference observations
p99.orig <- ref %>%
  ggplot() +
  geom_line(data=ref,aes(x=time,y=PRED),lwd=1) +
  geom_point(data=indivsam.ref.iiv.corr,      aes(x=time,y=OBS),col='red') +
  geom_point(data=indivsam.query.all.iiv.corr,aes(x=time,y=OBS),col='blue') +
  # geom_point(data=indivsam.ref.no.iiv.corr,      aes(x=time,y=OBS),col='orange') +
  # geom_point(data=indivsam.query.all.no.iiv.corr,aes(x=time,y=OBS),col='cyan') +
  coord_cartesian(xlim=c(0,max.x),ylim=c(0,max.y)) +
  labs(title="Final Vachette plot",
       caption=script)

# Vachette plot, both transformed query observations together with reference observations
p99 <- ref %>%
  ggplot() +
  geom_line(data=ref,aes(x=x,y=y),lwd=1) +
  # geom_point(data=indivsam.ref.iiv.corr,aes(x=x,y=y),col='red') +
  geom_point(data=indivsam.query.all.iiv.corr,aes(x=x.scaled,y=y.scaled),col='purple') +
  # geom_point(data=indivsam.ref.no.iiv.corr,aes(x=x,y=y),col='orange') +
  geom_point(data=indivsam.query.all.no.iiv.corr,aes(x=x.scaled,y=y.scaled),col='hotpink2') +
  coord_cartesian(xlim=c(0,max.x),ylim=c(0,max.y)) +
  labs(title="Final Vachette plot",
       caption=script)+
  render


# Vachette plot, both transformed query observations together with reference observations
p9.orig.cov.color <- ref %>%
  ggplot(aes(x=time,y=PRED,col=factor(round(ALB,1)))) +
  geom_line(data=ref,lwd=1) +
  geom_point(data=indivsam.ref,      aes(x=time,y=OBS),pch=19) +
  geom_point(data=indivsam.ref,      aes(x=time,y=OBS),pch=1,col='black') +
  geom_point(data=indivsam.query.all,aes(x=time,y=OBS),pch=19) +
  coord_cartesian(xlim=c(0,max.x),ylim=c(0,max.y)) +
  labs(title="Final Vachette plot",
       x='Time',
       y='Observations',
       caption=script)+
  guides(color=guide_legend(title="Albumin"))

# Vachette plot, both transformed query observations together with reference observations
p9.orig <- ref %>%
  ggplot() +
  geom_line(data=ref,aes(x=time,y=PRED),lwd=1) +
  geom_point(data=indivsam.ref,      aes(x=time,y=OBS),col='red') +
  geom_point(data=indivsam.query.all,aes(x=time,y=OBS),col='blue') +
  coord_cartesian(xlim=c(0,max.x),ylim=c(0,max.y)) +
  labs(title="Final Vachette plot",
       x='Time',
       y='Observations',
       caption=script)

# Vachette plot, both transformed query observations together with reference observations
if(IIV_CORR)  ylabel <- 'IIV-corrected and transformed query observations'
if(!IIV_CORR) ylabel <- 'Transformed observations (no IIV correction)'

p9.cov.color <- ref %>%
  ggplot(aes(x=time,y=y,col=factor(round(ALB,1)))) +
  geom_line(data=ref,lwd=1) +
  geom_point(data=indivsam.ref,      aes(x=time,y=y),pch=19) +
  geom_point(data=indivsam.ref,      aes(x=time,y=y),pch=1,col='black') +
  geom_point(data=indivsam.query.all,aes(x=x.scaled,y=y.scaled),pch=19) +
  coord_cartesian(xlim=c(0,max.x),ylim=c(0,max.y)) +
  labs(title="Final Vachette plot",
       x='Time',
       y=ylabel,
       caption=script)+
  guides(color=guide_legend(title="Albumin"))


p9 <- ref %>%
  ggplot() +
  geom_line(data=ref,aes(x=x,y=y),lwd=1) +
  geom_point(data=indivsam.ref,aes(x=x,y=y),col='red') +
  geom_point(data=indivsam.query.all,aes(x=x.scaled,y=y.scaled),col='purple') +
  coord_cartesian(xlim=c(0,max.x),ylim=c(0,max.y)) +
  labs(title="Final Vachette plot",
       x='Time',
       y=ylabel,
       caption=script)

# Vachette plot log scale, both transformed query observations together with reference observations
p10.orig.cov.color <- p9.orig.cov.color +
  scale_y_log10()+
  coord_cartesian(xlim=c(0,max.x),ylim=c(lowest.y,max.y))

p10.cov.color <- p9.cov.color +
  scale_y_log10()+
  coord_cartesian(xlim=c(0,max.x),ylim=c(lowest.y,max.y))

# Vachette plot log scale, both transformed query observations together with reference observations
p10.orig <- p9.orig +
  scale_y_log10()+
  coord_cartesian(xlim=c(0,max.x),ylim=c(lowest.y,max.y))

p10 <- p9 +
  scale_y_log10()+
  coord_cartesian(xlim=c(0,max.x),ylim=c(lowest.y,max.y))

# Vachette plot, both transformed query observations together with reference observations
# With original data point in grey
p11 <- ref %>%
  ggplot() +
  geom_line(data=ref,aes(x=x,y=y),lwd=1) +
  geom_point(data=indivsam.ref,aes(x=x,y=y),col='red') +
  geom_point(data=indivsam.query.all,aes(x=x,y=y),col='grey') +
  geom_point(data=indivsam.query.all,aes(x=x.scaled,y=y.scaled),col='purple') +
  coord_cartesian(xlim=c(0,max.x)) +
  labs(title="Final Vachette plot",
       caption=script)

# Vachette plot, both transformed query observations together with reference observations
# With original data point in grey
p11 <- ref %>%
  ggplot() +
  geom_line(data=ref,aes(x=x,y=y),lwd=1) +
  geom_point(data=indivsam.ref,aes(x=x,y=y),col='red') +
  geom_point(data=indivsam.query.all,aes(x=x,y=y),col='grey') +
  geom_point(data=indivsam.query.all,aes(x=x.scaled,y=y.scaled),col='purple') +
  coord_cartesian(xlim=c(0,max.x)) +
  labs(title="Final Vachette plot",
       caption=script)

# Vachette plot, both transformed query observations together with reference observations
# With original data point in grey w/o reference
p12 <- ref %>%
  ggplot() +
  geom_line(data=ref,aes(x=x,y=y),lwd=1) +
  geom_point(data=indivsam.query.all,aes(x=x,y=y),col='grey') +
  geom_point(data=indivsam.query.all,aes(x=x.scaled,y=y.scaled),col='purple') +
  coord_cartesian(xlim=c(0,max.x)) +
  labs(title="Final Vachette plot w/o reference observations",
       caption=script)


# Save plots:

pdf(paste0("../plots/vachette-main-",model,"-",my.covariate,"-",tag,extra,extra2,"-",setup,".pdf"))

plot(p0+render)
plot(p0+scale_y_log10()+render)
plot(p0.obs+render)
plot(p0.obs+scale_y_log10(limits=c(lowest.y,NA))+render)
plot(p0.lm+render)
plot(p0.lm+scale_y_log10()+render)
plot(p1+render)
plot(p1+scale_y_log10()+render)
plot(p7+render)
plot(p7+scale_y_log10(limits=c(lowest.y,NA))+render)
if (PROP_TR) plot(p8.log.y.diff+render)
if (ADD_TR)  plot(p8.y.diff+render)
plot(p9.orig+render)
plot(p9+render)
plot(p9.orig+render+scale_y_log10()+coord_cartesian(xlim=c(0,max.x),ylim=c(lowest.y,max.y)))
plot(p9+render+scale_y_log10()+coord_cartesian(xlim=c(0,max.x),ylim=c(lowest.y,max.y)))
plot(p9.orig.cov.color+render)
plot(p9.cov.color+render)
plot(p9.orig.cov.color+render+scale_y_log10()+coord_cartesian(xlim=c(0,max.x),ylim=c(lowest.y,max.y)))
plot(p9.cov.color+render+scale_y_log10()+coord_cartesian(xlim=c(0,max.x),ylim=c(lowest.y,max.y)))
plot(p11+render)
plot(p11+scale_y_log10(limits=c(lowest.y,NA))+render)
plot(p12+render)
plot(p12+scale_y_log10(limits=c(lowest.y,NA))+render)

dev.off()

# --- end ---

