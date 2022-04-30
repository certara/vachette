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
# v15 test use of simulation tolerances indirect response
#     Pembro model

# Clear memory
rm(list=ls())

tag      <- "v16"  
# model  <- "sigmoid"  
# model  <- "iv"  
model  <- "oral-absorption"
# model  <- "indirect-response"  
# model  <- "oral-two-dose"
# model <- "pembro"

script       <- paste0("vachette-main-",tag,".R")
setup        <- "A"     # choice of ref and query/queries. Choose B is ref/query should be swapped
my.covariate <- "none"  # To be changed for each model example

extra  <- "-delta01"      # extra info to add to output file names
mrgdelta <- 0.1

# Location of this script
home   <- dirname(rstudioapi::getActiveDocumentContext()$path) 
getwd()
setwd(home)

source('vachette-functions-v16.R')
source('models-v16.R')

set.seed(125)

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
# 
if (model == 'pembro')
{
  # my.ipred <- mcode("pembro", pembro.ipred)
  # my.pred  <- mcode("pembro", pembro.pred)
  
  # (Typical prediction)
  
  xstop         <- 14    # Actual observations end at xstop = 14 days
  xlast.user    <- 700    # User provided x value to max. simulate out - replaces expand.factor
  #expand.factor <- 10   # factor to expand simulated curve to find last X values for ref and query
  tolend        <- 0.001  # Max tolerance allowed minimizing difference open end ref and query
  tolnoise      <- 1e-8   # Max tolerance allowed between gridpoints to identify landmarks for stepsize=1
  step.x.factor <- 1.5  # Factor x-steps between which y-values are lower than the SCALED tolerance 
  ngrid.open.end<- 100  # Number of grid points to characterize open end curves
  

  # Select pembro covariate to assess:
  my.covariate <- "sex"
  
  # Original situation (query contraction to fit ref)
  # SEX covariate
  if(setup=="A") cov.ref        <- 1
  if(setup=="A") cov.query      <- 2
  # Alternative situation (query expansion to fit ref)
  if(setup=="B") cov.ref        <- 2
  if(setup=="B") cov.query      <- 1
  
  # ------
  
  # mod.read <- mread("../literature-examples/pembrolizumab/mrgsolve/pembro")
  mod.read <- mcode("pembro", pembro)
  
  
  #simulation settings
  dose <- 10                 #set pembrolizumab dose in mg/kg
  dose.interval <- c(14,21)  #set dosing interval in days
  dose.n <- 12
  inf.duration = 0.5/24
  maxtime = 25*7             #set end of simulation time interval
  subj.n = 2
  subj.n.iiv = 2188
  
  # sample.times1 <- c(0.5, 1, 2, 5, 10, 24, 48, 4*24, 8*24, 11*24, dose.interval[1]*24)/24
  # sample.times2 <- c(0:dose.n)*14
  # 
  # sample.times3 <- c(0.5, 1, 2, 5, 10, 24, 48, 4*24, 8*24, 11*24, 14*24, 17*24, dose.interval[2]*24)/24
  # sample.times4 <- c(0:dose.n)*21
  
  # samples <- merge(sample.times1, sample.times2) %>%
  #   mutate(times = x + y)
  # 
  # samples2 <- merge(sample.times3, sample.times4) %>%
  #   mutate(times = x + y)
  
  # ---- JL  start with single dose (single Q2W only)
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
  
  param(mod.read)
  # Model parameters (N=25):
  #   name     value   . name     value 
  # ADIAGN   1       | RESERR   -0.272
  # ALB      39.6    | SEX      1     
  # BECOG    1       | TVCL     0.22  
  # BSLD     89.6    | TVEXCLQ  0.595 
  # CLADIAGN 0.145   | TVEXV1V2 0.489 
  # CLALB    -0.907  | TVQ      0.795 
  # CLBECOGN -0.0739 | TVV1     3.48  
  # CLBSLD   0.0872  | TVV2     4.06  
  # CLEGFR   0.135   | V1ALB    -0.208
  # CLIPIP   0.14    | V1IPIP   0.0737
  # CLSEX    -0.152  | V1SEX    -0.134
  # EGFR     88.5    | WT       76.8  
  # IPIP     0       | .        .    
  
  #also update omega matrix to zero to simulate typical curves:
  mod     <- omat(mod.read, dmat(0,0))
  
  #individual covariate distributions (to match Table 2 in publication)
  sampling.n <- 1000
  ind.sex <- sample(c(1,2), size = sampling.n, replace = TRUE, prob = c(0.591, 0.409))
  ind.diag <- sample(c(1,2), size = sampling.n, replace = TRUE, prob = c(0.747, 0.253))
  ind.ecog <- sample(c(0,1), size = sampling.n, replace = TRUE, prob = c(0.576, 0.424))
  ind.ipi <- sample(c(0,1), size = sampling.n, replace = TRUE, prob = c(0.655, 0.345))
  ind.egfr <- exp(rnorm(sampling.n, mean = log(88.7), sd = 0.4))
  ind.alb <- exp(rnorm(sampling.n, mean = log(30), sd = 0.2))
  ind.bsld <- exp(rnorm(sampling.n, mean = log(86), sd = 0.5))
  
  #create simulation dataset:
  dosedat <- ev(ID = 1:subj.n, amt = dose, ii = dose.interval[1], addl = 49) %>%
    mutate(WT = 76.8,
           amt = dose * WT, 
           rate = amt/inf.duration,
           dose = dose)
  
  #typical covariate values
  tSEX = 1
  tIPIP = 0
  tALB = 39.6
  tEGFR = 88.47
  tBSLD = 89.60
  tBECOG = 1
  tADIAGN = 1
  
  covariates <- data.frame(tSEX, tIPIP, tALB, tEGFR, tBSLD, tBECOG, tADIAGN)
  if (my.covariate == 'sex')                icov <- 1
  if (my.covariate == 'IPI status')         icov <- 2
  if (my.covariate == 'albumin')            icov <- 3
  if (my.covariate == 'eGFR')               icov <- 4
  if (my.covariate == 'Baseline tumor burden') icov <- 5
  if (my.covariate == 'ECOG Status')        icov <- 6
  if (my.covariate == 'Cancer Type')        icov <- 7
  
  
  ind.dat.ori <- dosedat %>%
    dplyr::select(ID, WT) %>%
    mutate(SEX = tSEX,
           IPIP = tIPIP,
           ALB = tALB,
           EGFR = tEGFR,
           BSLD = tBSLD,
           BECOG = tBECOG,
           ADIAGN = tADIAGN)
  
  # ----- IIV model/dataset ----
  
  mod.iiv <- mod.read
  
  dosedat.iiv <- ev(ID = 1:sampling.n, amt = dose, ii = dose.interval, addl = 49) %>%
    mutate(WT = rnorm(sampling.n, 76.8, 8),
           amt = dose * WT,
           rate = amt/inf.duration,
           dose = dose)
  
  # Sex
  ind.dat.iiv <- dosedat.iiv %>%
    dplyr::select(ID, WT) %>%
    mutate(SEX = ind.sex,
           ADIAGN = tADIAGN,
           BECOG = tBECOG,
           IPIP = tIPIP,
           EGFR = tEGFR,
           ALB = tALB,
           BSLD = tBSLD)
  
  # -------------------------------
  
  cov.name <- str_replace(names(covariates)[icov], "t", "")
  
  if (cov.name == "SEX") {
    ind.data <- ind.dat.ori %>%
      mutate(SEX = range(ind.sex))
    colname = "Sex"
    header = "with different sexes"
  }
  if (cov.name == "IPIP") {
    ind.data <- ind.dat.ori %>%
      mutate(IPIP = range(ind.ipi))
    colname = "IPI Status"
    header = "with differing IPI status"
  }
  if (cov.name == "ALB") {
    ind.data <- ind.dat.ori %>%
      mutate(ALB = range(ind.alb))
    colname <- "Albumin (g/L)"
    header = "using different albumin levels"
  }
  if (cov.name == "EGFR") {
    ind.data <- ind.dat.ori %>%
      mutate(EGFR = range(ind.egfr))
    colname = "eGFR (ml/min/1.73 m2)"
    header = "using different eGFR values"
  }
  if (cov.name == "BSLD") {
    ind.data <- ind.dat.ori %>%
      mutate(BSLD = range(ind.bsld))
    colname = "Baseline tumor burden (mm)"
    header = "for differing tumor burden"
  }
  if (cov.name == "BECOG") {
    ind.data <- ind.dat.ori %>%
      mutate(BECOG = range(ind.ecog))
    colname = "ECOG Status"
    header = "for differing ECOG status"
  } 
  if (cov.name == "ADIAGN") {
    ind.data <- ind.dat.ori %>%
      mutate(ADIAGN = range(ind.diag))
    colname = "Cancer Type"
    header = "for different cancer types"
  }
  
  out.collect     <- NULL
  out.collect.iiv <- NULL
  
  for (iint in unique(dose.interval))
  {
    
    if (iint == 999) {
      times <- samples$times 
      dosedat <- dosedat %>%
        mutate(ii = 999)
    }
    
    if (iint == 21) {
      times <- samples2$times
      dosedat <- dosedat %>%
        mutate(ii = 21)
    }
    
    # ---- Typical -----
    sim <- mod %>%
      data_set(dosedat) %>%
      idata_set(ind.data) %>%
      carry_out(EVID) %>%
      Req(CONC, CONC_RUV, SEX, IPIP, ALB, EGFR, BSLD, BECOG, ADIAGN) %>%
      mrgsim(tad = TRUE, tgrid = times)
    
    output <- as.data.frame(sim) %>%
      filter(EVID == 0, time!=0) %>%   # time=0 causes problems at scaling
      dplyr::select(-EVID) %>%
      mutate(dosing.interval = paste0(iint, " days"))
    
    # Keep time=0 (to show issue when scaling)
    output0 <- as.data.frame(sim) %>%
      filter(EVID == 0) %>%   # time=0 causes problems at scaling
      dplyr::select(-EVID) %>%
      mutate(dosing.interval = paste0(iint, " days"))
    
    out.collect <- rbind(out.collect, output)
    
    # ---- Sampling, iiv ----
    sim.iiv <- mod.iiv %>%
      data_set(dosedat.iiv) %>%
      idata_set(ind.dat.iiv) %>%
      carry_out(EVID) %>%
      Req(CONC, CONC_RUV, SEX, IPIP, ALB, EGFR, BSLD, BECOG, ADIAGN) %>%
      mrgsim(tad = TRUE, tgrid = samples.iiv$times)
    
    output.iiv <- as.data.frame(sim.iiv) %>%
      filter(EVID == 0, time!=0) %>%   # time=0 causes problems at scaling
      dplyr::select(-EVID) 
    
    out.collect.iiv <- rbind(out.collect.iiv, output.iiv)
    
  }
  
  p.pembro1 <-  out.collect %>%
    filter(dosing.interval == "999 days") %>%
    ggplot() +
    geom_line(aes(x = time/7, y = CONC, color = factor(round(get(cov.name), 1)))) +
    labs(title = paste0("Pembrolizumab simulated typical PK profiles, ", header),
         x = "Time (weeks)",
         y = "Pembrolizumab Concentration (µg/mL)",
         color = colname) +
    facet_wrap(~dosing.interval, labeller = label_both, scales = "free") +
    theme_bw() +
    scale_y_log10()
  
  plot(p.pembro1)
  
  p.pembro2 <-  out.collect.iiv %>%
    ggplot() +
    geom_line(aes(x = time/7, y = CONC, group = ID, col=factor(SEX)))+
    labs(title = paste0("Pembrolizumab simulated typical PK profiles, ", header),
         x = "Time (weeks)",
         y = "Pembrolizumab Concentration (µg/mL)") +
    theme_bw() +
    scale_y_log10()

  # ---- Vachette ---
  # Start with one dose, SEX covariate only
  
  # Typical cuve for SEX
  long  <- as.data.frame(out.collect) %>%
    mutate(x=time, y=CONC, COV=SEX)                 # Standardize
  
  # all in standard parameters x,y
  longref   <- long %>% filter(COV==cov.ref)
  longquery <- long %>% filter(COV==cov.query)
  
  # (Representing observations)
  
  nsample <- 200
  
  # Sample nsample reference and nsample query data points:
  indivsam.query <- out.collect.iiv %>%
    filter(SEX==cov.query) %>%
    slice(sample(n(),nsample,replace=F)) %>% 
    mutate(x=time, y=CONC, COV=SEX)                 # Standardize
  
  indivsam.ref <- out.collect.iiv %>%
    filter(SEX==cov.ref) %>%
    slice(sample(n(),nsample,replace=F)) %>% 
    mutate(x=time, y=CONC, COV=SEX)                 # Standardize
  
  indivsam <- rbind(indivsam.query,indivsam.ref) %>%
    arrange(time) 
  
  # Plot typical + sample
  p.pembro3 <- longref %>%
    ggplot(aes(x=x,y=y)) + 
    geom_line(col='red') +
    geom_line(data=longquery,col='blue') +
    geom_point(data=indivsam.ref,col='red') +
    geom_point(data=indivsam.query,col='blue')
  
  # Plot typical + sample
  p.pembro4 <- output0 %>%
    ggplot(aes(x=time/7,y=CONC,col=factor(SEX))) + 
    geom_line()+
    labs(title = paste0("Pembrolizumab simulated typical PK profiles, ", header),
         x = "Time (weeks)",
         y = "Pembrolizumab Concentration (µg/mL)") +
    coord_cartesian(xlim=c(0,2)) +
    guides(color=guide_legend(title="SEX"))+
    render
  
  # Plot typical + sample
  p.pembro5 <- output %>%
    ggplot(aes(x=time/7,y=CONC,col=factor(SEX))) + 
    geom_line()+
    labs(title = paste0("Pembrolizumab simulated typical PK profiles, ", header),
         x = "Time (weeks)",
         y = "Pembrolizumab Concentration (µg/mL)") +
    coord_cartesian(xlim=c(0,2)) +
    guides(color=guide_legend(title="SEX"))+
    render
  
  pdf(paste0("../plots/pembro.pdf"))
  
  plot(p.pembro1)
  plot(p.pembro2)
  plot(p.pembro3)
  plot(p.pembro4)
  plot(p.pembro5)
  
  dev.off()
  
}  # if model=='pembro'

##########################################################
#                                                        #
#                    SIGMOID/EMAX                        #
#              DEFINE REFERENCE AND QUERY                #
#                                                        #
##########################################################

if (model == 'sigmoid')
{
  # (Typical prediction)
  
  xstop         <- 10**2.5  # Samples cannot be located "after" xstop
  xlast.user    <- 10**5  # User provided x value to max. simulate out - replaces expand.factor
  tolend        <- 0.001  # Max tolerance allowed minimizing difference open end ref and query
  tolnoise      <- 1e-8   # Max tolerance allowed between gridpoints to identify landmarks for stepsize=1
  step.x.factor <- 1.5    # Factor x-steps between which y-values are lower than the SCALED tolerance 
  ngrid.open.end<- 100    # Number of grid points to characterize open end curves
  # Original situation (query contraction to fit ref)
  if(setup=="A") cov.ref        <- 70
  if(setup=="A") cov.query      <- 30
  # Alternative situation (query expansion to fit ref)
  if(setup=="B") cov.ref        <- 30
  if(setup=="B") cov.query      <- 70
  
  my.covariate <- "weigth"
  
  # simulate pred
  # longref     <- longquery  <- data.frame(x = 10**(c(-40:40)/10))
  longref     <- longquery  <- data.frame(x = 10**(c((-10*log10(xlast.user)):(10*log10(xlast.user)))/10))
  emax        <- 10
  ic50.ref    <- ifelse(cov.ref==70,1,10)
  ic50.query  <- ifelse(cov.query==70,1,10)
  gamma       <- 1
  longref$y   <- sigmoid(longref$x, emax=emax, ic50=ic50.ref, gamma=gamma)
  longref$WT  <- cov.ref
  longref$COV <- cov.ref   # Standardize
  longquery$y <- sigmoid(longref$x, emax=emax, ic50=ic50.query, gamma=gamma)
  longquery$WT <- cov.query
  longquery$COV<- cov.query   # Standardize
  
  # Merge
  long <- rbind(longref,longquery)
  
  # (Representing observations)
  
  nsample <- 20
  # simulate ipred
  set.seed(20220131)
  xgrid.samples    <- 10**(c((-10*log10(xstop)):(10*log10(xstop)))/10)
  ngrid.samples    <- length(xgrid.samples)
  indivsam.ref     <- indivsam.ref    <- data.frame(x = xgrid.samples, WT=cov.ref)
  indivsam.query   <- indivsam.query  <- data.frame(x = xgrid.samples, WT=cov.query)
  indivsam         <- rbind(indivsam.ref, indivsam.query)
  indivsam$x       <- indivsam$x*rnorm((2*ngrid.samples),1,0.1)
  indivsam$iemax   <- 10*rnorm((2*ngrid.samples),1,0.1)
  indivsam$igamma  <- 1*rnorm((2*ngrid.samples),1,0.15)
  # Covariate effect:
  indivsam$iic50[indivsam$WT==cov.ref]   <- NA
  indivsam$iic50[indivsam$WT==cov.ref]   <- ic50.ref*rnorm(length(indivsam$iic50[indivsam$WT==cov.ref]),1,0.2)
  indivsam$iic50[indivsam$WT==cov.query] <- ic50.query*rnorm(length(indivsam$iic50[indivsam$WT==cov.query]),1,0.2)
  
  # Standardize:
  indivsam$COV <- indivsam$WT
  
  # Response
  indivsam$y   <- sigmoid(indivsam$x, emax=indivsam$iemax, ic50=indivsam$iic50, gamma=indivsam$igamma)
  
  # Sample nsample reference and nsample query data points:
  indivsam.query <- indivsam %>%
    filter(COV==cov.query) %>% 
    slice(sample(n(),nsample,replace=F))
  indivsam.ref   <- indivsam %>%
    filter(COV==cov.ref) %>% 
    slice(sample(n(),nsample,replace=F))
  indivsam <- rbind(indivsam.query,indivsam.ref) %>%
    arrange(x)
  
  # ---- V2ACHER transformation of x for comparison ----
  longquery <- longquery %>% 
    mutate(x.tr = ifelse(WT==cov.query,x/(ic50.query/ic50.ref),x))
  # No y-scaling needed in this case
  
  indivsam <- indivsam %>% 
    mutate(x.tr = ifelse(WT==cov.query,x/(ic50.query/ic50.ref),x))
  # No y-scaling needed in this case
  # ------------------------------------------
  
  p0 <- longref %>% 
    ggplot(aes(x=x,y=y))+
    geom_line(col='red')+
    geom_line(data=longquery,col='blue') +
    labs(title=paste0(model," Extrapolated Query (blue) and Reference (red) ",setup),
         caption=script)
  p0log <- longref %>% 
    ggplot(aes(x=x,y=y))+
    geom_line(col='red')+
    geom_line(data=longquery,col='blue') +
    scale_x_log10() +
    labs(title=paste0(model," Extrapolated Query (blue) and Reference (red) ",setup),
         caption=script)
  p0logV2acher <- longref %>% 
    ggplot(aes(x=x,y=y))+
    geom_line(col='red')+
    geom_line(data=longquery,aes(x=x.tr,y=y),col='blue',lty=2) +
    scale_x_log10() +
    labs(title=paste0(model," Extrapolated Query (blue) and Reference (red) ",setup),
         caption=script)
  p0samples <- longref %>% 
    ggplot(aes(x=x,y=y))+
    geom_line(col='red')+
    geom_line(data=longquery,col='blue') +
    geom_point(data=indivsam[indivsam$COV==cov.ref,],col='red') +
    geom_point(data=indivsam[indivsam$COV==cov.query,],col='blue') +
    scale_x_log10() +
    labs(title=paste0(model," Extrapolated Query (blue) and Reference (red) ",setup),
         caption=script)
  
}

##########################################################
#                                                        #
#                       IV                               #
#              DEFINE REFERENCE AND QUERY                #
#                                                        #
##########################################################


if (model == 'iv')
{
  my.ipred <- mcode("iv.ipred", iv.ipred)
  my.pred  <- mcode("iv.pred", iv.pred)
  
  # (Typical prediction)
  
  xstop         <- 12   # Actual observations end at xstop = 48 hr
  xlast.user    <- 72  # User provided x value to max. simulate out - replaces expand.factor
  #expand.factor <- 10   # factor to expand simulated curve to find last X values for ref and query
  tolend        <- 0.001  # Max tolerance allowed minimizing difference open end ref and query
  tolnoise      <- 1e-8   # Max tolerance allowed between gridpoints to identify landmarks for stepsize=1
  step.x.factor <- 1.5  # Factor x-steps between which y-values are lower than the SCALED tolerance 
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
  
  # simulate pred
  long.out <- my.pred %>%
    ev(amt=100,addl=0,ii=24,cmt=1) %>%
    idata_set(data) %>%
    mrgsim(delta=mrgdelta,end=step.x.factor*xlast.user) 
  
  long  <- as.data.frame(long.out) %>%
    filter(!(time==0 & CENT==0)) %>%     # No Dose rec
    mutate(x=time, y=CP, COV=WT)         # Standardize
  
  # all in standard parameters x,y
  longref   <- long %>% filter(COV==cov.ref)
  longquery <- long %>% filter(COV==cov.query)
  
  # (Representing observations)
  
  nsample <- 20
  
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
    filter(!(time==0 & CENT==0)) %>%     # No Dose rec
    mutate(x=time, y=CP, COV=WT)         # Standardize
  
  # Sample nsample reference and nsample query data points:
  indivsam.query <- indiv %>%
    filter(COV==cov.query) %>%
    slice(sample(n(),nsample,replace=F))
  indivsam.ref <- indiv %>%
    filter(COV==cov.ref) %>%
    slice(sample(n(),nsample,replace=F))
  indivsam <- rbind(indivsam.query,indivsam.ref) %>%
    arrange(time)
  
}

##########################################################
#                                                        #
#                    ORAL ABSORPTION                     #
#              DEFINE REFERENCE AND QUERY                #
#                                                        #
##########################################################


if (model == 'oral-absorption')
{
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
}

##########################################################
#                                                        #
#              ORAL ABSORPTION  - TWO DOSE               #
#              DEFINE REFERENCE AND QUERY                #
#                                                        #
##########################################################


if (model == 'oral-two-dose')
{
  my.ipred <- mcode("oral1cmt.ipred", oral1cmt.ipred)
  my.pred  <- mcode("oral1cmt.pred", oral1cmt.pred)
  
  # ---------------- SIMULATED TYPICAL/POPULATION CURVES -------------
  
  # (Typical prediction)
  
  xstop         <- 96   # Actual observations end at xstop = 48 hr
  xlast.user    <- 640  # User provided x value to max. simulate out - replaces expand.factor
  #expand.factor <- 10   # factor to expand simulated curve to find last X values for ref and query
  tolend        <- 0.001  # Max tolerance allowed minimizing difference open end ref and query
  tolnoise      <- 1e-8   # Max tolerance allowed between gridpoints to identify landmarks for stepsize=1
  step.x.factor <- 1.5  # Factor x-steps between which y-values are lower than the SCALED tolerance 
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
    ev(amt=100,addl=1,ii=48,cmt=1) %>%
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
    ev(amt=100,addl=1,ii=48,cmt=1) %>%
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
}

##########################################################
#                                                        #
#               INDIRECT RESPONSE EXAMPLE                #
#              DEFINE REFERENCE AND QUERY                #
#                                                        #
##########################################################


if (model == 'indirect-response')
{
  
  my.ipred <- mcode("indirect.response.ipred", indirect.response.ipred)
  my.pred  <- mcode("indirect.response.pred", indirect.response.pred)
  
  # ---------------- SIMULATED TYPICAL/POPULATION CURVES -------------
  
  # (Typical prediction)
  
  xstop         <- 48   # Actual curves end at xstop = 48 hr
  xlast.user    <- 96   # User provided x value to max. simulate out - replaces expand.factor
  #expand.factor <- 2.5  # factor to expand simulated curve to find last X values for ref and query
  tolend        <- 0.001  # Max tolerance allowed minimizing difference open end ref and query
  tolnoise      <- 1e-8   # Max tolerance allowed between gridpoints to identify landmarks for stepsize=1
  step.x.factor <- 1.5  # Factor x-steps between which y-values are lower than the SCALED tolerance 
  ngrid.open.end<- 100  # Number of grid points to characterize open end curves
  # Original situation (query contraction to fit ref)
  if(setup=="A") cov.ref        <- 70
  if(setup=="A") cov.query      <- 30
  # Alternative situation (query expansion to fit ref)
  if(setup=="B") cov.ref        <- 30
  if(setup=="B") cov.query      <- 70
  mrgdelta     <- 0.1
  
  my.covariate <- "weigth"
  
  # population data
  data <- expand.idata(ID=1,WT=c(cov.query,cov.ref))
  
  # simulate pred
  long.out <- my.pred %>%
    init(RESP = 10) %>%    # =TVBSLN
    ev(amt=100,addl=0,ii=24,cmt=1) %>%
    idata_set(data) %>%
    # mrgsim(delta=mrgdelta,end=step.x.factor*xlast.user) 
    mrgsim(delta=0.1,end=step.x.factor*xlast.user,atol=1E-20,rtol=1E-10)
  # Problems when delta=0.25 (small stepsize --> identifying f2)
  
  long  <- as.data.frame(long.out)%>%
    filter(!(time==0 & CENT==0)) %>%       # No Dose rec
    mutate(x=time, y=RESP, COV=WT)         # Standardize
  
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
    init(RESP = 10) %>%    # =TVBSLN
    ev(amt=100,addl=0,ii=24,cmt=1) %>%
    idata_set(idata) %>%
    mrgsim(delta=mrgdelta,end=xstop)
  
  indiv <- as.data.frame(imy.out) %>%
    filter(!(time==0 & CENT==0)) %>%       # No Dose rec
    mutate(x=time, y=RESP, COV=WT)         # Standardize
  
  # Sample nsample reference and nsample query data points:
  indivsam.query <- indiv %>%
    filter(WT==cov.query) %>%
    slice(sample(n(),nsample,replace=F))
  indivsam.ref <- indiv %>%
    filter(WT==cov.ref) %>%
    slice(sample(n(),nsample,replace=F))
  indivsam <- rbind(indivsam.query,indivsam.ref) %>%
    arrange(time)
}

# --------------------------------------------------------

p0 <- longref %>% 
  ggplot(aes(x=x,y=y))+
  geom_line(col='red')+
  geom_line(data=longquery,col='blue') +
  labs(title=paste0(model," Extrapolated Query (blue) and Reference (red) ",setup),
       caption=script)

# ------ explore smoothing method -----

# explore.sg(long[long$COV==cov.ref,]$x,long[long$COV==cov.ref,]$y,
#            xmin=42,xmax=56,ymin=-1,ymax=1,
#            w0=3,w1=5,w2=7,p0=0,p1=1,p2=2,tag='two-dose-1')
# explore.sg(long[long$COV==cov.ref,]$x,long[long$COV==cov.ref,]$y,
#            xmin=47,xmax=49,ymin=-1,ymax=1,
#            w0=3,w1=5,w2=7,p0=0,p1=1,p2=2,tag='two-dose-2')
# explore.sg(long[long$COV==cov.ref,]$x,long[long$COV==cov.ref,]$y,
#            xmin=47,xmax=49,ymin=-1,ymax=1,
#            w0=13,w1=17,w2=21,p0=0,p1=1,p2=2,tag='two-dose-3')

# --- Landmark detection -----

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

pdf(paste0("../plots/vachette-main-",model,"-",my.covariate,"-",tag,extra,"-",setup,".pdf"))

plot(p0+render)
plot(p0.lm.init+render)
plot(p0.lm+render)
plot(p.scaling1+render)
plot(p.scaling2+render)
plot(p1+render)
plot(p2+render)
plot(p3+render)
plot(p4+render)
plot(p5+render)
plot(p6+render)
plot(p7+render)
plot(p8+render)
plot(p8.log.y.diff+render)
plot(p9+render)
plot(p10+render)

if(model=='sigmoid')
{
  plot(p0+render+scale_x_log10())
  plot(p0+render+coord_cartesian(xlim=c(0,100)))
  plot(p7+render+coord_cartesian(xlim=c(0.0001,NA)) + scale_x_log10())
  plot(p8+render+coord_cartesian(xlim=c(0.0001,NA)) + scale_x_log10())
  plot(p9+render+coord_cartesian(xlim=c(0.0001,NA)) + scale_x_log10())
  plot(p10+render+coord_cartesian(xlim=c(0.0001,NA)) + scale_x_log10())
}
if(model=='indirect-response')
{
  plot(p0.lm+render+annotate("text",x=40,y=1,label=paste0("Delta=",mrgdelta))+coord_cartesian(xlim=c(0,48),ylim=c(0,11)))
  plot(p7+render+annotate("text",x=40,y=1,label=paste0("Delta=",mrgdelta))+guides(color="none"))
}
if(model=='pembro')
{
  plot(p0+scale_y_log10()+render)
  plot(p0.lm.init+scale_y_log10()+render)
  plot(p0.lm+scale_y_log10()+render)
  plot(p1+scale_y_log10()+render)
  plot(p2+scale_y_log10()+render)
  plot(p3+scale_y_log10()+render)
  plot(p4+scale_y_log10()+render)
  plot(p5+scale_y_log10()+render)
  plot(p6+scale_y_log10()+render)
  plot(p8+coord_cartesian(xlim=c(0,xstop),ylim=c(30,NA))+render)
  plot(p10+coord_cartesian(xlim=c(0,xstop),ylim=c(30,NA))+render)
}
if(model=='sigmoid') # Add V2acher scaling
{
  plot(p7+render+coord_cartesian(xlim=c(0.0001,NA)) + scale_x_log10() + 
         geom_point(data=indivsam[indivsam$COV==cov.query,],aes(x=x.tr,y=y), pch=1, size=4))
  plot(p8+render+coord_cartesian(xlim=c(0.0001,NA)) + scale_x_log10() + 
         geom_point(data=indivsam[indivsam$COV==cov.query,],aes(x=x.tr,y=y), pch=1, size=4))
  plot(p9+render+coord_cartesian(xlim=c(0.0001,NA)) + scale_x_log10() + 
         geom_point(data=indivsam[indivsam$COV==cov.query,],aes(x=x.tr,y=y), pch=1, size=4))
  plot(p10+render+coord_cartesian(xlim=c(0.0001,NA)) + scale_x_log10() + 
         geom_point(data=indivsam[indivsam$COV==cov.query,],aes(x=x.tr,y=y), pch=1, size=4))
}

dev.off()


# --- end ---

