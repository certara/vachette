##########################################################
#                                                        #
#          INDIRECT RESPONSE HIGH DOSE EXAMPLE           #
#                                                        #
##########################################################

#clear environment
rm(list = ls())

# ----------- To be added: ----------------
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyverse)

# Location of this script
home   <- dirname(rstudioapi::getActiveDocumentContext()$path)
getwd()
setwd(home)

filedir <- "flat-files"
ifelse(!dir.exists(filedir), dir.create(filedir), FALSE)

# Mrgsove models
source('vachette-models.R')

#set simulation seed
sim.seed <- 20230205

#define number of subjects to be simulated per replicate
nsim.indiv <- 100

#Note: Set to 1 below if you want to simulate individual subjects
# Observations
iiv        <- 1   # typ=0, obs=1
ruv        <- 1  # typ=0, obs=1

#Note: Define number of vpc replicates if vpc should be done, otherwise set to zero
# VPC
nvpc       <- 100   # typ=0, obs=0, vpc=100

# Other
SAVE       <- T #defines if results should be saved
PROP       <- T #defines if a proportional error is used

simtag   <- "indirect-response-high-dose" # script version number
my.model <- mcode("indirect.response", indirect.response)

# ----------- Model settings ------------------

set.seed(sim.seed)
obs.times     <- c(0.5,1,2,4,8,12,16,24,36,48)

# ---------------- SIMULATED TYPICAL/POPULATION CURVES -------------

xend.typ     <- 1440   # User provided x value to max. simulate out 

mrgdelta      <- 1

dose1 <- 1000000
dose2 <- 100

# (Low) Variability
# IIV
ECL   <- 0.0
EV    <- 0.0
# Low variability
EEC50 <- 0.01
EKIN  <- 0.01
EBSLN <- 0.01
# Residual
PROPVAL<- 0
ADDVAL <- 0
if(PROP)  PROPVAL<- 0.05
if(!PROP) ADDVAL <- 0.0

# ----------- Typical curves --------------

idata <- expand.idata(ID=1, DOSE = c(dose1,dose2))
doses <- expand.ev(amt = c(dose1,dose2), ii = 24, addl=0, cmt = 1)

# simulate pred 
long.out <- my.model %>%
  omat(IIV = dmat(0, 0, 0, 0, 0)) %>%
  smat(SGMA = dmat(0, 0)) %>%
  data_set(doses) %>%
  idata_set(idata) %>%
  carry_out(DOSE, EVID) %>%
  mrgsim(delta=mrgdelta,end=xend.typ)

output.typ  <- as.data.frame(long.out) %>%
  #mutate(ID=1) %>%
  filter(!(time==0 & CENT==0)) %>%      # No Dose rec
  mutate(x=time, y=RESP) %>%            # Standardize
  mutate(PRED=RESP) %>%
  mutate(vachette.cov1 = DOSE) %>%
  # Dose number
  mutate(dosenr = 1,
         ID = ID + 1000000) %>%
  dplyr::select(ID,time,EVID,PRED,DOSE)

# --------- Observations ---------------

# individual data
idata <- expand.idata(ID = c(1:nsim.indiv), DOSE = c(dose1,dose2))

doses <- data.frame(ID = idata$ID, amt = idata$DOSE, ii = 24, addl=0, cmt = 1, evid = 1, time = 0)

# simulate ipred
nid <- length(unique(idata$ID))

# Add variable obs.times
obs.times.var <- NULL
for(i in c(1:nid))
{
  add.times.var <- data.frame(time = c(obs.times + runif(length(obs.times),-0.2,0.2))) %>%
    mutate(ID = i)
  obs.times.var <- rbind(obs.times.var,add.times.var)
}

simdata <- obs.times.var %>%
  mutate(amt = 0, ii = 0, addl = 0, cmt = 2, evid = 0) %>%
  bind_rows(doses) %>%
  arrange(ID, time) %>%
  left_join(idata) 

simdata2 <- simdata %>%
  filter(evid == 1)
  
  # simulate
  imy.ipred <- as.data.frame(my.model %>%
                               omat(IIV = dmat(ECL, EV, EEC50, EKIN, EBSLN)) %>%
                               smat(SGMA = dmat(PROPVAL, ADDVAL)) %>%
                               #ev(amt=dose1,addl=0,ii=24,cmt=1) %>%
                               data_set(simdata) %>%
                               carry_out(DOSE, EVID) %>%
                               mrgsim())
                               #mrgsim(rtol = 1E-12))
  imy.pred <- as.data.frame(my.model %>%
                              omat(IIV = dmat(0, 0, 0, 0, 0)) %>%
                              smat(SGMA = dmat(0, 0)) %>%
                              #ev(amt=dose1,addl=0,ii=24,cmt=1) %>%
                              data_set(simdata) %>%
                              carry_out(DOSE, EVID) %>%
                              mrgsim())
                              #mrgsim(rtol = 1E-12)) 

imy.ipred.df <- as.data.frame(imy.ipred) %>%   mutate(vachette.cov1 = DOSE) 
imy.pred.df  <- as.data.frame(imy.pred)  %>%   mutate(vachette.cov1 = DOSE) 

indivsam.all <- imy.ipred.df %>%
  # Add PRED, keep required variables only
  left_join(imy.pred.df %>% dplyr::select(ID,time,RESP,vachette.cov1) %>% rename(PRED=RESP),
            by=c('ID','time','vachette.cov1')) %>%
  #filter(!(time==0)) %>%     # No Dose rec
  mutate(x=time, y=RESP) %>%             # Standardize
  mutate(IPRED = RESP) %>%
  mutate(DV = RESP_RUV) %>%
  # Dose number
  mutate(dosenr = 1) %>%
  dplyr::select(ID,time,EVID,DV,PRED,IPRED,DOSE)


# VPC, simulate nvpc times and pick same ID/timepoint combinations
if(nvpc>0)
{
  indivsam.vpc <- NULL
  for(ivpc in c(1:nvpc))
  {
    
      # high dose
    imy.ipred <- as.data.frame(my.model %>%
                                 omat(IIV = dmat(ECL, EV, EEC50, EKIN, EBSLN)) %>%
                                 smat(SGMA = dmat(PROPVAL, ADDVAL)) %>%
                                 #ev(amt=dose1,addl=0,ii=24,cmt=1) %>%
                                 data_set(simdata) %>%
                                 carry_out(DOSE, EVID) %>%
                                 mrgsim())
                                 #mrgsim(rtol = 1E-12))
    imy.pred <- as.data.frame(my.model %>%
                                omat(IIV = dmat(0, 0, 0, 0, 0)) %>%
                                smat(SGMA = dmat(0, 0)) %>%
                                #ev(amt=dose1,addl=0,ii=24,cmt=1) %>%
                                data_set(simdata) %>%
                                carry_out(DOSE, EVID) %>%
                                mrgsim())
                                #mrgsim(rtol = 1E-12))
    

    imy.ipred.df <- as.data.frame(imy.ipred) %>%   mutate(vachette.cov1 = DOSE) 
    imy.pred.df  <- as.data.frame(imy.pred)  %>%   mutate(vachette.cov1 = DOSE) 
    
    indivsam.ivpc <- imy.ipred.df %>%
      # Add PRED, keep required variables only
      left_join(imy.pred.df %>% dplyr::select(ID,time,RESP,vachette.cov1) %>% rename(PRED=RESP),
                by=c('ID','time','vachette.cov1')) %>%
      #filter(!(time==0)) %>%     # No Dose rec
      mutate(x=time, y=RESP) %>%             # Standardize
      mutate(IPRED = RESP) %>%
      mutate(DV = RESP_RUV) %>%
      # Dose number
      mutate(dosenr = 1) %>%
      dplyr::select(ID,time,EVID,PRED,IPRED,DV,DOSE,dosenr) %>%
      mutate(REP=ivpc)

    # Collect but keep required variables only
    indivsam.vpc <- rbindlist(list(indivsam.vpc,
                                   indivsam.ivpc %>% dplyr::select(REP,ID,time,EVID,DV,PRED,IPRED,DOSE)))
  }
}

if(SAVE)
{
  filetyp <- paste0(filedir, "/", simtag, "-typ.csv")
  fileobs <- paste0(filedir, "/", simtag, "-obs.csv")
  filevpc <- paste0(filedir, "/", simtag, "-vpc.csv")
  if(file.exists(filetyp) | file.exists(fileobs) | file.exists(filevpc)) {
    warning("At least one of the files you are attempting to write already exists and was not written, check file folder")
  }
  if(exists('output.typ') & !file.exists(filetyp))    write.csv(output.typ,file=filetyp,row.names=F)
  if(exists('indivsam.all') & !file.exists(fileobs))  write.csv(indivsam.all,file=fileobs,row.names=F)
  if(exists('indivsam.vpc') & !file.exists(filevpc))   write.csv(indivsam.vpc,file=filevpc,row.names=F)
}

if(exists('output.typ'))   print(dim(output.typ))
if(exists('indivsam.all')) print(dim(indivsam.all))
if(exists('indivsam.vpc'))  print(dim(indivsam.vpc))


# ---- end ----

