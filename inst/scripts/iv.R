##########################################################
#                                                        #
#                       IV                               #
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

home   <- dirname(rstudioapi::getActiveDocumentContext()$path)
getwd()
setwd(home)

filedir <- "flat-files"
ifelse(!dir.exists(filedir), dir.create(filedir), FALSE)

# Mrgsove models
source('vachette-models.R')

#set simulation seed
sim.seed <- 123456

#define number of subjects to be simulated per replicate
nsim.indiv <- 50

#Note: The typical profile is always simulated, so no need to simulate separately if iiv or ruv greater than zero
# Typical
# iiv        <- 0
# ruv        <- 0

#Note: Set to 1 below if you want to simulate individual subjects
# Observations
iiv        <- 1
ruv        <- 1

#Note: Define number of vpc replicates if vpc should be done, otherwise set to zero
# VPC
nvpc       <- 100
# Other
SAVE       <- T  #defines if results should be saved
PROP       <- T  #defines if a proportional error is used		  

  my.model <- mcode("iv", iv)
  simtag   <- "iv" # script version number
  
  # ----------- Model settings ------------------
  
  obs.times     <- c(0.25,0.5,1,2,3,4,8,12)
  
  # ---------------- SIMULATED TYPICAL/POPULATION CURVES -------------
  
  xend.typ    <- 108  # User provided x value to max. simulate out - replaces expand.factor

  mrgdelta      <- 0.1
  
  wt1 <- 30
  wt2 <- 70

  # Variability
  # IIV
  ECL <- 0.2
  EV  <- 0.2
  # Residual
  PROPVAL<- 0
  ADDVAL <- 0
  if(PROP)  PROPVAL<- 0.1
  if(!PROP) ADDVAL <- 0.01
  
  # ---------------- SIMULATED TYPICAL/POPULATION CURVES -------------
  
  # population input data
  data <- expand.idata(ID=1,WT=c(wt1,wt2))
  
  # simulate pred
  long.out <- my.model %>%
    omat(IIV = dmat(0, 0)) %>% 
    smat(SGMA = dmat(0, 0)) %>%
    ev(amt=100,addl=0,ii=24,cmt=1) %>%
    idata_set(data) %>%
    carry_out(EVID) %>%
    mrgsim(delta=mrgdelta,end=xend.typ) 
  
  output.typ  <- as.data.frame(long.out) %>%
    filter(!(time==0 & CENT==0)) %>%     # Only Dose rec at time=0
    mutate(x=time, y=CP) %>%            # Standardize
    mutate(PRED=CP) %>% 
    mutate(vachette.cov1 = WT) %>% 
    # Dose number
    mutate(dosenr = 1,
           ID = ID + 1000000) %>% 
    dplyr::select(ID,time,EVID,PRED,WT)
  
  # --------- Observations ---------------
  
  # individual data
  idata <- expand.idata(ID = c(1:nsim.indiv),
                        WT = c(wt1,wt2))
  
  nid <- length(unique(idata$ID))
  
  # Add variable obs.times
  obs.times.var <- NULL
  for(i in c(1:nid)) 
  {
    add.times.var <- t(data.frame(times = c(obs.times + runif(length(obs.times),-0.2,0.2))))
    obs.times.var <- rbind(obs.times.var,add.times.var)
  }
  
  # simulate ipred
  set.seed(sim.seed)
  
  imy.ipred <- NULL
  imy.pred  <- NULL
  for(iid in c(1:nid))
  {
    imy.ipred <- rbind(imy.ipred, as.data.frame(my.model %>%
                                                  omat(IIV = dmat(ECL, EV)) %>%    
                                                  smat(SGMA = dmat(PROPVAL, ADDVAL)) %>%    
                                                  ev(amt=100,addl=0,ii=24,cmt=1) %>%
                                                  idata_set(idata[iid,]) %>%
                                                  carry_out(EVID) %>%
                                                  mrgsim(tgrid=obs.times.var[iid,])))
    imy.pred <- rbind(imy.pred, as.data.frame(my.model %>%
                                                omat(IIV = dmat(0, 0)) %>% 
                                                smat(SGMA = dmat(0, 0)) %>% 
                                                ev(amt=100,addl=0,ii=24,cmt=1) %>%
                                                idata_set(idata[iid,]) %>%
                                                carry_out(EVID) %>%
                                                mrgsim(tgrid=obs.times.var[iid,])))
  }
  
  indivsam.all <- as.data.frame(imy.ipred) %>%
    # Add PRED, keep required variables only
    left_join(as.data.frame(imy.pred)[,c('ID','time','CP','WT')] %>% rename(PRED=CP),
              by=c('ID','time','WT')) %>% 
    #filter(!(time==0)) %>%     # No Dose rec
    mutate(x=time, y=CP) %>%             # Standardize
    mutate(vachette.cov1 = WT) %>% 
    mutate(IPRED = CP) %>% 
    mutate(DV = CP_RUV) %>% 
    # Dose number
    mutate(dosenr = 1) %>% 
    dplyr::select(ID,time,EVID,DV,PRED,WT)
  
  # VPC, simulate nvpc times and pick same ID/timepoint combinations
  if(nvpc>0)
  {
    indivsam.vpc <- NULL
    for(ivpc in c(1:nvpc))
    {
      imy.ipred <- NULL
      imy.pred  <- NULL
      for(iid in c(1:nid))
      {
        imy.ipred <- rbind(imy.ipred, as.data.frame(my.model %>%
                                                      omat(IIV = dmat(ECL, EV)) %>%    # ECL, EV
                                                      smat(SGMA = dmat(PROPVAL, ADDVAL)) %>%     # PROP, ADD 
                                                      ev(amt=100,addl=0,ii=24,cmt=1) %>%
                                                      idata_set(idata[iid,]) %>%
                                                      carry_out(EVID) %>%
                                                      mrgsim(tgrid=obs.times.var[iid,])))
        
        imy.pred <- rbind(imy.pred, as.data.frame(my.model %>%
                                                    omat(IIV = dmat(0, 0)) %>% 
                                                    smat(SGMA = dmat(0, 0)) %>% 
                                                    ev(amt=100,addl=0,ii=24,cmt=1) %>%
                                                    idata_set(idata[iid,]) %>%
                                                    carry_out(EVID) %>%
                                                    mrgsim(tgrid=obs.times.var[iid,])))
      }

      
      indivsam.ivpc <- as.data.frame(imy.ipred) %>%
        # Add PRED, keep required variables only
        left_join(as.data.frame(imy.pred)[,c('ID','time','CP','WT')] %>% rename(PRED=CP),
                  by=c('ID','time','WT')) %>% 
        #filter(!(time==0)) %>%     # No Dose rec
        mutate(x=time, y=CP) %>%             # Standardize
        mutate(vachette.cov1 = WT) %>% 
        mutate(IPRED = CP) %>% 
        mutate(DV = CP_RUV) %>% 
        # Dose number
        mutate(dosenr = 1) %>% 
        dplyr::select(ID,time,EVID,PRED,IPRED,DV,WT,dosenr) %>% 
        mutate(REP=ivpc)
      
      # Collect but keep required variables only
      indivsam.vpc <- rbindlist(list(indivsam.vpc, 
                            indivsam.ivpc %>% dplyr::select(REP,ID,time,EVID,DV,PRED,WT)))
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
   

