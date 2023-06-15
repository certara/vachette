##########################################################
#                                                        #
#                    PEMBROLIZUMAB MODEL                 #
#                                                        #
##########################################################
#clear environment
rm(list = ls())

# ----------- To be added: ----------------
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(data.table)

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

  my.model <- mcode("pembro3", pembro3)
  simtag   <- "pembro" # script version number
  
  # Note: Time=0 and trough already give max conc, so new dosenr start there
  
  # ----------- Model settings ------------------
  
  set.seed(sim.seed)
  
  obs.times.q3w     <- c(0.1,   0.2,  0.4,  1,  2,  4, 10, 21, 
                         21.1, 21.2, 21.4, 22, 23, 25, 31, 42) # Days, Q3W
  obs.times.q2w     <- c(0.1,   0.2 , 0.4,  1,  2,  5, 14,
                         14.1, 14.2, 14.4, 15, 16, 19, 28,
                         28.1, 28.2, 28.4, 29, 30, 33, 42) # Days, Q2W
  
  xend.typ     <- 800  # User provided x value to max. simulate out 

  mrgdelta      <- 0.1
  
  sched1 <- 'Q2W'
  sched2 <- 'Q3W'
  alb1   <- 16
  alb2   <- 53.5
  
  # Variability
  # IIV
  etaCLQ <- 0.067
  etaV   <- 0.0208
  # Ignoring correlation here
  # $OMEGA @annotated @block
  # //original: etaCLQ: 0.1341               : ETA on clearance and Q
  # //original: etaV:   0.0159127  0.0417008 : ETA on volumes
  # Reduced:
  # etaCLQ: 0.067         : ETA on clearance and Q
  # etaV:   0.008  0.0208 : ETA on volumes
  
  # Residual - proportional only
  # RESERR = -0.272147
  RESERR<- 0
  # ADDVAL <- 0
  if(PROP)  RESERR <- -0.272147
  if(!PROP) stop('No additive error in Pembro model implemented')
  
  # ---------------- SIMULATED TYPICAL/POPULATION CURVES -------------
  
  # population input data
  data <- expand.idata(ID=1,SCHED=c(sched1,sched2),ALB=c(alb1,alb2))
  
  # simulate pred
  long.out.q2w <- my.model %>%
    omat(IIV = dmat(0, 0)) %>% 
    # smat(SGMA = dmat(0, 0)) %>%
    ev(amt=100,addl=2,ii=14,cmt=1) %>% # Every 2 wks
    idata_set(data[data$SCHED=='Q2W',]) %>%
    carry_out(EVID, ADDL, II) %>%
    mrgsim(delta=mrgdelta,end=xend.typ) 
  long.out.q3w <- my.model %>%
    omat(IIV = dmat(0, 0)) %>% 
    # smat(SGMA = dmat(0, 0)) %>%
    ev(amt=100,addl=1,ii=21,cmt=1) %>% # Every 3 wks
    idata_set(data[data$SCHED=='Q3W',]) %>%
    carry_out(EVID, ADDL, II) %>%
    mrgsim(delta=mrgdelta,end=xend.typ) 
  
  output.typ  <- rbind(
    # Q2W
    as.data.frame(long.out.q2w) %>% 
      mutate(SCHED='Q2W') %>% 
      # dosenr
      mutate(dosenr = ifelse(time<42,floor(time/14)+1,3)),
    # Q3W
    as.data.frame(long.out.q3w) %>% 
      mutate(SCHED='Q3W') %>% 
    # dosenr
      mutate(dosenr = ifelse(time<42,floor(time/21)+1,2))) %>% 
    filter(!(time==0 & CONC==0)) %>%     # No Dose rec
    mutate(x=time, y=CONC) %>%       # Standardize
    mutate(PRED=CONC,
           ID = ID + 1000000) %>% 
    mutate(vachette.cov1 = SCHED) %>% 
    mutate(vachette.cov2 = ALB) %>%
    dplyr::select(ID, time, EVID, ADDL, II, PRED, SCHED, ALB)
  
  output.typ %>% ggplot(aes(x=x,y=y,group=ID))+
    geom_line()+
    facet_grid(SCHED~ALB)
  
  # --------- Observations / VPC prep ---------------
  
  # individual data
  idata <- expand.idata(ID   = c(1:nsim.indiv),
                        SCHED= c(sched1,sched2),
                        ALB  = c(alb1,alb2))
  
  # Q2W only:
  idata.q2w <- idata %>% filter(SCHED=='Q2W')
  nid.q2w <- length(unique(idata.q2w$ID))
  # Add variable obs.times
  obs.times.var.q2w <- NULL
  for(i in c(1:nid.q2w)) 
  {
    add.times.var.q2w <- t(data.frame(times = c(obs.times.q2w + runif(length(obs.times.q2w),-0.1,0.1))))
    obs.times.var.q2w <- rbind(obs.times.var.q2w,add.times.var.q2w)
  }
  
  # Q3W only:
  idata.q3w <- idata %>% filter(SCHED=='Q3W')
  nid.q3w <- length(unique(idata.q3w$ID))
  # Add variable obs.times
  obs.times.var.q3w <- NULL
  for(i in c(1:nid.q3w)) 
  {
    add.times.var.q3w <- t(data.frame(times = c(obs.times.q3w + runif(length(obs.times.q3w),-0.1,0.1))))
    obs.times.var.q3w <- rbind(obs.times.var.q3w,add.times.var.q3w)
  }
  
  # OBS, simulation Q2W
  imy.ipred.q2w <- NULL
  imy.pred.q2w  <- NULL
  for(iid in c(1:nid.q2w))
  {
    imy.ipred.q2w <- rbind(imy.ipred.q2w, as.data.frame(my.model %>%
                                                  omat(IIV = dmat(etaCLQ, etaV)) %>% 
                                                  smat(SGMA = dmat(1)) %>%  # Just prop err
                                                  ev(amt=100,addl=2,ii=14,cmt=1) %>% # Every 2 wks
                                                  idata_set(idata.q2w[iid,]) %>%
                                                  carry_out(EVID, ADDL, II) %>%
                                                  mrgsim(tgrid=obs.times.var.q2w[iid,])) %>% 
                             mutate(SCHED = 'Q2W') %>% 
                             # dosenr
                             mutate(dosenr = ifelse(time<42,floor(time/14)+1,3)))
    imy.pred.q2w <- rbind(imy.pred.q2w, as.data.frame(my.model %>%
                                                 omat(IIV = dmat(0,0)) %>% 
                                                 smat(SGMA = dmat(0)) %>%  # Just prop err
                                                 ev(amt=100,addl=2,ii=14,cmt=1) %>% # Every 2 wks
                                                 idata_set(idata.q2w[iid,]) %>%
                                                  carry_out(EVID, ADDL, II) %>%
                                                 mrgsim(tgrid=obs.times.var.q2w[iid,])))
  }
  
  # OBS, simulation Q3W
  imy.ipred.q3w <- NULL
  imy.pred.q3w  <- NULL
  for(iid in c(1:nid.q3w))
  {
    imy.ipred.q3w <- rbind(imy.ipred.q3w, as.data.frame(my.model %>%
                                                          omat(IIV = dmat(etaCLQ, etaV)) %>% 
                                                          smat(SGMA = dmat(1)) %>%  # Just prop err
                                                          ev(amt=100,addl=1,ii=21,cmt=1) %>% # Every 3 wks
                                                          idata_set(idata.q3w[iid,]) %>%
                                                          carry_out(EVID, ADDL, II) %>%
                                                          mrgsim(tgrid=obs.times.var.q3w[iid,])) %>% 
                             mutate(SCHED = 'Q3W') %>% 
                             # dosenr
                             mutate(dosenr = ifelse(time<42,floor(time/21)+1,2)))
    
    imy.pred.q3w <- rbind(imy.pred.q3w, as.data.frame(my.model %>%
                                                omat(IIV = dmat(0,0)) %>% 
                                                smat(SGMA = dmat(0)) %>%  # Just prop err
                                                ev(amt=100,addl=1,ii=21,cmt=1) %>% # Every 3 wks
                                                idata_set(idata.q3w[iid,]) %>%
                                                carry_out(EVID, ADDL, II) %>%
                                                mrgsim(tgrid=obs.times.var.q3w[iid,])))
    }
  
  indivsam.all <- rbind(
    imy.ipred.q2w %>%
    # Add PRED, keep required variables only
    left_join(imy.pred.q2w[,c('ID','time','CONC','ALB')] %>% rename(PRED=CONC),
              by=c('ID','time','ALB')),
    imy.ipred.q3w %>%
      # Add PRED, keep required variables only
      left_join(imy.pred.q3w[,c('ID','time','CONC','ALB')] %>% rename(PRED=CONC),
                by=c('ID','time','ALB'))) %>% 
    mutate(x=time, y=CONC) %>%       # Standardize
    mutate(vachette.cov1 = SCHED) %>% 
    mutate(vachette.cov2 = ALB) %>% 
    mutate(IPRED = CONC) %>% 
    mutate(DV = CONC_RUV) %>%
    dplyr::select(ID, time, EVID, ADDL, II, DV, PRED, IPRED, SCHED, ALB)
  
  # Check
  indivsam.all %>% 
    ggplot(aes(x=time))+
    geom_point(aes(y=DV))+
    geom_point(aes(y=IPRED),col='red')+
    geom_point(aes(y=PRED),col='cyan')+
    geom_line(aes(y=PRED))+
    facet_grid(SCHED~ALB)
  
  # VPC, simulate nvpc times and pick same ID/timepoint combinations
  if(nvpc>0)
  {
    indivsam.vpc <- NULL
    for(ivpc in c(1:nvpc))
    {
      
      # OBS, simulation Q2W
      imy.ipred.q2w <- NULL
      imy.pred.q2w  <- NULL
      for(iid in c(1:nid.q2w))
      {
        imy.ipred.q2w <- rbind(imy.ipred.q2w, as.data.frame(my.model %>%
                                                              omat(IIV = dmat(etaCLQ, etaV)) %>% 
                                                              smat(SGMA = dmat(1)) %>%  # Just prop err
                                                              ev(amt=100,addl=2,ii=14,cmt=1) %>% # Every 2 wks
                                                              idata_set(idata.q2w[iid,]) %>%
                                                              carry_out(EVID, ADDL, II) %>%
                                                              mrgsim(tgrid=obs.times.var.q2w[iid,])) %>% 
                                 mutate(SCHED = 'Q2W') %>% 
                                 # dosenr
                                 mutate(dosenr = ifelse(time<42,floor(time/14)+1,3)))
        
        imy.pred.q2w <- rbind(imy.pred.q2w, as.data.frame(my.model %>%
                                                            omat(IIV = dmat(0,0)) %>% 
                                                            smat(SGMA = dmat(0)) %>%  # Just prop err
                                                            ev(amt=100,addl=2,ii=14,cmt=1) %>% # Every 2 wks
                                                            idata_set(idata.q2w[iid,]) %>%
                                                            carry_out(EVID, ADDL, II) %>%
                                                            mrgsim(tgrid=obs.times.var.q2w[iid,])))
      }
      
      # OBS, simulation Q3W
      imy.ipred.q3w <- NULL
      imy.pred.q3w  <- NULL
      for(iid in c(1:nid.q3w))
      {
        imy.ipred.q3w <- rbind(imy.ipred.q3w, as.data.frame(my.model %>%
                                                              omat(IIV = dmat(etaCLQ, etaV)) %>% 
                                                              smat(SGMA = dmat(1)) %>%  # Just prop err
                                                              ev(amt=100,addl=1,ii=21,cmt=1) %>% # Every 3 wks
                                                              idata_set(idata.q3w[iid,]) %>%
                                                              carry_out(EVID, ADDL, II) %>%
                                                              mrgsim(tgrid=obs.times.var.q3w[iid,])) %>% 
                                 mutate(SCHED = 'Q3W') %>% 
                                 # dosenr
                                 mutate(dosenr = ifelse(time<42,floor(time/21)+1,2)))
        
        imy.pred.q3w <- rbind(imy.pred.q3w, as.data.frame(my.model %>%
                                                            omat(IIV = dmat(0,0)) %>% 
                                                            smat(SGMA = dmat(0)) %>%  # Just prop err
                                                            ev(amt=100,addl=1,ii=21,cmt=1) %>% # Every 3 wks
                                                            idata_set(idata.q3w[iid,]) %>%
                                                            carry_out(EVID, ADDL, II) %>%
                                                            mrgsim(tgrid=obs.times.var.q3w[iid,])))
      }
      
      indivsam.ivpc <- rbind(
        imy.ipred.q2w %>%
          # Add PRED, keep required variables only
          left_join(imy.pred.q2w[,c('ID','time','CONC','ALB')] %>% rename(PRED=CONC),
                    by=c('ID','time','ALB')),
        imy.ipred.q3w %>%
          # Add PRED, keep required variables only
          left_join(imy.pred.q3w[,c('ID','time','CONC','ALB')] %>% rename(PRED=CONC),
                    by=c('ID','time','ALB'))) %>% 
        mutate(x=time, y=CONC) %>%       # Standardize
        mutate(vachette.cov1 = SCHED) %>% 
        mutate(vachette.cov2 = ALB) %>% 
        mutate(IPRED = CONC) %>% 
        mutate(DV = CONC_RUV) %>% 
        mutate(REP = ivpc)
      
      # Collect but keep required variables only
      indivsam.vpc <- rbindlist(list(indivsam.vpc, 
                            indivsam.ivpc %>% dplyr::select(REP,ID,time,EVID,ADDL,II,DV,PRED,SCHED,ALB)))
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
  
  
