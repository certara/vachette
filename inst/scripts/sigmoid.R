##########################################################
#                                                        #
#             SIGMOID/EMAX, SPARSE SAMPLES               #
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
#source('vachette-models.R')

#set simulation seed
sim.seed <- 20230331

#define number of subjects to be simulated per replicate
nsim.indiv <- 100
nsim.samples <- 6

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
#sim.sigmoid <- function(nsim.indiv,iiv,ruv,nvpc,SAVE=F,PROP=T)
#{
  simtag   <- "sigmoid" # script version number
  
  # ----------- Model settings ------------------
  
  set.seed(20230206)
  
  obs.bmx1     <- c(0.003,0.01,0.1,1,10,100)
  obs.bmx2     <- c(0.01,0.03,0.3,3,30,300)
  xend.typ      <- 10**5  # User provided x value to max. simulate out 
  
  # ---------------- SIMULATED TYPICAL/POPULATION CURVES -------------
  
  
  wt1 <- 70
  wt2 <- 30
  
  # Typical
  ic50.wt1    <- 1
  ic50.wt2    <- 10
  
  emax.wt1     <- 10
  emax.wt2     <- 10
  
  gamma.wt1   <- 1
  gamma.wt2   <- 1
  
  # Variability (as stand.dev.)
  # IIV
  sd.ic50   <- 0.1
  sd.emax   <- 0.1
  sd.gamma  <- 0.0
  # Residual
  PROPVAL<- 0
  ADDVAL <- 0
  if(PROP)  PROPVAL<- 0.05  # As s.d.
  if(!PROP) ADDVAL <- 0.01  # As s.d.
  
  # ----------- Typical curves --------------
  
  # simulate pred
  output.typ  <- expand_grid(bmx = 10**(c((-10*log10(xend.typ)):(10*log10(xend.typ)))/10),
                             WT = c(wt1,wt2)) %>% 
    mutate(resp = ifelse(WT==wt1, emax.wt1*bmx**gamma.wt1/(bmx**gamma.wt1 + ic50.wt1**gamma.wt1),
                         ifelse(WT==wt2, emax.wt2*bmx**gamma.wt2/(bmx**gamma.wt2 + ic50.wt2**gamma.wt2),NA))) %>% 
    mutate(ID=1, x=bmx, y=resp) %>%            # Standardize
    mutate(PRED=resp) %>% 
    mutate(vachette.cov1 = WT) %>% 
    # Dose number
    mutate(dosenr = 1,
           ID = ifelse(WT == 30, 1000001, 1000002)) %>% 
    dplyr::select(ID,bmx,PRED,WT) %>%
    arrange(ID, bmx)
  
  # --------- Observations ---------------
  
  # simulate ipred
  nid <- nsim.indiv
										
  
  # Add variable obs.bmx (+/- 30% proportional)
  obs.bmx.var <- NULL
  obs.iiv.var <- NULL
  for (iweight in c(1,2)) {
  
  for(iid in c(1:nid)) 
  {
    tmp.samples <- sort(sample(obs.bmx1, nsim.samples))
    if (iweight == 2) {
      tmp.samples <- sort(sample(obs.bmx2, nsim.samples))
    }
    add.bmx.var <- t(data.frame(bmx = c(tmp.samples + tmp.samples*runif(length(tmp.samples),-0.30,0.30))))
    obs.bmx.var <- rbind(obs.bmx.var,add.bmx.var)
    
    add.iiv.var <- data.frame(ID = ifelse(iweight == 1, iid, iid + nsim.indiv),
                              WT = ifelse(iweight == 1, wt1, wt2),
                              ic50.wt1  = ic50.wt1*rnorm(1,1,sd.ic50),
                              ic50.wt2  = ic50.wt2*rnorm(1,1,sd.ic50),
                              emax.wt1  = emax.wt1*rnorm(1,1,sd.emax),
                              emax.wt2  = emax.wt2*rnorm(1,1,sd.emax),
                              gamma.wt1 = gamma.wt1*rnorm(1,1,sd.gamma),
                              gamma.wt2 = gamma.wt2*rnorm(1,1,sd.gamma))
    obs.iiv.var <- rbind(obs.iiv.var,add.iiv.var)
  }
  }
  
  # Could be more efficient, here follows same "structure" as other models 
  imy.ipred <- NULL
  imy.pred  <- NULL
  for(iid in unique(obs.iiv.var$ID))
  {
    # individual
    eic50.wt1   <- obs.iiv.var$ic50.wt1[iid]
    eic50.wt2   <- obs.iiv.var$ic50.wt2[iid]
    eemax.wt1   <- obs.iiv.var$emax.wt1[iid]
    eemax.wt2   <- obs.iiv.var$emax.wt2[iid]
    egamma.wt1  <- obs.iiv.var$gamma.wt1[iid]
    egamma.wt2  <- obs.iiv.var$gamma.wt2[iid]
    imy.ipred   <- rbind(imy.ipred, expand_grid(ID=iid,
                                              bmx = obs.bmx.var[iid,]) %>% 
                         mutate(WT = obs.iiv.var$WT[iid],
                                resp = ifelse(       WT==wt1, eemax.wt1*bmx**egamma.wt1/(bmx**egamma.wt1 + eic50.wt1**egamma.wt1),
                                                     ifelse(WT==wt2, eemax.wt2*bmx**egamma.wt2/(bmx**egamma.wt2 + eic50.wt2**egamma.wt2),NA))) %>% 
                         # Add residual variability
                         mutate(resp_ruv = resp*rnorm(1,1,PROPVAL) + rnorm(1,0,ADDVAL)))
    
    # Typical
    imy.pred <- rbind(imy.pred, expand_grid(ID=iid,
                                            bmx = obs.bmx.var[iid,]) %>% 
                        mutate(WT = obs.iiv.var$WT[iid],
                               resp = ifelse(       WT==wt1, emax.wt1*bmx**gamma.wt1/(bmx**gamma.wt1 + ic50.wt1**gamma.wt1),
                                                    ifelse(WT==wt2, emax.wt2*bmx**gamma.wt2/(bmx**gamma.wt2 + ic50.wt2**gamma.wt2),NA))))
  }
  
  indivsam.all <- as.data.frame(imy.ipred) %>%
    # Add PRED, keep required variables only
    left_join(as.data.frame(imy.pred)[,c('ID','bmx','resp','WT')] %>% rename(PRED=resp),
              by=c('ID','bmx','WT')) %>% 
    mutate(x=bmx, y=resp) %>%             # Standardize
    mutate(vachette.cov1 = WT) %>% 
    mutate(IPRED = resp) %>% 
    mutate(DV = resp_ruv) %>% 
    # Dose number
    mutate(dosenr = 1) %>% 
    dplyr::select(ID,bmx,DV,PRED,WT)
  
  
  # VPC, simulate nvpc times and pick same ID/timepoint combinations
  if(nvpc>0)
  {
    indivsam.vpc <- NULL
    for(ivpc in c(1:nvpc))
    {
      imy.ipred <- NULL
      imy.pred  <- NULL
      # Same bmx (independent variable), residual variability re-sampled
      for(iid in unique(obs.iiv.var$ID))
      {
        # individual 
        # eic50.wt1   <- obs.iiv.var$ic50.wt1[i]
        # eic50.wt2   <- obs.iiv.var$ic50.wt2[i]
        # eemax.wt1   <- obs.iiv.var$emax.wt1[i]
        # eemax.wt2   <- obs.iiv.var$emax.wt2[i]
        # egamma.wt1 <- obs.iiv.var$gamma.wt1[i]
        # egamma.wt2 <- obs.iiv.var$gamma.wt2[i]
        eic50.wt1  = ic50.wt1*rnorm(1,1,sd.ic50)
        eic50.wt2  = ic50.wt2*rnorm(1,1,sd.ic50)
        eemax.wt1  = emax.wt1*rnorm(1,1,sd.emax)
        eemax.wt2  = emax.wt2*rnorm(1,1,sd.emax)
        egamma.wt1 = gamma.wt1*rnorm(1,1,sd.gamma)
        egamma.wt2 = gamma.wt2*rnorm(1,1,sd.gamma)
        
        imy.ipred <- rbind(imy.ipred, expand_grid(ID=iid,
                                                  bmx = obs.bmx.var[iid,]) %>% 
                             mutate(WT = obs.iiv.var$WT[iid],
                                    resp = ifelse(       WT==wt1, eemax.wt1*bmx**egamma.wt1/(bmx**egamma.wt1 + eic50.wt1**egamma.wt1),
                                                         ifelse(WT==wt2, eemax.wt2*bmx**egamma.wt2/(bmx**egamma.wt2 + eic50.wt2**egamma.wt2),NA))) %>% 
                             # Add residual variability
                             mutate(resp_ruv = resp*rnorm(1,1,PROPVAL) + rnorm(1,0,ADDVAL)))
        
        # Typical at same bmx
        imy.pred <- rbind(imy.pred, expand_grid(ID=iid,
                                                bmx = obs.bmx.var[iid,]) %>% 
                            mutate(WT = obs.iiv.var$WT[iid],
                                   resp = ifelse(       WT==wt1, emax.wt1*bmx**gamma.wt1/(bmx**gamma.wt1 + ic50.wt1**gamma.wt1),
                                                        ifelse(WT==wt2, emax.wt2*bmx**gamma.wt2/(bmx**gamma.wt2 + ic50.wt2**gamma.wt2),NA))))
      }
    
    
    
    indivsam.ivpc <- as.data.frame(imy.ipred) %>%
      # Add PRED, keep required variables only
      left_join(as.data.frame(imy.pred)[,c('ID','bmx','resp','WT')] %>% rename(PRED=resp),
                by=c('ID','bmx','WT')) %>% 
      mutate(x=bmx, y=resp) %>%             # Standardize
      mutate(vachette.cov1 = WT) %>% 
      mutate(IPRED = resp) %>% 
      mutate(DV = resp_ruv) %>% 
      # Dose number
      mutate(dosenr = 1) %>% 
      dplyr::select(ID,bmx,PRED,IPRED,DV,WT,dosenr) %>% 
      mutate(REP=ivpc)
    
    # Collect but keep required variables only
    indivsam.vpc <- rbindlist(list(indivsam.vpc, 
                                   indivsam.ivpc %>% dplyr::select(REP,ID,bmx,DV,PRED,WT)))
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
  
#}

