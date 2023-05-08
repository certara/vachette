# ------------------------------------------------------------------------- 
#  Sponsor           : Merck 
#  Project Number    : MRK-FTE-VRSV-640
#  Project Root Path : Local
# ------------------------------------------------------------------------- 
#  Program : vachette-example-sigmoid-v30.R
#  Author  : Jos Lommerse - Certara 
#  Date    : 06 February 2022
#  Purpose : mrgsolve model codes
# ------------------------------------------------------------------------- 
#  Software : R version 4.1.2 (2021-11-01)
#  Platform : ThinkPad P15s Gen 2 Model 20W6008AUS
#  Environment : Windows-10, 11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz
# -------------------------------------------------------------------------

##########################################################
#                                                        #
#                    SIGMOID/EMAX                        #
#              DEFINE REFERENCE AND QUERY                #
#                                                        #
##########################################################

sim.sigmoid <- function(nsim.indiv,iiv,ruv,nvpc,SAVE=F,PROP=T)
{
  simtag   <- "sigmoid-v30" # script version number
  
  # ----------- Model settings ------------------
  
  set.seed(20230206)
  
  obs.bmx     <- c(0.003,0.01,0.03,0.1,0.3,1,3,10,30,100,300)
  
  # ---------------- SIMULATED TYPICAL/POPULATION CURVES -------------
  
  xstop         <- 10**2.5  # Samples cannot be located "after" xstop
  xlast.user    <- 10**5    # User provided x value to max. simulate out - replaces expand.factor
  
  step.x.factor <- 1.5
  mrgdelta      <- 0.1
  
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
  typ.data  <- expand_grid(bmx = 10**(c((-10*log10(xlast.user)):(10*log10(xlast.user)))/10),
                             WT = c(wt1,wt2)) %>% 
    mutate(resp = ifelse(WT==wt1, emax.wt1*bmx**gamma.wt1/(bmx**gamma.wt1 + ic50.wt1**gamma.wt1),
                         ifelse(WT==wt2, emax.wt2*bmx**gamma.wt2/(bmx**gamma.wt2 + ic50.wt2**gamma.wt2),NA))) %>% 
    mutate(ID=1, x=bmx, y=resp) %>%            # Standardize
    mutate(PRED=resp) %>% 
    mutate(vachette.cov1 = WT) %>% 
    # Dose number
    mutate(dosenr = 1)    
  
  # --------- Observations ---------------
  
  # simulate ipred
  nid <- nsim.indiv
  
  # Add variable obs.bmx (+/- 5% proportional)
  obs.bmx.var <- NULL
  obs.iiv.var <- NULL
  for(iid in c(1:nid)) 
  {
    add.bmx.var <- t(data.frame(bmx = c(obs.bmx + obs.bmx*runif(length(obs.bmx),-0.05,0.05))))
    obs.bmx.var <- rbind(obs.bmx.var,add.bmx.var)
    
    add.iiv.var <- data.frame(ID = iid,
                              ic50.wt1  = ic50.wt1*rnorm(1,1,sd.ic50),
                              ic50.wt2  = ic50.wt2*rnorm(1,1,sd.ic50),
                              emax.wt1  = emax.wt1*rnorm(1,1,sd.emax),
                              emax.wt2  = emax.wt2*rnorm(1,1,sd.emax),
                              gamma.wt1 = gamma.wt1*rnorm(1,1,sd.gamma),
                              gamma.wt2 = gamma.wt2*rnorm(1,1,sd.gamma))
    obs.iiv.var <- rbind(obs.iiv.var,add.iiv.var)
  }
  
  # Could be more efficient, here follows same "structure" as other models 
  imy.ipred <- NULL
  imy.pred  <- NULL
  for(iid in c(1:nid))
  {
    # individual
    eic50.wt1   <- obs.iiv.var$ic50.wt1[iid]
    eic50.wt2   <- obs.iiv.var$ic50.wt2[iid]
    eemax.wt1   <- obs.iiv.var$emax.wt1[iid]
    eemax.wt2   <- obs.iiv.var$emax.wt2[iid]
    egamma.wt1  <- obs.iiv.var$gamma.wt1[iid]
    egamma.wt2  <- obs.iiv.var$gamma.wt2[iid]
    imy.ipred   <- rbind(imy.ipred, expand_grid(ID=iid,
                                              bmx = obs.bmx.var[iid,],
                                              WT = c(wt1,wt2)) %>% 
                         mutate(resp = ifelse(       WT==wt1, eemax.wt1*bmx**egamma.wt1/(bmx**egamma.wt1 + eic50.wt1**egamma.wt1),
                                                     ifelse(WT==wt2, eemax.wt2*bmx**egamma.wt2/(bmx**egamma.wt2 + eic50.wt2**egamma.wt2),NA))) %>% 
                         # Add residual variability
                         mutate(resp_ruv = resp*rnorm(1,1,PROPVAL) + rnorm(1,0,ADDVAL)))
    
    # Typical
    imy.pred <- rbind(imy.pred, expand_grid(ID=iid,
                                            bmx = obs.bmx.var[iid,],
                                            WT = c(wt1,wt2)) %>% 
                        mutate(resp = ifelse(       WT==wt1, emax.wt1*bmx**gamma.wt1/(bmx**gamma.wt1 + ic50.wt1**gamma.wt1),
                                                    ifelse(WT==wt2, emax.wt2*bmx**gamma.wt2/(bmx**gamma.wt2 + ic50.wt2**gamma.wt2),NA))))
  }
  
  indivsam.all <- as.data.frame(imy.ipred) %>%
    # Add PRED, keep required variables only
    left_join(as.data.frame(imy.pred)[,c('ID','bmx','resp','WT')] %>% rename(PRED=resp),
              by=c('ID','bmx','WT')) %>% 
    mutate(x=bmx, y=resp) %>%             # Standardize
    mutate(vachette.cov1 = WT) %>% 
    mutate(IPRED = resp) %>% 
    mutate(OBS = resp_ruv) %>% 
    # Dose number
    mutate(dosenr = 1) %>% 
    dplyr::select(ID,x,PRED,IPRED,OBS,vachette.cov1,dosenr)
  
  
  # VPC, simulate nvpc times and pick same ID/timepoint combinations
  if(nvpc>0)
  {
    sim.data <- NULL
    for(ivpc in c(1:nvpc))
    {
      imy.ipred <- NULL
      imy.pred  <- NULL
      # Same bmx (independent variable), residual variability re-sampled
      for(iid in c(1:nid))
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
                                                  bmx = obs.bmx.var[iid,],
                                                  WT = c(wt1,wt2)) %>% 
                             mutate(resp = ifelse(       WT==wt1, eemax.wt1*bmx**egamma.wt1/(bmx**egamma.wt1 + eic50.wt1**egamma.wt1),
                                                         ifelse(WT==wt2, eemax.wt2*bmx**egamma.wt2/(bmx**egamma.wt2 + eic50.wt2**egamma.wt2),NA))) %>% 
                             # Add residual variability
                             mutate(resp_ruv = resp*rnorm(1,1,PROPVAL) + rnorm(1,0,ADDVAL)))
        
        # Typical at same bmx
        imy.pred <- rbind(imy.pred, expand_grid(ID=iid,
                                                bmx = obs.bmx.var[iid,],
                                                WT = c(wt1,wt2)) %>% 
                            mutate(resp = ifelse(       WT==wt1, emax.wt1*bmx**gamma.wt1/(bmx**gamma.wt1 + ic50.wt1**gamma.wt1),
                                                        ifelse(WT==wt2, emax.wt2*bmx**gamma.wt2/(bmx**gamma.wt2 + ic50.wt2**gamma.wt2),NA))))
      }
    }
    
    
    indivsam.ivpc <- as.data.frame(imy.ipred) %>%
      # Add PRED, keep required variables only
      left_join(as.data.frame(imy.pred)[,c('ID','bmx','resp','WT')] %>% rename(PRED=resp),
                by=c('ID','bmx','WT')) %>% 
      mutate(x=bmx, y=resp) %>%             # Standardize
      mutate(vachette.cov1 = WT) %>% 
      mutate(IPRED = resp) %>% 
      mutate(OBS = resp_ruv) %>% 
      # Dose number
      mutate(dosenr = 1) %>% 
      dplyr::select(ID,x,PRED,IPRED,OBS,vachette.cov1,dosenr) %>% 
      mutate(isim=ivpc)
    
    # Collect but keep required variables only
    sim.data <- rbindlist(list(sim.data, 
                                   indivsam.ivpc %>% dplyr::select(isim,ID,x,PRED,IPRED,OBS,vachette.cov1,dosenr)))
  }
  
  if(SAVE)
  {
    filetyp <- paste0("../flat-files/vachette-example-", simtag, "-typ.csv")
    fileobs <- paste0("../flat-files/vachette-example-", simtag, "-obs.csv")
    filevpc <- paste0("../flat-files/vachette-example-", simtag, "-vpc.csv")
    if(iiv==0 & ruv==0 & nvpc==0 & !file.exists(filetyp))    write.csv(typ.data,file=filetyp,row.names=F)
    if((iiv!=0 | ruv!=0) & nvpc==0 & !file.exists(fileobs))  write.csv(indivsam.all,file=fileobs,row.names=F)
    if((iiv!=0 | ruv!=0) & nvpc>0 & !file.exists(filevpc))   write.csv(sim.data,file=filevpc,row.names=F)
  }
  
  if(iiv==0 & ruv==0 & nvpc==0)   print(dim(typ.data))
  if((iiv!=0 | ruv!=0) & nvpc==0) print(dim(indivsam.all))
  if((iiv!=0 | ruv!=0) & nvpc>0)  print(dim(sim.data))
  
  if(iiv==0 & ruv==0 & nvpc==0)   return(typ.data)
  if((iiv!=0 | ruv!=0) & nvpc==0) return(indivsam.all)
  if((iiv!=0 | ruv!=0) & nvpc>0)  return(sim.data)
}

