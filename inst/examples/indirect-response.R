# ------------------------------------------------------------------------- 
#  Sponsor           : Merck 
#  Project Number    : MRK-FTE-VRSV-640
#  Project Root Path : Local
# ------------------------------------------------------------------------- 
#  Program : vachette-indirect-response-v30.R
#  Author  : Jos Lommerse - Certara 
#  Date    : 6 February 2023
#  Purpose : Model simulations using mrgsolve
# ------------------------------------------------------------------------- 
#  Software : R version 4.2.2 (2022-10-31 ucrt)
#  Platform : ThinkPad P15s Gen 2 Model 20W6008AUS
#  Environment : Windows-10, 11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz
# -------------------------------------------------------------------------


##########################################################
#                                                        #
#               INDIRECT RESPONSE EXAMPLE                #
#              DEFINE REFERENCE AND QUERY                #
#                                                        #
##########################################################


sim.indirect.response <- function(nsim.indiv,iiv,ruv,nvpc,SAVE=F,PROP=T)
{

  simtag   <- "indirect-response-v30" # script version number
  
  my.model <- mcode("indirect.response", indirect.response)
  
  # ----------- Model settings ------------------
  
  set.seed(20230206)
  obs.times     <- c(0.5,1,2,4,8,12,16,24,36,48)
  
  # ---------------- SIMULATED TYPICAL/POPULATION CURVES -------------
  
  xstop         <- 48   # Actual curves end at xstop = 48 hr
  xlast.user    <- 96   # User provided x value to max. simulate out - replaces expand.factor

  step.x.factor <- 1.5
  mrgdelta      <- 0.1
  
  wt1 <- 30
  wt2 <- 70
  
  # Variability
  # IIV
  ECL   <- 0.0
  EV    <- 0.0
  EEC50 <- 0.1
  EKIN  <- 0.1
  EBSLN <- 0.1
  # Residual
  PROPVAL<- 0
  ADDVAL <- 0
  if(PROP)  PROPVAL<- 0.1
  if(!PROP) ADDVAL <- 0.01
  
  # ----------- Typical curves --------------
  
  data <- expand.idata(ID=1,WT=c(wt1,wt2))

  # simulate pred
  long.out <- my.model %>%
    omat(IIV = dmat(0, 0, 0, 0, 0)) %>% 
    smat(SGMA = dmat(0, 0)) %>%
    ev(amt=100,addl=0,ii=24,cmt=1) %>%
    idata_set(data) %>%
    mrgsim(delta=mrgdelta,end=step.x.factor*xlast.user) 
  
  typ.data  <- as.data.frame(long.out) %>%
    filter(!(time==0 & CENT==0)) %>%      # No Dose rec
    mutate(x=time, y=RESP) %>%            # Standardize
    mutate(PRED=RESP) %>% 
    mutate(vachette.cov1 = WT) %>% 
    # Dose number
    mutate(dosenr = 1) %>% 
    dplyr::select(ID,x,PRED,vachette.cov1,dosenr)
  
  # --------- Observations ---------------
  
  # individual data
  idata <- expand.idata(ID = c(1:nsim.indiv),
                        WT = c(wt1,wt2))
  
  # simulate ipred
  nid <- length(unique(idata$ID))
  
  # Add variable obs.times
  obs.times.var <- NULL
  for(i in c(1:nid)) 
  {
    add.times.var <- t(data.frame(times = c(obs.times + runif(length(obs.times),-0.2,0.2))))
    obs.times.var <- rbind(obs.times.var,add.times.var)
  }

  imy.ipred <- NULL
  imy.pred  <- NULL
  for(iid in c(1:nid))
  {
    imy.ipred <- rbind(imy.ipred, as.data.frame(my.model %>%
                                                  omat(IIV = dmat(ECL, EV, EEC50, EKIN, EBSLN)) %>%    
                                                  smat(SGMA = dmat(PROPVAL, ADDVAL)) %>%    
                                                  ev(amt=100,addl=0,ii=24,cmt=1) %>%
                                                  idata_set(idata[iid,]) %>%
                                                  mrgsim(tgrid=obs.times.var[iid,])))
    imy.pred <- rbind(imy.pred, as.data.frame(my.model %>%
                                                omat(IIV = dmat(0, 0, 0, 0, 0)) %>% 
                                                smat(SGMA = dmat(0, 0)) %>% 
                                                ev(amt=100,addl=0,ii=24,cmt=1) %>%
                                                idata_set(idata[iid,]) %>%
                                                mrgsim(tgrid=obs.times.var[iid,])))
  }
  
  indivsam.all <- as.data.frame(imy.ipred) %>%
    # Add PRED, keep required variables only
    left_join(as.data.frame(imy.pred)[,c('ID','time','RESP','WT')] %>% rename(PRED=RESP),
              by=c('ID','time','WT')) %>% 
    filter(!(time==0)) %>%     # No Dose rec
    mutate(x=time, y=RESP) %>%             # Standardize
    mutate(vachette.cov1 = WT) %>% 
    mutate(IPRED = RESP) %>% 
    mutate(OBS = RESP_RUV) %>% 
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
      for(iid in c(1:nid))
      {
        imy.ipred <- rbind(imy.ipred, as.data.frame(my.model %>%
                                                      omat(IIV = dmat(ECL, EV, EEC50, EKIN, EBSLN)) %>%    
                                                      smat(SGMA = dmat(PROPVAL, ADDVAL)) %>%     # PROP, ADD 
                                                      ev(amt=100,addl=0,ii=24,cmt=1) %>%
                                                      idata_set(idata[iid,]) %>%
                                                      mrgsim(tgrid=obs.times.var[iid,])))
        
        imy.pred <- rbind(imy.pred, as.data.frame(my.model %>%
                                                    omat(IIV = dmat(0, 0, 0, 0, 0)) %>% 
                                                    smat(SGMA = dmat(0, 0)) %>% 
                                                    ev(amt=100,addl=0,ii=24,cmt=1) %>%
                                                    idata_set(idata[iid,]) %>%
                                                    mrgsim(tgrid=obs.times.var[iid,])))
      }
      
      
      indivsam.ivpc <- as.data.frame(imy.ipred) %>%
        # Add PRED, keep required variables only
        left_join(as.data.frame(imy.pred)[,c('ID','time','RESP','WT')] %>% rename(PRED=RESP),
                  by=c('ID','time','WT')) %>% 
        filter(!(time==0)) %>%     # No Dose rec
        mutate(x=time, y=RESP) %>%             # Standardize
        mutate(vachette.cov1 = WT) %>% 
        mutate(IPRED = RESP) %>% 
        mutate(OBS = RESP_RUV) %>% 
        # Dose number
        mutate(dosenr = 1) %>% 
        dplyr::select(ID,x,PRED,IPRED,OBS,vachette.cov1,dosenr) %>% 
        mutate(isim=ivpc)
      
      # Collect but keep required variables only
      sim.data <- rbindlist(list(sim.data, 
                                     indivsam.ivpc %>% dplyr::select(isim,ID,x,PRED,IPRED,OBS,vachette.cov1,dosenr)))
    }
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

