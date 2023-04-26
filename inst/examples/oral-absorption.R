# ------------------------------------------------------------------------- 
#  Sponsor           : Merck 
#  Project Number    : MRK-FTE-VRSV-640
#  Project Root Path : Local
# ------------------------------------------------------------------------- 
#  Program : vachette-example-oral-absorption-v38-james.R
#  Author  : Jos Lommerse - Certara 
#  Date    : 31 January 2023
#  Purpose : Model simulations using mrgsolve
# ------------------------------------------------------------------------- 
#  Software : R version 4.2.2 (2022-10-31 ucrt)
#  Platform : ThinkPad P15s Gen 2 Model 20W6008AUS
#  Environment : Windows-10, 11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz
# -------------------------------------------------------------------------


##########################################################
#                                                        #
#                    ORAL ABSORPTION                     #
#              DEFINE REFERENCE AND QUERY                #
#                                                        #
##########################################################


sim.oral.absorption <- function(nsim.indiv,iiv,ruv,nvpc,SAVE=F,PROP=T)
{
  my.model <- mcode("oral1cmt", oral1cmt)
  simtag   <- "oral-absorption-v38-james" # script version number
  
  # ----------- Model settings ------------------
  
  obs.times     <- c(0.5,1,2,4,6,8,12,18,24,36,48)
  
  xstop         <- 48   # Actual observations end at xstop = 48 hr
  xlast.user    <- 500  # User provided x value to max. simulate out - replaces expand.factor
  
  step.x.factor <- 1.5
  mrgdelta      <- 0.1
  
  wt1 <- 30
  wt2 <- 70
  
  # Variability
  # IIV
  ECL <- 0.1
  EV  <- 0.1
  EKA <- 0.1
  # Residual
  PROPVAL<- 0
  ADDVAL <- 0
  if(PROP)  PROPVAL<- 0.05
  if(!PROP) ADDVAL <- 0.0001
  
  # ---------------- SIMULATED TYPICAL/POPULATION CURVES -------------
  
  # population input data
  data <- expand.idata(ID=1,WT=c(wt1,wt2))
  
  # simulate pred
  long.out <- my.model %>%
    omat(IIV = dmat(0, 0, 0)) %>% 
    smat(SGMA = dmat(0, 0)) %>%
    ev(amt=100,addl=0,ii=24,cmt=1) %>%
    idata_set(data) %>%
    mrgsim(delta=mrgdelta,end=step.x.factor*xlast.user) 
  
  typ.data  <- as.data.frame(long.out) %>%
    filter(!(time==0 & GUT==0)) %>%     # No Dose rec
    mutate(x=time, y=CP) %>%       # Standardize
    mutate(PRED=CP) %>% 
    mutate(vachette.cov1 = WT) %>% 
    # Dose number
    mutate(dosenr = 1)
  
  # --------- Observations ---------------
  
  # individual data
  idata <- expand.idata(ID = c(1:nsim.indiv),
                        WT = c(wt1,wt2))
  nid <- length(unique(idata$ID))
  
  # Add variable obs.times
  obs.times.var <- NULL
  for(i in c(1:nid)) 
  {
    add.times.var <- t(data.frame(times = c(obs.times + runif(length(obs.times),-0.1,0.1))))
    obs.times.var <- rbind(obs.times.var,add.times.var)
  }
  
  # simulate ipred
  set.seed(20230126)
  
  imy.ipred <- NULL
  imy.pred  <- NULL
  for(iid in c(1:nid))
  {
    imy.ipred <- rbind(imy.ipred, as.data.frame(my.model %>%
                         omat(IIV = dmat(ECL, EV, EKA)) %>%        # ECL, EV, EKA
                         smat(SGMA = dmat(PROPVAL, ADDVAL)) %>%          # PROP, ADD
                         ev(amt=100,addl=0,ii=24,cmt=1) %>%
                         idata_set(idata[iid,]) %>%
                         mrgsim(tgrid=obs.times.var[iid,])))
    imy.pred <- rbind(imy.pred, as.data.frame(my.model %>%
                         omat(IIV = dmat(0, 0, 0)) %>%        # ECL, EV, EKA
                         smat(SGMA = dmat(0, 0)) %>%          # PROP, ADD 
                         ev(amt=100,addl=0,ii=24,cmt=1) %>%
                         idata_set(idata[iid,]) %>%
                         mrgsim(tgrid=obs.times.var[iid,])))
  }
  
  indivsam.all <- imy.ipred %>%
    # Add PRED, keep required variables only
    left_join(imy.pred[,c('ID','time','CP','WT')] %>% rename(PRED=CP),
              by=c('ID','time','WT')) %>% 
    filter(!(time==0 & GUT==0)) %>%     # No Dose rec
    mutate(x=time, y=CP)%>%       # Standardize
    mutate(vachette.cov1 = WT) %>% 
    mutate(IPRED = CP) %>% 
    mutate(OBS = CP_RUV) %>% 
    # Dose number
    mutate(dosenr = 1)

  
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
                                                      omat(IIV = dmat(ECL, EV, EKA)) %>%        # ECL, EV, EKA
                                                      smat(SGMA = dmat(PROPVAL, ADDVAL)) %>%          # PROP, ADD
                                                      ev(amt=100,addl=0,ii=24,cmt=1) %>%
                                                      idata_set(idata[iid,]) %>%
                                                      mrgsim(tgrid=obs.times.var[iid,])))
        imy.pred <- rbind(imy.pred, as.data.frame(my.model %>%
                                                    omat(IIV = dmat(0, 0, 0)) %>%        # ECL, EV, EKA
                                                    smat(SGMA = dmat(0, 0)) %>%          # PROP, ADD 
                                                    ev(amt=100,addl=0,ii=24,cmt=1) %>%
                                                    idata_set(idata[iid,]) %>%
                                                    mrgsim(tgrid=obs.times.var[iid,])))
      }
      
      indivsam.ivpc <- as.data.frame(imy.ipred) %>%
        # Add PRED, keep required variables only
        left_join(as.data.frame(imy.pred)[,c('ID','time','CP','WT')] %>% rename(PRED=CP),
                  by=c('ID','time','WT')) %>% 
        filter(!(time==0 & GUT==0)) %>%     # No Dose rec
        mutate(x=time, y=CP)%>%       # Standardize
        mutate(vachette.cov1 = WT) %>% 
        mutate(IPRED = CP) %>% 
        mutate(OBS = CP_RUV) %>% 
        # Dose number
        mutate(dosenr = 1) %>% 
        mutate(isim = ivpc)
      
      # Collect but keep required variables only
      sim.data <- rbindlist(list(sim.data, 
                            indivsam.ivpc %>% dplyr::select(isim,ID,x,PRED,IPRED,OBS,vachette.cov1,dosenr)))
    }
  }
  
  if(SAVE)
  {
    filetyp <- paste0("../flat-files-james/vachette-example-", simtag, "-typ.csv")
    fileobs <- paste0("../flat-files-james/vachette-example-", simtag, "-obs.csv")
    filevpc <- paste0("../flat-files-james/vachette-example-", simtag, "-vpc.csv")
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
