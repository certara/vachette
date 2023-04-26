# ------------------------------------------------------------------------- 
#  Sponsor           : Merck 
#  Project Number    : MRK-FTE-VRSV-640
#  Project Root Path : Local
# ------------------------------------------------------------------------- 
#  Program : vachette-example-oral-absorption-v38.R
#  Author  : Jos Lommerse - Certara 
#  Date    : 13 Feb 2023
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


sim.pembro <- function(nsim.indiv,iiv,ruv,nvpc,SAVE=F,PROP=T)
{
  my.model <- mcode("pembro3", pembro3)
  simtag   <- "pembro-v1" # script version number
  
  # Note: Time=0 and through already give max conc, so new dosenr start there
  
  # ----------- Model settings ------------------
  
  set.seed(20230126)
  
  obs.times.q3w     <- c(0.1,   0.2,  0.4,  1,  2,  4, 10, 21, 
                         21.1, 21.2, 21.4, 22, 23, 25, 31, 42) # Days, Q3W
  obs.times.q2w     <- c(0.1,   0.2 , 0.4,  1,  2,  5, 14,
                         14.1, 14.2, 14.4, 15, 16, 19, 28,
                         28.1, 28.2, 28.4, 29, 30, 33, 42) # Days, Q2W
  
  xstop         <- 42   # Actual observations end at xstop = 48 hr
  xlast.user    <- 400  # User provided x value to max. simulate out - replaces expand.factor
  
  step.x.factor <- 1.5
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
    mrgsim(delta=mrgdelta,end=step.x.factor*xlast.user) 
  long.out.q3w <- my.model %>%
    omat(IIV = dmat(0, 0)) %>% 
    # smat(SGMA = dmat(0, 0)) %>%
    ev(amt=100,addl=1,ii=21,cmt=1) %>% # Every 2 wks
    idata_set(data[data$SCHED=='Q3W',]) %>%
    mrgsim(delta=mrgdelta,end=step.x.factor*xlast.user) 
  
  typ.data  <- rbind(
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
    mutate(PRED=CONC) %>% 
    mutate(vachette.cov1 = SCHED) %>% 
    mutate(vachette.cov2 = ALB)
  
  typ.data %>% ggplot(aes(x=x,y=y,group=ID,col=factor(dosenr)))+
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
                                                  mrgsim(tgrid=obs.times.var.q2w[iid,])) %>% 
                             mutate(SCHED = 'Q2W') %>% 
                             # dosenr
                             mutate(dosenr = ifelse(time<42,floor(time/14)+1,3)))
    imy.pred.q2w <- rbind(imy.pred.q2w, as.data.frame(my.model %>%
                                                 omat(IIV = dmat(0,0)) %>% 
                                                 smat(SGMA = dmat(0)) %>%  # Just prop err
                                                 ev(amt=100,addl=2,ii=14,cmt=1) %>% # Every 2 wks
                                                 idata_set(idata.q2w[iid,]) %>%
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
                                                          mrgsim(tgrid=obs.times.var.q3w[iid,])) %>% 
                             mutate(SCHED = 'Q3W') %>% 
                             # dosenr
                             mutate(dosenr = ifelse(time<42,floor(time/21)+1,2)))
    
    imy.pred.q3w <- rbind(imy.pred.q3w, as.data.frame(my.model %>%
                                                omat(IIV = dmat(0,0)) %>% 
                                                smat(SGMA = dmat(0)) %>%  # Just prop err
                                                ev(amt=100,addl=1,ii=21,cmt=1) %>% # Every 3 wks
                                                idata_set(idata.q3w[iid,]) %>%
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
    mutate(OBS = CONC_RUV) 
  
  # Check
  indivsam.all %>% 
    ggplot(aes(x=x))+
    geom_point(aes(y=OBS))+
    geom_point(aes(y=IPRED),col='red')+
    geom_point(aes(y=PRED),col='cyan')+
    geom_line(aes(y=PRED))+
    facet_grid(SCHED~ALB)
  # Prima
  
  # VPC, simulate nvpc times and pick same ID/timepoint combinations
  if(nvpc>0)
  {
    sim.data <- NULL
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
                                                              mrgsim(tgrid=obs.times.var.q2w[iid,])) %>% 
                                 mutate(SCHED = 'Q2W') %>% 
                                 # dosenr
                                 mutate(dosenr = ifelse(time<42,floor(time/14)+1,3)))
        
        imy.pred.q2w <- rbind(imy.pred.q2w, as.data.frame(my.model %>%
                                                            omat(IIV = dmat(0,0)) %>% 
                                                            smat(SGMA = dmat(0)) %>%  # Just prop err
                                                            ev(amt=100,addl=2,ii=14,cmt=1) %>% # Every 2 wks
                                                            idata_set(idata.q2w[iid,]) %>%
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
                                                              mrgsim(tgrid=obs.times.var.q3w[iid,])) %>% 
                                 mutate(SCHED = 'Q3W') %>% 
                                 # dosenr
                                 mutate(dosenr = ifelse(time<42,floor(time/21)+1,2)))
        
        imy.pred.q3w <- rbind(imy.pred.q3w, as.data.frame(my.model %>%
                                                            omat(IIV = dmat(0,0)) %>% 
                                                            smat(SGMA = dmat(0)) %>%  # Just prop err
                                                            ev(amt=100,addl=1,ii=21,cmt=1) %>% # Every 3 wks
                                                            idata_set(idata.q3w[iid,]) %>%
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
        mutate(OBS = CONC_RUV) %>% 
        mutate(isim = ivpc)
      
      # Collect but keep required variables only
      sim.data <- rbindlist(list(sim.data, 
                            indivsam.ivpc %>% dplyr::select(isim,ID,x,PRED,IPRED,OBS,vachette.cov1,vachette.cov2,dosenr)))
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
