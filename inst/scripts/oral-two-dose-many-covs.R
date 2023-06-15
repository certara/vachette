##########################################################
#                                                        #
#               ORAL ABSORPTION TWO DOSE WITH            #
#               MANY DIFFERENT COVARIATE VALUES          #
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
set.seed(sim.seed)

#define number of subjects to be simulated per replicate
nsim.indiv <- 100

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

# sim.oral.absorption <- function(nsim.indiv,iiv,ruv,nvpc,SAVE=F,PROP=T)
# {
  my.model <- mcode("oral1cmt", oral1cmt)
  simtag   <- "oral-two-dose-many-covs" # script version number
  
  # ----------- Model settings ------------------
  
  obs.times     <- c(0.5,1,2,4,6,8,12,18,24,36,48,
                     48.5,49,50,52,54,56,60,66,72,84,96)
  
  xend.typ      <- 600  # User provided x value to max. simulate out 
  
  mrgdelta      <- 0.1
  
  wt <- exp(rnorm(nsim.indiv, mean = log(73), sd = 0.2)) 
  
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
  inidata <- expand.idata(ID=1:nsim.indiv) %>%
    mutate(WT = wt) %>%
    mutate(amt = 0, addl=0, ii=0, cmt=2, evid=0)
  
  idata <- reframe((inidata %>% group_by(ID) %>% slice(1))[,1], time= obs.times) %>%
    left_join(inidata)
  
  doses <- expand.ev(amt=100,addl=1,ii=48,cmt=1) %>%
    dplyr::select(-ID)
  
  doses <- merge(idata[,c(1,3)] %>% group_by(ID) %>% slice(1), doses) 
  
  data <- rbind(doses, idata) %>%
    arrange(ID, time)
  
  # simulate pred
  long.out <- my.model %>%
    omat(IIV = dmat(0, 0, 0)) %>% 
    smat(SGMA = dmat(0, 0)) %>%
    idata_set(data %>% group_by(ID) %>% slice(1)) %>%
    ev(amt=100,addl=1,ii=48,cmt=1) %>%
    carry_out(EVID, ADDL, II) %>%
    mrgsim_df(start=0.1,delta=mrgdelta,end=xend.typ) 
  
  long.out %>%
    ggplot() +
    geom_line(aes(x = time, y = CP, group = ID)) +
    scale_x_continuous(limits = c(0,100))
  
  output.typ  <- long.out %>%
    #filter(!(time==0 & GUT==0)) %>%     # Only Dose rec
    mutate(x=time, y=CP) %>%       # Standardize
    rename(PRED=CP) %>% 
    mutate(vachette.cov1 = WT) %>% 
    # Dose number
    mutate(dosenr = ifelse(time <= 48,1,2),
           ID = ID + 1000000) %>%
    dplyr::select(ID,time,EVID,ADDL,II,PRED,WT)
  
  # --------- Observations ---------------
  
    # simulate ipred
  
  imy.ipred <- as.data.frame(my.model %>%
                               omat(IIV = dmat(ECL, EV, EKA)) %>%        # ECL, EV, EKA
                               smat(SGMA = dmat(PROPVAL, ADDVAL)) %>%          # PROP, ADD
                               #ev(amt=100,addl=1,ii=48,cmt=1) %>%
                               data_set(data) %>%
                               carry_out(EVID, ADDL, II) %>%
                               mrgsim())
  imy.pred <- as.data.frame(my.model %>%
                              omat(IIV = dmat(0, 0, 0)) %>%        # ECL, EV, EKA
                              smat(SGMA = dmat(0, 0)) %>%          # PROP, ADD 
                              #ev(amt=100,addl=1,ii=48,cmt=1) %>%
                              data_set(data) %>%
                              carry_out(EVID, ADDL, II) %>%
                              mrgsim())

  indivsam.all <- imy.ipred %>%
    dplyr::select(-GUT) %>%
    # Add PRED, keep required variables only
    left_join(imy.pred[,c('ID','time','CP','WT')] %>% rename(PRED=CP),
              by=c('ID','time','WT')) %>% 
    #filter(!(time==0 & GUT==0)) %>%     # No Dose rec
    mutate(x=time, y=CP)%>%       # Standardize
    mutate(vachette.cov1 = WT) %>% 
    rename(IPRED = CP) %>% 
    rename(DV = CP_RUV) %>% 
    # Dose number
    mutate(dosenr = ifelse(time <= 48,1,2)) %>%
    dplyr::select(ID,time,EVID,ADDL,II,DV,PRED,WT)

  
  # VPC, simulate nvpc times and pick same ID/timepoint combinations
  if(nvpc>0)
  {
    indivsam.vpc <- NULL
    for(ivpc in c(1:nvpc))
    {
      imy.ipred <- as.data.frame(my.model %>%
                                   omat(IIV = dmat(ECL, EV, EKA)) %>%        # ECL, EV, EKA
                                   smat(SGMA = dmat(PROPVAL, ADDVAL)) %>%          # PROP, ADD
                                   #ev(amt=100,addl=1,ii=48,cmt=1) %>%
                                   data_set(data) %>%
                                   carry_out(EVID, ADDL, II) %>%
                                   mrgsim())
      imy.pred <- as.data.frame(my.model %>%
                                  omat(IIV = dmat(0, 0, 0)) %>%        # ECL, EV, EKA
                                  smat(SGMA = dmat(0, 0)) %>%          # PROP, ADD 
                                  #ev(amt=100,addl=1,ii=48,cmt=1) %>%
                                  data_set(data) %>%
                                  carry_out(EVID, ADDL, II) %>%
                                  mrgsim())
      
      
      indivsam.ivpc <- as.data.frame(imy.ipred) %>%
        dplyr::select(-GUT) %>%
        # Add PRED, keep required variables only
        left_join(as.data.frame(imy.pred)[,c('ID','time','CP','WT')] %>% rename(PRED=CP)) %>% 
        #filter(!(time==0 & GUT==0)) %>%     # No Dose rec
        mutate(x=time, y=CP)%>%       # Standardize
        mutate(vachette.cov1 = WT) %>% 
        rename(IPRED = CP) %>% 
        rename(DV = CP_RUV) %>% 
        # Dose number
        mutate(dosenr = ifelse(time <= 48,1,2)) %>% 
        mutate(REP = ivpc)
      
      # Collect but keep required variables only
      indivsam.vpc <- rbindlist(list(indivsam.vpc, 
                            indivsam.ivpc %>% dplyr::select(REP,ID,time,EVID,ADDL,II,DV,PRED,WT)))
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
  
# }
