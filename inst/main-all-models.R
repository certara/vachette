# -------------------------------------------------------------------------
#  Sponsor           : Merck
#  Project Number    : MRK-FTE-VRSV-640
#  Project Root Path : Local
# -------------------------------------------------------------------------
#  Program : vachette-main-v29.R
#  Author  : Jos Lommerse - Certara
#  Date    : 31 January 2023
#  Purpose : Vachette method for Visualization of Analyses with Covariates
# -------------------------------------------------------------------------
#  Software : R version 4.1.2 (2021-11-01)
#  Platform : ThinkPad P15s Gen 2 Model 20W6008AUS
#  Environment : Windows-10, 11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz
# -------------------------------------------------------------------------

# v53-james Reduced version for James, i.v. and oral absorption only

# Clear memory
rm(list=ls())

# List of possible examples (will be ~8 in total)
models <- c("iv",
            "oral-absorption")

# Model example to run
run.model  <- "iv"

tag        <- "v53-james"
script     <- paste0("vachette-main-all-models-",tag,".R")

# ------------ Settings ----------

# ----- Simulations -----

# Simulate observations or take observations from flat file (located in "../flat-files/ folder")
SIM_OBS           <- F  # Simulate typical curves/observation (else flat files will be read)
SIM_SAVE          <- T  # Save simulated data, note that flat files will not be overwritten
# Carry out VPC simulations too?
VVPC              <- T
# Simulate model with proportional or additve error?
PROP_SIM          <- T  # Else simulation will be using additive error

# Number of indiviudals per covariate
nsim.indiv       <- 15          # nsim individuals (nsim.indiv indivs per covariate combo) with each obs.times observations
# Number of simulated replicates
nsim.vpc         <- 100         # nsim.vpc number of reps --> ntimepoints*nsim.indiv*nsim.vpc data points

# ----- Vachette transformation -----

# Before Vachette transformation, carry out IIV correction? (For VPC, IIV_CORR should be F)
IIV_CORR          <- F
if(VVPC) IIV_CORR <- F

# Proportional or additive error for Vachette transformation
PROP_TR           <- T
ADD_TR            <- ifelse(PROP_TR,F,T)

# Carry out landmark position refinement step?
LM_REFINE         <- F  # Landmark refinement

# ------ General defaults (can be adjusted by user if required) ----------
tolend        <- 0.001  # Max tolerance to determine last.x (allowed minimizing difference open end ref and query)
tolnoise      <- 1e-8   # Max tolerance allowed between gridpoints to identify landmarks for stepsize=1 (first and second derivative extremes)
step.x.factor <- 1.5    # Factor x-steps between which y-values are lower than the SCALED tolerance
ngrid.open.end<- 100    # Number of grid points to characterize open end curves

# =================================================================

# Location of this script
home   <- dirname(rstudioapi::getActiveDocumentContext()$path)
getwd()
setwd(home)

library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
# should not use entire tidyverse dependency, parse out using exact dependencies.
library(tidyverse)
library(dtw)
library(rootSolve)
library(pracma)
library(gtools)
library(Rfast)
library(prospectr)
library(data.table)
library(Hmisc)

# Vachette functions
source('vachette-functions-v39-james.R')

# Mrgsove models
#source('vachette-models-v33-james.R')    # Slightly reduced variability for oral abs

# PMx examples
source('vachette-example-iv-v35-james.R')
source('vachette-example-oral-absorption-v38-james.R') # with SAVE option

# Utility functions
source('vachette-utils-v30-james.R')

set.seed(121)

# -------------- Graph titles of models ---------------------------

models.txt = c("iv" = "Intravenous (i.v.)",
               "oral-absorption" = "Oral absorption")
data.frame(models.txt) %>% remove_rownames()

# ------------------------------------------------------

model  <- run.model

###########################################################################
#                                                                         #
#                                USER INPUT                               #
#                                                                         #
###########################################################################

# Tell Vachette which covariate to choose as reference

# ---------------- IV ----------------------------
if(model=='iv')
{
  ref.cov1   <- 70    # e.g. (kg)
  ref.dose   <- 1     # First (and only) dose
}

# ---------------- SINGLE DOSE ORAL ----------------------------
if(model=='oral-absorption')
{
  ref.cov1   <- 70
  ref.dose   <- 1     # First (and only) dose
}

###########################################################################
#                                                                         #
#   GENERATE SIMULATED TYPICAL CURVES, OBSERVATIONS AND VPC REPLICATES    #
#                                                                         #
###########################################################################

# Region is the Vachette terminology for the time between two dose administrations
ref.region <- ref.dose

# ------ Simulate ----

# @James: simulation of data can be removed altogether, not a part of the R-package
if(SIM_OBS)
{
  # iv
  if(model=='iv')        output.typ   <- sim.iv(nsim.indiv, iiv=0, ruv=0, nvpc=0, SAVE=SIM_SAVE,PROP=PROP_SIM)
  if(model=='iv')        indivsam.obs <- sim.iv(nsim.indiv, iiv=1, ruv=1, nvpc=0, SAVE=SIM_SAVE,PROP=PROP_SIM)
  if(model=='iv' & VVPC) indivsam.vpc <- sim.iv(nsim.indiv, iiv=1, ruv=1, nvpc=nsim.vpc, SAVE=SIM_SAVE,PROP=PROP_SIM)
  # oral
  if(model=='oral-absorption')        output.typ   <- sim.oral.absorption(nsim.indiv, iiv=0, ruv=0, nvpc=0, SAVE=SIM_SAVE,PROP=PROP_SIM)
  if(model=='oral-absorption')        indivsam.obs <- sim.oral.absorption(nsim.indiv, iiv=1, ruv=1, nvpc=0, SAVE=SIM_SAVE,PROP=PROP_SIM)
  if(model=='oral-absorption' & VVPC) indivsam.vpc <- sim.oral.absorption(nsim.indiv, iiv=1, ruv=1, nvpc=nsim.vpc, SAVE=SIM_SAVE,PROP=PROP_SIM)
}

# ------ Retrieve from flat files ----

# @ James: to keep. Files to be stored in some sort of a library of example data
if(!SIM_OBS) # Read flat file
{
  # iv
  if(model=='iv')        output.typ   <- read.csv("../flat-files-james/vachette-example-iv-v35-james-typ.csv",stringsAsFactors = F)
  if(model=='iv')        indivsam.obs <- read.csv("../flat-files-james/vachette-example-iv-v35-james-obs.csv",stringsAsFactors = F)
  if(model=='iv' & VVPC) indivsam.vpc <- read.csv("../flat-files-james/vachette-example-iv-v35-james-vpc.csv",stringsAsFactors = F)
  # oral
  if(model=='oral-absorption')        output.typ   <- read.csv("../flat-files-james/vachette-example-oral-absorption-v38-james-typ.csv",stringsAsFactors = F)
  if(model=='oral-absorption')        indivsam.obs <- read.csv("../flat-files-james/vachette-example-oral-absorption-v38-james-obs.csv",stringsAsFactors = F)
  if(model=='oral-absorption' & VVPC) indivsam.vpc <- read.csv("../flat-files-james/vachette-example-oral-absorption-v38-james-vpc.csv",stringsAsFactors = F)
}

# -------------------------------------------------

# ---- For the script flow ----

# "Dummy" observations replicate number
indivsam.obs$isim      <- 1
# "Dummy" vpc simulated observations dataset
# @James: also fine to add if-statements throughout the code up to "VACHETTE TRANSFORMATION STEPS"
if(!VVPC) indivsam.vpc <- indivsam.obs

# Keep required data fileds only:
indivsam.obs  <- indivsam.obs %>% dplyr::select(isim,ID,x,PRED,IPRED,OBS,vachette.cov1,dosenr)
indivsam.vpc  <- indivsam.vpc %>% dplyr::select(isim,ID,x,PRED,IPRED,OBS,vachette.cov1,dosenr)

###########################################################################
#                                                                         #
#                    VACHETTE PREPARATION                                 #
#                       - UNIQUE COVARIATE COMBINATIONS                   #
#                                                                         #
###########################################################################

# Extract covariates/regimens (cov1, cov2, ....)
vachette.covs  <- names(output.typ)[grepl("vachette.cov",names(output.typ))]
nvachette.covs <- length(vachette.covs)

# Extract last observed x (look for max x/time in obs and sim data)
xstop <- max(indivsam.obs$x,indivsam.vpc$x)

# Define unique covariate combination: To be improved
if(nvachette.covs==1)
{
  obs.orig   <- indivsam.obs %>%
    mutate(COV = paste(paste(vachette.cov1)))
  sim.orig   <- indivsam.vpc %>%
    mutate(COV = paste(paste(vachette.cov1)))
  output.typ   <- output.typ %>%
    mutate(COV = paste(paste(vachette.cov1)))
}
if(nvachette.covs==2)
{
  obs.orig   <- indivsam.obs %>%
    mutate(COV = paste(paste(vachette.cov1,vachette.cov2)))
  sim.orig   <- indivsam.vpc %>%
    mutate(COV = paste(paste(vachette.cov1,vachette.cov2)))
  output.typ <- output.typ  %>%
    mutate(COV = paste(paste(vachette.cov1,vachette.cov2)))
}
# Etc. for 3, 4 5 etc different covariates/covariate values
# @James: please improve, possibly via new function

# ----- Define regions and provide a number to all combined covariate effects ----

# Region number
output.typ$region  <- NA
obs.orig$region    <- NA
sim.orig$region    <- NA

# Region type (open/closed)
output.typ$region.type  <- NA
obs.orig$region.type    <- NA
sim.orig$region.type    <- NA

# Covariate combinations number
output.typ$ucov  <- NA
obs.orig$ucov    <- NA
sim.orig$ucov    <- NA

n.ucov    <- 0
tab.ucov  <- NULL # Table with ucov properties

# All covariate combinations:
if(nvachette.covs==1) comb.ucov <- output.typ %>% expand(vachette.cov1)
if(nvachette.covs==2) comb.ucov <- output.typ %>% expand(vachette.cov1,vachette.cov2)

# Collect unique covariate/dose combinations and define region.type
# Add to info to data frames
if(nvachette.covs==1)
{
  for(i in c(1:dim(comb.ucov)[1])) #nrow()
  {
    # Get number of doses for each covariate combination:
    z <- output.typ %>% filter(vachette.cov1 == comb.ucov$vachette.cov1[i])
    ndose <- length(unique(z$dosenr))

    for(idose in c(1:ndose))
    {
      if (idose == ndose) region.type <- 'open'   # Last dose
      if (idose <  ndose) region.type <- 'closed' # Any dose before last dose

      # New unique combination
      n.ucov <- n.ucov + 1
      add <- data.frame(ucov = n.ucov,
                        vachette.cov1 = comb.ucov$vachette.cov1[i],
                        region = idose,
                        region.type = region.type)
      tab.ucov <- rbind(tab.ucov, add)

      # Selection of new unique combination
      sel.typ  <- output.typ$vachette.cov1 == comb.ucov$vachette.cov1[i] &
        output.typ$dosenr == idose
      sel.obs  <- obs.orig$vachette.cov1 == comb.ucov$vachette.cov1[i] &
        obs.orig$dosenr == idose
      sel.sim  <- sim.orig$vachette.cov1 == comb.ucov$vachette.cov1[i] &
        sim.orig$dosenr == idose

      # Add new unique combination info
      output.typ$region[sel.typ]   <- idose
      obs.orig$region[sel.obs]     <- idose
      sim.orig$region[sel.sim]     <- idose

      output.typ$ucov[sel.typ]   <- n.ucov
      obs.orig$ucov[sel.obs]     <- n.ucov
      sim.orig$ucov[sel.sim]     <- n.ucov

      output.typ$region.type[sel.typ]   <- region.type
      obs.orig$region.type[sel.obs]     <- region.type
      sim.orig$region.type[sel.sim]     <- region.type
    }
  }
}

if(nvachette.covs==2)
{
  for(i in c(1:dim(comb.ucov)[1]))
  {
    # Get number of doses for covariate combination:
    z <- output.typ %>% filter(vachette.cov1 == comb.ucov$vachette.cov1[i] &
                                 vachette.cov2 == comb.ucov$vachette.cov2[i])
    ndose <- length(unique(z$dosenr))

    for(idose in c(1:ndose))
    {
      if (idose == ndose) region.type <- 'open'
      if (idose <  ndose) region.type <- 'closed'

      # New unique combination
      n.ucov <- n.ucov + 1
      add <- data.frame(ucov = n.ucov,
                        vachette.cov1 = comb.ucov$vachette.cov1[i],
                        vachette.cov2 = comb.ucov$vachette.cov2[i],
                        region = idose,
                        region.type = region.type)
      tab.ucov <- rbind(tab.ucov, add)

      # Selection of new unique combination
      sel.typ  <- output.typ$vachette.cov1 == comb.ucov$vachette.cov1[i] &
        output.typ$vachette.cov2 == comb.ucov$vachette.cov2[i] &
        output.typ$dosenr == idose
      sel.obs  <- obs.orig$vachette.cov1 == comb.ucov$vachette.cov1[i] &
        obs.orig$vachette.cov2 == comb.ucov$vachette.cov2[i] &
        obs.orig$dosenr == idose
      sel.sim  <- sim.orig$vachette.cov1 == comb.ucov$vachette.cov1[i] &
        sim.orig$vachette.cov2 == comb.ucov$vachette.cov2[i] &
        sim.orig$dosenr == idose

      # Add new unique combination info
      output.typ$region[sel.typ]   <- idose
      obs.orig$region[sel.obs]     <- idose
      sim.orig$region[sel.sim]     <- idose

      output.typ$ucov[sel.typ]   <- n.ucov
      obs.orig$ucov[sel.obs]     <- n.ucov
      sim.orig$ucov[sel.sim]     <- n.ucov

      output.typ$region.type[sel.typ]   <- region.type
      obs.orig$region.type[sel.obs]     <- region.type
      sim.orig$region.type[sel.sim]     <- region.type
    }
  }
}
# @James: to be improve for any number of covariates. Create function?

# Add flag for reference
tab.ucov$ref <- "No"
for(i.ucov in c(1:dim(tab.ucov)[1]))
{
  if(nvachette.covs==1)
    if(tab.ucov$vachette.cov1[i.ucov] == ref.cov1 & tab.ucov$region[i.ucov] == ref.region)
      tab.ucov$ref[i.ucov] <- "Yes"
  if(nvachette.covs==2)
    if(tab.ucov$vachette.cov1[i.ucov] == ref.cov1 & tab.ucov$vachette.cov2[i.ucov] == ref.cov2 & tab.ucov$region[i.ucov] == ref.region)
      tab.ucov$ref[i.ucov] <- "Yes"
}

# Add reference flag to output.typ
output.typ <- output.typ %>%
  left_join(tab.ucov[,c('ucov','ref')],by='ucov')

###########################################################################
#                                                                         #
#                    VACHETTE TRANSFORMATION STEPS                        #
#                                                                         #
###########################################################################


# ----- Step 1: Done (modeling+observations) -----

# ----- Step 2: Correct for IIV (or omit) -----

# Determine dependent variable value y's
output.typ   <- output.typ   %>% mutate(y=PRED)

# iiv correction for both for query and reference
if(IIV_CORR) obs.orig$y <- obs.orig$OBS  - obs.orig$IPRED + obs.orig$PRED
if(IIV_CORR) sim.orig$y <- sim.orig$OBS  - sim.orig$IPRED + sim.orig$PRED

# omit iiv correction:
if(!IIV_CORR) obs.orig$y <- obs.orig$OBS
if(!IIV_CORR) sim.orig$y <- sim.orig$OBS

# ---------- Preparation step 3+4 (Vachette x,y-transformations) ----------

# @James: possibly this (large) piece of code into a function?

# Collect all Vachette query curves and original/transformed observations (incl reference)
curves.all    <- NULL
obs.all       <- NULL

my.ref.lm.all   <- NULL
my.query.lm.all <- NULL
lm.all <- NULL

# --------- Repeat for each covariate combination --------------------------

# # Loop for all combinations of covariates
for(i.ucov in c(1:dim(tab.ucov)[1]))
{
  # i.ucov <- 1

  # A. ----- Define reference and query typical curves and observations to transform  ----------

  # Ref: may change (by extensions), so define (again) every new combination
  if(nvachette.covs==1) ref <- output.typ %>% filter(!is.na(region) & vachette.cov1 == ref.cov1 & region == ref.region)
  if(nvachette.covs==2) ref <- output.typ %>% filter(!is.na(region) & vachette.cov1 == ref.cov1 & vachette.cov2 == ref.cov2 & region == ref.region)
  ref.ucov <- unique(ref$ucov)
  if(length(ref.ucov) != 1) stop("Error: length ref.ucov != 1")     # A single ref only possible
  ref.region.type <- tab.ucov$region.type[tab.ucov$ucov==ref.ucov]
  # ref = Typical reference curve

  # Query: may be the same as reference
  query.cov1         <- tab.ucov$vachette.cov1[tab.ucov$ucov==i.ucov]
  if(nvachette.covs==2)
    query.cov2         <- tab.ucov$vachette.cov2[tab.ucov$ucov==i.ucov]
  query.region       <- tab.ucov$region[tab.ucov$ucov==i.ucov]
  query.region.type  <- tab.ucov$region.type[tab.ucov$ucov==i.ucov]

  # Typical query curve
  query <- output.typ %>%
    filter(ucov == i.ucov) %>%
    mutate(ref = tab.ucov$ref[i.ucov])    # Flag for reference

  # We have to run Vachette twice if observation AND simulated replicates have to be generated (for VPC)
  # @James: you may be able to do this more efficient - leverage tidyvpc here?

  # Observations
  if (!VVPC)
  {
    obs.query   <- obs.orig %>%
      filter(ucov == i.ucov) %>%
      mutate(ref = tab.ucov$ref[i.ucov])  # Flag for reference data
  }
  # Simulated observations for VPC
  if (VVPC)
  {
    obs.query   <- sim.orig %>%
      filter(ucov == i.ucov)  %>%
      mutate(ref = tab.ucov$ref[i.ucov])  # Flag for reference data
  }

  #  B. -------- Get landmarks approximate positions ----

  # Tolerance to apply for finding maximums, minimums and inflection points
  # Small stepsize -> small tolerance
  # Large stepsize -> large tol
  tolapply <- tolnoise*(output.typ$x[2]-output.typ$x[1])   # Grid stepsize

  my.ref.lm.init    <- get.x.multi.landmarks(ref$x,ref$y,w=17,tol=tolapply)
  my.ref.lm.init$y  <- approx(ref$x,ref$y, xout=my.ref.lm.init$x)$y

  my.query.lm.init    <- get.x.multi.landmarks(query$x,query$y,w=17,tol=tolapply)
  my.query.lm.init$y  <- approx(query$x,query$y, xout=my.query.lm.init$x)$y

  if(!LM_REFINE)
  {
    my.ref.lm.refined      <- my.ref.lm.init
    my.query.lm.refined    <- my.query.lm.init
  }
  if(LM_REFINE)
  {
    my.ref.lm.refined    <- refine.x.multi.landmarks(x=ref$x,y=ref$y,lm=my.ref.lm.init,tol=0.01*tolapply)
    my.ref.lm.refined$y  <- approx(ref$x, ref$y, xout=my.ref.lm.refined$x)$y

    my.query.lm.refined    <- refine.x.multi.landmarks(query$x,query$y,lm=my.query.lm.init,tol=0.01*tolapply)
    my.query.lm.refined$y  <- approx(query$x,query$y, xout=my.query.lm.refined$x)$y
  }

  # Check if query and ref have same landmarks in same order
  if(dim(my.ref.lm.refined)[1] != dim(my.query.lm.refined)[1] |
     !sum(my.ref.lm.refined$type == my.query.lm.refined$type)>0)
    stop("No matching landmarks between reference and query")

  # C. ---------- Get last.x for reference and query typical curves ------------------

  scaling <- 'linear'

  # ------ Landmark last.x determination depends on segment type ---------

  # 1. ref and query both open ends ("default")
  # 2. ref is open end, query is closed --> get.query.open.end() for last x **ref** only
  # 3. ref is closed, query is open end --> get.query.open.end() for last x query only
  # 4. ref and query both closed        --> get.query.open.end() for last x query only

  query.region.type <- unique(query$region.type)
  if(length(query.region.type) != 1) stop("Error: length query.region.type != 1")

  # ---- Reference and query last x -----

  # 1. Find ref last x, Fit query last x
  if(ref.region.type == 'open' & query.region.type == 'open')
  {
    my.ref.lm   <- get.ref.x.open.end(ref$x,ref$y,my.ref.lm.refined,step.x.factor=step.x.factor,tol=tolend)
    my.query.lm <- get.query.x.open.end(ref,query,my.ref.lm,my.query.lm.refined,ngrid=ngrid.open.end,scaling=scaling)
  }
  # 2. Fit ref last x, Fix query last x
  if(ref.region.type == 'open' & query.region.type == 'closed')
  {
    my.ref.lm   <- get.query.x.open.end(query,ref,my.query.lm.refined,my.ref.lm.refined,ngrid=ngrid.open.end,scaling=scaling)
    my.query.lm <- my.query.lm.refined
  }
  # 3. Fix ref last x, Fit query last x
  if(ref.region.type == 'closed' & query.region.type == 'open')
  {
    my.ref.lm   <- my.ref.lm.refined
    my.query.lm <- get.query.x.open.end(ref,query,my.ref.lm,my.query.lm.refined,ngrid=ngrid.open.end,scaling=scaling)
  }
  # 4. Fit ref last x, Fix query last x
  if(ref.region.type == 'closed' & query.region.type == 'closed')
  {
    my.ref.lm   <- get.query.x.open.end(query,ref,my.query.lm.refined,my.ref.lm.refined,ngrid=ngrid.open.end,scaling=scaling)
    my.query.lm <- my.query.lm.refined
  }

  # @James recent developments
  # Overrides above (above needs reprogramming)
  # Always FIX ref and FIT query to find best fitting x ("open"). Then extrapolate if required
  if(ref.region.type == 'closed')
  {
    my.ref.lm   <- my.ref.lm.refined
    my.query.lm <- get.query.x.open.end(ref,query,my.ref.lm,my.query.lm.refined,ngrid=ngrid.open.end,scaling=scaling)
  }

  # Recalc landmark y's
  my.query.lm$y  <- approx(query$x,query$y, xout=my.query.lm$x)$y
  my.ref.lm$y    <- approx(ref$x,ref$y, xout=my.ref.lm$x)$y

  # Collect all landmarks
  my.ref.lm.all   <- rbind(my.ref.lm.all,   my.ref.lm   %>% mutate(i.ucov = i.ucov))
  my.query.lm.all <- rbind(my.query.lm.all, my.query.lm %>% mutate(i.ucov = i.ucov))

  lm.all <- rbind(lm.all, my.query.lm %>% mutate(ucov = i.ucov))

  # D. ----- Define and map segments -----------

  nseg     <- dim(my.ref.lm)[1]-1

  # Initial assignment segments to typical curves, collect all query curves
  ref$seg   <- NA
  query$seg <- NA
  for(iseg in c(1:nseg))
  {
    if(iseg==1) # includes first datapoint
    {
      ref <- ref %>%
        mutate(seg = ifelse(x>=my.ref.lm$x[iseg] & x<=my.ref.lm$x[iseg+1],iseg,seg))
      query <- query %>%
        mutate(seg = ifelse(x>=my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1],iseg,seg))
    }
    if(iseg>1)
    {
      ref <- ref %>%
        mutate(seg = ifelse(x>my.ref.lm$x[iseg] & x<=my.ref.lm$x[iseg+1],iseg,seg))
      query <- query %>%
        mutate(seg = ifelse(x>my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1],iseg,seg))
    }
  }

  # E. ---- x-scaled query segments -----------

  # New query.scaled data frame
  query.scaled <- NULL
  for(iseg in c(1:nseg))
  {
    # Ref segment length
    ref.seg.length <- my.ref.lm$x[iseg+1] - my.ref.lm$x[iseg]
    # Query segment length
    query.seg.length <- my.query.lm$x[iseg+1] - my.query.lm$x[iseg]
    # Scaling factor
    my.query.add <- query %>%
      filter((nseg==1 & x>=my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1]) |
               (nseg>1 & x>my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1])) %>%
      mutate(seg=iseg) %>%
      mutate(y.scaling     = 1, y.scaled =1) %>%                # dummies, To remove
      mutate(x.shift.ref   = 0 - my.ref.lm$x[iseg]) %>%
      mutate(x.shift.query = 0 - my.query.lm$x[iseg]) %>%
      mutate(x.trans       = x + x.shift.query) %>%             # May be redefine to my.ref.lm$x[iseg]
      mutate(x.scaling     = ifelse(query.seg.length>0,
                                    ref.seg.length/query.seg.length,
                                    NA)) %>%                    # NA causes troubles?
      mutate(x.scaled      = x.scaling*x.trans - x.shift.ref)   # Back to corresponding ref position

    query.scaled <- rbind(query.scaled,my.query.add)
  }

  # F. ---- Check for need of extended reference curve by extrapolation -----------

  # Default no extension for extrapolation
  tab.extension <- NULL

  # Check if there are observation at all
  OBSERV <- ifelse(dim(obs.query)[1]>0,T,F)

  # If all observations are not before query last x, then extrapolate and add to reference region curve
  if (OBSERV)
  {
    if(max(my.query.lm$x) < max(obs.query$x))
    {

      # @James: I have to check this part of the code using a suitable example

      message('---------------------------------------------------')
      message('EXTENSION REF CURVE >> extrapolation needs checking')
      message('---------------------------------------------------')
      print(paste0("**** EXTEND REFERENCE CURVE FOR query i.ucov ",i.ucov," (Region: ",tab.ucov$region[tab.ucov$ucov==i.ucov],") *****"))

      max.obs.x <- max(obs.query$x)
      # On the typical query curve:
      max.obs.y <- approx(query$x,query$y,xout=max.obs.x)$y

      # Extrapolate (linear) max.obs.x in translated x of query curve
      xmax.scaling <- approxExtrap(x=query.scaled$x.trans, query.scaled$x.scaling, xout = (max.obs.x+query.scaled$x.shift.query[1]))$y
      ymax.scaling <- approxExtrap(x=query.scaled$x.trans, query.scaled$y.scaling, xout = (max.obs.x+query.scaled$x.shift.query[1]))$y
      xmax.scaled  <- xmax.scaling * (max.obs.x+query.scaled$x.shift.query[1])
      ymax.scaled  <- ymax.scaling * max.obs.y

      # Add to scaling data frame
      # Last known query segment "continues"
      cur.query.seg <- query$seg[dim(query[!is.na(query$seg),])[1]]
      query.scaled.add <- query.scaled %>%
        slice(1) %>%
        mutate(x=max.obs.x, # By default will be last observation x too, so OK with segment assignment
               y=max.obs.y,
               x.trans  = max.obs.x+query.scaled$x.shift.query[1],
               x.scaled = xmax.scaled,
               y.scaled = ymax.scaled,
               x.scaling = xmax.scaling,
               y.scaling = ymax.scaling) %>%
        mutate(seg=cur.query.seg)
      query.scaled <- rbind(query.scaled,query.scaled.add)

      tab.extension <- rbind(tab.extension,query.scaled.add)

      # The method ensures (?) that query last x superimposes on reference last x after scaling
      # So in principle we just have to add the additional query x.scaled/y.scaled to ref

      # Last ref point is extrapolated, so segment "continues"
      cur.ref.seg <- ref$seg[dim(ref[!is.na(ref$seg),])[1]]
      ref.add <- ref %>%
        slice(1) %>%
        mutate(x=xmax.scaled,
               y=ymax.scaled) %>%
        mutate(seg=cur.ref.seg)

      ref <- rbind(ref,ref.add)

      # Check:
      if(cur.ref.seg != cur.query.seg) stop("Error segment number assignment in reference extrapolation block")

      print("**** END EXTENSION *****")
    }
  }

  # Collect all typical curves with scaling factors and scaled x,y values
  query.scaled$y.scaled <- approx(ref$x,ref$y,xout=query.scaled$x.scaled)$y
  curves.all            <- rbind(curves.all,query.scaled)

  #  Carry out x,y, transformation of observation only if there are query observations for this covariate combination
  if(dim(obs.query)[1]>0)
  {
    # Steps
    # a. (Re-)assign segments
    # b. Determine obs new x-position
    # c. Determine obs y-scaling ADDITIVE OR PROPORTIONAL

    # a. (Re-)assign segments, x.shift.ref, x.shift.query and x.scaling to observations and query curve (ref already done)
    obs.query$seg <- NA
    for(iseg in c(1:nseg))
    {
      extended.query.x <- ifelse(length(tab.extension$x[tab.extension$seg==iseg])>0,tab.extension$x[tab.extension$seg==iseg],NA)
      first.query.x    <- my.query.lm$x[iseg]
      last.query.x     <- ifelse(!is.na(extended.query.x),extended.query.x,my.query.lm$x[iseg+1])

      obs.query <- obs.query %>% mutate(seg = ifelse((nseg==1 & x>=first.query.x & x<=last.query.x) |
                                                       (nseg>1 & x>first.query.x & x<=last.query.x),iseg, seg))
      query     <- query     %>% mutate(seg = ifelse((nseg==1 & x>=first.query.x & x<=last.query.x) |
                                                       (nseg>1 & x>first.query.x & x<=last.query.x),iseg, seg))
    }

    # ----- Steps 3: Vachette x transformation -----

    # b. Determine obs new x-position
    obs.query.orig <- obs.query
    obs.query      <- NULL
    for(iseg in c(1:nseg))
    {
      # Ref segment length
      ref.seg.length <- my.ref.lm$x[iseg+1] - my.ref.lm$x[iseg]
      # Query segment length
      query.seg.length <- my.query.lm$x[iseg+1] - my.query.lm$x[iseg]
      # Scaling factor
      obs.query.add <- obs.query.orig %>%
        filter((nseg==1 & x>=my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1]) |
                 (nseg>1 & x>my.query.lm$x[iseg] & x<=my.query.lm$x[iseg+1])) %>%
        mutate(seg=iseg) %>%
        mutate(y.scaling     = 1, y.scaled =1) %>%               # dummies, To remove
        mutate(x.shift.ref   = 0 - my.ref.lm$x[iseg]) %>%
        mutate(x.shift.query = 0 - my.query.lm$x[iseg]) %>%
        mutate(x.trans       = x + x.shift.query) %>%            # May be redefine to my.ref.lm$x[iseg]
        mutate(x.scaling     = ref.seg.length/query.seg.length) %>%
        mutate(x.scaled      = x.scaling*x.trans - x.shift.ref)  # Back to corresponding ref position

      obs.query <- rbind(obs.query,obs.query.add)
    }

    # ----- Steps 4: Vachette y transformation -----

    # c. Determine obs y-scaling ADDITIVE OR PROPORTIONAL

    # In case of proportional error model
    if(PROP_TR)
    {
      # PROPORTIONAL - log scale addition. Take full ref and curve curves ....
      ref            <- ref   %>% mutate(ylog = ifelse(y>0,log(y),NA))
      query          <- query %>% mutate(ylog = ifelse(y>0,log(y),NA))
      obs.query      <- obs.query %>% mutate(ylog = ifelse(y>0,log(y),NA))

      # Position on query curve
      obs.query$ylog.query       <- ifelse(!is.na(obs.query$ylog),
                                           approx(x=query$x, y=query$ylog, xout = obs.query$x)$y,
                                           NA)
      # log difference obs.query to query curve
      obs.query$ylog.diff        <- obs.query$ylog - obs.query$ylog.query

      # Position on ref curve (position at x.scaled!!)
      obs.query$ylog.ref         <- ifelse(!is.na(obs.query$ylog),
                                           approx(x=ref$x, y=ref$ylog, xout = obs.query$x.scaled)$y,
                                           NA)
      # New position
      obs.query$ylog.scaled      <- obs.query$ylog.diff + obs.query$ylog.ref
      obs.query$y.scaled         <- exp(obs.query$ylog.scaled)
    }

    # In case of additive error model
    if(ADD_TR)
    {
      # ADDITIVE - normal scale addition

      # Position on query curve
      obs.query$y.query       <- ifelse(!is.na(obs.query$y),
                                        approx(x=query$x, y=query$y, xout = obs.query$x)$y,
                                        NA)
      # difference obs.query to query curve
      obs.query$y.diff        <- obs.query$y - obs.query$y.query

      # Position on ref curve (position at x.scaled!!)
      obs.query$y.ref         <- ifelse(!is.na(obs.query$y),
                                        approx(x=ref$x, y=ref$y, xout = obs.query$x.scaled)$y,
                                        NA)
      # New position
      obs.query$y.scaled      <- obs.query$y.diff + obs.query$y.ref
    }

    obs.all <- rbind(obs.all,obs.query)

  } # If there are observations to transform

} # All ucov

###########################################################################
#                                                                         #
#                    VACHETTE OUTPUT                                      #
#                                                                         #
###########################################################################

# Save Vachette transformed data
write.table(output.typ,paste0("../tables-james/vachette-curves-",model,"-",tag,".csv"),
            row.names=F,col.names=T,sep=',',na='.',quote=F)
if (!VVPC) write.table(obs.all,paste0("../tables-james/vachette-obs-query-",model,"-",tag,".csv"),
                       row.names=F,col.names=T,sep=',',na='.',quote=F)
if (VVPC)  write.table(obs.all,paste0("../tables-james/vachette-sim-query-",model,"-",tag,".csv"),
                       row.names=F,col.names=T,sep=',',na='.',quote=F)

###########################################################################
#                                                                         #
#                    VACHETTE DIAGNOSTICS                                 #
#                                                                         #
###########################################################################

# --------- Distances to curves -----------

for(i.ucov in c(1:n.ucov))
{
  typ.curve <- output.typ %>% filter(ucov==i.ucov)
  ref.curve <- output.typ %>% filter(ucov==ref.ucov)
  obs       <- obs.all    %>% filter(ucov==i.ucov)

  # Additive distances to original curves
  obs.all$dist.add.orig[obs.all$ucov==i.ucov] <- approx(typ.curve$x,typ.curve$y,xout=obs$x)$y - obs$y
  # Vachette transformed - distances to ref curve
  obs.all$dist.add.transformed[obs.all$ucov==i.ucov] <- approx(ref.curve$x,ref.curve$y,xout=obs$x.scaled)$y - obs$y.scaled

  # Proportional distances to original curves (assume all y > 0)
  # obs.all$dist.prop.orig[obs.all$ucov==i.ucov]        <- log(approx(typ.curve$x,typ.curve$y,xout=obs$x)$y) - log(obs$y)
  obs.all$dist.prop.orig[obs.all$ucov==i.ucov]        <- approx(typ.curve$x,log(typ.curve$y),xout=obs$x)$y - log(obs$y)
  # Vachette transformed - distances to ref curve
  # obs.all$dist.prop.transformed[obs.all$ucov==i.ucov] <- log(approx(ref.curve$x,ref.curve$y,xout=obs$x.scaled)$y) - log(obs$y.scaled)
  obs.all$dist.prop.transformed[obs.all$ucov==i.ucov] <- approx(ref.curve$x,log(ref.curve$y),xout=obs$x.scaled)$y - log(obs$y.scaled)
}


###########################################################################
#                                                                         #
#                    VACHETTE PLOTTING                                    #
#                                                                         #
###########################################################################


# ---------- Diagnostic plots -----------------

# Dummy:
lm.all$seg <- NA

p.scaled.typical.curves.landmarks <- curves.all %>%
  ggplot(aes(x=x,y=y,group=ucov,
             col = factor(seg)))+
  geom_line(lwd=1)+
  geom_line(data=curves.all %>% filter(ref=="Yes"),col='black',lty=2,lwd=1,alpha=0.5)+
  # Add landmark positions
  geom_point(data=lm.all,pch=3,size=4,col='black',stroke = 2)+
  coord_cartesian(xlim=c(0,xstop)) +
  labs(title=paste0(models.txt[model]," - Typical curve segments"),
       subtitle = "Dashed line: Reference typical curve",
       caption=script,
       col="Segment")+
  render

# Full curves with landmarks
p.scaled.typical.full.curves.landmarks <- curves.all %>%
  ggplot(aes(x=x,y=y,group=ucov,
             col = factor(seg)))+
  geom_line(lwd=1)+
  geom_line(data=curves.all %>% filter(ref=="Yes"),col='black',lty=2,lwd=1,alpha=0.5)+
  # Add landmark positions
  geom_point(data=lm.all,pch=3,size=4,col='black',stroke = 2)+
  labs(title=paste0(models.txt[model]," - Typical curve segments"),
       subtitle = "Dashed line: Reference typical curve",
       caption=script,
       col="Segment")+
  render

# Scaling factors
p.scaling.factor <- curves.all %>%
  mutate(mycurve = ifelse(ref=='Yes','Reference','Query')) %>%
  ggplot(aes(x=x,y=x.scaling,col=factor(seg)))+
  geom_line(lwd=1)+
  facet_wrap(~paste(mycurve," covariate=",COV))+
  coord_cartesian(xlim=c(NA,max(obs.all$x,obs.all$x.scaled)),
                  ylim=c(0,max(curves.all$x.scaling[curves.all$x<=max(obs.all$x,obs.all$x.scaled)])))+
  labs(title=paste0(models.txt[model]," - x-scaling factors"),
       subtitle = paste0(if(ADD_TR) "Additive Error",if(PROP_TR) "Proportional Error"),
       caption=script,
       col="Segment")+
  render


p.scaled.typical.curves <- curves.all %>%
  ggplot(aes(x=x,y=y,group=ucov))+
  # geom_line(lwd=1.5,alpha=0.5)+
  geom_line(data=curves.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Yes'),lwd=1.5,alpha=0.60)+
  geom_line(data=curves.all %>% filter(ref=='No'),aes(x=x,y=y,col='No'),lwd=1.5)+
  geom_line(aes(x=x.scaled, y=y.scaled), col='black',lwd=0.75,lty=2)+
  # scale_color_manual(values=c('red','blue')) +
  # scale_linetype_manual(values = c(3, 1),
  #                       labels = c("Covariate", "Reference")) +
  coord_cartesian(xlim=c(0,xstop)) +
  scale_color_manual(name='Reference\ncurve',
                     breaks=c('No', 'Yes'),
                     values=c('No'='blue',
                              'Yes'='red'))+
  labs(title=paste0(models.txt[model]," - Covariate typical curves"),
       subtitle = "Dashed lines: Covariate typical curves after Vachette transformation",
       caption=script,
       col="Covariate value\n(Reference)")+
  render

# Copy ref curves to each ID
myids          <- unique(obs.all$ID)
curves.all.ids <- NULL
for(iid in c(1:length(myids))) {curves.all.ids <- rbind(curves.all.ids,curves.all %>% mutate(ID=myids[iid]))}

# Overlay observation curves per ID
p.scaled.observation.curves <- obs.all %>%
  ggplot(aes(x=x,y=y,group=paste(ID,COV)))+
  geom_line(data=curves.all.ids %>% filter(ref=="Yes"),aes(x=x,y=y,col='Typical reference'),lwd=1)+
  geom_line(data=obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference observations'))+
  geom_point(data=obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference observations'))+
  geom_line(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query observations'))+
  geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query observations'))+
  geom_line(data=obs.all %>% filter(ref=='No'),aes(x=x.scaled, y=y.scaled, col='Query transformed'),lty=2)+
  geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x.scaled, y=y.scaled, col='Query transformed'))+
  # scale_color_manual(values=c('red','blue')) +
  # scale_linetype_manual(values = c(3, 1),
  #                       labels = c("Covariate", "Reference")) +
  scale_color_manual(name='Data type',
                     breaks=c('Query transformed',
                              'Query observations',
                              'Reference observations',
                              'Typical reference'),
                     values=c('Query transformed'='purple',
                              'Query observations'='blue',
                              'Reference observations'='red',
                              'Typical reference'='black'))+
  coord_cartesian(xlim=c(0,xstop)) +
  labs(title=paste0(models.txt[model]," - Individual observation curves"),
       subtitle = "Dashed lines: Individual observation curves after Vachette transformation",
       caption=script,
       col="Covariate value\n(Reference)")+
  render

# Observation curve for each ID
for(iid in c(1:length(myids))) {curves.all.ids <- rbind(curves.all.ids,curves.all %>% mutate(ID=myids[iid]))}
p.scaled.observation.curves.by.id <- obs.all %>%
  ggplot(aes(x=x,y=y,group=paste(ID,COV)))+
  geom_line(data=curves.all.ids %>% filter(ref=="Yes"),aes(x=x,y=y,col='Typical reference'),lwd=1)+
  geom_line(data=obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference observations'))+
  geom_point(data=obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference observations'))+
  geom_line(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query observations'))+
  geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query observations'))+
  geom_line(data=obs.all %>% filter(ref=='No'),aes(x=x.scaled, y=y.scaled, col='Query transformed'),lty=2)+
  geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x.scaled, y=y.scaled, col='Query transformed'))+
  # scale_color_manual(values=c('red','blue')) +
  # scale_linetype_manual(values = c(3, 1),
  #                       labels = c("Covariate", "Reference")) +
  scale_color_manual(name='Data type',
                     breaks=c('Query transformed',
                              'Query observations',
                              'Reference observations',
                              'Typical reference'),
                     values=c('Query transformed'='purple',
                              'Query observations'='blue',
                              'Reference observations'='red',
                              'Typical reference'='black'))+
  facet_wrap(~ID)+
  coord_cartesian(xlim=c(0,xstop)) +
  labs(title=paste0(models.txt[model]," - Individual observation curves"),
       subtitle = "Dashed lines: Individual observation curves after Vachette transformation",
       caption=script,
       col="Covariate value\n(Reference)")+
  render

# Transformation steps
# Select first or second segment as example and first region
p.seg    <- ifelse(nseg>1,2,1)
p.region <- 1
p.transformation <- obs.all %>%
  ggplot(aes(x=x,y=y,col=factor(seg),pch=factor(paste(vachette.cov1)))) +
  # geom_point(data = obs.all %>% filter(seg==2),size=2,alpha=0.25) +
  # geom_point(data = obs.all %>% filter(ref=="No",seg==2),aes(x=x.scaled,y=y.scaled),col='purple',size=2) +
  # Transformation arrows
  # geom_segment(aes(x=x,y=y,xend=x.scaled,yend=y.scaled),
  #              arrow = arrow(length = unit(0.2, "cm")),col='grey') +
  geom_line(data=curves.all %>% filter(ref!="No"),lwd=1,col='grey') +
  geom_line(data=curves.all %>% filter(ref=="No"),lwd=1,col='grey') +
  geom_line(data=curves.all %>% filter(ref!="No",seg==p.seg),lwd=1,col='red') +
  geom_line(data=curves.all %>% filter(ref=="No",seg==p.seg,region==p.region),lwd=1,col='blue') +
  geom_line(data=curves.all %>% filter(ref=="No",seg==p.seg,region==p.region),aes(x=x.scaled),lwd=1,col='cyan') +
  geom_line(data=curves.all %>% filter(ref=="No",seg==p.seg,region==p.region),aes(x=x.scaled,y=y.scaled),lwd=1,col='purple',lty=2) +
  coord_cartesian(xlim=c(0,max(obs.all$x,obs.all$x.scaled)))+
  labs(title=paste0(models.txt[model]," - Example Vachette segment-",p.seg," transformation"),
       subtitle = paste0(if(ADD_TR) "Additive Error",if(PROP_TR) "Proportional Error"),
       caption=script,
       col="Segm",
       shape="Covs")+
  render

# Additive distances before and after transformation
p.add.distances <- obs.all %>%
  ggplot(aes(x=dist.add.orig,y=dist.add.transformed,col=factor(seg)))+
  geom_abline(slope=1)+
  geom_point()+
  labs(title=paste0(models.txt[model]," - Normal distances: original and after transformation"),
       subtitle = paste0(if(ADD_TR) "Additive Error",if(PROP_TR) "Proportional Error"),
       caption=script,
       x = 'Original distance',
       x = 'Distance after transformation',
       col="Segm.")+
  render

# Proportional distances before and after transformation
p.prop.distances <- obs.all %>%
  ggplot(aes(x=dist.prop.orig,y=dist.prop.transformed,col=factor(seg)))+
  geom_abline(slope=1)+
  geom_point()+
  labs(title=paste0(models.txt[model]," - Proportional distances: original and after transformation"),
       subtitle = paste0(if(ADD_TR) "Additive Error",if(PROP_TR) "Proportional Error"),
       caption=script,
       x = 'Original distance on log scale',
       x = 'Distance on log scale after transformation',
       col="Segm.")+
  render

# -------------- Transformation etc. ---------------------------

if(nvachette.covs==1)
{
  p.obs.ref.query <- obs.all %>%
    ggplot(aes(x=x,y=y)) +
    geom_line(data=curves.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query'),lwd=1) +
    geom_line(data=curves.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference'),lwd=1) +
    geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query'),pch=19) +
    geom_point(data=obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference'),pch=19) +

    scale_color_manual(name='Data type',
                       breaks=c('Query',
                                'Reference'),
                       values=c('Query'='blue',
                                'Reference'='red'))+
    coord_cartesian(xlim=c(0,max(obs.all$x))) +
    labs(title=paste0(models.txt[model]," - Observation + typical curves"),
         subtitle = paste0(if(ADD_TR) "Additive Error",if(PROP_TR) "Proportional Error"),
         caption=script)+
    render

  p.obs.cov <- obs.all %>%
    ggplot(aes(x=x,y=y)) +
    geom_line(data=curves.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference'),lwd=1) +
    geom_line(data=curves.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query'),lwd=1) +
    geom_point(data=obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference'),pch=19) +
    geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query'),pch=19) +

    scale_color_manual(name='Data type',
                       breaks=c('Query',
                                'Reference'),
                       values=c('Query'='blue',
                                'Reference'='red'))+
    facet_wrap(~paste("ucov",ucov)) +
    coord_cartesian(xlim=c(0,max(obs.all$x))) +
    labs(title=paste0(models.txt[model]," - Observation + typical curves"),
         subtitle = paste0(if(ADD_TR) "Additive Error",if(PROP_TR) "Proportional Error"),
         caption=script)+
    render

  p.vachette.arrow <- obs.all %>%
    ggplot(aes(x=x,y=y)) +
    # Transformation arrows
    geom_segment(aes(x=x,y=y,xend=x.scaled,yend=y.scaled),
                 arrow = arrow(length = unit(0.2, "cm")),col='grey') +

    geom_line(data=curves.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference'),lwd=1) +
    geom_line(data=curves.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query'),lwd=1) +
    geom_point(data=obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference'),pch=19,alpha=0.25) +
    geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x,y=y,col='Query'),pch=19,alpha=0.25) +
    geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x.scaled,y=y.scaled,col='Transformed'),pch=19) +

    scale_color_manual(name='Data type',
                       breaks=c('Query',
                                'Reference',
                                'Transformed'),
                       values=c('Query'='blue',
                                'Reference'='red',
                                'Transformed' = 'purple'))+

    coord_cartesian(xlim=c(0,max(obs.all$x,obs.all$x.scaled)))+
    labs(title=paste0(models.txt[model]," - Observations + transformations"),
         subtitle = paste0(if(ADD_TR) "Additive Error",if(PROP_TR) "Proportional Error"),
         caption=script)+
    render

  # Vachette final
  p.vachette <- obs.all %>%
    ggplot(aes(x=x,y=y)) +

    geom_line(data=curves.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference'),lwd=1) +
    geom_point(data=obs.all %>% filter(ref=='Yes'),aes(x=x,y=y,col='Reference'),pch=19) +
    geom_point(data=obs.all %>% filter(ref=='No'),aes(x=x.scaled,y=y.scaled,col='Query\nTransformed'),pch=19) +

    scale_color_manual(name='Data type',
                       breaks=c('Query',
                                'Reference',
                                'Query\nTransformed'),
                       values=c('Query'='blue',
                                'Reference'='red',
                                'Query\nTransformed' = 'purple'))+

    coord_cartesian(xlim=c(0,max(obs.all$x,obs.all$x.scaled)))+
    labs(title=paste0(models.txt[model]," - Observations + transformations"),
         subtitle = paste0(if(ADD_TR) "Additive Error",if(PROP_TR) "Proportional Error"),
         caption=script)+
    render
}

###########################################################################
#                                                                         #
#                             SAVE PLOTS                                  #
#                                                                         #
###########################################################################

# Save plots in pdf-file:
if (!VVPC) pdf(paste0("../plots-james/vachette-obs-",model,"-",tag,".pdf"))
if (VVPC)  pdf(paste0("../plots-james/vachette-sim-",model,"-",tag,".pdf"))
# scope VVPC argument, in the context of VPC options
if(ADD_TR)  plot(p.add.distances)
if(PROP_TR) plot(p.prop.distances)
plot(p.scaled.typical.curves.landmarks)
plot(p.scaled.typical.full.curves.landmarks)
plot(p.scaling.factor)
plot(p.scaled.typical.curves)
plot(p.scaled.observation.curves)
plot(p.scaled.observation.curves.by.id)

plot(p.obs.ref.query)
plot(p.obs.cov)
plot(p.vachette.arrow)
plot(p.vachette)

dev.off()

plot(p.vachette)

###########################################################################
#                                                                         #
#                            SUMMARY SETTINGS                             #
#                                                                         #
###########################################################################

print(' ')
print("----- Summary ------")
print(' ')

print(paste0("MODEL"))
print(paste0("  ",run.model))
print(' ')

print("SIMULATIONS")
if(SIM_OBS)  print(paste0("  Observations of ",nsim.indiv," individuals per covariate have been simulated"))
if(SIM_OBS & PROP_SIM)  print("  Simulation of proportional residual error only")
if(SIM_OBS & !PROP_SIM)  print("  Simulation of additive residual error only")
if(SIM_OBS & VVPC) print(paste0("  Simulation of ",nsim.vpc," replicates for Vachette VPC done"))
if(!SIM_OBS) print("  Observations have been read from '../flat-files-james/' folder")
print(' ')

print("UNIQUE COVARIATES")
print(tab.ucov)
print(' ')


print("VACHETTE TRANSFORMATION")
if(IIV_CORR)   print(paste0("  Vachette IIV correction applied"))
if(!IIV_CORR)  print(paste0("  Vachette IIV correction not applied"))
if(PROP_TR)    print("  Vachette transformations based on preserving log scale distances (proportional)")
if(ADD_TR)     print("  Vachette transformations based on preserving normal scale distances (additive)")
if(LM_REFINE)  print("  Landmark position refinement carried out")
print(paste0("  Tolerance determination last.x: ",tolend))
print(paste0("  Tolerance identification maximums, minimums and inflection points: ",tolnoise))
print(paste0("  Factor to extend curve to determine last.x: ",step.x.factor))
print(paste0("  Number of gridpoints used for minimizing shape-difference open segments: ",ngrid.open.end))

print(' ')
print("----------------------")

# --- end ---


