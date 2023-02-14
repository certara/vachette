# ------------------------------------------------------------------------- 
#  Sponsor           : Merck 
#  Project Number    : MRK-FTE-VRSV-640
#  Project Root Path : Local
# ------------------------------------------------------------------------- 
#  Program : vachette-vpc-v11-james.R
#  Author  : Jos Lommerse - Certara 
#  Date    : 31 January 2023
#  Purpose : VPC, pcVPC and V2PC (Vachette-VPC)
# ------------------------------------------------------------------------- 
#  Software : R version 4.2.2 (2022-10-31 ucrt)
#  Platform : ThinkPad P15s Gen 2 Model 20W6008AUS
#  Environment : Windows-10, 11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz
# -------------------------------------------------------------------------

# v51-james Check-VPC-script for James

# @James: perhaps somehow you see ways to more efficiently make use of tidyvpc?
#         the tidyvpc part should be part of the Vachette package if people want to run VPC
#         So generating some default plot like below: plot(p.vpc.vachette + render)

# Clear memory
rm(list=ls())

tag <- 'v11-james' 

# Required packages
library(dplyr)
library(tidyr)
library(tidyverse)
library(tidyvpc)
library(ggplot2)
library(gginnards)

source("vachette-utils-v30-james.R")

# "Observed" and simulated dataset tag
tagfile <- "v53-james"
model   <- 'iv'

# Select prediction interval (from "80perc", "90perc")
pred.int <- "90perc"

# Reference definition:
# ref.alb      <- "16"
# ref.sched    <- "Q3W"
# ref.region   <- 4  # First dose Q3W

# VPC settings:
if(model=='iv')              all.breaks   <- c(0.25,0.5,1,2,3,4,8,12)
if(model=='iv')              ref.breaks   <- c(0.25,0.5,1,2,3,4,8,12)
if(model=='oral-absorption') all.breaks   <- c(0.5,1,2,4,6,8,12,18,24,36,48) 
if(model=='oral-absorption') ref.breaks   <- c(0.5,1,2,4,6,8,12,18,24,36,48)

# ----- Locations (using R-Studio interface) ----------------

# Location of this script
home   <- dirname(rstudioapi::getActiveDocumentContext()$path)
getwd()
# Home
setwd(home)

################################################################################
#                                                                              #
#                     READ SIMULATION OUTPUT                                   #
#                                                                              #
################################################################################

# Read observations from file 
typ.all <- read.csv(paste0("../tables-james/vachette-curves-",model,"-",tagfile,".csv"), head=T)

# Read observations from file 
obs.all <- read.csv(paste0("../tables-james/vachette-obs-query-",model,"-",tagfile,".csv"), head=T) 

# Read simulated data from file 
sim.all <- read.csv(paste0("../tables-james/vachette-sim-query-",model,"-",tagfile,".csv"), head=T) %>% 
  arrange(isim)

# Characterize data
obs.all %>% 
  group_by(ucov) %>% 
  dplyr::summarise(nindiv = length(unique(ID)),
                   nobs=n())
sim.all %>% 
  group_by(ucov) %>% 
  dplyr::summarise(nindiv = length(unique(ID)),
                   nobs=n())


# Quick plot - correct data?
obs.all %>% ggplot(aes(x=x,y=y,col=factor(region)))+geom_point()
sim.all %>% ggplot(aes(x=x,y=y,col=factor(region)))+geom_point()

head(obs.all)
head(sim.all)
# OK y=OBS (no IIV correction)

# Make numeric and remove NA's
for(j in c(1:dim(obs.all)[2])) obs.all[,j] <- as.numeric(as.character(obs.all[,j]))
for(j in c(1:dim(sim.all)[2])) sim.all[,j] <- as.numeric(as.character(sim.all[,j]))

y.max <- max(c(obs.all$y, obs.all$y.scaled, sim.all$y, sim.all$y.scaled),na.rm=T)
# y.max<-100
# obs.all$COV = paste(obs.all$SCHED,obs.all$ALB)
# sim.all$COV = paste(sim.all$SCHED,sim.all$ALB)

# ----------------- VPC's -----------------------
#define prediction interval based on selection above
interval <- c(0.05, 0.5, 0.95)
if (pred.int == "80perc") {
  interval <- c(0.1, 0.5, 0.9)
}

# plot numbers
pnr <- 0

# ------------------- Just obs data -----------------

pnr <- pnr + 1
p.data.unstratified <- sim.all %>% 
  ggplot(aes(x=x,y=y,col=factor(ucov)))+
  geom_vline(xintercept = all.breaks, col='gray80')+
  geom_point()+
  geom_point(data=obs.all,col='black',size=1)+
  scale_y_continuous(limits=c(0,y.max))+
  ggtitle(paste0(pnr,". All data, not stratified"))+
  guides(color=guide_legend(title="Covariate"))

# Normal VPC (first test ref only)
vpc.unstratified <- observed(obs.all, x=x, y=y) %>%
  simulated(sim.all, x=x, y=y) %>%
  #stratify(~COV) %>%
  binning(bin = "centers", centers=ref.breaks) %>%
  vpcstats(qpred = interval)
vpc.unstratified$stats

# Normal VPC - ref only
pnr <- pnr + 1
p.vpc.unstratified <- ggplot(vpc.unstratified$stats, aes(x = xbin)) + 
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = qname, col = qname, group = qname),alpha = 0.25, col = NA) + 
  geom_line(aes(y = md, col = qname, group = qname)) +
  geom_line(aes(y = y, linetype = qname), size = 1) +
  #facet_wrap(~COV)+
  scale_y_continuous(limits=c(0,y.max))+
  scale_color_manual(values = c("blue", "red", "blue")) +
  scale_fill_manual(values = c("blue", "red", "blue")) +
  labs(title=paste0(pnr,". Normal VPC. Not Stratified"),
       color="Simulated",
       fill="CI95",
       linetype="Observed") 

# -----------  General stratified --------------------

pnr <- pnr + 1
p.data.stratified <- sim.all %>% 
  ggplot(aes(x=x,y=y,col=factor(ucov)))+
  geom_vline(xintercept = all.breaks, col='gray80')+
  geom_point()+
  geom_point(data=obs.all,col='black',size=1)+
  facet_wrap(~COV)+
  scale_y_continuous(limits=c(0,y.max))+
  ggtitle(paste0(pnr,". Stratified"))+
  guides(color=guide_legend(title="Covariate"))

# Normal VPC (first test ref only)
vpc.stratified <- observed(obs.all, x=x, y=y) %>%
  simulated(sim.all, x=x, y=y) %>%
  stratify(~ucov) %>%
  #binning(bin = "breaks", breaks=all.breaks) %>%
  binning(bin = "centers", centers=ref.breaks) %>%
  vpcstats(qpred = interval)
vpc.stratified$stats

# Normal VPC - ref only
pnr <- pnr + 1
p.vpc.stratified <- ggplot(vpc.stratified$stats, aes(x = xbin)) + 
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = qname, col = qname, group = qname),alpha = 0.25, col = NA) + 
  geom_line(aes(y = md, col = qname, group = qname)) +
  geom_line(aes(y = y, linetype = qname), size = 1) +
  facet_wrap(~ucov)+
  scale_y_continuous(limits=c(0,y.max))+
  scale_color_manual(values = c("blue", "red", "blue")) +
  scale_fill_manual(values = c("blue", "red", "blue")) +
  labs(title=paste0(pnr,". Normal VPC. Stratified"),
       color="Simulated",
       fill="CI95",
       linetype="Observed")

# ------------- pcVPC -------------------------------------

# pcVPC (first test ref only)
pnr <- pnr + 1
pcvpc.unstratified <- observed(obs.all, x=x, y=y) %>%
  simulated(sim.all, x=x, y=y) %>%
  binning(bin = "centers", centers=ref.breaks) %>%
  predcorrect(pred=PRED) %>%
  vpcstats()
pcvpc.unstratified$stats

#plot of prediction-corrected data points
pcvpc.data.stratified <- pcvpc.unstratified$sim %>%
  mutate(ucov = rep(obs.all[,"ucov"], max(sim.all$isim))) %>%
  ggplot(aes(x=x,y=ypc,col=factor(ucov)))+
  geom_vline(xintercept = ref.breaks, col='gray80')+
  geom_point()+
  geom_point(data=pcvpc.unstratified$obs,col='black',size=1)+
  scale_y_continuous(limits=c(0,y.max))+
  ggtitle(paste0(pnr,". Prediction-corrected data (all)"))+
  guides(color=guide_legend(title="Covariate"))

pnr <- pnr + 1
# pcVPC - ref only
pc.vpc.unstratified <- ggplot(pcvpc.unstratified$stats, aes(x = xbin)) + 
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = qname, col = qname, group = qname),alpha = 0.25, col = NA) + 
  geom_line(aes(y = md, col = qname, group = qname)) +
  geom_line(aes(y = y, linetype = qname), size = 1) +
  scale_y_continuous(limits=c(0,y.max))+
  scale_color_manual(values = c("blue", "red", "blue")) +
  scale_fill_manual(values = c("blue", "red", "blue")) +
  labs(title=paste0(pnr,". pcVPC. Unstratified"),
       color="Simulated",
       fill="CI95",
       linetype="Observed") 


# ----------- All data Vachette transformed ----------------

pnr <- pnr + 1
p.data.vachette <- sim.all %>% 
  ggplot(aes(x=x.scaled,y=y.scaled,col=factor(region)))+
  geom_vline(xintercept = all.breaks, col='gray80')+
  geom_point()+
  geom_point(data=obs.all,col='black',size=2)+
  guides(color=guide_legend(title="Region"))+
  scale_y_continuous(limits=c(0,y.max))+
  ggtitle(paste0(pnr,". Vachette transformed data (all) to true reference"))

vachette.vpc <- observed(obs.all,x=x.scaled,y=y.scaled) %>%
  simulated(sim.all,y=y.scaled) %>%
  #binning(bin = x.scaled) %>%
  binning(bin = "centers", centers=all.breaks) %>%
  vpcstats(qpred = interval)
vachette.vpc$stats


# Normal VPC - ref only
pnr <- pnr + 1
p.vpc.vachette <- ggplot(vachette.vpc$stats, aes(x = xbin)) + 
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = qname, col = qname, group = qname),alpha = 0.25, col = NA) + 
  geom_line(aes(y = md, col = qname, group = qname)) +
  geom_line(aes(y = y, linetype = qname), size = 1) +
  scale_y_continuous(limits=c(0,y.max))+
  scale_color_manual(values = c("blue", "red", "blue")) +
  scale_fill_manual(values = c("blue", "red", "blue")) +
  labs(title=paste0(pnr,". Vachette-VPC"),
       color="Simulated",
       fill="CI95",
       linetype="Observed")
  

# -------------------------------------------------------

# Save plots in pdf-file:
pdf(paste0("../plots-james/vachette-vpc-",model,"-",tag,".pdf"))

plot(p.data.unstratified + render)
plot(p.vpc.unstratified + render)

plot(p.data.stratified + render)
plot(p.vpc.stratified + render)

plot(pcvpc.data.stratified + render)
plot(pc.vpc.unstratified + render)

plot(p.data.vachette + render)
plot(p.vpc.vachette + render)

plot(p.data.unstratified + render+scale_y_log10())
plot(p.vpc.unstratified + render+scale_y_log10())

plot(p.data.stratified + render+scale_y_log10())
plot(p.vpc.stratified + render+scale_y_log10())

plot(pcvpc.data.stratified + render +scale_y_log10())
plot(pc.vpc.unstratified + render+scale_y_log10())

plot(p.data.vachette + render+scale_y_log10())
plot(p.vpc.vachette + render+scale_y_log10())

dev.off()

# --- end ---
