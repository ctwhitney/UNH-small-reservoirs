#Analysis and figures for 'Small Dams' manuscript
#Created 2019-11-19

#Sections:
#1 - Read in grab sample data




rm(list=ls())
#Set Working Directory
#setwd("C:/Users/ctw1/Box/Data/Analysis/Dams/")
Sys.Date<-as.character(Sys.Date())

library(gridExtra)
#library(reshape2)
#library(ggpubr)
library(tidyverse)


##### Define some functions #####

ProcessUSGSdaily <- function(x) {
  x <- x %>% mutate(
    "Date" = as.Date(strptime(x[,3],format="%Y-%m-%d")),
    "Year" = as.numeric(format(as.Date(strptime(x[,3],format="%Y-%m-%d")),"%Y")),
    "Discharge.m3s" =  x[,4]/35.3147,
    "Discharge.Ls" = x[,4]/35.3147/1000,
    "Date2" = lubridate::floor_date(as.Date(strptime(x[,3],format="%Y-%m-%d")),"month"))
  x <- x[,c(6,8,9,7,10)]
}


get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

lm_eqn2 = function(dx, y){
  m = df %>% select(y_var = !!enquo(y), x_var = !!enquo(x)) %>%
    lm(y_var ~ x_var, data = .);
  eq <- substitute(atop(italic(y) == a + b %.% italic(x),italic(r)^2~"="~r2), 
                   list(a = format(coef(summary(m))[1,1], digits = 2), 
                        b = format(coef(summary(m))[2,1], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 2)))
  as.character(as.expression(eq));                 
}


################################################################################
## Import data ##
################################################################################

## Read in Grab sample data
GRABS<-read.csv(file="Data/Reservoir_Grab_2022-05-06.csv",header=TRUE,stringsAsFactors=FALSE)
GRABS$DateTime<-as.POSIXct(strptime(paste(GRABS$Date,GRABS$Time),format="%m/%d/%Y %H:%M",tz="America/New_York"))
attributes(GRABS$DateTime)$tzone<-"EST"
GRABS$Date<-as.POSIXct(strptime(GRABS$Date,format="%m/%d/%Y"))
GRABS$Year<-as.numeric(format(GRABS$Date,"%Y"))
GRABS$DateTime15<-as.POSIXct(round(as.numeric(GRABS$DateTime)/900)*900,origin='1970-01-01',tz="EST")

## Remove Location == NA, Location == "IP" and Storm data from GRABS
GRABS<-GRABS[!is.na(GRABS$Location),]
GRABS<-GRABS[!grepl("_Storm",GRABS$Sample),]
GRABS<-GRABS[!GRABS$Sample=="SMD_BD",]
GRABS<-GRABS[!(GRABS$Site=="CCBP" | GRABS$Site=="DBBP" | GRABS$Site=="IP"),]

## Add in estimated catchment area for SMD_SIDE (closest pixel in StreamStats)
GRABS[GRABS$Sample=="SMD_SIDE",]$CArea.km2<-113.1825

# Calculate N2:Ar disequilibrium where available
GRABS$NArDisEq<-GRABS$NArcalc.mean-GRABS$NAr_sat

# Calculate DOC/DON ratio
GRABS$DOCtoDON<-GRABS$DOC.mgL/GRABS$DON.mgL
# Calculate DOC/TDN ratio
GRABS$DOCtoTDN<-GRABS$DOC.mgL/GRABS$DON.mgL

#####
## Create Seitzinger et al. 2006 and David et al 2006 R vs HL relationships
#####
LitRemoval<-data.frame(
  Hl=seq(1.1,10000,length.out = 100)
)
LitRemoval$Removal_Seit<-0.88453*LitRemoval$Hl^-0.3677
LitRemoval$Removal_Seit100 <- -88.453*LitRemoval$Hl^-0.3677
LitRemoval$Removal_Seit<-LitRemoval$Removal_Seit+rnorm(length(LitRemoval$Removal_Seit),sd=0.01)
LitRemoval$Removal_Seit100<-LitRemoval$Removal_Seit100+rnorm(length(LitRemoval$Removal_Seit100),sd=0.01)

Seit <- nls(Removal_Seit100~a*Hl^b, LitRemoval, start=list(a=-88, b=-0.3))

LitRemoval$Removal_David<-2.43*LitRemoval$Hl^-0.5632
LitRemoval$Removal_David<-LitRemoval$Removal_David+rnorm(length(LitRemoval$Removal_David),sd=0.01)

LitRemoval <- LitRemoval %>% pivot_longer(!Hl)

#LitRemoval<-reshape2::melt(LitRemoval,id.vars=c("Hl"))
LitRemoval$Eq<-NA
LitRemoval[LitRemoval$name=="Removal_Seit",]$Eq<-"Seitzinger et al."
LitRemoval[LitRemoval$name=="Removal_David",]$Eq<-"David et al."

Seit <- nls(value~a*Hl^b,LitRemoval[LitRemoval$name == "Removal_Seit",], start=list(a=100,b=-1))
David <- nls(value~a*Hl^b,LitRemoval[LitRemoval$name == "Removal_David",], start=list(a=100,b=-1))

#####
################ Read in discharge data from USGS website ################
#####

##USGS has a package "dataRetrieval" that can download data - look into this

## Can skip the below downloading/.processing steps and import data from file:

Q_Daily <- read.csv("Data/Q_Daily_2022-06-08.csv", stringsAsFactors = FALSE, header = TRUE)
Q_Daily$Date <- as.Date(Q_Daily$Date, format = "%Y-%m-%d")

## Start and end dates of all grab sample data
Startdate<-format(GRABS$Date[1],"%Y-%m-%d")
Enddate<-format(GRABS$Date[length(GRABS$Date)],"%Y-%m-%d")

## Assign gage numbers to gages
SMDGage<-"01101500"
IDGage<-"01102000"
PDGage<-"01101000"
OYSGage<-"01073000"
LMPGage<-"01073500"

## Download daily discharge data by site
SMD_Q_Daily<-read.table(file=paste0("https://waterdata.usgs.gov/nwis/dv?cb_00060=on&format=rdb&site_no=",SMDGage,"&referred_module=sw&period=&begin_date=",Startdate,"&end_date=",Enddate),skip=29,sep="\t",stringsAsFactors=FALSE)
ID_Q_Daily<-read.table(file=paste0("https://waterdata.usgs.gov/nwis/dv?cb_00060=on&format=rdb&site_no=",IDGage,"&referred_module=sw&period=&begin_date=",Startdate,"&end_date=",Enddate),skip=29,sep="\t",stringsAsFactors=FALSE)
PD_Q_Daily<-read.table(file=paste0("https://waterdata.usgs.gov/nwis/dv?cb_00060=on&format=rdb&site_no=",PDGage,"&referred_module=sw&period=&begin_date=",Startdate,"&end_date=",Enddate),skip=29,sep="\t",stringsAsFactors=FALSE)
OYS_Q_Daily<-read.table(file=paste0("https://waterdata.usgs.gov/nwis/dv?cb_00060=on&format=rdb&site_no=",OYSGage,"&referred_module=sw&period=&begin_date=",Startdate,"&end_date=",Enddate),skip=29,sep="\t",stringsAsFactors=FALSE)
LMP_Q_Daily<-read.table(file=paste0("https://waterdata.usgs.gov/nwis/dv?cb_00060=on&format=rdb&site_no=",LMPGage,"&referred_module=sw&period=&begin_date=",Startdate,"&end_date=",Enddate),skip=29,sep="\t",stringsAsFactors=FALSE)

## Process daily USGS data
SMD_Q_Daily <- ProcessUSGSdaily(SMD_Q_Daily)
ID_Q_Daily <- ProcessUSGSdaily(ID_Q_Daily)
PD_Q_Daily <- ProcessUSGSdaily(PD_Q_Daily)
OYS_Q_Daily <- ProcessUSGSdaily(OYS_Q_Daily)
LMP_Q_Daily <- ProcessUSGSdaily(LMP_Q_Daily)

## Scale USGS gage data to individual sampling locations

# SMD_UP, SMD_OUT
SMD_Q<-data.frame(Date=SMD_Q_Daily$Date,Discharge.m3s=((111.9/115.2545)*SMD_Q_Daily$Discharge.m3s),Site="SMD", Sample="SMD_UP")
SMD_Q<-rbind(SMD_Q,data.frame(Date=SMD_Q_Daily$Date,Discharge.m3s=((113.8/115.2545)*SMD_Q_Daily$Discharge.m3s),Site="SMD", Sample="SMD_OUT"))
# ID_UP, MR, ID
ID_Q<-data.frame(Date=ID_Q_Daily$Date,Discharge.m3s=((337.8/323.749)*ID_Q_Daily$Discharge.m3s),Site="ID",Sample="ID_UP")
ID_Q<-rbind(ID_Q,data.frame(Date=ID_Q_Daily$Date,Discharge.m3s=((43.3/323.749)*ID_Q_Daily$Discharge.m3s),Site="ID",Sample="MR"))
ID_Q<-rbind(ID_Q,data.frame(Date=ID_Q_Daily$Date,Discharge.m3s=((387.5/323.749)*ID_Q_Daily$Discharge.m3s),Site="ID",Sample="ID"))
# PD_UP, WB, PD
PD_Q<-data.frame(Date=PD_Q_Daily$Date,Discharge.m3s=((55.7/55.16675)*PD_Q_Daily$Discharge.m3s),Site="PD",Sample="PD_UP")
PD_Q<-rbind(PD_Q,data.frame(Date=PD_Q_Daily$Date,Discharge.m3s=((6/55.16675)*PD_Q_Daily$Discharge.m3s),Site="PD",Sample="WB"))
PD_Q<-rbind(PD_Q,data.frame(Date=PD_Q_Daily$Date,Discharge.m3s=((63.9/55.16675)*PD_Q_Daily$Discharge.m3s),Site="PD",Sample="PD"))
# LMPD_UP, LMPD
LMPD_Q<-data.frame(Date=LMP_Q_Daily$Date,Discharge.m3s=((471.2/479.148)*LMP_Q_Daily$Discharge.m3s),Site="LMPD",Sample="LMPD_IN")
LMPD_Q<-rbind(LMPD_Q,data.frame(Date=LMP_Q_Daily$Date,Discharge.m3s=((475.3/479.148)*LMP_Q_Daily$Discharge.m3s),Site="LMPD",Sample="LMPD"))
#MCLN_IN, Sidetrib_IN, Picassic, MCLN_OUT
MCLN_Q<-data.frame(Date=LMP_Q_Daily$Date,Discharge.m3s=((480.038/479.148)*LMP_Q_Daily$Discharge.m3s),Site="MCLN",Sample="MCLN_IN")
MCLN_Q<-rbind(MCLN_Q,data.frame(Date=LMP_Q_Daily$Date,Discharge.m3s=((1.735233/479.148)*LMP_Q_Daily$Discharge.m3s),Site="MCLN",Sample="Sidetrib_IN"))
MCLN_Q<-rbind(MCLN_Q,data.frame(Date=LMP_Q_Daily$Date,Discharge.m3s=((56.9778/479.148)*LMP_Q_Daily$Discharge.m3s),Site="MCLN",Sample="Piscassic"))
MCLN_Q<-rbind(MCLN_Q,data.frame(Date=LMP_Q_Daily$Date,Discharge.m3s=((549.0588/479.148)*LMP_Q_Daily$Discharge.m3s),Site="MCLN",Sample="MCLN_OUT"))
#LH_IN, LH_OUT
LH_Q<-data.frame(Date=OYS_Q_Daily$Date,Discharge.m3s=((1.061859/31.33886)*OYS_Q_Daily$Discharge.m3s),Site="LH",Sample="LH_IN")
LH_Q<-rbind(LH_Q,data.frame(Date=OYS_Q_Daily$Date,Discharge.m3s=((1.113657/31.33886)*OYS_Q_Daily$Discharge.m3s),Site="LH",Sample="LH_OUT"))
#BRDS_UP, PTEB, BRDS_OUT
BRDS_Q<-data.frame(Date=OYS_Q_Daily$Date,Discharge.m3s=((4.9014/31.33886)*OYS_Q_Daily$Discharge.m3s),Site="BRDS",Sample="BRDS")
BRDS_Q<-rbind(BRDS_Q,data.frame(Date=OYS_Q_Daily$Date,Discharge.m3s=((2.51/31.33886)*OYS_Q_Daily$Discharge.m3s),Site="BRDS",Sample="PTEB"))
BRDS_Q<-rbind(BRDS_Q,data.frame(Date=OYS_Q_Daily$Date,Discharge.m3s=((8.32/31.33886)*OYS_Q_Daily$Discharge.m3s),Site="BRDS",Sample="BRDS_OUT"))
#OMPD_IN, CLGB, HAM
OMPD_Q<-data.frame(Date=OYS_Q_Daily$Date,Discharge.m3s=((43.17363/31.33886)*OYS_Q_Daily$Discharge.m3s),Site="OMPD",Sample="OMPD_IN")
OMPD_Q<-rbind(OMPD_Q,data.frame(Date=OYS_Q_Daily$Date,Discharge.m3s=((2.006/31.33886)*OYS_Q_Daily$Discharge.m3s),Site="OMPD",Sample="CLGB"))
OMPD_Q<-rbind(OMPD_Q,data.frame(Date=OYS_Q_Daily$Date,Discharge.m3s=((1.7364/31.33886)*OYS_Q_Daily$Discharge.m3s),Site="OMPD",Sample="HAM"))
OMPD_Q<-rbind(OMPD_Q,data.frame(Date=OYS_Q_Daily$Date,Discharge.m3s=((50.59/31.33886)*OYS_Q_Daily$Discharge.m3s),Site="OMPD",Sample="OMPD"))


## Rbind individual site discharge records into one record
Q_Daily<-rbind(SMD_Q,ID_Q,PD_Q,LMPD_Q,MCLN_Q,LH_Q,BRDS_Q,OMPD_Q)

## Remove individual USGS discharge data
rm(ID_Q_Daily,SMD_Q_Daily,PD_Q_Daily,OYS_Q_Daily,LMP_Q_Daily)

## Remove site-specific discharge data
rm(SMD_Q,ID_Q,PD_Q,LMPD_Q,MCLN_Q,LH_Q,BRDS_Q,OMPD_Q)

## Write Q_Daily to file so it can be accessed later without re-downloading from USGS
write.csv(Q_Daily, file="Data/Q_Daily_2022-06-08.csv", row.names = FALSE)

# Merge GRABS with daily Q data
GRABS <- GRABS %>% left_join(., Q_Daily, by=c("Date", "Site", "Sample"))


#############################################################
##### Effective Discharge Analysis #####
#############################################################

## This may be able to go here in the future but for now this is a standalone script




############# Other processing on GRABS data ################

#SMD - 60809 including buffered channel above impoundment #Was 60751
#CCBP - 6900
#ID - 129610 including all wetland areas #Was 129493
#IP - 20288 #Was 15000 based on estimates form Google Maps
#PD - 162161 including all wetland areas #Was 64077
#BRDS - 61656 #Was 43363
#LH - 3607 #Was 3583
#OMPD - 91525 #Was 79655
#LMPD - #194560 including wetland area, 102773 without wetland #Was 121400
#MCLN - 526548 #Was 459829


GRABS$SArea.km2<-NA
#GRABS[GRABS$Site=="IP",]$SArea.km2<-20288/1e6
GRABS[GRABS$Site=="SMD",]$SArea.km2<-60809/1e6
GRABS[GRABS$Site=="ID",]$SArea.km2<-129610/1e6
GRABS[GRABS$Site=="PD",]$SArea.km2<-162161/1e6
## GRABS[GRABS$Site=="CCBP",]$SArea.km2<-6900/1e6
GRABS[GRABS$Site=="MCLN",]$SArea.km2<-526548/1e6
GRABS[GRABS$Site=="LMPD",]$SArea.km2<-194560/1e6
GRABS[GRABS$Site=="OMPD",]$SArea.km2<-91525/1e6
GRABS[GRABS$Site=="BRDS",]$SArea.km2<-61656/1e6
GRABS[GRABS$Site=="LH",]$SArea.km2<-3607/1e6

#Weight nutrient concentrations by catchment area for calculating removal
GRABS$Sp_CondCA<-GRABS$Sp_Cond.uScm*GRABS$CArea.km2
GRABS$TempCA<-GRABS$Temp.C*GRABS$CArea.km2
GRABS$DO.mgLCA<-GRABS$DO.mgL*GRABS$CArea.km2
GRABS$DO.pctCA<-GRABS$DO.pct*GRABS$CArea.km2
GRABS$TSSCA<-GRABS$TSS.mgL*GRABS$CArea.km2
GRABS$ChlCA<-GRABS$Chl.mgL*GRABS$CArea.km2
GRABS$PO4CA<-GRABS$PO4.ugL*GRABS$CArea.km2
GRABS$NH4CA<-GRABS$NH4.ugL*GRABS$CArea.km2
GRABS$ClCA<-GRABS$Cl.mgL*GRABS$CArea.km2
GRABS$NO3CA<-GRABS$NO3.mgL*GRABS$CArea.km2
GRABS$SO4CA<-GRABS$SO4.mgL*GRABS$CArea.km2
GRABS$BrCA<-GRABS$Br.mgL*GRABS$CArea.km2
GRABS$DOCCA<-GRABS$DOC.mgL*GRABS$CArea.km2
GRABS$TDNCA<-GRABS$TDN.mgL*GRABS$CArea.km2
GRABS$DONCA<-GRABS$DON.mgL*GRABS$CArea.km2
GRABS$DINCA<-GRABS$DIN.mgL*GRABS$CArea.km2
GRABS$TDPCA<-GRABS$TDP.mgL*GRABS$CArea.km2
GRABS$N2OCA<-GRABS$N2O.uatm*GRABS$CArea.km2
GRABS$CH4CA<-GRABS$CH4.uatm*GRABS$CArea.km2
#GRABS$CO2CA<-GRABS$CO2.ppm*GRABS$CArea.km2
GRABS$AFDM.mgLCA<-GRABS$AFDM.mgL*GRABS$CArea.km2
GRABS$Ash.mgLCA<-GRABS$Ash.mgL*GRABS$CArea.km2
GRABS$PCCA<-GRABS$PC.mgL*GRABS$CArea.km2
GRABS$PNCA<-GRABS$PN.mgL*GRABS$CArea.km2
GRABS$TNCA<-GRABS$TN.mgL*GRABS$CArea.km2
GRABS$NArCA<-GRABS$NArcalc.mean*GRABS$CArea.km2
GRABS$N2CA<-GRABS$N2.mM.mean*GRABS$CArea.km2
#GRABS$DOCtoDONCA<-GRABS$DOCtoDON*GRABS$CArea.km2 ## Don't need to weight the ratios themselves
#GRABS$DOCtoTDNCA<-GRABS$DOCtoTDN*GRABS$CArea.km2

# Flow weight nutrient concentrations (flux - kg/km2/day)
GRABS$TSS.kgkm2Day<-(GRABS$TSS.mgL*(GRABS$Discharge.m3s*1000)*86400/1e+6)
GRABS$Chl.kgkm2Day<-(GRABS$Chl.mgL*(GRABS$Discharge.m3s*1000)*86400/1e+6)
GRABS$PO4.kgkm2Day<-((GRABS$PO4.ugL*1000)*(GRABS$Discharge.m3s*1000)*86400/1e+6)
GRABS$NH4.kgkm2Day<-((GRABS$NH4.ugL*1000)*(GRABS$Discharge.m3s*1000)*86400/1e+6)
GRABS$Cl.kgkm2Day<-(GRABS$Cl.mgL*(GRABS$Discharge.m3s*1000)*86400/1e+6)
GRABS$NO3.kgkm2Day<-(GRABS$NO3.mgL*(GRABS$Discharge.m3s*1000)*86400/1e+6)
GRABS$SO4.kgkm2Day<-(GRABS$SO4.mgL*(GRABS$Discharge.m3s*1000)*86400/1e+6)
GRABS$Br.kgkm2Day<-(GRABS$Br.mgL*(GRABS$Discharge.m3s*1000)*86400/1e+6)
GRABS$DOC.kgkm2Day<-(GRABS$DOC.mgL*(GRABS$Discharge.m3s*1000)*86400/1e+6)
GRABS$TDN.kgkm2Day<-(GRABS$TDN.mgL*(GRABS$Discharge.m3s*1000)*86400/1e+6)
GRABS$DON.kgkm2Day<-(GRABS$DON.mgL*(GRABS$Discharge.m3s*1000)*86400/1e+6)
GRABS$DIN.kgkm2Day<-(GRABS$DIN.mgL*(GRABS$Discharge.m3s*1000)*86400/1e+6)
GRABS$TDP.kgkm2Day<-(GRABS$TDP.mgL*(GRABS$Discharge.m3s*1000)*86400/1e+6)
#GRABS$N2O.kgkm2Day<-(GRABS$N2O.uatm*(GRABS$Discharge.m3s*1000)*86400/1e+6)
#GRABS$CH4.kgkm2Day<-(GRABS$CH4.uatm*(GRABS$Discharge.m3s*1000)*86400/1e+6)
#GRABS$CO2CA<-(GRABS$CO2.ppm*(GRABS$Discharge.m3s*1000)*86400/1e+6)
GRABS$AFDM.kgkm2Day<-(GRABS$AFDM.mgL*(GRABS$Discharge.m3s*1000)*86400/1e+6)
GRABS$Ash.kgkm2Day<-(GRABS$Ash.mgL*(GRABS$Discharge.m3s*1000)*86400/1e+6)
GRABS$PC.kgkm2Day<-(GRABS$PC.mgL*(GRABS$Discharge.m3s*1000)*86400/1e+6)
GRABS$PN.kgkm2Day<-(GRABS$PN.mgL*(GRABS$Discharge.m3s*1000)*86400/1e+6)
GRABS$TN.kgkm2Day<-(GRABS$TN.mgL*(GRABS$Discharge.m3s*1000)*86400/1e+6)
#GRABS$DOCtoDON.kgkm2Day <- GRABS$DOC.kgkm2Day/GRABS$DON.kgkm2Day
#GRABS$DOCtoTDN.kgkm2Day <- GRABS$DOC.kgkm2Day/GRABS$TDN.kgkm2Day
#GRABS$NArCA<-(GRABS$NArcalc.mean*(GRABS$Discharge.m3s*1000)*86400/1e+6)
#GRABS$N2CA<-(GRABS$N2.mM.mean*(GRABS$Discharge.m3s*1000)*86400/1e+6)
#GRABS$DOCtoDON.kgkm2Day2<-(GRABS$DOCtoDON*(GRABS$Discharge.m3s*1000)*86400/1e+6)
#GRABS$DOCtoTDN.kgkm2Day<-(GRABS$DOCtoTDN*(GRABS$Discharge.m3s*1000)*86400/1e+6)



## Separate data by site ##

Removal<-GRABS[grep("UP",GRABS$Location),] #c(1:80,83)

## Something probably changed here so the original columns selected here no longer work??
#Removal<-aggregate(Removal[,c(8:14,16:26,28:32,34:41,45,49,54:106)], by=list(Removal$Site,Removal$Date),FUN=sum,na.rm=TRUE) #c(8:14,16:26,28:32,34:41,45,49,54:79)
## ALSO!!! Changed code so DOCtoDONCA and DOCtoDON.kgkm2Day were not calculated. **54:106 -> 54:102**
## Can't do these things with ratios - calculate ratios from DOCCA, DONCA, DOC.kgkm2Day, and DON.kgkm2Day
## ALSO, got rid of DOCtoDON and DOCtoTDN in Removal df (cols 55 and 56) - can't sum ratios!!

Removal<-aggregate(Removal[,c(6,8:14,16:26,28:32,34:41,45,49,54,57:102)], by=list(Removal$Site,Removal$Date),FUN=sum) #c(8:14,16:26,28:32,34:41,45,49,54:79)
colnames(Removal)[1:2]<-c("Site","Date")
#Removal$Location<-'UP'
Removal<-merge(Removal,GRABS[GRABS$Location == "OUT",],by=c("Date","Site"),suffixes=c(".UP",".OUT"),all=TRUE) #c(1:80,83)

## Add DOC:DON into Removal df
Removal$DOCtoDON.UP <- Removal$DOCCA.UP / Removal$DONCA.UP
Removal$DOCtoDON.OUT <- Removal$DOCCA.OUT / Removal$DONCA.OUT
Removal$DOCtoTDN.UP <- Removal$DOCCA.UP / Removal$TDNCA.UP
Removal$DOCtoTDN.OUT <- Removal$DOCCA.OUT / Removal$TDNCA.OUT

##### Calculate removal for all variables of interest #####

## Removal dataframe ##
Removal$RSpCond<-with(Removal,((Sp_CondCA.UP-Sp_CondCA.OUT)/Sp_CondCA.UP))
Removal$RTemp<-with(Removal,((TempCA.UP-TempCA.OUT)/TempCA.UP))
Removal$RDO.mgL<-with(Removal,((DO.mgLCA.UP-DO.mgLCA.OUT)/DO.mgLCA.UP))
Removal$RDO.pct<-with(Removal,((DO.pctCA.UP-DO.pctCA.OUT)/DO.pctCA.UP))
Removal$RTSS<-with(Removal,((TSSCA.UP-TSSCA.OUT)/TSSCA.UP))
Removal$RChl<-with(Removal,((ChlCA.UP-ChlCA.OUT)/ChlCA.UP))
Removal$RPO4<-with(Removal,((PO4CA.UP-PO4CA.OUT)/PO4CA.UP))
Removal$RNH4<-with(Removal,((NH4CA.UP-NH4CA.OUT)/NH4CA.UP))
Removal$RCl<-with(Removal,((ClCA.UP-ClCA.OUT)/ClCA.UP))
Removal$RNO3<-with(Removal,((NO3CA.UP-NO3CA.OUT)/NO3CA.UP))
Removal$RSO4<-with(Removal,((SO4CA.UP-SO4CA.OUT)/SO4CA.UP))
Removal$RBr<-with(Removal,((BrCA.UP-BrCA.OUT)/BrCA.UP))
Removal$RDOC<-with(Removal,((DOCCA.UP-DOCCA.OUT)/DOCCA.UP))
Removal$RTDN<-with(Removal,((TDNCA.UP-TDNCA.OUT)/TDNCA.UP))
Removal$RDON<-with(Removal,((DONCA.UP-DONCA.OUT)/DONCA.UP))
Removal$RDIN<-with(Removal,((DINCA.UP-DINCA.OUT)/DINCA.UP))
#Removal$RTDN2<-with(Removal,((TDNCA.OUT-TDNCA.UP)/TDNCA.OUT)) #CHANGED THIS
#Removal$RDON2<-with(Removal,((DONCA.OUT-DONCA.UP)/DONCA.OUT)) #CHANGED THIS
#Removal$RDIN2<-with(Removal,((DINCA.OUT-DINCA.UP)/DINCA.OUT)) #CHANGED THIS
Removal$RTDP<-with(Removal,((TDPCA.UP-TDPCA.OUT)/TDPCA.UP))
Removal$RN2O<-with(Removal,((N2OCA.UP-N2OCA.OUT)/N2OCA.UP))
Removal$RCH4<-with(Removal,((CH4CA.UP-CH4CA.OUT)/CH4CA.UP))
Removal$RAFDM<-with(Removal,((AFDM.mgLCA.UP-AFDM.mgLCA.OUT)/AFDM.mgLCA.UP))
Removal$RAsh<-with(Removal,((Ash.mgLCA.UP-Ash.mgLCA.OUT)/Ash.mgLCA.UP))
Removal$RPC<-with(Removal,((PCCA.UP-PCCA.OUT)/PCCA.UP))
Removal$RPN<-with(Removal,((PNCA.UP-PNCA.OUT)/PNCA.UP))
Removal$RTN<-with(Removal,((TNCA.UP-TNCA.OUT)/TNCA.UP))
Removal$RNAr<-with(Removal,((NArCA.UP-NArCA.OUT)/NArCA.UP))
Removal$RN2<-with(Removal,((N2CA.UP-N2CA.OUT)/N2CA.UP))
Removal$RDOCtoDON<-with(Removal,((DOCtoDON.UP-DOCtoDON.OUT)/DOCtoDON.UP))
Removal$RDOCtoTDN<-with(Removal,((DOCtoTDN.UP-DOCtoTDN.OUT)/DOCtoTDN.UP))


## Create flag for chloride mass balance ##
Removal$nInputs<-NA
Removal[Removal$Site=="SMD",]$nInputs<-1
Removal[Removal$Site=="ID",]$nInputs<-2
Removal[Removal$Site=="PD",]$nInputs<-2
Removal[Removal$Site=="LMPD",]$nInputs<-1
Removal[Removal$Site=="MCLN",]$nInputs<-3
Removal[Removal$Site=="LH",]$nInputs<-1
Removal[Removal$Site=="BRDS",]$nInputs<-2
Removal[Removal$Site=="OMPD",]$nInputs<-3

Removal$RFlag<-ifelse((Removal$RCl>(Removal$nInputs/10) | Removal$RCl<(-Removal$nInputs/10)),1,0)
# Create second RFlag set at ±20%
Removal$RFlag2<-ifelse((Removal$RCl>(0.2) | Removal$RCl<(-0.2)),1,0)

## Manually flag certain values
Removal[87,c("RFlag","RFlag2")]<-2 # DIN Removal = -9.38, no flow over dam
Removal[93,c("RFlag","RFlag2")]<-2 # No flow over dam
Removal[105,c("RFlag","RFlag2")]<-2 # No flow over dam
Removal[126,c("RFlag","RFlag2")]<-2 # Massive NH4 production (RNH4 = -6.11) - not sure why
#Removal[146,c("RFlag","RFlag2")]<-2 # ID 2018-05-08 DIN production > 100% - not flagged because DON, TDN are OK
Removal[147,c("RFlag","RFlag2")]<-2 # NO3, NH4 production (inputs BDL, output above DL); little change in TDN, DON
Removal[134,c("RFlag","RFlag2")]<-3 # ID_OUT does not exist (ice-no sample) and causing problems


#Remove infinite values
Removal<-do.call(data.frame,lapply(Removal, function(x) replace(x, is.infinite(x),NA)))

# Calculate change in flux or flow-weighted concentrations
# Flow weight nutrient concentrations (flux - kg/day)
Removal$Delta_TSS.kgkm2Day<-Removal$TSS.kgkm2Day.UP-Removal$TSS.kgkm2Day.OUT
Removal$Delta_Chl.kgkm2Day<-Removal$Chl.kgkm2Day.UP-Removal$Chl.kgkm2Day.OUT
Removal$Delta_PO4.kgkm2Day<-Removal$PO4.kgkm2Day.UP-Removal$PO4.kgkm2Day.OUT
Removal$Delta_NH4.kgkm2Day<-Removal$NH4.kgkm2Day.UP-Removal$NH4.kgkm2Day.OUT
Removal$Delta_Cl.kgkm2Day<-Removal$Cl.kgkm2Day.UP-Removal$Cl.kgkm2Day.OUT
Removal$Delta_NO3.kgkm2Day<-Removal$NO3.kgkm2Day.UP-Removal$NO3.kgkm2Day.OUT
Removal$Delta_SO4.kgkm2Day<-Removal$SO4.kgkm2Day.UP-Removal$SO4.kgkm2Day.OUT
Removal$Delta_Br.kgkm2Day<-Removal$Br.kgkm2Day.UP-Removal$Br.kgkm2Day.OUT
Removal$Delta_DOC.kgkm2Day<-Removal$DOC.kgkm2Day.UP-Removal$DOC.kgkm2Day.OUT
Removal$Delta_TDN.kgkm2Day<-Removal$TDN.kgkm2Day.UP-Removal$TDN.kgkm2Day.OUT
Removal$Delta_DON.kgkm2Day<-Removal$DON.kgkm2Day.UP-Removal$DON.kgkm2Day.OUT
Removal$Delta_DIN.kgkm2Day<-Removal$DIN.kgkm2Day.UP-Removal$DIN.kgkm2Day.OUT
Removal$Delta_TDP.kgkm2Day<-Removal$TDP.kgkm2Day.UP-Removal$TDP.kgkm2Day.OUT
#Removal$Delta_N2O.kgkm2Day<-
#Removal$Delta_CH4.kgkm2Day<-
#Removal$Delta_CO2CA<-
Removal$Delta_AFDM.kgkm2Day<-Removal$AFDM.kgkm2Day.UP-Removal$AFDM.kgkm2Day.OUT
Removal$Delta_Ash.kgkm2Day<-Removal$Ash.kgkm2Day.UP-Removal$Ash.kgkm2Day.OUT
Removal$Delta_PC.kgkm2Day<-Removal$PC.kgkm2Day.UP-Removal$PC.kgkm2Day.OUT
Removal$Delta_PN.kgkm2Day<-Removal$PN.kgkm2Day.UP-Removal$PN.kgkm2Day.OUT
Removal$Delta_TN.kgkm2Day<-Removal$TN.kgkm2Day.UP-Removal$TN.kgkm2Day.OUT
#Removal$Delta_NArCA<-
#Removal$Delta_N2CA<-
#Removal$Delta_DOCtoDON.kgkm2Day<-Removal$DOCtoDON.kgkm2Day.UP-Removal$DOCtoDON.kgkm2Day.OUT
#Removal$Delta_DOCtoTDN.kgkm2Day<-Removal$DOCtoTDN.kgkm2Day.UP-Removal$DOCtoTDN.kgkm2Day.OUT

## Calculate percent removal of flow-weighted fluxes
Removal$RDIN.kgkm2Day <- with(Removal,(DIN.kgkm2Day.UP-DIN.kgkm2Day.OUT)/DIN.kgkm2Day.UP)
Removal$RDON.kgkm2Day <- with(Removal,(DON.kgkm2Day.UP-DON.kgkm2Day.OUT)/DON.kgkm2Day.UP)
Removal$RTDN.kgkm2Day <- with(Removal,(TDN.kgkm2Day.UP-TDN.kgkm2Day.OUT)/TDN.kgkm2Day.UP)


## Calculate hydraulic load for each site/date ##
Removal$HL.myr<-(Removal$Discharge.m3s.OUT/(Removal$SArea.km2.OUT*1e+6))*86400*365

## Add some other grouping variables to the data frame
# Month
Removal$Month<-format(Removal$Date,"%B")

# Season
for(i in 1:nrow(Removal)){
  Removal$Season[as.numeric(format(Removal$Date,"%m"))%in%c(1,2,12)]<-"Winter"
  Removal$Season[as.numeric(format(Removal$Date,"%m"))%in%c(3:5)]<-"Spring"
  Removal$Season[as.numeric(format(Removal$Date,"%m"))%in%c(6:8)]<-"Summer"
  Removal$Season[as.numeric(format(Removal$Date,"%m"))%in%c(9:11)]<-"Autumn"
}

# Day of year (numeric)
Removal$DOY<-as.POSIXlt(Removal$Date)$yday+1

##### Calculate mean annual RDIN #####
# Discharge values are averages of mean annual flow from 2015-2019 for each gage (Q.m3s_5avg),
# scaled to each site using catchment area (similar to above). Surface areas for calculating
# hydraulic load are taken from GRABS file.

MA_RDIN<-data.frame(Q.m3s_2019=c(2.4942,2.7177,8.6325,1.6433,10.8341,12.5154,0.0263,0.1962,1.1931),
                    Q.m3s_5avg=c(1.7547,1.9119,5.8605,1.1099,7.7499,8.9525,0.0189,0.1409,0.8566),
                    SA.km2=c(0.020288,0.060809,0.12961,0.162161,0.19456,0.526548,0.003607,0.061656,0.091525),
                    Site=c("IP","SMD","ID","PD","LMPD","MCLN","LH","BRDS","OMPD"))

MA_RDIN$HL.myr_2019<-(MA_RDIN$Q.m3s_2019/(MA_RDIN$SA.km2*1e+6)*86400*365)
MA_RDIN$HL.myr_5avg<-(MA_RDIN$Q.m3s_5avg/(MA_RDIN$SA.km2*1e+6)*86400*365)
RDINfit<-nls(RNO3~SSlogis(log10(HL.myr), Asym, xmid, scal),data=Removal[Removal$RFlag==0 & !Removal$Site=="IP",])
RDINfit<-nls(RDIN~a*HL.myr^b,Removal[Removal$HL.myr >= 1 & Removal$RFlag==0 & ! Removal$Site=="IP",],start=list(a=1,b=-1))
MA_RDIN$RDIN_2019<-predict(RDINfit,data.frame(HL.myr=MA_RDIN$HL.myr_2019))
MA_RDIN$RDIN_5avg<-predict(RDINfit,data.frame(HL.myr=MA_RDIN$HL.myr_5avg))
MA_RDIN$RDIN_Seit_2019<-(88*MA_RDIN$HL.myr_2019^-0.368)/100
MA_RDIN$RDIN_Seit_5avg<-(0.88*MA_RDIN$HL.myr_5avg^-0.368)


###################### Mean annual change ######################
MAChange<-data.frame(Site=c("MCLN","LMPD","ID","SMD","PD","OMPD","BRDS","LH"),
                     MAQ.m3s=c(9.43,8.16,6.57,1.88,1.25,0.92,0.15,0.02),
                     SA.km2=c(0.526548,0.194560,0.129610,0.060809,0.162161,0.091525,0.061656,0.003607),
                     CA.km2=c(549.1,475.3,387.5,113.8,63.9,50.6,8.3,1.1))
MAChange$HL.myr<-MAChange$MAQ.m3s/(MAChange$SA.km2*1e+6)*86400*365
RDINfitfpl<-nls(-RDIN*100~SSfpl(log10(HL.myr), A, B, xmid, scal),data=Removal[Removal$RFlag==0 & Removal$RDIN > -2,])
RDINfitlog<-lm(-RDIN*100~log10(HL.myr),Removal[Removal$RFlag==0 & Removal$RDIN > -2,])
RDINfit2<-lm(-RDIN~log10(HL.myr),Removal[Removal$RFlag==0,])
RDONfit<-lm(RDON*100~log10(HL.myr),Removal[Removal$RFlag==0,])
RTDNfit<-lm(-RTDN*100~log10(HL.myr),Removal[Removal$RFlag==0,])
MAChange$RDIN<-predict(RDINfitfpl,data.frame(HL.myr=MAChange$HL.myr))
MAChange$RDIN2<-predict(RDINfit2,data.frame(HL.myr=MAChange$HL.myr))
MAChange$RDIN3<-(MAChange$HL.myr*-0.024604+29.870802)
MAChange$RDON<-predict(RDONfit,data.frame(HL.myr=MAChange$HL.myr))
MAChange$RTDN<-predict(RTDNfit,data.frame(HL.myr=MAChange$HL.myr))

## Export MA Change data for Wil
write.csv(MAChange,file=paste0("Data/MA_Removal_Export_",Sys.Date(),".csv"),row.names=FALSE)


## RDIN fit comparison
RDINfitfpl<-nls(-RDIN*100~SSfpl(log10(HL.myr), A, B, xmid, scal),data=Removal[Removal$RFlag2==0 & Removal$RDIN > -2,])
RDINfitlog<-lm(-RDIN*100~log10(HL.myr),Removal[Removal$RFlag2==0 & Removal$RDIN > -2,])
RDINfitpow<-nls(-RDIN*100~a*log10(HL.myr)^b,Removal[Removal$RFlag2==0 & Removal$RDIN > -2 & Removal$HL.myr > 0.5,],start=list(a=1,b=-1))

AIC(RDINfitfpl,RDINfitlog,RDINfitpow)
BIC(RDINfitfpl,RDINfitlog,RDINfitpow)



###################### Mean removal/production df for synthesis ####################
Means<-Removal[Removal$RFlag==0,]
Means<-aggregate(Means[,c(182:226,135,227,86,90,91,137)],by=list(Means$Site),FUN=median,na.rm=TRUE)
names(Means)[1]<-"Site"
Means$Depth.m<-c(3.48,2.85,0.75,3.16,4.07,1.42,0.76,1.27)
Means$Age.yrs<-c(68,121,70,110,134,86,46,68)
Means$DAtoSA<-Means$CArea.km2/Means$SArea.km2.OUT
Means$Vol<-c(215859,370044,3000,616740,2146255,130749,123348,77709.2)
Means$Dev<-c(49.10,34.62,51.11,7.89,8.50,16.62,18.33,53.53)
Means$For<-c(37.84,35.11,42.88,70.34,68.22,59.16,50.09,22.58)
Means$Wet<-c(4.83,23.14,4.70,11.68,13.09,11.74,23.29,20.76)
Means$SAtoVol<-(Means$SArea.km2.OUT*1e+6)/Means$Vol
Means$ResT<-Means$Discharge.m3s.OUT/Means$Vol*86400
Means$Runoff.mmDay<-Means$Discharge.m3s.OUT/(Means$CArea.km2*1e+6)*86400*1000
#Not sure if this belongs here:
Means$MAF<-c(0.15,6.57,0.02,8.16,9.43,0.92,1.25,1.88)
Means$MAHL.myr<-(Means$MAF/(Means$SArea.km2.OUT*1e+6))*86400*365

Means$RDINcalc<-predict(MARDINfit,data.frame(HL.myr=Means$MAHL.myr))

Means$VfDIN.myr<-((-Means$HL.myr)*log(1-Means$RNO3))



###################### Final figures for the manuscript ###########
# Starting at Figure 2 because Figure 1 is the study site/map

# Figure 2 - DIN/DON proportions of TDN at all sites for inflow and outflow
# First, need to create melted dataframe from GRABS data
# Add new 'TN.mgL' column that uses TDN for samples where there is no PN

## Updated on 2021-08-10:: Changed sum to mean in aggregate calls not using proportion of TDN - 
## now using absolute concentrations with proportions of DIN/DON

GRABS$TN.mgL2<-rowSums(GRABS[,c("TDN.mgL","PN.mgL")],na.rm=TRUE)
GRABSmelt<-aggregate(GRABS[grep("UP",GRABS$Location),c(23:25)],list(GRABS[grep("UP",GRABS$Location),]$Site), mean, na.rm=TRUE) #Also previously included: ,39,40 (102)
GRABSmelt<-merge(GRABSmelt,aggregate(GRABS[GRABS$Location=="OUT",c(23:25)],list(GRABS[GRABS$Location=="OUT",]$Site), mean, na.rm=TRUE),by="Group.1",suffixes=c(".UP",".OUT")) #Also previously included: ,39,102
names(GRABSmelt)[1]<-"Site"
GRABSmelt$PDIN.UP<-GRABSmelt$DIN.mgL.UP/GRABSmelt$TDN.mgL.UP
GRABSmelt$PDON.UP<-GRABSmelt$DON.mgL.UP/GRABSmelt$TDN.mgL.UP
#GRABSmelt$PPN.UP<-GRABSmelt$PN.mgL.UP/GRABSmelt$TDN.mgL.UP
GRABSmelt$PDIN.OUT<-GRABSmelt$DIN.mgL.OUT/GRABSmelt$TDN.mgL.OUT
GRABSmelt$PDON.OUT<-GRABSmelt$DON.mgL.OUT/GRABSmelt$TDN.mgL.OUT
#GRABSmelt$PPN.OUT<-GRABSmelt$PN.mgL.OUT/GRABSmelt$TDN.mgL.OUT
GRABSmelt<-melt(GRABSmelt,id.vars="Site")
GRABSmelt$Location<-NA
GRABSmelt[grep("UP",GRABSmelt$variable),]$Location<-"Inflow"
GRABSmelt[grep("OUT",GRABSmelt$variable),]$Location<-"Outflow"
GRABSmelt$Location_f<-factor(GRABSmelt$Location,levels=c("Inflow","Outflow"))
GRABSmelt$Species<-NA
#GRABSmelt[grep("TDN",GRABSmelt$variable),]$Species<-"TDN"
GRABSmelt[grep("DIN",GRABSmelt$variable),]$Species<-"DIN"
GRABSmelt[grep("DON",GRABSmelt$variable),]$Species<-"DON"
#GRABSmelt[grep("PN",GRABSmelt$variable),]$Species<-"PN"
#GRABSmelt<-GRABSmelt[-grep(c("TN|mgL"),GRABSmelt$variable),]
# Sorted by decreasing surface area
GRABSmelt$Site_f<-factor(GRABSmelt$Site,levels=c("MCLN","LMPD","PD","ID","OMPD","BRDS","SMD","LH","IP"))
# Sorted by decreasing drainage area
GRABSmelt$Site_f2<-factor(GRABSmelt$Site,levels=c("MCLN","LMPD","ID","SMD","PD","OMPD","BRDS","LH","IP"))
# Sorted by increasing drainage area
GRABSmelt$Site_f3<-factor(GRABSmelt$Site,levels=c("IP","LH","BRDS","OMPD","PD","SMD","ID","LMPD","MCLN"))


## Test for significant difference in TDN between >65 km2 and <65km2
t <- t.test(Removal[Removal$CArea.km2.OUT >= 65,]$TDN.mgL.UP,Removal[Removal$CArea.km2.OUT <= 65,]$TDN.mgL.UP)


## Mean RTDN at all sites
Rdf <- Removal %>% filter(RFlag2==0) %>% group_by(., Site) %>% summarize(RTDN = -mean(RTDN, na.rm=TRUE) * 100) %>% 
  mutate(Removal %>% filter(RFlag2==0) %>% group_by(., Site) %>% summarize(RDIN = -mean(RDIN, na.rm=TRUE) * 100)) %>%
  mutate(Removal %>% filter(RFlag2==0) %>% group_by(., Site) %>% summarize(RDON = -mean(RDON, na.rm=TRUE) * 100)) %>%
  as.data.frame()


## Figure 2 - DIN/DON proportions of TDN
quartz(file=paste0("Plots/Fig_2_",Sys.Date(),".pdf"),type='pdf')
cairo_pdf(file=paste0("Plots/Fig_2_",Sys.Date(),".pdf"),family="mono")
#png(filename = paste0("Plots/Fig_2_",Sys.Date(),".png"),type="windows")
pdf(file=paste0("Plots/Fig_2_",Sys.Date(),".pdf"),family="Calibri")

Fig2<-ggplot()+
  #geom_col(data=GRABSmelt[! GRABSmelt$Site=="IP",],aes(x=Site_f2,y=value,fill=Species))+
  geom_col(data=GRABSmelt[GRABSmelt$variable=="DIN.mgL.UP" | GRABSmelt$variable=="DIN.mgL.OUT" | GRABSmelt$variable=="DON.mgL.UP" | GRABSmelt$variable=="DON.mgL.OUT",],aes(x=Site_f3,y=value,fill=Species))+
  #geom_hline(yintercept=0.59932520)+
  geom_vline(xintercept=4.5,size=1.2)+
  annotate("segment",x=1, xend=2, y=-0.11, yend=-0.11)+
  #annotate("text",x=2,y=0.9,label=('bold("DA" > "65 km"^2)'),parse=TRUE,size=3)+
  #annotate("text",x=7,y=0.9,label=('bold("DA" < "65 km"^2)'),parse=TRUE,size=3)+
  scale_fill_manual("",values=c("DIN"="orange","DON"="purple"))+
  facet_wrap(~Location_f,nrow=2)+
  #labs(x=expression(bold("\u2014 Decreasing Drainage Area \u2192")),y=expression(bold(Proportion~of~TDN~("%"))))+
  #labs(x=expression(bold("\u2014"~"Decreasing Drainage Area"~"\u279E")),y=expression(bold(Proportion~of~TDN~("%"))))+
  labs(x=expression(bold("<"~"65"~km^"2"~~~~~~~~~~"\u2014"~"Increasing Drainage Area"~"\u2192"~~~~~~~~~~~">"~"65"~km^"2")),y=expression(bold(Mean~TDN~(mg~L^"-1"))))+
  #labs(x=expression(bold("<"~"65"~km^"2"~~~~~~~~~~"\u2014"~"Increasing Drainage Area" %->% ~~~~~~~~~~">"~"65"~km^"2")),y=expression(bold(Mean~TDN~(mg~L^"-1"))))+
  #labs(x=expression("Farts" %->% "Turds"),y=expression("Balls"))+
  #labs(x=expression(bold(Site)),y=expression(bold(Proportion~of~TDN~("%"))))+
  theme_bw()+
  theme(panel.border = element_rect(color="black",fill=NA,size=1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12,color="black"),#axis.text.x=element_text(family="Arial"),
        axis.title=element_text(face="bold",size=13), axis.title.y=element_text(vjust=2),
        axis.title.x=element_text(vjust=0),legend.title = element_blank(),legend.key=element_blank(),
        plot.title=element_text(size=13),legend.text=element_text(size=9),legend.position="bottom",
        strip.background = element_blank(),strip.text.x = element_text(face="bold",size=12),
        plot.margin=unit(c(5.5,5.5,0,5.5),"pt"),legend.box.margin=margin(t=-0.3,unit='cm'))
Fig2
dev.off()

ggsave(file=paste0("Plots/Fig_2_",Sys.Date(),".pdf"),plot=Fig2,device="pdf")
ggsave(file=paste0("~/Desktop/Dissertation/Manuscripts/CH1_Small_Reservoirs/Figures/Fig2_",
                   Sys.Date(),".pdf"),plot=Fig2)


## Chloride mass balance
Fig3<-ggplot()+
  geom_point(data=Removal[Removal$RFlag2==0 & !is.na(Removal$RFlag2),],aes(x=HL.myr,y=-RCl*100,color="Sample Included",fill="Sample Included"),size=3)+ #,fill="black",color="black"
  geom_point(data=Removal[Removal$RFlag2==1 & !is.na(Removal$RFlag2),],aes(x=HL.myr,y=-RCl*100,color="Sample Omitted",fill="Sample Omitted"),size=3)+ #,fill="red",color="red"
  #geom_hline(yintercept=c(-30,-20,-10,10,20,30))+
  geom_hline(yintercept=c(-20,20))+
  #annotate("text",x=10000,y=-13,label="-10%",size=2.5)+
  annotate("text",x=10000,y=-23,label="-20%",size=2.5)+
  #annotate("text",x=10000,y=-33,label="-30%",size=2.5)+
  #annotate("text",x=10000,y=13,label="+10%",size=2.5)+
  annotate("text",x=10000,y=23,label="+20%",size=2.5)+
  #annotate("text",x=10000,y=33,label="+30%",size=2.5)+
  #geom_hline(yintercept=0)+
  scale_x_continuous(limits=c(0.06,10000),breaks=c(0,1,10,100,1000,10000),trans="log10")+ #limits=c(0.5,10000),
  scale_y_continuous(limits=c(-100,100))+
  scale_shape_manual("# of Inputs",values=c("1"=21,"2"=22,"3"=23),na.translate=FALSE)+
  scale_color_manual("",values=c("Sample Included"="black","Sample Omitted"="red"))+
  scale_fill_manual("",values=c("Sample Included"="black","Sample Omitted"="red"))+
  labs(x=expression(bold(HL~(m~yr^-1))),y=expression(bold(Delta~Cl~("%"))))+
  theme_bw()+
  theme(panel.border = element_rect(color="black",fill=NA,size=1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12,color="black"),
        axis.title=element_text(face="bold",size=13), axis.title.y=element_text(vjust=0),
        axis.title.x=element_text(vjust=0),legend.title = element_text(),legend.key=element_blank(),
        plot.title=element_text(size=13),legend.text=element_text(size=9),legend.position=c(0.25,0.15),
        strip.background = element_blank(),strip.text.x = element_text(face="bold",size=12),legend.box="horizontal",
        plot.margin=unit(c(5.5,5.5,0,5.5),"pt"),legend.spacing.y=unit(-0.1,"cm"),legend.margin=margin(0,0,0,0,unit="cm"))+ #legend.box.margin=margin(t=-0.3,unit='cm')
  guides(shape=guide_legend(order=1),color=guide_legend(order=2),fill=guide_legend(order=2))
Fig3

ggsave(paste0("Plots/Fig_3_",Sys.Date(),".pdf"),plot=Fig3)



# Figure 3 - DIN % removal vs. hydraulic load

Fig4a<-ggplot()+ 
  geom_point(data=Removal[Removal$RFlag2==0 & Removal$RDIN > -2,],aes(x=HL.myr,y=-RDIN*100,shape=Site,color=Temp.C.OUT),size=1.5,stroke=0.75,show.legend=TRUE)+ #,stroke=1
  #Four-Parameter Logistic
  geom_smooth(data=Removal[Removal$RFlag2==0 & Removal$RDIN > -2,],aes(x=HL.myr,y=-RDIN*100,linetype="Logistic Fit - This Study"),method="nls",formula=y~SSfpl(x, A, B, xmid, scal),
              method.args=list(start=list(A=-84,B=25.69,xmid=1.90,scal=0.44)),se=FALSE,color="black",size=0.75)+
  # Power fit
  # This only works if the x-axis has a lower limit of 0.5 m yr
  #geom_smooth(data=Removal[Removal$RFlag2==0 & Removal$RDIN > -2,], aes(x=HL.myr, y=-RDIN*100, linetype="Power Fit - This Study"),
  #            method = 'nls', formula = 'y~a*x^b', method.args = list(start=list(a=-100, b=-1)), se=FALSE, color="black", size=0.75)+
  # Power Fit from equation
  # This works with the same x-axis as the FPL
  stat_function(data=Removal[Removal$RFlag2 == 0 & Removal$RDIN > -2,], fun = function(HL.myr) -66.43064*HL.myr^-0.23992, #-66.43064*HL.myr^-0.23992
                aes(linetype = "Power Fit - This Study"),size=0.75, color="black")+
  #Logarithmic fit
  #geom_smooth(data=Removal[Removal$RFlag2==0 & Removal$RDIN > -2,],
  #            aes(x=HL.myr,y=-RDIN*100,linetype="Logarithmic Fit - This Study"),method='lm',formula="y~log10(x)",se=FALSE,color="black",size=0.75)+
  # Logarithmic Fit
  stat_function(data=Removal[Removal$RFlag2==0 & Removal$RDIN > -2,], fun = function(HL.myr) -85.962+log10(HL.myr)*31.024,
                aes(linetype = "Logarithmic Fit - This Study"), size=0.75, color="black")+
  # Mean annual RDIN for verification purposes
  #geom_point(data=MAChange,aes(x=HL.myr,y=-RDIN*100,shape=Site),color="red",size=2,stroke=1.5)+
  #Seitzinger et al. 2006 #0.5
  #NOT THIS ONE geom_smooth(data=LitRemoval,aes(x=Hl,y=(-Removal_Seit),linetype="Seitzinger et al. 2006"),method='nls',formula='y~a*x^b',method.args=list(start=list(a=88.453,b=-0.3677)),se=FALSE,fullrange=FALSE,color="red")+
  stat_function(data=data.frame(x=seq(0.5,10000,by=1)),fun= function(x) -88.453*x^-0.3677,aes(linetype="Seitzinger et al. 2006"),size=0.75,color="black",show.legend=FALSE)+
  #David et al. 2006
  #stat_function(data=data.frame(x=seq(0.5,10000,by=1)),fun= function(x) -243*x^-0.5632,aes(linetype="David et al. 2006"),size=0.75)+
  #annotate("text",x=0.07,y=85,label='bold("A")',parse=TRUE,size=4)+
  geom_hline(yintercept=c(0))+
  geom_hline(yintercept=c(-20,20),linetype="dashed")+
  #geom_vline(xintercept=290.953)+ #320
  annotate("text",x=8500,y=24,label="+20%",size=2.2)+
  annotate("text",x=8500,y=-24,label="-20%",size=2.2)+
  scale_x_continuous(limits=c(0.06,10000),breaks=c(0,1,10,100,1000,10000),trans="log10")+ #limits=c(0.5,10000),
  #scale_x_continuous(limits=c(0.5,10000),breaks=c(0,1,10,100,1000,10000),trans="log10")+ #limits=c(0.5,10000),
  scale_y_continuous(limits=c(-100,NA))+
  scale_shape_manual("",values=c("BRDS"=0,"ID"=1,"LH"=3,"LMPD"=4,"MCLN"=5,"OMPD"=6,"PD"=2,"SMD"=8),na.translate=FALSE)+
  # DOY scale_color_gradientn("",colors=c("darkblue","blue","green","yellowgreen","orangered","orangered2","orange","blue","darkblue"),values=c(0,0.25,0.5,0.75,1))+ #,values=c(0,0.2,0.4,0.6,0.8,1)
  scale_color_gradient2(name=expression(bold(Temp~(degree~C))),low="darkblue",mid="green",high="red",midpoint=16)+
  #scale_color_manual("",values=c("Winter"="darkblue","Spring"="green","Summer"="red","Autumn"="orange"),breaks=c("Winter","Spring","Summer","Autumn"))+
  scale_linetype_manual("",values=c("Logistic Fit - This Study"="solid", "Power Fit - This Study" = "dotted", "Logarithmic Fit - This Study" = "dotdash", "Seitzinger et al. 2006"="dashed","David et al. 2006"="dotted"),breaks=c("Logistic Fit - This Study", "Power Fit - This Study" ,"Logarithmic Fit - This Study","Seitzinger et al. 2006"))+ #, "David et al. 2006"
  labs(x=expression(bold(HL~(m~yr^-1))),y=expression(bold(Delta~DIN~("%"))))+
  #facet_wrap(~Site,scales="free_y")+ #,scales="free"
  theme_bw()+
  theme(panel.border = element_rect(color="black",fill=NA,size=1), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=8,color="black"),axis.text.x=element_text(hjust=0.5),
        axis.title=element_text(face="bold",size=9), axis.title.y=element_text(vjust=0),
        axis.title.x=element_text(vjust=-0.6),legend.title = element_text(face="bold",size=8),
        legend.key=element_blank(), legend.key.width=unit(2,"lines"), plot.title=element_text(size=13),
        legend.text=element_text(face="bold",size=6),legend.position=c(0.2,0.45),legend.box="vertical",
        legend.margin=margin(-0.5,0,-0.5,0,unit="cm"))+ #,legend.box="vertical"
  guides(linetype=guide_legend(nrow=4,byrow=TRUE,override.aes=list(shape=NA),keyheight=0.5),color="none",
         shape="none") #,linetype=guide_legend(override.                                                                                                                                           es=list(size=2))
# guides(linetype=FALSE,color=guide_colorbar(title.position="top",title.hjust=0.5),
#        shape=guide_legend(label.theme=element_text(size=8,angle=0),override.aes=list(color="black",stroke=1,size=1))) #,linetype=guide_legend(override.aes=list(size=2))
# guides(linetype=guide_legend(nrow=2,byrow=TRUE,override.aes=list(shape=NA,color=c("black","black","black"))), color=FALSE,
#        shape=FALSE) #,linetype=guide_legend(override.aes=list(size=2))
# guides(linetype=guide_legend(nrow=2,byrow=TRUE,override.aes=list(shape=NA,color=c("black","black","black"))), color=guide_colorbar(title.position="top",title.hjust=0.5),
#        shape=guide_legend(label.theme=element_text(size=8,angle=0),override.aes=list(color="black",stroke=1,size=1))) #,linetype=guide_legend(override.aes=list(size=2))
Fig4a



ggsave(file=paste0("Plots/Fig_4a_",Sys.Date(),".png"),plot=Fig4a) #,width=10,height=6 more recently: ,width=7.29,height=6
ggsave(file=paste0("~/Desktop/Dissertation/Manuscripts/CH1_Small_Reservoirs/Figures/Fig3a_",
                   Sys.Date(),".pdf"),plot=Fig3a)

summary(lm(-RDIN*100~log10(HL.myr),Removal[Removal$RFlag==0,]))

## Four parameter logistic
summary(nls(-RDIN*100~SSfpl(log10(HL.myr), A, B, xmid, scal),data=Removal[Removal$RFlag==0,]))

#A: maximum value that can be obtained
#B: minimum value that can be obtained
#xmid: 



# Figure 3b - DON % removal vs hydraulic load
Fig4b<-ggplot()+
  geom_point(data=Removal[Removal$RFlag2==0 & ! Removal$Site=="IP",],aes(x=HL.myr,y=-RDON*100,shape=Site,color=Temp.C.OUT),size=1.5,stroke=0.75)+
  geom_smooth(data=Removal[Removal$RFlag2==0,],aes(x=HL.myr,y=-RDON*100,linetype="Logarithmic Fit"),method='lm',formula='y~x',se=FALSE,color="black",size=0.75)+
  geom_segment(aes(x=0.06,y=100,xend=0.4,yend=100),size=1)+
  annotate("text",x=3.5,y=100,label='bold("Logarithmic Fit")',parse=TRUE,size=2.2)+
  ## Four aprameter logistic fit
  #geom_smooth(data=Removal[Removal$RFlag==0,],aes(x=HL.myr,y=-RDON*100),method="nls",formula=y~SSfpl(x, A, B, xmid, scal),method.args=list(start=list(A=24.59,B=0.63,xmid=2.12,scal=0.05)),se=FALSE,color="red",size=0.75)+
  #geom_smooth(data=Removal[Removal$HL.myr >-1 & Removal$RFlag==0 & ! Removal$Site=="IP",],aes(x=HL.myr,y=RDON*100,linetype="Power Fit"),method='nls',formula='y~a*x^b',method.args=list(start=list(a=-34.7611,b=-0.2277)),se=FALSE,color="black")+
  #stat_function(data=data.frame(x=seq(1,10000,by=1)),fun= function(x) -0.34581*x^-0.22204,aes(linetype="Power Fit"),size=0.75)+
  #annotate("text",x=0.07,y=100,label='bold("B")',parse=TRUE,size=4)+
  geom_hline(yintercept=c(0))+
  geom_hline(yintercept=c(-20,20),linetype="dashed")+
  #geom_vline(xintercept = 3165.018)+
  annotate("text",x=8500,y=24,label="+20%",size=2.2)+
  annotate("text",x=8500,y=-24,label="-20%",size=2.2)+
  #scale_y_continuous(limits=c(-100,100))+
  scale_x_continuous(limits=c(0.06,10000),breaks=c(0,1,10,100,1000,10000),trans="log10")+ #,trans="log10" limits=c(0.5,10000),
  scale_shape_manual("",values=c("BRDS"=0,"ID"=1,"LH"=3,"LMPD"=4,"MCLN"=5,"OMPD"=6,"PD"=2,"SMD"=8))+
  scale_color_gradient2(name=expression(bold(Temp~(degree~C))),low="darkblue",mid="green",high="red",midpoint=16)+
  scale_linetype_manual("",values=c("Logarithmic Fit"="solid","Seitzinger et al. 2006"="dotted","David et al. 2006"="dotdash"),breaks=c("Logarithmic Fit","Seitzinger et al. 2006", "David et al. 2006"))+
  labs(x=expression(bold(HL~(m~yr^-1))),y=expression(bold(Delta~DON~("%"))))+
  #facet_wrap(~Site,scales="free_y")+ #,scales="free"
  theme_bw()+
  theme(panel.border = element_rect(color="black",fill=NA,size=1), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=8,color="black"),axis.text.x=element_text(hjust=0.5),
        axis.title=element_text(face="bold",size=9), axis.title.y=element_text(vjust=0),
        axis.title.x=element_text(vjust=-0.6),legend.title = element_text(face="bold",size=8),legend.key=element_blank(),
        plot.title=element_text(size=13),legend.text=element_text(face="bold",size=10),legend.position="bottom")+
  guides(linetype=FALSE, color=guide_colorbar(title.position="top",title.hjust=0.5),
         shape=guide_legend(label.theme=element_text(size=8,angle=0),override.aes=list(color="black",stroke=1.5))) #,linetype=guide_legend(override.aes=list(size=2))
# guides(linetype=guide_legend(nrow=2,override.aes=list(shape=NA,color=c("black"))), color=guide_colorbar(title.position="top",title.hjust=0.5),
#        shape=guide_legend(label.theme=element_text(size=8,angle=0),override.aes=list(color="black",stroke=1.5))) #,linetype=guide_legend(override.aes=list(size=2))
Fig4b

ggsave(file=paste0("Plots/Fig_3b_4PL_",Sys.Date(),".pdf"),plot=Fig3b,width=7.29,height=6) #,width=10,height=6

summary(lm(RDON~log10(HL.myr),Removal[Removal$RFlag==0 & Removal$RDON > -1,]))
summary(lm(RDON~HL.myr,Removal[Removal$RFlag==0 & Removal$RDON > -1,]))


# Figure 3c - TDN % removal vs hydraulic load

Fig4c<-ggplot()+
  geom_point(data=Removal[Removal$RFlag2==0 & !is.na(Removal$RTDN),],aes(x=HL.myr,y=-RTDN*100,shape=Site,color=Temp.C.OUT),size=1.5,stroke=0.75)+
  geom_smooth(data=Removal[Removal$RFlag2==0,],aes(x=HL.myr,y=-RTDN*100,linetype="Logarithmic Fit"),method='lm',formula='y~x',se=FALSE,color="black",size=0.75)+
  geom_segment(aes(x=0.06,y=50,xend=0.4,yend=50),size=1)+
  annotate("text",x=3.5,y=50,label='bold("Logarithmic Fit")',parse=TRUE,size=2.2)+
  #geom_smooth(data=Removal[Removal$RFlag==0,],aes(x=HL.myr,y=-RTDN*100),method='lm',formula='y~x',se=FALSE,color="black",size=0.75)+
  #geom_smooth(data=Removal[Removal$RFlag==0,],aes(x=HL.myr,y=-RTDN*100),method="nls",formula=y~a*x^b,method.args=list(start=list(a=-19,b=-0.4)),se=FALSE)+
  #geom_smooth(aes(x=HL.myr,y=RTDN,linetype="Logistic fit"),method="nls",formula=y~Asym/(1+exp((xmid-x)/scal)),method.args=list(start=list(Asym=0.55853,xmid=0.08192,scal=-0.29219)),se=FALSE,color="black")+
  #geom_smooth(data=LitRemoval,aes(x=Hl,y=(Removal_Seit),linetype="Seitzinger et al. 2006"),method='nls',formula='y~a*x^b',method.args=list(start=list(a=1,b=-0.368)),se=FALSE,fullrange=FALSE,color="black")+
  #annotate("text",x=0.07,y=85,label='bold("C")',parse=TRUE,size=4)+ #85
  geom_hline(yintercept=0)+
  geom_hline(yintercept=c(-20,20),linetype="dashed")+
  #geom_vline(xintercept=305.0057)+
  annotate("text",x=8500,y=24,label="+20%",size=2.2)+
  annotate("text",x=8500,y=-24,label="-20%",size=2.2)+
  #scale_y_continuous(limits=c(-100,100))+
  scale_x_continuous(limits=c(0.06,10000),breaks=c(0,1,10,100,1000,10000),trans="log10")+ #,trans="log10" ,labels=scales::format_format(big.mark="",decimal_mark="",scientific=FALSE)
  scale_shape_manual("",values=c("BRDS"=0,"ID"=1,"LH"=3,"LMPD"=4,"MCLN"=5,"OMPD"=6,"PD"=2,"SMD"=8))+
  scale_color_gradient2(name=expression(bold(Temp~(degree~C))),low="darkblue",mid="green",high="red",midpoint=16)+
  scale_linetype_manual("",values=c("Logarithmic Fit"="solid","Seitzinger et al. 2006"="dotted"))+
  labs(x=expression(bold(HL~(m~yr^-1))),y=expression(bold(Delta~TDN~("%"))))+
  #facet_wrap(~Site,scales="free_y")+ #,scales="free"
  theme_bw()+
  theme(panel.border = element_rect(color="black",fill=NA,size=1), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=8,color="black"),axis.text.x=element_text(hjust=0.5),
        axis.title=element_text(face="bold",size=9), axis.title.y=element_text(vjust=0),
        axis.title.x=element_text(vjust=-0.6),legend.title = element_text(face="bold",size=8),legend.key=element_blank(),
        plot.title=element_text(size=13),legend.text=element_text(face="bold",size=8),legend.position="bottom",
        legend.box="vertical")+
  guides(linetype=FALSE, color=guide_colorbar(title.position="top",title.hjust=0.5, order=1),
         shape=guide_legend(label.theme=element_text(size=8,angle=0),override.aes=list(color="black",stroke=1.5), order=2)) #,linetype=guide_legend(override.aes=list(size=2))
# guides(linetype=guide_legend(nrow=2,override.aes=list(shape=NA,color=c("black"))), color=guide_colorbar(title.position="top",title.hjust=0.5),
#        shape=guide_legend(label.theme=element_text(size=8,angle=0),override.aes=list(color="black",stroke=1.5))) #,linetype=guide_legend(override.aes=list(size=2))
Fig4c

ggsave(file=paste0("Plots/Fig_3c_",Sys.Date(),".pdf"),plot=Fig3c,width=7.29,height=6) #,width=10,height=6

summary(lm(-RTDN*100~log10(HL.myr),Removal[Removal$RFlag==0,]))
summary(lm(-RTDN*100~HL.myr,Removal[Removal$RFlag==0,]))
summary(nls(-RTDN*100~a*HL.myr^b,Removal[Removal$RFlag==0,],start=list(a=100,b=-1)))


## Figure 4 combination plot
# Fig4comb<-grid.arrange(ggplotGrob(Fig4a+theme(axis.title.x=element_blank())), #+theme(axis.title.x=element_blank())
#                        ggplotGrob(Fig4b+theme(legend.position="none",axis.title.x=element_text(hjust=-3))), #,axis.title.x=element_blank()
#                        ggplotGrob(Fig4c+theme(legend.position="none")),
#                        get_legend(Fig4c),
#                        nrow=2,ncol=2) #,heights=c(2,1.7),widths=c(2,1)
# ggsave(paste0("Plots/Fig_4_V2_",Sys.Date(),".pdf"),plot=Fig4comb,width=6,height=6) #,width=6,height=6

# Fig4comb3<-grid:::grid.draw(rbind(ggplotGrob(Fig4a+theme(axis.title.x=element_blank())),ggplotGrob(Fig4b+theme(legend.position="none",axis.title.x=element_blank())),ggplotGrob(Fig4c),size="last"))
# 
# Fig4comb4<-grid:::grid.draw(rbind(cbind(ggplotGrob(Fig4a+theme(axis.title.x=element_blank(),legend.position=c(0.45,0.935))),
#                                         ggplotGrob(Fig4b+theme(legend.position="none",axis.title.x=element_text(hjust=0.4)))),
#                                   cbind(ggplotGrob(Fig4c+theme(legend.position="none")),
#                                         get_legend(Fig4c))))
# grid::grid.draw(cbind(ggplotGrob(Fig4a+theme(axis.title.x=element_blank(),legend.position=c(0.45,0.935))),
#                       ggplotGrob(Fig4b+theme(legend.position="none",axis.title.x=element_text(hjust=0.4)))))
# 
# grid::grid.draw(cbind(ggplotGrob(Fig4c+theme(legend.position="none")),
#                       cowplot::get_legend(Fig4c)))
# 
# ggsave(file=paste0("~/Desktop/Dissertation/Manuscripts/CH1_Small_Reservoirs/Figures/Fig3_",
#                    Sys.Date(),".pdf"),plot=Fig3comb)

Fig4comb2<-cowplot::plot_grid(cowplot::plot_grid(
  Fig4a+theme(axis.title.x=element_blank(),legend.position=c(0.35,0.84)),
  Fig4b+theme(legend.position="none",axis.title.x=element_text(hjust=-0.4)), #,axis.title.x=element_blank() 
  ncol=2,
  align='h',
  labels=c("A","B"),
  hjust=-5,
  vjust=2),
  cowplot::plot_grid(
    Fig4c+theme(legend.position="none"), #,plot.margin=margin("r"=-1,unit="cm") axis.text.y=element_text(hjust=-0.5)
    get_legend(Fig4c),
    #nrow=2,
    ncol=2,
    labels="C",
    hjust=-5,
    vjust=2),
  #rel_widths=1),
  #align='v'),
  nrow=2,
  ncol=1
)
print(Fig4comb2)

ggsave(paste0("Plots/Fig_4_V2_",Sys.Date(),".pdf"),plot=Fig4comb2,width=6,height=6) #,width=6,height=6

print(ggarrange(Fig4a,Fig4b,Fig4c,common.legend = TRUE))
grid.arrange(Fig3a,Fig3b,Fig3c)


## Figure 4 model fits
## RDIN fit comparison (copied from above)
RDINfitfpl<-nls(-RDIN*100~SSfpl(log10(HL.myr), A, B, xmid, scal),data=Removal[Removal$RFlag2==0 & Removal$RDIN > -2,],na.action = "na.exclude") #,na.action = "na.exclude"
RDINfitlog<-lm(-RDIN*100~log10(HL.myr),Removal[Removal$RFlag2==0 & Removal$RDIN > -2,],na.action = "na.exclude")
RDINfitpow<-nls(RDIN~a*HL.myr^b,Removal[Removal$RFlag2==0 & Removal$RDIN > -2,],start=list(a=1,b=-1),na.action = "na.exclude") # & Removal$HL.myr > 0.5

#{Used this to try to use with as.lm.nls function
## This regression omits the negative and '*100' on RDIN
## Can change that by adding Removal$RDIN2<-(-Removal$RDIN*100)
Removal$RDIN2 <- (-Removal$RDIN * 100)
RDINfitfpl2<-nls(RDIN2~SSfpl(log10(HL.myr), A, B, xmid, scal),data=Removal[Removal$RFlag2==0 & Removal$RDIN > -2,],na.action = "na.omit")
RDINfitpow2<-nls(RDIN2~a*log10(HL.myr)^b,Removal[Removal$RFlag2==0 & Removal$RDIN > -2 & Removal$HL.myr > 0.5,],start=list(a=1,b=-1),na.action = "na.omit")
## Using predictNLS function (in 'R' directory)
#summary(predictNLS(RDINfitfpl2))
##Using as.lm.nls function (in 'R' directory})
summary(as.lm.nls(RDINfitfpl2))
summary(as.lm.nls(RDINfitpow2))
##Values are different but are an approximation due to linearization of nls model:
## adj, R squared = 0.6743; p value = < 2.2e-16; F-stat = 53.28 on 4 & 97 DF

ggplot()+
  geom_point(data=Removal[Removal$RFlag2==0 & Removal$RDIN > -2,],aes(x=HL.myr, y=-RDIN*100))+
  geom_smooth(data=Removal[Removal$RFlag2==0 & Removal$RDIN > -2,],aes(x=HL.myr,y=-RDIN*100),method="nls",formula=y~SSfpl(x, A, B, xmid, scal),
              method.args=list(start=list(A=-84,B=25.69,xmid=1.90,scal=0.44)),se=FALSE,color="black",size=0.75)+
  geom_smooth(data=Removal[Removal$RFlag2==0 & Removal$RDIN > -2 & Removal$HL.myr > 0.5,],aes(x=HL.myr,y=-RDIN*100),method='nls',
              formula='y~a*x^b',method.args=list(start=list(a=-66.43064,b=-0.23992)),se=FALSE,color="black")+
  scale_x_continuous(limits=c(0.06,10000),breaks=c(0,1,10,100,1000,10000),trans="log10")


# ## Create predicted RDIN from 4pl
# Removal$RDINpred4pl<-predict(RDINfitfpl,newdata=data.frame(HL.myr=Removal$HL.myr))
#
# ## RDIN-obs vs RDINpred to "get" significance
# summary(lm(RDIN~RDINpred4pl,Removal[Removal$RFlag==0 & Removal$RDIN > -2,]))

## RDON fit
RDONfitlog<-lm(-RDON*100~log10(HL.myr),Removal[Removal$RFlag2==0,],na.action=na.exclude)

## RTDN fit
RTDNfitlog<-lm(-RTDN*100~log10(HL.myr),Removal[Removal$RFlag2==0,],na.action=na.exclude)

## Calculate RMSE

RMSE = function(obs, fit){
  sqrt(mean((obs - fit)^2, na.rm = TRUE))
}

RMSE(Removal[Removal$RFlag2==0 & Removal$RDIN > -2,]$RDIN, predict(Seit,newdata=Removal[Removal$RFlag2==0 & Removal$RDIN > -2,]$HL.myr))


RMSE(Removal[Removal$RFlag2==0 & Removal$RDIN > -2,]$RDIN, fitted(RDINfitfpl))

RMSE(Removal[Removal$RFlag2==0 & Removal$RDIN > -2 & Removal$HL.myr > 0.5,]$RDIN,fitted(RDINfitpow))

RMSE(Removal[Removal$RFlag2==0 & Removal$RDIN > -2,]$RDIN, fitted(RDINfitlog))


sqrt(mean((Removal[Removal$RFlag2==0 & Removal$RDIN > -2,]$RDIN - predict(RDINfitfpl,newdata=Removal[Removal$RFlag2==0 & Removal$RDIN > -2,]))^2, na.rm=TRUE))


tmp<-data.frame(v1 = predict(RDINfitfpl, newdata = Removal[Removal$RFlag2==0 & Removal$RDIN > -2,]))
tmp <- tmp[!is.na(tmp$v1)]
tmp$v2 <- fitted(RDINfitfpl)

tmp$v2 <- unlist(RDINfitfpl$m$fitted())

## Comparew model fits
## Trying again because using negative and * 100 and log10(HL.myr) was not producing RDIN similar to RDIN

RDINfitfpl<-nls(-RDIN*100~SSfpl(log10(HL.myr), A, B, xmid, scal),data=Removal[Removal$RFlag2==0 & Removal$RDIN > -2,],na.action = "na.exclude") #,na.action = "na.exclude"
RDINfitlog<-lm(-RDIN*100~log10(HL.myr),Removal[Removal$RFlag2==0 & Removal$RDIN > -2,],na.action = "na.exclude")
RDINfitpow<-nls(-RDIN*100~a*HL.myr^b,Removal[Removal$RFlag2==0 & Removal$RDIN > -2,],start=list(a=100,b=-1),na.action = "na.exclude") # & Removal$HL.myr > 0.5

RDINfitfpl2<-nls(RDIN2~SSfpl(log10(HL.myr), A, B, xmid, scal),data=Removal[Removal$RFlag2==0 & Removal$RDIN > -2,],na.action = "na.omit")
RDINfitlog2<-lm(RDIN~log10(HL.myr),Removal[Removal$RFlag2==0 & Removal$RDIN > -2,],na.action = "na.omit")
RDINfitpow2<-nls(RDIN~a*HL.myr^b,Removal[Removal$RFlag2==0 & Removal$RDIN > -2,],start=list(a=1,b=-1),na.action = "na.omit") # & Removal$HL.myr > 0.5


Removal$RDINfpl <- predict(RDINfitfpl, newdata=data.frame(HL.myr = Removal$HL.myr))
Removal$RDINpow <- predict(RDINfitpow, newdata=data.frame(HL.myr = Removal$HL.myr))
Removal$RDINlog <- predict(RDINfitlog, newdata=data.frame(HL.myr = Removal$HL.myr))

Removal$RDINSeit2 <- 0.88453*Removal$HL.myr^-0.3677
Removal$RDINSeit2 <- Removal$RDINSeit2+rnorm(length(Removal$RDINSeit2),sd=0.01)


Removal$RDINSeit <- predict(Seit, newdata=data.frame(Hl = Removal$HL.myr))
Removal$RDINSeit <- Removal$RDINSeit+rnorm(length(Removal$RDINSeit),sd=0.01)
Removal$RDINDavid <- predict(David, newdata=data.frame(Hl = Removal$HL.myr))
Removal$RDINDavid <- Removal$RDINDavid+rnorm(length(Removal$RDINDavid),sd=0.01)


RMSE(Removal[Removal$RFlag2==0 & Removal$RDIN > -2,]$RDIN, fitted(RDINfitfpl2))
RMSE(Removal[Removal$RFlag2==0 & Removal$RDIN > -2,]$RDIN, fitted(RDINfitpow))
RMSE(Removal[Removal$RFlag2==0 & Removal$RDIN > -2,]$RDIN, fitted(RDINfitlog))
RMSE(-Removal[Removal$RFlag2==0 & Removal$RDIN > -2,]$RDIN*100, 
     predict(Seit, newdata=data.frame(Hl = Removal$HL.myr))) #Seitzinger model fit to my HL
RMSE(Removal[Removal$RFlag2==0 & Removal$RDIN > -2,]$RDIN, 
     predict(David, newdata=data.frame(Hl = Removal$HL.myr))) #David model fit to my HL

Seit4 <- nls(RDINSeit2~a*HL.myr^b, Removal[Removal$RFlag2 == 0 & Removal$RDIN > -2,], start=list(a=0.88, b=-0.3),na.action="na.exclude")
David2 <- nls(RDINDavid~a*HL.myr^b, Removal[Removal$RFlag2 == 0 & Removal$RDIN > -2,], start=list(a=2.43, b=-0.568))


nls(RDIN~B + ((A-B)/(1 + exp((HL.myr-xmid)/scal))), data=Removal[Removal$RFlag2==0 & Removal$RDIN > -2,],
    start=list(A = 6.3, B = 0.75, xmid = 62.9, scal = 2),nls.control(maxiter=15000, minFactor = 6e-10))

AIC(RDINfitfpl,RDINfitpow,RDINfitlog)
AIC(RDINfitfpl2,RDINfitpow2,RDINfitlog2)

BIC(RDINfitfpl,RDINfitpow,RDINfitlog)
BIC(RDINfitfpl2,RDINfitpow2,RDINfitlog2)



