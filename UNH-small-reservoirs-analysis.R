#Analysis and figures for 'Small Dams' manuscript
#Created 2019-11-19

#Sections:
#1 - Read in grab sample data




rm(list=ls())
#Set Working Directory
#setwd("C:/Users/ctw1/Box/Data/Analysis/Dams/")
Sys.Date<-as.character(Sys.Date())

library(gridExtra)
library(reshape2)
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
LitRemoval$Removal_Seit<-88*LitRemoval$Hl^-0.368
LitRemoval$Removal_Seit<-LitRemoval$Removal_Seit+rnorm(length(LitRemoval$Removal_Seit),sd=0.01)

LitRemoval$Removal_David<-243*LitRemoval$Hl^-0.5632
LitRemoval$Removal_David<-LitRemoval$Removal_David+rnorm(length(LitRemoval$Removal_David),sd=0.01)

LitRemoval<-melt(LitRemoval,id.vars=c("Hl"))
LitRemoval$Eq<-NA
LitRemoval[LitRemoval$variable=="Removal_Seit",]$Eq<-"Seitzinger et al."
LitRemoval[LitRemoval$variable=="Removal_David",]$Eq<-"David et al."


#####
################ Read in discharge data from USGS website ################
#####

##USGS has a package "dataRetrieval" that can download data - look into this

## Can skip the below downloading/.processing steps and import data from file:

Q_Daily <- read.csv("Data/Q_Daily_2022-06-08.csv", stringsAsFactors = FALSE, header = TRUE)

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





















