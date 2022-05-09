#Analysis and figures for 'Small Dams' manuscript
#Created 2019-11-19

#Sections:
#1 - Read in grab sample data




rm(list=ls())
#Set Working Directory
setwd("C:/Users/ctw1/Box/Data/Analysis/Dams/")
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
GRABS$DateTime<-as.POSIXct(strptime(paste(GRABS$Date,GRABS$Time),format="%m/%d/%y %H:%M",tz="America/New_York"))
attributes(GRABS$DateTime)$tzone<-"EST"
GRABS$Date<-as.POSIXct(strptime(GRABS$Date,format="%m/%d/%y"))
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










