





Q <- do.call(dplyr::bind_rows, lapply(list.files(path=("C:/Users/ctw1/Box/Data/Analysis/LTER_Update_2021/Data/CC/Discharge/"),pattern="*.csv",full=TRUE),
                                 read.csv, skip=0, header=TRUE, stringsAsFactors=FALSE))


Q$DateTime <- as.POSIXct(strptime(paste(Q$Date, Q$Time), format = "%m/%d/%Y %H:%M:%S", tz = "EST"))
Q$Date <- as.Date(Q$Date,format = "%m/%d/%Y")


Qdaily <- Q %>% group_by(Date) %>% select(-c(Time, Depth_Device, DateTime)) %>% summarize_all(mean ,na.rm=TRUE)

write.csv(Qdaily, file= "C:/Users/ctw1/Box/Data/Analysis/LTER_Update_2021/Data/CC_Q_Daily_2022-05-09.csv",row.names=FALSE)


OYS_Q_Daily<-read.table(file=paste0("https://waterdata.usgs.gov/nwis/dv?cb_00060=on&format=rdb&site_no=01073000&referred_module=sw&period=&begin_date=",Startdate,"&end_date=",Enddate),skip=30,sep="\t",stringsAsFactors=FALSE)

OYS_Q_Daily$Date <- as.Date(strptime(OYS_Q_Daily[,3],format="%Y-%m-%d"))
OYS_Q_Daily$Year <- as.numeric(format(OYS_Q_Daily$Date,"%Y"))
OYS_Q_Daily$Discharge.m3s <- OYS_Q_Daily[,4]/35.3147
OYS_Q_Daily$Date2 <- lubridate::floor_date(as.Date(strptime(OYS_Q_Daily[,3],format="%Y-%m-%d")),"month")
OYS_Q_Daily2 <- OYS_Q_Daily[,c(6,8,7)]

write.csv(OYS_Q_Daily2, file= "C:/Users/ctw1/Box/Data/Analysis/Dams/Data/Oyster_Q_Daily_2022-05-09.csv", row.names=FALSE)



ProcessUSGSdaily(OYS_Q_Daily)
ProcessUSGSdaily(ID_Q_Daily)



ProcessUSGSdaily <- function(x) {
  x <- x %>% mutate(
    "Date" = as.Date(strptime(x[,3],format="%Y-%m-%d")),
    "Year" = as.numeric(format(as.Date(strptime(x[,3],format="%Y-%m-%d")),"%Y")),
    "Discharge.m3s" =  x[,4]/35.3147,
    "Discharge.Ls" = x[,4]/35.3147/1000,
    "Date2" = lubridate::floor_date(as.Date(strptime(x[,3],format="%Y-%m-%d")),"month"))
  x <- x[,c(6,8,9,7,10)]
}
