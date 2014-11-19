##############wind correction##################\
library(ncdf4)
library(proj4)
library(rgdal)
library(raster)
library(shapefiles)
require(spsurvey)
require(geosphere)
require(TTR)
require(hydroTSM)
require(qmap)
require(ggplot2)
require(grid)
require(gridExtra)


Model='CanESM2'
##change the directory for each watershed #######
watsed=c("A1","A2","B1","B21","B22","C1","C21","C22")
w=1
setwd(paste("E:\\USU_Research_work\\TOPNET PROJECT\\MODEL COMPARISON\\",watsed[w],"_watershed_final\\",watsed[w],"_CC_data\\",Model,sep=""))##for A1:


##change rain gauge shape file name

lat_lon_rg=read.dbf(paste("rg_",watsed[w],".dbf",sep="")) 
lat_lon_rg=do.call(cbind,lat_lon_rg)
lat_lon_rg1=albersgeod(lat_lon_rg[,1], lat_lon_rg[,2], sph="GRS80", clon=-96, clat=23, sp1=29.5, sp2=45.5)
lon_ws=signif(lat_lon_rg1[,1]+360,7)
lat_ws=signif(lat_lon_rg1[,2],7)

setwd(paste("F:\\USU_Research_work_update_March_30_2014\\TOPNET PROJECT\\MODEL COMPARISON\\CLIMATE_CHANGE_ALLDATA\\",Model,sep=""))
sfc_file=list.files(path =paste("F:\\USU_Research_work_update_March_30_2014\\TOPNET PROJECT\\MODEL COMPARISON\\CLIMATE_CHANGE_ALLDATA\\",Model,sep=""),pattern =paste("*sfcWind_day_",Model,sep=""))
#get index of the nearest lat,lon of the grid ###
sf.nc = nc_open(sfc_file[1])
sa_latall=ncvar_get(sf.nc,'lat')
sa_lonall=ncvar_get(sf.nc,'lon')
len=length(lon_ws)
sa_lon_ID=matrix(NA,nrow=len,ncol=1)
sa_lat_ID=matrix(NA,nrow=len,ncol=1)
for ( i in 1:len){
  sa_lon_ID[i,1]=which.min(abs(sa_lonall-lon_ws[i])) 
  sa_lat_ID[i,1]=which.min(abs(sa_latall-lat_ws[i])) 
}
sa_lat_lon_ID=cbind(sa_lon_ID,sa_lat_ID)
sa_ID=unique(sa_lat_lon_ID)



##get index of the nearest lat,lon of the

wind <- function(file ){
  
  w.nc = nc_open(file)
  w_time=as.Date(ncvar_get(w.nc,'time'),origin="1850-01-01")
  #strDates= as.character(hur_time)
  #ua=matrix(NA,nrow=length(ua_time),ncol=length(ua_ID[,1]))
  w=matrix(NA,nrow=length(w_time),ncol=1)
  
  for( i in 1:1){
    start=c(sa_ID[1,1],sa_ID[1,2],1)
    count=c(1,1,length(w_time))
    w[,i]=ncvar_get(w.nc,'sfcWind',start,count)
  }
  wdata=data.frame(data=cbind(w, w_time-0.5))
  return(wdata)
}




ws=data.frame(matrix(NA,nrow=1,ncol=1+1))

for ( m in 1:length(sfc_file)){
  wd=wind(sfc_file[m])
  ne=colnames(wd)
  colnames(ws)[]=ne
  ws=rbind(ws,wd)
}



all_wind=data.frame(windsp=ws)
all_wind=data.frame(all_wind[-1,])
wind_data=cbind(all_wind[1:(nrow(all_wind)/2),],all_wind[((nrow(all_wind)/2)+1):nrow(all_wind),])

dat_1=seq(as.Date("2056/1/1"), as.Date("2065/12/31"), "day")
dat_2=seq(as.Date("2090/1/1"), as.Date("2099/12/31"), "day")

match_w1=match(dat_1,wind_data[,2])
match_w2=match(dat_2,wind_data[,2])
wind_data=data.frame(wind_data)
wind_rcp_data=rbind(wind_data[match_w1,],wind_data[match_w2,])
wind_rcp_data=wind_rcp_data[,c(1,3)]
############WRITING WIND DATA###################
strDates=rbind(data.frame(date=dat_1),data.frame(date=dat_2))
gh=gsub("-","", strDates$date, fixed=TRUE)
day=data.frame(time=gh)
hour=rep.int(240000, nrow(day))
rcp=2
rcp_name_w=c("rcp4.5","rcp8.5")
for (j in 1:rcp){
  wind=data.frame(wind_rcp_data[,j])
  wind[] <- lapply(wind, function(.col){ if (is.numeric(.col)) return(sprintf("%8.2f",.col))else return(.col)})
  wind=data.frame(wind,day,hour)
  file=paste("wind_",rcp_name_w[j],".dat",sep="")
  sink(file)
  cat(sprintf("This file provides daily values of wind speed for each wind station"),file=file,append=TRUE)
  cat("\n", file=file, append=TRUE)
  cat(sprintf("Wind speed is provided in m/sec"),file=file,append=TRUE)
  cat("\n", file=file, append=TRUE)
  sites=seq(1,1,1) ## need to automate this one
  cat(sprintf("%s %d ", "ver2",1),file=file,append=TRUE) 
  cat(sprintf( "%d", sites),(sprintf( "%s", "Date Hour","\n")),file=file,append=TRUE) 
  cat("\n", file=file, append=TRUE)
  write.table(wind, file =file,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
  sink()
}








