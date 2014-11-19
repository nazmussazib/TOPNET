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


Model='CCSM4'
setwd(paste("F:\\USU_Research_work_update_March_30_2014\\TOPNET PROJECT\\MODEL COMPARISON\\CLIMATE_CHANGE_ALLDATA\\",Model,sep=""))
pr.nc = nc_open("BCCAv2_0.125deg_pr_day_CCSM4_rcp45_r2i1p1_20560101-20651231.nc")
latall=ncvar_get(pr.nc ,'latitude')
lonall=ncvar_get(pr.nc ,'longitude')


###change the directory for each watershed #########
watsed=c("A1","A2","B1","B21","B22","C1","C21","C22")
w=2
setwd(paste("E:\\USU_Research_work\\TOPNET PROJECT\\MODEL COMPARISON\\",watsed[w],"_watershed_final\\",watsed[w],"_CC_data",sep=""))##for A1:


##change rain gauge shape file name

lat_lon_rg=read.dbf(paste("rg_",watsed[w],".dbf",sep="")) 
lat_lon_rg=do.call(cbind,lat_lon_rg)
lat_lon_rg1=albersgeod(lat_lon_rg[,1], lat_lon_rg[,2], sph="GRS80", clon=-96, clat=23, sp1=29.5, sp2=45.5)

lon_all=as.numeric(lonall)
lon_ws=signif(lat_lon_rg1[,1]+360,7)
lon_ID=match(lon_ws,lon_all,)

lat_all=as.numeric(latall)
lat_ws=signif(lat_lon_rg1[,2],7)
lat_ID=match(lat_ws,lat_all)
lat_lon_ID=data.frame(lon_ID=lon_ID,lat_ID=lat_ID)
len=length(lat_ID)

setwd(paste("F:\\USU_Research_work_update_March_30_2014\\TOPNET PROJECT\\MODEL COMPARISON\\CLIMATE_CHANGE_ALLDATA\\",Model,sep=""))
ua_file=list.files(path =paste("F:\\USU_Research_work_update_March_30_2014\\TOPNET PROJECT\\MODEL COMPARISON\\CLIMATE_CHANGE_ALLDATA\\",Model,sep=""),pattern =paste("*ua_day_",Model,"_rcp",sep=""))
va_file=list.files(path =paste("F:\\USU_Research_work_update_March_30_2014\\TOPNET PROJECT\\MODEL COMPARISON\\CLIMATE_CHANGE_ALLDATA\\",Model,sep=""),pattern =paste("*va_day_",Model,"_rcp",sep=""))
##get index of the nearest lat,lon of the grid ###
ua.nc = nc_open(ua_file[1])
ua_latall=ncvar_get(ua.nc,'lat')
ua_lonall=ncvar_get(ua.nc,'lon')
len=length(lon_ws)
ua_lon_ID=matrix(NA,nrow=len,ncol=1)
ua_lat_ID=matrix(NA,nrow=len,ncol=1)
for ( i in 1:len){
  ua_lon_ID[i,1]=which.min(abs(ua_lonall-lon_ws[i])) 
  ua_lat_ID[i,1]=which.min(abs(ua_latall-lat_ws[i])) 
}
ua_lat_lon_ID=cbind(ua_lon_ID,ua_lat_ID)
ua_ID=unique(ua_lat_lon_ID)



##get index of the nearest lat,lon of the
uawind <- function(file ){
  
  ua.nc = nc_open(file)
  ua_time=as.Date(ncvar_get(ua.nc,'time'),origin="2005-01-01")
  #strDates= as.character(hur_time)
  #ua=matrix(NA,nrow=length(ua_time),ncol=length(ua_ID[,1]))
  ua=matrix(NA,nrow=length(ua_time),ncol=1)
  
  for( i in 1:1){
    start=c(ua_ID[i,1],ua_ID[i,2],3,1)
    count=c(1,1,1,length(ua_time))
    ua[,i]=ncvar_get(ua.nc,'ua',start,count)
  }
  uadata=data.frame(data=cbind(ua, ua_time-0.5))
  return(uadata)
}




uad=data.frame(matrix(NA,nrow=1,ncol=1+1))

for ( m in 1:length(ua_file)){
  ud=uawind(ua_file[m])
  ne=colnames(ud)
  colnames(uad)[]=ne
  uad=rbind(uad,ud)
}


##get index of the nearest lat,lon of the grid ###
va.nc = nc_open(va_file[1])
va_latall=ncvar_get(va.nc,'lat')
va_lonall=ncvar_get(va.nc,'lon')
len=length(lon_ws)
va_lon_ID=matrix(NA,nrow=len,ncol=1)
va_lat_ID=matrix(NA,nrow=len,ncol=1)
for ( i in 1:len){
  va_lon_ID[i,1]=which.min(abs(va_lonall-lon_ws[i])) 
  va_lat_ID[i,1]=which.min(abs(va_latall-lat_ws[i])) 
}
va_lat_lon_ID=cbind(va_lon_ID,va_lat_ID)
va_ID=unique(va_lat_lon_ID)

##get index of the nearest lat,lon of the
vawind <- function(file ){
  
  va.nc = nc_open(file)
  va_time=as.Date(ncvar_get(va.nc,'time'),origin="2005-01-01")
  #strDates= as.character(hur_time)
  va=matrix(NA,nrow=length(va_time),ncol=1)
  for( i in 1:1){
    start=c(va_ID[i,1],va_ID[i,2],3,1)
    count=c(1,1,1,length(va_time))
    va[,i]=ncvar_get(va.nc,'va',start,count)
  }
  vadata=data.frame(data=cbind(va, va_time-0.5))
  return(vadata)
}

###it seems that file based on th rcp###
vad=data.frame(matrix(NA,nrow=1,ncol=1+1))
for ( m in 1:length(va_file)){
  vd=vawind(va_file[m])
  ne=colnames(vd)
  colnames(vad)[]=ne
  vad=rbind(vad,vd)
}

all_wind=sqrt(((uad[,1]))^2+((vad[,1]))^2)
all_wind=data.frame(wind=all_wind)

wind_data=cbind(all_wind[2:9491,],all_wind[9492:18981,],uad[2:9491,2])
dat_1=seq(as.Date("2056/1/1"), as.Date("2065/12/31"), "day")
dat_2=seq(as.Date("2090/1/1"), as.Date("2099/12/31"), "day")
match_w1=match(dat_1,wind_data[,3])
match_w2=match(dat_2,wind_data[,3])
wind_data=data.frame(wind_data)
wind_rcp_data=rbind(wind_data[match_w1,],wind_data[match_w2,])

##convert at 2m height
##we took wind speed at 85000 pa level that means 1.41km avove the ground

##Reference

##But we need wind speed at 2m avpve ground

##equation taken form http://www.fao.org/docrep/x0490e/x0490e07.htm
wind_rcp_data_mod=wind_rcp_data[,1:2]*4.87/(log(67.8*3.6*1000-5.42)) ##1.2km*(100000-850000)=3.6km=3.6*1000 m


############WRITING WIND DATA###################
strDates=rbind(data.frame(date=dat_1),data.frame(date=dat_2))
gh=gsub("-","", strDates$date, fixed=TRUE)
day=data.frame(time=gh)
hour=rep.int(240000, nrow(day))
rcp=2
rcp_name_w=c("rcp4.5","rcp8.5")
for (j in 1:rcp){
  wind=data.frame(wind_rcp_data_mod[,j])
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










