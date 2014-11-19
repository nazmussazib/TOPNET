###########tmaxtmintdewcorrection#############
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

Mod=c('CCSM4','CanESM2','GFDL-ESM2G','MPI-ESM-MR','MRI-CGCM3','CNRM-CM5','MIROC5')
Model=Mod[1]
setwd(paste("F:\\USU_Research_work_update_March_30_2014\\TOPNET PROJECT\\MODEL COMPARISON\\CLIMATE_CHANGE_ALLDATA\\",Model,sep=""))
pr.nc = nc_open("BCCA_0.125deg_tasmax_day_CCSM4_rcp45_r2i1p1_20060101-20151231.nc")
latall=ncvar_get(pr.nc ,'latitude')
lonall=ncvar_get(pr.nc ,'longitude')


###change the directory for each watershed #########
watsed=c("A1","A2","B1","B21","B22","C1","C21","C22")
w=2
setwd(paste("E:\\USU_Research_work\\TOPNET PROJECT\\MODEL COMPARISON\\",watsed[w],"_watershed_final\\",watsed[w],"_CC_data\\",Model,sep=""))##for A1:


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

#####DAYMET data for each of the pixel #####################


xx="file"
nc=len
for (i in 1:nc){
  xx[i]=paste("filename",i,".csv",sep="")
}
rainfall <- function(file ){
  mydata2 = read.csv(file,skip=6)
  my.data.frame <- mydata2[(mydata2$year%%4==0) & (mydata2$yday ==365), ]
  my.data.frame$yday=my.data.frame$yday+1
  total <- rbind( my.data.frame,mydata2)
  gh= total[with(total, order(year,yday)), ]
  rain=data.frame(rainmm=gh$prcp..mm.day., tmaxcel=gh$tmax..deg.c.,tmincel=gh$tmin..deg.c.,vppa= gh$vp..Pa.)
  
  return(rain)
}

dates_day=seq(as.Date("1980/1/1"), as.Date("2012/12/31"), "day") 
dss=data.frame(date=dates_day)

for ( m in 1:length(xx)){
  hjk=rainfall(xx[m])
  dss=cbind(dss,hjk)
}

dewT_d=matrix(ncol =nc, nrow = nrow(dss))
avgT_d=matrix(ncol = nc, nrow = nrow(dss))
satVp_d=matrix(ncol = nc, nrow = nrow(dss))
actVp_d=matrix(ncol = nc, nrow = nrow(dss))
precipitation_d=matrix(ncol = nc, nrow = nrow(dss))
tmaxtmintdew_d=matrix(ncol = 1, nrow = nrow(dss))


for (j in 1:(nc)){for (i in 1:nrow(dss))
{
  avgT_d[i,j]=0.5*(dss[i,4*j-1]+dss[i,4*j])
  satVp_d[i,j]=6.112*exp(17.67*avgT_d[i,j]/(avgT_d[i,j]+243.5)) # unit is hecto-pascal
  actVp_d[i,j]=satVp_d[i,j]-dss[i,4*j+1]/10000 # convert vapor pressure deficit pa to hpa
  
  dewT_d[i,j]=(log(actVp_d[i,j]/(6.11)))*243.5/(17.67-log(actVp_d[i,j]/(6.11)))
  
}
}

for ( k in 1:(nc)){
  precipitation_d[,k]=dss[,4*k-2]
  tmaxtmintdew_d=cbind(tmaxtmintdew_d,dss[,4*k-1],dss[,4*k],dewT_d[,k])
}


Elevation<- function(file ){
  station_ele = as.matrix(read.csv(file))[3]
  ele = (as.numeric((unlist(strsplit(station_ele, split=' ', fixed=TRUE))[2])))
  return(ele)
}

Ele_station=matrix(ncol =1, nrow = nc)

for ( m in 1:length(xx)){
  
  Ele_station[m,]=Elevation(xx[m])
  
}


station_monthly<- data.frame(monthlyfunction(avgT_d, FUN=mean, na.rm=TRUE,dates=dates_day))
station_annual<- t(as.matrix(annualfunction(avgT_d, FUN=mean, na.rm=TRUE,dates=dates_day)))
sta_ID=seq(1,nc,1)
sta_ele=rep(-999,nc)
station_monthly=cbind(sta_ID,lat_lon_rg1[,2],lat_lon_rg1[,1],sta_ele,station_monthly)
lapse_mat=cbind(Ele_station,station_annual)
dddd=lm(lapse_mat[,2]~lapse_mat[,1])
station_monthly[,4]=lapse_mat[,1]
station_monthly[] <- lapply(station_monthly, function(.col){ if (is.numeric(.col)) return(sprintf("%8.2f",.col))else return(.col)}) 

#writng clipar.dat
sink("clipar.dat") 
cat(sprintf("This file provides meta-data for each temperature station"),file='clipar.dat',append=TRUE)
cat("\n", file="clipar.dat", append=TRUE)
cat(sprintf("Temperature ranges are in Degrees Celcius"),file='clipar.dat',append=TRUE)
cat("\n", file="clipar.dat", append=TRUE)
sites=seq(1,(nc),1)
cat(sprintf("%d %s ", nc,"! Number of temperature stations"),file='clipar.dat',append=TRUE) 
cat("\n", file="clipar.dat", append=TRUE)
##please change standard longitude according to watershed location

std_longitude=-120.00;
cat(sprintf("%3.2f %s ",std_longitude,"!Standard Longitude of time zone used in local time calculations"),file='clipar.dat',append=TRUE) 
cat("\n", file="clipar.dat", append=TRUE)
cat(sprintf("'Station_id  Lat    Lon     Elev_m   Jan    Feb    Mar     Apr       May    Jun    Jul      Aug     Sep    Oct      Nov    Dec"),file='clipar.dat',append=TRUE)
cat("\n", file="clipar.dat", append=TRUE)
write.table(station_monthly, file = "clipar.dat",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
sink()




temper_day=data.frame(tmaxtmintdew_d[,-1])
dates_bc=seq(as.Date("1980/1/1"), as.Date("2005/12/31"), "day") 
date_overlap_day=match(dates_bc,dates_day)
#####processiing temperature  data #############

setwd(paste("F:\\USU_Research_work_update_March_30_2014\\TOPNET PROJECT\\MODEL COMPARISON\\CLIMATE_CHANGE_ALLDATA\\",Model,sep=""))

toNumerics <- function(Date) {
  stopifnot(inherits(Date, c("Date", "POSIXt")))
  day <- as.numeric(strftime(Date, format = "%d"))
  month <- as.numeric(strftime(Date, format = "%m"))
  year <- as.numeric(strftime(Date, format = "%Y"))
  list(year = year, month = month, day = day)
}

dff=t(matrix(unlist(lapply(dates_bc, toNumerics)),nrow=3,ncol=length(dates_bc)))



############Analysis of Rainfall Results compare to base period 2003-2012 ##########

dates_day=seq(as.Date("1980/1/1"), as.Date("2012/12/31"), "day") 
dates_base=seq(as.Date("2003/1/1"), as.Date("2012/12/31"), "day") ##base line period for 2003-2012

date_overlap_base=match(dates_base,dates_day) 



##################TEMPERATURE DATA PROCESSING#########################

dates_obs=seq(as.Date("1980/1/1"), as.Date("2005/12/31"), "day") 
date_overlap_obs=match(dates_obs,dates_day)
tmax_obs=(tmaxtmintdew_d[date_overlap_obs,seq(2,dim(tmaxtmintdew_d)[2],3)])
tmin_obs=(tmaxtmintdew_d[date_overlap_obs,seq(3,dim(tmaxtmintdew_d)[2],3)])
tdew_obs=(tmaxtmintdew_d[date_overlap_obs,seq(4,dim(tmaxtmintdew_d)[2],3)])


tmax <- function(file ){
  
  tmax.nc = nc_open(file)
  tmax_time=as.Date(ncvar_get(tmax.nc,'time'),origin="1950-01-01")
  #strDates= as.character(rain_time)
  tmax=matrix(NA,nrow=length(tmax_time),ncol=len)
  for (i in 1:len){
    start=c(lat_lon_ID[i,1],lat_lon_ID[i,2],1)
    count=c(1,1,length( tmax_time))
    tmax[,i]=ncvar_get(tmax.nc,'tasmax',start,count)
  }
  tmaxdata=data.frame(data=cbind(tmax, tmax_time-0.5))
  
  return(tmaxdata)
}


tmin <- function(file ){
  
  tmin.nc = nc_open(file)
  tmin_time=as.Date(ncvar_get(tmin.nc,'time'),origin="1950-01-01")
  #strDates= as.character(rain_time)
  tmin=matrix(NA,nrow=length(tmin_time),ncol=len)
  for (i in 1:len){
    start=c(lat_lon_ID[i,1],lat_lon_ID[i,2],1)
    count=c(1,1,length( tmin_time))
    tmin[,i]=ncvar_get(tmin.nc,'tasmin',start,count)
  }
  tmindata=data.frame(data=cbind(tmin, tmin_time-0.5))
  
  return(tmindata)
}



#############Reading projected data for##########

hist_tmax_file=list.files(path =paste("F:\\USU_Research_work_update_March_30_2014\\TOPNET PROJECT\\MODEL COMPARISON\\CLIMATE_CHANGE_ALLDATA\\",Model,sep=""),pattern =paste("*_tasmax_day_",Model,"_historical_",sep=""))
hist_tmin_file=list.files(path =paste("F:\\USU_Research_work_update_March_30_2014\\TOPNET PROJECT\\MODEL COMPARISON\\CLIMATE_CHANGE_ALLDATA\\",Model,sep=""),pattern =paste("*_tasmin_day_",Model,"_historical_",sep=""))

proj_tmax_file=list.files(path =paste("F:\\USU_Research_work_update_March_30_2014\\TOPNET PROJECT\\MODEL COMPARISON\\CLIMATE_CHANGE_ALLDATA\\",Model,sep=""),pattern =paste("*_tasmax_day_",Model,"_rcp",sep=""))
proj_tmin_file=list.files(path =paste("F:\\USU_Research_work_update_March_30_2014\\TOPNET PROJECT\\MODEL COMPARISON\\CLIMATE_CHANGE_ALLDATA\\",Model,sep=""),pattern =paste("*_tasmin_day_",Model,"_rcp",sep=""))

his_tmaximum=data.frame(matrix(NA,nrow=1,ncol=len+1))
his_tminimum=data.frame(matrix(NA,nrow=1,ncol=len+1))
proj_tmaximum=data.frame(matrix(NA,nrow=1,ncol=len+1))
proj_tminimum=data.frame(matrix(NA,nrow=1,ncol=len+1))



for ( m in 1:length(hist_tmax_file)){
  tmx=tmax(hist_tmax_file[m])
  ne=colnames(tmx)
  colnames(his_tmaximum)[]=ne
  his_tmaximum=rbind(his_tmaximum,tmx)
  tmn=tmin(hist_tmin_file[m])
  nen=colnames(tmn)
  colnames(his_tminimum)[]=nen
  his_tminimum=rbind(his_tminimum,tmn)
  
}


for ( m in 1:length(proj_tmax_file)){
  tmx=tmax(proj_tmax_file[m])
  ne=colnames(tmx)
  colnames(proj_tmaximum)[]=ne
  proj_tmaximum=rbind(proj_tmaximum,tmx)
  tmn=tmin(proj_tmin_file[m])
  nen=colnames(tmn)
  colnames(proj_tminimum)[]=nen
  proj_tminimum=rbind(proj_tminimum,tmn)
  
}


his_tmaximum1=as.matrix((his_tmaximum[-1,1:len]))
his_tminimum1=as.matrix((his_tminimum[-1,1:len]))

##only for gddl#####
gfdl_bc=seq(as.Date("1976/1/1"), as.Date("2005/12/31"), "day") 
match_gfdl=match(dates_bc,gfdl_bc)
his_tmaximum1=his_tmaximum1[match_gfdl,]
his_tminimum1=his_tminimum1[match_gfdl,]
###########

proj_tmaximum1=as.matrix((proj_tmaximum[-1,1:len]))
proj_tminimum1=as.matrix((proj_tminimum[-1,1:len]))

tmax_rcp_all=cbind(proj_tmaximum1[1:12783,],proj_tmaximum1[12784:25566,],proj_tmaximum[2:12784,len+1])
tmin_rcp_all=cbind(proj_tminimum1[1:12783,],proj_tminimum1[12784:25566,],proj_tminimum[2:12784,len+1])
##for MPI
tmax_rcp_all=cbind(proj_tmaximum1[1:(nrow(proj_tmaximum1)/2),],proj_tmaximum1[(1+(nrow(proj_tmaximum1)/2)):nrow(proj_tmaximum1),],proj_tmaximum[2:(1+(nrow(proj_tmaximum1)/2)),len+1])
tmin_rcp_all=cbind(proj_tminimum1[1:(nrow(proj_tminimum1)/2),],proj_tminimum1[(1+(nrow(proj_tminimum1)/2)):nrow(proj_tminimum1),],proj_tminimum[2:(1+(nrow(proj_tminimum1)/2)),len+1])


rcp=2
dat_1=seq(as.Date("2006/1/1"), as.Date("2012/12/31"), "day")
dat_2=seq(as.Date("2056/1/1"), as.Date("2065/12/31"), "day")
dat_3=seq(as.Date("2090/1/1"), as.Date("2099/12/31"), "day")
proj_dff1=t(matrix(unlist(lapply(dat_1, toNumerics)),nrow=3,ncol=length(dat_1)))
proj_dff2=t(matrix(unlist(lapply(dat_2, toNumerics)),nrow=3,ncol=length(dat_2)))
proj_dff3=t(matrix(unlist(lapply(dat_3, toNumerics)),nrow=3,ncol=length(dat_3)))

match_tmax1=match(dat_1,tmax_rcp_all[,rcp*len+1])
match_tmax2=match(dat_2,tmax_rcp_all[,rcp*len+1])
match_tmax3=match(dat_3,tmax_rcp_all[,rcp*len+1])

###Order of data
#rcp 26/....................rcp26 rcp45.....................rcp45 rcp60...............rcp60 rcp85..................rcp85
tmax_rcp_final=rbind(tmax_rcp_all[match_tmax1,],tmax_rcp_all[match_tmax2,],tmax_rcp_all[match_tmax3,])
tmin_rcp_final=rbind(tmin_rcp_all[match_tmax1,],tmin_rcp_all[match_tmax2,],tmin_rcp_all[match_tmax3,])


######for MPI
tmax_rcp_final=rbind(tmax_rcp_all[match_tmax2,],tmax_rcp_all[match_tmax3,])
tmin_rcp_final=rbind(tmin_rcp_all[match_tmax2,],tmin_rcp_all[match_tmax3,])


#########BIAS CORRECTION FOR TEMPERATURE#############
ex=rcp*len+1
proj_tmax1=cbind(proj_dff1,tmax_rcp_final[1:length(dat_1),-ex])
proj_tmax2=cbind(proj_dff2,tmax_rcp_final[(1+length(dat_1)):(length(dat_2)+length(dat_1)),-ex])
proj_tmax3=cbind(proj_dff3,tmax_rcp_final[(1+length(dat_2)+length(dat_1)):(length(dat_2)+length(dat_1)+length(dat_3)),-ex])

proj_tmin1=cbind(proj_dff1,tmin_rcp_final[1:length(dat_1),-ex])
proj_tmin2=cbind(proj_dff2,tmin_rcp_final[(1+length(dat_1)):(length(dat_2)+length(dat_1)),-ex])
proj_tmin3=cbind(proj_dff3,tmin_rcp_final[(1+length(dat_2)+length(dat_1)):(length(dat_2)+length(dat_1)+length(dat_3)),-ex])

###MPI
proj_tmax2=cbind(proj_dff2,tmax_rcp_final[1:length(dat_2),-ex])
proj_tmax3=cbind(proj_dff3,tmax_rcp_final[(1+length(dat_2)):(length(dat_2)+length(dat_3)),-ex])

proj_tmin2=cbind(proj_dff2,tmin_rcp_final[1:length(dat_2),-ex])
proj_tmin3=cbind(proj_dff3,tmin_rcp_final[(1+length(dat_2)):(length(dat_2)+length(dat_3)),-ex])



proj_tmax1=data.frame(proj_tmax1)
proj_tmax2=data.frame(proj_tmax2)
proj_tmax3=data.frame(proj_tmax3)

proj_tmin1=data.frame(proj_tmin1)
proj_tmin2=data.frame(proj_tmin2)
proj_tmin3=data.frame(proj_tmin3)

proj_bc_cu_tmax=proj_tmax1 ##2006-2012
proj_bc_cu_tmin=proj_tmin1 ##2006-2012

proj_bc_fu_tmax=rbind(proj_tmax2,proj_tmax3) 
proj_bc_fu_tmax_raw=rbind(proj_tmax2,proj_tmax3)
proj_bc_fu_tmin=rbind(proj_tmin2,proj_tmin3) 
proj_bc_fu_tmin_raw=rbind(proj_tmin2,proj_tmin3)

his_dym=t(matrix(unlist(lapply(dates_bc, toNumerics)),nrow=3,ncol=length(dates_bc)))
fut_dym=proj_bc_fu_tmax[,2:3]
dym=his_dym[367:731,2:3]

tmax_obs=data.frame(tmax_obs)
tmax_obs=data.frame(rep(tmax_obs,2))
tmax_sim=data.frame(his_tmaximum1)
tmax_sim=data.frame(rep(tmax_sim,2))


tmin_obs=data.frame(tmin_obs)
tmin_obs=data.frame(rep(tmin_obs,2))
tmin_sim=data.frame(his_tminimum1)
tmin_sim=data.frame(rep(tmin_sim,2))


dtr_obs=tmax_obs-tmin_obs
dtr_sim=tmax_sim-tmin_sim
dtr_proj=proj_bc_fu_tmax-proj_bc_fu_tmin
dtr_proj1=dtr_proj[1:3653,]
dtr_proj2=dtr_proj[3654:7305,]


for(i in 1:(rcp*len)){
  
  
  for ( j in 1:365){
    
    
    id_his=which(his_dym[,2]==dym[j,1] & his_dym[,3]==dym[j,2])
    obs=tmax_obs[id_his,i]
    sim=tmax_sim[id_his,i] 
    
    
    id_proj1=which(fut_dym[1:3653,1]==dym[j,1] & fut_dym[1:3653,2]==dym[j,2])
    id_proj2=which(fut_dym[3654:7305,1]==dym[j,1] & fut_dym[3654:7305,2]==dym[j,2])
    
    
    proj1=sort(proj_tmax2[id_proj1,i+3])
    proj2=sort(proj_tmax3[id_proj2,i+3])
    
    
    tmaxfit1 <- fitQmapRQUANT(sim,proj1,qstep=0.01,nboot=1,wet.day=FALSE)
    corr1 <- doQmapRQUANT(proj1,tmaxfit1,type="tricub")
    
    tmaxfit2 <- fitQmapRQUANT(obs,proj1,qstep=0.01,nboot=1,wet.day=FALSE)
    corr2<- doQmapRQUANT(proj1,tmaxfit2,type="tricub")
    
    dx1=corr2-corr1
    proj_tmax2[id_proj1,i+3]=proj_tmax2[id_proj1,i+3]+dx1
    
    
    tmaxfit3 <- fitQmapRQUANT(sim,proj2,qstep=0.01,nboot=1,wet.day=FALSE)
    corr3 <- doQmapRQUANT(proj2,tmaxfit3,type="tricub")
    
    tmaxfit4 <- fitQmapRQUANT(obs,proj2,qstep=0.01,nboot=1,wet.day=FALSE)
    corr4<- doQmapRQUANT(proj2,tmaxfit4,type="tricub")
    
    dx2=corr4-corr3
    proj_tmax3[id_proj2,i+3]=proj_tmax3[id_proj2,i+3]+dx2
    
    
    
  }
  print(i)
}

tmax_final=rbind(proj_tmax2,proj_tmax3)


for(i in 1:(rcp*len)){
  
  for ( j in 1:365){
    
    
    
    id_his=which(his_dym[,2]==dym[j,1] & his_dym[,3]==dym[j,2])
    obs=dtr_obs[id_his,i]
    sim=dtr_sim[id_his,i] 
    
    
    id_proj1=which(fut_dym[1:3653,1]==dym[j,1] & fut_dym[1:3653,2]==dym[j,2])
    id_proj2=which(fut_dym[3654:7305,1]==dym[j,1] & fut_dym[3654:7305,2]==dym[j,2])
    
    
    proj1=sort(dtr_proj1[id_proj1,i+3])
    proj2=sort(dtr_proj2[id_proj2,i+3])
    
    
    dtrfit1 <- fitQmapRQUANT(sim,proj1,qstep=0.01,nboot=1,wet.day=FALSE)
    corr1 <- doQmapRQUANT(proj1,dtrfit1,type="tricub")
    
    dtrfit2 <- fitQmapRQUANT(obs,proj1,qstep=0.01,nboot=1,wet.day=FALSE)
    corr2<- doQmapRQUANT(proj1,dtrfit2,type="tricub")
    
    dx1=corr2-corr1
    dtr_proj1[id_proj1,i+3]=dtr_proj1[id_proj1,i+3]+dx1
    
    
    dtrfit3 <- fitQmapRQUANT(sim,proj2,qstep=0.01,nboot=1,wet.day=FALSE)
    corr3 <- doQmapRQUANT(proj2,dtrfit3,type="tricub")
    
    dtrfit4 <- fitQmapRQUANT(obs,proj2,qstep=0.01,nboot=1,wet.day=FALSE)
    corr4<- doQmapRQUANT(proj2,dtrfit4,type="tricub")
    
    dx2=corr4-corr3
    dtr_proj2[id_proj2,i+3]=dtr_proj2[id_proj2,i+3]+dx2
    
    
  }
  print(i)
} 


dtr_final=rbind(dtr_proj1,dtr_proj2)
#
#plot(dtr_final[3654:7305,4])
#plot(dtr_proj[3654:7305,4])

tmin_final=tmax_final-dtr_final
#tmin_final=dtr_final


###temperature analysis######

tmax_base=(tmaxtmintdew_d[date_overlap_base,seq(2,dim(tmaxtmintdew_d)[2],3)])
tmin_base=(tmaxtmintdew_d[date_overlap_base,seq(3,dim(tmaxtmintdew_d)[2],3)])
tdew_base=(tmaxtmintdew_d[date_overlap_base,seq(4,dim(tmaxtmintdew_d)[2],3)])

tmax_dm_mon<- data.frame(monthlyfunction(tmax_base, FUN=mean, na.rm=TRUE,dates=dates_base))
tmin_dm_mon<- data.frame(monthlyfunction(tmin_base, FUN=mean, na.rm=TRUE,dates=dates_base))
tdew_dm_mon<- data.frame(monthlyfunction(tdew_base, FUN=mean, na.rm=TRUE,dates=dates_base))

ng=len
tmax_bcca_mon_corr_56_65<- data.frame(monthlyfunction(tmax_final[1:3653,seq(4,rcp*ng+3,1)], FUN=mean, na.rm=TRUE,dates=dat_2))
tmax_bcca_mon_corr_90_99<- data.frame(monthlyfunction(tmax_final[3654:7305,seq(4,rcp*ng+3,1)], FUN=mean, na.rm=TRUE,dates=dat_3))


tmin_bcca_mon_corr_56_65<- data.frame(monthlyfunction(tmin_final[1:3653,seq(4,rcp*ng+3,1)], FUN=mean, na.rm=TRUE,dates=dat_2))
tmin_bcca_mon_corr_90_99<- data.frame(monthlyfunction(tmin_final[3654:7305,seq(4,rcp*ng+3,1)], FUN=mean, na.rm=TRUE,dates=dat_3))


par(mfrow=c(2,2))

plot(colMeans(tmax_dm_mon),ylim=c(-5,35))
lines(colMeans(tmax_bcca_mon_corr_56_65[1:10,]))

plot(colMeans(tmin_dm_mon),ylim=c(-8,20))
lines(colMeans(tmin_bcca_mon_corr_56_65[1:10,]))


plot(colMeans(tmax_dm_mon),ylim=c(-5,35))
lines(colMeans(tmax_bcca_mon_corr_56_65[11:20,]))

plot(colMeans(tmin_dm_mon),ylim=c(-10,20))
lines(colMeans(tmin_bcca_mon_corr_56_65[11:20,]))



par(mfrow=c(2,2))
barplot(colMeans(tmax_bcca_mon_corr_56_65[(ng*1-ng+1):(ng*1),]),ylab="Monthly Rainfall(mm)", main="RCP2.6")
barplot(colMeans(tmax_bcca_mon_corr_56_65[(ng*2-ng+1):(ng*2),]),main="RCP4.5")
barplot((colMeans(tmax_bcca_mon_corr_56_65[(ng*1-ng+1):(ng*1),])-colMeans(tmax_dm_mon)),ylab="Rainfall difference(mm)",main="RCP2.6")
barplot(colMeans(tmax_bcca_mon_corr_56_65[(ng*2-ng+1):(ng*2),])-colMeans(tmax_dm_mon),main="RCP4.5")

par(mfrow=c(2,2))
barplot(colMeans(tmax_bcca_mon_corr_90_99[(ng*1-ng+1):(ng*1),]),ylab="Monthly Rainfall(mm)", main="RCP2.6")
barplot(colMeans(tmax_bcca_mon_corr_90_99[(ng*2-ng+1):(ng*2),]),main="RCP4.5")
barplot((colMeans(tmax_bcca_mon_corr_90_99[(ng*1-ng+1):(ng*1),])-colMeans(tmax_dm_mon)),ylab="Rainfall difference(mm)",main="RCP2.6")
barplot(colMeans(tmax_bcca_mon_corr_90_99[(ng*2-ng+1):(ng*2),])-colMeans(tmax_dm_mon),main="RCP4.5")

par(mfrow=c(2,2))
barplot((0.5*(colMeans(tmax_bcca_mon_corr_90_99[(ng*1-ng+1):(ng*1),])+colMeans(tmin_bcca_mon_corr_90_99[(ng*1-ng+1):(ng*1),]))),ylab=expression(paste("Avg Temperature " ^ 0,"C",sep="")),main="RCP4.5",ylim=c(-5,25))
barplot((0.5*(colMeans(tmax_bcca_mon_corr_90_99[(ng*2-ng+1):(ng*2),])+colMeans(tmin_bcca_mon_corr_90_99[(ng*2-ng+1):(ng*2),]))),main="RCP8.5",ylim=c(-5,25))



avgT_56_65_4=(0.5*(colMeans(tmax_bcca_mon_corr_56_65[(ng*1-ng+1):(ng*1),])+colMeans(tmin_bcca_mon_corr_56_65[(ng*1-ng+1):(ng*1),])))
avgT_56_65_8=(0.5*(colMeans(tmax_bcca_mon_corr_56_65[(ng*2-ng+1):(ng*2),])+colMeans(tmin_bcca_mon_corr_56_65[(ng*2-ng+1):(ng*2),])))

avgT_90_99_4=(0.5*(colMeans(tmax_bcca_mon_corr_90_99[(ng*1-ng+1):(ng*1),])+colMeans(tmin_bcca_mon_corr_90_99[(ng*1-ng+1):(ng*1),])))
avgT_90_99_8=(0.5*(colMeans(tmax_bcca_mon_corr_90_99[(ng*2-ng+1):(ng*2),])+colMeans(tmin_bcca_mon_corr_90_99[(ng*2-ng+1):(ng*2),])))


avgT_dm=0.5*(colMeans(tmax_dm_mon)+colMeans(tmin_dm_mon))


par(mfrow=c(2,2))
barplot((avgT_56_65_4-avgT_dm),ylab="change of Temperature(c)", main="RCP4.5 56-65")
barplot((avgT_56_65_8-avgT_dm),ylab="change of Temperature(c)", main="RCP8.5 56-65")
barplot((avgT_90_99_4-avgT_dm),ylab="change of Temperature(c)", main="RCP4.5 90-99")
barplot((avgT_90_99_8-avgT_dm),ylab="change of Temperature(c)", main="RCP8.5 90-99")


#####findign relation between tdew and tmin ###########

tdewobs=rowMeans(tdew_obs)
tminobs=rowMeans(tmin_obs)
tmaxobs=rowMeans(tmax_obs)
x=matrix(nrow=41,ncol=1)
z=matrix(nrow=41,ncol=1)

for ( i in 1:41){
  #plot(tdew_obs[,i],tmin_obs[,i])
  rltdewtmin=lm(tdew_obs[,i]~tmin_obs[,i])
  x[i]=rltdewtmin$coefficients[1]
  z[i]=rltdewtmin$coefficients[2]
  
}

rltdewtmin=lm(tdewobs ~ tmaxobs)

strDates=rbind(data.frame(date=dat_2),data.frame(date=dat_3))
gh=gsub("-","", strDates$date, fixed=TRUE)
day=data.frame(time=gh)
hour=rep.int(240000, nrow(day))

rcp_name=c("rcp4.5","rcp8.5")
nc=ng
tdew=tmin_final[,]*1.137+6.309
tdew1=tdew
tdew1[tdew>tmax_final]=tmin_final[tdew>tmax_final]


tmaxtmintdew=data.frame(matrix(NA,nrow=1,ncol=1))
for ( k in 1:(rcp*len)){
  tmaxtmintdew=cbind(tmaxtmintdew,tmax_final[,k+3],tmin_final[,k+3]-1,tmin_final[,k+3]*1.137+7)
}

temper=data.frame(tmaxtmintdew[,-1])

nc=ng
for ( j in 1:rcp){
  
  temperpp=data.frame(temper[,seq((nc*j*3-nc*3+1),(nc*j*3),1)])
  temperpp[] <- lapply(temperpp, function(.col){ if (is.numeric(.col)) return(sprintf("%8.2f",.col))else return(.col)}) 
  temperpp=data.frame(temperpp,day,hour)
  
  
  ##writing rain.dat file
  file=paste("tmaxtmintdew_",rcp_name[j],".dat",sep="")
  sink(file) 
  cat(sprintf("This file provides daily values of Tmax/Tmin/Tdew for each temperature station"),file=file,append=TRUE)
  cat("\n", file=file, append=TRUE)
  cat(sprintf("Temperature is provided in degrees Celcius"),file=file,append=TRUE)
  cat("\n", file=file, append=TRUE)
  sites=seq(1,(nc),1)
  cat(sprintf("%s %d ", "ver2",nc),file=file,append=TRUE) 
  cat(sprintf( "%d", sites),(sprintf( "%s", "Date Hour","\n")),file=file,append=TRUE) 
  cat("\n", file=file, append=TRUE)
  write.table(temperpp, file =file,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
  sink()
  
}




































