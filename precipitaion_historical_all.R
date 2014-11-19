#############precipitaion correction$##############
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
pr.nc = nc_open("BCCA_0.125deg_tasmax_day_CCSM4_historical_r2i1p1_19800101-19891231.nc")
latall=ncvar_get(pr.nc ,'latitude')
lonall=ncvar_get(pr.nc ,'longitude')


###change the directory for each watershed #########
watsed=c("A1","A2","B1","B21","B22","C1","C21","C22")
w=1
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


######server is really slow and problemw with downloading data so extract information from the 
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

hist_dat=seq(as.Date("1980/1/1"), as.Date("1995/12/31"), "day") 
proj_dat=seq(as.Date("1996/1/1"), as.Date("2005/12/31"), "day") 
dates_bc=seq(as.Date("1980/1/1"), as.Date("2005/12/31"), "day") 
precip_day=data.frame(precipitation_d)
pr_obs_his=precip_day[1:length(hist_dat),]
pr_obs_base=precip_day[(1+length(hist_dat)):length(dates_bc),]


tmax_obs_his=(tmaxtmintdew_d[1:length(hist_dat),seq(2,dim(tmaxtmintdew_d)[2],3)])
tmax_obs_base=(tmaxtmintdew_d[(1+length(hist_dat)):length(dates_bc),seq(2,dim(tmaxtmintdew_d)[2],3)])

tmin_obs_his=(tmaxtmintdew_d[1:length(hist_dat),seq(3,dim(tmaxtmintdew_d)[2],3)])
tmin_obs_base=(tmaxtmintdew_d[(1+length(hist_dat)):length(dates_bc),seq(3,dim(tmaxtmintdew_d)[2],3)])

tdew_obs_his=(tmaxtmintdew_d[1:length(hist_dat),seq(4,dim(tmaxtmintdew_d)[2],3)])
tdew_obs_base=(tmaxtmintdew_d[(1+length(hist_dat)):length(dates_bc),seq(4,dim(tmaxtmintdew_d)[2],3)])

#####processiing historical rainfall  data #############

setwd(paste("F:\\USU_Research_work_update_March_30_2014\\TOPNET PROJECT\\MODEL COMPARISON\\CLIMATE_CHANGE_ALLDATA\\",Model,sep=""))

rain <- function(file ){
  
  rain.nc = nc_open(file)
  rain_time=as.Date(ncvar_get(rain.nc,'time'),origin="1950-01-01")
  #strDates= as.character(rain_time)
  ra=matrix(NA,nrow=length( rain_time),ncol=len)
  for (i in 1:len){
    start=c(lat_lon_ID[i,1],lat_lon_ID[i,2],1)
    count=c(1,1,length( rain_time))
    ra[,i]=ncvar_get(rain.nc,'pr',start,count)
  }
  radata=data.frame(data=cbind(ra, rain_time-0.5))
  
  return(radata)
}

ng=len
his_rain_file=list.files(path =paste("F:\\USU_Research_work_update_March_30_2014\\TOPNET PROJECT\\MODEL COMPARISON\\CLIMATE_CHANGE_ALLDATA\\",Model,sep=""),pattern =paste("*_pr_day_",Model,"_historical_",sep=""))
his_precp=data.frame(matrix(NA,nrow=1,ncol=ng+1))
for ( m in 1:length(his_rain_file)){
  pp=rain(his_rain_file[m])
  ne=colnames(pp)
  colnames(his_precp)[]=ne
  his_precp=rbind(his_precp,pp)
}


his_precp1=as.matrix((his_precp[-1,1:len]))
pr_sim_his=his_precp1[1:length(hist_dat),]
pr_sim_base=his_precp1[(1+length(hist_dat)):length(dates_bc),]

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

his_tmaximum=data.frame(matrix(NA,nrow=1,ncol=len+1))
his_tminimum=data.frame(matrix(NA,nrow=1,ncol=len+1))

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
his_tmaximum1=as.matrix((his_tmaximum[-1,1:len]))
his_tminimum1=as.matrix((his_tminimum[-1,1:len]))

tmax_sim_his=his_tmaximum1[1:length(hist_dat),]
tmax_sim_base=his_tmaximum1[(1+length(hist_dat)):length(dates_bc),]

tmin_sim_his=his_tminimum1[1:length(hist_dat),]
tmin_sim_base=his_tminimum1[(1+length(hist_dat)):length(dates_bc),]


######quantile mapping bias correction method##########

toNumerics <- function(Date) {
  stopifnot(inherits(Date, c("Date", "POSIXt")))
  day <- as.numeric(strftime(Date, format = "%d"))
  month <- as.numeric(strftime(Date, format = "%m"))
  year <- as.numeric(strftime(Date, format = "%Y"))
  list(year = year, month = month, day = day)
}
proj_dff1=t(matrix(unlist(lapply(proj_dat, toNumerics)),nrow=3,ncol=length(proj_dat)))
proj_rain1=cbind(proj_dff1,pr_sim_base)
proj_rain1=data.frame(proj_rain1)

proj_tmax1=cbind(proj_dff1,tmax_sim_base)
proj_tmax1=data.frame(proj_tmax1)

proj_tmin1=cbind(proj_dff1,tmin_sim_base)
proj_tmin1=data.frame(proj_tmin1)
proj_dtr1=proj_tmax1-proj_tmin1

his_dym=t(matrix(unlist(lapply(hist_dat, toNumerics)),nrow=3,ncol=length(hist_dat)))
fut_dym=proj_rain1[,2:3]
dym=his_dym[367:731,2:3]

rainobs=data.frame(pr_obs_his)
rainsim=data.frame(pr_sim_his)

tmaxobs=data.frame(tmax_obs_his)
tmaxsim=data.frame(tmax_sim_his)

tminobs=data.frame(tmin_obs_his)
tminsim=data.frame(tmin_sim_his)

dtrobs=tmaxobs-tminobs
dtrsim=tmaxsim-tminsim


###rainfall correction###########

for(i in 1:(len)){
  
  for ( j in 1:365){
    
    
    
    id_his=which(his_dym[,2]==dym[j,1] & his_dym[,3]==dym[j,2])
    obs1=rainobs[id_his,i]
    sim1=rainsim[id_his,i] 
    id_proj1=which(fut_dym[,1]==dym[j,1] & fut_dym[,2]==dym[j,2])
    obs1=sort(obs1)
    sim1=sort(sim1)
   
    
    proj1=proj_rain1[id_proj1,i+3]
    proj=sort(proj1)
    
    if((length(which(sim1>0))>7) & (length(which(proj>0))>5) & (length(which(obs1>0))>7)){
      rain_qm.fit1 <- fitQmapRQUANT(sim1,proj,qstep=0.1,nboot=1,wet.day=0)
      his_corr1 <- doQmapRQUANT(proj1,rain_qm.fit1,type="tricub")
      
      rain_qm.fit2 <- fitQmapRQUANT(obs1,proj,qstep=0.1,nboot=1,wet.day=0)
      his_corr2<- doQmapRQUANT(proj1,rain_qm.fit2,type="tricub")
      
      dx1=his_corr2-his_corr1
      proj_rain1[id_proj1,i+3]=proj_rain1[id_proj1,i+3]+dx1}else {proj_rain1[id_proj1,i+3]=proj_rain1[id_proj1,i+3]}
   
    
  }
  print(i)
}

his_rain1=cbind(proj_dff1,pr_sim_base)
his_rain1=data.frame(his_rain1)
proj_bc_fu1=proj_rain1
proj_bc_fu1[proj_bc_fu1<0]=his_rain1[proj_bc_fu1<0]

###cumsum
plot(cumsum(rowMeans(pr_obs_base)))
lines(cumsum(rowMeans(proj_bc_fu1[,4:16])),col='green')
lines(cumsum(rowMeans(pr_sim_base)),col='red')


###temperature correction############
for(i in 1:(len)){
  
  
  for ( j in 1:365){
    
    
    id_his=which(his_dym[,2]==dym[j,1] & his_dym[,3]==dym[j,2])
    obs1=sort(tmaxobs[id_his,i])
    sim1=sort(tmaxsim[id_his,i])
    obs2=sort(dtrobs[id_his,i])
    sim2=sort(dtrsim[id_his,i])
    
    
    id_proj1=which(fut_dym[,1]==dym[j,1] & fut_dym[,2]==dym[j,2])
    tmproj=proj_tmax1[id_proj1,i+3]
    dtrproj=proj_dtr1[id_proj1,i+3]
    tmaxproj1=sort(tmproj)
    dtrproj1=sort(dtrproj)
    
    
    tmaxfit1 <- fitQmapRQUANT(sim1,tmaxproj1,qstep=0.01,nboot=1,wet.day=FALSE)
    corr1 <- doQmapRQUANT(tmproj,tmaxfit1,type="tricub")
    
    tmaxfit2 <- fitQmapRQUANT(obs1,tmaxproj1,qstep=0.01,nboot=1,wet.day=FALSE)
    corr2<- doQmapRQUANT(tmproj,tmaxfit2,type="tricub")
    
    dx1=corr2-corr1
    proj_tmax1[id_proj1,i+3]=proj_tmax1[id_proj1,i+3]+dx1
    
    
    tmaxfit3 <- fitQmapRQUANT(sim2,dtrproj1,qstep=0.01,nboot=1,wet.day=FALSE)
    corr3 <- doQmapRQUANT(dtrproj,tmaxfit3,type="tricub")
    
    tmaxfit4 <- fitQmapRQUANT(obs2,dtrproj1,qstep=0.01,nboot=1,wet.day=FALSE)
    corr4<- doQmapRQUANT(dtrproj,tmaxfit4,type="tricub")
    
    dx2=corr4-corr3
    proj_dtr1[id_proj1,i+3]=proj_dtr1[id_proj1,i+3]+dx2
    
    
    
  }
  print(i)
}

tmax_final=proj_tmax1
dtr_final=proj_dtr1
tmin_final=tmax_final-dtr_final


###estimating tdew from tmin #########


tmax_dm_mon<- data.frame(monthlyfunction(tmax_obs_base, FUN=mean, na.rm=TRUE,dates=proj_dat))
tmin_dm_mon<- data.frame(monthlyfunction(tmin_obs_base, FUN=mean, na.rm=TRUE,dates=proj_dat))
tdew_dm_mon<- data.frame(monthlyfunction(tdew_obs_base, FUN=mean, na.rm=TRUE,dates=proj_dat))

tmax_sim_mon<- data.frame(monthlyfunction(tmax_sim_base, FUN=mean, na.rm=TRUE,dates=proj_dat))
tmin_sim_mon<- data.frame(monthlyfunction(tmin_sim_base, FUN=mean, na.rm=TRUE,dates=proj_dat))

tmax_sim_mon_cor<- data.frame(monthlyfunction(tmax_final[,4:16], FUN=mean, na.rm=TRUE,dates=proj_dat))
tmin_sim_mon_cor<- data.frame(monthlyfunction(tmin_final[,4:16], FUN=mean, na.rm=TRUE,dates=proj_dat))






par(mfrow=c(2,2))

plot(colMeans(tmax_dm_mon),ylim=c(-5,35))
lines(colMeans(tmax_sim_mon))
lines(colMeans(tmax_sim_mon_cor),col='red')

plot(colMeans(tmin_dm_mon),ylim=c(-15,25))
lines(colMeans(tmin_sim_mon))
lines(colMeans(tmin_sim_mon_cor),col='red')


######correcting tdew########
plot(tmin_obs_his ~tdew_obs_his)

x1=matrix(nrow=13,ncol=1)
x2=matrix(nrow=13,ncol=1)
for ( i in 1:13) {
rltmd=lm(tdew_obs_his[,i] ~ tmin_obs_his[,i])
x1[i]=rltmd$coefficients[1]
x2[i]=rltmd$coefficients[2]
}

tdew_sim_base=tmin_final[,4:16]*x2+x1




tmaxtmintdew=data.frame(matrix(NA,nrow=1,ncol=1))
for ( k in 1:len){
  tmaxtmintdew=cbind(tmaxtmintdew,tmax_final[,k+3],tmin_final[,k+3],tdew_sim_base)
}

temper=data.frame(tmaxtmintdew[,-1])

precip=data.frame(proj_bc_fu1[,4:16])
temper=data.frame(tmaxtmintdew[,-1])

precip[] <- lapply(precip, function(.col){ if (is.numeric(.col)) return(sprintf("%8.2f",.col))else return(.col)}) 
temper[] <- lapply(temper, function(.col){ if (is.numeric(.col)) return(sprintf("%8.2f",.col))else return(.col)})
strDates=rbind(data.frame(date=proj_dat))
gh=gsub("-","", strDates$date, fixed=TRUE)
day=data.frame(time=gh)
hour=rep.int(240000, nrow(day))
precip=data.frame(precip,day,hour)
temper=data.frame(temper,day,hour)


##writing rain.dat file

sink("rain.dat") 
cat(sprintf("This file provides daily precipitation rate for each rain station"),file='rain.dat',append=TRUE)
cat("\n", file="rain.dat", append=TRUE)
cat(sprintf("Precipitation is provided in mm/day"),file='rain.dat',append=TRUE)
cat("\n", file="rain.dat", append=TRUE)
sites=seq(1,(nc),1)
cat(sprintf("%s %d ", "ver2",nc),file='rain.dat',append=TRUE) 
cat(sprintf( "%d", sites),(sprintf( "%s", "Date Hour","\n")),file='rain.dat',append=TRUE) 
cat("\n", file="rain.dat", append=TRUE)
write.table(precip, file = "rain.dat",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
sink()






##writing tmaxtmintdew.dat file

sink('tmaxtmintdew.dat')
cat(sprintf("This file provides daily values of Tmax/Tmin/Tdew for each temperature station"),file='tmaxtmintdew.dat',append=TRUE)
cat("\n", file="tmaxtmintdew.dat", append=TRUE)
cat(sprintf("Temperature is provided in degrees Celcius"),file='tmaxtmintdew.dat',append=TRUE)
cat("\n", file="tmaxtmintdew.dat", append=TRUE)
sites=seq(1,(nc),1)
cat(sprintf("%s %d ", "ver2",nc),file='tmaxtmintdew.dat',append=TRUE) 
cat(sprintf( "%d", sites),(sprintf( "%s", "Date Hour","\n")),file='tmaxtmintdew.dat',append=TRUE) 
cat("\n", file="tmaxtmintdew.dat", append=TRUE)
write.table(temper, file = "tmaxtmintdew.dat",row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
sink()
























































