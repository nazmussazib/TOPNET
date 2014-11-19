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
Model=Mod[2]
setwd(paste("F:\\USU_Research_work_update_March_30_2014\\TOPNET PROJECT\\MODEL COMPARISON\\CLIMATE_CHANGE_ALLDATA\\",Model,sep=""))
pr.nc = nc_open("BCCA_0.125deg_tasmax_day_CNRM-CM5_rcp45_r1i1p1_20060101-20101231.nc")
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

for(i in 1:len){
  cat("filename",i,".csv", file="latlon.txt", sep="",append=TRUE)
  cat(",", file="latlon.txt", append=TRUE)
  cat(lat_lon_rg1[i,2], lat_lon_rg1[i,1], file="latlon.txt", sep=",",append=TRUE)
  cat(",", file="latlon.txt", append=TRUE)
  cat("ignore stuff",i, file="latlon.txt", sep="",append=TRUE)
  cat("\n", file="latlon.txt", append=TRUE)}

system("java -Xms512m -Xmx1024m -jar daymet_multiple_extraction.jar latlon.txt")


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


precipitation_d=matrix(ncol = nc, nrow = nrow(dss))
for ( k in 1:(nc)){
  precipitation_d[,k]=dss[,4*k-2]
  
}

precip_day=data.frame(precipitation_d)
dates_bc=seq(as.Date("1980/1/1"), as.Date("2005/12/31"), "day") 
date_overlap_day=match(dates_bc,dates_day)
pr_DAY=(precip_day[date_overlap_day,]) ##daymet rainfall for 1996-2005 for each station


#####processiing projected rainfall  data #############

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
rcp=2
proj_rain_file=list.files(path =paste("F:\\USU_Research_work_update_March_30_2014\\TOPNET PROJECT\\MODEL COMPARISON\\CLIMATE_CHANGE_ALLDATA\\",Model,sep=""),pattern =paste("*_pr_day_",Model,"_rcp",sep=""))
his_rain_file=list.files(path =paste("F:\\USU_Research_work_update_March_30_2014\\TOPNET PROJECT\\MODEL COMPARISON\\CLIMATE_CHANGE_ALLDATA\\",Model,sep=""),pattern =paste("*_pr_day_",Model,"_historical_",sep=""))


his_precp=data.frame(matrix(NA,nrow=1,ncol=ng+1))
for ( m in 1:length(his_rain_file)){
  pp=rain(his_rain_file[m])
  ne=colnames(pp)
  colnames(his_precp)[]=ne
  his_precp=rbind(his_precp,pp)
}

precp=data.frame(matrix(NA,nrow=1,ncol=ng+1))
for ( m in 1:length(proj_rain_file)){
  pp=rain(proj_rain_file[m])
  ne=colnames(pp)
  colnames(precp)[]=ne
  precp=rbind(precp,pp)
}


#######precp means projected rainfall and his_precp means hitorical rainfall #######
##reshape the data

##for ccsm4,canESm2 and MPI_ESMMR

his_precp1=as.matrix((his_precp[-1,1:len]))
his_precp_all=his_precp1[1:length(dates_bc),]

##for gfdl only

his_precp1=as.matrix((his_precp[-1,1:len]))
gfdl_bc=seq(as.Date("1976/1/1"), as.Date("2005/12/31"), "day") 
match_gfdl=match(dates_bc,gfdl_bc)
his_precp_all=his_precp1[match_gfdl,]


##for ccsm4,canESm2 
precp1=as.matrix((precp[-1,1:len]))
precp_rcp_all=cbind(precp1[1:(nrow(precp1)/2),],precp1[(1+nrow(precp1)/2):nrow(precp1),],precp[2:(1+nrow(precp1)/2),len+1])

##for MPI_ESMMR
precp1=as.matrix((precp[-1,1:len]))
precp_rcp_all=cbind(precp1[1:(nrow(precp1)/2),],precp1[(1+nrow(precp1)/2):nrow(precp1),],precp[2:(1+nrow(precp1)/2),len+1])

dat_1=seq(as.Date("2006/1/1"), as.Date("2012/12/31"), "day")
dat_2=seq(as.Date("2056/1/1"), as.Date("2065/12/31"), "day")
dat_3=seq(as.Date("2090/1/1"), as.Date("2099/12/31"), "day")

match_precp1=match(dat_1,precp_rcp_all[,rcp*len+1])
match_precp2=match(dat_2,precp_rcp_all[,rcp*len+1])
match_precp3=match(dat_3,precp_rcp_all[,rcp*len+1])

###Order of data
#rcp 26/....................rcp26 rcp45.....................rcp45 rcp60...............rcp60 rcp85..................rcp85
precp_rcp_final=rbind(precp_rcp_all[match_precp1,],precp_rcp_all[match_precp2,],precp_rcp_all[match_precp3,])
##for MPI_ESMMR
precp_rcp_final=rbind(precp_rcp_all[match_precp2,],precp_rcp_all[match_precp3,])

######quantile mapping bias correction method##########

toNumerics <- function(Date) {
  stopifnot(inherits(Date, c("Date", "POSIXt")))
  day <- as.numeric(strftime(Date, format = "%d"))
  month <- as.numeric(strftime(Date, format = "%m"))
  year <- as.numeric(strftime(Date, format = "%Y"))
  list(year = year, month = month, day = day)
}

dff=t(matrix(unlist(lapply(dates_bc, toNumerics)),nrow=3,ncol=length(dates_bc)))

proj_dff1=t(matrix(unlist(lapply(dat_1, toNumerics)),nrow=3,ncol=length(dat_1)))
proj_dff2=t(matrix(unlist(lapply(dat_2, toNumerics)),nrow=3,ncol=length(dat_2)))
proj_dff3=t(matrix(unlist(lapply(dat_3, toNumerics)),nrow=3,ncol=length(dat_3)))

ex=rcp*len+1
proj_rain1=cbind(proj_dff1,precp_rcp_final[1:length(dat_1),-ex])
proj_rain2=cbind(proj_dff2,precp_rcp_final[(1+length(dat_1)):(length(dat_2)+length(dat_1)),-ex])
proj_rain3=cbind(proj_dff3,precp_rcp_final[(1+length(dat_2)+length(dat_1)):(length(dat_2)+length(dat_1)+length(dat_3)),-ex])

##for MPI
proj_rain2=cbind(proj_dff2,precp_rcp_final[1:length(dat_2),-ex])
proj_rain3=cbind(proj_dff3,precp_rcp_final[(1+length(dat_2)):(length(dat_2)+length(dat_3)),-ex])

proj_rain1=data.frame(proj_rain1) ## for year 2006-2012
proj_rain2=data.frame(proj_rain2) ## for year 2056-2065
proj_rain3=data.frame(proj_rain3) ## for year 2090-2099

proj_bc_cu=proj_rain1 ##2006-2012
proj_bc_fu=rbind(proj_rain2,proj_rain3)          ## 2056-65 and 2090-99
#proj_bc_fu1=rbind(proj_rain2,proj_rain3) 

his_dym=t(matrix(unlist(lapply(dates_bc, toNumerics)),nrow=3,ncol=length(dates_bc)))
fut_dym=proj_bc_fu[,2:3]
dym=his_dym[367:731,2:3]

obs=data.frame(pr_DAY)
obs=data.frame(rep(obs,2))
sim=data.frame(his_precp_all)
sim=data.frame(rep(sim,2))


for(i in 1:(rcp*len)){
  
  for ( j in 1:365){
    
    
    
    id_his=which(his_dym[,2]==dym[j,1] & his_dym[,3]==dym[j,2])
    obs1=obs[id_his,i]
    sim1=sim[id_his,i] 
    id_proj1=which(fut_dym[1:3653,1]==dym[j,1] & fut_dym[1:3653,2]==dym[j,2])
    id_proj2=which(fut_dym[3654:7305,1]==dym[j,1] & fut_dym[3654:7305,2]==dym[j,2])
    
    proj1=proj_rain2[id_proj1,i+3]
    proj2=proj_rain3[id_proj2,i+3]
    
    if((length(which(sim1>0))>10) & (length(which(proj1>0))>3) & (length(which(obs1>0))>10)){
        rain_qm.fit1 <- fitQmapRQUANT(sim1,proj1,qstep=0.1,nboot=10,wet.day=0)
        his_corr1 <- doQmapRQUANT(proj1,rain_qm.fit1,type="tricub")
        
        rain_qm.fit2 <- fitQmapRQUANT(obs1,proj1,qstep=0.1,nboot=10,wet.day=0)
        his_corr2<- doQmapRQUANT(proj1,rain_qm.fit2,type="tricub")
        
        dx1=his_corr2/his_corr1
        proj_rain2[id_proj1,i+3]=proj_rain2[id_proj1,i+3]*dx1}else {proj_rain2[id_proj1,i+3]=proj_rain2[id_proj1,i+3]}
    
    if((length(which(sim1>0))>10) & (length(which(proj2>0))>3) & (length(which(obs1>0))>10)){
    rain_qm.fit3 <- fitQmapRQUANT(sim1,proj2,qstep=0.1,nboot=10,wet.day=0)
    his_corr3 <- doQmapRQUANT(proj2,rain_qm.fit3,type="tricub")
    
    rain_qm.fit4 <- fitQmapRQUANT(obs1,proj2,qstep=0.1,nboot=10,wet.day=0)
    his_corr4<- doQmapRQUANT(proj2,rain_qm.fit4,type="tricub")
    
    dx2=his_corr4/his_corr3
    proj_rain3[id_proj2,i+3]=proj_rain3[id_proj2,i+3]*dx2}else {proj_rain3[id_proj2,i+3]=proj_rain3[id_proj2,i+3]}
  
          
  }
print(i)
}

proj_bc_fu1=rbind(proj_rain2,proj_rain3)
proj_bc_fu1[proj_bc_fu1==Inf]=0
proj_bc_fu1[is.na(proj_bc_fu1)]=proj_bc_fu[is.na(proj_bc_fu1)]


############Analysis of Rainfall Results compare to base period 2003-2012 ##########

dates_day=seq(as.Date("1980/1/1"), as.Date("2012/12/31"), "day") 
dates_base=seq(as.Date("2003/1/1"), as.Date("2012/12/31"), "day") ##base line period for 2003-2012

date_overlap_base=match(dates_base,dates_day) 
pr_DAY_base=(precip_day[date_overlap_base,])
pr_daymet_mean_base=rowMeans(pr_DAY_base)
pr_daymet_mean_base=matrix(pr_daymet_mean_base)
rain_dm_mon_base<- (data.frame(monthlyfunction(pr_daymet_mean_base, FUN=sum, na.rm=TRUE,dates=dates_base)))/10  ##as base period is 10 years

rain_bcca_mon_corr_56_65<- data.frame(monthlyfunction(proj_bc_fu1[1:3653,seq(4,rcp*len+3,1)], FUN=sum, na.rm=TRUE,dates=dat_2))/10
rain_bcca_mon_corr_90_99<- data.frame(monthlyfunction(proj_bc_fu1[3654:7305,seq(4,rcp*len+3,1)], FUN=sum, na.rm=TRUE,dates=dat_3))/10

par(mfrow=c(2,2))
barplot(colMeans(rain_bcca_mon_corr_56_65[(ng*1-ng+1):(ng*1),]),ylab="Monthly Rainfall(mm)", main="RCP4.5 56-65",ylim=c(0,200))
barplot(colMeans(rain_bcca_mon_corr_56_65[(ng*2-ng+1):(ng*2),]),main="RCP8.5 56-65",ylim=c(0,200))
barplot(colMeans(rain_bcca_mon_corr_90_99[(ng*1-ng+1):(ng*1),]),ylab="Monthly Rainfall(mm)", main="RCP4.5 90-99",ylim=c(0,200))
barplot(colMeans(rain_bcca_mon_corr_90_99[(ng*2-ng+1):(ng*2),]),main="RCP8.5 90-99",ylim=c(0,200))


par(mfrow=c(2,2))
barplot((colMeans(rain_bcca_mon_corr_56_65[(ng*1-ng+1):(ng*1),])-data.matrix(rain_dm_mon_base)),ylab="Rainfall difference(mm)",main="RCP4.5 56-65",ylim=c(-25,50))
barplot(colMeans(rain_bcca_mon_corr_56_65[(ng*2-ng+1):(ng*2),])-data.matrix(rain_dm_mon_base),main="RCP8.5 56-65",ylim=c(-25,50))
barplot((colMeans(rain_bcca_mon_corr_90_99[(ng*1-ng+1):(ng*1),])-data.matrix(rain_dm_mon_base)),ylab="Rainfall difference(mm)",main="RCP4.5 90-99",ylim=c(-25,50))
barplot(colMeans(rain_bcca_mon_corr_90_99[(ng*2-ng+1):(ng*2),])-data.matrix(rain_dm_mon_base),main="RCP8.5 90-99",ylim=c(-25,50))





#############WRITING RAIN.DAT for DIFFERNECT RCP SCENARIO#################

strDates=rbind(data.frame(date=dat_2),data.frame(date=dat_3))
gh=gsub("-","", strDates$date, fixed=TRUE)
day=data.frame(time=gh)
hour=rep.int(240000, nrow(day))

rcp_name=c("rcp4.5","rcp8.5")
nc=ng

###change the directory for each watershed #########
setwd(paste("E:\\USU_Research_work\\TOPNET PROJECT\\MODEL COMPARISON\\",watsed[w],"_watershed_final\\",watsed[w],"_CC_data\\",Model,sep=""))##for A1:


precp_rcp_bias_corr=proj_bc_fu1[,seq(4,rcp*len+3,1)]
for ( j in 1:rcp){
  precipp=data.frame(precp_rcp_bias_corr[,seq((nc*j-nc+1),nc*j,1)])
  
  #precipp=data.frame(precp_rcp_final[,seq((nc*j-nc+1),nc*j,1)])
  
  ##change here as date of 2056-2065 and 2090-2099
  
  precipp[] <- lapply(precipp, function(.col){ if (is.numeric(.col)) return(sprintf("%8.2f",.col))else return(.col)}) 
  precipp=data.frame(precipp,day,hour)
  
  
  ##writing rain.dat file
  file=paste("rain_",rcp_name[j],".dat",sep="")
  sink(file) 
  cat(sprintf("This file provides daily precipitation rate for each rain station"),file=file,append=TRUE)
  cat("\n", file=file, append=TRUE)
  cat(sprintf("Precipitation is provided in mm/day"),file=file,append=TRUE)
  cat("\n", file=file, append=TRUE)
  sites=seq(1,(nc),1)
  cat(sprintf("%s %d ", "ver2",nc),file=file,append=TRUE) 
  cat(sprintf( "%d", sites),(sprintf( "%s", "Date Hour","\n")),file=file,append=TRUE) 
  cat("\n", file=file, append=TRUE)
  write.table(precipp, file =file,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
  sink()
  
}


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




















