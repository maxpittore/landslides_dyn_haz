library(ncdf4)
library(reshape2)
library(dplyr)
library(ggplot2)
require(ggmap)
library(viridis)
library(gganimate)
library(raster)
#library(rasterVis)
library(maptools)
#library(stpp)
library(readr)
library(raster)
#library(stpp)
#library(spatstat)
#library(beepr)
library(pracma)


#raster::plot(nc_brick[[1]], main = "ERA-5 Reanalysis Precipitation")
#maps::map("world", add = TRUE)

daily_sum<-function(rasterbrick)
{
  indices<-as.numeric(format(as.Date(names(rasterbrick), format = "X%Y.%m.%d.%H"), format = "%Y%m%d"))
  indices<-indices-indices[1]+1
  rasterbrick_daily = raster::stackApply(rasterbrick,indices,fun=sum,na.rm=TRUE)
  rasterbrick_daily<-setNames(rasterbrick_daily,unique(format(as.Date(names(rasterbrick), format = "X%Y.%m.%d.%H"), format = "%Y-%m-%d")))
  return (rasterbrick_daily)
}

#rast<-nc_runoff_daily
#rast<-nc_precip_daily
#cell_row<-19
#cell_col<-9
#t_start<-733881
#t_end<-734674
proc_cell<-function(rast, cell_row,cell_col,thr,t_start=-1,t_end=-1) {
  #get all layers
  t <-getValuesBlock(rast,cell_row,1,cell_col,1,format='v')
  dm <-as.Date(dimnames(t)[[2]],format="X%Y.%m.%d")
  names(t)<-dm
  #filter based on time range
  if (t_start >= 0 & t_end >= 0) {
    #dm <-as.Date(dimnames(t)[[2]],format="X%Y.%m.%d")
    t<-t[(dm >= as.Date(t_start,origin='0-01-01')) &
            (dm <= as.Date(t_end,origin='0-01-01'))]
    #names(t)<-dm[(dm >= as.Date(t_start,origin='0-01-01')) &
    #              (dm <= as.Date(t_end,origin='0-01-01'))]
  }
  #t<-t[,600:700]
  t_mask<-as.numeric(t>thr)
  #compute difference
  dtt<-diff(t_mask)
  #get transition up
  dttup<-which(dtt>0)
  if (isempty(dttup)) {
    return()
  }
  dtu_df<-as.data.frame(dttup,make.names=TRUE)
  dtu_df$dir<-1
  names(dtu_df)<-c("ind","dir")
  dtu_df$date<-rownames(dtu_df)
  #get transition down
  dttdown<-which(dtt<0)
  dtd_df<-as.data.frame(dttdown)
  dtd_df$dir<--1
  names(dtd_df)<-c("ind","dir")
  dtd_df$date<-rownames(dtd_df)
  #merge transitions 
  mm<-merge(dtu_df,dtd_df,all=TRUE)
  #get starting index of the sequences over threshold 
  start<-which(mm$dir==1)
  #get a copy of the array values
  ttcopy<-t
  #init a dataframe for the wet spell transitions
  ntrans<-size(unique(mm$date))[2] #number of transitions
  trans_df<-data.frame(index=numeric(ntrans),start=numeric(ntrans),
                       end=numeric(ntrans),duration=numeric(ntrans),cumval=numeric(ntrans))
  transind<-1
  #substitute values with cumulative sums
  for (i in start) {
    #i<-start[1]
    if (i>=nrow(mm)) {
      break
    }
    ttcopy[(mm[i,]$ind+1):mm[i+1,]$ind]<-cumsum(t[(mm[i,]$ind+1):mm[i+1,]$ind])
    #print(cumsum(tt[(mm[i,]$ind+1):mm[i+1,]$ind]))
    trans_df$index[transind]<-transind
    trans_df$start[transind]<-mm[i,]$ind+1
    trans_df$end[transind]<-mm[i+1,]$ind
    trans_df$duration[transind]<-mm[i+1,]$ind-mm[i,]$ind+1
    trans_df$cumval[transind]<-sum(t[(mm[i,]$ind+1):mm[i+1,]$ind])
    transind<-transind+1
  }
  #apply mask
  ttcopy_masked<-ttcopy*t_mask
  t_df<-data.frame("date"=as.Date(names(ttcopy_masked)), 
                   "value"=as.numeric(ttcopy_masked), row.names=NULL)
#  t_df<-data.frame("date"=as.Date(dimnames(ttcopy_masked)[[2]],format='X%Y.%m.%d'), 
#                   "value"=as.numeric(ttcopy_masked), row.names=NULL)
  
  #return masked cumsum and list of transitions
  return(list("daily_vals"=t_df,"transitions"=trans_df))
}


proc_cell(nc_precip_daily,1,1,0.1)

#for a given ROI (extent) compute the "events" (wet spell in case of precipitation,
#melt-spell in case of snowmelt) and populate a dataframe of events
process_roi<-function(rois,roi_ind,daily_rast,min_thr=0.1) {
  #rois<-as_3_3_4
  #roi_ind<-1
  #daily_rast<-nc_precip_daily
  #min_thr=0.1
  cells<-cellsFromExtent(daily_rast,
                         extent(rois[roi_ind,]$xmin,rois[roi_ind,]$xmax,rois[roi_ind,]$ymin,rois[roi_ind,]$ymax))
  rc_roi<-rowColFromCell(daily_rast,cells)
  joint_trans<-proc_cell(daily_rast, rc_roi[1,1],rc_roi[1,2],min_thr,rois[roi_ind,]$tmin,rois[roi_ind,]$tmax)$transitions
  for (i in 2:nrow(rc_roi)) {
    res_p<-proc_cell(daily_rast, rc_roi[i,1],rc_roi[i,2],min_thr,rois[roi_ind,]$tmin,rois[roi_ind,]$tmax)
    joint_trans<-rbind(joint_trans,res_p$transitions)
  }
  if (isempty(joint_trans)) {
    return()
  }
  joint_trans$roi<-roi_ind
  return(joint_trans)
}
###########################################################################################

proj_utm43<-"+proj=utm +zone=43 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

#inpath = '/Users/pittore/Documents/workspace/ERA5/data'
ncfile_precip = '/Users/pittore/Documents/workspace/ERA5/data/2009-2018_totalprecipitation.nc'
nc_precip_daily <- daily_sum(raster::brick(ncfile_precip))
curproj<-projection(nc_precip_daily)
#nc_precip_daily<-projectRaster(nc_precip_daily, crs=proj_utm43)
raster::plot(nc_precip_daily[[1]], main = "ERA-5 Reanalysis Precipitation")

#nc<-nc_open(ncfile_precip)
#print(nc)

ncfile_snowmelt = '/Users/pittore/Documents/workspace/ERA5/data/2009-2018_snowmelt.nc'
nc_snowmelt_daily <- daily_sum(raster::brick(ncfile_snowmelt))
raster::plot(nc_snowmelt_daily[[5]], main = "ERA-5 Reanalysis Snow melt")

ncfile_runoff = '/Users/pittore/Documents/workspace/ERA5/data/2009-2018_surfacerunoff.nc'
nc_runoff_daily <- daily_sum(raster::brick(ncfile_runoff))
raster::plot(nc_runoff_daily[[100]], main = "ERA-5 Reanalysis Runoff")

ncfile_subsurfrunoff = '/Users/pittore/Documents/workspace/ERA5/data/2009-2018_subsurfacerunoff.nc'
nc_subsurfrunoff_daily <- daily_sum(raster::brick(ncfile_subsurfrunoff))
raster::plot(nc_subsurfrunoff_daily[[100]], main = "ERA-5 Reanalysis Subsurface Runoff")

'testp<-"/Users/pittore/Documents/PROJECTS/WorldBank/hazard/vs30/deliverable/vs30_usgs_utm43.tif"
testp_r<-raster(testp)
str(testp_r)
projection(testp_r)

testp2<-"/Users/pittore/Documents/workspace/ERA5/data/test_prec_utm43N.tif"
testp2_r<-raster(testp2)
str(testp2_r)
raster::plot(testp2_r, main = "ERA-5 Reanalysis Precipitation")
projection(testp2_r)
'
nc_precip_daily_df<-data.frame(rasterToPoints(nc_precip_daily[[100]]))
nc_precip_daily_df<-setNames(nc_precip_daily_df,c("x","y","value"))

#### PROCESS AREAS NO 1 #####################################################
#read area sets from csv file
as_3_3_1<-read_csv("/Users/pittore/Documents/ACTIVITIES/landslides/CATENA/areaset_3x3x1_wgs84.csv")
str(as_3_3_1)

#roi_ind<-7
sps<-SpatialPoints(cbind(c(as_3_3_1[roi_ind,]$xmin,as_3_3_1[roi_ind,]$xmax),
      c(as_3_3_1[roi_ind,]$ymin,as_3_3_1[roi_ind,]$ymax)),CRS(proj_utm43))
extent(sps)
#sps2<-spTransform(sps,curproj)
sps_df<-data.frame(x1=c(extent(sps)@xmin),x2=c(extent(sps)@xmax),y1=c(extent(sps)@ymin),y2=c(extent(sps)@ymax))
#test_rast<-crop(nc_precip_daily[[1]],testext)

#plot raster and ROIs
ggplot()+geom_tile(data=nc_precip_daily_df,aes(x=x,y=y,fill=value))+
 geom_rect(data=as_3_3_1,mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=min_count),colour='red',alpha=0.5)

#plot raster and ROIs - alternative
ggplot(data=as_3_3_1)+geom_rect(mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=min_count),colour='red',alpha=0.5) + 
  geom_tile(nc_precip_daily_df,mapping=aes(x=x,y=y,fill=value),alpha=0.5)

#---- precipitations -----
pr_thr<-0.1
joint_trans_precip<-process_roi(as_3_3_1,1,nc_precip_daily,min_thr = pr_thr)
for (roi_ind in 2:nrow(as_3_3_1)) {
  print(sprintf("processing roi %d...",roi_ind))
  tmp_df<-process_roi(as_3_3_1,roi_ind,nc_precip_daily,min_thr = pr_thr)
  joint_trans_precip<-rbind(joint_trans_precip,tmp_df)
}
joint_trans_precip$type<-'wetspell'
ggplot(joint_trans_precip,aes(x=duration,y=cumval))+
  geom_point(alpha=0.5)+
  facet_wrap(~roi)+geom_smooth(method='lm')

# ---- snow ------
sm_thr<-0.1
joint_trans_snow<-process_roi(as_3_3_1,1,nc_snowmelt_daily,min_thr = sm_thr)
for (roi_ind in 2:nrow(as_3_3_1)) {
  print(sprintf("processing roi %d...",roi_ind))
  tmp_df<-process_roi(as_3_3_1,roi_ind,nc_snowmelt_daily,min_thr = sm_thr)
  joint_trans_snow<-rbind(joint_trans_snow,tmp_df)
}
joint_trans_snow$type<-'snowmelt'
ggplot(joint_trans_snow,aes(x=duration,y=cumval))+
  geom_point(alpha=0.5)+
  facet_wrap(~roi)+geom_smooth()

# --- surface runoff ---
sr_thr<-0.1
joint_trans_runoff<-process_roi(as_3_3_1,1,nc_runoff_daily,min_thr = sr_thr)
for (roi_ind in 2:nrow(as_3_3_1)) {
  print(sprintf("processing roi %d...",roi_ind))
  tmp_df<-process_roi(as_3_3_1,roi_ind,nc_runoff_daily,min_thr = sr_thr)
  joint_trans_runoff<-rbind(joint_trans_runoff,tmp_df)
}
joint_trans_runoff$type<-'runoff'
ggplot(joint_trans_runoff,aes(x=duration,y=cumval))+
  geom_point(alpha=0.5)+
  facet_wrap(~roi)+geom_smooth()

# --- sub-surface runoff ---
ssr_thr<-0.1
joint_trans_subsurfrunoff<-process_roi(as_3_3_1,1,nc_subsurfrunoff_daily,min_thr = ssr_thr)
for (roi_ind in 2:nrow(as_3_3_1)) {
  print(sprintf("processing roi %d...",roi_ind))
  #roi_ind<-5
  tmp_df<-process_roi(as_3_3_1,roi_ind,nc_subsurfrunoff_daily,min_thr = ssr_thr)
  joint_trans_subsurfrunoff<-rbind(joint_trans_subsurfrunoff,tmp_df)
}
joint_trans_subsurfrunoff$type<-'subsurf-runoff'
ggplot(joint_trans_subsurfrunoff,aes(x=duration,y=cumval))+
  geom_point(alpha=0.5)+
  facet_wrap(~roi)+geom_smooth()


joint_df<-joint_trans_snow
joint_df<-rbind(joint_df,joint_trans_precip)
joint_df<-rbind(joint_df,joint_trans_runoff)
joint_df<-rbind(joint_df,joint_trans_subsurfrunoff)

ggplot(joint_df,aes(x=duration,y=cumval,colour=type))+
  geom_point(alpha=0.5)+
  facet_wrap(~roi)#+geom_smooth()

#### PROCESS AREAS - TEST NO 2 #####################################################
#read area sets from csv file
as_3_3_4<-read_csv("/Users/pittore/Documents/ACTIVITIES/landslides/CATENA/areaset_3x3x4_wgs84.csv")
str(as_3_3_4)

#roi_ind<-7
sps<-SpatialPoints(cbind(c(as_3_3_4[roi_ind,]$xmin,as_3_3_4[roi_ind,]$xmax),
                         c(as_3_3_4[roi_ind,]$ymin,as_3_3_4[roi_ind,]$ymax)),CRS(proj_utm43))
extent(sps)
#sps2<-spTransform(sps,curproj)
sps_df<-data.frame(x1=c(extent(sps)@xmin),x2=c(extent(sps)@xmax),y1=c(extent(sps)@ymin),y2=c(extent(sps)@ymax))
#test_rast<-crop(nc_precip_daily[[1]],testext)

#plot raster and ROIs
ggplot()+geom_tile(data=nc_precip_daily_df,aes(x=x,y=y,fill=value))+
  geom_rect(data=as_3_3_4,mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=min_count),colour='red',alpha=0.5)

#plot raster and ROIs - alternative
ggplot(data=as_3_3_4)+geom_rect(mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=min_count),colour='red',alpha=0.5) + 
  geom_tile(nc_precip_daily_df,mapping=aes(x=x,y=y,fill=value),alpha=0.5)

#---- precipitations -----
pr_thr<-0.1
joint_trans_precip2<-process_roi(as_3_3_4,1,nc_precip_daily,min_thr = pr_thr)
for (roi_ind in 2:nrow(as_3_3_4)) {
  print(sprintf("processing roi %d...",roi_ind))
  tmp_df<-process_roi(as_3_3_4,roi_ind,nc_precip_daily,min_thr = pr_thr)
  joint_trans_precip2<-rbind(joint_trans_precip2,tmp_df)
}
joint_trans_precip2$type<-'wetspell'
ggplot(joint_trans_precip2,aes(x=duration,y=cumval))+
  geom_point(alpha=0.5)+
  facet_wrap(~roi,ncol=4)+geom_smooth(method='lm')

#ggplot(joint_trans_precip2)+ geom_histogram(aes(x=cumval))

# ---- snow ------
sm_thr<-0.1
joint_trans_snow2<-process_roi(as_3_3_4,1,nc_snowmelt_daily,min_thr = sm_thr)
for (roi_ind in 2:nrow(as_3_3_4)) {
  print(sprintf("processing roi %d...",roi_ind))
  tmp_df<-process_roi(as_3_3_4,roi_ind,nc_snowmelt_daily,min_thr = sm_thr)
  joint_trans_snow2<-rbind(joint_trans_snow2,tmp_df)
}
joint_trans_snow2$type<-'snowmelt'
ggplot(joint_trans_snow2,aes(x=duration,y=cumval))+
  geom_point(alpha=0.5)+
  facet_wrap(~roi,ncol=4)+geom_smooth()

# --- surface runoff ---
sr_thr<-0.1
joint_trans_runoff2<-process_roi(as_3_3_4,1,nc_runoff_daily,min_thr = sr_thr)
for (roi_ind in 2:nrow(as_3_3_4)) {
  print(sprintf("processing roi %d...",roi_ind))
  tmp_df<-process_roi(as_3_3_4,roi_ind,nc_runoff_daily,min_thr = sr_thr)
  joint_trans_runoff2<-rbind(joint_trans_runoff2,tmp_df)
}
joint_trans_runoff2$type<-'runoff'
ggplot(joint_trans_runoff2,aes(x=duration,y=cumval))+
  geom_point(alpha=0.5)+
  facet_wrap(~roi,ncol=4)+geom_smooth()

# --- sub-surface runoff ---
ssr_thr<-0.1
joint_trans_subsurfrunoff2<-process_roi(as_3_3_4,1,nc_subsurfrunoff_daily,min_thr = ssr_thr)
for (roi_ind in 2:nrow(as_3_3_4)) {
  print(sprintf("processing roi %d...",roi_ind))
  #roi_ind<-5
  tmp_df<-process_roi(as_3_3_4,roi_ind,nc_subsurfrunoff_daily,min_thr = ssr_thr)
  joint_trans_subsurfrunoff2<-rbind(joint_trans_subsurfrunoff2,tmp_df)
}
joint_trans_subsurfrunoff2$type<-'subsurf-runoff'
ggplot(joint_trans_subsurfrunoff2,aes(x=duration,y=cumval))+
  geom_point(alpha=0.5)+
  facet_wrap(~roi,ncol=4)+geom_smooth()


joint_df2<-joint_trans_snow2
joint_df2<-rbind(joint_df2,joint_trans_precip2)
joint_df2<-rbind(joint_df2,joint_trans_runoff2)
joint_df2<-rbind(joint_df2,joint_trans_subsurfrunoff2)

ggplot(joint_df2,aes(x=duration,y=cumval,colour=type))+
  geom_point(alpha=0.5)+
  facet_wrap(~roi,ncol=4)#+geom_smooth()

ggplot(joint_df2[joint_df2$roi %in% 33:36,],aes(x=duration,y=cumval,colour=type))+
  geom_point(alpha=0.5)+
  facet_wrap(~roi,ncol=4)#+geom_smooth()

outpath <-"/Users/pittore/Documents/workspace/ERA5/Results"
outfile <-"spells_min0.1.csv"
write_csv(joint_df2,paste(outpath,outfile,sep='/'))

