library(ncdf4)
library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
library(gganimate)
library(raster)
#library(rasterVis)
library(maptools)
library(stpp)
library(readr)
library(zoo)

######################################################################################
#         FUNCTIONS
######################################################################################

# return the x and y rage of the entire dataset (in UTM43 proj)
get_xyrange<-function(ls_df) {
  x_range = range(ls_df$E_main)
  y_range = range(ls_df$N_main)
  return (c(x_range,y_range))
}

# return the x and y rage of the landslide centroids 
# in a dataset (in UTM43 proj)
# note: it expects the fields centroid_x,centroid_y
get_xyrange_centroids<-function(ls_df) {
  x_range = range(ls_df$centroid_x)
  y_range = range(ls_df$centroid_y)
  return (c(x_range,y_range))
}

#retun overall time range as desumed from the observation ref dates
get_timerange<-function(ls_df) {
  first_day = min(ls_df$FirstRefDate)
  last_day = max(ls_df$LastRefDate)
  obs_period = range(first_day,last_day)
  return(obs_period)
}

#select only the Behling landslides, group by landslide event and count pixels, ref dates and ref period
get_behling_events<-function(ls_df) {
  lss_behling_dates <- ls_df %>% 
    filter (LSno_Behling>0) %>%
    #select(E_main, N_main,LSno_Behling,FirstRefDate,LastRefDate) %>% 
    group_by(LSno_Behling) %>% 
    summarise(
      size=n(),
      FirstRefDate=first(FirstRefDate),
      LastRefDate=first(LastRefDate),
      centroid_x = mean(E_main), 
      centroid_y = mean(N_main)) %>%
    mutate(refperiod=LastRefDate-FirstRefDate)
 return (lss_behling_dates)
}

# compute the cdf of the landslide count (minimum and maximum), optionally over a given ROI
# the ladnslides can be filtered according to size (area, in terms of number of pixels)
# NOTE: reference projection for N_main and E_main is UTM43
# NOTE: each refdate refer to the number of days elapsed since January 1st 0000 (annum zero) <- CHECK IT
compute_cdfcount_behling<-function(ls_df,x_range=NULL,y_range=NULL,size_range=NULL) {
  
  x_range <- as.numeric(x_range)
  y_range <- as.numeric(y_range)
  
  if (is.null(x_range)) x_range = range(ls_df$E_main)
  if (is.null(y_range)) y_range = range(ls_df$N_main)
  if (is.null(size_range)) size_range = range(0,1e10)
  
  #select only the Behling landslides, group by landslide event and count pixels, ref dates and ref period
  lss_behling_dates <- get_behling_events(lss)
  
  lss_behling_dates_filt <- lss_behling_dates %>%
    filter (centroid_x >= x_range[1] & centroid_x <= x_range[2]) %>%
    filter (centroid_y >= y_range[1] & centroid_y <= y_range[2]) %>%
    filter (size >= size_range[1] & size <= size_range[2])
  
  if (nrow(lss_behling_dates_filt)<1) return(NULL)
  
  #get the overall observation period
  first_day = min(lss_behling_dates_filt$FirstRefDate)
  last_day = max(lss_behling_dates_filt$LastRefDate)
  obs_period = range(first_day,last_day)
  
  #initialize a vector of days from the first to the last observations
  ls_counts_days = data.frame("days"=seq(first_day,last_day),"min_count"=0,"max_count"=0)
  
  #fill up the CDF (cumulative function) of the yearly minimum and maximum count
  for(i in 1:nrow(lss_behling_dates_filt)) {
    start=lss_behling_dates_filt[i,]$FirstRefDate
    end=lss_behling_dates_filt[i,]$LastRefDate
    ls_counts_days[ls_counts_days$days >= end,]$min_count = 
      ls_counts_days[ls_counts_days$days >= end,]$min_count +1
    ls_counts_days[ls_counts_days$days >= start,]$max_count = 
      ls_counts_days[ls_counts_days$days >= start,]$max_count +1
  }
  #add a date
  ls_counts_days$date = as.Date(ls_counts_days$days, origin = "0-01-01")
  return(ls_counts_days)
}

#convert days into dates
days2date<-function(days){
  return(as.Date(days, origin = "0-01-01"))
}

#compute the max and minimal landslides number in the temp_win centered on each day
#temp_win = 60 #temporal wíndow in days
lscnt_twin<-function(ls_cdf,twin){
  ls_sum <- ls_cdf %>% 
    mutate( min_count_twin = lead(min_count,n=twin/2)-lag(min_count,n=twin/2))  %>% 
    mutate( max_count_twin = lead(max_count,n=twin/2)-lag(min_count,n=twin/2))
  return(ls_sum)
}
#test###########################################
#tst<-lscnt_twin(ls_counts[ls_counts$date < '2016-01-01',],90)
#tst
#and plot it
#ggplot(tst,aes(date,min_count_twin))+geom_line(colour='blue')+
#  geom_line(aes(date,max_count_twin),colour='red')+
#  geom_line(aes(date,min_count),colour='grey')+
#  geom_line(aes(date,max_count),colour='green')
################################################

#compute the minimum and maximum number of landslides in 
#a given time range given a landslide count CDF
get_minmax_counts<-function(ls_cdf,trange){
  ls_trange<-range(ls_cdf$days)
  #round up days in trange to the nearest integer
  trange<-as.numeric(round(trange))
  
  if (trange[1] %in% ls_cdf$days) {
    t1<-ls_cdf[ls_cdf$days==trange[1],]
  } else {
    if (trange[1]<ls_trange[2]) {
      t1<-ls_cdf[ls_cdf$days==ls_trange[1],] #use the first time ref
    } else {
        return (c(0,0)) # the time range is incompatible 
      }
  }
  
  if (trange[2] %in% ls_cdf$days) {
    t2<-ls_cdf[ls_cdf$days==trange[2],]
  } else {
    if (trange[2]>=ls_trange[1]) {
      t2<-ls_cdf[ls_cdf$days==ls_trange[2],]# use the last time ref
    } else {
      return (c(0,0)) # the time range is incompatible 
    }
  }
  #print (c(t1,t2))
  maxc<-t2$max_count-t1$min_count
  minc<-t2$min_count-t1$min_count
  return(c(minc,maxc))
}
#get_minmax_counts(ls_cnt,roi_df[2,6:7])

#generate a ROI of size x =size[1],y=size[2] from the range
#if bRandomSize is TRUE than size is interpreted as minsizex,maxsizex,minsizey,maxsizey
get_random_roi<-function(range,size,bRandomSize=FALSE){
  #sample the north-western / topleft corner of the roi
  #the granularity is one meter, since the proj is UTM43
  xseed<-base::sample(range[1]:(range[2]-size[1]),1)
  yseed<-base::sample((range[3]+size[2]):range[4],1)
  #define the roi
  roi<-c(xseed,xseed+size[1],yseed-size[2],yseed)
  return(roi)
}

#generate nx x ny rectangular ROIs covering the range
#each roi is defined by a spatial range (xmin,xmax,ymin,ymax) and area in squared meters
#and a time range (daymin, daymax) from 01.01.1970 
get_regular_roi<-function(srange,trange,nx,ny,nt) {
  dx <- (srange[2]-srange[1])/nx
  dy <- (srange[4]-srange[3])/ny
  dt <- (trange[2]-trange[1])/nt
  roi_df<-data.frame(index=numeric(nx*ny*nt),xmin=numeric(nx*ny*nt),xmax=numeric(nx*ny*nt),
                     ymin=numeric(nx*ny*nt),ymax=numeric(nx*ny*nt),tmin=numeric(nx*ny*nt),
                     tmax=numeric(nx*ny*nt),area=numeric(nx*ny))
  i<-1
  for (iy in 0:(ny-1)) {
   for (ix in 0:(nx-1)) {
     for (it in  0:(nt-1)) {
       roi_df$index[i]<-i
       roi_df$xmin[i]<-srange[1]+ix*dx
       roi_df$xmax[i]<-srange[1]+(ix+1)*dx
       roi_df$ymin[i]<-srange[3]+iy*dy
       roi_df$ymax[i]<-srange[3]+(iy+1)*dy
       roi_df$area[i]<-dx*dy
       roi_df$tmin[i]<-floor(trange[1]+it*dt)
       #print(trange[1]+it*dt)
       #print(floor((trange[1]+it*dt)))
       #print (trange[1]+(it+1)*dt-1)
       #floor (floor((trange[1]+(it+1)*dt-1)))
       roi_df$tmax[i]<-floor(trange[1]+(it+1)*dt-1) #round to the preceding int
       if (it==nt-1) roi_df$tmax[i]<-floor(trange[1]+(it+1)*dt)
       i<-i+1
    }
   }  
  }
  return(roi_df)
}
#get_regular_roi(xyrange,timerange,1,1,4)
#roi_df

#project a xy range from the projection source_proj 
#to the projection target_proj
project_xyrange<-function(input_range, source_proj, target_proj) {
  input_range<-xyrange
  source_proj<-proj_utm43
  target_proj<-proj_wgs84
  roi_llcorners_pro1<-SpatialPoints(cbind(input_range[1],input_range[3]),CRS(source_proj))
  roi_urcorners_pro1<-SpatialPoints(cbind(input_range[2],input_range[4]),CRS(source_proj))
  roi_llcorners_pro2<-spTransform(roi_llcorners_pro1,target_proj)
  roi_urcorners_pro2<-spTransform(roi_urcorners_pro1,target_proj)
  return(unname(c(roi_llcorners_pro2@coords[,1],
           roi_urcorners_pro2@coords[,1],
           roi_llcorners_pro2@coords[,2],
           roi_urcorners_pro2@coords[,2])))
}

#project a set of rois (in a dataframe) from the projection source_proj 
#to the projection target_proj
project_rois<-function(rois,source_proj,target_proj) {
  roi_llcorners_pro1<-SpatialPoints(cbind(rois$xmin,rois$ymin),CRS(source_proj))
  roi_urcorners_pro1<-SpatialPoints(cbind(rois$xmax,rois$ymax),CRS(source_proj))
  roi_llcorners_pro2<-spTransform(roi_llcorners_pro1,target_proj)
  roi_urcorners_pro2<-spTransform(roi_urcorners_pro1,target_proj)
  
  res_df<-data.frame("xmin"= roi_llcorners_pro2@coords[,1],"xmax"=roi_urcorners_pro2@coords[,1],
                     "ymin"= roi_llcorners_pro2@coords[,2],"ymax"= roi_urcorners_pro2@coords[,2])
  return(res_df)
}
######################################################################################

setwd("/Users/pittore/Documents/ACTIVITIES/landslides/CATENA")

#NOTE: the following takes a while ! 
lss <- read_csv("/Users/pittore/Documents/ACTIVITIES/landslides/CATENA/LSSusceptibility_v5.csv")
str(lss)

#get the info on the behling dataset
lss_behling<-get_behling_events(lss)
#and plot it
ggplot(lss_behling,aes(centroid_x,centroid_y,size=size,colour=refperiod))+
  geom_point(alpha=0.5) + scale_color_gradient(low="blue", high="red")

#get the extent of the study area in space (in UTM43) 
#and in time (in days from 01.01.1970)
xyrange <- get_xyrange_centroids(lss_behling) #[xmin,xmax,ymin,ymax]
timerange<-get_timerange(lss_behling) # daymin,daymax

#get a dataframe with a regular tessellation of the range
#including the time range
roi_df<-get_regular_roi(xyrange,timerange,3,3,4)
roi_df

proj_utm43<-"+proj=utm +zone=43 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj_wgs84<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

#convert xyrange into WGS84 and obtain a new set of ROIs
xyrange_wgs84<-project_xyrange(xyrange,proj_utm43,proj_wgs84)
roi_df_wgs84<-get_regular_roi(xyrange_wgs84,timerange,3,3,4)

#caution, the new rois might overlap ! 
#roi_df_wgs84<-project_rois(roi_df,proj_utm43,proj_wgs84)
roi_df_wgs84

#get the min/max CDF of the events count in the defined rois
nrois<-nrow(roi_df)
roi_df$min_count<-rep(0,nrois)
roi_df$max_count<-rep(0,nrois)
for (roi_ind in 1:nrois) {
  print(sprintf("processing roi %d of %d ...",roi_ind,nrois))
  ls_cnt <- compute_cdfcount_behling(lss_behling,x_range=roi_df[roi_ind,2:3],y_range=roi_df[roi_ind,4:5])
  lcntinterval<-get_minmax_counts(ls_cnt,roi_df[roi_ind,6:7])
  roi_df_wgs84$min_count[roi_ind]<-lcntinterval[1]
  roi_df_wgs84$max_count[roi_ind]<-lcntinterval[2]
}
roi_df_wgs84$area<-roi_df$area/1.e6
roi_df_wgs84
sum(roi_df_wgs84$min_count)
sum(roi_df_wgs84$max_count)

write_csv(roi_df_wgs84,"/Users/pittore/Documents/ACTIVITIES/landslides/CATENA/areaset_3x3x4_wgs84.csv")

# add dates
ls_dates_file<-"/Users/pittore/Documents/ACTIVITIES/landslides/CATENA/LSSusceptibility_Dates_v2.csv"
ls_dates<-read_csv(ls_dates_file)
lss['firstrefdate'] <- ls_dates$FirstRefDate
lss['lastrefdate'] <- ls_dates$LastRefDate

str(robdf_ls)

get_dates<-function(e,n)
{
 # e<- robdf_ls$E_mean[1]
 # n<- robdf_ls$N_mean[1]
 res<- which.min((lss$E_main-e)^2+(lss$N_main-n)^2)
 return(res[1])
}

res1<-lapply(seq(1,nrow(robdf_ls)),function(x){get_dates(robdf_ls[i,]$E_mean,robdf_ls[i,]$N_mean)})


#lss2 <- read_csv("/Users/pittore/Documents/ACTIVITIES/landslides/CATENA/Susceptibility/LSSusceptibility_v2.csv")
#str(lss2)


#dataset from Robert Behling
robdf<-read_csv("RobertLSData.csv")

#show all pixels, colormapped according to slope
p0<-ggplot(robdf)+geom_raster(aes(x=E_main,y=N_main,fill=Slope))+
  scale_fill_viridis(begin=1,end=0)
show(p0)

#filter out non-lanslide pixels, group by landslide id
robdf_prc<- robdf %>%
  filter(LSno > 0) %>%
  group_by(LSno)

#cmpute centroid for each landslide
robdf_ls<- robdf_prc %>%
  summarise(
    count = n(),
    E_mean=mean(E_main,na.rm=TRUE),
    N_mean=mean(N_main,na.rm=TRUE),
    area_min=min(LSArea,na.rm=TRUE),
    area_max=max(LSArea,na.rm=TRUE),
    slope_min=min(Slope,na.rm=TRUE),
    slope_max=max(Slope,na.rm=TRUE),
    slope_mean=mean(Slope,na.rm=TRUE)
    )

size(robdf_ls$area_max)
cdf <-ecdf(robdf_ls$area_max)
#data.frame(cdf)
plot(cdf, verticals = TRUE, do.points = TRUE,log='x')

ggplot(NULL, aes(x=robdf_ls$area_max))+geom_step(stat='ecdf')+scale_x_log10()


#check distribution of landslides area
pa<-ggplot(data=robdf_ls)+geom_histogram(aes(x=area_max))+scale_x_log10()+scale_y_log10()
show(pa)

#p2<-ggplot(robdf_ls)+geom_point(aes(x=E_mean,y=N_mean,size=count),alpha=0.1,col='brown') 
p2<-ggplot()+geom_point(data=robdf_ls,mapping=aes(x=E_mean,y=N_mean,size=area_max,col=slope_max),alpha=0.6) 
show(p2)
ggsave("Behling_ls.png",dpi=150)



########## HAVENITH DATA ######################################

#dtaset from Hans HAvenith
havdf<-read_csv("HavenithLSData.csv")

#filter out non-lanslide pixels, group by landslide id
havdf_prc<- havdf %>%
  filter(LSno > 0) %>%
  group_by(LSno)

#cmpute centroid for each landslide
havdf_ls<- havdf_prc %>%
  summarise(
    count = n(),
    E_mean=mean(E_main,na.rm=TRUE),
    N_mean=mean(N_main,na.rm=TRUE)
  )

#ph<-ggplot(havdf_ls)+geom_point(aes(x=E_mean,y=N_mean,size=count),alpha=0.1,col='red') 
ph<-p0+geom_point(data=havdf_ls,mapping=aes(x=E_mean,y=N_mean,size=count),alpha=0.1,col='red')+
  labs(title='Havenith dataset',x='easting',y='northing')
show(ph)
ggsave("havenith_ls.png",dpi=150)

#estimate the bandwith ?
bandwidth.nrd(havdf_ls$E_mean)

#plot the 2d kernel density of the landslides
pj<-ggplot()+geom_point(data=havdf_ls,mapping=aes(x=E_mean,y=N_mean,size=count),alpha=0.5,col='red')+
  geom_density_2d(data=havdf_ls,aes(x=E_mean,y=N_mean),h=c(15000,15000),col='gray',alpha=0.9)+
  labs(title='Havenith dataset (density)',x='easting',y='northing',size='area')
show(pj)
ggsave('havenith_density.png',dpi=150)

pj<-ggplot()+geom_point(data=robdf_ls,mapping=aes(x=E_mean,y=N_mean,size=count),alpha=0.1,col='red')+
  geom_density_2d(data=robdf_ls,aes(x=E_mean,y=N_mean),h=c(15000,15000),col='gray',alpha=0.9)+
  labs(title='Behling dataset (density)',x='easting',y='northing',size='area')
show(pj)
ggsave('behling_density.png',dpi=150)


library(MASS)
#density plot
hav_dens<-kde2d(havdf_ls$E_mean,havdf_ls$N_mean,n=100)

image(hav_dens)
contour(hav_dens, add=T)


size(robdf_ls$count)

#library(R.matlab)
#test<-readMat('/Users/pittore/Documents/ACTIVITIES/landslides/CATENA/Susceptibility/Everything_v190621.mat')



