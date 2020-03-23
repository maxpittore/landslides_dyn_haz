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
library(stpp)
library(readr)
library(raster)
library(stpp)
library(spatstat)
library(beepr)
library(pracma)

setwd("/Users/pittore/Documents/ACTIVITIES/landslides/CATENA/Susceptibility")

#read data. Can take a while !
filename <- 'LSSusceptibility.csv'
T <- read_csv(filename)
dim(T)

#pointmap(as.points(T[,1:2]))
ls_mask <- T %>% select('E_main','N_main','Mask_Behling') %>%
  rename('x'='E_main','y'='N_main','m'='Mask_Behling')
ls_mask_r <- rasterFromXYZ(ls_mask)
projection(ls_mask_r) <- "+proj=utm +zone=43 +datum=WGS84"
ls_mask_r[ls_mask_r < 0.5] <- NA
plot(ls_mask_r)

#pointmap(as.points(T[,1:2]))
ls_behling_df <- T %>% select('E_main','N_main','Suscebtibility_Behling') %>%
                      filter('Mask_Behling'>0.5) %>%
                       rename('x'='E_main','y'='N_main','ls'='Suscebtibility_Behling')
#create a 2D raster
ls_behling_r <- rasterFromXYZ(ls_behling_df)
projection(ls_behling_r) <- "+proj=utm +zone=43 +datum=WGS84"

#apply mask
ls_behling_rm<-mask(ls_behling_r,ls_mask_r)
region = extent(ls_behling_r)
plot(ls_behling_rm)
#save raster to geotiff
writeRaster(ls_behling_rm,'ls_behling_rm.tif')

#aggregate to decrease resolution 
#original resolution: 30m x 30m - target resolution: 1km x 1km
ls_behling_rm_agg<-aggregate(ls_behling_rm,fact= 33.33333, fun=mean)
plot(ls_behling_rm_agg)
#writeRaster(ls_behling_rm_agg,'ls_behling_rm_lowres.tif')

#reproject to wgs84
newcrs<-'+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
ls_behling_r_wgs<-projectRaster(ls_behling_rm_agg,crs=newcrs)
plot(ls_behling_r_wgs)
writeRaster(ls_behling_r_wgs,'ls_behling_rm_lowres_wgs.tif')

region_wgs<-extent(ls_behling_r_wgs)

p0<-ggplot() +geom_raster(data=ls_behling_df, aes(x=x, y=y, fill=ls), alpha=0.8)+
  scale_fill_viridis()
show(p0)

################ precipitation data #############################
nc<-nc_open("CHIRPSstacked.nc")
ncs<-nc_open("CHIRPSstacked_swap.nc")
#print(paste("The file has",nc$nvars,"variables,",nc$ndims,"dimensions and",nc$natts,"NetCDF attributes"))
#print(nc)
#attributes(nc$var)$names
#ncatt_get(nc, attributes(nc$dim)$names[1])

ncs_lat <- ncvar_get( ncs, attributes(ncs$dim)$names[3])
nc_lon <- ncvar_get( nc, attributes(nc$dim)$names[2])
#nc_prcp<- ncvar_get( nc, attributes(nc$var)$names[3])

chirp.b<-brick("CHIRPSstacked.nc")
#transpose to fix problem of x/y swap in the raster brick loading...
chirp.b <- t(chirp.b)
projection(chirp.b)<-'+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
extent(chirp.b)

#raster::animate(har.b,pause=0.05)
res(chirp.b)
crs(chirp.b)
#################################################################
#harmonize extents
common_extent <- extent(72.49999,73.89842,40.09528,41.42152)

ls_behling_r_wgs_crop <- crop(ls_behling_r_wgs,common_extent)
chirp.b_crop <- crop(chirp.b,common_extent)

extent(ls_behling_r_wgs_crop)
extent(chirp.b_crop)
#################################################################
# harmonize_resolution
s<-Sys.time()
chirp.b_crop<-resample(chirp.b_crop,ls_behling_r_wgs_crop)
e<-Sys.time()
e-s

plot(chirp.b_crop[[10]])
writeRaster(chirp.b_crop[[10]],'chirp_b_crop_layer10.tif')

#test<-aggregate(chirp.b_crop,fact=c(1,1,10),fun=sum)
#raster::animate(test)

#################################################################
## moving cumulative map
compute_movsum<-function(brickdata,cum_days)
{
  nl<-nlayers(brickdata)
  mw<-cum_days #number of days
  
  movsum<-function(x)
  {
    start<-x-mw
    end<-x
    sub <-brickdata[[start:end]]
    return (sum(sub))
  }
  #compute for each day the sum of the precipitation in the previous mw (10) days
  #and turn it into a raster brick
  mov_sum<-brick(lapply((mw+1):nl,movsum))
  return (mov_sum)
}

mov_sum_10d<-compute_movsum(chirp.b_crop,10)

plot(mov_sum_10d[[110]])
writeRaster(mov_sum_10d[[110]],'chirp_movsum_10d_layer110.tif')

raster::animate(mov_sum_10d[[100:150]],n=1)

#har_spdf<-as(har.b[[1]], "SpatialPixelsDataFrame")
chirp_df <- as.data.frame(chirp.b_crop[[1]], xy = TRUE) 
chirp_df2 <- as.data.frame(chirp.b, xy = TRUE) 

p<-ggplot() +geom_raster(data=chirp_df, aes(x=x, y=y, fill=layer.1), alpha=0.8)+
  scale_fill_viridis()
show(p)

#gganimate(p)

########### TEST STATIC POISSON PROCESS #############

#generate inhomogeneous distribution of points
#according to the input raster, as intensity of a 
#Poisson Point Process
#return a set of spatial points
PPPoints<-function(rast,coef)
{
  mbox<-bbox(rast)
  mat<-as.matrix(flip(rast,'y'))
  ndens<-as.im(mat*coef)
  Z<-rpoispp(ndens)
  #rescale the points on the raster extent
  x<-Z$x*((mbox[1,2]-mbox[1,1])/ndens$xrange[2])+mbox[1,1]
  y<-(Z$y)*((mbox[2,2]-mbox[2,1])/ndens$yrange[2])+mbox[2,1]
  return(SpatialPoints(cbind(x,y)))
}

# generate stochastic realization after the static susceptibility
#single realization
pts_dens <- PPPoints(ls_behling_r_wgs_crop,1)
plot(ls_behling_r_wgs_crop)
plot(pts_dens,add=TRUE)

#save to csv file
pts_dens_df<-as.data.frame(pts_dens)
write.csv(pts_dens_df, file = "ls_behling_r_wgs_crop_realisation.csv")

########### CHECK  check for consistency with count number ###############
#multiple realizations - check for consistency with count number
pts_list<-sapply(1:100,function(x) {length(PPPoints(ls_behling_r_wgs,1))})
plot(pts_list)
print(mean(pts_list))
print(cellStats(ls_behling_r_wgs,sum,na.rm=TRUE))
########

#init map for plotting
mp<- get_stamenmap(bbox = c(left = region_wgs@xmin, bottom = region_wgs@ymin, right =
                              region_wgs@xmax, top = region_wgs@ymax), zoom = 10, maptype = "terrain") 

plot1 <- mp %>% ggmap() + stat_density_2d(data = as.data.frame(pts_dens),
                                          aes(x = x,
                                              y = y,
                                              fill = stat(level)),
                                          alpha = .04,
                                          bins = 25,
                                          geom = "polygon") + 
  geom_density_2d(data = as.data.frame(pts_dens),
                  aes(x = x,
                      y = y,
                      fill = stat(level)),
                  alpha = .8,
                  bins = 5)
show(plot1)

test_spdf <- as(ls_behling_r_wgs_crop, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
plot01 <- mp %>% ggmap() + geom_tile(data=test_df,aes(x=x,y=y,fill=factor(ls)))
show(plot01)

pts <- PPPoints(ls_behling_r_wgs,1e-4)
pts
plot(pts,add=TRUE)

plot2 <- plot1 + geom_point(data=as.data.frame(pts), aes(x=x,y=y,col='red'),alpha=0.9) #+
# geom_point(data=ad_stations,aes(x=lon,y=lat),shape=24,fill="red") +geom_text_repel(data=ad_stations,aes(label=stat_code)) + labs(title="all events")
show(plot2)

################## TEST DYNAMIC POISSON PROCESS ###################

#base model - only static susceptibility
ntimes<-10
dynamic0<-brick(lapply(1:ntimes,function(x) ls_behling_r_wgs_crop))
plot(dynamic0[[1]])
#image(dynamic0[[1]])

reshape_model<-function(model)
{
  tmp1<-raster::as.array(flip(model,direction='y'))
  model_array<-aperm(tmp1, c(2,1,3),resize=TRUE)
  return(model_array)
}

dynamic0_arr <- reshape_model(dynamic0)
image(dynamic0_arr[,,1])

#average total number of points expected
print(sum(dynamic0_arr,na.rm=TRUE))

ext<-extent(dynamic0)
roi<-matrix(c(ext@xmin,ext@xmax,ext@xmax,ext@xmin,ext@ymin,ext@ymin,ext@ymax,ext@ymax),c(4,2))
#roi<-matrix(c(73.0,73.5,73.5,73.0,40.2,40.2,40.8,40.8),c(4,2))
polymap(roi)
roi.area<-(ext@xmax-ext@xmin)*(ext@ymax-ext@ymin)


##### TEST WITH TIME- CONSTANT LAMBDA
sp_sim<-rpp(lambda=200,s.region<-roi,t.region=c(1,10),discrete.time=TRUE)
str(sp_sim)


##### TEST WITH  TIME- CONSTANT LAMBDA

lambda_test<-dynamic0_arr*1000
lambda_test[is.na(lambda_test)]<-0
lambda_test[is.nan(lambda_test)]<-0

#compute 2D integral of the density over the grid
get_surf<-function(x,y)
{
  x1<-linspace(0,1,149)
  y1<-linspace(0,1,120)
  Z<-lambda_test[,,1]
  interp2(x=x1,y=y1,Z,xp=x,yp=y,method='linear')
}
mean_evts<-integral2(get_surf, 0, 1, 0, 1, reltol = 1e-10)

#sp_sim0<-rpp(lambda=lambda_test,discrete.time=FALSE)
#str(sp_sim0)
#sp_sim0.pts<-as.data.frame(as.points(sp_sim0$xyt[,1:3]))
#plot(ls_behling_r_wgs_crop)

#lambda_test1<-array(10,c(120,149,10))
#sum(lambda_test1[,,1])/length(which(lambda_test1[,,1]>1e-6))#*roi.area*100

simulate_0<-function(dens)
{
  sp_s<-rpp(lambda=dens,discrete.time=FALSE)
  return(length(sp_s$xyt)/3)
}

res<-sapply(1:200,function (x) simulate_0(lambda_test))
beep(4)
mean(res)

############ ONLY PRECIPITATION ######################

dynamic1<-mov_sum_10d[[20:50]]
raster::animate(dynamic1,n=1)
dynamic1_arr <- reshape_model(dynamic1)
image(dynamic1_arr[,,1])

lambda_test1<-dynamic1_arr*10
lambda_test1[is.na(lambda_test1)]<-0
lambda_test1[is.nan(lambda_test1)]<-0

#compute 2D integral of the density over the grid
get_surf<-function(x,y)
{
  x1<-linspace(0,1,149)
  y1<-linspace(0,1,120)
  Z<-lambda_test1[,,30]
  interp2(x=x1,y=y1,Z,xp=x,yp=y,method='linear')
}
mean_evts<-integral2(get_surf, 0, 1, 0, 1, reltol = 1e-10)

sp_sim1<-rpp(lambda=lambda_test1,t.region=c(1,30),discrete.time=FALSE)
str(sp_sim1)
plot(sp_sim1$xyt,cex=0.5)
sp_sim0.pts<-as.data.frame(as.points(sp_sim0$xyt[,1:3]))

############ SUSCEPTIBILITY AND PRECIPITATION ######################

#generate a dynamic susceptibility by adding a multiplicative 
#component based on the moving cumulative precipitation
#weighted by a suitable coefficient

par1<-10
dynamic2<-ls_behling_r_wgs_crop*(1+mov_sum_10d[[20:50]]*par1)
#dynamic1<-ls_behling_r_wgs_crop*(1+mov_sum_10d*par1)
raster::animate(dynamic2,n=1)

dynamic2_arr <- reshape_model(dynamic2)
image(dynamic2_arr[,,10])

lambda_test2<-dynamic2_arr*0.1
lambda_test2[is.na(lambda_test2)]<-0
lambda_test2[is.nan(lambda_test2)]<-0

#compute 2D integral of the density over the grid
get_surf<-function(x,y)
{
  x1<-linspace(0,1,149)
  y1<-linspace(0,1,120)
  Z<-lambda_test2[,,1]
  interp2(x=x1,y=y1,Z,xp=x,yp=y,method='linear')
}
mean_evts<-integral2(get_surf, 0, 1, 0, 1, reltol = 1e-10)

sp_sim2<-rpp(lambda=lambda_test2,s.region=roi,t.region=c(1,30),discrete.time=FALSE)
str(sp_sim2)
plot(sp_sim2$xyt,cex=0.5)

sp_sim2.pts<-as.data.frame(as.points(sp_sim2$xyt[,1:3]))
write.csv(sp_sim2.pts, file = "ls_behling_precip_30_50_realiz.csv")


#############################################################


## analysis
nl <- dim(lambda_test2)[3]
st_sum <- sapply(1:nl,FUN =function(x) sum(lambda_test2[,,x],na.rm = TRUE))
st_mean <- sapply(1:nl,FUN =function(x) mean(lambda_test2[,,x],na.rm = TRUE))
st_max <- sapply(1:nl,FUN =function(x) max(lambda_test2[,,x],na.rm = TRUE))
plot(st_sum,type='l')
plot(st_mean,type='l')
plot(st_max,type='l')

st_mean_df<-as.data.frame(st_mean)*2
st_mean_df['t']<-1:31

ph<-ggplot(sp_sim2.pts) + geom_histogram(aes(t),colour='lightgray')+theme_bw()+
  geom_line(data=st_mean_df,aes(x=t,y=st_mean),colour='blue')+
  xlab('time')+scale_y_continuous('counts',
                                  sec.axis = sec_axis(~ .,name='mean cumulated precipitation'))
show(ph)
ggsave('hist_counts_precip.png',dpi='retina',device='png')

#which(stat>400)
plot(dynamic2[[1]])

###################################################
#save plots to create animations


for (n_layer in 1:nl) 
{
  #n_layer=1
  print(n_layer)
  
  pts_dens <- PPPoints(dynamic2[[n_layer]],1)
  pts_dens_df<-as.data.frame(pts_dens)
  
  plt <- mp %>% ggmap() + 
    stat_density_2d(data = pts_dens_df,
                    aes(x = x,
                        y = y,
                        fill = stat(level)),
                    alpha = .04,
                    bins = 25,
                    geom = "polygon") + labs(fill='proc. density') +
  geom_density_2d(data = pts_dens_df,
                    aes(x = x,
                        y = y,
                        alpha = .8,
                    bins = 5)) + 
  geom_point(data = sp_sim2.pts[sp_sim2.pts$t<n_layer+1,],
                    aes(x = x,
                        y = y),
                        size=1.5,
                        alpha = .8,
                        colour = 'yellow')
  show(plt)
  filename=sprintf('lambda_test2_pts2_%2.2d.png',n_layer)
  ggsave(filename,dpi='retina',device='png')
#show(plt)
}
beep(4)
########################################################

#ext<-extent(dynamic1)
#roi<-matrix(c(ext@xmin,ext@xmax,ext@xmax,ext@xmin,ext@ymin,ext@ymin,ext@ymax,ext@ymax),c(4,2))
#polymap(roi)

sp_sim<-rpp(npoints=1000,lambda=lambda,s.region<-roi,t.region=c(50,150),discrete.time=TRUE)
beep(4)
#stpp::animation(test$xyt,runtime = 10)

#extracts points from spatio-temporal simulation sp_sim
sp_sim.pts<-as.data.frame(as.points(sp_sim$xyt[,1:3]))
hist(sp_sim.pts$t, breaks=seq(50,150,1))


plot3 <- plot1 + geom_point(data=sp_sim.pts, aes(x=x,y=y,col=t),alpha=0.9) #+
# geom_point(data=ad_stations,aes(x=lon,y=lat),shape=24,fill="red") +geom_text_repel(data=ad_stations,aes(label=stat_code)) + labs(title="all events")
show(plot3)

plot(st_sum[50:250],type='l',add=TRUE)

simulate<-function(s,e,npts)
{
  sp_sim<-rpp(npoints=npts,lambda=lambda,s.region<-roi,t.region=c(s,e),discrete.time=TRUE)
  sp_sim.pts<-as.data.frame(as.points(sp_sim$xyt[,1:3]))
  sp_sim.pts<-as.data.frame(as.points(sp_sim$xyt[,1:3]))
  res<-hist(sp_sim.pts$t, breaks=seq(s,e,1))
  #plot(res)
  print('ok')
  return(res$counts)
}

s = 1
e = 250
a<-sapply(1:20,FUN=function(x) simulate(s,e,200))
beep(4)
df<-data.frame(a)
df['ind']<-s:(e-1)
ggplot(data=melt(df,id.vars = 'ind'))+geom_line(aes(x=ind,y=value,col=variable))+
  geom_smooth(aes(x=ind,y=value,col=variable))
