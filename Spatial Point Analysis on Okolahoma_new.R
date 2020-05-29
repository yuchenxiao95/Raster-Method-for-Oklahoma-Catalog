library(spatstat)                              # Spatial statistical analysis by many authors
library(dplyr)  # A fast, consistent tool for working with data frame like objects by Hadley Wickam et al.
library(tidyverse)
library(smerc)                                 #Implements statistical methods for analyzing the counts of areal data, with a 
library(maptools)
library(sf)
library(raster)
library(rasterVis)                        # The rasterVis package complements the raster package, providing a set of methods for enhanced                                       # visualization and interaction. 
library(dismo)
library(spatial)
library(sp)
library(maxstat)
library(rgdal)
library(RColorBrewer)
library(LaplacesDemon)
library(imager)
library(pracma)
library("rgeos")
library(maps)
theme_set(theme_bw())
library("rnaturalearth")
library("rnaturalearthdata")
library(ggspatial)
library(ggplot2)
library(wesanderson)
library(mltools)
library(gganimate)                           # create anamination
library(tools)
library(maps)
library(zoo)
library(lubridate)
library(ggnewscale)                          #use multiple colour and fill scales in ggplot2
# Set the working directory, I always like to do this so I don't lose files and to simplify subsequent read and writes
setwd("C:/Users/12819/Dropbox/My PC (DESKTOP-9VC04LJ)/Documents/Research/Induced Seismicity/InducedSeismicity/Iason Oklahoma/OKL 2006_2018")                           # choose your local working directory


##################################################################################################################
# read in comma delimited data file for SWD
EOR_data_temp = read.csv("homog_may2020_AOIp_AllFormations_2006-2018_EPA_KNS_EOR_corrtp2.csv")
head(EOR_data_temp)
EOR_data <- EOR_data_temp %>%
  dplyr::select(API, year,latitude,longitude,month_1,month_2,month_3,month_4,month_5, month_6,month_7,month_8,month_9,month_10,month_11,month_12) %>%
  rename_at(vars(starts_with("month_")),funs(str_replace(.,"month_",""))) %>%
  gather(key = 'month',value = 'monthly_injection', '1','2','3','4','5','6','7','8','9','10','11','12') %>%
  as.data.frame() %>%
  arrange(API,year)
head(EOR_data)
# change int and chr to numeric
EOR_data$month <- as.numeric(EOR_data$month)
EOR_data$year <- as.numeric(EOR_data$year)
# combine year and month
EOR_data$Date <- with(EOR_data, sprintf("%d-%02d", year, month))
head(EOR_data)

# create a simplified dataframe to capture cmmulative injection volume for each SWD
EOR_cumInj <- EOR_data %>%
  dplyr::select(API,latitude,longitude,`monthly_injection`) %>%
  group_by(API,longitude,latitude) %>%
  summarise(total_inj = sum(monthly_injection))
head(EOR_cumInj)

# create a spatial dataframe for SWD with correct projection
EOR_df <- data.frame(ID=1:nrow(EOR_cumInj), data= EOR_cumInj[,4])
EOR_points <- SpatialPoints(EOR_cumInj[,2:3])
crs(EOR_points) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
new_crs <- CRS("+proj=lcc +lat_1=35.56666666666667 +lat_2=36.76666666666667 +lat_0=35 +lon_0=-98 +x_0=609601.2192024384 +y_0=0 +ellps=clrk66 +datum=NAD27 +to_meter=0.3048006096012192 +no_defs ")
EOR_points <- spTransform(EOR_points,new_crs)
extent <- extent(EOR_points)
EOR_spadf <- SpatialPointsDataFrame(EOR_points,data=EOR_df)
EOR_regdf <- as.data.frame(EOR_spadf)
head(EOR_regdf)
########################################################################################################
# read in comma delimited data file for SWD
SWD_data_temp = read.csv("homog_may2020_AOIp_AllFormations_2006-2018_EPA_KNS_SWD_corrtp2.csv")
SWD_data <- SWD_data_temp %>%
  dplyr::select(API, year,latitude,longitude,month_1,month_2,month_3,month_4,month_5, month_6,month_7,month_8,month_9,month_10,month_11,month_12) %>%
  rename_at(vars(starts_with("month_")),funs(str_replace(.,"month_",""))) %>%
  gather(key = 'month',value = 'monthly_injection', '1','2','3','4','5','6','7','8','9','10','11','12') %>%
  mutate(year=as.numeric(year), month=as.numeric(month),API=as.numeric(API)) %>%
  mutate(Date=make_datetime(SWD_data$year, SWD_data$month)) %>%
  group_by(API) %>%
  mutate(cum_inj=cumsum(monthly_injection)) %>%
  mutate(tol_inj=sum(monthly_injection)) %>%
  as.data.frame() 
head(SWD_data)
str(SWD_data)


# SWD_data$tol_inj<- SWD_data%>%
#   dplyr::select(API,monthly_injection) %>%
#   group_by(API) %>%
#   mutate(cumsum(monthly_injection))

# create a simplified dataframe to capture cmmulative injection volume for each SWD
SWD_cumInj <- SWD_data %>%
  dplyr::select(API,latitude,longitude,`monthly_injection`) %>%
  group_by(API,longitude,latitude) %>%
  summarise(total_inj = sum(monthly_injection))
head(SWD_cumInj)

# create a spatial dataframe for SWD with correct projection
SWD_df <- data.frame(ID=1:nrow(SWD_cumInj), data= SWD_cumInj[,4])
SWD_points <- SpatialPoints(SWD_cumInj[,2:3])
crs(SWD_points) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
new_crs <- CRS("+proj=lcc +lat_1=35.56666666666667 +lat_2=36.76666666666667 +lat_0=35 +lon_0=-98 +x_0=609601.2192024384 +y_0=0 +ellps=clrk66 +datum=NAD27 +to_meter=0.3048006096012192 +no_defs ")
SWD_points_new <- spTransform(SWD_points,new_crs)
SWD_extent <- extent(SWD_points)
SWD_lat_spadf<-SpatialPointsDataFrame(SWD_points,data=SWD_df)
SWD_lat_regdf <- as.data.frame(SWD_lat_spadf)
SWD_spadf <- SpatialPointsDataFrame(SWD_points_new,data=SWD_df)
SWD_regdf <- as.data.frame(SWD_spadf)
head(SWD_lat_regdf)

static <- ggplot(data=world)+
  geom_sf() + 
  geom_sf(data=states) +
  geom_sf(data=counties, fill = NA, color = gray(.5)) +
  geom_text(data=counties_points, aes(x=X,y=Y, label=name),
            color='darkblue', size=2.5)+
  coord_sf(xlim = c(-99.85, -96.10), ylim = c(34.97, 37.04), expand = FALSE) +
  geom_point(data = SWD_data, aes(x = longitude, y = latitude,size=cum_inj, color=cum_inj, alpha=0.2))+
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.35, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering)+
  #facet_wrap("year")+
  scale_colour_gradient(low='blue',high='red')+
  labs(color='Total Injected Volume (bbl)',
       size = "SWD Volume") +
  xlab("Longitude (degree)") + ylab("Latitude (degree)")+
  ggtitle("SWD Map in Okolahoma") 

static

################################################################################################################################
# read in comma delimited data file for earthquake
EOR_data <- EOR_data_temp %>%
  dplyr::select(API, year,latitude,longitude,month_1,month_2,month_3,month_4,month_5, month_6,month_7,month_8,month_9,month_10,month_11,month_12) %>%
  rename_at(vars(starts_with("month_")),funs(str_replace(.,"month_",""))) %>%
  gather(key = 'month',value = 'monthly_injection', '1','2','3','4','5','6','7','8','9','10','11','12') %>%
  as.data.frame() %>%
  arrange(API,year)

Earthquake_data = read.csv("OKL_jun2019_v6b_reasen_2006-2019_Mc2.5_Mmin2.5.csv") 
colnames(Earthquake_data) <- c("Lon","Lat","Year","Month","Day","Magnitude", "Depth_km", "Hour","Minute","Second","Random")
Earthquake_data$Date <- with(Earthquake_data, sprintf("%d-%02d", Year, Month))
head(Earthquake_data)
# Create a data.frame with the same number of rows as the geometries
Earthquake_df <- data.frame(ID=1:nrow(Earthquake_data), Earthquake_data[,c('Magnitude','Depth_km','Year','Month','Date')])
head(Earthquake_df)

Earthquake_points <- SpatialPoints(Earthquake_data[,1:2])
crs(Earthquake_points) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
new_crs <- CRS("+proj=lcc +lat_1=35.56666666666667 +lat_2=36.76666666666667 +lat_0=35 +lon_0=-98 +x_0=609601.2192024384 +y_0=0 +ellps=clrk66 +datum=NAD27 +to_meter=0.3048006096012192 +no_defs ")
Earthquake_points_new <- spTransform(Earthquake_points,new_crs)
Earthquake_extent <- extent(Earthquake_points)

# Combine the SpatialPoints with the data frame
Earthquake_spadf <- SpatialPointsDataFrame(Earthquake_points_new,data=Earthquake_df)
Earthquake_lat_df <- SpatialPointsDataFrame(Earthquake_points,data=Earthquake_df)
Earthquake_regdf <- as.data.frame(Earthquake_spadf)
Earthquake_lat_regdf <- as.data.frame(Earthquake_lat_df)
# Let's visualize the first several rows of our data so we can make sure we successfully loaded it
head(Earthquake_lat_regdf) # show the first several rows of a data table in the console

# Check out the summary statistics for each column
str(Earthquake_regdf) 

#crs(Earthquake_spadf) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
#new_crs <- CRS("+proj=utm +zone=14 +ellps=GRS80 +units=m +no_defs")
#Earthquake_spadf <- spTransform(Earthquake_spadf,new_crs)

# Make a raster with variable nrow and ncol and specify the extent
extent <- c(xmn=1420000, xmx=2600000, ymn=1000, ymx=730000)
Earthquake_raster <- raster(Earthquake_spadf, nrow=20,ncol=20)
extent(Earthquake_raster) <- Earthquake_extent
# Rastreize: transfer values associated with spatial points to raster cells, with function that sum the number of points in cell
Earthquake_rasterize <- rasterize(Earthquake_points, Earthquake_raster, fun='count', na.rm='T')

# Assign the spatial extents 
xlim = c(1420000,2600000); ylim = c(1000,730000)

##############################################################################################################################
str(Earthquake_regdf)
world <- ne_countries(scale = "medium", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
states <- cbind(states, st_coordinates(st_centroid(states)))
states$ID <- toTitleCase(as.character(states$ID))
counties <- st_as_sf(map("county", plot = FALSE, fill = TRUE))
counties <- subset(counties, grepl("oklahoma", counties$ID))
counties$name <- as.character(counties$ID) %>% substring(first=10)
counties_points <- st_centroid(counties)
counties_points <- cbind(counties, st_coordinates(st_centroid(counties$geom)))
counties$area <- as.numeric(st_area(counties))

head(counties)

static <- ggplot(data=world)+
  geom_sf() + 
  geom_sf(data=states) +
  geom_sf(data=counties, fill = NA, color = gray(.5)) +
  geom_text(data=counties_points, aes(x=X,y=Y, label=name),
            color='darkblue', size=2.5)+
  coord_sf(xlim = c(-99.85, -96.10), ylim = c(34.97, 37.04), expand = FALSE) +
  geom_point(data = Earthquake_lat_regdf, aes(x = Lon, y = Lat, size = Magnitude,color = Magnitude), alpha=0.2)+
  scale_colour_gradient(low = "black", high = "red")+
  #geom_point(data = SWD_data, aes(x = longitude, y = latitude, size=cum_inj,color=cum_inj),shape=23,alpha=0.2)+
  #scale_colour_gradient(low='light blue',high='purple') +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.35, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  #facet_wrap("Year")+
  xlab("Longitude (degree)") + ylab("Latitude (degree)")+
  ggtitle("Earthquake Map in Okolahoma") 


# plot earthquake 
static <- ggplot(data = Earthquake_lat_regdf, aes(x = Lon, y = Lat, color = Magnitude, size = Magnitude))+
  geom_point(alpha = 0.2)+
  #facet_wrap("Year")+
  scale_colour_gradient(low = "black", high = "red")+
  xlab("Longitude (degree)") + ylab("Latitude (degree)")+
  ggtitle("Earthquake Map in Okolahoma") 

static
ggsave('Earthquake and SWD in Oklahoma.png',static,device='png',width=40, height=40,units='cm',path = "C:/Users/12819/Dropbox/My PC (DESKTOP-9VC04LJ)/Documents/Research/Induced Seismicity/InducedSeismicity/Iason Oklahoma/OKL 2006_2018/Plots", dpi=600)


anim<- static + transition_time(Year) + 
  enter_fade() +
  exit_disappear() +
  labs(title = "Year: {frame_time}")
anim

Earthquake_anim<- animate(anim, width=20, height=20,units='cm', device='png')
Earthquake_anim

anim_save(filename = "Earthquake_animation.gif", Earthquake_anim,path = "C:/Users/12819/Dropbox/My PC (DESKTOP-9VC04LJ)/Documents/Research/Induced Seismicity/InducedSeismicity/Iason Oklahoma/OKL 2006_2018/Plots",dpi=1000)

##########################################################################################################################
# Begin spatial point pattern analysis on earthquakes Mc 2.5
# set owin
est_win = owin(xlim,ylim)
#creates an object that representing a point pattern dataset in the 2D plane for event and Pore Pressure data
pp_event <- ppp(Earthquake_regdf$Lon, Earthquake_regdf$Lat, window = est_win, unitname=c("meters","meters"))
# Notice that we have just pased the locations (x,y) and the extents of the box aligned with the x and y coordinates and units

# Let's check for and remove duplicates
pp_unique_event = unique(pp_event,warn=TRUE)
pp_unique_event
summary(pp_unique_event)
plot(pp_unique_event, cex = 0.5,main = "Events and Rectangle Hull", asp=1.0)

####################################################################################################
#NEAREST NEIGHBOR AND THE G FUNCTION
par(mfrow=c(1,1)) 
# Calculate the distribution of nearest neighbours
# nndist determines the nearest neighbor for a set of observations within a certain radius
dist_nearest_neighbor_event <- nndist.ppp(pp_unique_event)
hist(dist_nearest_neighbor_event,main="Distance to Nearest Event",xlab="Distance (metre)",nclass = 20) # histogram
plot(ecdf(dist_nearest_neighbor_event),main="Distance to Nearest Event",xlab="Distance (ft)",ylab="Cumulative Probability",col="red",xlim = c(0,10000))

#Nearest neighbour distance function G
#In exploratory analysis, the estimate of G is a useful statistic summarising one aspect of the clustering of points
#For inference purposes, the estimate of G is usually compared to the true value of G for a complete random (possion) point process
#rs: the reduced sample or border correction estimator of G
#km: the spatial Kaplan-Meier estimator of G
#han: the Hanisch correction estimator of G
nearest_neighbour_event <- Gest(pp_unique_event,correction=c("rs", "km", "han"), domain=NULL)
plot(nearest_neighbour_event, main="Nearest Neighbor Distance Function G",xlim = c(0,5000))
#### plot alstats by hand
par(mfrow=c(2,2))
Gmult_e <- envelope(pp_unique_event, fun=Gest, nsim=500, correction='border', nlarge = 3000)
plot(Gmult_e,main="G Functon")
Fmult_e <- envelope(pp_unique_event, fun=Fest, nsim=500, correction='border', nlarge = 3000)
plot(Fmult_e,main="F Function", xlim=c(0,10000))
legend("topleft",legend=a$meaning, lty=a$lty, col=a$col)
Lmult_e <- envelope(pp_unique_event, fun=Lest, nsim=500, correction='border', nlarge = 3000)
plot(Lmult_e,main="L Function", xlim=c(0,10000))
legend("topleft",legend=a$meaning, lty=a$lty, col=a$col)
Jmult_e <- envelope(pp_unique_event, fun=Jest, nsim=500, correction='rs', nlarge = 3000)
plot(Jmult_e,main="J Function", xlim=c(0,10000),ylim = c(0,1))


#### Anisotropy and the nearest-neighbor orientation
par(mfrow=c(1,1))
rose<- rose(nnorient(pp_unique_event), col='grey',start="N", clockwise=TRUE, main= "Rose Diagram for Earthquakes")

#### Morisita Index plot 
miplot(pp_unique_event, xlim=c(0,600000),ylim=c(0,50))

#### Fry plot
#plot(frypoints(pp_unique_event))

###############################################################
par(mfrow=c(1,1))
# Let's start by checking the regular Ripley's K on the events
# Estimates Ripley's reduced second moment function K(r) from a point pattern in a window of arbitrary shape
K = Kest(pp_unique_event, correction="isotropic", nlarge=3000, domain=NULL, var.approx=TRUE, ratio=FALSE)
d <- as.data.frame(K)
a <- plot(K, legend=FALSE)
legend("topleft",legend=a$meaning, lty=a$lty, col=a$col)
Kmult_e <- envelope(pp_unique_event, fun=Kest, nsim=100, correction='isotropic', nlarge = 3000)
plot(Kmult_e,main="Event Ripley K with Confidence Interval")
###Estimating the pair correlation function
###It contains contributions only from interpoint distance less than or equal to r
### There are two ways to compute the pcf, one is to directly compute the pcf, another is to compute pcf from K function
#(1)

#################################################################
#### PCF for Earthquake
par(mfrow=c(2,1)) 
gmult_e <- envelope(pp_unique_event, fun=dg, nsim=500)
plot(gmult_e,main="pcf")

gk <- pcf(K)
plot(gk, main = "Earthquake Pair Correlation Function with Correction", xlim = c(0, 3000))

#####################################

### Use allstats to eficiently computs all of the F,G,J, anf K functions
#plot(allstats(pp_unique_event), correction = c('best'),main = "Function Summary for All Events")

#Calculate an estimate of the L-functio (Besag's trnasformation of Ripley's K-function) for a spatial point pattern
## Since the expected K function grows by the squared distance, both the expected ad observed values of the K function become large as 
#### d increases. As a result, small differences can be difficult to see between expected and observed values in the data. 
#### The L-function as a expected values of 0 for a random point pattern
### K(r) is proportional to r^2 in the plannar case but L(r) is always proportional to r. Thus in practice, when a pattern is assessed for complete 
###spatial randomness, the graph of its L-function is compared to a line, whereas that of the K-function is compared to a parabolic curve. 
###The second-order behavior of a point process can be visualised and interpreted more easily based on the L-function than based on the K-function.
### Statistical experience shows that the fluctuations of estimated K-functions increase with increasing r. The root transformation stabilises these
###fluctuations (both means and variances) and can even make them independent of distance r. 
par(mfrow = c(1,1))
L = Lest(pp_unique_event, correction="isotropic", nlarge=3000, domain=NULL, var.approx=T, ratio=FALSE)
Lmult_e <- envelope(pp_unique_event, fun=Lest, nsim=500, correction='isotropic', nlarge = 3000)
L
plot(Lmult_e, .-r ~ r,main="Earthquake L-Function (Besag's transformation of Ripley's K-function) with Confidence Interval")

#####################################################################################################################
# Spatial Analysis in Injection Only

# Plots the points and calculate and plot the convex hull (clock wise)
par(mfrow=c(1,1)) 


# ppp creates an object of class "ppp" representing a point pattern dataset in the two-dimensional plane
pp_inj <- ppp(SWD_regdf$longitude, SWD_regdf$latitude,window = est_win, unitname=c("meters","meters"))
# Notice that we have just pased the locations (x,y) and the extents of the box aligned with the x and y coordinates and units

# Let's see what we have made
pp_inj

# Let's check for and remove duplicates
pp_unique_inj = unique(pp_inj,warn=TRUE)
pp_unique_inj
summary(pp_unique_inj)
clark = clarkevans(pp_unique_inj, correction=c("guard"), clipregion=est_win)
clark
##########################
###ploting all summary function of injections
#plot(allstats(pp_unique_inj), main = "Function Summary for All Injections")

par(mfrow=c(2,2))
Gmult_i <- envelope(pp_unique_inj, fun=Gest, nsim=1000, correction='best', nlarge = 3000)
plot(Gmult_i,main="G Functon")
Fmult_i <- envelope(pp_unique_inj, fun=Fest, nsim=1000, correction='best', nlarge = 3000)
plot(Fmult_i,main="F Function",xlim=c(0,10000))
Lmult_i <- envelope(pp_unique_inj, fun=Lest, nsim=1000, correction='isotropic', nlarge = 3000)
plot(Lmult_i,main="L Function", xlim=c(0,10000))
Jmult_i <- envelope(pp_unique_inj, fun=Jest, nsim=1000, correction='best', nlarge = 3000)
plot(Jmult_i,main="J Function",ylim = c(0,1), xlim=c(0,10000))


#########################################

par(mfrow=c(2,1)) 
# Calculate the distribution of nearest neighbours
dist_nearest_neighbor_inj <- nndist(pp_unique_inj)
hist(dist_nearest_neighbor_inj,main="Distance to Nearest Inj",xlab="Distance (ft)",nclass = 20) # histogram
plot(ecdf(dist_nearest_neighbor_inj),main="Distance to Nearest Injection",xlab="Distance (ft)",ylab="Cumulative Probability",col="blue",xlim = c(0,60000))

# Nearest neighbour distance function G
nearest_neighbour_inj <- Gest(pp_unique_inj,correction=c("rs", "km", "han"), domain=NULL)
plot(nearest_neighbour_inj, main="Nearest Neighbor Distance Function G for Injection",xlim = c(0,40000))

######################################

# Let's start by checking the regular Ripley's K on the injections
par(mfrow=c(1,1))
K_inj = Kest(pp_unique_inj, correction="isotropic", nlarge=3000, domain=NULL, var.approx=FALSE, ratio=FALSE)
Kmult_i <- envelope(pp_unique_inj, fun=Kest, nsim=1000, correction='isotropic', nlarge = 3000)

plot(Kmult_i,main="Injection Ripley K Function with Confidence Interval")
############################
###Estimating the pair correlation function
###It contains contributions only from interpoint distance less than or equal to r
### There are two ways to compute the pcf, one is to directly compute the pcf, another is to compute pcf from K function
#(1)
par(mfrow=c(2,1))
g_inj <- pcf(pp_unique_inj)
plot(g_inj, main = "Injection Pair Correlation Function without Correction for Small r", xlim = c(0,5000))
gk <- pcf(K_inj, spar = 0.5)
plot(gk, main = "Injection Pair Correction Function from K function", xlim= c(0,5000))
################################

#the L-function for Injections
par(mfrow = c(1,1))
L = Lest(pp_unique_inj, correction="isotropic", nlarge=3000, domain=NULL, var.approx=FALSE, ratio=FALSE)
Lmult_i <- envelope(pp_unique_inj, fun=Lest, nsim=1000, correction='isotropic', nlarge = 3000)
plot(Lmult_i, .-r ~ r,main="SWD L-Function with Confidence Interval")


#### Spatial Analysis on Events and Injections

SWD_temp = SWD_regdf %>%
  dplyr::select(longitude,latitude) %>%
  rename(Lon=longitude,Lat=latitude) %>%
  mutate(TYPE =1)
  
Earthquake_temp = Earthquake_regdf %>%
  dplyr::select(Lon,Lat) %>%
  mutate(TYPE =0)

mydata<- rbind(SWD_temp, Earthquake_temp)

# Plots the points and calculate and plot the convex hull (clock wise)
plot(mydata$Lon, mydata$Lat, cex = 0.5,xlab = 'X (ft)',ylab ='Y(ft)',xlim = xlim, ylim = ylim, main='All and Convex Hull',asp=1.0,col = rgb(0,0,0,alpha=1.0))
hpts = ripras(mydata$SURFX, mydata$SURFY, shape="convex")
plot(hpts,col= rgb(0,0,0,alpha=0.2),add=TRUE)

#plot(mydata_inj$SURFX, mydata_inj$SURFY, col = rgb(1,0,0,alpha=1.0),add=TRUE)
plot(hpts_event,col = rgb(1,0,0,alpha=0.2),add=TRUE)

#Points is a generic function to draw a sequence of points at the specified coordinates
points(mydata_event$SURFX, mydata_event$SURFY,col = rgb(1,0,0,alpha=0.2))
plot(hpts_inj,col = rgb(0,0,1,alpha=0.2),add=TRUE)
points(mydata_inj$SURFX, mydata_inj$SURFY,col = rgb(0,0,1,alpha=0.2))

##Create a point pattern dataset for both events and injection
pp_unmarked <- ppp(mydata$Lon,mydata$Lat,window = est_win, unitname=c("meters","meters"))
# Notice that we have just pased the locations (x,y) and the extents of the box aligned with the x and y coordinates and units

# Let's see what we have made
pp_unmarked

# Let's check for and remove duplicates
pp_unique_unmarked = unique(pp_unmarked,warn=TRUE)
pp_unique_unmarked

# Let's start by checking the regular Ripley's K of all the available sample data; Mixing both events and injection
K = Kest(pp_unique_unmarked, correction="isotropic", nlarge=3000, domain=NULL, var.approx=FALSE, ratio=FALSE)
K
plot(K)
# We use the "isotropic" boundary correction as it is best for rectangular windows (recal our study area is a simple 4000m x 4000m box) 
# We can observe that our well spatial point pattern is ""identifical"" to a random (Poisson) process over all isotropic spatial 
# length scales.  This is a reasonable result as the orginal synthetic data set is based on random sampling.


# Let's perform Cross-K function
marker <- mydata$TYPE

# Now we are ready to make a point set, with locations, bounding box and units
pp <- ppp(mydata$Lon,mydata$Lat, window = est_win, marks = marker, unitname=c("metre","metres"))
# We can also add the marks, we convert our 1 / 2 values as doubles to factors to be recognized as levels 
# pp creates a vector of length sum(y) of zeros with a one at the end of each uncensored time interval for use with 'ehr'
## the function "factor" is used to encode a vector as a factor
####use the binary operator "%mark%" to add marks
pp <- pp %mark% factor(marker)

####################################################################################

# Let's confirm that we have a planar point pattern with 109 points, multitype with levels 1 and 2
pp
# This can be confirmed in the output below

# Let's check for and remove duplicates
pp_unique = unique(pp,warn=TRUE)
par(mfrow=c(1,1))
plot(pp_unique, main="Convex Hull for Both Events and Injection", asp = 1.0)

# Now we can calculate the cross K function for the high values (2) relative to the low values (1)
#Kcross: for a multitype point pattern, estimate the multitype K function which counts the expected number of points of type J within in given 
#distance of a point of type I
K <- Kcross(pp_unique, "0", "1",correction = c("isotropic"))
Kmultevent <- envelope(pp_unique, fun=Kcross, i ="0", j="1", nsim=500, correction='isotropic', ratio = TRUE)
# We can use the built in plot feature to observe the result
par(mfrow=c(1,1))
plot(Kmultevent,xlim=c(0,120000),ylim = c(0,4e10), main = "Cross K Function with Confidence Interval Center on Events")
#######################################################################################

Kinj <- Kcross(pp_unique, "1", "0",correction = c("isotropic"))
Kmultinj <- envelope(pp_unique, fun=Kcross,i="1",j="0", nsim=500, correction='isotropic', ratio = TRUE)

# We can use the built in plot feature to observe the result
par(mfrow=c(1,1))
plot(Kmultinj,xlim=c(0,120000),ylim = c(0,4e10), main = "Cross K Function with Confidence Interval Cente on Injection")
# Over all scales we can see a lower density of high values relative to low values for porosity than predicted by
# random.  This suggestions some type of repulsion, or "anti-cross clustering".  This makes sense as the original
# data came from a simulation with spatial correlation (clusters or highs and lows).  The lows and high are repulsing
# eachother by design.  This occurs consistently acorss all scales observed.

#Cross L function centers on events 
Lmultevent <- envelope(pp_unique, fun = Lcross, i="0",j="1", nsim=500, correction="isotropic",ratio=T)
par(mfrow=c(1,1))
plot(Lmultevent,xlim=c(0,120000),ylim = c(0,150000), main = "Cross L Function with Confidence Interval Center on Events")

#Cross L function centers on injections
Lmultinj <- envelope(pp_unique, fun = Lcross, i="1",j="0", nsim=500, correction="isotropic",ratio=T)
par(mfrow=c(1,1))
plot(Lmultinj, xlim=c(0,120000),ylim = c(0,150000), main = "Cross L Function with Confidence Interval Center on Injections")

#Cross PCF function centers on events
p <- pcfcross(pp_unique, i = "0",j = "1", correction = "isotropic", stoyan = 0.1, main="Cross PCF")
plot(p, main="Cross PCF")

par(oma=c(1,7,0,0),mfrow=c(2,1))
crosspcfevent <- envelope(pp_unique, fun = pcfcross, i="0",j="1", nsim=500, correction="isotropic",ratio=T)
plot(crosspcfevent,main = "Cross PCF Function on Earthquakes",xlim=c(0,150000),ylim = c(0,3),legendpos = "right")
#Cross PCF function centers on injections
crosspcfinj <- envelope(pp_unique, fun = pcfcross, i="1",j="0", nsim=500, correction="isotropic",ratio=T)
plot(crosspcfinj,xlim=c(0,150000),ylim = c(0,3), main = "Cross PCF Function on SWDs",legendpos = "right")



############################################################################################################
## Spatial Point Pattern Analysis on EOR and Earthquakes
par(mfrow=c(1,1)) 
plot(EOR_regdf$longitude, EOR_regdf$latitude, cex = 0.5,xlab = 'X (meter)',ylab='Y(meter)', main='Injection and Convex Hull',asp=1.0,col='red')
extent = extent(EOR_points)
plot(hpts_inj,add=TRUE)

# Assign the spatial extents for eor
xlim = c(1420000,2600000); ylim = c(1000,730000)
est_win = owin(xlim, ylim)

# ppp creates an object of class "ppp" representing a point pattern dataset in the two-dimensional plane
pp_eor <- ppp(EOR_regdf$longitude, EOR_regdf$latitude,window = est_win, unitname=c("meters","meters"))
# Notice that we have just pased the locations (x,y) and the extents of the box aligned with the x and y coordinates and units

# Let's see what we have made
pp_eor

# Let's check for and remove duplicates
pp_unique_eor = unique(pp_eor,warn=TRUE)
pp_unique_eor
summary(pp_unique_eor)

##########################
###ploting all summary function of injections
plot(allstats(pp_unique_eor), main = "Function Summary for All EOR")


par(mfrow=c(2,2))
Gmult_e <- envelope(pp_unique_eor, fun=Gest, nsim=1000, correction='best', nlarge = 3000)
plot(Gmult_e,main="G Functon" )
Fmult_e <- envelope(pp_unique_eor, fun=Fest, nsim=1000, correction='best', nlarge = 3000)
plot(Fmult_e,main="F Function",xlim=c(0,10000))
Lmult_e <- envelope(pp_unique_eor, fun=Lest, nsim=1000, correction='isotropic', nlarge = 3000)
plot(Lmult_e,main="L Function",xlim=c(0,10000))
Jmult_e <- envelope(pp_unique_eor, fun=Jest, nsim=1000, correction='best', nlarge = 3000)
plot(Jmult_e,main="J Function",ylim = c(0,1), xlim=c(0,20000))

#########################################
# Let's start by checking the regular Ripley's K on the eor
K_eor = Kest(pp_unique_eor, correction="isotropic", nlarge=3000, domain=NULL, var.approx=FALSE, ratio=FALSE)

###Estimating the pair correlation function
###It contains contributions only from interpoint distance less than or equal to r
### There are two ways to compute the pcf, one is to directly compute the pcf, another is to compute pcf from K function
#(1)
par(mfrow=c(2,1))
g_eor <- pcf(pp_unique_eor)
plot(g_eor, main = "EOR Pair Correlation Function without Correction for Small r", xlim = c(0,5000))
gk <- pcf(K_eor, spar = 0.5)
plot(gk, main = "EOR Pair Correction Function from K function", xlim= c(0,15000))
################################

#the L-function for EOR
par(mfrow = c(1,1))
L = Lest(pp_unique_eor, correction="isotropic", nlarge=3000, domain=NULL, var.approx=FALSE, ratio=FALSE)
Lmult_e <- envelope(pp_unique_eor, fun=Lest, nsim=1000, correction='isotropic', nlarge = 3000)
plot(Lmult_e, .-r ~ r,main="EOR L-Function with Confidence Interval")


### Spatial Analysis on Events and EOR

EOR_temp = EOR_regdf %>%
  dplyr::select(longitude,latitude) %>%
  rename(Lon=longitude,Lat=latitude) %>%
  mutate(TYPE =1)

Earthquake_temp = Earthquake_regdf %>%
  dplyr::select(Lon,Lat) %>%
  mutate(TYPE =0)

mydata<- rbind(EOR_temp, Earthquake_temp)

# Plots the points and calculate and plot the convex hull (clock wise)
plot(mydata$Lon, mydata$Lat, cex = 0.5,xlab = 'X (ft)',ylab ='Y(ft)',xlim = xlim, ylim = ylim, main='All and Convex Hull',asp=1.0,col = rgb(0,0,0,alpha=1.0))

##Create a point pattern dataset for both events and injection
pp_unmarked <- ppp(mydata$Lon,mydata$Lat,window = est_win, unitname=c("meters","meters"))
# Notice that we have just pased the locations (x,y) and the extents of the box aligned with the x and y coordinates and units

# Let's see what we have made
pp_unmarked
# Let's check for and remove duplicates
pp_unique_unmarked = unique(pp_unmarked,warn=TRUE)
pp_unique_unmarked

# Let's start by checking the regular Ripley's K of all the available sample data; Mixing both events and injection
K = Kest(pp_unique_unmarked, correction="isotropic", nlarge=3000, domain=NULL, var.approx=FALSE, ratio=FALSE)
plot(K)
# We use the "isotropic" boundary correction as it is best for rectangular windows (recal our study area is a simple 4000m x 4000m box) 
# We can observe that our well spatial point pattern is ""identifical"" to a random (Poisson) process over all isotropic spatial 
# length scales.  This is a reasonable result as the orginal synthetic data set is based on random sampling.

# Let's perform Cross-K function
marker <- mydata$TYPE

# Now we are ready to make a point set, with locations, bounding box and units
pp <- ppp(mydata$Lon,mydata$Lat, window = est_win, marks = marker, unitname=c("metre","metres"))
# We can also add the marks, we convert our 1 / 2 values as doubles to factors to be recognized as levels 
# pp creates a vector of length sum(y) of zeros with a one at the end of each uncensored time interval for use with 'ehr'
## the function "factor" is used to encode a vector as a factor
####use the binary operator "%mark%" to add marks
pp <- pp %mark% factor(marker)
####################################################################################

# Let's confirm that we have a planar point pattern with 109 points, multitype with levels 1 and 2
pp
# This can be confirmed in the output below

# Let's check for and remove duplicates
pp_unique = unique(pp,warn=TRUE)
par(mfrow=c(1,1))
plot(pp_unique, main="Convex Hull for Both Events and Injection", asp = 1.0)

# Now we can calculate the cross K function for the high values (2) relative to the low values (1)
#Kcross: for a multitype point pattern, estimate the multitype K function which counts the expected number of points of type J within in given 
#distance of a point of type I
K <- Kcross(pp_unique, "0", "1",correction = c("isotropic"))
Kmultevent <- envelope(pp_unique, fun=Kcross, i ="0", j="1", nsim=500, correction='isotropic', ratio = TRUE)
# We can use the built in plot feature to observe the result
par(mfrow=c(1,1))
plot(Kmultevent,xlim=c(0,120000),ylim = c(0,4e10), main = "Cross K Function Center on Events with EOR")
#######################################################################################

Keor <- Kcross(pp_unique, "1", "0",correction = c("isotropic"))
Kmulteor <- envelope(pp_unique, fun=Kcross,i="1",j="0", nsim=500, correction='isotropic', ratio = TRUE)

# We can use the built in plot feature to observe the result
par(mfrow=c(1,1))
plot(Kmulteor,xlim=c(0,120000),ylim = c(0,4e10), main = "Cross K Function Center on EOR")
# Over all scales we can see a lower density of high values relative to low values for porosity than predicted by
# random.  This suggestions some type of repulsion, or "anti-cross clustering".  This makes sense as the original
# data came from a simulation with spatial correlation (clusters or highs and lows).  The lows and high are repulsing
# eachother by design.  This occurs consistently acorss all scales observed.

#Cross L function centers on events 
Lmultevent <- envelope(pp_unique, fun = Lcross, i="0",j="1", nsim=500, correction="isotropic",ratio=T)
par(mfrow=c(1,1))
plot(Lmultevent,xlim=c(0,120000),ylim = c(0,150000), main = "Cross L Function with Confidence Interval Center on Events")

#Cross L function centers on injections
Lmulteor <- envelope(pp_unique, fun = Lcross, i="1",j="0", nsim=500, correction="isotropic",ratio=T)
par(mfrow=c(1,1))
plot(Lmulteor, xlim=c(0,120000),ylim = c(0,150000), main = "Cross L Function with Confidence Interval Center on Injections")

#Cross PCF function centers on events
p_eor <- pcfcross(pp_unique, i = "0",j = "1", correction = "isotropic", stoyan = 0.1, main="Cross PCF")
plot(p_eor, main="Cross PCF")

par(oma=c(1,7,0,0),mfrow=c(2,1))
crosspcfevent <- envelope(pp_unique, fun = pcfcross, i="0",j="1", nsim=500, correction="isotropic",ratio=T)
plot(crosspcfevent,main = "Cross PCF Function on Earthquakes with EOR",xlim=c(0,150000),ylim = c(0,3),legendpos = "right")
#Cross PCF function centers on injections
crosspcfinj <- envelope(pp_unique, fun = pcfcross, i="1",j="0", nsim=500, correction="isotropic",ratio=T)
plot(crosspcfinj,xlim=c(0,150000),ylim = c(0,3), main = "Cross PCF Function on ERO ",legendpos = "right")
#########################################################################################################################################

#Cross J function on event

crossJevn <- envelope(pp_unique, fun = Jcross, i ="0", j="1", nsim=500, correction= "rs")
plot(crossJevn, main = "Cross J Function with Confidence Interval Center on Earthquake")




# Now we can now reverse and calculate the cross K function for the low values (1) relative to the high values (1)
K <- Kcross(pp, "1", "0",correction = "isotropic")

# Again, we can use the built in plot feature to observe the result
plot(K)

###scan.stat: Spatial scan statistic
#calculates the spatial scan statistic for a zone(a set of spatial regions). 


##### event and injection, K and L functions 
par(mfrow=c(2,2))

#plot(Kmultevent,xlim=c(0,120000),ylim = c(0,4e10), main = "Cross K Function on Earthquakes")

plot(Kmultinj,xlim=c(0,120000),ylim = c(0,4e10), main = "Cross K Function on SWDs")
plot(Kmulteor,xlim=c(0,120000),ylim = c(0,4e10), main = "Cross K Function on EOR")

plot(Lmulteor,.-r ~ r, main = "Cross L Function on EORs")
plot(Lmultinj,.-r ~ r, main = "Cross L Function on SWDs")




par(mfrow=c(2,2))

plot(Kmult_e, main=" Ripley's K Function on Earthquakes")

plot(Kmult_i,main="Ripley's K Function on SWDs")

plot(Lmult_e, .-r ~ r,main="L-Function on Earthquakes")

plot(Lmult_i, .-r ~ r,main="L-Function on SWDs")




par(mfrow=c(1,3))


gk_e <- pcf(K, spar = 0.5)
plot(gk_e, main = "PCF on Earthquakes with Correction", xlim = c(0, 30000))
gk_inj <- pcf(K_inj, spar = 0.5)
plot(gk_inj, main = "PCF on SWDs with Correction ", xlim= c(0,5000))
gk_eor <- pcf(K_eor, spar = 0.5)
plot(gk_inj, main = "PCF on EORs with Correction ", xlim= c(0,5000))



