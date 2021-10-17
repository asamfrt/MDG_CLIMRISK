#rm(ls=list())
#setwd("D:/R - HoMe/MDG_CLIMRISK/")

#load required libraries
library(sp)
library(rgdal)
library(raster)
library(tidyverse)

extrafont::loadfonts(device = "win")
theme_sleek<-function(base_size = 11, base_family = "Helvetica") {
  half_line <- base_size/2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "grey30"),
      strip.text.y = element_text(colour = "grey30"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      legend.title = element_text(colour = "grey30", size = rel(0.9)),
      panel.border = element_rect(fill = NA, colour = "grey70", size = 0.5),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), colour = "grey30"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(colour = "grey30", size = rel(1)),
      plot.subtitle = element_text(colour = "grey30", size = rel(.85))
    )
}

mn_pts<-readOGR("./shps/mangrovePoints.shp")
mn_pts<-spTransform(mn_pts, CRS="+proj=moll")
dim(mn_pts)

mn_layers<-stack(list.files("./raster", pattern = ".tif$", full.names = TRUE))
plot(mn_layers[[1]]);plot(mn_pts, add = TRUE, col = "red")

mn_layers<-projectRaster(mn_layers, crs = "+proj=moll", res = 24000)
res(mn_layers)

mdg_extData<-extract(mn_layers, mn_pts,
                     method='simple', 
                     buffer=24000, 
                     sp = TRUE,
                     dp = TRUE,
                     na.rm = TRUE,
                     fun=mean)
mdg_extData<-as.data.frame(mdg_extData)
dplyr::glimpse(mdg_extData)


mdg_data<-mdg_extData #create a copy of the data
names(mdg_data)
coordinates(mdg_data) <- ~coords.x1+coords.x2
r1 <- rasterize(mdg_data, mn_layers[[1]], 'INS_HIST2000_MDG', fun=mean)
r2 <- rasterize(mdg_data, mn_layers[[2]], 'INS_SSP52050_MDG', fun=mean)
r3 <- rasterize(mdg_data, mn_layers[[3]], 'INS_SSP52100_MDG', fun=mean)

r4 <- rasterize(mdg_data, mn_layers[[4]], 'VEL_HIST2000_MDG', fun=mean)
r5 <- rasterize(mdg_data, mn_layers[[5]], 'VEL_SSP52050_MDG', fun=mean)
r6 <- rasterize(mdg_data, mn_layers[[6]], 'VEL_SSP52100_MDG', fun=mean)

r<-stack(r1, r2, r3, r4, r5, r6)
crs(r)<-"+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
rr<-projectRaster(r, crs = "+proj=longlat")


plot(rr)
mdg_data<-as.data.frame(rr, xy = TRUE, na.rm=T)
dplyr::glimpse(mdg_data)
names(mdg_data)<-c("x", "y", 
                   'INS_HIST2000_MDG',
                   'INS_SSP52050_MDG',
                   'INS_SSP52100_MDG',
                   'VEL_HIST2000_MDG',
                   'VEL_SSP52050_MDG',
                   'VEL_SSP52100_MDG')

library(Hmisc)
mdg_data$clim_brk<-as.numeric(as.character(Hmisc::cut2(mdg_data$VEL_SSP52050_MDG, g = 100, levels.mean = TRUE)))
mdg_data$luse_brk<-as.numeric(as.character(Hmisc::cut2(mdg_data$INS_SSP52050_MDG, g = 100, levels.mean = TRUE)))
mdg_data$clim_brk<-scales::rescale(log(mdg_data$clim_brk), to = c(0.01, 1))
mdg_data$luse_brk<-scales::rescale(log(mdg_data$luse_brk+0.01), to = c(0.01, 1))

mdg_data<-mdg_data%>%
  filter(is.finite(clim_brk), is.finite(luse_brk))%>%
  mutate(mix1 = rgb(red = round(luse_brk,2), green = round(clim_brk, 2), blue = 0.5))


#Load administrative boundaries
focusISO3<-readOGR("./shps/iso3.shp")%>%spTransform(CRS="+proj=longlat")%>%st_as_sf(focusISO3)

#load sectors
mgt_units<-readOGR("./shps/wio_sectors.shp")%>%st_as_sf(mgt_units)
crs(mgt_units)

#Plot
ggplot()+
  geom_sf(data = focusISO3, aes(fill = NA), lwd=.5, colour = "grey70")+
  geom_sf(data = mgt_units, aes(fill = NA), lwd=.5, colour = "grey50")+
  geom_tile(data = mdg_data, aes(x= x, y = y, fill = mix1))+
  scale_fill_identity()+theme_sleek(base_size = 18)+
  labs(x="", y="")+
  coord_sf(xlim = c(25, 55), ylim = c(0, -30), expand = TRUE)
ggsave("mgdmain.png", dpi = 600, width = 7.03, height = 8.29)


ggplot()+
  geom_sf(data = focusISO3, aes(fill = NA), lwd=.5, colour = "grey70")+
  geom_tile(data = mdg_data, aes(x= x, y = y, fill = mix1))+
  scale_fill_identity()+theme_sleek(base_size = 18)+
  labs(x="", y="")+
  coord_sf(xlim = c(30, 50), ylim = c(-10, -28), expand = TRUE)
ggsave("mgd1.png", dpi = 600, width = 7.03, height = 8.29)


ggplot()+
  geom_sf(data = focusISO3, aes(fill = NA), lwd=.5, colour = "grey70")+
  geom_tile(data = mdg_data, aes(x= x, y = y, fill = mix1))+
  scale_fill_identity()+theme_sleek(base_size = 18)+
  scale_x_continuous(breaks = c(40, 45))+
  scale_y_continuous(breaks = c(0, -5, -10, -15))+
  labs(x="", y="")+coord_sf(xlim = c(38, 45), ylim = c(-10, 0), expand = T)
ggsave("mgd2.png", dpi = 600, width = 4, height = 8.29)


#Plot legend
all_choro <- expand.grid(luhP=1:100, cchP=1:100)
all_choro <- all_choro %>%
  select(everything()) %>%
  mutate(mix1 = rgb(red = round(luhP/100, 2), green = round(cchP/100, 2), blue = 0.4))

ggplot(all_choro, aes(x = luhP, y = cchP)) + 
  geom_raster(aes(fill = mix1)) + scale_fill_identity() +
  theme_sleek(base_size = 14)+
  scale_x_continuous(expand = c(0,0), breaks = c(1,100), labels = c(0, 100))+
  scale_y_continuous(expand = c(0,0), breaks = c(1,100), labels = c(0, 100))+
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(angle = 90),
        plot.background = element_rect(fill = NA))+
  labs(x = "Land-use instability (km/yr, log-scaled)", y = "Climate velocity (km/yr, log-scaled)")
ggsave(here::here("lgd.tiff"), dpi = 600, width = 5, height = 5)



#Prepare data for Bleaching Event Plots
#create only spatial polygon to mask the study area for further analysis
#https://julyonr.wordpress.com/2015/10/23/create-a-rectangle-spatialpolygonsdataframe-from-scratch/
x1 = 25
x2 = 70
y1 = 0
y2 = -30

myPolygon = Polygon(cbind(c(x1,x1,x2,x2,x1),c(y1,y2,y2,y1,y1)))
myPolygons = Polygons(list(myPolygon), ID = "A")
SpPolygon = SpatialPolygons(list(myPolygons))
plot(SpPolygon)
proj4string(SpPolygon) = "+proj=longlat +datum=WGS84 +no_defs"
df = matrix(data = c(0))
rownames(df) = "A"
spp = SpatialPolygonsDataFrame(SpPolygon, data= as.data.frame(df))


mdg_rasta<-raster("C:/Users/45019738/Dropbox/MDG_CVA/RiskSlides/asb_year_0_ssp585.tif")
mdg_rasta<-mask(mdg_rasta, spp)
mdg_rasta<-as.data.frame(mdg_rasta, xy = TRUE, na.rm=TRUE)
#mdg_rasta$asb_year_0_ssp585<-as.numeric(as.character(mdg_rasta$asb_year_0_ssp585))
str(mdg_rasta)
ggplot()+
  geom_sf(data = focusISO3, aes(), fill = NA, lwd=.5, colour = "grey70")+
  geom_tile(data = mdg_rasta, aes(x= x, y = y, fill = asb_year_0_ssp585))+
  scale_fill_viridis_c(name = "Year", option = "plasma", 
                       breaks = c(2015, 2020, 2025, 2030, 2035, 2040, 2045, 2050,
                                  2055, 2060, 2065)) +
  theme_sleek(base_size = 18)+labs(x="", y="")+
  theme(legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(.5, 'cm'))+
  coord_sf(xlim = c(25, 60), ylim = c(0, -30), expand = TRUE)
ggsave("mgdYR.png", dpi = 600, width = 8.5, height = 9)



mdg_rasta2<-raster("C:/Users/45019738/Dropbox/MDG_CVA/RiskSlides/diff_ssp245_ssp585.tif")
mdg_rasta2<-mask(mdg_rasta2, spp)
mdg_rasta2<-as.data.frame(mdg_rasta2, xy = TRUE, na.rm=TRUE)
ggplot()+
  geom_sf(data = focusISO3, aes(), fill = NA, lwd=.5, colour = "grey70")+
  geom_tile(data = mdg_rasta2, aes(x= x, y = y, fill = diff_ssp245_ssp585))+
  scale_fill_viridis_c(name = "Difference\n(years)", option = "D") +
  theme_sleek(base_size = 18)+labs(x="", y="")+
  theme(legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(.5, 'cm'))+
  coord_sf(xlim = c(25, 60), ylim = c(0, -30), expand = TRUE)
ggsave("mgddfiffYR.png", dpi = 600, width = 8.5, height = 9)

