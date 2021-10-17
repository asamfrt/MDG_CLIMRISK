#Examine risks to PAs of East Africa
mdsPnts<-readOGR("./shps/mgd_wpa_pnts.shp")
mdgraster<-stack(list.files("./raster", pattern = ".tif$", full.names = TRUE))
plot(mdgraster[[1]])

#extract velocity values to points
mdgExtract<-extract(mdgraster, mdsPnts, sp=TRUE, fun=mean)
mdgExtract.df<-as.data.frame(mdgExtract)
glimpse(na.omit(mdgExtract.df))

#Compute and store deviation of projected estimates (2021-2050) from global median of baseline (1971-200) for both metrics
data_fnl <- mdgExtract.df %>% select(everything())%>%
  mutate(clm852050 = log(VEL_SSP52050_MDG/median(VEL_HIST2000_MDG, na.rm = TRUE)), 
         lnd852050 = log(INS_SSP52050_MDG/median(INS_HIST2000_MDG, na.rm = TRUE)))

#Check data distribution 
hist(data_fnl$clm852050, breaks = 100)
hist(data_fnl$lnd852050, breaks = 100)
with(data_fnl, cor.test(lnd852050, clm852050, method = "spearman")) #-0.18173; S = 2.0554e+15, p-value < 2.2e-16

#Normalised estimates
data_fnl$resLUI<-scales::rescale(data_fnl$lnd852050, to = c(0.01, 1))
data_fnl$resCCV<-scales::rescale(data_fnl$clm852050, to = c(0.01, 1)) 

dataFig3 <- data_fnl %>% select(everything())%>%filter(is.finite(resLUI), is.finite(resCCV))%>%
  mutate(catLUI = ifelse(resLUI > .5, 0.8, .2), catCCV = ifelse(resCCV > .5, 0.8, .2),
         biVar_mix1 = rgb(catLUI, catCCV, 0.4))

#Import PA polygons
wdpaPoly<-readOGR("./shps/mgd_wpa_poly.shp")%>%st_as_sf(wdpaPoly)
wdpaPolyMerge<-merge(wdpaPoly,dataFig3[, c("WDPAID", "biVar_mix1")], by= "WDPAID")

(fig3a<-ggplot()+
    geom_sf(data = focusISO3, aes(fill = NA), lwd=.5, colour = "grey70")+
    geom_sf(data = wdpaPolyMerge, aes(fill = biVar_mix1), lwd=0.001, colour = "grey50")+
    scale_fill_identity()+theme_sleek(base_size = 18)+
    labs(x="", y="")+
    coord_sf(xlim = c(25, 55), ylim = c(10, -35), expand = TRUE))
ggsave("plts/mgdMaina.tiff", dpi = 600, width = 7.03, height = 8.29)

#Estimate percentage of pixels under each bivariate space grouped by IUCN categories
Obs<-nrow(dataFig3) #29345
(dataFig3 %>% select(catLUI, catCCV, biVar_mix1)%>% 
    group_by(catLUI, catCCV, biVar_mix1)%>%
    summarise(N = n(), Percent = round(100*(N/Obs), 0))%>%ungroup()) 

##Mix colours for legend
dchoro <- expand.grid(luhP=1:100/100, cchP=1:100/100)%>%
  mutate(catLUI = ifelse(luhP > .5, 0.8, .2), catCCV = ifelse(cchP > .5, 0.8, .2),
         biVar_mix1 = rgb(catLUI, catCCV, 0.4))

legFig3 <- ggplot(data = NULL) + 
  geom_raster(data=dchoro, aes(x = luhP, y = cchP, fill = biVar_mix1)) + scale_fill_identity() +
  theme_ncc(base_size = 12)+
  scale_x_continuous(expand = c(0, 0), breaks = c(0.01, 1), labels = c(0, 100))+
  scale_y_continuous(expand = c(0, 0), breaks = c(0.01, 1), labels = c(0, 100))+
  theme(axis.text.y = element_text(size = 10),axis.text.x = element_text(size = 10),axis.title.y = element_text(angle = 90))+
  labs(x = "Land-use instabbility relative to \nglobal median of 1971-2000 \n(log-rescale)", 
       y = "Climate velocity relative to \nglobal median of 1971-2000 \n(log-rescale)")+
  geom_point(data = dataFig3, aes(x = resLUI, y = resCCV),
             size=.1, colour = "grey10",shape = 1, alpha= .3) +
  scale_colour_identity()+
  theme(legend.position = "none", plot.margin = unit(c(-5.5,-5.5,4,3), "mm"))+
  annotate("text", label = 'atop(bold("BL"), bold("13%"))', x = .10, y = .10, size = 4.5, colour = "black", parse = T)+
  annotate("text", label = 'atop(bold("TL"), bold("15%"))', x = .10, y = .90, size = 4.5, colour = "black", parse = T)+
  annotate("text", label = 'atop(bold("BR"), bold("53%"))', x = .90, y = .10, size = 4.5, colour = "black", parse = T)+
  annotate("text", label = 'atop(bold("TR"), bold("19%"))', x = .90, y = .90, size = 4.5, colour = "black", parse = T)

#legNtake<-ggExtra::ggMarginal(legNtake, type="histogram", col="grey30", fill="white", size = 3)
densLUI <- ggplot(data = dataFig3,aes(x = resLUI)) + 
  geom_histogram(color = "grey10", fill = "white", size = .5) + 
  scale_x_continuous(expand = c(0,0))+theme_void()+theme(plot.margin = unit(c(0,0,0,-5.5), "mm"))

densCCV <- ggplot(data = dataFig3,aes(x = resCCV)) +
  geom_histogram(color = "grey10", fill = "white", size = .5) + 
  theme_void() + scale_x_continuous(expand = c(0,0))+coord_flip()+theme(plot.margin = unit(c(0,0,0,-5.5), "mm"))

library(patchwork)
(legxx <- densLUI + plot_spacer() + legFig3 + densCCV + plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4)))
ggsave("legend.tiff", dpi = 600, width = 4, height = 4)

#Combine plots and Save
legxx<-cowplot::plot_grid(ggEmpty, legxx, ggEmpty, ncol = 3, rel_widths = c(.2, 2.5, .2))
(plt3<-cowplot::plot_grid(fig3a, legxx, ncol = 1, labels = c("(a)", "(b)"), rel_heights = c(1.5, 2), label_size = 20))
