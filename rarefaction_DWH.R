## import packages
library(iNEXT)
library(ggplot2)
library(dplyr)
library(vegan)
library(directlabels)

theme_craig <- function () { 
  theme_bw(base_size=12) %+replace% 
    theme(
      # change stuff here
      axis.line = element_line(colour = "darkgrey"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      legend.position="none",
      plot.title = element_text(lineheight=.8, face="bold", hjust = 0))
}


#get sample data
setwd("~/Dropbox/Return of the Woodfall/Data")
DWH <- data.frame(read.csv("transectbytaxa_refined.csv", header=TRUE, stringsAsFactors=FALSE))
glimpse(DWH)

DWH$Group<-as.factor(DWH$Group)
DWH$Group2<-as.factor(DWH$Group2)

DWH_nozero <- DWH %>%
  filter(Total.Abundance>0) %>%
  group_by(Group2) %>% #groups by variable site
  select(Group2:Enypniastes_eximia) %>% #drops the columns group and sub.group 
  summarise_all(sum) %>% #for each species column takes the sum
  select(-Group2) 
  


#bit of code because iNext requires code in ridiculous format with everything switched
DWH_t <- data.frame(t(DWH_nozero)) #transpose and then data frame
colnames(DWH_t) <- c("Background 1", "Background 2", "Background 3", "Background 4", 
                         "DWH 2010 2000-S", "DWH 2010 500-N", 
                         "DWH 2017 2000-S", "DWH 2017 500-N") #relabel the column heads
str(DWH_t)

#rarefaction analysis and plot
out1 <- iNEXT(DWH_t, q=0, datatype="abundance")
ggiNEXT(out1, type=1) + theme_bw(base_size=10) 

#get everything on the same plot with even more tedious and ridicolous code
#this is not worth explaining other than it extracts out the bits of data needed from the
#mess that iNext spits out so we can plot it how would like.

df_deep <- fortify(out1, type=1)

df.point_deep <- df_deep[which(df_deep$method=="observed"),]
df.line_deep <- df_deep[which(df_deep$method!="observed"),]
df.line_deep$method <- factor(df.line_deep$method, 
                              c("interpolated", "extrapolated"),
                              c("interpolation", "extrapolation"))



#designate the color pallette
study_color_pallette <- scale_color_manual(values = c("DWH 2010 500-N" = "indianred2" , 
                                                      "DWH 2010 2000-S" = "indianred3",
                                                      "DWH 2017 500-N" = "orange2", 
                                                      "DWH 2017 2000-S" = "orange3", 
                                                      "Background 1" = "steelblue1", 
                                                      "Background 2" = "steelblue2", 
                                                      "Background 3" = "steelblue3", 
                                                      "Background 4" = "steelblue4"))
study_color_pallette2 <- scale_fill_manual(values = c("DWH 2010 500-N" = "indianred2" , 
                                                      "DWH 2010 2000-S" = "indianred3",
                                                      "DWH 2017 500-N" = "orange2", 
                                                      "DWH 2017 2000-S" = "orange3", 
                                                      "Background 1" = "steelblue1", 
                                                      "Background 2" = "steelblue2", 
                                                      "Background 3" = "steelblue3", 
                                                      "Background 4" = "steelblue4"))

fig1 <- ggplot(df_deep, aes(x=x, y=y, color=site, fill=site)) + 
  study_color_pallette+ #says we want to use our own colors for outlines and points
  study_color_pallette2+ #use our own colors for fills 
  #geom_point(aes(shape=site), size=5, data=df.point_deep) + #plots the end points of the curves
  geom_line(aes(linetype=method), lwd=1.5, data=df.line_deep) + #plots acutal lines
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                  fill=site, colour=NULL), alpha=0.2) + #plots envelopes
  labs(x="Number of Individuals", y="Expected Number of Species") + #labelling the axis
  xlim(0,550)+
  ylim(1,18)+
  geom_dl(aes(label = site), method = list(dl.combine("last.points"), cex = 0.8))+
  theme_craig() #swithcing the the graph theme to something plain

setwd("~/Desktop")
ggsave("figure1.pdf", width = 6, height = 4)

