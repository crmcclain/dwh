################necessary packages################
  require(adespatial)
  require(dplyr)
  require(ggplot2)
  require(vegan)
  require(DescTools)
  library(tibble)
  require(ggrepel)
  require(gridExtra)
  source("~/Dropbox/DWH/Final Paper/evplot.R")
  source("~/Dropbox/DWH/Final Paper/clean.plot.R")

################get data################
  setwd("~/Dropbox/Return of the Woodfall/Data")
  DWH <- data.frame(read.csv("transectbytaxa_refined.csv", header=TRUE, stringsAsFactors=FALSE))
  glimpse(DWH)
  
  DWH$Group<-as.factor(DWH$Group)
  DWH$Group2<-as.factor(DWH$Group2)

  DWH_taxa <- DWH %>% 
    dplyr::select(Acanthonus_armatus:Enypniastes_eximia)
  
################refine data################
  #remove stations that have zero total abundance
  DWH_nozero <- DWH %>%
    filter(Total.Abundance>0) 
  
  #view dataset
  glimpse(DWH_nozero)
  levels(DWH_nozero$Group)
  levels(DWH_nozero$Group2)
  
  summary(DWH_nozero$Group)
  summary(DWH_nozero$Group2)
  
################data analysis################
  
  #create data frame to hold beta diversity estimates
  Beta_results <- data.frame(c("DWH2017", "DWH2010", "DWH2017 500m", "DWH 2017 2000m", 
                         "DWH 2010 500m" , "DWH 2010 2000m", 
                         "Background 1", "Background 2", "Background 3", "Background 4"))

  #DWH2017
    #all 2017 DWH transects
    DWH2017 <-  DWH_nozero %>%
      filter(Group=="DWH2017") %>%
      dplyr::select(Acanthonus_armatus:Enypniastes_eximia)
    #beta div
    DWH2017_beta <-beta.div(DWH2017, method="hellinger" ,nperm = 999)
    Beta_results[1,2] <- DWH2017_beta$beta[2]
    Beta_results[1,3] <- length(DWH2017_beta$LCBD)
    colnames(Beta_results) <- c("Group","BDTotal", "Transects")
  
  #DWH2010
  
    #all 2010 DWH transects
    DWH2010 <-  DWH_nozero %>%
      filter(Group=="DWH2010") %>%
      dplyr::select(Acanthonus_armatus:Enypniastes_eximia)
    #beta div
    DWH2010_beta <-beta.div(DWH2010, method="hellinger" ,nperm = 999)
    Beta_results[2,2] <- DWH2010_beta$beta[2]
    Beta_results[2,3] <- length(DWH2010_beta$LCBD)
  
  #DWH2017 500m
    #get transects
    DWH2017_500 <-  DWH_nozero %>%
      filter(Group2=="DWH2017 500m") %>%
      dplyr::select(Acanthonus_armatus:Enypniastes_eximia)
    #beta div
    DWH2017_500_beta <-beta.div(DWH2017_500, method="hellinger" ,nperm = 999)
    Beta_results[3,2] <- DWH2017_500_beta$beta[2]
    Beta_results[3,3] <- length(DWH2017_500_beta$LCBD)
  
  #DWH 2017 2000m
    #get transects
    DWH2017_2000 <-  DWH_nozero %>%
      filter(Group2=="DWH 2017 2000m") %>%
      dplyr::select(Acanthonus_armatus:Enypniastes_eximia)
    #beta div
    DWH2017_2000_beta <-beta.div(DWH2017_2000, method="hellinger" ,nperm = 999)
    Beta_results[4,2] <- DWH2017_2000_beta$beta[2]
    Beta_results[4,3] <- length(DWH2017_2000_beta$LCBD)
    
  #DWH 2010 500m
    #get transects
    DWH2010_500 <-  DWH_nozero %>%
      filter(Group2=="DWH 2010 500m") %>%
      dplyr::select(Acanthonus_armatus:Enypniastes_eximia)
    #beta div
    DWH2010_500_beta <-beta.div(DWH2010_500, method="hellinger" ,nperm = 999)
    Beta_results[5,2] <- DWH2010_500_beta$beta[2]
    Beta_results[5,3] <- length(DWH2010_500_beta$LCBD)
  
  #DWH 2010 2000m
    DWH2010_2000 <-  DWH_nozero %>%
      filter(Group2=="DWH 2010 2000m") %>%
      dplyr::select(Acanthonus_armatus:Enypniastes_eximia)
    #beta div
    DWH2010_2000_beta <-beta.div(DWH2010_2000, method="hellinger" ,nperm = 999)
    Beta_results[6,2] <- DWH2010_2000_beta$beta[2]
    Beta_results[6,3] <- length(DWH2010_2000_beta$LCBD)
    
  #Background 1
    WF1 <-  DWH_nozero %>%
      filter(Group2=="Background 1") %>%
      dplyr::select(Acanthonus_armatus:Enypniastes_eximia)
    #beta div
    WF1_beta <-beta.div(WF1, method="hellinger" ,nperm = 999)
    Beta_results[7,2] <- WF1_beta$beta[2]
    Beta_results[7,3] <- length(WF1_beta$LCBD)
  
  #Background 2
    WF2 <-  DWH_nozero %>%
      filter(Group2=="Background 2") %>%
      dplyr::select(Acanthonus_armatus:Enypniastes_eximia)
    #beta div
    WF2_beta <-beta.div(WF2, method="hellinger" ,nperm = 999)
    Beta_results[8,2] <- WF2_beta$beta[2]
    Beta_results[8,3] <- length(WF2_beta$LCBD)
    
  #Background 3
    WF3 <-  DWH_nozero %>%
      filter(Group2=="Background 3") %>%
      dplyr::select(Acanthonus_armatus:Enypniastes_eximia)
    #beta div
    WF3_beta <-beta.div(WF3, method="hellinger" ,nperm = 999)
    Beta_results[9,2] <- WF3_beta$beta[2]
    Beta_results[9,3] <- length(WF3_beta$LCBD)
  
  #Background 4
    WF4 <-  DWH_nozero %>%
      filter(Group2=="Background 4") %>%
      dplyr::select(Acanthonus_armatus:Enypniastes_eximia)
    #beta div
    WF4_beta <-beta.div(WF4, method="hellinger" ,nperm = 999)
    Beta_results[10,2] <- WF4_beta$beta[2]
    Beta_results[10,3] <- length(WF4_beta$LCBD)
  

################multivariate################
    DWH_nozero_taxa <- DWH_nozero %>%
      dplyr::select(Acanthonus_armatus:Enypniastes_eximia)
    rownames(DWH_nozero_taxa)<-DWH_nozero$Transect.Number
    
    #Hellinger pre-transformation of the species matrix
    dwh.h <- decostand (DWH_nozero_taxa, "hellinger")
    dwh.h.pca <- rda(dwh.h)
    dwh.h.pca.summ <- summary(dwh.h.pca)

    #plot eigenvlues and % variance for each axis
    ev <- dwh.h.pca$CA$eig
    evplot(ev)
    
    #PCA biplots
    cleanplot.pca(dwh.h.pca, ahead=0)
    
    #pull coordinates and loadings
    df1  <- data.frame(dwh.h.pca.summ$sites[,1:2])
    DWH_nozero$PC1 <- df1[,1]
    DWH_nozero$PC2 <- df1[,2]
    
    #pull loadings
    DWH.loadings  <- data.frame(dwh.h.pca.summ$species[,1:2])
    
    #filtering out top species contributions
    DWH.loadings.top <- DWH.loadings %>%
      rownames_to_column('name') %>% #need to preserve labels
      filter(PC1<=-0.25 | PC1>=0.25 | PC2<=-0.25 | PC2>=0.25) %>%
      column_to_rownames('name')
    
    #analyses
    dwh.h.dl <- dist(dwh.h)
    pc_test<- betadisper(dwh.h.dl, DWH_nozero$Group)
    anova(pc_test)
    permutest(pc_test, pairwise = TRUE)
    
################alpha diversity and abundance ################  
    require(multcomp)
    require(lme4)
    require(car)
    require(WordR)

    
    DWH$H <- diversity(DWH_taxa, index = "shannon")
    DWH$ES5 <- rarefy(DWH_taxa, 5)
    DWH$ES5 <- rarefy(DWH_taxa, 5)
    
    ####
    DWH_H_Group2<- lmer(H  ~Group +  + (1|Site), data = DWH, na.action=na.fail)
    setwd("~/Desktop")
    sink("H anova.txt")
    summary(DWH_H_Group2)
    Anova(DWH_H_Group2, type = "II")
    summary(glht(DWH_H_Group2, mcp(Group="Tukey")))
    sink()
    ####
    
    DWH_abundance_Group<- lmer(Total.Abundance  ~Group +  + (1|Site), data = DWH, na.action=na.fail)
    setwd("~/Desktop")
    sink("abundance anova.txt")
    DWH_abundance_Group<- lmer(Total.Abundance  ~Group +  + (1|Site), data = DWH, na.action=na.fail)
    summary(DWH_abundance_Group)
    Anova(DWH_abundance_Group, type = "II")
    summary(glht(DWH_abundance_Group, mcp(Group="Tukey")))
    sink()
    
    DWH_abundance_Group2<- lmer(Total.Abundace2  ~Group +  + (1|Site), data = DWH, na.action=na.fail)
    setwd("~/Desktop")
    sink("abundance2 anova.txt")
    summary(DWH_abundance_Group2)
    Anova(DWH_abundance_Group2, type = "II")
    summary(glht(DWH_abundance_Group2, mcp(Group="Tukey")))
    sink()
    
################plots################

    #need to reorder the levels because R package sucks    
    levels(Beta_results$Group)
    Beta_results$Group <- reorder(Beta_results$Group, 
                                new.order=c("DWH2010","DWH 2010 500m" , "DWH 2010 2000m",
                                            "DWH2017", "DWH2017 500m", "DWH 2017 2000m", 
                                            "Background 1", "Background 2", "Background 3", "Background 4"))
    DWH_nozero$Group2 = factor(DWH_nozero$Group2,levels(DWH_nozero$Group2)
                              [c(6,5,8,7,1:4)])
    DWH$Group2 = factor(DWH$Group2,levels(DWH$Group2)
                               [c(6,5,8,7,1:4)])
    
    
    #designate the color pallette
    study_color_pallette <- scale_color_manual(values = c("DWH2010" = "indianred1",
                                                      "DWH 2010 500m" = "indianred2" , 
                                                      "DWH 2010 2000m" = "indianred3",
                                                      "DWH2017" = "orange1", 
                                                      "DWH2017 500m" = "orange2", 
                                                      "DWH 2017 2000m" = "orange3", 
                                                      "Background 1" = "steelblue1", 
                                                      "Background 2" = "steelblue2", 
                                                      "Background 3" = "steelblue3", 
                                                      "Background 4" = "steelblue4",
                                                      "2017Background" = "steelblue1"))
    study_color_pallette2 <- scale_fill_manual(values = c("DWH2010" = "indianred1",
                                                          "DWH 2010 500m" = "indianred2" , 
                                                          "DWH 2010 2000m" = "indianred3",
                                                          "DWH2017" = "orange1", 
                                                          "DWH2017 500m" = "orange2", 
                                                          "DWH 2017 2000m" = "orange3", 
                                                          "Background 1" = "steelblue1", 
                                                          "Background 2" = "steelblue2", 
                                                          "Background 3" = "steelblue3", 
                                                          "Background 4" = "steelblue4",
                                                          "2017Background" = "steelblue1"))
   #Principale Coordinates plot 
    
   p12 <- ggplot(DWH_nozero, aes(x=PC1, y=PC2, fill=Group, color=Group)) + 
      study_color_pallette+
      study_color_pallette2+
      geom_text_repel(aes(label=Transect.Number),size=4) +
      geom_hline(yintercept=0, linetype="dotted") +
      geom_vline(xintercept=0, linetype="dotted") +
      geom_text(data=DWH.loadings.top, aes(x=PC1, y=PC2,label=rownames(DWH.loadings.top)),
                                    color="grey30",
                                    inherit.aes = FALSE)+
     ylab('PC 1')+
     xlab('PC 2')+
     ggtitle("B")+
     theme_bw(base_size=12)+
     theme(axis.line = element_line(colour = "grey70"),
           panel.grid.minor = element_blank(),
           panel.grid.major = element_blank(),
           legend.position="none",
           plot.title = element_text(lineheight=.8, face="bold", hjust = 0))

    
  #plot of group beta diversity       
  p13 <- ggplot(Beta_results, aes(y=BDTotal, x=Group, fill=Group, color=Group, label=Transects))+
      geom_point(cex=3) + 
      geom_text(hjust=-1)+
      geom_segment(aes(x=Group, 
                       xend=Group, 
                       y=0, 
                       yend=BDTotal))+
      study_color_pallette+
      study_color_pallette2+
      ylab('Total Within Group Beta-Diversity')+
      xlab('Location')+
      ggtitle("C")+
      theme_bw(base_size=12)+
      theme(axis.line = element_line(colour = "darkgrey"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            legend.position="none",
            axis.text.x = element_text(angle = -90),
            plot.title = element_text(lineheight=.8, face="bold", hjust = 0))
  
  #boxplot of total abundances 
  Transect2017 <- DWH %>%
    filter(Group=="DWH2017")
  
  p14 <- ggplot(DWH, aes(y=Total.Abundance, x=Group2,color=Group2, fill=Group2))+
    geom_boxplot(alpha=.5)+
    geom_jitter()+
    geom_boxplot(data=Transect2017, aes(y=Total.Abundace2, x=Group2,color=Group2, fill=Group2),
                 alpha=.5, inherit.aes = FALSE)+
    ggtitle("D")+
    study_color_pallette+
    study_color_pallette2+
    ylab('Total Abundance')+
    xlab('Location')+
    theme_bw(base_size=12)+
    theme(axis.line = element_line(colour = "darkgrey"),
          panel.grid.minor = element_blank(),
          legend.position="none",
          axis.text.x = element_text(angle = -90),
          plot.title = element_text(lineheight=.8, face="bold", hjust = 0))
  
  #boxplot of H' 
  p11 <- ggplot(DWH, aes(y=H, x=Group2,color=Group2, fill=Group2))+
    geom_boxplot(alpha=.5)+
    geom_jitter()+
    ggtitle("A")+
    study_color_pallette+
    study_color_pallette2+
    ylab('Shannon H')+
    xlab('Location')+
    theme_bw(base_size=12)+
    theme(axis.line = element_line(colour = "darkgrey"),
          panel.grid.minor = element_blank(),
          legend.position="none",
          axis.text.x = element_text(angle = -90),
          plot.title = element_text(lineheight=.8, face="bold", hjust = 0))
  

  setwd("~/Desktop")
  pdf(file="figure1.pdf",width=11, height=8.5, useDingbats=FALSE)
  par(mar=c(5,3,2,2)+0.1) #removes space from around edges of pdf
  grid.arrange(p11,p12,p13,p14, ncol=2)
  dev.off()
  
  ###########
  
  citation("lme4")
  citation("car")
  citation("multcomp")
  
  
