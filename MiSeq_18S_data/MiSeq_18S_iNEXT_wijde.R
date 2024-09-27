## MiSeq iNEXT - Wijde

library(httr)
library(iNEXT)
library(ggplot2)
library(readxl)
library(tidyverse)
library(dplyr)

## set wd
wd <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/")


WijdeF02 <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_Wijde.xlsx", sheet = "Wijde_F02")
WijdeF3 <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_Wijde.xlsx", sheet = "Wijde_F3")
WijdeNP <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_Wijde.xlsx", sheet = "Wijde_NP")

#str(WijdeF02)

WijdeF02 <- as.data.frame(WijdeF02)
WijdeF3 <- as.data.frame(WijdeF3)
WijdeNP <- as.data.frame(WijdeNP)


q123_WijdeF02 <- iNEXT(WijdeF02, q = c(0,1,2), datatype = "abundance", endpoint = 30000)
q123_WijdeF3 <- iNEXT(WijdeF3, q = c(0,1,2), datatype = "abundance", endpoint = 6000000)
q123_WijdeNP <- iNEXT(WijdeNP, q = c(0,1,2), datatype = "abundance", endpoint = 15000)



p_F02 <- ggiNEXT(q123_WijdeF02, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#ffe75b", "#e9e370", "#d3df84", "#bddb99", "#a7d7ad",
                                                                                                                                 "#90d2c2", "#7aced6", "#64caeb", "#4ec6ff")) +
  labs(title = "Hill numbers - Wijdefjorden - MiSeq 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "F02 filters \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                           axis.title=element_text(size=19), #change font size of axis titles
                           plot.title=element_text(size=26), #change font size of plot title
                           plot.subtitle = element_text(size=22), 
                           legend.text=element_text(size=16), #change font size of legend text
                           legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("MiSeq_Wijde_F02.png", p_F02, path = wd, width = 15, height = 15, dpi = 600)


p_F3 <- ggiNEXT(q123_WijdeF3, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#e87eff", "#eb86f3", "#ee8de8", "#f195dc", "#f49cd1",
                                                                                                                               "#f6a4c5", "#f9abb9", "#fcb3ae", "#ffbaa2")) + 
  labs(title = "Hill numbers - Wijdefjorden - MiSeq 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "F3 filters \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )


ggsave("MiSeq_Wijde_F3.png", p_F3, path = wd, width = 15, height = 15, dpi = 600)

p_NP <- ggiNEXT(q123_WijdeNP, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#6fa940", "#81ae41", "#93b342", "#a5b842", "#b7bd43",
                                                                                                                               "#c9c244", "#bdc745", "#edcc45", "#ffd146")) + 
  labs(title = "Hill numbers - Wijdefjorden - MiSeq 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "Pump samples \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )


ggsave("MiSeq_Wijde_NP.png", p_NP, path = wd, width = 15, height = 15, dpi = 600)








