## MiSeq iNEXT - Rijp

library(httr)
library(iNEXT)
library(ggplot2)
library(readxl)
library(tidyverse)
library(dplyr)

## set wd
wd <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/")


RijpF02 <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_Rijp.xlsx", sheet = "Rijp_F02")
RijpF3 <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_Rijp.xlsx", sheet = "Rijp_F3")
RijpNP <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_Rijp.xlsx", sheet = "Rijp_P")

str(RijpF02)

RijpF02 <- as.data.frame(RijpF02)
RijpF3 <- as.data.frame(RijpF3)
RijpNP <- as.data.frame(RijpNP)


q123_RijpF02 <- iNEXT(RijpF02, q = c(0,1,2), datatype = "abundance", endpoint = 450000)
q123_RijpF3 <- iNEXT(RijpF3, q = c(0,1,2), datatype = "abundance", endpoint = 100000)
q123_RijpNP <- iNEXT(RijpNP, q = c(0,1,2), datatype = "abundance", endpoint = 1000)


head(q123_RijpF02)
head(q123_RijpF3)
head(q123_RijpNP)



p_F02 <- ggiNEXT(q123_RijpF02, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#0c7edc", "#2a68d6", "#4952d0", "#673cca", "#8626c4", "#a410be")) +
                                                                                                                                   #geom_line(aes(linetype = "lty"), linewidth = 1.5) + 
  labs(title = "Hill numbers - Rijpfjorden - MiSeq 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "F02 filters \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("MiSeq_Rijp_F02_2.png", p_F02, path = wd, width = 15, height = 15, dpi = 600)


p_F3 <- ggiNEXT(q123_RijpF3, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#dca27d", "#b9a874", "#96ad6b", "#73b363", "#50b85a", "#2dbe51")) + 
  labs(title = "Hill numbers - Rijpfjorden - MiSeq 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "F3 filters \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("MiSeq_Rijp_F3.png", p_F3, path = wd, width = 15, height = 15, dpi = 600)

p_NP <- ggiNEXT(q123_RijpNP, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#3ddcba", "#64be9c")) + 
   labs(title = "Hill numbers - Rijpfjorden - MiSeq 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "Pump samples \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("MiSeq_Rijp_NP.png", p_NP, path = wd, width = 15, height = 15, dpi = 600)








