## MiSeq iNEXT - van Mijenfjorden

library(httr)
library(iNEXT)
library(ggplot2)
library(readxl)
library(tidyverse)
library(dplyr)

## set wd
wd <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/")

vMijenF02 <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_vMijen.xlsx", sheet = "vMijen_F02")
vMijenF3 <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_vMijen.xlsx", sheet = "vMijen_F3")
vMijenNP <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_vMijen.xlsx", sheet = "vMijen_NP")

vMijenF02 <- as.data.frame(vMijenF02)
vMijenF3 <- as.data.frame(vMijenF3)
vMijenNP <- as.data.frame(vMijenNP)

q123_vMijenF02 <- iNEXT(vMijenF02, q = c(0,1,2), datatype = "abundance", endpoint = 40000)
q123_vMijenF3 <- iNEXT(vMijenF3, q = c(0,1,2), datatype = "abundance", endpoint = 1100000)
q123_vMijenNP <- iNEXT(vMijenNP, q = c(0,1,2), datatype = "abundance", endpoint = 100000)                          


p_F02 <- ggiNEXT(q123_vMijenF02, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#71d2e4", "#6cc5e1", "#66b9df", "#61acdc", "#5ca0da", "#5793d7",
                                                                                                                                  "#5187d5", "#4c7ad2", "#476ed0", "#4261cd", "#3c55cb", "#3748c8")) + 
  labs(title = "Hill numbers - van Mijenfjorden - MiSeq 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "F02 filters \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("MiSeq_vMijen_F02_2.jpeg", p_F02, path = wd, width = 15, height = 15, dpi = 600)


p_F3 <- ggiNEXT(q123_vMijenF3, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#34e4a1", "#40d2aa", "#4cc1b2", "#59afbb", "#659ec3", "#718ccc", 
                                                                                                                                "#7d7bd4", "#8969dd", "#9558e5", "#a246ee", "#ae35f6", "#ba23ff")) + 
  labs(title = "Hill numbers - van Mijenfjorden - MiSeq 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "F3 filters \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("MiSeq_vMijen_F3.jpeg", p_F3, path = wd, width = 15, height = 15, dpi = 600)


p_NP <- ggiNEXT(q123_vMijenNP, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#e4bf6f", "#e5b67b", "#e5ad87", "#e6a493", "#e79b9f", "#e792ab", "#e889b7",
                                                                                                                                "#e980c3", "#e977cf", "#ea6edb", "#eb65e7", "#eb5cf3", "#ec53ff")) + 
  labs(title = "Hill numbers - van Mijenfjorden - MiSeq 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "Net and pump samples \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("MiSeq_vMijen_NP_2.jpeg", p_NP, path = wd, width = 15, height = 15, dpi = 600)




