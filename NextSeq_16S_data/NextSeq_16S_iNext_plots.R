## iNEXT Plots NextSeq 16S Prokaryotes
##

library(httr)
library(iNEXT)
library(ggplot2)
library(readxl)
library(tidyverse)
library(dplyr)

## 
wd <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/iNext_analysis/")

vMijen16S <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/iNext_analysis/NextSeq_16S_ampli_cleaned.xlsx", sheet = "vMijenfjorden")
Kongs16S <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/iNext_analysis/NextSeq_16S_ampli_cleaned.xlsx", sheet = "Kongsfjorden")
Wijde16S <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/iNext_analysis/NextSeq_16S_ampli_cleaned.xlsx", sheet = "Wijdefjorden")
Rijp16S <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/iNext_analysis/NextSeq_16S_ampli_cleaned.xlsx", sheet = "Rijpfjorden")

vMijen16S <- as.data.frame(vMijen16S)
Kongs16S <- as.data.frame(Kongs16S)
Wijde16S <- as.data.frame(Wijde16S)
Rijp16S <- as.data.frame(Rijp16S)

q123vMijen16S <- iNEXT(vMijen16S, q=c(0,1,2), datatype = "abundance", endpoint = 750000)
q123Kongs16S <- iNEXT(Kongs16S, q = c(0,1,2), datatype = "abundance", endpoint = 150000)
q123Wijde16S <- iNEXT(Wijde16S, q = c(0,1,2), datatype = "abundance", endpoint = 200000)
q123Rijp16S <- iNEXT(Rijp16S, q = c(0,1,2), datatype = "abundance", endpoint = 200000)

## Plots

p_Mi16 <- ggiNEXT(q123vMijen16S, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#fc35ff", "#fd5bde", "#fd71cc", "#fd82bd", "#fe91b0", "#fe9fa4",
                                                                                                                                  "#feab99", "#feb88e", "#fec384", "#ffce7b", "#ffd971", "#ffe468")) +
  labs(title = "Hill numbers - van Mijenfjorden", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "NextSeq 16S \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq_16S_vMijen_F02.jpeg", p_Mi16, path = wd, width = 15, height = 15, dpi = 600)


p_Ko16 <- ggiNEXT(q123Kongs16S, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#ff6009", "#e56d0a", "#cc7a0b", "#b2870c", "#99940d", "#7fa103",
                                                                                                                                 "#65ae0f", "#4cbb10", "#32c811")) +
  labs(title = "Hill numbers - Kongsfjorden", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "NextSeq 16S \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq_16S_Kongs_F02.jpeg", p_Ko16, path = wd, width = 15, height = 15, dpi = 600)


p_Wi16 <- ggiNEXT(q123Wijde16S, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#6346e4", "#5a65c9", "#547bb6", "#4f8da6", "#4a9e98", "#45ae8a",
                                                                                                                                 "#41be7c", "#3ccd6f", "#38dc62")) +
  labs(title = "Hill numbers - Wijdefjorden", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "NextSeq 16S \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq_16S_Wijde_F02.jpeg", p_Wi16, path = wd, width = 15, height = 15, dpi = 600)


p_Ri16 <- ggiNEXT(q123Rijp16S, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#e41145", "#e2394b", "#e16251", "#df8a58", "#f6c03d", "#eee038")) +
   labs(title = "Hill numbers - Rijpfjorden", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "NextSeq 16S \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq_16S_Rijp_F02.jpeg", p_Ri16, path = wd, width = 15, height = 15, dpi = 600)
