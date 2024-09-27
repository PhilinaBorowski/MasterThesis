## NextSeq2 - iNEXT - Plots
##

## iNEXT Plots NextSeq 16S Prokaryotes
##

library(httr)
library(iNEXT)
library(ggplot2)
library(readxl)
library(tidyverse)
library(dplyr)

## 
wd <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/iNext_analysis/")


# read in files

vMijenDeep <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/iNext_analysis/NextSeq2_18S_wo_metazoa.xlsx", sheet = "vMijenDeep")
vMijenF02 <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/iNext_analysis/NextSeq2_18S_wo_metazoa.xlsx", sheet = "vMijenF02")
vMijenF3 <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/iNext_analysis/NextSeq2_18S_wo_metazoa.xlsx", sheet = "vMijenF3")
vMijenNP <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/iNext_analysis/NextSeq2_18S_wo_metazoa.xlsx", sheet = "vMijenNP")

KongsDeep <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/iNext_analysis/NextSeq2_18S_wo_metazoa.xlsx", sheet = "KongsDeep")
KongsF02 <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/iNext_analysis/NextSeq2_18S_wo_metazoa.xlsx", sheet = "KongsF02")
KongsF3 <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/iNext_analysis/NextSeq2_18S_wo_metazoa.xlsx", sheet = "KongsF3")
KongsNP <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/iNext_analysis/NextSeq2_18S_wo_metazoa.xlsx", sheet = "KongsNP")

WijdeDeep <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/iNext_analysis/NextSeq2_18S_wo_metazoa.xlsx", sheet = "WijdeDeep")
WijdeF02 <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/iNext_analysis/NextSeq2_18S_wo_metazoa.xlsx", sheet = "WijdeF02")
WijdeF3 <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/iNext_analysis/NextSeq2_18S_wo_metazoa.xlsx", sheet = "WijdeF3")
WijdeNP <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/iNext_analysis/NextSeq2_18S_wo_metazoa.xlsx", sheet = "WijdeNP")

RijpDeep <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/iNext_analysis/NextSeq2_18S_wo_metazoa.xlsx", sheet = "RijpDeep")
RijpF02 <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/iNext_analysis/NextSeq2_18S_wo_metazoa.xlsx", sheet = "RijpF02")
RijpF3 <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/iNext_analysis/NextSeq2_18S_wo_metazoa.xlsx", sheet = "RijpF3")
RijpNP <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/iNext_analysis/NextSeq2_18S_wo_metazoa.xlsx", sheet = "RijpNP")

# bring it into right format

vMijenDeep <- as.data.frame(vMijenDeep)
vMijenF02 <- as.data.frame(vMijenF02)
vMijenF3 <- as.data.frame(vMijenF3)
vMijenNP <- as.data.frame(vMijenNP)

KongsDeep <- as.data.frame(KongsDeep)
KongsF02 <- as.data.frame(KongsF02)
KongsF3 <- as.data.frame(KongsF3)
KongsNP <- as.data.frame(KongsNP)

WijdeDeep <- as.data.frame(WijdeDeep)
WijdeF02 <- as.data.frame(WijdeF02)
WijdeF3 <- as.data.frame(WijdeF3)
WijdeNP <- as.data.frame(WijdeNP)

RijpDeep <- as.data.frame(RijpDeep)
RijpF02 <- as.data.frame(RijpF02)
RijpF3 <- as.data.frame(RijpF3)
RijpNP <- as.data.frame(RijpNP)

# iNEXT

q123_vMijenDeep <- iNEXT(vMijenDeep, q = c(0,1,2), datatype = "abundance", endpoint = 120000)
q123_vMijenF02 <- iNEXT(vMijenF02, q = c(0,1,2), datatype = "abundance", endpoint = 600000)
q123_vMijenF3 <- iNEXT(vMijenF3, q = c(0,1,2), datatype = "abundance", endpoint = 700000)
q123_vMijenNP <- iNEXT(vMijenNP, q = c(0,1,2), datatype = "abundance", endpoint = 250000)

q123_KongsDeep <- iNEXT(KongsDeep, q = c(0,1,2), datatype = "abundance", endpoint = 160000)
q123_KongsF02 <- iNEXT(KongsF02, q = c(0,1,2), datatype = "abundance", endpoint = 800000)      ##
q123_KongsF3 <- iNEXT(KongsF3, q = c(0,1,2), datatype = "abundance", endpoint = 1200000)        ##
q123_KongsNP <- iNEXT(KongsNP, q = c(0,1,2), datatype = "abundance", endpoint = 170000)

q123_WijdeDeep <- iNEXT(WijdeDeep, q = c(0,1,2), datatype = "abundance", endpoint = 80000)
q123_WijdeF02 <- iNEXT(WijdeF02, q = c(0,1,2), datatype = "abundance", endpoint = 300000)
q123_WijdeF3 <- iNEXT(WijdeF3, q = c(0,1,2), datatype = "abundance", endpoint = 17500000)   ##
q123_WijdeNP <- iNEXT(WijdeNP, q = c(0,1,2), datatype = "abundance", endpoint = 125000)

q123_RijpDeep <- iNEXT(RijpDeep, q = c(0,1,2), datatype = "abundance", endpoint = 100000)
q123_RijpF02 <- iNEXT(RijpF02, q = c(0,1,2), datatype = "abundance", endpoint = 300000)
q123_RijpF3 <- iNEXT(RijpF3, q = c(0,1,2), datatype = "abundance", endpoint = 17000000)
q123_RijpNP <- iNEXT(RijpNP, q = c(0,1,2), datatype = "abundance", endpoint = 120000)


## Plots

# vMijen
p_vMDeep <- ggiNEXT(q123_vMijenDeep, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#a0d0d7", "#7eb1d5", "#5b92d4", "#3973d2")) + 
   labs(title = "Hill numbers - van Mijenfjorden - NextSeq - 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "Deep DNA \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq2_18S_vMijen_deep.jpeg", p_vMDeep, path = wd, width = 15, height = 15, dpi = 600)

p_vM02 <- ggiNEXT(q123_vMijenF02, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#67379a", "#bf3eff", "#6046ab", "#9A32CD", "#5a54bb", "#9932CC", 
                                                                                                                                   "#00008b", "#506ad5", "#0000FF", "#4a78e5", "#0000ee", "#4387f6")) + 
  labs(title = "Hill numbers - van Mijenfjorden - NextSeq - 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "F02 filters \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq2_18S_vMijen_F02.jpeg", p_vM02, path = wd, width = 15, height = 15, dpi = 600)

p_vMF3 <- ggiNEXT(q123_vMijenF3, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#da1365", "#cd2572", "#c0377f", "#b24a8d", "#a55c9a", "#986ea7", 
                                                                                                                                  "#8b80b4", "#7e92c1", "#71a4ce", "#63b7dc", "#56c9e9", "#49dbf6")) + 
  labs(title = "Hill numbers - van Mijenfjorden - NextSeq - 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "F3 filters \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq2_18S_vMijen_F3_2.png", p_vMF3, path = wd, width = 15, height = 15, dpi = 600)

p_vMNP <- ggiNEXT(q123_vMijenNP, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#1f3495", "#413c90", "#58458b", "#6a4e86", "#795781", "#88617b", "#966b76", "#a37670", 
                                                                                                                                  "#af8069", "#bc8b62", "#c8965a", "#d3a151", "#dfac47", "#ebb83a", "#f6c329")) + 
  labs(title = "Hill numbers - van Mijenfjorden - NextSeq - 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "Net and pump samples \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq2_18S_vMijen_NP.jpeg", p_vMNP, path = wd, width = 15, height = 15, dpi = 600)


# Kongs
p_KoDeep <- ggiNEXT(q123_KongsDeep, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#e184ff", "#eeab7b", "#f6c329")) + 
  labs(title = "Hill numbers - Kongsfjorden - NextSeq - 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "Deep DNA \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq2_18S_Kongs_deep.jpeg", p_KoDeep, path = wd, width = 15, height = 15, dpi = 600)

p_Ko02 <- ggiNEXT(q123_KongsF02, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#d47cf1", "#5D478B", "#AB82FF", "#9ca5c4", "#7CCD7C", "#76c0a7",
                                                                                                                                  "#64ce98", "#90EE90", "#3ee97a")) + 
  labs(title = "Hill numbers - Kongsfjorden - NextSeq - 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "F02 filters \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq2_18S_Kongs_F02.jpeg", p_Ko02, path = wd, width = 15, height = 15, dpi = 600)

p_KoF3 <- ggiNEXT(q123_KongsF3, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#48e4be", "#44cca7", "#40b391", "#3c9b7a", "#388363",
                                                                                                                                 "#336a4c", "#2f5236", "#2b391f", "#272108")) + 
  labs(title = "Hill numbers - Kongsfjorden - NextSeq - 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "F3 filters \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq2_18S_Kongs_F3.jpeg", p_KoF3, path = wd, width = 15, height = 15, dpi = 600)

p_KoNP <- ggiNEXT(q123_KongsNP, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#7245b9", "#BA55D3", "#8a3897", "#8968CD", "#a22c76", "#8B4789", "#ae2665",
                                                                                                                                 "#CD00CD", "#c61943", "#FF00FF", "#CD69C9", "#f60000")) + 
  labs(title = "Hill numbers - Kongsfjorden - NextSeq - 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "Net and pump samples \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq2_18S_Kongs_NP.jpeg", p_KoNP, path = wd, width = 15, height = 15, dpi = 600)


# wijde

p_WijdeDeep <- ggiNEXT(q123_WijdeDeep, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#9d60ff", "#ca3ab0", "#f61261")) + 
  labs(title = "Hill numbers - Wijdefjorden - NextSeq - 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "Deep DNA \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq2_18S_Wijde_deep.jpeg", p_WijdeDeep, path = wd, width = 15, height = 15, dpi = 600)

p_Wijde02 <- ggiNEXT(q123_WijdeF02, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#2528ff", "#6346eb", "#8061d6", "#927ac1", "#9e93ab", 
                                                                                                                                     "#a5ac94", "#a9c47b", "#a8dd5d", "#a3f630")) + 
  labs(title = "Hill numbers - Wijdefjorden - NextSeq - 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "F02 filters \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq2_18S_Wijde_F02.jpeg", p_Wijde02, path = wd, width = 15, height = 15, dpi = 600)

p_WijdeF3 <- ggiNEXT(q123_WijdeF3, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#1a50a7", "#2866aa", "#367cac", "#4492af", "#53a8b2", "#61bdb4",
                                                                                                                                    "#6fd3b7", "#7de9b9", "#8bffbc")) + 
  labs(title = "Hill numbers - Wijdefjorden - NextSeq - 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "F3 filters \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq2_18S_Wijde_F3.jpeg", p_WijdeF3, path = wd, width = 15, height = 15, dpi = 600)

p_WijdeNP <- ggiNEXT(q123_WijdeNP, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#58ff82", "#56f191", "#55e394", "#53d598", "#51c79b", "#4fb99e", "#4eaca1",
                                                                                                                                    "#4c9ea4", "#4a90a7", "#4882ab", "#4774ae", "#4566b1", "#4358b4")) + 
  labs(title = "Hill numbers - Wijdefjorden - NextSeq - 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "Net and pump samples \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq2_18S_Wijde_NP.jpeg", p_WijdeNP, path = wd, width = 15, height = 15, dpi = 600)


# rijp
p_RijpDeep <- ggiNEXT(q123_RijpDeep, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#58ff8e", "#4eaca1", "#4358b4")) + 
  labs(title = "Hill numbers - Rijpfjorden - NextSeq - 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "Deep DNA \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq2_18S_Rijp_deep.jpeg", p_RijpDeep, path = wd, width = 15, height = 15, dpi = 600)

p_Rijp02 <- ggiNEXT(q123_RijpF02, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#ff20fb", "#d92bed", "#b436df", "#8e42d0", "#CD2990", "#4358b4")) + 
  labs(title = "Hill numbers - Rijpfjorde - NextSeq - 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "F02 filters \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq2_18S_Rijp_F02.jpeg", p_Rijp02, path = wd, width = 15, height = 15, dpi = 600)

p_RijpF3 <- ggiNEXT(q123_RijpF3, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#52ff23", "#70d15b", "#77a579", "#72798f", "#5f4ca3", "#3417b4")) + 
  labs(title = "Hill numbers - Rijpfjorden - NextSeq - 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "F3 filters \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq2_18S_Rijp_F3.jpeg", p_RijpF3, path = wd, width = 15, height = 15, dpi = 600)

p_RijpNP <- ggiNEXT(q123_RijpNP, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#611ad2", "#AB82FF", "#8238d3", "#5D478B", "#a255d5", "#FF34B3",
                                                                                                                                  "#c373d6", "#CD2990")) + 
  labs(title = "Hill numbers - Rijpfjorden - NextSeq - 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "Net and pump samples \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )

ggsave("NextSeq2_18S_Rijp_NP.jpeg", p_RijpNP, path = wd, width = 15, height = 15, dpi = 600)






