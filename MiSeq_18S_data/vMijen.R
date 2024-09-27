## MiSeq iNEXT - van Mijen

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
vMijenNP <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_vMijen.xlsx", sheet = "vMijenNP")

str(vMijenF02)

vMijenF02 <- as.data.frame(vMijenF02)
vMijenF3 <- as.data.frame(vMijenF3)
vMijenNP <- as.data.frame(vMijenNP)


q123_vMijenF02 <- iNEXT(vMijenF02, q = c(0,1,2), datatype = "abundance", endpoint = 500000)
q123_vMijenF3 <- iNEXT(vMijenF3, q = c(0,1,2), datatype = "abundance", endpoint = 1200000)
q123_vMijenNP <- iNEXT(vMijenNP, q = c(0,1,2), datatype = "abundance", endpoint = 500000)


head(q123_vMijenF02)
head(q123_vMijenF3)
head(q123_vMijenNP)



p_F02 <- ggiNEXT(q123_vMijenF02, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#b3e0a6", "#a1d991", "#92D282", "#86ca78", "#7bc16e",
                                                                                                                          "#71b966", "#68b05d", "#60a855", "#559f52", "#4b974f",
                                                                                                                          "#418e4d", "#358747")) +
     labs(title = "Hill numbers - van Mijenfjord", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "F02 filters \n q0 = Species Richness \n q1 = Shannon diversiry \n q2 = Simpson diversity")

ggsave("MiSeq_vMijen_F02.jpeg", p_F02, path = wd, width = 30, height = 30, dpi = 600)


p_F3 <- ggiNEXT(q123_vMijenF3, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#aca4e2", "#a0a7e2", "#93abe1", "#86aedf", "#78b1dc", "#6ab4d8",
                                                                                                                        "#5cb7d3", "#4fb9cd", "#3fbcc3", "#38bdbc", "#38beb4", "#3dbeac")) + 
     labs(title = "Hill numbers - van Mijenfjord", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "F3 filters \n q0 = Species Richness \n q1 = Shannon diversiry \n q2 = Simpson diversity")

ggsave("MiSeq_vMijen_F3.jpeg", p_F3, path = wd, width = 30, height = 30, dpi = 600)

p_NP <- ggiNEXT(q123_vMijenF3, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#abb065", "#b3ae64", "#bbab66", "#c3a869", "#caa66e", "#d0a374", "#d5a07b"),
                                                                                                                        "#dc9c87", "#df9a90", "#e29899", "#e396a2", "#e494aa", "#e494b3") + 
     labs(title = "Hill numbers - van Mijenfjord", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "Pump and Net samples \n q0 = Species Richness \n q1 = Shannon diversiry \n q2 = Simpson diversity")

ggsave("MiSeq_vMijen_NP.jpeg", p_NP, path = wd, width = 30, height = 30, dpi = 600)








