## MiSeq iNEXT - Kongsfjorden

library(httr)
library(iNEXT)
library(ggplot2)
library(readxl)
library(tidyverse)
library(dplyr)

## set wd
wd <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/")

kongs <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_metabarcoding_wo_metazoa.xlsx", sheet = "Kongsfjorden")

#str(kongs)
#kongs[1:18] <- lapply(kongs[1:18], as.numeric)

#str(kongs)
#kongs <- as.data.frame(kongs)


#hill_kongs <- iNEXT(kongs, q = 0, datatype = "abundance", se = T, conf = 0.95)
#head(hill_kongs$AsyEst)

#cleaned_kongs <- hill_kongs$AsyEst %>% select(-one_of("s.e.", "LCL", "UCL", "Estimator")) 
#head(cleaned_kongs)
#wide_kongs <- spread(cleaned_kongs, key = Diversity, value = Observed)
#head(wide_kongs)

#se_kongs <- hill_kongs$AsyEst %>% select(-one_of("LCL", "UCL", "Estimator", "Observed"))
#head(se_kongs)
#wide_se_kongs <- spread(se_kongs, key = Diversity, value = s.e.)
#head(wide_se_kongs)
#colnames(wide_se_kongs) <- c("Assemblage","Shannon_se","Simpson_se", "Richness_se")
#head(wide_se_kongs)

#wide_kongs_se <- merge(wide_kongs, wide_se_kongs, all = F)
#head(wide_kongs_se)

#wide_kongs_se <- wide_kongs_se[, c(1,2,5,3,6,4,7)]

#write.csv(hill_kongs$AsyEst, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_kongs_Hill.csv")
#write.csv(wide_kongs_se, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_kongs_Hill_short.csv")


## to se wether it makes a different in Hill if q=0 or q=012

#hill_kongs_123 <- iNEXT(kongs, q = c(0,1,2), datatype = "abundance", se = T, conf = 0.95)
#head(hill_kongs_123$AsyEst)

#cleaned_kongs_123 <- hill_kongs_123$AsyEst %>% select(-one_of("s.e.", "LCL", "UCL", "Estimator")) 
#head(cleaned_kongs_123)
#wide_kongs_123 <- spread(cleaned_kongs_123, key = Diversity, value = Observed)
#head(wide_kongs_123)

#se_kongs_123 <- hill_kongs_123$AsyEst %>% select(-one_of("LCL", "UCL", "Estimator", "Observed"))
#head(se_kongs_123)
#wide_se_kongs_123 <- spread(se_kongs_123, key = Diversity, value = s.e.)
#head(wide_se_kongs_123)
#colnames(wide_se_kongs_123) <- c("Assemblage","Shannon_se","Simpson_se", "Richness_se")
#head(wide_se_kongs_123)
#
#wide_kongs_se_123 <- merge(wide_kongs_123, wide_se_kongs_123, all = F)
#head(wide_kongs_se_123)
#
#wide_kongs_se_123 <- wide_kongs_se_123[, c(1,2,5,3,6,4,7)]
#
#write.csv(hill_kongs_123$AsyEst, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_kongs_Hill123.csv")
#write.csv(wide_kongs_se_123, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_kongs_Hill_short123.csv")


# Plots

KongsF02 <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_Kongs.xlsx", sheet = "Kongs_F02")
KongsF3 <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_Kongs.xlsx", sheet = "Kongs_F3")


KongsF02 <- as.data.frame(KongsF02)
KongsF3 <- as.data.frame(KongsF3)


q123_KongsF02 <- iNEXT(KongsF02, q = c(0,1,2), datatype = "abundance", endpoint = 1600000)
q123_KongsF3 <- iNEXT(KongsF3, q = c(0,1,2), datatype = "abundance", endpoint = 60000)



p_F02 <- ggiNEXT(q123_KongsF02, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#93ff20", "#84ff3b", "#76ff55", "#67ff70", "#59ff8a", 
                                                                                                                                 "#4affa5", "#3bffbf", "#2dffda", "#1efff4")) + 
  labs(title = "Hill numbers - Kongsfjorden - MiSeq 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "F02 filters \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )


ggsave("MiSeq_Kongs_F02_2.png", p_F02, path = wd, width = 15, height = 15, dpi = 600)


p_F3 <- ggiNEXT(q123_KongsF3, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#ff84fd", "#eb86fd", "#d888fe", "#c48afe", "#b18cfe",
                                                                                                                               "#9d8dfe", "#898fff", "#7691ff", "#6293ff")) + 
  labs(title = "Hill numbers - Kongsfjorden - MiSeq 18S", xlab = "Number of individuals", ylab = "Species diversity", subtitle = "F3 filters \n q0 = Species Richness \n q1 = Shannon diversity \n q2 = Simpson diversity") +
  theme_light() +  theme(axis.text=element_text(size=14), #change font size of axis text
                         axis.title=element_text(size=19), #change font size of axis titles
                         plot.title=element_text(size=26), #change font size of plot title
                         plot.subtitle = element_text(size=22), 
                         legend.text=element_text(size=16), #change font size of legend text
                         legend.title=element_text(size=20)) #change font size of legend title   )


ggsave("MiSeq_Kongs_F3.png", p_F3, path = wd, width = 15, height = 15, dpi = 600)




