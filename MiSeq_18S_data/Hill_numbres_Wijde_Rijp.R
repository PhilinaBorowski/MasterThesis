## take a look at euk MiSeq data
##
#
#install.packages("readxl")

library(httr)
library(iNEXT)
library(ggplot2)
library(readxl)
library(tidyverse)
library(dplyr)

## set wd
wd <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/")


Mi18S <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_wo_metazoa_wo_taxa_seq_w_sheets.xlsx", sheet = 1)
str(Mi18S)
Mi18S[1:96] <- lapply(Mi18S[1:96], as.numeric)
str(Mi18S)

#van_Mijen <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_wo_metazoa_wo_taxa_seq_w_sheets.xlsx", sheet = "van_Mijenfjorden")
#kongs <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_wo_metazoa_wo_taxa_seq_w_sheets.xlsx", sheet = "Kongsfjorden")
wijde <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_wo_metazoa_wo_taxa_seq_w_sheets.xlsx", sheet = "Wijdefjorden")
rijp <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_wo_metazoa_wo_taxa_seq_w_sheets.xlsx", sheet = "Rijpfjorden")

#str(van_Mijen)
#van_Mijen[1:37] <- lapply(van_Mijen[1:37], as.numeric)
#kongs[1:18] <- lapply(kongs[1:18], as.numeric)
wijde[1:27] <- lapply(wijde[1:27], as.numeric)
rijp[1:14] <- lapply(rijp[1:14], as.numeric)
#37+18+27+14
#96

#str(van_Mijen)

#van_Mijen <- as.data.frame(van_Mijen)
#kongs <- as.data.frame(kongs)
wijde <- as.data.frame(wijde)
rijp <- as.data.frame(rijp)



## Wijdefjorden 

hill_wijde <- iNEXT(wijde, q = 0, datatype = "abundance", se = T, conf = 0.95)
head(hill_wijde$AsyEst)

cleaned_wijde <- hill_wijde$AsyEst %>% select(-one_of("s.e.", "LCL", "UCL", "Estimator")) 
head(cleaned_wijde)
wide_wijde <- spread(cleaned_wijde, key = Diversity, value = Observed)
head(wide_wijde)

se_wijde <- hill_wijde$AsyEst %>% select(-one_of("LCL", "UCL", "Estimator", "Observed"))
head(se_wijde)
wide_se_wijde <- spread(se_wijde, key = Diversity, value = s.e.)
head(wide_se_wijde)
colnames(wide_se_wijde) <- c("Assemblage","Shannon_se","Simpson_se", "Richness_se")
head(wide_se_wijde)

wide_wijde_se <- merge(wide_wijde, wide_se_wijde, all = F)
head(wide_wijde_se)

wide_wijde_se <- wide_wijde_se[, c(1,2,5,3,6,4,7)]


write.csv(hill_wijde$AsyEst, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_wijde_Hill.csv")
write.csv(wide_wijde_se, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_wijde_Hill_short.csv")


## Rijpfjorden 

hill_rijp <- iNEXT(rijp, q = 0, datatype = "abundance", se = T, conf = 0.95)
head(hill_rijp$AsyEst)

cleaned_rijp <- hill_rijp$AsyEst %>% select(-one_of("s.e.", "LCL", "UCL", "Estimator")) 
head(cleaned_rijp)
wide_rijp <- spread(cleaned_rijp, key = Diversity, value = Observed)
head(wide_rijp)

se_rijp <- hill_rijp$AsyEst %>% select(-one_of("LCL", "UCL", "Estimator", "Observed"))
head(se_rijp)
wide_se_rijp <- spread(se_rijp, key = Diversity, value = s.e.)
head(wide_se_rijp)
colnames(wide_se_rijp) <- c("Assemblage","Shannon_se","Simpson_se", "Richness_se")
head(wide_se_rijp)

wide_rijp_se <- merge(wide_rijp, wide_se_rijp, all = F)
head(wide_rijp_se)

wide_rijp_se <- wide_rijp_se[, c(1,2,5,3,6,4,7)]


write.csv(hill_rijp$AsyEst, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_rijp_Hill.csv")
write.csv(wide_rijp_se, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_rijp_Hill_short.csv")
