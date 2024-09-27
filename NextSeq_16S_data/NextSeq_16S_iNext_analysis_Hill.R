## Hill numbers NextSeq 16S
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
wd <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/iNext_analysis")


N16S <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/iNext_analysis/NextSeq_16S_ampli_cleaned.xlsx", sheet = 1)
str(N16S)
N16S[1:36] <- lapply(N16S[1:36], as.numeric)
str(N16S)

van_Mijen <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/iNext_analysis/NextSeq_16S_ampli_cleaned.xlsx", sheet = "vMijenfjorden")
kongs <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/iNext_analysis/NextSeq_16S_ampli_cleaned.xlsx", sheet = "Kongsfjorden")
wijde <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/iNext_analysis/NextSeq_16S_ampli_cleaned.xlsx", sheet = "Wijdefjorden")
rijp <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/iNext_analysis/NextSeq_16S_ampli_cleaned.xlsx", sheet = "Rijpfjorden")

str(van_Mijen)
van_Mijen[1:12] <- lapply(van_Mijen[1:12], as.numeric)
kongs[1:9] <- lapply(kongs[1:9], as.numeric)
wijde[1:9] <- lapply(wijde[1:9], as.numeric)
rijp[1:6] <- lapply(rijp[1:6], as.numeric)
12+9+9+6
#36

str(van_Mijen)

van_Mijen <- as.data.frame(van_Mijen)
kongs <- as.data.frame(kongs)
wijde <- as.data.frame(wijde)
rijp <- as.data.frame(rijp)

## get Hill numbers
## & save in table

## van Mijenfjorden 

hill_van_Mijen <- iNEXT(van_Mijen, q = 0, datatype = "abundance", se = T, conf = 0.95)
#SC = the estimated sample coverage for a sample of size m.
#m = sample size for which diversity estimates of order q are computed; by default setting (in the left hand side of the screen), m represents the sample size for each of the 40 knots between 1 and the default endpoint (double the reference sample size). Under 'General Settings', you can also either specify the endpoint and number of knots or specify the samples sizes for which you would like to calculate diversity estimates.
#Method = Rarefaction, Observed, or Extrapolation, depending on whether the size m is less than, equal to, or greater than the reference sample size.
#Order.q = the diversity order of q you selected in the 'General Settings' on the left hand side of the screen.
#qD = the estimated diversity of order q for a sample of size m.
#qD.LCL, qD.UCL = the bootstrap lower and upper confidence limits for the diversity of order q at the specified level in the settings (with a default value of 0.95).

head(hill_van_Mijen$DataInfo)   # summarizing data information
head(hill_van_Mijen$iNextEst)   # showing diversity estimates along with related statistics for a series of rarefied and extrapolated samples 
head(hill_van_Mijen$AsyEst)     # showing asymptotic diversity estimates along with related statistics

cleaned_vMijen <- hill_van_Mijen$AsyEst %>% select(-one_of("s.e.", "LCL", "UCL", "Estimator")) 
head(cleaned_vMijen)
wide_van_Mijen <- spread(cleaned_vMijen, key = Diversity, value = Observed)
head(wide_van_Mijen)

se_vMijen <- hill_van_Mijen$AsyEst %>% select(-one_of("LCL", "UCL", "Estimator", "Observed"))
head(se_vMijen)
wide_se_vMijen <- spread(se_vMijen, key = Diversity, value = s.e.)
head(wide_se_vMijen)
colnames(wide_se_vMijen) <- c("Assemblage","Shannon_se","Simpson_se", "Richness_se")
head(wide_se_vMijen)

wide_van_Mijen_se <- merge(wide_van_Mijen, wide_se_vMijen, all = F)
head(wide_van_Mijen_se)

wide_van_Mijen_se <- wide_van_Mijen_se[, c(1,2,5,3,6,4,7)]

write.csv(hill_van_Mijen$AsyEst, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/iNext_analysis/NextSeq_16S_vMijen_all_Hill.csv")
write.csv(wide_van_Mijen_se, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/iNext_analysis/NextSeq_16S_vMijen_Hill_short.csv")


## Kongsfjorden

hill_kongs <- iNEXT(kongs, q = 0, datatype = "abundance", se = T, conf = 0.95)
head(hill_kongs$AsyEst)

cleaned_kongs <- hill_kongs$AsyEst %>% select(-one_of("s.e.", "LCL", "UCL", "Estimator")) 
head(cleaned_kongs)
wide_kongs <- spread(cleaned_kongs, key = Diversity, value = Observed)
head(wide_kongs)

se_kongs <- hill_kongs$AsyEst %>% select(-one_of("LCL", "UCL", "Estimator", "Observed"))
head(se_kongs)
wide_se_kongs <- spread(se_kongs, key = Diversity, value = s.e.)
head(wide_se_kongs)
colnames(wide_se_kongs) <- c("Assemblage","Shannon_se","Simpson_se", "Richness_se")
head(wide_se_kongs)

wide_kongs_se <- merge(wide_kongs, wide_se_kongs, all = F)
head(wide_kongs_se)

wide_kongs_se <- wide_kongs_se[, c(1,2,5,3,6,4,7)]

write.csv(hill_kongs$AsyEst, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/iNext_analysis/NextSeq_16S_kongs_Hill.csv")
write.csv(wide_kongs_se, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/iNext_analysis/NextSeq_16S_kongs_Hill_short.csv")


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


write.csv(hill_wijde$AsyEst, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/iNext_analysis/NextSeq_16S_wijde_Hill.csv")
write.csv(wide_wijde_se, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/iNext_analysis/NextSeq_16S_wijde_Hill_short.csv")


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


write.csv(hill_rijp$AsyEst, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/iNext_analysis/NextSeq_16S_rijp_Hill.csv")
write.csv(wide_rijp_se, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/iNext_analysis/NextSeq_16S_rijp_Hill_short.csv")

