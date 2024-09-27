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


#Mi18S <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_wo_metazoa_wo_taxa_seq_w_sheets.xlsx", sheet = 1)
#str(Mi18S)
#Mi18S[1:96] <- lapply(Mi18S[1:96], as.numeric)
#str(Mi18S)

van_Mijen <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_metabarcoding_wo_metazoa.xlsx", sheet = "vMijenfjorden")
kongs <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_metabarcoding_wo_metazoa.xlsx", sheet = "Kongsfjorden")
wijde <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_metabarcoding_wo_metazoa.xlsx", sheet = "Wijdefjorden")
rijp <- read_xlsx("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_metabarcoding_wo_metazoa.xlsx", sheet = "Rijpfjorden")

str(van_Mijen)
van_Mijen[1:37] <- lapply(van_Mijen[1:37], as.numeric)
kongs[1:18] <- lapply(kongs[1:18], as.numeric)
wijde[1:27] <- lapply(wijde[1:27], as.numeric)
rijp[1:14] <- lapply(rijp[1:14], as.numeric)
37+18+27+14
#96

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
# size? endpoint? knots? nboot? 
# knots: user may also specify the number of knots in the range of sample size between 1 and the endpoint. If you choose a large number of knots, 
#then it may take a long time to obtain the output due to the time‐consuming bootstrap method

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

write.csv(hill_van_Mijen$AsyEst, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_van_Mijenf_all_Hill.csv")
write.csv(wide_van_Mijen_se, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_van_Mijenf_Hill_short.csv")

## can we also add the s.e. ?????????? somehow????????????
## pls???

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

write.csv(hill_kongs$AsyEst, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_kongs_Hill.csv")
write.csv(wide_kongs_se, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/iNext_analysis/MiSeq_18S_kongs_Hill_short.csv")


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



#?ggiNEXT()
# sample-size-based rarefaction/extrapolation curve (type = 1) 
# sample completeness curve (type = 2)
# coverage-based rarefaction/extrapolation curve (type = 3)
# se = T or F, to see intrapolation of curves

#q123_hill_van_Mijen <- iNEXT(van_Mijen, q = c(0,1,2), datatype = "abundance") ##

#p1 <- ggiNEXT(hill_van_Mijen, type = 1, se = TRUE, facet.var="Assemblage") + scale_colour_manual(values = c("#8A2BE2", "#6495ED", "#00008B", "#008b8b", "#a9a9a9",
#                                                                                                 "#006400", "#8b008b", "#8b0000", "#e9967a", "#483d8b",
#                                                                                                 "#2f4f4f", "#00ced1", "#228b22", "#adff2f", "#f08080",
#                                                                                                 "#90ee90", "#20b2aa", "#8470ff", "#778899", "#66cdaa",
#                                                                                                 "#0000cd", "#48d1cc", "#7b68ee", "#00fa9a", "#c71585",
#                                                                                                 "#ffe4b5", "#cd853f", "#8b4513", "#ee82ee", "#d02090",
#                                                                                                 "#36628d", "#00cd66", "#ffe1ff", "#ff6347", "#ff46ff",
#                                                                                                 "#ee9a49", "#8b8989"))


#p1_wrap <- p1 + facet_wrap(~Assemblage, ncol = 5)
#ggsave("plot1.jpeg", p1, path = wd, width = 30, height = 30, dpi = 600)
#ggsave("plot1_wrap.jpeg", p1_wrap, path = wd, width = 30, height = 30, dpi = 600)
#?ggsave()

## maybe try to change the length of x axis? down to 1500000

#q123_hill_van_Mijen <-  iNEXT(van_Mijen, q = c(0,1,2), datatype = "abundance", se = T, conf = 0.95)
#head(q123_hill_van_Mijen)
#p2 <- ggiNEXT(q123_hill_van_Mijen, type = 1, facet.var="Order", color.var = "Assemblage") + scale_colour_manual(values = c("#8A2BE2", "#6495ED", "#00008B", "#00008B", "#008b8b", "#a9a9a9",
#                                                                                                                     "#006400", "#8b008b", "#8b0000", "#e9967a", "#483d8b",
#                                                                                                                     "#2f4f4f", "#00ced1", "#228b22", "#adff2f", "#f08080",
#                                                                                                                     "#90ee90", "#20b2aa", "#8470ff", "#778899", "#66cdaa",
#                                                                                                                     "#0000cd", "#48d1cc", "#7b68ee", "#00fa9a", "#c71585",
#                                                                                                                     "#ffe4b5", "#cd853f", "#8b4513", "#ee82ee", "#d02090",
#                                                                                                                     "#36628d", "#00cd66", "#ffe1ff", "#ff6347", "#ff46ff",
#                                                                                                                     "#ee9a49", "#8b8989")) +
#              ggtitle("Hill numbers - van Mijenfjord") + xlab("Number of individuals") + ylab("Species diversity") #+ scale_x_continuous(limits = c(0, 15000)) + scale_y_continuous(limits = c(0, 400))


#p2_wrap <- p2 + facet_wrap(~Assemblage, ncol = 5)
#ggsave("plot2_wrap.jpeg", p2_wrap, path = wd, width = 30, height = 30, dpi = 600)
#ggsave("plot2_wrap_closer.jpeg", p2_wrap, path = wd, width = 30, height = 30, dpi = 600)
#ggsave("plot_q123.jpeg", p2, path = wd, width = 30, height = 30, dpi = 600)



## statistical tests?
## t-test?
## aov? anova? manova? manova for multiple outcome variables












