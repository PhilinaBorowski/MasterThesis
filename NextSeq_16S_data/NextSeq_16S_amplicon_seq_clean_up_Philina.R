## NextSeq - 16S - amplicon clean-up
#

### set wd ###
setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S")

### packages

library(dplyr)

### read in files
## change names
X16S_ASVTable <- read.csv("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/NextSeq_seqtab_all_16S.csv", row.names = 1)
X16S_taxonomy <- read.csv("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/NextSeq_16S_assigned_all_Philina.csv", row.names = 1)

str(X16S_taxonomy) # 2501 obs
###

# Remove mitochondria + chloroplasts
# check in Excel before! 
# how many chloroplasts?     # 271   in Order
# how many mitochondroia?    # 186   in Family


chloro <- X16S_taxonomy %>% dplyr::filter(Order == "Chloroplast") 
str(chloro) # 271 obs. of 43

wo_chloro_16 <- X16S_taxonomy%>%dplyr::filter(!Order %in% chloro$Order)
# 2230

2501 - 271
#2230

mitoch <- X16S_taxonomy %>% dplyr::filter(Family == "Mitochondria") 
str(mitoch) # 186 obs. of 43

wo_chloro_mitoch <- wo_chloro_16 %>% dplyr::filter(!Family %in% mitoch$Family)

2501 - 271 - 186
# 2044

write.csv(wo_chloro_mitoch, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/NextSeq_16S_ampli_cleaned.csv", quote=FALSE )

