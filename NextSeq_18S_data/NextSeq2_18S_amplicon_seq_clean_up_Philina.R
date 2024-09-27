## NextSeq II - 18S - Amplicon clean-up

### set wd ###
setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2")

### packages

library(dplyr)

### read in files
X18S_taxonomy <- read.csv("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/NextSeq2_18S_merged_assigned_all.csv", row.names = 1)

str(X18S_taxonomy)


###
#Remove metazoans

## first check how many in Excel and which column 
## 685 ASVs as Metazoa / in Division

metazoa <- X18S_taxonomy %>% dplyr::filter(Division == "Metazoa") 
str(metazoa) # 685 obs. of 142 var

other18S <- X18S_taxonomy%>%dplyr::filter(!Division %in% metazoa$Division)

# 16423 ASVs before 
16423 - 685 
# 15738

write.csv(other18S, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/NextSeq2_18S_wo_metazoa.csv", quote=FALSE )
