## merging both NextSeq 18S datasets
##
## 

#install.packages("tidyverse")
#library(tidyverse)
library(dada2)

nextseq1 <- read.csv("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq/analysis_of_the_other_11_samples/NextSeq1_18S_11samples_assigned_all.csv", header = T, row.names = 1) # copied into NextSeq2 folder
#5278 ASVs
nextseq2 <- read.csv("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/NextSeq2_18S_assigned_all_Philina.csv", header = T, row.names = 1)
#13830 ASVs

nextseq1_t <- t(nextseq1)
nextseq2_t <- t(nextseq2)

str(nextseq1)

class(nextseq1); mode(nextseq1); dim(nextseq1)  # data.frame, list, 5278 21
class(nextseq2); mode(nextseq2); dim(nextseq2)  # data.frame, list, 1380 132

dada2:::is.sequence.table(nextseq1)
dada2:::is.sequence.table(nextseq2)

rownames(nextseq1)
colnames(nextseq1)
any(is.na(nextseq1))

nextseq1 <- as.matrix(nextseq1)

head(nextseq1)

merged <- mergeSequenceTables(table1 = nextseq1, table2 = nextseq2, repeats = error)
merged <- mergeSequenceTables(nextseq1_t, nextseq2_t, repeats = error)

write.csv(t(merged), "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/NextSeq_18S_merged_assigned_all.csv")
# 16424 ASVs in the end 