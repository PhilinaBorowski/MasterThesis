## assign taxa NextSeq merged

library(dada2)
library(ShortRead)
library(Biostrings)
library(stringr)
library(R.utils)
library(plyr)


seqtab <- read.csv("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/NextSeq_18S_merged_seqtab_all.csv", row.names = 1)


seqtab_t <- t(seqtab)
db_dir <- file.path("/isibhv/projects/p_bioinf2/dbs/dada2_dbs")
taxDB <- "pr2_version_5.0.0_SSU_dada2.fasta"
taxLvs <- c("Kingdom", "Supergroup","Phylum","Division", "Class", "Order", "Family", "Genus", "Species")

taxa <- assignTaxonomy(seqtab_t, file.path(db_dir,taxDB), taxLevels = taxLvs, multithread=TRUE, minBoot = 70)

write.csv(taxa, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/NextSeq_18S_merged_taxa.csv")

taxa <- read.csv("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/NextSeq_18S_merged_taxa.csv")

bind <- cbind((seqtab), taxa)

write.csv(bind, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/NextSeq2_18S_merged_assigned_all.csv", quote=FALSE ) # CHANGE ME to output directory.

