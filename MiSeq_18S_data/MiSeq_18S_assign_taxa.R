## MiSeq1 18S - assign taxonomy


library(dada2)
library(ShortRead)
library(Biostrings)
library(stringr)
library(R.utils)
library(plyr)


### set wd ###
setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/MiSeq_results/")

seqtab.nochim <- read.csv("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/MiSeq_results/MiSeq_18S_seqtab_all.csv", row.names = 1)

#seqtab.nochim_t <- t(seqtab.nochim)
#db_dir <- file.path("/isibhv/projects/p_bioinf2/dbs/dada2_dbs")
#taxDB <- "pr2_version_5.0.0_SSU_dada2.fasta"
#taxLvs <- c("Kingdom", "Supergroup","Phylum","Division", "Class", "Order", "Family", "Genus", "Species")

#taxa <- assignTaxonomy(seqtab.nochim_t, file.path(db_dir,taxDB), taxLevels = taxLvs ,multithread=TRUE, minBoot=70)

#write.csv(taxa, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/MiSeq_results/MiSeq_18S_taxa.csv")

taxa <- read.csv("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/MiSeq_results/MiSeq_18S_taxa.csv")

bind <- cbind(seqtab.nochim, taxa)

write.csv(bind, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/MiSeq_results/MiSeq_18S_assigned_all.csv", quote=FALSE ) # CHANGE ME to output directory.
