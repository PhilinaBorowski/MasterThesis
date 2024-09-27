## NextSeq1 18S - 11 samples - assign taxonomy


library(dada2)
library(ShortRead)
library(Biostrings)
library(stringr)
library(R.utils)
library(plyr)


### set wd ###
setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq/analysis_of_the_other_11_samples")

#seqtab.nochim <- read.csv("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq/analysis_of_the_other_11_samples/NextSeq1_18S_11samples_seqtab_all.csv", row.names = 1)

seqtab.nochim_t <- t(seqtab.nochim)
db_dir <- file.path("/isibhv/projects/p_bioinf2/dbs/dada2_dbs")
taxDB <- "pr2_version_5.0.0_SSU_dada2.fasta"
taxLvs <- c("Kingdom", "Supergroup","Phylum","Division", "Class", "Order", "Family", "Genus", "Species")

taxa <- assignTaxonomy(seqtab.nochim_t, file.path(db_dir,taxDB), taxLevels = taxLvs ,multithread=TRUE, minBoot=70)

write.csv(taxa, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq/analysis_of_the_other_11_samples/NextSeq1_18S_11samples_taxa.csv")

#taxa <- read.csv("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq/analysis_of_the_other_11_samples/NextSeq1_18S_11samples_taxa.csv")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
bind <- cbind(seqtab.nochim, taxa)

write.csv(bind, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq/analysis_of_the_other_11_samples/NextSeq1_18S_11samples_assigned_all_18S.csv", quote=FALSE ) # CHANGE ME to output directory.
