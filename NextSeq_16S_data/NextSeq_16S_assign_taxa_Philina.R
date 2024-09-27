## NextSeq 2 - Pro - assign taxa

library(dada2)
library(ShortRead)
library(Biostrings)
library(stringr)
library(R.utils)
library(plyr)


## set wd
setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S")

seqtab.nochim <- read.csv("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/NextSeq_seqtab_all_16S.csv", row.names = 1)

seqtab.nochim_t <- t(seqtab.nochim)
#ASSIGN TAXONOMY
# for 16S use Silva !!
db_dir <- file.path("/isibhv/projects/p_bioinf2/dbs/dada2_dbs")
taxDB <- "silva_nr99_v138.1_wSpecies_train_set_Pros.fa"
taxLvs <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

taxa <- assignTaxonomy(seqtab.nochim_t, file.path(db_dir,taxDB), taxLevels = taxLvs ,multithread=TRUE, minBoot=50)

write.csv(taxa, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/NextSeq2_16S_taxa.csv")
taxa <- read.csv("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/NextSeq2_16S_taxa.csv")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
bind <- cbind((seqtab.nochim), taxa)

write.csv(bind, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/NextSeq_16S_assigned_all_Philina.csv", quote=FALSE ) # CHANGE ME to output directory.

