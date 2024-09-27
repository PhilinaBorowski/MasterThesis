### DADA2 NextSeq - Prokaryotes

library(dada2)
library(ShortRead)
library(Biostrings)
library(stringr)
library(R.utils)
library(plyr)


## set wd
setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S")
wd <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S")

raw_dir <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/FastQ_files_unzip")
raw_dir <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/FastQ_files_unzip")# run 2x
#View(raw_dir)
list.files(raw_dir)


###
### create folders for later analysis
###

preFilt_dir <- file.path(wd,"preFilt_16S") # vorgefilterte 
primerCut5_dir <- file.path(wd,"primerCut5_16S") # primer 5' cuts
primerCut3_dir <- file.path(wd,"primerCut3_16S") # primer 3' cuts
qualFiltTrim_dir <- file.path(wd,"qualFiltTrim_16S") # quality filter trim 


###
### "CONSTRUCT NEEDED FILE LISTS"
### fn = full names

fnFs.raw <- sort(list.files(raw_dir, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs.raw <- sort(list.files(raw_dir, pattern = "_R2_001.fastq", full.names = TRUE))
fnFs.preFilt <- file.path(preFilt_dir,basename(fnFs.raw))
fnRs.preFilt <- file.path(preFilt_dir,basename(fnRs.raw))
fnFs.primerCut5 <- file.path(primerCut5_dir,basename(fnFs.raw))
fnRs.primerCut5 <- file.path(primerCut5_dir,basename(fnRs.raw))
fnFs.primerCut3 <- file.path(primerCut3_dir,basename(fnFs.raw))
fnRs.primerCut3 <- file.path(primerCut3_dir,basename(fnRs.raw))
fnFs.qualFiltTrim <- file.path(qualFiltTrim_dir,basename(fnFs.raw))
fnRs.qualFiltTrim <- file.path(qualFiltTrim_dir,basename(fnRs.raw))


###
#GET proper SAMPLE NAMES; shorten names
basename(fnFs.raw)
sample.names <- str_remove(basename(fnFs.raw),"_R1_001.fastq")
#sample.names <- str_remove(basename(fnFs.raw), "HE627-")
head(sample.names)

basename(fnRs.raw)
basename(sample.names)


#PREFILTERING
### doesnt work on windows haha -> have to put "multithread = FALSE" instead of "TRUE"
filterAndTrim(fnFs.raw,fnFs.preFilt,fnRs.raw,fnRs.preFilt,truncQ=2,minQ=2,minLen=50,maxN=0,multithread = FALSE, compress = FALSE)
# filterAndTrim(input, output)
# truncQ = Truncate reads at the first instance of a quality score less than or equal to truncQ
# minQ = After truncation, reads contain a quality score less than minQ will be discarded.
# minLen = Remove reads with length less than minLen. minLen is enforced after trimming and truncation.
# maxN = After truncation, sequences with more than maxN Ns will be discarded; DADA2 requires N to be 0!

## no applicable method for 'depth' applied to an object of class "NULL"




### 
#IDENTIFY PRIMERs

FWD_PRIMER <- "GTGCCAGCMGCCGCGGTAA"                          # primer Sequence (5’ - 3')
REV_PRIMER <- "GGACTACHVGGGTWTCTAAT" 

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))
}


FWD_PRIMER.orients <- allOrients(FWD_PRIMER)
REV_PRIMER.orients <- allOrients(REV_PRIMER)

FWD_PRIMER.orients   # what is RevComp?
REV_PRIMER.orients


primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

#check one sample (1st sample - St1-F02-1)
rbind(FWD_PRIMER.ForwardReads = sapply(FWD_PRIMER.orients, primerHits, fn = fnFs.preFilt[[1]]),
      FWD_PRIMER.ReverseReads = sapply(FWD_PRIMER.orients, primerHits, fn = fnRs.preFilt[[1]]),
      REV_PRIMER.ForwardReads = sapply(REV_PRIMER.orients, primerHits, fn = fnFs.preFilt[[1]]),
      REV_PRIMER.ReverseReads = sapply(REV_PRIMER.orients, primerHits, fn = fnRs.preFilt[[1]]))


#                          Forward Complement Reverse RevComp
#FWD_PRIMER.ForwardReads  185363          0       0       7
#FWD_PRIMER.ReverseReads       9          0       0  167124
#REV_PRIMER.ForwardReads       8          0       0  204197
#REV_PRIMER.ReverseReads  206598          0       0       8



system2("cutadapt", args = "--version") # Run shell commands from R
# 2.8



### 
###
#create output dirs
if(!dir.exists(primerCut5_dir)) dir.create(primerCut5_dir)
if(!dir.exists(primerCut3_dir)) dir.create(primerCut3_dir)


FWD_PRIMER.RC <- dada2:::rc(FWD_PRIMER)
REV_PRIMER.RC <- dada2:::rc(REV_PRIMER)

# Run Cutadapt both primer at 5' end
# min_overlap primer bases -1
# reverse: 20 bases. so 19
# forward: 19, so 18

for(i in seq_along(fnFs.preFilt)) {
  system2("cutadapt", args = c("-g", paste("\"",FWD_PRIMER,";min_overlap=18;max_error_rate=0.15","\"",sep=""),
                             "-G", paste("\"",REV_PRIMER,";min_overlap=19;max_error_rate=0.15","\"",sep=""),
                             "--discard-untrimmed", "-o", fnFs.primerCut5[i], "-p", fnRs.primerCut5[i],
                             fnFs.preFilt[i], fnRs.preFilt[i]))
  
}


###
### check if it worked; should be 0, I guess
#check one sample - 1st sample - St1-F02-1
rbind(FWD_PRIMER.ForwardReads = sapply(FWD_PRIMER.orients, primerHits, fn = fnFs.primerCut5[[1]]),
      FWD_PRIMER.ReverseReads = sapply(FWD_PRIMER.orients, primerHits, fn = fnRs.primerCut5[[1]]),
      REV_PRIMER.ForwardReads = sapply(REV_PRIMER.orients, primerHits, fn = fnFs.primerCut5[[1]]),
      REV_PRIMER.ReverseReads = sapply(REV_PRIMER.orients, primerHits, fn = fnRs.primerCut5[[1]]))

#                         Forward Complement Reverse RevComp
#FWD_PRIMER.ForwardReads       1          0       0       0
#FWD_PRIMER.ReverseReads       0          0       0  167036
#REV_PRIMER.ForwardReads       0          0       0  202945
#REV_PRIMER.ReverseReads       2          0       0       0



###
# Run Cutadapt for primer at 3' end
# min_overlap primer bases -1
# reverse: 20 bases. so 19
# forward: 19, so 18
for(i in seq_along(fnFs.primerCut5)) {
  system2("cutadapt", args = c("-a", paste("\"",REV_PRIMER.RC,";min_overlap=19;max_error_rate=0.15","\"",sep=""),
                             "-A", paste("\"",FWD_PRIMER.RC,";min_overlap=18;max_error_rate=0.15","\"",sep=""),
                             "-o", fnFs.primerCut3[i], "-p", fnRs.primerCut3[i],
                             fnFs.primerCut5[i], fnRs.primerCut5[i]))
  
}


    
###
### check if it worked with 1st sample (should be 0)
rbind(FWD_PRIMER.ForwardReads = sapply(FWD_PRIMER.orients, primerHits, fn = fnFs.primerCut3[[1]]),
      FWD_PRIMER.ReverseReads = sapply(FWD_PRIMER.orients, primerHits, fn = fnRs.primerCut3[[1]]),
      REV_PRIMER.ForwardReads = sapply(REV_PRIMER.orients, primerHits, fn = fnFs.primerCut3[[1]]),
      REV_PRIMER.ReverseReads = sapply(REV_PRIMER.orients, primerHits, fn = fnRs.primerCut3[[1]]))

## results
#                         Forward Complement Reverse RevComp
#FWD_PRIMER.ForwardReads       1          0       0       0
#FWD_PRIMER.ReverseReads       0          0       0       0
#REV_PRIMER.ForwardReads       0          0       0       0
#REV_PRIMER.ReverseReads       0          0       0       0



### QUality Control
###
### Quality filtering

filterOut <- filterAndTrim(fnFs.primerCut3,fnFs.qualFiltTrim,
                           fnRs.primerCut3,fnRs.qualFiltTrim,
                           maxN=0,maxEE=c(2.7,2.2),
                           truncLen=c(240,160),         # 16 S !!!
                           verbose = TRUE, rm.phix = TRUE,
                           compress = TRUE, multithread = FALSE)

print(filterOut)

#                                    reads.in     reads.out
#HE627-pro-St1-F02-1_S1_R1_001.fastq     212886    194261
#HE627-pro-St1-F02-2_S2_R1_001.fastq      96210     88293
#HE627-pro-St1-F02-3_S3_R1_001.fastq     190817    174464
#HE627-pro-St10-F02-1_S25_R1_001.fastq   188152    173895
#HE627-pro-St10-F02-2_S26_R1_001.fastq   179254    166464
#HE627-pro-St10-F02-3_S27_R1_001.fastq   185633    168844
#HE627-pro-St11-F02-1_S28_R1_001.fastq   152416    138823
#HE627-pro-St11-F02-2_S29_R1_001.fastq   179922    167684
#HE627-pro-St11-F02-3_S30_R1_001.fastq   156757    146694
#HE627-pro-St12-F02-1_S31_R1_001.fastq   194472    181558
#HE627-pro-St12-F02-2_S32_R1_001.fastq   206788    194329
#HE627-pro-St12-F02-3_S33_R1_001.fastq   154597    142436
#HE627-pro-St13-F02-1_S34_R1_001.fastq   183987    170142
#HE627-pro-St13-F02-2_S35_R1_001.fastq   137097    126839
#HE627-pro-St13-F02-3_S36_R1_001.fastq   165505    152840
#HE627-pro-St2-F02-1_S4_R1_001.fastq     702424    644103
#HE627-pro-St2-F02-2_S5_R1_001.fastq     212979    196206
#HE627-pro-St2-F02-3_S6_R1_001.fastq     166042    153985
#HE627-pro-St3-F02-1_S7_R1_001.fastq     192832    178829
#HE627-pro-St3-F02-2_S8_R1_001.fastq     171730    157450
#HE627-pro-St3-F02-3_S9_R1_001.fastq     148715    137577
#HE627-pro-St4-F02-1_S10_R1_001.fastq    192479    177908
#HE627-pro-St4-F02-2_S11_R1_001.fastq    171200    156986
#HE627-pro-St4-F02-3_S12_R1_001.fastq    196101    179954
#HE627-pro-St5-F02-1_S13_R1_001.fastq    179773    164747
#HE627-pro-St5-F02-2_S14_R1_001.fastq    184437    169829
#HE627-pro-St5-F02-3_S15_R1_001.fastq    208053    191119
#HE627-pro-St6-F02-1_S16_R1_001.fastq    190061    176299
#HE627-pro-St6-F02-2_S17_R1_001.fastq    204907    191674
#HE627-pro-St6-F02-3_S18_R1_001.fastq    188951    175231
#HE627-pro-St7-F02-1_S19_R1_001.fastq    151897    138255
#HE627-pro-St7-F02-2_S20_R1_001.fastq    191043    177679
#HE627-pro-St7-F02-3_S21_R1_001.fastq    207484    192475
#HE627-pro-St8-F02-1_S22_R1_001.fastq    161446    148343
#HE627-pro-St8-F02-2_S23_R1_001.fastq    181153    167143
#HE627-pro-St8-F02-3_S24_R1_001.fastq    204499    190292




#DE-REPLICATE and keep going only with existing files

exists <- file.exists(fnFs.qualFiltTrim)
fnFs.deRep <- derepFastq(fnFs.qualFiltTrim[exists], verbose=TRUE)
fnRs.deRep <- derepFastq(fnRs.qualFiltTrim[exists], verbose=TRUE)
names(fnFs.deRep) <- sample.names[exists]
names(fnRs.deRep) <- sample.names[exists]


#LEARN ERRORS AND PLOT

errF <- learnErrors(fnFs.deRep, multithread=10,randomize=TRUE, nbases = 1e8)
#124879920 total bases in 520333 reads from 3 samples will be used for learning the error rates.
errR <- learnErrors(fnRs.deRep, multithread=10,randomize=TRUE, nbases = 1e8)
#107428000 total bases in 671425 reads from 4 samples will be used for learning the error rates.

errF
errR



#SAMPLE INFERENCE

dadaFs <- dada(fnFs.deRep, err=errF, multithread=TRUE)
#Sample 1 - 194261 reads in 16683 unique sequences.
#Sample 2 - 88293 reads in 8711 unique sequences.
#Sample 3 - 174464 reads in 14083 unique sequences.
#Sample 4 - 173895 reads in 14886 unique sequences.
#Sample 5 - 166464 reads in 13464 unique sequences.
#Sample 6 - 168844 reads in 14236 unique sequences.
#Sample 7 - 138823 reads in 13589 unique sequences.
#Sample 8 - 167684 reads in 14104 unique sequences.
#Sample 9 - 146694 reads in 13043 unique sequences.
#Sample 10 - 181558 reads in 14746 unique sequences.
#Sample 11 - 194329 reads in 15881 unique sequences.
#Sample 12 - 142436 reads in 13184 unique sequences.
#Sample 13 - 170142 reads in 14679 unique sequences.
#Sample 14 - 126839 reads in 11502 unique sequences.
#Sample 15 - 152840 reads in 13204 unique sequences.
#Sample 16 - 644103 reads in 31415 unique sequences.
#Sample 17 - 196206 reads in 14940 unique sequences.
#Sample 18 - 153985 reads in 11979 unique sequences.
#Sample 19 - 178829 reads in 15892 unique sequences.
#Sample 20 - 157450 reads in 13333 unique sequences.
#Sample 21 - 137577 reads in 11604 unique sequences.
#Sample 22 - 177908 reads in 14438 unique sequences.
#Sample 23 - 156986 reads in 12824 unique sequences.
#Sample 24 - 179954 reads in 14267 unique sequences.
#Sample 25 - 164747 reads in 14219 unique sequences.
#Sample 26 - 169829 reads in 14263 unique sequences.
#Sample 27 - 191119 reads in 14230 unique sequences.
#Sample 28 - 176299 reads in 14473 unique sequences.
#Sample 29 - 191674 reads in 15611 unique sequences.
#Sample 30 - 175231 reads in 14521 unique sequences.
#Sample 31 - 138255 reads in 11913 unique sequences.
#Sample 32 - 177679 reads in 13965 unique sequences.
#Sample 33 - 192475 reads in 15350 unique sequences.
#Sample 34 - 148343 reads in 12269 unique sequences.
#Sample 35 - 167143 reads in 13119 unique sequences.
#Sample 36 - 190292 reads in 14471 unique sequences.

dadaRs <- dada(fnRs.deRep, err=errR, multithread=TRUE)
#Sample 1 - 194261 reads in 10866 unique sequences.
#Sample 2 - 88293 reads in 7686 unique sequences.
#Sample 3 - 174464 reads in 10474 unique sequences.
#Sample 4 - 173895 reads in 10734 unique sequences.
#Sample 5 - 166464 reads in 9700 unique sequences.
#Sample 6 - 168844 reads in 9958 unique sequences.
#Sample 7 - 138823 reads in 8816 unique sequences.
#Sample 8 - 167684 reads in 11589 unique sequences.
#Sample 9 - 146694 reads in 9806 unique sequences.
#Sample 10 - 181558 reads in 10585 unique sequences.
#Sample 11 - 194329 reads in 12413 unique sequences.
#Sample 12 - 142436 reads in 10023 unique sequences.
#Sample 13 - 170142 reads in 12426 unique sequences.
#Sample 14 - 126839 reads in 8768 unique sequences.
#Sample 15 - 152840 reads in 9771 unique sequences.
#Sample 16 - 644103 reads in 23106 unique sequences.
#Sample 17 - 196206 reads in 10423 unique sequences.
#Sample 18 - 153985 reads in 8910 unique sequences.
#Sample 19 - 178829 reads in 11793 unique sequences.
#Sample 20 - 157450 reads in 10869 unique sequences.
#Sample 21 - 137577 reads in 9092 unique sequences.
#Sample 22 - 177908 reads in 11379 unique sequences.
#Sample 23 - 156986 reads in 9273 unique sequences.
#Sample 24 - 179954 reads in 10914 unique sequences.
#Sample 25 - 164747 reads in 11174 unique sequences.
#Sample 26 - 169829 reads in 10081 unique sequences.
#Sample 27 - 191119 reads in 10654 unique sequences.
#Sample 28 - 176299 reads in 11367 unique sequences.
#Sample 29 - 191674 reads in 11916 unique sequences.
#Sample 30 - 175231 reads in 11246 unique sequences.
#Sample 31 - 138255 reads in 9244 unique sequences.
#Sample 32 - 177679 reads in 10322 unique sequences.
#Sample 33 - 192475 reads in 11919 unique sequences.
#Sample 34 - 148343 reads in 9623 unique sequences.
#Sample 35 - 167143 reads in 9617 unique sequences.
#Sample 36 - 190292 reads in 11331 unique sequences.

#inspect obtained ASVs
dadaFs[[1]]  #593 sequence variants were inferred from 16683 input unique sequences.
dadaRs[[1]]  #516 sequence variants were inferred from 10866 input unique sequences.

#MERGE PAIRED ENDS
mergers <- mergePairs(dadaFs, fnFs.deRep, dadaRs, fnRs.deRep, minOverlap=20,verbose=TRUE) ## min overlap maybe wrong? 
#Inspect the merger data.frame from the first sample
head(mergers[[1]])



#CONSTRUCT SEQUENCE TABLE
seqtab <- makeSequenceTable(mergers)
dim(seqtab)  # 36 2896

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


#REMOVE CHIMERAS AND EXPORT ASV TABLE
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 546 bimeras out of 3130 input sequences.
dim(seqtab.nochim)
#36 2501
write.csv(t(seqtab.nochim), "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/NextSeq_seqtab_all_16S.csv", quote=FALSE ) # CHANGE ME to output directory.
saveRDS(filterOut, file = "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/NextSeq_seqtab_nochim_16S.RDA") # CHANGE ME to output directory.

#obtain the frequencies of chimeras in the dataset
sum(seqtab.nochim)/sum(seqtab)
# 0.995767



#TRACK READS
getN <- function(x) sum(getUniques(x))
#bind columns with same length
track1 <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
rownames(filterOut) <- sample.names
#merge with columns with different length and add NA when there is no value because there were no output sequences from filtering
track_reads <-merge (filterOut, track1, by = 0, all = TRUE)
colnames(track_reads) <- c("sample", "input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
#change NAs by zeros and add rownames
track_reads[is.na(track_reads)] <- 0
rownames(track_reads) <- track_reads$sample
track_reads <- track_reads[,-1]

write.table(track_reads, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/NextSeq_16S_tracked_reads_Philina.txt") # CHANGE ME to output directory.



#ASSIGN TAXONOMY
# for 16S use Silva !!
db_dir <- file.path("/isibhv/projects/p_bioinf2/dbs/dada2_dbs")
taxDB <- "silva_nr99_v138.1_wSpecies_train_set_Pros.fa"
taxLvs <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

taxa <- assignTaxonomy(seqtab, file.path(db_dir,taxDB), taxLevels = taxLvs ,multithread=TRUE, minBoot=50)

write.csv(taxa, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/NextSeq_16S_taxa.csv")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
bind <- cbind(t(seqtab.nochim), taxa)

write.csv(bind, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S//NextSeq_16S_assigned_all_Philina.csv", quote=FALSE ) # CHANGE ME to output directory.


##
### bzw. über zweites Script, um nicht nochmal alles laufen zu lassen: 

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

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
bind <- cbind(t(seqtab.nochim), taxa)  # funktioniert nicht 

write.csv(bind, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/NextSeq_16S_assigned_all_Philina.csv", quote=FALSE ) # CHANGE ME to output directory.





