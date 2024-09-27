### MiSew Analysis - FastQ Files ###

### Install packages ###

#install.packages("dada2")
#install.packages("ShortRead")
#install.packages("Biostrings")
#install.packages("stringr")
#install.packages("R.utils")
#install.packages("Rcpp")

library(dada2)
library(ShortRead)
library(Biostrings)
library(stringr)
library(R.utils)
library(plyr)


### set wd ###
setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/MiSeq_results/")
wd <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/MiSeq_results/")

raw_dir <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/MiSeq_results/FastQ_files_unzip/")
raw_dir <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/MiSeq_results/FastQ_files_unzip/")# run 2x
#View(raw_dir)
list.files(raw_dir)
#?list.files()

###
### create folders for later analysis
###

preFilt_dir <- file.path(wd,"preFilt_18S") # vorgefilterte 
primerCut5_dir <- file.path(wd,"primerCut5_18S") # primer 5' cuts
primerCut3_dir <- file.path(wd,"primerCut3_18S") # primer 3' cuts
qualFiltTrim_dir <- file.path(wd,"qualFiltTrim_18S") # quality filter trim 


###
### "CONSTRUCT NEEDED FILE LISTS"
### fn = full names

fnFs.raw <- sort(list.files(raw_dir, pattern = "_L001_R1_001.fastq", full.names = TRUE))
fnRs.raw <- sort(list.files(raw_dir, pattern = "_L001_R2_001.fastq", full.names = TRUE))
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
sample.names <- str_remove(basename(fnFs.raw),"_L001_R1_001.fastq")
head(sample.names)

basename(fnRs.raw)
basename(sample.names)

###
#PREFILTERING
### doesnt work on windows haha -> have to put "multithread = FALSE" instead of "TRUE"
filterAndTrim(fnFs.raw, fnFs.preFilt, fnRs.raw, fnRs.preFilt, truncQ = 2, minQ = 2,
              minLen = 50, maxN = 0, multithread = FALSE, compress = FALSE)
# filterAndTrim(input, output)
# truncQ = Truncate reads at the first instance of a quality score less than or equal to truncQ
# minQ = After truncation, reads contain a quality score less than minQ will be discarded.
# minLen = Remove reads with length less than minLen. minLen is enforced after trimming and truncation.
# maxN = After truncation, sequences with more than maxN Ns will be discarded; DADA2 requires N to be 0!

## took about 50-55 min

### 
#IDENTIFY PRIMER
#FWD_PRIMER="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCAGCASCYGCGGTAATTCC"  # 5' - 3' + illumina tail
#REV_PRIMER="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGACTTTCGTTCTTGAT"      # 5' - 3' + illumina tail
#dont use these cause they are already gone after demultiplexing


FWD_PRIMER <- "CCAGCASCYGCGGTAATTCC"                          # primer Sequence (5’ - 3')
REV_PRIMER <- "ACTTTCGTTCTTGAT"                               # primer Sequence (5’ - 3')

### 
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


#                         Forward  Complement  Reverse  RevComp
#FWD_PRIMER.ForwardReads   20876          0       0       0
#FWD_PRIMER.ReverseReads       0          0       0      33
#REV_PRIMER.ForwardReads       0          0       0       3
#REV_PRIMER.ReverseReads   20904          0       0       0


###
### REMOVE PRIMERS
###

#cutadapt <- "C:/Users/phili/Downloads/cutadapt.exe" # CHANGE ME to the cutadapt path on your machine
system2("cutadapt", args = "--version") # Run shell commands from R
#2.8


### 
###
#create output dirs
if(!dir.exists(primerCut5_dir)) dir.create(primerCut5_dir)
if(!dir.exists(primerCut3_dir)) dir.create(primerCut3_dir)


FWD_PRIMER.RC <- dada2:::rc(FWD_PRIMER)
REV_PRIMER.RC <- dada2:::rc(REV_PRIMER)

# Run Cutadapt both primer at 5' end
# min_overlap primer bases -1
# reverse: 15 bases. so 14
# forward: 20, so 19

for(i in seq_along(fnFs.preFilt)) {
  system2("cutadapt", args = c("-g", paste("\"",FWD_PRIMER,";min_overlap=19;max_error_rate=0.15","\"",sep=""),
                             "-G", paste("\"",REV_PRIMER,";min_overlap=14;max_error_rate=0.15","\"",sep=""),
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
#FWD_PRIMER.ForwardReads       0          0       0       0
#FWD_PRIMER.ReverseReads       0          0       0       6
#REV_PRIMER.ForwardReads       0          0       0       1
#REV_PRIMER.ReverseReads       0          0       0       0



###
# Run Cutadapt for primer at 3' end
# min_overlap primer bases -1
# reverse: 15 bases. so 14
# forward: 20, so 19
for(i in seq_along(fnFs.primerCut5)) {
  system2("cutadapt", args = c("-a", paste("\"",REV_PRIMER.RC,";min_overlap=14;max_error_rate=0.15","\"",sep=""),
                             "-A", paste("\"",FWD_PRIMER.RC,";min_overlap=19;max_error_rate=0.15","\"",sep=""),
                             "-o", fnFs.primerCut3[i], "-p", fnRs.primerCut3[i],
                             fnFs.primerCut5[i], fnRs.primerCut5[i]))
  
}

###
### check if it worked with 1st sample (should be 0)
rbind(FWD_PRIMER.ForwardReads = sapply(FWD_PRIMER.orients, primerHits, fn = fnFs.primerCut3[[1]]),
      FWD_PRIMER.ReverseReads = sapply(FWD_PRIMER.orients, primerHits, fn = fnRs.primerCut3[[1]]),
      REV_PRIMER.ForwardReads = sapply(REV_PRIMER.orients, primerHits, fn = fnFs.primerCut3[[1]]),
      REV_PRIMER.ReverseReads = sapply(REV_PRIMER.orients, primerHits, fn = fnRs.primerCut3[[1]]))

#                           Forward Complement Reverse RevComp
#FWD_PRIMER.ForwardReads       0          0       0       0
#FWD_PRIMER.ReverseReads       0          0       0       0
#REV_PRIMER.ForwardReads       0          0       0       0
#REV_PRIMER.ReverseReads       0          0       0       0



### QUality Control
###
### Quality filtering

filterOut <- filterAndTrim(fnFs.primerCut3,fnFs.qualFiltTrim,
                           fnRs.primerCut3,fnRs.qualFiltTrim,
                           maxN=0,maxEE=c(2.7,2.2),
                           truncLen=c(270,220),   # 18 S !!
                           verbose = TRUE, rm.phix = TRUE,
                           compress = TRUE, multithread = FALSE)

print(filterOut)

###
#DE-REPLICATE and keep going only with existing files

exists <- file.exists(fnFs.qualFiltTrim)
fnFs.deRep <- derepFastq(fnFs.qualFiltTrim[exists], verbose=TRUE)
fnRs.deRep <- derepFastq(fnRs.qualFiltTrim[exists], verbose=TRUE)
names(fnFs.deRep) <- sample.names[exists]
names(fnRs.deRep) <- sample.names[exists]


#LEARN ERRORS AND PLOT

errF <- learnErrors(fnFs.deRep, multithread=10,randomize=TRUE, nbases = 1e8)
errR <- learnErrors(fnRs.deRep, multithread=10,randomize=TRUE, nbases = 1e8)

errF
errR



#SAMPLE INFERENCE

dadaFs <- dada(fnFs.deRep, err=errF, multithread=TRUE)
dadaRs <- dada(fnRs.deRep, err=errR, multithread=TRUE)
#inspect obtained ASVs
dadaFs[[1]]
dadaRs[[1]]

#MERGE PAIRED ENDS
mergers <- mergePairs(dadaFs, fnFs.deRep, dadaRs, fnRs.deRep, minOverlap=20,verbose=TRUE)
#Inspect the merger data.frame from the first sample
head(mergers[[1]])



#CONSTRUCT SEQUENCE TABLE
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


#REMOVE CHIMERAS AND EXPORT ASV TABLE
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#96 4602
write.csv(t(seqtab.nochim), "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/MiSeq_results/MiSeq_18S_seqtab_all.csv", quote=FALSE ) # CHANGE ME to output directory.
saveRDS(filterOut, file = "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/MiSeq_results/MiSeq_18S_seqtab_nochim.RDA") # CHANGE ME to output directory.
#obtain the frequencies of chimeras in the dataset
sum(seqtab.nochim)/sum(seqtab)
# 0.9564204




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

write.table(track_reads, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/MiSeq_results/MiSeq_18S_tracked_reads.txt") # CHANGE ME to output directory.

#ASSIGN TAXONOMY
seqtab.nochim_t <- t(seqtab.nochim)
db_dir <- file.path("/isibhv/projects/p_bioinf2/dbs/dada2_dbs")
taxDB <- "pr2_version_5.0.0_SSU_dada2.fasta"
taxLvs <- c("Kingdom", "Supergroup","Phylum","Division", "Class", "Order", "Family", "Genus", "Species")

taxa <- assignTaxonomy(seqtab.nochim_t, file.path(db_dir,taxDB), taxLevels = taxLvs ,multithread=TRUE, minBoot=70)

write.csv(taxa, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/MiSeq_results/MiSeq_18S_taxa.csv")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
bind <- cbind(seqtab.nochim, taxa)

write.csv(bind, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/MiSeq_results/MiSeq_18S_assigned_all_18S.csv", quote=FALSE ) # CHANGE ME to output directory.





