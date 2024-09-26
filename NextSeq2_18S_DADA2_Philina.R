### DADA2 NextSeq 2 - Eukaryotes

## sample St7_F3_3 is empty/zero

library(dada2)
library(ShortRead)
library(Biostrings)
library(stringr)
library(R.utils)
library(plyr)


## set wd
setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2")
wd <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2")

raw_dir <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/fastq_files_unzip")
raw_dir <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/fastq_files_unzip")# run 2x
#View(raw_dir)
list.files(raw_dir)


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
filterAndTrim(fnFs.raw, fnFs.preFilt, fnRs.raw, fnRs.preFilt, truncQ = 2, minQ = 2,
              minLen = 50, maxN = 0, multithread = FALSE, compress = FALSE)
# filterAndTrim(input, output)
# truncQ = Truncate reads at the first instance of a quality score less than or equal to truncQ
# minQ = After truncation, reads contain a quality score less than minQ will be discarded.
# minLen = Remove reads with length less than minLen. minLen is enforced after trimming and truncation.
# maxN = After truncation, sequences with more than maxN Ns will be discarded; DADA2 requires N to be 0!
# compress = T for fastq.gz files, compress = F for fastq files 



### 
#IDENTIFY PRIMERs

FWD_PRIMER <- "CCAGCASCYGCGGTAATTCC"                          # primer Sequence (5’ - 3')
REV_PRIMER <- "ACTTTCGTTCTTGAT" 

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
#FWD_PRIMER.ForwardReads  241018          0       0       0
#FWD_PRIMER.ReverseReads       4          0       0      51
#REV_PRIMER.ForwardReads       4          0       0      42
#REV_PRIMER.ReverseReads  246846          0       0       0


#cutadapt <- "/isibhv/projects/AG_John/Expeditions/HE627_Data/silva+pr2+cutadapt/cutadapt.exe" # CHANGE ME to the cutadapt path on your machine
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
# reverse: 15 bases. so 14
# forward: 20, so 19

for(i in seq_along(fnFs.preFilt)) {
  system2("cutadapt", args = c("-g", paste("\"",FWD_PRIMER,";min_overlap=19;max_error_rate=0.15","\"",sep=""),
                             "-G", paste("\"",REV_PRIMER,";min_overlap=14;max_error_rate=0.15","\"",sep=""),
                             "--discard-untrimmed", "-o", fnFs.primerCut5[i], "-p", fnRs.primerCut5[i],
                             fnFs.preFilt[i], fnRs.preFilt[i]))
  
}

###
### check if it worked; should be 0
#check one sample - 1st sample - St1-F02-1
rbind(FWD_PRIMER.ForwardReads = sapply(FWD_PRIMER.orients, primerHits, fn = fnFs.primerCut5[[1]]),
      FWD_PRIMER.ReverseReads = sapply(FWD_PRIMER.orients, primerHits, fn = fnRs.primerCut5[[1]]),
      REV_PRIMER.ForwardReads = sapply(REV_PRIMER.orients, primerHits, fn = fnFs.primerCut5[[1]]),
      REV_PRIMER.ReverseReads = sapply(REV_PRIMER.orients, primerHits, fn = fnRs.primerCut5[[1]]))

#                         Forward Complement Reverse RevComp
#FWD_PRIMER.ForwardReads       0          0       0       0
#FWD_PRIMER.ReverseReads       0          0       0      29
#REV_PRIMER.ForwardReads       0          0       0      26
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

## results
#                         Forward Complement Reverse RevComp
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
                           truncLen=c(270,220),         # 18 S !!!
                           verbose = TRUE, rm.phix = TRUE,
                           compress = TRUE, multithread = FALSE)

print(filterOut)    # output ~ 76 - 80 % 

#                                    reads.in     reads.out
#HE627-euk-St1-dDNA_S111_R1_001.fastq    248970    198914
#HE627-euk-St1-F02-1_S1_R1_001.fastq     306846    242328
#HE627-euk-St1-F02-2_S2_R1_001.fastq     316426    249555
#HE627-euk-St1-F02-3_S3_R1_001.fastq     318080    254029
#HE627-euk-St1-F3-1_S4_R1_001.fastq      274060    221416
#HE627-euk-St1-F3-2_S5_R1_001.fastq      331018    266165
#HE627-euk-St1-F3-3_S6_R1_001.fastq      294383    237858
#HE627-euk-St1-N_S8_R1_001.fastq         207599    167194
#HE627-euk-St1-P1_S7_R1_001.fastq        313508    252136
#HE627-euk-St10-dDNA_S120_R1_001.fastq   242752    195772
#HE627-euk-St10-F02-1_S29_R1_001.fastq   347996    267561
#HE627-euk-St10-F02-2_S30_R1_001.fastq   292847    231601
#HE627-euk-St10-F02-3_S31_R1_001.fastq   251235    192349
#HE627-euk-St10-F3-1_S56_R1_001.fastq    271350    216445
#HE627-euk-St10-F3-2_S57_R1_001.fastq    314060    249992
#HE627-euk-St10-F3-3_S58_R1_001.fastq    263166    209720
#HE627-euk-St10-N_S106_R1_001.fastq      286132    229852
#HE627-euk-St10-P-1_S87_R1_001.fastq     318099    254867
#HE627-euk-St10-P-2_S88_R1_001.fastq     289194    232692
#HE627-euk-St10-P-3_S89_R1_001.fastq     262080    209561
#HE627-euk-St11-dDNA_S121_R1_001.fastq   269088    215021
#HE627-euk-St11-F02-1_S32_R1_001.fastq   245750    190264
#HE627-euk-St11-F02-2_S33_R1_001.fastq   258951    201649
#HE627-euk-St11-F02-3_S34_R1_001.fastq   328325    260192
#HE627-euk-St11-F3-1_S59_R1_001.fastq    305669    245803
#HE627-euk-St11-F3-2_S60_R1_001.fastq    266336    213622
#HE627-euk-St11-N_S107_R1_001.fastq      356973    285111
#HE627-euk-St11-P-1_S90_R1_001.fastq     311391    251546
#HE627-euk-St11-P-2_S91_R1_001.fastq     257079    201733
#HE627-euk-St11-P-3_S92_R1_001.fastq     276445    226038
#HE627-euk-St12-dDNA_S122_R1_001.fastq   336220    269955
#HE627-euk-St12-F02-1_S35_R1_001.fastq   277591    218514
#HE627-euk-St12-F02-2_S36_R1_001.fastq   397760    315563
#HE627-euk-St12-F02-3_S37_R1_001.fastq   260635    208265
#HE627-euk-St12-F3-1_S61_R1_001.fastq    274606    221049
#HE627-euk-St12-F3-2_S62_R1_001.fastq    260913    209414
#HE627-euk-St12-F3-3_S63_R1_001.fastq    191532    154007
#HE627-euk-St12-N_S108_R1_001.fastq      284522    228568
#HE627-euk-St12-P-1_S93_R1_001.fastq     243022    195439
#HE627-euk-St12-P-2_S94_R1_001.fastq     300265    239659
#HE627-euk-St12-P-3_S95_R1_001.fastq     288531    228708
#HE627-euk-St13-dDNA_S123_R1_001.fastq   284942    228169
#HE627-euk-St13-F02-1_S38_R1_001.fastq   348433    275286
#HE627-euk-St13-F02-3_S39_R1_001.fastq   227198    182925
#HE627-euk-St13-F3-1_S64_R1_001.fastq    227341    181727
#HE627-euk-St13-F3-2_S65_R1_001.fastq    245003    198814
#HE627-euk-St13-F3-3_S66_R1_001.fastq     71876     57642
#HE627-euk-St13-N_S109_R1_001.fastq      282995    227661
#HE627-euk-St13-P-1_S96_R1_001.fastq     312611    252159
#HE627-euk-St13-P-2_S97_R1_001.fastq     274745    218448
#HE627-euk-St13-P-3_S98_R1_001.fastq     266258    213581
#HE627-euk-St2-dDNA_S112_R1_001.fastq    296743    237685
#HE627-euk-St2-F02-1_S9_R1_001.fastq     215456    171933
#HE627-euk-St2-F02-2_S10_R1_001.fastq    318945    252865
#HE627-euk-St2-F02-3_S11_R1_001.fastq    238149    191447
#HE627-euk-St2-F3-2_S40_R1_001.fastq     257780    208082
#HE627-euk-St2-F3-3_S41_R1_001.fastq     314581    252486
#HE627-euk-St2-N_S99_R1_001.fastq        301019    242403
#HE627-euk-St2-P-1_S67_R1_001.fastq      291094    236034
#HE627-euk-St2-P-2_S68_R1_001.fastq      273809    219922
#HE627-euk-St2-P-3_S69_R1_001.fastq      125734    101131
#HE627-euk-St3-dDNA_S113_R1_001.fastq    315133    253208
#HE627-euk-St3-F02-2_S12_R1_001.fastq    377188    301577
#HE627-euk-St3-F02-3_S13_R1_001.fastq    219407    174756
#HE627-euk-St3-F3-1_S42_R1_001.fastq     276371    221066
#HE627-euk-St3-F3-2_S43_R1_001.fastq     218933    174964
#HE627-euk-St3-F3-3_S44_R1_001.fastq     280196    224395
#HE627-euk-St3-N_S100_R1_001.fastq       270942    216943
#HE627-euk-St3-P-1_S70_R1_001.fastq      279207    224891
#HE627-euk-St3-P-2_S71_R1_001.fastq      282246    225332
#HE627-euk-St3-P-3_S72_R1_001.fastq      268155    214741
#HE627-euk-St4-dDNA_S114_R1_001.fastq    282944    227371
#HE627-euk-St4-F02-1_S14_R1_001.fastq    396115    311991
#HE627-euk-St4-F02-2_S15_R1_001.fastq    271652    215716
#HE627-euk-St4-F02-3_S16_R1_001.fastq    294767    231812
#HE627-euk-St4-F3-1_S45_R1_001.fastq     234723    188622
#HE627-euk-St4-F3-2_S46_R1_001.fastq     226965    184025
#HE627-euk-St4-F3-3_S47_R1_001.fastq     344262    274675
#HE627-euk-St4-N_S101_R1_001.fastq       278139    222411
#HE627-euk-St4-P-1_S73_R1_001.fastq      286235    231360
#HE627-euk-St4-P-2_S74_R1_001.fastq      281412    226664
#HE627-euk-St4-P-3_S75_R1_001.fastq      258320    207050
#HE627-euk-St5-dDNA_S115_R1_001.fastq    280804    223378
#HE627-euk-St5-F02-1_S17_R1_001.fastq    303934    242762
#HE627-euk-St5-F02-2_S18_R1_001.fastq    276071    215562
#HE627-euk-St5-F02-3_S19_R1_001.fastq    232192    180837
#HE627-euk-St5-F3-3_S48_R1_001.fastq     248905    199712
#HE627-euk-St5-N_S102_R1_001.fastq       307116    245214
#HE627-euk-St5-N-150_S110_R1_001.fastq   278967    222375
#HE627-euk-St5-P-1_S76_R1_001.fastq      305367    244927
#HE627-euk-St5-P-2_S77_R1_001.fastq      270039    214780
#HE627-euk-St5-P-3_S78_R1_001.fastq      322285    258093
#HE627-euk-St6-dDNA_S116_R1_001.fastq    283717    226074
#HE627-euk-St6-F02-1_S20_R1_001.fastq    265402    202909
#HE627-euk-St6-F02-2_S21_R1_001.fastq    264687    204380
#HE627-euk-St6-F02-3_S22_R1_001.fastq    275679    217069
#HE627-euk-St6-F3-1_S49_R1_001.fastq     248137    193961
#HE627-euk-St6-F3-2_S50_R1_001.fastq     224228    177891
#HE627-euk-St6-N_S103_R1_001.fastq       273963    214184
#HE627-euk-St6-P-1_S79_R1_001.fastq      289620    228537
#HE627-euk-St6-P-3_S80_R1_001.fastq      296779    238355
#HE627-euk-St7-dDNA_S117_R1_001.fastq    294315    233504
#HE627-euk-St7-F02-1_S23_R1_001.fastq    275812    209645
#HE627-euk-St7-F02-2_S24_R1_001.fastq    358614    277952
#HE627-euk-St7-F02-3_S25_R1_001.fastq    292455    226002
#HE627-euk-St7-F3-2_S51_R1_001.fastq   21021088  16681861
#HE627-euk-St7-N_S104_R1_001.fastq       270179    215908
#HE627-euk-St7-P-1_S81_R1_001.fastq      255218    205865
#HE627-euk-St7-P-2_S82_R1_001.fastq      223212    179050
#HE627-euk-St7-P-3_S83_R1_001.fastq      210699    165809
#HE627-euk-St8-dDNA_S118_R1_001.fastq    248947    199109
#HE627-euk-St8-F02-1_S26_R1_001.fastq    294440    226829
#HE627-euk-St8-F02-2_S27_R1_001.fastq    306937    244269
#HE627-euk-St8-F02-3_S28_R1_001.fastq    220165    173457
#HE627-euk-St8-F3-1_S53_R1_001.fastq     263086    211298
#HE627-euk-St8-F3-2_S54_R1_001.fastq     269260    211180
#HE627-euk-St8-F3-3_S55_R1_001.fastq   22867284  17871017
#HE627-euk-St8-N_S105_R1_001.fastq       288441    229397
#HE627-euk-St8-P-1_S84_R1_001.fastq      124551    101208
#HE627-euk-St8-P-2_S85_R1_001.fastq      189956    151105
#HE627-euk-St8-P-3_S86_R1_001.fastq      281619    223649
#HE627-euk-St9-dDNA_S119_R1_001.fastq    220578    174326



#DE-REPLICATE and keep going only with existing files

exists <- file.exists(fnFs.qualFiltTrim)
fnFs.deRep <- derepFastq(fnFs.qualFiltTrim[exists], verbose=TRUE)
fnRs.deRep <- derepFastq(fnRs.qualFiltTrim[exists], verbose=TRUE)
names(fnFs.deRep) <- sample.names[exists]
names(fnRs.deRep) <- sample.names[exists]


#LEARN ERRORS AND PLOT

errF <- learnErrors(fnFs.deRep, multithread=10,randomize=TRUE, nbases = 1e8)
#146769840 total bases in 543592 reads from 2 samples will be used for learning the error rates.

errR <- learnErrors(fnRs.deRep, multithread=10,randomize=TRUE, nbases = 1e8)
#131803760 total bases in 599108 reads from 3 samples will be used for learning the error rates.


errF
errR



#SAMPLE INFERENCE

dadaFs <- dada(fnFs.deRep, err=errF, multithread=TRUE)
#Warning messages:
#1: In rval[, 1:ncol(tt)] + tt : NAs produced by integer overflow
#2: In rval[, 1:ncol(tt)] + tt : NAs produced by integer overflow

## but ran all 122 samples

dadaRs <- dada(fnRs.deRep, err=errR, multithread=TRUE)
#Warning messages:
#1: In rval[, 1:ncol(tt)] + tt : NAs produced by integer overflow
#2: In rval[, 1:ncol(tt)] + tt : NAs produced by integer overflow

## but ran all 122 samples

#inspect obtained ASVs
dadaFs[[1]]  #585 sequence variants were inferred from 9033 input unique sequences. OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16
dadaRs[[1]]  #437 sequence variants were inferred from 10364 input unique sequences. OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

#MERGE PAIRED ENDS
mergers <- mergePairs(dadaFs, fnFs.deRep, dadaRs, fnRs.deRep, minOverlap=20,verbose=TRUE)
#Inspect the merger data.frame from the first sample
head(mergers[[1]])



#CONSTRUCT SEQUENCE TABLE
seqtab <- makeSequenceTable(mergers)
dim(seqtab)  # 122 36001

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


#REMOVE CHIMERAS AND EXPORT ASV TABLE
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 21930 bimeras out of 36001 input sequences.
dim(seqtab.nochim)
# 122 14071
write.csv(t(seqtab.nochim), "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/NextSeq2_18S_seqtab_all.csv", quote=FALSE ) # CHANGE ME to output directory.
saveRDS(filterOut, file = "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/NextSeq2_18S_seqtab_nochim.RDA") # CHANGE ME to output directory.
#obtain the frequencies of chimeras in the dataset
sum(seqtab.nochim)/sum(seqtab)
# 0.9757045




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

write.table(track_reads, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/NextSeq2_18S_tracked_reads.txt") # CHANGE ME to output directory.




## ASSIGN TAXONOMY
# pr2 for euks
seqtab.nochim_t <- t(seqtab.nochim)
db_dir <- file.path("/isibhv/projects/p_bioinf2/dbs/dada2_dbs")
taxDB <- "pr2_version_5.0.0_SSU_dada2.fasta"
taxLvs <- c("Kingdom", "Supergroup","Phylum","Division", "Class", "Order", "Family", "Genus", "Species")

taxa <- assignTaxonomy(seqtab.nochim_t, file.path(db_dir,taxDB), taxLevels = taxLvs ,multithread=TRUE, minBoot=70)

write.csv(taxa, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/NextSeq2_18S_taxa.csv")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
bind <- cbind(t(seqtab.nochim), taxa)

write.csv(bind, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/NextSeq2_18S_assigned_all_Philina.csv", quote=FALSE ) # CHANGE ME to output directory.





