### MiSeq Plots ###


library(dada2)
library(ShortRead)
library(Biostrings)
library(stringr)
library(R.utils)
library(plyr)


### set wd ###
setwd("C:/Users/phili/OneDrive/Desktop/Master_Thesis")


raw_dir <- setwd("C:/Users/phili/OneDrive/Desktop/Master_Thesis/MiSeq results/FastQ files_unzip") # run 2x
list.files(raw_dir)

preFilt_dir <- file.path(wd,"preFilt_18S") # vorgefilterte 
primerCut5_dir <- file.path(wd,"primerCut5_18S") # primer 5' cuts
primerCut3_dir <- file.path(wd,"primerCut3_18S") # primer 3' cuts
qualFiltTrim_dir <- file.path(wd,"qualFiltTrim_18S") # quality filter trim 


plotQualityProfile(fnRs.raw[1:2]) #600:450
plotQualityProfile(fnRs.raw) #2000:3000

plotQualityProfile(fnFs.preFilt[1:2])
plotQualityProfile(fnFs.preFilt)

plotQualityProfile(fnRs.preFilt[1:2])
plotQualityProfile(fnRs.preFilt)

#plotQualityProfile(fnFs.primerCut3[1:2])
#plotQualityProfile(fnFs.primerCut3)

#plotQualityProfile(fnRs.primerCut3[1:2])
#plotQualityProfile(fnRs.primerCut3)

#plotQualityProfile(fnFs.primerCut5[1:2])
#plotQualityProfile(fnRF.primerCut5)

#plotQualityProfile(fnRs.primerCut5[1:2])
#plotQualityProfile(fnRs.primerCut5)

#plotQualityProfile(fnFs.qualFiltTrim[1:2])
#plotQualityProfile(fnFs.qualFiltTrim)

#plotQualityProfile(fnRs.qualFiltTrim[1:2])
#plotQualityProfile(fnRs.qualFiltTrim)

