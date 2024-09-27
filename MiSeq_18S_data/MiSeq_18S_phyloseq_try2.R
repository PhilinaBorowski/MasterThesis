## Phyloseq MiSeq 18S - new 

library(dplyr)
library(ggplot2)
library(phyloseq)
library(tibble)
library(microbiomeMarker)


## set wd
setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/Phyloseq_analysis/")

## read in data
## but reduced data 
# so here instead og 96 samples: only 92 samples 

tax <- read.csv("MiSeq_metabarcoding_wo_metazoa_92.csv")

str(tax)

# abundance as num
# rest as factor? 

