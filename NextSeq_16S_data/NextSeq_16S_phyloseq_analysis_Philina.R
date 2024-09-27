## NextSeq - 16S - phyloseq analysis
#

## read in packages

library(dplyr)
library(ggplot2)
library(phyloseq)
library(tibble)
library(microbiomeMarker)


## set wd
setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/phyloseq_analysis/")

## read in files

asvs <- read.csv("ASV_abundances_18S_miseq.csv")
taxa <- read.csv("ASV_table_18S_miseq.csv")                 ###ASV0001 - ASV4143
meta_samples <- read.csv("meta_phyloseq_18S_miseq.csv")  

## add row names

asvs <- asvs %>% tibble::column_to_rownames("otu")
taxa <- taxa %>% tibble::column_to_rownames("otu")
meta_samples <- meta_samples %>% tibble::column_to_rownames("sample")


asvs <- as.matrix(asvs)
taxa <- as.matrix(taxa)

str(asvs)

OTU = otu_table(asvs, taxa_are_rows = T)
TAX = tax_table(taxa)
samples = sample_data(meta_samples)

## got a mismatch in names
sample_names(OTU)            # all "-" got converted into "."
sample_names(samples)        # converted all "-" into "." manually

MiSeq_18S <- phyloseq(OTU, TAX, samples)
MiSeq_18S

sample_names(MiSeq_18S)
rank_names(MiSeq_18S)
sample_variables(MiSeq_18S)


## didnt do that
## should I?
#carbom <- subset_taxa(carbom, Division %in% c("Chlorophyta", "Dinophyta", "Cryptophyta", 
#                                              "Haptophyta", "Ochrophyta", "Cercozoa"))
#carbom <- subset_taxa(carbom, !(Class %in% c("Syndiniales", "Sarcomonadea")))
#carbom


plot_bar(MiSeq_18S, fill = "Family")   # takes long, cant see shit
plot_bar(MiSeq_18S, fill = "Class")    # takes long, cant see shit 
plot_bar(MiSeq_18S, fill = "Phylum")

## remove sample St7.F3.3, St13.F02.1, St2.F3.1
MiSeq_18S_cleaned <- subset_samples(MiSeq_18S, sample_names(MiSeq_18S) !="St7.F3.3" & 
                                      sample_names(MiSeq_18S) != "St13.F02.1" & 
                                      sample_names(MiSeq_18S) != "St2.F3.1")



#### OR 
## do the thing Uwe suggested
## but
## idk what




# check bar plot again
plot_bar(MiSeq_18S_cleaned, fill = "Class") 



# normalise it? 
MiSeq_18S_norm <- normalize(MiSeq_18S_cleaned, method = "CLR")  # package microbiomeMarker



plot_bar(MiSeq_18S_norm, fill = "Class") + geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")
plot_bar(MiSeq_18S_norm, fill = "Phylum")

MiSeq_18S_norm_fraction <- merge_samples(MiSeq_18S_norm, "fraction")
## doesnt work; but maybe can add it manually? what is pico and what is nano? 



plot_heatmap(MiSeq_18S_norm, method = "NMDS", distance = "bray")


# plot a diversity
plot_richness(MiSeq_18S_norm, measures = "Shannon")



## ordination

MiSeq_18S_ord <- ordinate(MiSeq_18S, "NMDS", "bray")
plot_ordination(MiSeq_18S, MiSeq_18S_ord, type = "taxa", color = "Order", shape = "Class")
plot_ordination(MiSeq_18S, MiSeq_18S_ord, type = "taxa", color = "Order", shape = "Phylum")
plot_ordination(MiSeq_18S, MiSeq_18S_ord, type = "taxa", color = "Order", shape = "Kingdom")


MiSeq_18S_abund <- filter_taxa(MiSeq_18S, function(x) sum(x > total*0.2) > 0, TRUE)

plot_net(MiSeq_18S_abund, distance = "(A+B-2*J)/(A+B)", type = "taxa", maxdist = 0.6, color = "Order")


