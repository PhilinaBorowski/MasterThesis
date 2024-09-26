## NextSeq 18S - look at filters and net/pump seperately

library(viridis)         
library(dplyr)
library(ggplot2)
library(phyloseq)
library(tibble)
library(ape)
#library(vegan)
library(microViz)
library(tidyr)
library(zCompositions)
library(compositions)
library(patchwork)
library(fantaxtic)
library(corrplot)


## set wd
setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/Phyloseq_analysis/")
wd <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/Phyloseq_analysis/")

# remove all of the NAs in annotation
# like NA in Kingdom

#data <- read.csv("NextSeq_18S_reduced_no_1rowsum_133.csv")  # 10425 ASVs

#new_df <- subset(data, Kingdom != "NA")                     # 10197 ASVs

#write.csv(new_df, "NextSeq_18S_reduced_no_1rowsum_KingdomNAs_133.csv")


# prepare ASV + tax tables 
# separate filter + pump/net/dDNA samples? 


## read in data

asvsf_n <- read.csv("NextSeq_18S_abundances_filters.csv")           ## need the ones wo metazoa + singletons !!!
asvsnp_n <- read.csv("NextSeq_18S_abundances_pump_net_deepDNA.csv")
taxa_n <- read.csv("NextSeq_18S_taxa.csv")                 ###ASV0001 - ASV2502
meta_samples_nf <- read.csv("metadata_phyloseq_NextSeq_18S_filters.csv")  
meta_samples_npd <- read.csv("metadata_phyloseq_NextSeq_18S_NPdeepDNA.csv")

asvsf_n <- asvsf_n %>% tibble::column_to_rownames("otu")
asvsnp_n <- asvsnp_n %>% tibble::column_to_rownames("otu")
taxa_n <- taxa_n %>% tibble::column_to_rownames("otu")
meta_samples_nf <- meta_samples_nf %>% tibble::column_to_rownames("sample")
meta_samples_npd <- meta_samples_npd %>% tibble::column_to_rownames("sample")


## prep for CLR + CLR 
#d.czm <- cmultRepl(t(d.1),  label=0, method="CZM", frac= 0.99, adjust=FALSE, z.warning = 0.99)
#d.clr <- t(apply(d.czm, 1, function(x){log(x) - mean(log(x))})) #centered log ratio transformation
# first only for filter; check the other samples and whether they need that

d.czmf_n <- cmultRepl(t(asvsf_n),  label=0, method="CZM", frac= 0.999, adjust=FALSE, z.warning = 0.999)
#Column no. 605 containing >99.9% zeros/unobserved values deleted (see arguments z.warning and z.delete).
## ans many more
# but normal, right? 
d.clrf_n <- (apply(d.czmf_n, 1, function(x){log(x) - mean(log(x))})) #centered log ratio transformation

#move all into positive values
asvs_normf_n <- apply(d.clrf_n, 2, function(x) x - min(x[x < 0]))

#make it relative
asvs_normf_n_P <- apply(asvs_normf_n, 2, function(x) {x/sum(x)}) 

# check later whether i need to clr it
asvs_normnp_n_P <- apply(asvsnp_n, 2, function(x) {x/sum(x)})


##
## creating phyloseq objects - for filters
## 

asvs_normf_n <- as.matrix(asvs_normf_n)
ncol(asvs_normf_n)
# 72

taxa_n <- as.matrix(taxa_n)
#taxa <- as.data.frame(taxa)
str(asvs_normf_n)

OTUf_n = otu_table(asvs_normf_n, taxa_are_rows = T)
ncol(OTUf_n) # 72
TAX_n = tax_table(taxa_n)
samples_n = sample_data(meta_samples_nf)


## do the names match? 
sample_names(OTUf_n)           
sample_names(samples_n)        

Next_18Sf <- phyloseq(OTUf_n, TAX_n, samples_n)
Next_18Sf
# 72 samples
# 8980 taxa



## sort from van Mijen to Rijfjorden, from outher to inner fjord
Next_18Sf@otu_table

Next_18Sf <- ps_reorder(Next_18Sf, sample_order = c("N.St1.F02.1", "N.St1.F02.2", "N.St1.F02.3",	"N.St1.F3.1",	"N.St1.F3.2",	"N.St1.F3.3",                  # 6
                                                      "N.St2.F02.1",	"N.St2.F02.2", "N.St2.F02.3",	"N.St2.F3.1",	"N.St2.F3.2",	"N.St2.F3.3",	               # 6
                                                      "N.St4.F02.1",	"N.St4.F02.2",	"N.St4.F02.3",	"N.St4.F3.1",	"N.St4.F3.2",	"N.St4.F3.3",              # 6
                                                      "N.St3.F02.1",	"N.St3.F02.2",	"N.St3.F02.3",	"N.St3.F3.1",	"N.St3.F3.2",	"N.St3.F3.3",              # 6
                                                      "N.St12.F02.1",	"N.St12.F02.2",	"N.St12.F02.3",	"N.St12.F3.1",	"N.St12.F3.2",	"N.St12.F3.3",         # 6	
                                                      "N.St13.F02.1",	"N.St13.F02.2",	"N.St13.F02.3",	"N.St13.F3.1",	"N.St13.F3.2",	"N.St13.F3.3",         # 6
                                                      "N.St11.F02.1",	"N.St11.F02.2",	"N.St11.F02.3",	"N.St11.F3.1",	"N.St11.F3.2",	"N.St11.F3.3.N",       # 6	
                                                      "N.St5.F02.1",	"N.St5.F02.2",	"N.St5.F02.3",	"N.St5.F3.1",	"N.St5.F3.2",	"N.St5.F3.3",              # 6	
                                                      "N.St7.F02.1",	"N.St7.F02.2",	"N.St7.F02.3",	"N.St7.F3.1",	"N.St7.F3.2.N",	"N.St7.F3.3.N",          # 6
                                                      "N.St6.F02.1",	"N.St6.F02.2",	"N.St6.F02.3",	"N.St6.F3.1",	"N.St6.F3.2",	"N.St6.F3.3",              # 6	
                                                      "N.St10.F02.1",	"N.St10.F02.2",	"N.St10.F02.3",	"N.St10.F3.1",	"N.St10.F3.2",	"N.St10.F3.3",         # 6	
                                                      "N.St8.F02.1",	"N.St8.F02.2",	"N.St8.F02.3",	"N.St8.F3.1",	"N.St8.F3.2",	"N.St8.F3.3.N"))           # 6 

Next_18Sf@otu_table

sample_data(Next_18Sf)$sample_order <- factor(sample_names(Next_18Sf))
Next_18Sf@sam_data


# for percentage

asvs_normf_n_P <- as.matrix(asvs_normf_n_P)
str(asvs_normf_n_P)

OTUfP = otu_table(asvs_normf_n_P, taxa_are_rows = T)

NextSeq_18SfP <- phyloseq(OTUfP, TAX_n, samples_n)
NextSeq_18SfP

# 72 samples w 8980 taxa

NextSeq_18SfP <- ps_reorder(NextSeq_18SfP, sample_order = c("N.St1.F02.1", "N.St1.F02.2", "N.St1.F02.3",	"N.St1.F3.1",	"N.St1.F3.2",	"N.St1.F3.3",                # 6
                                                            "N.St2.F02.1",	"N.St2.F02.2", "N.St2.F02.3",	"N.St2.F3.1",	"N.St2.F3.2",	"N.St2.F3.3",	               # 6
                                                            "N.St4.F02.1",	"N.St4.F02.2",	"N.St4.F02.3",	"N.St4.F3.1",	"N.St4.F3.2",	"N.St4.F3.3",              # 6
                                                            "N.St3.F02.1",	"N.St3.F02.2",	"N.St3.F02.3",	"N.St3.F3.1",	"N.St3.F3.2",	"N.St3.F3.3",              # 6
                                                            "N.St12.F02.1",	"N.St12.F02.2",	"N.St12.F02.3",	"N.St12.F3.1",	"N.St12.F3.2",	"N.St12.F3.3",         # 6	
                                                            "N.St13.F02.1",	"N.St13.F02.2",	"N.St13.F02.3",	"N.St13.F3.1",	"N.St13.F3.2",	"N.St13.F3.3",         # 6
                                                            "N.St11.F02.1",	"N.St11.F02.2",	"N.St11.F02.3",	"N.St11.F3.1",	"N.St11.F3.2",	"N.St11.F3.3.N",       # 6	
                                                            "N.St5.F02.1",	"N.St5.F02.2",	"N.St5.F02.3",	"N.St5.F3.1",	"N.St5.F3.2",	"N.St5.F3.3",              # 6	
                                                            "N.St7.F02.1",	"N.St7.F02.2",	"N.St7.F02.3",	"N.St7.F3.1",	"N.St7.F3.2.N",	"N.St7.F3.3.N",          # 6
                                                            "N.St6.F02.1",	"N.St6.F02.2",	"N.St6.F02.3",	"N.St6.F3.1",	"N.St6.F3.2",	"N.St6.F3.3",              # 6	
                                                            "N.St10.F02.1",	"N.St10.F02.2",	"N.St10.F02.3",	"N.St10.F3.1",	"N.St10.F3.2",	"N.St10.F3.3",         # 6	
                                                            "N.St8.F02.1",	"N.St8.F02.2",	"N.St8.F02.3",	"N.St8.F3.1",	"N.St8.F3.2",	"N.St8.F3.3.N"))           # 6 

NextSeq_18SfP@otu_table

sample_data(NextSeq_18SfP)$sample_order <- factor(sample_names(NextSeq_18SfP))
NextSeq_18SfP@sam_data


## --  NP + deep DNA

asvsnp_n <- as.matrix(asvsnp_n)
str(asvsnp_n)

OTUnp = otu_table(asvsnp_n, taxa_are_rows = T)
TAX = tax_table(taxa_n)
samples_npd = sample_data(meta_samples_npd)

## do the names match? 
sample_names(OTUnp)           # all - are now a . 
sample_names(samples_npd)        

NextSeq_18Snp <- phyloseq(OTUnp, TAX, samples_npd)
NextSeq_18Snp
# 61 samples with 10425 samples

sample_names(NextSeq_18Snp)

NextSeq_18Snp <- ps_reorder(NextSeq_18Snp, sample_order = c("N.St1.dDNA",	"N.St1.N",	"N.St1.P1", "N.St1.P3",                          # 4
                                                        "N.St2.dDNA",	"N.St2.N",	"N.St2.P.1", "N.St2.P.2", "N.St2.P.3",               # 5
                                                        "N.St4.dDNA",	"N.St4.N",	"N.St4.P.1", "N.St4.P.2", "N.St4.P.3",               # 5
                                                        "N.St3.dDNA",	"N.St3.N",	"N.St3.P.1", "N.St3.P.2", "N.St3.P.3",               # 5
                                                        "N.St12.dDNA", "N.St12.N",	"N.St12.P.1", "N.St12.P.2", "N.St12.P.3",          # 5
                                                        "N.St13.dDNA", "N.St13.N",	"N.St13.P.1", "N.St13.P.2", "N.St13.P.3",          # 5
                                                        "N.St11.dDNA", "N.St11.N",	"N.St11.P.1", "N.St11.P.2", "N.St11.P.3",          # 5
                                                        "N.St5.dDNA",	"N.St5.N", "N.St5.N.150", "N.St5.P.1", "N.St5.P.2", "N.St5.P.3", # 6
                                                        "N.St7.dDNA",	"N.St7.N", "N.St7.P.1", "N.St7.P.2", "N.St7.P.3",                # 5
                                                        "N.St6.dDNA",	"N.St6.N", "N.St6.P.1", "N.St6.P.2", "N.St6.P.3",                # 5
                                                        "N.St10.dDNA",	"N.St10.N", "N.St10.P.1", "N.St10.P.2", "N.St10.P.3",          # 5
                                                        "N.St8.dDNA",	"N.St8.N", "N.St8.P.1", "N.St8.P.2", "N.St8.P.3",
                                                        "N.St9.dDNA"))               # 5

NextSeq_18Snp@otu_table # right order
NextSeq_18Snp@sam_data

sample_data(NextSeq_18Snp)$sample_order <- factor(sample_names(NextSeq_18Snp))
NextSeq_18Snp@sam_data



# percentage npdeppDNA

asvs_normnp_n_P <- as.matrix(asvs_normnp_n_P)
str(asvs_normnp_n_P)

OTUnpdP = otu_table(asvs_normnp_n_P, taxa_are_rows = T)

NextSeq_18SnpP <- phyloseq(OTUnpdP, TAX_n, samples_npd)
NextSeq_18SnpP
# 61 samples and 10425 taxa

NextSeq_18SnpP <- ps_reorder(NextSeq_18SnpP, sample_order = c("N.St1.dDNA",	"N.St1.N",	"N.St1.P1", "N.St1.P3",                              # 4
                                                              "N.St2.dDNA",	"N.St2.N",	"N.St2.P.1", "N.St2.P.2", "N.St2.P.3",               # 5
                                                              "N.St4.dDNA",	"N.St4.N",	"N.St4.P.1", "N.St4.P.2", "N.St4.P.3",               # 5
                                                              "N.St3.dDNA",	"N.St3.N",	"N.St3.P.1", "N.St3.P.2", "N.St3.P.3",               # 5
                                                              "N.St12.dDNA", "N.St12.N",	"N.St12.P.1", "N.St12.P.2", "N.St12.P.3",          # 5
                                                              "N.St13.dDNA", "N.St13.N",	"N.St13.P.1", "N.St13.P.2", "N.St13.P.3",          # 5
                                                              "N.St11.dDNA", "N.St11.N",	"N.St11.P.1", "N.St11.P.2", "N.St11.P.3",          # 5
                                                              "N.St5.dDNA",	"N.St5.N", "N.St5.N.150", "N.St5.P.1", "N.St5.P.2", "N.St5.P.3", # 6
                                                              "N.St7.dDNA",	"N.St7.N", "N.St7.P.1", "N.St7.P.2", "N.St7.P.3",                # 5
                                                              "N.St6.dDNA",	"N.St6.N", "N.St6.P.1", "N.St6.P.2", "N.St6.P.3",                # 5
                                                              "N.St10.dDNA",	"N.St10.N", "N.St10.P.1", "N.St10.P.2", "N.St10.P.3",          # 5
                                                              "N.St8.dDNA",	"N.St8.N", "N.St8.P.1", "N.St8.P.2", "N.St8.P.3"))               # 5


NextSeq_18SnpP@otu_table # right order
NextSeq_18SnpP@sam_data

sample_data(NextSeq_18SnpP)$sample_order <- factor(sample_names(NextSeq_18SnpP))
NextSeq_18SnpP@sam_data


### ----- 
## plots

# filters

## try to merge and only look at phylum
nf.phy = tax_glom(Next_18Sf, taxrank="Phylum")

Next_18Sf@sam_data
nf.phy@sam_data

sample_data(nf.phy)$NewID <- factor(sample_names(nf.phy))
sample_data(nf.phy)$NewID <- factor(sample_data(nf.phy)$NewID, levels = (sample_data(nf.phy)$sample_order))

str(nf.phy@sam_data)


png("NextSeq_filter_phylum_try.png", width = 1500, height = 800)
plot_bar(nf.phy, x = "NewID", fill = "Phylum") + labs(title = "NextSeq 18S - Phylum per filter sample - CLR transformed") + xlab("Sample") +
  theme_light() + scale_fill_viridis(discrete = T, option = "C") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                                         axis.text=element_text(size=14), #change font size of axis text
                                                                         axis.title=element_text(size=19), #change font size of axis titles
                                                                         plot.title=element_text(size=26), #change font size of plot title
                                                                         legend.text=element_text(size=16), #change font size of legend text
                                                                         legend.title=element_text(size=20), #change font size of legend title   )
                                                                         legend.key.size = unit(1, "line"))
dev.off()


## try to merge and only look at class
n18.gen.class = tax_glom(Next_18Sf, taxrank="Class", NArm = F)
n18.gen.class@tax_table

Next_18Sf@tax_table

sample_data(n18.gen.class)$NewID <- factor(sample_names(n18.gen.class))
sample_data(n18.gen.class)$NewID <- factor(sample_data(n18.gen.class)$NewID, levels = (sample_data(n18.gen.class)$sample_order))

TOPS18 = names(sort(taxa_sums(n18.gen.class), TRUE)[1:21]) # change 1:X till I have top20
TOPS18
TOPS118 = prune_taxa(TOPS18, n18.gen.class)
TOPS118 # 20 taxa, 72 samples

png("NextSeq_18S_class_asvs_top20.png", width = 1500, height = 800)
plot_bar(TOPS118, x = "NewID", fill = "Class") + labs(title = "NextSeq 18S - Class per filter sample - Top 20 - CLR transformed") + xlab("Sample") +
  theme_light() + scale_fill_viridis(discrete = T, option = "C") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                                         axis.text=element_text(size=14), #change font size of axis text
                                                                         axis.title=element_text(size=19), #change font size of axis titles
                                                                         plot.title=element_text(size=26), #change font size of plot title
                                                                         legend.text=element_text(size=16), #change font size of legend text
                                                                         legend.title=element_text(size=20)) #change font size of legend title   )
dev.off()

#try2 - no nas in plot
TOPS1_218 <- subset_taxa(TOPS118, !is.na(Class) & !Class %in% c("", "uncharacterized"))


png("NextSeq_18S_class_asvs_top20_noNAs_rainbow.png", width = 1500, height = 800)
plot_bar(TOPS1_218, !is.na("Class"), x = "NewID", y = "Abundance", fill = "Class") + labs(title = "NextSeq 18S - Class per filter sample - Top 20 - CLR transformed") + xlab("Sample") +
  theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                        axis.text=element_text(size=14), #change font size of axis text
                        axis.title=element_text(size=19), #change font size of axis titles
                        plot.title=element_text(size=26), #change font size of plot title
                        legend.text=element_text(size=16), #change font size of legend text
                        legend.title=element_text(size=20)) #change font size of legend title   )
dev.off()


rhg_cols <- c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F", "#8CD17D", "#B6992D", "#F1CE63", "#499894", "#86BCB6",
                       "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", "#D37295", "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6")
                       
png("NextSeq_18S_class_asvs_top20_noNAs_other_colour.png", width = 1500, height = 800)
plot_bar(TOPS1_218, !is.na("Class"), x = "NewID", y = "Abundance", fill = "Class") + labs(title = "NextSeq 18S - Class per filter sample - Top 20 - CLR transformed") + xlab("Sample") +
  theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                        axis.text=element_text(size=14), #change font size of axis text
                        axis.title=element_text(size=19), #change font size of axis titles
                        plot.title=element_text(size=26), #change font size of plot title
                        legend.text=element_text(size=16), #change font size of legend text
                        legend.title=element_text(size=20)) +  #change font size of legend title   ) 
  scale_fill_manual(values = rhg_cols)
dev.off()

# only dinos, diatoms + haptos
bacdinoprym18 <- subset_taxa(TOPS1_218, Class %in% c("Bacillariophyceae", "Dinophyceae", "Prymnesiophyceae"))

png("NextSeq_18S_class_only_DiatDinosPrymn.png", width = 1500, height = 800)
plot_bar(bacdinoprym18, x = "NewID", y = "Abundance", fill = "Class") + labs(title = "NextSeq 18S - Bacillariophyceae, Dinophyceae, Prymnesiophyceae (Class)- CLR transformed") + xlab("Sample") +
  theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                        axis.text=element_text(size=14), #change font size of axis text
                        axis.title=element_text(size=19), #change font size of axis titles
                        plot.title=element_text(size=26), #change font size of plot title
                        legend.text=element_text(size=16), #change font size of legend text
                        legend.title=element_text(size=20)) + #change font size of legend title   )
  scale_fill_manual(values = c("#3e8fff", "#ff6cd0", "#ff2a58"))
dev.off()


# -----

# create a phyloseq w only dinos, diatoms + prymnesio? 

tax_try <- read.csv("NextSeq_18S_taxa.csv")     

tax_try <- tax_try %>% tibble::column_to_rownames("otu")
TAX_n_diadinoprym <- subset(tax_try, Class == "Bacillariophyceae" | Class == "Dinophyceae" | Class == "Prymnesiophyceae")



TAX_n_diadinoprym <- as.matrix(TAX_n_diadinoprym)
TAX_n_diadinoprym_phy = tax_table(TAX_n_diadinoprym)

phylo_dinodiaprym <- phyloseq(OTUf_n, TAX_n_diadinoprym_phy, samples_n)
phylo_dinodiaprym
# 72 samples, 3820 taxa
phylo_dinodiaprym@sam_data

phylo_dinodiaprym <- ps_reorder(phylo_dinodiaprym, sample_order = c("N.St1.F02.1", "N.St1.F02.2", "N.St1.F02.3",	"N.St1.F3.1",	"N.St1.F3.2",	"N.St1.F3.3",                # 6
                                                            "N.St2.F02.1",	"N.St2.F02.2", "N.St2.F02.3",	"N.St2.F3.1",	"N.St2.F3.2",	"N.St2.F3.3",	               # 6
                                                            "N.St4.F02.1",	"N.St4.F02.2",	"N.St4.F02.3",	"N.St4.F3.1",	"N.St4.F3.2",	"N.St4.F3.3",              # 6
                                                            "N.St3.F02.1",	"N.St3.F02.2",	"N.St3.F02.3",	"N.St3.F3.1",	"N.St3.F3.2",	"N.St3.F3.3",              # 6
                                                            "N.St12.F02.1",	"N.St12.F02.2",	"N.St12.F02.3",	"N.St12.F3.1",	"N.St12.F3.2",	"N.St12.F3.3",         # 6	
                                                            "N.St13.F02.1",	"N.St13.F02.2",	"N.St13.F02.3",	"N.St13.F3.1",	"N.St13.F3.2",	"N.St13.F3.3",         # 6
                                                            "N.St11.F02.1",	"N.St11.F02.2",	"N.St11.F02.3",	"N.St11.F3.1",	"N.St11.F3.2",	"N.St11.F3.3.N",       # 6	
                                                            "N.St5.F02.1",	"N.St5.F02.2",	"N.St5.F02.3",	"N.St5.F3.1",	"N.St5.F3.2",	"N.St5.F3.3",              # 6	
                                                            "N.St7.F02.1",	"N.St7.F02.2",	"N.St7.F02.3",	"N.St7.F3.1",	"N.St7.F3.2.N",	"N.St7.F3.3.N",          # 6
                                                            "N.St6.F02.1",	"N.St6.F02.2",	"N.St6.F02.3",	"N.St6.F3.1",	"N.St6.F3.2",	"N.St6.F3.3",              # 6	
                                                            "N.St10.F02.1",	"N.St10.F02.2",	"N.St10.F02.3",	"N.St10.F3.1",	"N.St10.F3.2",	"N.St10.F3.3",         # 6	
                                                            "N.St8.F02.1",	"N.St8.F02.2",	"N.St8.F02.3",	"N.St8.F3.1",	"N.St8.F3.2",	"N.St8.F3.3.N"))           # 6 


sample_data(phylo_dinodiaprym)$sample_order <- factor(sample_names(phylo_dinodiaprym))
phylo_dinodiaprym@sam_data

sample_data(phylo_dinodiaprym)$NewID <- factor(sample_names(phylo_dinodiaprym))
sample_data(phylo_dinodiaprym)$NewID <- factor(sample_data(phylo_dinodiaprym)$NewID, levels = (sample_data(phylo_dinodiaprym)$sample_order))

phylo_dinodiaprym@tax_table

phylo_dinodiaprym.glom = tax_glom(phylo_dinodiaprym, taxrank="Order")

png("NextSeq_18S_class_only_DiatDinosPrymn_try.png", width = 1500, height = 800)
plot_bar(phylo_dinodiaprym.glom, x = "NewID", fill = "Order") + labs(title = "NextSeq 18S - Bacillariophyceae, Dinophyceae, Prymnesiophyceae (Class) - CLR transformed") + xlab("Sample") +
  theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                        axis.text=element_text(size=14), #change font size of axis text
                        axis.title=element_text(size=19), #change font size of axis titles
                        plot.title=element_text(size=26), #change font size of plot title
                        legend.text=element_text(size=16), #change font size of legend text
                        legend.title=element_text(size=20))
dev.off()


# only diatoms

tax_try <- read.csv("NextSeq_18S_taxa.csv")     

tax_try <- tax_try %>% tibble::column_to_rownames("otu")
TAX_n_dia <- subset(tax_try, Class == "Bacillariophyceae")



TAX_n_dia <- as.matrix(TAX_n_dia)
TAX_n_dia_phy = tax_table(TAX_n_dia)

phylo_dia <- phyloseq(OTUf_n, TAX_n_dia_phy, samples_n)
phylo_dia
# 72 samples, 75 taxa
phylo_dia@sam_data
phylo_dia@tax_table

phylo_dia <- ps_reorder(phylo_dia, sample_order = c("N.St1.F02.1", "N.St1.F02.2", "N.St1.F02.3",	"N.St1.F3.1",	"N.St1.F3.2",	"N.St1.F3.3",                # 6
                                                                    "N.St2.F02.1",	"N.St2.F02.2", "N.St2.F02.3",	"N.St2.F3.1",	"N.St2.F3.2",	"N.St2.F3.3",	               # 6
                                                                    "N.St4.F02.1",	"N.St4.F02.2",	"N.St4.F02.3",	"N.St4.F3.1",	"N.St4.F3.2",	"N.St4.F3.3",              # 6
                                                                    "N.St3.F02.1",	"N.St3.F02.2",	"N.St3.F02.3",	"N.St3.F3.1",	"N.St3.F3.2",	"N.St3.F3.3",              # 6
                                                                    "N.St12.F02.1",	"N.St12.F02.2",	"N.St12.F02.3",	"N.St12.F3.1",	"N.St12.F3.2",	"N.St12.F3.3",         # 6	
                                                                    "N.St13.F02.1",	"N.St13.F02.2",	"N.St13.F02.3",	"N.St13.F3.1",	"N.St13.F3.2",	"N.St13.F3.3",         # 6
                                                                    "N.St11.F02.1",	"N.St11.F02.2",	"N.St11.F02.3",	"N.St11.F3.1",	"N.St11.F3.2",	"N.St11.F3.3.N",       # 6	
                                                                    "N.St5.F02.1",	"N.St5.F02.2",	"N.St5.F02.3",	"N.St5.F3.1",	"N.St5.F3.2",	"N.St5.F3.3",              # 6	
                                                                    "N.St7.F02.1",	"N.St7.F02.2",	"N.St7.F02.3",	"N.St7.F3.1",	"N.St7.F3.2.N",	"N.St7.F3.3.N",          # 6
                                                                    "N.St6.F02.1",	"N.St6.F02.2",	"N.St6.F02.3",	"N.St6.F3.1",	"N.St6.F3.2",	"N.St6.F3.3",              # 6	
                                                                    "N.St10.F02.1",	"N.St10.F02.2",	"N.St10.F02.3",	"N.St10.F3.1",	"N.St10.F3.2",	"N.St10.F3.3",         # 6	
                                                                    "N.St8.F02.1",	"N.St8.F02.2",	"N.St8.F02.3",	"N.St8.F3.1",	"N.St8.F3.2",	"N.St8.F3.3.N"))           # 6 


sample_data(phylo_dia)$sample_order <- factor(sample_names(phylo_dia))

sample_data(phylo_dia)$NewID <- factor(sample_names(phylo_dia))
sample_data(phylo_dia)$NewID <- factor(sample_data(phylo_dia)$NewID, levels = (sample_data(phylo_dia)$sample_order))


phylo_dia.glom = tax_glom(phylo_dia, taxrank="Genus")


rhg_cols_15 <- c("#231942", "#5e548e", "#9f86c0", "#be95c4", "#e0b1cb", "#22577a", "#38a3a5", "#57cc99",
                       "#80ed99", "#c7f9cc", "#004949", "#009292", "#FF6DB6", "#5fa8d3", "#2ec4b6")
                       

png("NextSeq_18S_genus_diatoms_try.png", width = 1500, height = 800)
plot_bar(phylo_dia.glom, x = "NewID", fill = "Genus") + labs(title = "NextSeq 18S - Bacillariophyceae - Top 15 Genera - CLR transformed") + xlab("Sample") +
  theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                        axis.text=element_text(size=14), #change font size of axis text
                        axis.title=element_text(size=19), #change font size of axis titles
                        plot.title=element_text(size=26), #change font size of plot title
                        legend.text=element_text(size=16), #change font size of legend text
                        legend.title=element_text(size=20)) + scale_fill_manual(values = rhg_cols_15)
dev.off()


# only dinos

tax_try <- read.csv("NextSeq_18S_taxa.csv")     

tax_try <- tax_try %>% tibble::column_to_rownames("otu")
TAX_n_dino <- subset(tax_try, Class == "Dinophyceae")



TAX_n_dino <- as.matrix(TAX_n_dino)
TAX_n_dino_phy = tax_table(TAX_n_dino)

phylo_dino <- phyloseq(OTUf_n, TAX_n_dino_phy, samples_n)
phylo_dino
# 72 samples, 3467 taxa


phylo_dino <- ps_reorder(phylo_dino, sample_order = c("N.St1.F02.1", "N.St1.F02.2", "N.St1.F02.3",	"N.St1.F3.1",	"N.St1.F3.2",	"N.St1.F3.3",                # 6
                                                    "N.St2.F02.1",	"N.St2.F02.2", "N.St2.F02.3",	"N.St2.F3.1",	"N.St2.F3.2",	"N.St2.F3.3",	               # 6
                                                    "N.St4.F02.1",	"N.St4.F02.2",	"N.St4.F02.3",	"N.St4.F3.1",	"N.St4.F3.2",	"N.St4.F3.3",              # 6
                                                    "N.St3.F02.1",	"N.St3.F02.2",	"N.St3.F02.3",	"N.St3.F3.1",	"N.St3.F3.2",	"N.St3.F3.3",              # 6
                                                    "N.St12.F02.1",	"N.St12.F02.2",	"N.St12.F02.3",	"N.St12.F3.1",	"N.St12.F3.2",	"N.St12.F3.3",         # 6	
                                                    "N.St13.F02.1",	"N.St13.F02.2",	"N.St13.F02.3",	"N.St13.F3.1",	"N.St13.F3.2",	"N.St13.F3.3",         # 6
                                                    "N.St11.F02.1",	"N.St11.F02.2",	"N.St11.F02.3",	"N.St11.F3.1",	"N.St11.F3.2",	"N.St11.F3.3.N",       # 6	
                                                    "N.St5.F02.1",	"N.St5.F02.2",	"N.St5.F02.3",	"N.St5.F3.1",	"N.St5.F3.2",	"N.St5.F3.3",              # 6	
                                                    "N.St7.F02.1",	"N.St7.F02.2",	"N.St7.F02.3",	"N.St7.F3.1",	"N.St7.F3.2.N",	"N.St7.F3.3.N",          # 6
                                                    "N.St6.F02.1",	"N.St6.F02.2",	"N.St6.F02.3",	"N.St6.F3.1",	"N.St6.F3.2",	"N.St6.F3.3",              # 6	
                                                    "N.St10.F02.1",	"N.St10.F02.2",	"N.St10.F02.3",	"N.St10.F3.1",	"N.St10.F3.2",	"N.St10.F3.3",         # 6	
                                                    "N.St8.F02.1",	"N.St8.F02.2",	"N.St8.F02.3",	"N.St8.F3.1",	"N.St8.F3.2",	"N.St8.F3.3.N"))           # 6 


sample_data(phylo_dino)$sample_order <- factor(sample_names(phylo_dino))

sample_data(phylo_dino)$NewID <- factor(sample_names(phylo_dino))
sample_data(phylo_dino)$NewID <- factor(sample_data(phylo_dino)$NewID, levels = (sample_data(phylo_dino)$sample_order))

phylo_dino.glom = tax_glom(phylo_dino, taxrank="Genus")

TOPSdino = names(sort(taxa_sums(phylo_dino.glom), TRUE)[1:15]) # change 1:X till I have top20
TOPSdino
TOPSdinoss = prune_taxa(TOPSdino, phylo_dino.glom)
TOPSdinoss # 15 taxa, 72 samples


png("NextSeq_18S_genus_dinos_try.png", width = 1500, height = 800)
plot_bar(TOPSdinoss, x = "NewID", fill = "Genus") + labs(title = "NextSeq 18S - Dinophyceae - Top 15 Genera - CLR transformed") + xlab("Sample") +
  theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                        axis.text=element_text(size=14), #change font size of axis text
                        axis.title=element_text(size=19), #change font size of axis titles
                        plot.title=element_text(size=26), #change font size of plot title
                        legend.text=element_text(size=16), #change font size of legend text
                        legend.title=element_text(size=20))  + scale_fill_manual(values = rhg_cols_15)
dev.off()


# only haptos

tax_try <- read.csv("NextSeq_18S_taxa.csv")     

tax_try <- tax_try %>% tibble::column_to_rownames("otu")
TAX_n_hapto <- subset(tax_try, Phylum == "Haptophyta")



TAX_n_hapto <- as.matrix(TAX_n_hapto)
TAX_n_hapto_phy = tax_table(TAX_n_hapto)

phylo_hapto <- phyloseq(OTUf_n, TAX_n_hapto_phy, samples_n)
phylo_hapto
# 72 samples, 300 taxa


phylo_hapto <- ps_reorder(phylo_hapto, sample_order = c("N.St1.F02.1", "N.St1.F02.2", "N.St1.F02.3",	"N.St1.F3.1",	"N.St1.F3.2",	"N.St1.F3.3",                # 6
                                                      "N.St2.F02.1",	"N.St2.F02.2", "N.St2.F02.3",	"N.St2.F3.1",	"N.St2.F3.2",	"N.St2.F3.3",	               # 6
                                                      "N.St4.F02.1",	"N.St4.F02.2",	"N.St4.F02.3",	"N.St4.F3.1",	"N.St4.F3.2",	"N.St4.F3.3",              # 6
                                                      "N.St3.F02.1",	"N.St3.F02.2",	"N.St3.F02.3",	"N.St3.F3.1",	"N.St3.F3.2",	"N.St3.F3.3",              # 6
                                                      "N.St12.F02.1",	"N.St12.F02.2",	"N.St12.F02.3",	"N.St12.F3.1",	"N.St12.F3.2",	"N.St12.F3.3",         # 6	
                                                      "N.St13.F02.1",	"N.St13.F02.2",	"N.St13.F02.3",	"N.St13.F3.1",	"N.St13.F3.2",	"N.St13.F3.3",         # 6
                                                      "N.St11.F02.1",	"N.St11.F02.2",	"N.St11.F02.3",	"N.St11.F3.1",	"N.St11.F3.2",	"N.St11.F3.3.N",       # 6	
                                                      "N.St5.F02.1",	"N.St5.F02.2",	"N.St5.F02.3",	"N.St5.F3.1",	"N.St5.F3.2",	"N.St5.F3.3",              # 6	
                                                      "N.St7.F02.1",	"N.St7.F02.2",	"N.St7.F02.3",	"N.St7.F3.1",	"N.St7.F3.2.N",	"N.St7.F3.3.N",          # 6
                                                      "N.St6.F02.1",	"N.St6.F02.2",	"N.St6.F02.3",	"N.St6.F3.1",	"N.St6.F3.2",	"N.St6.F3.3",              # 6	
                                                      "N.St10.F02.1",	"N.St10.F02.2",	"N.St10.F02.3",	"N.St10.F3.1",	"N.St10.F3.2",	"N.St10.F3.3",         # 6	
                                                      "N.St8.F02.1",	"N.St8.F02.2",	"N.St8.F02.3",	"N.St8.F3.1",	"N.St8.F3.2",	"N.St8.F3.3.N"))           # 6 


sample_data(phylo_hapto)$sample_order <- factor(sample_names(phylo_hapto))

sample_data(phylo_hapto)$NewID <- factor(sample_names(phylo_hapto))
sample_data(phylo_hapto)$NewID <- factor(sample_data(phylo_hapto)$NewID, levels = (sample_data(phylo_hapto)$sample_order))
phylo_hapto@sam_data

phylo_hapto.glom = tax_glom(phylo_hapto, taxrank="Genus")

TOPShapto = names(sort(taxa_sums(phylo_hapto.glom), TRUE)[1:15]) # change 1:X till I have top20
TOPShapto
TOPShaptoss = prune_taxa(TOPShapto, phylo_hapto.glom)
TOPShaptoss # 15 taxa, 72 samples

png("NextSeq_18S_genus_hapto_try.png", width = 1500, height = 800)
plot_bar(TOPShaptoss, x = "NewID", fill = "Genus") + labs(title = "NextSeq 18S - Haptophyta - Top 15 Genera - CLR transformed") + xlab("Sample") +
  theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                        axis.text=element_text(size=14), #change font size of axis text
                        axis.title=element_text(size=19), #change font size of axis titles
                        plot.title=element_text(size=26), #change font size of plot title
                        legend.text=element_text(size=16), #change font size of legend text
                        legend.title=element_text(size=20))  + scale_fill_manual(values = rhg_cols_15)
dev.off()


## --- 


nfP.phy = tax_glom(NextSeq_18SfP, taxrank="Phylum", NArm=FALSE)

nfP.phy@sam_data

sample_data(nfP.phy)$NewID <- factor(sample_names(nfP.phy))
sample_data(nfP.phy)$NewID <- factor(sample_data(nfP.phy)$NewID, levels = (sample_data(nfP.phy)$sample_order))

nfP.phy@sam_data

png("NextSeq_filter_phylum_try_P.png", width = 1500, height = 800)
plot_bar(nfP.phy, x = "NewID", fill = "Phylum") + labs(title = "NextSeq 18S - Phylum per filter sample - Relative - CLR transformed") + xlab("Sample") + 
  theme_light() + scale_fill_viridis(discrete = T, option = "C") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                                          axis.text=element_text(size=14), #change font size of axis text
                                                                          axis.title=element_text(size=19), #change font size of axis titles
                                                                          plot.title=element_text(size=26), #change font size of plot title
                                                                          legend.text=element_text(size=16), #change font size of legend text
                                                                          legend.title=element_text(size=20)) #change font size of legend title   )
dev.off()

# net/pump

mnp.phy = tax_glom(NextSeq_18Snp, taxrank="Phylum", NArm=FALSE)

mnp.phy@sam_data

sample_data(mnp.phy)$NewID <- factor(sample_names(mnp.phy))
sample_data(mnp.phy)$NewID <- factor(sample_data(mnp.phy)$NewID, levels = (sample_data(mnp.phy)$sample_order))

mnp.phy@sam_data

png("NextSeq_netpump_phylum_try.png", width = 1500, height = 800)
plot_bar(mnp.phy, x = "NewID", fill = "Phylum") + labs(title = "NextSeq 18S - Phylum per net/pump sample") + xlab("Sample") +
  theme_light() + scale_fill_viridis(discrete = T, option = "C") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                                         axis.text=element_text(size=14), #change font size of axis text
                                                                         axis.title=element_text(size=19), #change font size of axis titles
                                                                         plot.title=element_text(size=26), #change font size of plot title
                                                                         legend.text=element_text(size=16), #change font size of legend text
                                                                         legend.title=element_text(size=20)) #change font size of legend title   )
dev.off()

# net/pump %

mnpP.phy = tax_glom(NextSeq_18SnpP, taxrank="Phylum", NArm=FALSE)

mnpP.phy@sam_data

sample_data(mnpP.phy)$NewID <- factor(sample_names(mnpP.phy))
sample_data(mnpP.phy)$NewID <- factor(sample_data(mnpP.phy)$NewID, levels = (sample_data(mnpP.phy)$sample_order))

mnpP.phy@sam_data

png("NextSeq_netpump_phylum_try_P.png", width = 1500, height = 800)
plot_bar(mnpP.phy, x = "NewID", fill = "Phylum") + labs(title = "NextSeq 18S - Phylum per net/pump sample - Relative") + xlab("Sample") +
  theme_light() + scale_fill_viridis(discrete = T, option = "C") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                                          axis.text=element_text(size=14), #change font size of axis text
                                                                          axis.title=element_text(size=19), #change font size of axis titles
                                                                          plot.title=element_text(size=26), #change font size of plot title
                                                                          legend.text=element_text(size=16), #change font size of legend text
                                                                          legend.title=element_text(size=20)) #change font size of legend title   )
## + geom_bar(aes(fill=Phylum), stat="identity", position="stack") 
## for a smooth overlay of colours
dev.off()





# ---------------------
# PCA

### ordination
# PCA
# beta diversity !!

# filter samples !!

Next_18Sf@sam_data
# create a new phyloseq wo/ clr transf before

asvsf_n <- as.matrix(asvsf_n)

taxa_n <- as.matrix(taxa_n)
#taxa <- as.data.frame(taxa)
str(asvsf_n)

OTUfn_n = otu_table(asvsf_n, taxa_are_rows = T)
TAXfn_n = tax_table(taxa_n)
samplesfn_n = sample_data(meta_samples_nf)

## do the names match? 
sample_names(OTUfn_n)           # all - are now a . 
sample_names(samplesfn_n)        

Next_18Sf_default <- phyloseq(OTUfn_n, TAXfn_n, samplesfn_n)
Next_18Sf_default
# 72 samples
# 10425 taxa

Next_18Sf_default@tax_table
# some NAs

tax_fix_interactive(Next_18Sf_default)

NextSeq_fixed <- Next_18Sf_default %>%
  tax_fix(
    min_length = 4,
    unknowns = c("NA"),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified")


str(NextSeq_fixed@sam_data)
NextSeq_fixed@sam_data

NextSeq_fixed@sam_data$stat <- as.factor(NextSeq_fixed@sam_data$stat)
Next_18Sf_default@sam_data$stat <- as.factor(Next_18Sf_default@sam_data$stat)

# add fjords as a variable!!

## asvs - the nicest ones cause lot of info
test1 <- NextSeq_fixed %>% 
  tax_transform("clr", rank ="unique") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "fjord", size = 3) + theme_bw() +
  scale_colour_manual(values = c("van Mijenfjorden" = "#ff00ad", "Kongsfjorden" = "#ab55c7", "Wijdefjorden" = "#5e91e4", "Rijpfjorden"= "#0edaff"), 
                      breaks = c('van Mijenfjorden', 'Kongsfjorden', 'Wijdefjorden', "Rijpfjorden")) +
  labs(title = "PCA Analysis - NextSeq 18S ASVs", color = "Fjord") + theme(legend.text=element_text(size=10), 
                                                                           legend.title=element_text(size=12),
                                                                           plot.title=element_text(size=16)) +
  stat_ellipse(aes(colour = fjord))

ggsave("pca_NextSeq_18S_try_asvs1_try2.png", width = 6.5, height = 6)

subset_ord_plot(test1) # ??? to get the pca values? 

# axis 1 + 3
NextSeq_fixed %>% tax_transform("clr", rank ="unique") %>%
  ord_calc() %>% 
  ord_plot(axes = c(1,3), color = "fjord", size = 3) + theme_bw() +
  scale_colour_manual(values = c("van Mijenfjorden" = "#ff00ad", "Kongsfjorden" = "#ab55c7", "Wijdefjorden" = "#5e91e4", "Rijpfjorden"= "#0edaff"), 
                      breaks = c('van Mijenfjorden', 'Kongsfjorden', 'Wijdefjorden', "Rijpfjorden")) +
  labs(title = "PCA Analysis - NextSeq 18S ASVs - Axis 1 + 3", color = "Fjord") + theme(legend.text=element_text(size=10), 
                                                                                      legend.title=element_text(size=12),
                                                                                      plot.title=element_text(size=16)) +
  stat_ellipse(aes(colour = fjord))



ggsave("pca_NextSeq_18S_try_asvs_2.png", width = 6.5, height = 6)



NextSeq_18Sf_default %>%
  tax_transform("identity", rank = "unique") %>% # don't transform!
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_get() %>%
  phyloseq::plot_scree() + theme(axis.text.x = element_text(size = 6)) + labs(title = "NextSeq 18S - PCA eigenvalues")

ggsave("pca_NextSeq_18S_axes_asvs.png", width = 10, height= 7)


#stats



## RDA

NextSeq_fixed %>% 
  tax_transform("clr", rank = "unique") %>%
  ord_calc(constraints = c("temp", "sal", "dens", "oxyg", "chla", "lat", "long"),
           scale_cc = F, method = "RDA") %>%
  ord_plot(auto_caption = NA, color = "fjord", size = 4.5, plot_taxa = 1:10, constraint_lab_style = constraint_lab_style(size = 4), tax_lab_style = tax_lab_style(size = 4)) + theme_bw() + 
  scale_colour_manual(values = c("van Mijenfjorden" = "#ff00ad", "Kongsfjorden" = "#ab55c7", "Wijdefjorden" = "#5e91e4", "Rijpfjorden"= "#0edaff"), 
                      breaks = c('van Mijenfjorden', 'Kongsfjorden', 'Wijdefjorden', "Rijpfjorden")) +
  labs(title = "RDA Analysis - NextSeq 18S - Top 10 ASVs", color = "Fjord") + theme(legend.text=element_text(size=12), 
                                                                                    legend.title=element_text(size=16),
                                                                                    plot.title=element_text(size=20),
                                                                                    axis.title = element_text(size=14)) 
ggsave("RDA_NextSeq_18S_try_asvs.png", width = 13, height = 8)

#stats

aa <- NextSeq_fixed %>% 
  #tax_transform("clr", rank = "unique") %>%
  dist_calc(dist = "aitchison")

set.seed(111)
PERM <- aa%>%dist_permanova(seed = 1, variables = c("temp", "sal", "dens", "oxyg", "chla", "lat", "long"), n_processes = 1, n_perms = 999)
PERM

PERM2 <- aa%>%dist_permanova(seed = 1, variables = "fjord", n_processes = 1, n_perms = 999)
PERM2


# try corrplot

meta_samples_nf_red <- meta_samples_nf[, c(2,3,7,8,11,13,14)]

meta_samples_nf_red_CORR <- cor(meta_samples_nf_red, method = "pearson")

png("corr_test_18S.png", width = 600, height = 600)
corrplot(meta_samples_nf_red_CORR, method = 'color', addCoef.col = 'black'
         #title = "Correlation Plot - Meta data used for RDA analysis - NextSeq 18S",
         )
dev.off()

png("corr_test_mixed_18S.png", width = 600, height = 600)
corrplot.mixed(meta_samples_nf_red_CORR, order = 'AOE')
dev.off()

meta_samples_nbr_only <- meta_samples_nf[, c(2:14)] 
meta_samples_nbr_only__CORR <- cor(meta_samples_nbr_only, method = "pearson")

png("corr_test_all_vari_18S.png", width = 600, height = 600)
corrplot(meta_samples_nbr_only__CORR, method = 'color', addCoef.col = 'black')
dev.off()



## ---

NextSeq_fixed@sam_data$Temperature <- NextSeq_fixed@sam_data$temp
NextSeq_fixed@sam_data$Salinity <- NextSeq_fixed@sam_data$sal
NextSeq_fixed@sam_data$Density <- NextSeq_fixed@sam_data$dens
NextSeq_fixed@sam_data$Oxygen <- NextSeq_fixed@sam_data$oxyg
#NextSeq_fixed@sam_data$Attenuation <- NextSeq_fixed@sam_data$atten
NextSeq_fixed@sam_data$Chlorophyll_a <- NextSeq_fixed@sam_data$chla
NextSeq_fixed@sam_data$Latitude <- NextSeq_fixed@sam_data$lat
NextSeq_fixed@sam_data$Longitude <- NextSeq_fixed@sam_data$long

NextSeq_fixed %>% 
  tax_transform("clr", rank = "unique") %>%
  ord_calc(constraints = c("Temperature", "Salinity", "Density", "Oxygen", "Chlorophyll_a", "Latitude", "Longitude"),
           scale_cc = F, method = "RDA") %>%
  ord_plot(auto_caption = NA, color = "fjord", size = 4.5, plot_taxa = 1:10, constraint_lab_style = constraint_lab_style(size = 4), tax_lab_style = tax_lab_style(size = 4)) + theme_bw() + 
  scale_colour_manual(values = c("van Mijenfjorden" = "#ff00ad", "Kongsfjorden" = "#ab55c7", "Wijdefjorden" = "#5e91e4", "Rijpfjorden"= "#0edaff"), 
                      breaks = c('van Mijenfjorden', 'Kongsfjorden', 'Wijdefjorden', "Rijpfjorden")) +
  labs(title = "RDA Analysis - NextSeq 18S - Top 10 ASVs", color = "Fjord") + theme(legend.text=element_text(size=12), 
                                                                                    legend.title=element_text(size=16),
                                                                                    plot.title=element_text(size=20),
                                                                                    axis.title = element_text(size=14)) 

ggsave("RDA_NextSeq_18S_try_asvs_2.png", width = 13, height = 8)
