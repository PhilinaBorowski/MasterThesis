## miseq test - look at filters and net/pump seperately
devtools::session_info()
# to see which packages were used?

library(ggside)
library(viridis)           # Load

library(dplyr)
library(ggplot2)
library(phyloseq)
library(tibble)
#library(microbiomeMarker)
library(ape)
library(vegan)
library(microViz)
library(tidyr)
library(zCompositions)
library(compositions)

## set wd
setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/Phyloseq_analysis/")
wd <- ("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/Phyloseq_analysis/")



asvsf <- read.csv("MiSeq_test_abundance_only_filters.csv")           ## need the ones wo metazoa + singletons !!!
asvsnp <- read.csv("MiSeq_test_abundance_only_NP.csv")
taxa <- read.csv("MiSeq_18S_ASV_taxonomy_wo_singl_metaz_for_phyloseq_analysis.csv")                 ###ASV0001 - ASV2502
meta_samples <- read.csv("meta_phyloseq_18S_miseq.csv")  

asvsf <- asvsf %>% tibble::column_to_rownames("otu")
asvsnp <- asvsnp %>% tibble::column_to_rownames("otu")
taxa <- taxa %>% tibble::column_to_rownames("otu")
meta_samples <- meta_samples %>% tibble::column_to_rownames("sample")


## prep for CLR + CLR 
#d.czm <- cmultRepl(t(d.1),  label=0, method="CZM", frac= 0.99, adjust=FALSE, z.warning = 0.99)
#d.clr <- t(apply(d.czm, 1, function(x){log(x) - mean(log(x))})) #centered log ratio transformation

d.czmf <- cmultRepl(t(asvsf),  label=0, method="CZM", frac= 0.999, adjust=FALSE, z.warning = 0.999)
d.clrf <- (apply(d.czmf, 1, function(x){log(x) - mean(log(x))})) #centered log ratio transformation

#move all into positive values
asvs_normf <- apply(d.clrf, 2, function(x) x - min(x[x < 0]))

#make it relative
asvs_normf_P <- apply(asvs_normf, 2, function(x) {x/sum(x)}) 

asvs_normnp_P <- apply(asvsnp, 2, function(x) {x/sum(x)})



## creating phyloseq objects - for filters

asvs_normf <- as.matrix(asvs_normf)

taxa <- as.matrix(taxa)
#taxa <- as.data.frame(taxa)
str(asvs_normf)

OTUf = otu_table(asvs_normf, taxa_are_rows = T)
TAX = tax_table(taxa)
samples = sample_data(meta_samples)

## do the names match? 
sample_names(OTUf)           # all - are now a . 
sample_names(samples)        

MiSeq_18Sf <- phyloseq(OTUf, TAX, samples)
MiSeq_18Sf
# 71 samples
# 2468 taxa

MiSeq_18Sf@otu_table

## sort from van Mijen to Rijfjorden, from outher to inner fjord

MiSeq_18Sf <- ps_reorder(MiSeq_18Sf, sample_order = c("M.St1.F02.1", "M.St1.F02.2", "M.St1.F02.3",	"M.St1.F3.1",	"M.St1.F3.2",	"M.St1.F3.3",
                                            "M.St2.F02.1",	"M.St2.F02.2",	"M.St2.F3.1",	"M.St2.F3.2",	"M.St2.F3.3",	
                                            "M.St4.F02.1",	"M.St4.F02.2",	"M.St4.F02.3",	"M.St4.F3.1",	"M.St4.F3.2",	"M.St4.F3.3",
                                            "M.St3.F02.1",	"M.St3.F02.2",	"M.St3.F02.3",	"M.St3.F3.1",	"M.St3.F3.2",	"M.St3.F3.3",
                                            "M.St12.F02.1",	"M.St12.F02.2",	"M.St12.F02.3",	"M.St12.F3.1",	"M.St12.F3.2",	"M.St12.F3.3",	
                                            "M.St13.F02.1_N",	"M.St13.F02.2",	"M.St13.F02.3",	"M.St13.F3.1",	"M.St13.F3.2",	"M.St13.F3.3",
                                            "M.St11.F02.1",	"M.St11.F02.2",	"M.St11.F02.3",	"M.St11.F3.1",	"M.St11.F3.2",	"M.St11.F3.3",	
                                            "M.St5.F02.1",	"M.St5.F02.2",	"M.St5.F02.3",	"M.St5.F3.1",	"M.St5.F3.2",	"M.St5.F3.3",	
                                            "M.St7.F02.1",	"M.St7.F02.2",	"M.St7.F02.3",	"M.St7.F3.1",	"M.St7.F3.2",	"M.St7.F3.3_N",
                                            "M.St6.F02.1",	"M.St6.F02.2",	"M.St6.F02.3",	"M.St6.F3.1",	"M.St6.F3.2",	"M.St6.F3.3",	
                                            "M.St10.F02.1",	"M.St10.F02.2",	"M.St10.F02.3",	"M.St10.F3.1",	"M.St10.F3.2",	"M.St10.F3.3",	
                                            "M.St8.F02.1",	"M.St8.F02.2",	"M.St8.F02.3",	"M.St8.F3.1",	"M.St8.F3.2",	"M.St8.F3.3"))

MiSeq_18Sf@sam_data

MiSeq_18Sf@otu_table # right order 
View(MiSeq_18Sf@otu_table)


# percentage 

asvs_normf_P <- as.matrix(asvs_normf_P)
str(asvs_normf_P)

OTUfP = otu_table(asvs_normf_P, taxa_are_rows = T)

MiSeq_18SfP <- phyloseq(OTUfP, TAX, samples)
MiSeq_18SfP

MiSeq_18SfP <- ps_reorder(MiSeq_18SfP, sample_order = c("M.St1.F02.1", "M.St1.F02.2", "M.St1.F02.3",	"M.St1.F3.1",	"M.St1.F3.2",	"M.St1.F3.3",
                                                      "M.St2.F02.1",	"M.St2.F02.2",	"M.St2.F3.1",	"M.St2.F3.2",	"M.St2.F3.3",	
                                                      "M.St4.F02.1",	"M.St4.F02.2",	"M.St4.F02.3",	"M.St4.F3.1",	"M.St4.F3.2",	"M.St4.F3.3",
                                                      "M.St3.F02.1",	"M.St3.F02.2",	"M.St3.F02.3",	"M.St3.F3.1",	"M.St3.F3.2",	"M.St3.F3.3",
                                                      "M.St12.F02.1",	"M.St12.F02.2",	"M.St12.F02.3",	"M.St12.F3.1",	"M.St12.F3.2",	"M.St12.F3.3",	
                                                      "M.St13.F02.1_N",	"M.St13.F02.2",	"M.St13.F02.3",	"M.St13.F3.1",	"M.St13.F3.2",	"M.St13.F3.3",
                                                      "M.St11.F02.1",	"M.St11.F02.2",	"M.St11.F02.3",	"M.St11.F3.1",	"M.St11.F3.2",	"M.St11.F3.3",	
                                                      "M.St5.F02.1",	"M.St5.F02.2",	"M.St5.F02.3",	"M.St5.F3.1",	"M.St5.F3.2",	"M.St5.F3.3",	
                                                      "M.St7.F02.1",	"M.St7.F02.2",	"M.St7.F02.3",	"M.St7.F3.1",	"M.St7.F3.2",	"M.St7.F3.3_N",
                                                      "M.St6.F02.1",	"M.St6.F02.2",	"M.St6.F02.3",	"M.St6.F3.1",	"M.St6.F3.2",	"M.St6.F3.3",	
                                                      "M.St10.F02.1",	"M.St10.F02.2",	"M.St10.F02.3",	"M.St10.F3.1",	"M.St10.F3.2",	"M.St10.F3.3",	
                                                      "M.St8.F02.1",	"M.St8.F02.2",	"M.St8.F02.3",	"M.St8.F3.1",	"M.St8.F3.2",	"M.St8.F3.3"))

sample_data(MiSeq_18SfP)$sample_order <- factor(sample_names(MiSeq_18SfP))

MiSeq_18SfP@otu_table # right order
MiSeq_18SfP@sam_data

## -- only NP

asvsnp <- as.matrix(asvsnp)

taxa <- as.matrix(taxa)
#taxa <- as.data.frame(taxa)
str(asvsnp)

OTUnp = otu_table(asvsnp, taxa_are_rows = T)
TAX = tax_table(taxa)
samples = sample_data(meta_samples)

## do the names match? 
sample_names(OTUnp)           # all - are now a . 
sample_names(samples)        

MiSeq_18Snp <- phyloseq(OTUnp, TAX, samples)
MiSeq_18Snp
# 23 samples

MiSeq_18Snp <- ps_reorder(MiSeq_18Snp, sample_order = c("M.St1.N",	"M.St1.P1",	"M.St1.P3",
                                                        "M.St2.P1",	"M.St2.P2",	"M.St2.P3",
                                                        "M.St4.P1",	"M.St4.P2",	"M.St4.P3",
                                                        "M.St3.P1",	"M.St3.P2",	"M.St3.P3",
                                                        "M.St5.P1",	"M.St5.P2",	"M.St5.P3",
                                                        "M.St7.P1",	"M.St7.P2",	"M.St7.P3",
                                                        "M.St6.P1", "M.St6.P2",	"M.St6.P3",
                                                         "M.St8.P1",	"M.St8.P2"))

sample_data(MiSeq_18Snp)$sample_order <- factor(sample_names(MiSeq_18Snp))

MiSeq_18Snp@otu_table # right order
MiSeq_18Snp@sam_data

# percentage np

asvs_normnp_P <- as.matrix(asvs_normnp_P)
str(asvs_normnp_P)

OTUnpP = otu_table(asvs_normnp_P, taxa_are_rows = T)

MiSeq_18SnpP <- phyloseq(OTUnpP, TAX, samples)
MiSeq_18SnpP

MiSeq_18SnpP <- ps_reorder(MiSeq_18SnpP, sample_order = c("M.St1.N",	"M.St1.P1",	"M.St1.P3",
                                                        "M.St2.P1",	"M.St2.P2",	"M.St2.P3",
                                                        "M.St4.P1",	"M.St4.P2",	"M.St4.P3",
                                                        "M.St3.P1",	"M.St3.P2",	"M.St3.P3",
                                                        "M.St5.P1",	"M.St5.P2",	"M.St5.P3",
                                                        "M.St7.P1",	"M.St7.P2",	"M.St7.P3",
                                                        "M.St6.P1", "M.St6.P2",	"M.St6.P3",
                                                        "M.St8.P1",	"M.St8.P2"))


sample_data(MiSeq_18SnpP)$sample_order <- factor(sample_names(MiSeq_18SnpP))

MiSeq_18SnpP@otu_table # right order
MiSeq_18SnpP@sam_data





### ----- 
## plots

# filters

## try to merge and only look at phylum
mf.phy = tax_glom(MiSeq_18Sf, taxrank="Phylum")

MiSeq_18Sf@sam_data
mf.phy@sam_data

sample_data(mf.phy)$NewID <- factor(sample_names(mf.phy))
sample_data(mf.phy)$NewID <- factor(sample_data(mf.phy)$NewID, levels = (sample_data(mf.phy)$sample_order))

mf.phy@sam_data


png("MiSeq_filter_phylum_try.png", width = 1500, height = 800)
plot_bar(mf.phy, x = "NewID", fill = "Phylum") + labs(title = "MiSeq 18S - Phylum per filter sample - CLR transformed") + xlab("Sample") +
  theme_light() + scale_fill_viridis(discrete = T, option = "C") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                                         axis.text=element_text(size=14), #change font size of axis text
                                                                         axis.title=element_text(size=19), #change font size of axis titles
                                                                         plot.title=element_text(size=26), #change font size of plot title
                                                                         legend.text=element_text(size=16), #change font size of legend text
                                                                         legend.title=element_text(size=20)) #change font size of legend title   )
dev.off()



## --- 


mfP.phy = tax_glom(MiSeq_18SfP, taxrank="Phylum", NArm=FALSE)

mfP.phy@sam_data

sample_data(mfP.phy)$NewID <- factor(sample_names(mfP.phy))
sample_data(mfP.phy)$NewID <- factor(sample_data(mfP.phy)$NewID, levels = (sample_data(mfP.phy)$sample_order))

mfP.phy@sam_data

png("MiSeq_filter_phylum_try_P.png", width = 1500, height = 800)
plot_bar(mfP.phy, x = "NewID", fill = "Phylum") + labs(title = "MiSeq 18S - Phylum per filter sample - Relative - CLR transformed") + xlab("Sample") + 
  theme_light() + scale_fill_viridis(discrete = T, option = "C") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                                          axis.text=element_text(size=14), #change font size of axis text
                                                                          axis.title=element_text(size=19), #change font size of axis titles
                                                                          plot.title=element_text(size=26), #change font size of plot title
                                                                          legend.text=element_text(size=16), #change font size of legend text
                                                                          legend.title=element_text(size=20)) #change font size of legend title   )
dev.off()

# net/pump

mnp.phy = tax_glom(MiSeq_18Snp, taxrank="Phylum", NArm=FALSE)

mnp.phy@sam_data

sample_data(mnp.phy)$NewID <- factor(sample_names(mnp.phy))
sample_data(mnp.phy)$NewID <- factor(sample_data(mnp.phy)$NewID, levels = (sample_data(mnp.phy)$sample_order))

mnp.phy@sam_data

png("MiSeq_netpump_phylum_try.png", width = 1500, height = 800)
plot_bar(mnp.phy, x = "NewID", fill = "Phylum") + labs(title = "MiSeq 18S - Phylum per net/pump sample") + xlab("Sample") +
  theme_light() + scale_fill_viridis(discrete = T, option = "C") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                                         axis.text=element_text(size=14), #change font size of axis text
                                                                         axis.title=element_text(size=19), #change font size of axis titles
                                                                         plot.title=element_text(size=26), #change font size of plot title
                                                                         legend.text=element_text(size=16), #change font size of legend text
                                                                         legend.title=element_text(size=20)) #change font size of legend title   )
dev.off()

# net/pump %

mnpP.phy = tax_glom(MiSeq_18SnpP, taxrank="Phylum", NArm=FALSE)

mnpP.phy@sam_data

sample_data(mnpP.phy)$NewID <- factor(sample_names(mnpP.phy))
sample_data(mnpP.phy)$NewID <- factor(sample_data(mnpP.phy)$NewID, levels = (sample_data(mnpP.phy)$sample_order))

mnpP.phy@sam_data

png("MiSeq_netpump_phylum_try_P.png", width = 1500, height = 800)
plot_bar(mnpP.phy, x = "NewID", fill = "Phylum") + labs(title = "MiSeq 18S - Phylum per net/pump sample - Relative") + xlab("Sample") +
  theme_light() + scale_fill_viridis(discrete = T, option = "C") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                                          axis.text=element_text(size=14), #change font size of axis text
                                                                          axis.title=element_text(size=19), #change font size of axis titles
                                                                          plot.title=element_text(size=26), #change font size of plot title
                                                                          legend.text=element_text(size=16), #change font size of legend text
                                                                          legend.title=element_text(size=20)) #change font size of legend title   )
## + geom_bar(aes(fill=Phylum), stat="identity", position="stack") 
## for a smooth overlay of colours
dev.off()








## 
## try to merge triplicates together? and look at them together? 
## try mean of sample values

Mi_18f_merged <- merge_samples(mf.phy, "sample_merg", fun = mean)
#In asMethod(object) : NAs introduced by coercion

sample_names(Mi_18f_merged)
sample_variables(Mi_18f_merged)

Mi_18f_merged # 24 samples

png("Mi_filter_merged.png", width = 1000)
plot_bar(Mi_18f_merged, fill="Phylum", title = "MiSeq 18S - Merged filter samples (mean of replicates)") +
  theme_light() + scale_fill_viridis(discrete = T, option = "C") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                                          axis.text=element_text(size=14), #change font size of axis text
                                                                          axis.title=element_text(size=19), #change font size of axis titles
                                                                          plot.title=element_text(size=26), #change font size of plot title
                                                                          legend.text=element_text(size=14), #change font size of legend text
                                                                          legend.title=element_text(size=18)) #change font size of legend title   )
dev.off()



## ???????????????????????????????
#Mi_18f_merged_sum <- merge_samples(mf.phy, "sample_merg", fun = sum)
#In asMethod(object) : NAs introduced by coercion

#png("Mi_filter_merged_sum.png", width = 1000)
#plot_bar(Mi_18f_merged_sum, fill="Phylum", title = "MiSeq 18S - Merged samples (sum of replicates)") +
#  theme_light() + scale_fill_viridis(discrete = T, option = "C") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
#                                                                          axis.text=element_text(size=14), #change font size of axis text
#                                                                          axis.title=element_text(size=19), #change font size of axis titles
#                                                                          plot.title=element_text(size=26), #change font size of plot title
#                                                                          legend.text=element_text(size=16), #change font size of legend text
#                                                                          legend.title=element_text(size=20)) #change font size of legend title   )
#dev.off()





Mi_18fP_merged <- merge_samples(mfP.phy, "sample_merg", fun = mean)
#In asMethod(object) : NAs introduced by coercion

png("Mi_filter_mergedP.png", width = 1000)
plot_bar(Mi_18fP_merged, fill="Phylum", title = "MiSeq 18S - Merged filter samples (mean of replicates) - Relatives") +
  theme_light() + scale_fill_viridis(discrete = T, option = "C") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                                          axis.text=element_text(size=14), #change font size of axis text
                                                                          axis.title=element_text(size=19), #change font size of axis titles
                                                                          plot.title=element_text(size=26), #change font size of plot title
                                                                          legend.text=element_text(size=14), #change font size of legend text
                                                                          legend.title=element_text(size=18)) #change font size of legend title   )
dev.off()


# net/pumo

Mi_18np_merged <- merge_samples(mnp.phy, "sample_merg", fun = mean)
#In asMethod(object) : NAs introduced by coercion

png("Mi_netpump_merged.png", width = 1000)
plot_bar(Mi_18np_merged, fill="Phylum", title = "MiSeq 18S - Merged net/pump samples (mean of replicates)") +
  theme_light() + scale_fill_viridis(discrete = T, option = "C") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                                          axis.text=element_text(size=14), #change font size of axis text
                                                                          axis.title=element_text(size=19), #change font size of axis titles
                                                                          plot.title=element_text(size=26), #change font size of plot title
                                                                          legend.text=element_text(size=14), #change font size of legend text
                                                                          legend.title=element_text(size=18)) #change font size of legend title   )
dev.off()



Mi_18npP_merged <- merge_samples(mnpP.phy, "sample_merg", fun = mean)
#In asMethod(object) : NAs introduced by coercion

png("Mi_netpump_mergedP.png", width = 1000)
plot_bar(Mi_18npP_merged, fill="Phylum", title = "MiSeq 18S - Merged net/pump samples (mean of replicates) - Relatives") +
  theme_light() + scale_fill_viridis(discrete = T, option = "C") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                                          axis.text=element_text(size=14), #change font size of axis text
                                                                          axis.title=element_text(size=19), #change font size of axis titles
                                                                          plot.title=element_text(size=26), #change font size of plot title
                                                                          legend.text=element_text(size=14), #change font size of legend text
                                                                          legend.title=element_text(size=18)) #change font size of legend title   )
dev.off()





### ordination
# PCA
# beta diversity !!



# filter samples !!


MiSeq_18Sf@sam_data
# create a new phyloseq wo/ clr transf before

asvsf <- as.matrix(asvsf)

taxa <- as.matrix(taxa)
#taxa <- as.data.frame(taxa)
str(asvsf)

OTUn = otu_table(asvsf, taxa_are_rows = T)
TAX = tax_table(taxa)
samples = sample_data(meta_samples)

## do the names match? 
sample_names(OTUn)           # all - are now a . 
sample_names(samples)        

MiSeq_18Sf_default <- phyloseq(OTUn, TAX, samples)
MiSeq_18Sf_default
# 71 samples
# 2502 taxa

MiSeq_18Sf_default@tax_table
# some NAs

tax_fix_interactive(MiSeq_18Sf_default)

MiSeq_fixed <- MiSeq_18Sf_default %>%
  tax_fix(
    min_length = 4,
    unknowns = c("NA"),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified")


str(MiSeq_fixed@sam_data)
MiSeq_fixed@sam_data

MiSeq_fixed@sam_data$stat <- as.factor(MiSeq_fixed@sam_data$stat)
MiSeq_18Sf_default@sam_data$stat <- as.factor(MiSeq_18Sf_default@sam_data$stat)



#MiSeq_fixed %>% tax_transform("clr", rank ="Class") %>%
#                ord_calc() %>% 
#                ord_plot(color = "stat", shape = "fjord") + theme_bw() + labs(title = "MiSeq 18S - PCA Class")
#ggsave("pca_try_class.png")



## asvs - the nicest ones cause lot of info
MiSeq_fixed %>% tax_transform("clr", rank ="unique") %>%
  ord_calc() %>% 
  ord_plot(color = "fjord") + theme_bw() +
  scale_colour_manual(values = c("van Mijenfjorden" = "#ff00ad", "Kongsfjorden" = "#ab55c7", "Wijdefjorden" = "#5e91e4", "Rijpfjorden"= "#0edaff"), 
                      breaks = c('van Mijenfjorden', 'Kongsfjorden', 'Wijdefjorden', "Rijpfjorden")) +
  labs(title = "PCA Analysis - MiSeq 18S ASVs", color = "Fjord") + theme(legend.text=element_text(size=10), 
                                                                                      legend.title=element_text(size=12),
                                                                                      plot.title=element_text(size=16)) 

ggsave("pca_try_asvs1.png", width = 8, height = 7)


# axis 1 + 3
MiSeq_fixed %>% tax_transform("clr", rank ="unique") %>%
  ord_calc() %>% 
  ord_plot(axes = c(1,3), color = "fjord") + theme_bw() +
  scale_colour_manual(values = c("van Mijenfjorden" = "#ff00ad", "Kongsfjorden" = "#ab55c7", "Wijdefjorden" = "#5e91e4", "Rijpfjorden"= "#0edaff"), 
                      breaks = c('van Mijenfjorden', 'Kongsfjorden', 'Wijdefjorden', "Rijpfjorden")) +
  labs(title = "PCA Analysis - MiSeq 18S ASVs - Axis 1 + 3", color = "Fjord") + theme(legend.text=element_text(size=10), 
                                                                    legend.title=element_text(size=12),
                                                                    plot.title=element_text(size=16)) 



ggsave("pca_try_asvs_2.png", width = 8, height = 7)



MiSeq_18Sf_default %>%
  tax_transform("identity", rank = "unique") %>% # don't transform!
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_get() %>%
  phyloseq::plot_scree() + theme(axis.text.x = element_text(size = 6)) + labs(title = "MiSeq 18S - PCA eigenvalues")

ggsave("pca_axes_asvs.png", width = 10, height= 7)






## ordination pump + net samples !!!





## try RDA
# station
# t
# s 
#MiSeq_18Sf_default@sam_data$fjord <- as.factor(MiSeq_18Sf_default@sam_data$fjord)
MiSeq_18Sf_default@sam_data$stat <- as.factor(MiSeq_18Sf_default@sam_data$stat)


#MiSeq_18Sf_default %>% 
#  tax_transform("clr", rank = "unique") %>%
#  ord_calc(constraints = c("temp", "sal"),
#           scale_cc = F, method = "RDA") %>%
#  ord_plot(color = "stat", shape = "fjord") + theme_bw()

#ggsave("rda_try.png")


MiSeq_18Sf_default %>% 
  tax_transform("clr", rank = "unique") %>%
  ord_calc(constraints = c("temp", "sal", "dens", "oxyg", "atten", "chla", "nobs", "lat", "long"),
           scale_cc = F, method = "RDA") %>%
  ord_plot(color = "fjord", plot_taxa = 1:10) + theme_bw() + 
  scale_colour_manual(values = c("van Mijenfjorden" = "#ff00ad", "Kongsfjorden" = "#ab55c7", "Wijdefjorden" = "#5e91e4", "Rijpfjorden"= "#0edaff"), 
                      breaks = c('van Mijenfjorden', 'Kongsfjorden', 'Wijdefjorden', "Rijpfjorden")) +
  labs(title = "RDA Analysis - MiSeq 18S - Top 10 ASVs", color = "Fjord") + theme(legend.text=element_text(size=10), 
                                                   legend.title=element_text(size=12),
                                                   plot.title=element_text(size=16)) 

ggsave("MiSeq_rda.png", width= 10, height = 8)


MiSeq_18Sf_default %>% 
  tax_transform("clr", rank = "unique") %>%
  ord_calc(constraints = c("temp", "sal", "dens", "oxyg", "atten", "chla", "nobs", "lat", "long"),
           scale_cc = F, method = "RDA") %>%
  ord_plot(color = "fjord") + theme_bw() + 
  scale_colour_manual(values = c("van Mijenfjorden" = "#ff00ad", "Kongsfjorden" = "#ab55c7", "Wijdefjorden" = "#5e91e4", "Rijpfjorden"= "#0edaff"), 
                      breaks = c('van Mijenfjorden', 'Kongsfjorden', 'Wijdefjorden', "Rijpfjorden")) +
  labs(title = "RDA Analysis - MiSeq 18S", color = "Fjord") + theme(legend.text=element_text(size=10), 
                                                                    legend.title=element_text(size=12),
                                                                    plot.title=element_text(size=16)) 

ggsave("MiSeq_rda_same_but_wo_taxa.png", width= 10, height = 8)




MiSeq_18Sf_default@sam_data$fjord <- as.integer(MiSeq_18Sf_default@sam_data$fjord)











