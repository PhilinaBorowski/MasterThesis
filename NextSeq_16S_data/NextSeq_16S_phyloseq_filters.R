## NextSeq - 16S - only filter samples
devtools::session_info()
# to see which packages were used?

library(ggside)
library(viridis)           # 
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
library(paletteer)


## set wd
setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/phyloseq_analysis/")
wd <- ("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/phyloseq_analysis/")


asvs <- read.csv("NextSeq_16S_asvs.csv")           ## need the ones wo mitochondira/chloroplasts + singletons !!!
taxa <- read.csv("NextSeq_16S_taxa.csv")                 ###ASV0001 - ASV2044
meta_samples <- read.csv("NextSeq_16S_metadata.csv")  

asvs <- asvs %>% tibble::column_to_rownames("otu")
taxa <- taxa %>% tibble::column_to_rownames("otu")
meta_samples <- meta_samples %>% tibble::column_to_rownames("sample")


## prep for CLR + CLR 
#d.czm <- cmultRepl(t(d.1),  label=0, method="CZM", frac= 0.99, adjust=FALSE, z.warning = 0.99)
#d.clr <- t(apply(d.czm, 1, function(x){log(x) - mean(log(x))})) #centered log ratio transformation

d.czmf <- cmultRepl(t(asvs),  label=0, method="CZM", frac= 0.999, adjust=FALSE, z.warning = 0.999)
d.clrf <- (apply(d.czmf, 1, function(x){log(x) - mean(log(x))})) #centered log ratio transformation

#move all into positive values
asvs_normf <- apply(d.clrf, 2, function(x) x - min(x[x < 0]))

#make it relative
asvs_normf_P <- apply(asvs_normf, 2, function(x) {x/sum(x)}) 



## creating phyloseq object

asvs_normf <- as.matrix(asvs_normf)

taxa <- as.matrix(taxa)
#taxa <- as.data.frame(taxa)
str(asvs_normf)

OTUf = otu_table(asvs_normf, taxa_are_rows = T)
TAX = tax_table(taxa)
samples = sample_data(meta_samples)

## do the names match? 
sample_names(OTUf)           
sample_names(samples)        

NextSeq_16Sf <- phyloseq(OTUf, TAX, samples)
NextSeq_16Sf
# 36 samples
# 2044 taxa

NextSeq_16Sf@otu_table

## sort from van Mijen to Rijfjorden, from outher to inner fjord

NextSeq_16Sf <- ps_reorder(NextSeq_16Sf, sample_order = c("Pro_St1.F02.1", "Pro_St1.F02.2", "Pro_St1.F02.3",   #3
                                                      "Pro_St2.F02.1",	"Pro_St2.F02.2", "Pro_St2.F02.3",      #3
                                                      "Pro_St4.F02.1",	"Pro_St4.F02.2",	"Pro_St4.F02.3",     #3
                                                      "Pro_St3.F02.1",	"Pro_St3.F02.2",	"Pro_St3.F02.3",	   #3
                                                      "Pro_St12.F02.1",	"Pro_St12.F02.2",	"Pro_St12.F02.3",	   #3
                                                      "Pro_St13.F02.1",	"Pro_St13.F02.2",	"Pro_St13.F02.3",	   #3
                                                      "Pro_St11.F02.1",	"Pro_St11.F02.2",	"Pro_St11.F02.3",    #3
                                                      "Pro_St5.F02.1",	"Pro_St5.F02.2",	"Pro_St5.F02.3",	   #3
                                                      "Pro_St7.F02.1",	"Pro_St7.F02.2",	"Pro_St7.F02.3",     #3
                                                      "Pro_St6.F02.1",	"Pro_St6.F02.2",	"Pro_St6.F02.3",	   #3
                                                      "Pro_St10.F02.1",	"Pro_St10.F02.2",	"Pro_St10.F02.3",	   #3
                                                      "Pro_St8.F02.1",	"Pro_St8.F02.2",	"Pro_St8.F02.3"))    #3

sample_data(NextSeq_16Sf)$sample_order <- factor(sample_names(NextSeq_16Sf))

NextSeq_16Sf@sam_data

NextSeq_16Sf@otu_table # right order 
View(NextSeq_16Sf@otu_table)


# percentage 

asvs_normf_P <- as.matrix(asvs_normf_P)
str(asvs_normf_P)

OTUfP = otu_table(asvs_normf_P, taxa_are_rows = T)

NextSeq_16SfP <- phyloseq(OTUfP, TAX, samples)
NextSeq_16SfP

NextSeq_16SfP <- ps_reorder(NextSeq_16SfP, sample_order = c("Pro_St1.F02.1", "Pro_St1.F02.2", "Pro_St1.F02.3",       #3
                                                            "Pro_St2.F02.1",	"Pro_St2.F02.2", "Pro_St2.F02.3",      #3
                                                            "Pro_St4.F02.1",	"Pro_St4.F02.2",	"Pro_St4.F02.3",     #3
                                                            "Pro_St3.F02.1",	"Pro_St3.F02.2",	"Pro_St3.F02.3",	   #3
                                                            "Pro_St12.F02.1",	"Pro_St12.F02.2",	"Pro_St12.F02.3",	   #3
                                                            "Pro_St13.F02.1",	"Pro_St13.F02.2",	"Pro_St13.F02.3",	   #3
                                                            "Pro_St11.F02.1",	"Pro_St11.F02.2",	"Pro_St11.F02.3",    #3
                                                            "Pro_St5.F02.1",	"Pro_St5.F02.2",	"Pro_St5.F02.3",	   #3
                                                            "Pro_St7.F02.1",	"Pro_St7.F02.2",	"Pro_St7.F02.3",     #3
                                                            "Pro_St6.F02.1",	"Pro_St6.F02.2",	"Pro_St6.F02.3",	   #3
                                                            "Pro_St10.F02.1",	"Pro_St10.F02.2",	"Pro_St10.F02.3",	   #3
                                                            "Pro_St8.F02.1",	"Pro_St8.F02.2",	"Pro_St8.F02.3"))    #3

sample_data(NextSeq_16SfP)$sample_order <- factor(sample_names(NextSeq_16SfP))

NextSeq_16SfP@otu_table # right order
NextSeq_16SfP@sam_data







### ----- 
## plots

# filters

## try to merge and only look at phylum
n16.phy = tax_glom(NextSeq_16Sf, taxrank="Phylum", NArm = F)

NextSeq_16Sf@sam_data
n16.phy@sam_data

sample_data(n16.phy)$NewID <- factor(sample_names(n16.phy))
sample_data(n16.phy)$NewID <- factor(sample_data(n16.phy)$NewID, levels = (sample_data(n16.phy)$sample_order))
# dont have a sample order here? 

n16.phy@sam_data


png("NextSeq_16S_phylum_asvs.png", width = 1500, height = 800)
plot_bar(n16.phy, x = "NewID", fill = "Phylum") + labs(title = "NextSeq 16S - Phylum per filter sample - CLR transformed") + xlab("Sample") +
  theme_light() + scale_fill_viridis(discrete = T, option = "D") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                                         axis.text=element_text(size=14), #change font size of axis text
                                                                         axis.title=element_text(size=19), #change font size of axis titles
                                                                         plot.title=element_text(size=26), #change font size of plot title
                                                                         legend.text=element_text(size=16), #change font size of legend text
                                                                         legend.title=element_text(size=20)) #change font size of legend title   )
dev.off()



## try to merge and only look at genus
n16.gen.genus = tax_glom(NextSeq_16Sf, taxrank="Genus", NArm = F)

sample_data(n16.gen.genus)$NewID <- factor(sample_names(n16.gen.genus))
sample_data(n16.gen.genus)$NewID <- factor(sample_data(n16.gen.genus)$NewID, levels = (sample_data(n16.gen.genus)$sample_order))

TOPS = names(sort(taxa_sums(n16.gen.genus), TRUE)[1:38]) # change 1:X till I have top20
TOPS
TOPS1 = prune_taxa(TOPS, n16.gen.genus)
TOPS1

png("NextSeq_16S_genus_asvs_top20.png", width = 1500, height = 800)
plot_bar(TOPS1, x = "NewID", fill = "Genus") + labs(title = "NextSeq 16S - Genus per filter sample - Top 20 - CLR transformed") + xlab("Sample") +
  theme_light() + scale_fill_viridis(discrete = T, option = "D") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                                         axis.text=element_text(size=14), #change font size of axis text
                                                                         axis.title=element_text(size=19), #change font size of axis titles
                                                                         plot.title=element_text(size=26), #change font size of plot title
                                                                         legend.text=element_text(size=16), #change font size of legend text
                                                                         legend.title=element_text(size=20)) #change font size of legend title   )
dev.off()

#try2
drop_na(TOPS1)

TOPS1_2 <- subset_taxa(TOPS1, !is.na(Genus) & !Genus %in% c("", "uncharacterized"))


png("NextSeq_16S_genus_asvs_top20_noNAs_rainbow.png", width = 1500, height = 800)
plot_bar(TOPS1_2, !is.na("Genus"), x = "NewID", y = "Abundance", fill = "Genus") + labs(title = "NextSeq 16S - Genus per filter sample - Top 20 - CLR transformed") + xlab("Sample") +
  theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                                         axis.text=element_text(size=14), #change font size of axis text
                                                                         axis.title=element_text(size=19), #change font size of axis titles
                                                                         plot.title=element_text(size=26), #change font size of plot title
                                                                         legend.text=element_text(size=16), #change font size of legend text
                                                                         legend.title=element_text(size=20)) #change font size of legend title   )
dev.off()


rhg_cols <- c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F", "#8CD17D", "#B6992D", "#F1CE63", "#499894", "#86BCB6",
                       "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", "#D37295", "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6")

png("NextSeq_16S_genus_asvs_top20_noNAs_other_colour.png", width = 1500, height = 800)
plot_bar(TOPS1_2, !is.na("Genus"), x = "NewID", y = "Abundance", fill = "Genus") + labs(title = "NextSeq 16S - Genus per filter sample - Top 20 - CLR transformed") + xlab("Sample") +
  theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                        axis.text=element_text(size=14), #change font size of axis text
                        axis.title=element_text(size=19), #change font size of axis titles
                        plot.title=element_text(size=26), #change font size of plot title
                        legend.text=element_text(size=16), #change font size of legend text
                        legend.title=element_text(size=20)) +  #change font size of legend title   ) 
scale_fill_manual(values = rhg_cols)
dev.off()


## --- 


n16P.phy = tax_glom(NextSeq_16SfP, taxrank="Phylum", NArm=FALSE)

n16P.phy@sam_data

sample_data(n16P.phy)$NewID <- factor(sample_names(n16P.phy))
sample_data(n16P.phy)$NewID <- factor(sample_data(n16P.phy)$NewID, levels = (sample_data(n16P.phy)$sample_order))

n16P.phy@sam_data

png("NextSeq_16S_phylum_asvs_P.png", width = 1500, height = 800)
plot_bar(n16P.phy, x = "NewID", fill = "Phylum") + labs(title = "NextSeq 16S - Phylum per filter sample - Relative - CLR transformed") + xlab("Sample") + 
  theme_light() + scale_fill_viridis(discrete = T, option = "D") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                                          axis.text=element_text(size=14), #change font size of axis text
                                                                          axis.title=element_text(size=19), #change font size of axis titles
                                                                          plot.title=element_text(size=26), #change font size of plot title
                                                                          legend.text=element_text(size=16), #change font size of legend text
                                                                          legend.title=element_text(size=20)) #change font size of legend title   )
dev.off()





## 
## try to merge triplicates together? and look at them together? 
## try mean of sample values

N_18f_merged <- merge_samples(n16.phy, "sample_merg", fun = mean)
#In asMethod(object) : NAs introduced by coercion

sample_names(N_18f_merged)
sample_variables(N_18f_merged)

N_18f_merged # 12 samples

png("N_filter_merged.png", width = 1000)
plot_bar(N_18f_merged, fill="Phylum", title = "NextSeq 16S - Merged filter samples (mean of replicates)") +
  theme_light() + scale_fill_viridis(discrete = T, option = "C") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                                          axis.text=element_text(size=14), #change font size of axis text
                                                                          axis.title=element_text(size=19), #change font size of axis titles
                                                                          plot.title=element_text(size=26), #change font size of plot title
                                                                          legend.text=element_text(size=14), #change font size of legend text
                                                                          legend.title=element_text(size=18)) #change font size of legend title   )
dev.off()






N_18fP_merged <- merge_samples(n16P.phy, "sample_merg", fun = mean)
#In asMethod(object) : NAs introduced by coercion

png("N_filter_mergedP.png", width = 1000)
plot_bar(N_18fP_merged, fill="Phylum", title = "NextSeq 16S - Merged filter samples (mean of replicates) - Relatives") +
  theme_light() + scale_fill_viridis(discrete = T, option = "C") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                                                          axis.text=element_text(size=14), #change font size of axis text
                                                                          axis.title=element_text(size=19), #change font size of axis titles
                                                                          plot.title=element_text(size=26), #change font size of plot title
                                                                          legend.text=element_text(size=14), #change font size of legend text
                                                                          legend.title=element_text(size=18)) #change font size of legend title   )
dev.off()




### ordination
# PCA
# beta diversity of the filter samples


NextSeq_16Sf@sam_data
# create a new phyloseq wo/ clr transf before

asvs <- as.matrix(asvs)

taxa <- as.matrix(taxa)
#taxa <- as.data.frame(taxa)
str(asvs)

OTUn = otu_table(asvs, taxa_are_rows = T)
TAX = tax_table(taxa)
samples = sample_data(meta_samples)

## do the names match? 
sample_names(OTUn)           # all - are now a . 
sample_names(samples)        

NextSeq_16Sf_default <- phyloseq(OTUn, TAX, samples)
NextSeq_16Sf_default
# 36 samples
# 2044 taxa

NextSeq_16Sf_default@tax_table
# some NAs

tax_fix_interactive(NextSeq_16Sf_default)

NextSeq_fixed <- NextSeq_16Sf_default %>%
  tax_fix(
    min_length = 4,
    unknowns = c("NA"),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified")


str(NextSeq_fixed@sam_data)
NextSeq_fixed@sam_data

NextSeq_fixed@sam_data$stat <- as.factor(NextSeq_fixed@sam_data$stat)
NextSeq_16Sf_default@sam_data$stat <- as.factor(NextSeq_16Sf_default@sam_data$stat)


## asvs - the nicest ones cause lot of info
NextSeq_fixed %>% tax_transform("clr", rank ="unique") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "fjord", size = 3) + theme_bw() +
  scale_colour_manual(values = c("van Mijenfjorden" = "#ff00ad", "Kongsfjorden" = "#ab55c7", "Wijdefjorden" = "#5e91e4", "Rijpfjorden"= "#0edaff"), 
                      breaks = c('van Mijenfjorden', 'Kongsfjorden', 'Wijdefjorden', "Rijpfjorden")) +
  labs(title = "PCA Analysis - NextSeq 16S ASVs", color = "Fjord") + theme(legend.text=element_text(size=10), 
                                                                         legend.title=element_text(size=12),
                                                                         plot.title=element_text(size=16)) +
  stat_ellipse(aes(colour = fjord))

ggsave("NextSeq_16S_pca_asvs1_try2.png", width = 6.5, height = 6)


# axis 1 + 3
NextSeq_fixed %>% tax_transform("clr", rank ="unique") %>%
  ord_calc() %>% 
  ord_plot(axes = c(1,3), color = "fjord", size = 3) + theme_bw() +
  scale_colour_manual(values = c("van Mijenfjorden" = "#ff00ad", "Kongsfjorden" = "#ab55c7", "Wijdefjorden" = "#5e91e4", "Rijpfjorden"= "#0edaff"), 
                      breaks = c('van Mijenfjorden', 'Kongsfjorden', 'Wijdefjorden', "Rijpfjorden")) +
  labs(title = "PCA Analysis - NextSeq 16S ASVs - Axis 1 + 3", color = "Fjord") + theme(legend.text=element_text(size=10), 
                                                                                      legend.title=element_text(size=12),
                                                                                      plot.title=element_text(size=16)) +
  stat_ellipse(aes(colour = fjord))



ggsave("NextSeq_16S_pca_try_asvs_try2.png", width = 6.5, height = 6)



NextSeq_16Sf_default %>%
  tax_transform("identity", rank = "unique") %>% # don't transform!
  dist_calc("aitchison") %>%
  ord_calc("PCoA") %>%
  ord_get() %>%
  phyloseq::plot_scree() + theme(axis.text.x = element_text(size = 6)) + labs(title = "NextSeq 16S - PCA eigenvalues")

ggsave("NextSeq_16S_pca_axes_asvs_try2.png", width = 10, height = 7)




# try pca with circles around
#NextSeq_fixed %>% tax_transform("clr", rank ="unique") %>%
#  ord_calc(method = "PCA") %>% 
#  ord_plot(color = "fjord", size = 3) + theme_bw() +
#  scale_colour_manual(values = c("van Mijenfjorden" = "#ff00ad", "Kongsfjorden" = "#ab55c7", "Wijdefjorden" = "#5e91e4", "Rijpfjorden"= "#0edaff"), 
#                      breaks = c('van Mijenfjorden', 'Kongsfjorden', 'Wijdefjorden', "Rijpfjorden")) +
#  labs(title = "PCA Analysis - NextSeq 16S ASVs", color = "Fjord") + theme(legend.text=element_text(size=10), 
#                                                                           legend.title=element_text(size=12),
#                                                                           plot.title=element_text(size=16)) +
#  stat_ellipse(aes(colour = fjord))

#ggsave("NextSeq_16S_pca_asvs1_try2.png", width = 6.5, height = 6)




## try RDA
# station
# t
# s 
NextSeq_16Sf_default@sam_data$stat <- as.factor(NextSeq_16Sf_default@sam_data$stat)


#NextSeq_16Sf_default %>% 
#  tax_transform("clr", rank = "unique") %>%
#  ord_calc(constraints = c("temp", "sal"),
#           scale_cc = F, method = "RDA") %>%
#  ord_plot(color = "stat", shape = "fjord") + theme_bw()

#ggsave("rda_try.png")

NextSeq_16Sf_default@sam_data$Temperature <- NextSeq_16Sf_default@sam_data$temp
NextSeq_16Sf_default@sam_data$Salinity <- NextSeq_16Sf_default@sam_data$sal
NextSeq_16Sf_default@sam_data$Density <- NextSeq_16Sf_default@sam_data$dens
NextSeq_16Sf_default@sam_data$Oxygen <- NextSeq_16Sf_default@sam_data$oxyg
NextSeq_16Sf_default@sam_data$Attenuation <- NextSeq_16Sf_default@sam_data$atten
NextSeq_16Sf_default@sam_data$Chlorophyll_a <- NextSeq_16Sf_default@sam_data$chla
NextSeq_16Sf_default@sam_data$Latitude <- NextSeq_16Sf_default@sam_data$lat
NextSeq_16Sf_default@sam_data$Longitude <- NextSeq_16Sf_default@sam_data$long

NextSeq_16Sf_default@sam_data

NextSeq_16Sf_default %>% 
  tax_transform("clr", rank = "unique") %>%
  ord_calc(constraints = c("temp", "sal", "dens", "oxyg", "atten", "chla", "lat", "long"),
           scale_cc = F, method = "RDA") %>%
  ord_plot(auto_caption = NA, color = "fjord", size = 4.5, plot_taxa = 1:10, constraint_lab_style = constraint_lab_style(size = 4), tax_lab_style = tax_lab_style(size = 4)) + theme_bw() + 
  scale_colour_manual(values = c("van Mijenfjorden" = "#ff00ad", "Kongsfjorden" = "#ab55c7", "Wijdefjorden" = "#5e91e4", "Rijpfjorden"= "#0edaff"), 
                      breaks = c('van Mijenfjorden', 'Kongsfjorden', 'Wijdefjorden', "Rijpfjorden")) +
  labs(title = "RDA Analysis - NextSeq 16S - Top 10 ASVs", color = "Fjord") + theme(legend.text=element_text(size=12), 
                                                                    legend.title=element_text(size=16),
                                                                    plot.title=element_text(size=20),
                                                                    axis.title = element_text(size=14)) 

ggsave("NextSeq_16S_rda_top10_asvs_try2.png", width= 13, height = 8)



ab <- NextSeq_16Sf_default %>% 
  #tax_transform("clr", rank = "unique") %>%
  dist_calc(dist = "aitchison")

set.seed(111)
PERMn16 <- ab%>%dist_permanova(seed = 1, variables = c("temp", "sal", "dens", "oxyg", "chla", "lat", "long"), n_processes = 1, n_perms = 999)
PERMn16

PERM216 <- ab%>%dist_permanova(seed = 1, variables = "fjord", n_processes = 1, n_perms = 999)
PERM216



# again but w longer meta names
NextSeq_16Sf_default %>% 
  tax_transform("clr", rank = "unique") %>%
  ord_calc(constraints = c("Temperature", "Salinity", "Density", "Oxygen", "Chlorophyll_a", "Latitude", "Longitude"),
           scale_cc = F, method = "RDA") %>%
  ord_plot(auto_caption = NA, color = "fjord", size = 4.5, plot_taxa = 1:10, constraint_lab_style = constraint_lab_style(size = 4), tax_lab_style = tax_lab_style(size = 4)) + theme_bw() + 
  scale_colour_manual(values = c("van Mijenfjorden" = "#ff00ad", "Kongsfjorden" = "#ab55c7", "Wijdefjorden" = "#5e91e4", "Rijpfjorden"= "#0edaff"), 
                      breaks = c('van Mijenfjorden', 'Kongsfjorden', 'Wijdefjorden', "Rijpfjorden")) +
  labs(title = "RDA Analysis - NextSeq 16S - Top 10 ASVs", color = "Fjord") + theme(legend.text=element_text(size=12), 
                                                                                    legend.title=element_text(size=16),
                                                                                    plot.title=element_text(size=20),
                                                                                    axis.title = element_text(size=14)) 

ggsave("NextSeq_16S_rda_top10_asvs_try2_2.png", width= 13, height = 8)



### not imporatant? but nice to only see fjords 


NextSeq_16Sf_default %>% 
  tax_transform("clr", rank = "unique") %>%
  ord_calc(constraints = c("temp", "sal", "dens", "oxyg", "atten", "chla", "nobs", "lat", "long"),
           scale_cc = F, method = "RDA") %>%
  ord_plot(color = "fjord", size = 3) + theme_bw() + 
  scale_colour_manual(values = c("van Mijenfjorden" = "#ff00ad", "Kongsfjorden" = "#ab55c7", "Wijdefjorden" = "#5e91e4", "Rijpfjorden"= "#0edaff"), 
                      breaks = c('van Mijenfjorden', 'Kongsfjorden', 'Wijdefjorden', "Rijpfjorden")) +
  labs(title = "RDA Analysis - NextSeq 16S", color = "Fjord") + theme(legend.text=element_text(size=10), 
                                                                    legend.title=element_text(size=12),
                                                                    plot.title=element_text(size=16)) 

ggsave("NextSeq_16S_rda_try2_same_but_wo_taxa_try2.png", width= 12, height = 6)











