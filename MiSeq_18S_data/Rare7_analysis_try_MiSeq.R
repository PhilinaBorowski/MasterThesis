## Rare7 - MiSeq values

library(Rare7)
library(dplyr)
library(tidyr)
library(microbiome)
library(ggplot2)

## set wd
setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/Rare7")

# ----------------------------------

#### PRE-work
## start filtering only ASVs with species level 

#Mi_data <- read.csv("MiSeq_18S_reduced_no_1rowsum.csv", row.names = 1)

#Mi_spec_no_NAs <- Mi_data %>% drop_na()
#str(Mi_spec_no_NAs)
## total of 1545 ASVs left

#write.csv(Mi_spec_no_NAs, "MiSeq_18S_reduced_no_NAs_in_spec.csv")


# ------------------------------------------------------------
# write a data frame w/ lat + long + species + numb + habitat

#'data.frame': 5 obs. of 5 variables:
#$ specie : Factor w/ 5 levels "sp.1","sp.2",..: 1 2 3 4 5
#$ lat : num -4.7 -14.7 -14.7 -24.7 -24.7
#$ long : num -42.9 -39.9 -45.9 -45.1 -45.9
#$ NumIndiv: int 2 1 3 1 32
#$ habitat : Factor w/ 2 levels "forest","savanna": 2 1 2 1 2


# create phyloseq object

asvsm <- read.csv("MiSeq_abundance_only_specs.csv")
taxam <- read.csv("MiSeq_taxa_only_specs.csv")                 
meta_samplesm <- read.csv("meta_phyloseq_18S_miseq.csv")  

asvsm <- asvsm %>% tibble::column_to_rownames("otu")
taxam <- taxam %>% tibble::column_to_rownames("otu")
meta_samplesm <- meta_samplesm %>% tibble::column_to_rownames("sample")

asvsm <- as.matrix(asvsm)
taxam <- as.matrix(taxam)

OTUm = otu_table(asvsm, taxa_are_rows = T)
TAXm = tax_table(taxam)
samplesm = sample_data(meta_samplesm)

# make sure sample names match
species_rarem <- phyloseq(OTUm, TAXm, samplesm)
species_rarem
# 1545 taxa
# 94 samples

## ---------------------------
## delete ? 

## 1st: merge same species 

#samem <- tax_glom(species_rarem, taxrank = "Species", NArm = F)
#samem # 398 taxa left

#samem@tax_table


#all_specsm <- samem@tax_table
# save all species in a list? 
# now tax tables
#all_speciesm <- tax_table(samem)[, "Species"]
#species_listm <- as(all_speciesm, "vector")


## 2nd: merge samples

same_mergedm <- merge_samples(species_rarem, "sample_merg", fun = sum)
# In asMethod(object) : NAs introduced by coercion
# NAs in Seuqncer, DNA; depth, sampling_method, sample_merg 
# but why? cause they're characters? 

same_mergedm # 33 samples left w/ 398 taxa

same_mergedm@otu_table <- t(same_mergedm@otu_table)
same_mergedm@otu_table            # st1 f02 for asv0001 : 4484
species_rarem@otu_table           # st1 f02 for asv0001 : 2080 + 36 + 2368 = 4484
### ---


#merged_specm <- otu_table(same_mergedm)
#merged_spec_listm <- as(merged_specm, "vector")
#dont think it makes sense 


## ----------------------------

# maybe merge first and then tax_glom? 
# makes more sense

same_merged_specsm <- tax_glom(same_mergedm, taxrank = "Species", NArm = F)
same_merged_specsm
# 33 samples w 398 taxa, same as above? 


same_merged_specsm@otu_table



# Extract abundance matrix from the phyloseq object
OTU1m = as(otu_table(same_merged_specsm), "matrix")
# transpose if necessary
if(taxa_are_rows(same_merged_specsm)){OTU1m <- t(OTU1m)}
# Coerce to data.frame
OTUdfm = as.data.frame(OTU1m)

write.csv((OTUdfm), "M18_abundances_merged_spec.csv")

# for taxa
TAX1m = as(tax_table(same_merged_specsm), "matrix")
# transpose if necessary
#if(taxa_are_rows(same_merged_specsm)){TAX1m <- t(TAX1m)}
# Coerce to data.frame
TAXdfm = as.data.frame(TAX1m)

write.csv((TAXdfm), "Mi_taxnomomy_merged_taxa_specs.csv")  


# delete everything but species
TAXdf_specm <- subset(TAXdfm, select = -c(1:8))
TAXdf_specm


# merge stations together ? 
# only look at filter/pump/net or all together? 

# try only looking at F02 + F3 cause only N for st1 and P for some stations
OTUdfm <- t(OTUdfm)
OTUdfm <- as.data.frame(OTUdfm)
str(OTUdfm)
head(OTUdfm)

OTUdfm$St1 <- OTUdfm$St1.F02 + OTUdfm$St1.F3
OTUdfm$St2 <- OTUdfm$St2.F02 + OTUdfm$St2.F3
OTUdfm$St3 <- OTUdfm$St3.F02 + OTUdfm$St3.F3
OTUdfm$St4 <- OTUdfm$St4.F02 + OTUdfm$St4.F3
OTUdfm$St5 <- OTUdfm$St5.F02 + OTUdfm$St5.F3
OTUdfm$St6 <- OTUdfm$St6.F02 + OTUdfm$St6.F3
OTUdfm$St7 <- OTUdfm$St7.F02 + OTUdfm$St7.F3
OTUdfm$St8 <- OTUdfm$St8.F02 + OTUdfm$St8.F3

OTUdfm$St10 <- OTUdfm$St10.F02 + OTUdfm$St10.F3
OTUdfm$St11 <- OTUdfm$St11.F02 + OTUdfm$St11.F3
OTUdfm$St12 <- OTUdfm$St12.F02 + OTUdfm$St12.F3
OTUdfm$St13 <- OTUdfm$St13.F02 + OTUdfm$St.13.F3

str(OTUdfm)

OTU113 <- subset(OTUdfm, select = -c(1:33)) # save the  merged stations
OTU113


# merge taxa + abundance datasets 

merge(x = OTU113, y = TAXdf_specm, by = "row.names") # works!!!!
merged_tax_abund <- merge(x = OTU113, y = TAXdf_specm, by = "row.names")
str(merged_tax_abund)


##check the cafeteria

#cafe <- merged_tax_abund%>%filter(Species == "Cafeteria_burkhardae")
#taxa_t <- as.data.frame(taxa)
#cafe_t <- taxa_t%>%filter(Species == "Cafeteria_burkhardae")
#asvs_t <- as.data.frame(asvs)
#cafe_OTU <- asvs_t[row.names(asvs_t) == "ASV1468",]
# but need stat + count in long format 

merged_tax_abund_long <- merged_tax_abund %>%
                                   pivot_longer(cols=c("St1", "St2", "St3", "St4", "St5", "St6", "St7", "St8", "St10", "St11", "St12", "St13"),
                                   names_to = "Station", values_to = "Abundance")
# 4776 rows
# 398 obs for 12 station
398*12
# = 4776

merged_tax_abund_long
merged_tax_abund_long$Station <- as.factor(merged_tax_abund_long$Station)
merged_tax_abund_long$Abundance <- as.numeric(merged_tax_abund_long$Abundance)
str(merged_tax_abund_long)



# add long-/latitude 
# somehow??

same_merged_specsm
same_merged_specsm@sam_data

#samem
#samem@sam_data
#meta_same <- samem@sam_data

merged_tax_abund_long$CodSite <- rep(c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13), 398) # added stat to maybe then add lat/long? 

#merged_tax_abund_long$lat <- if("stat" = 1, "lat" = meta_same$lat)
  
  
#merge(merged_tax_abund_long, meta_same, by = "stat")
#merged_tax_abund_long
#meta_same

#merge(merged_tax_abund_long, meta_same, by = NULL)
# doesnt match 

#merge(merged_tax_abund_long, meta_same, all = T)


# try to add manually? but make sure that lat/long is according to stat
# merged_tax_abund_long$stat <- rep(c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13), 398)

merged_tax_abund_long$lat <- rep(c(77.7116, #1
                                   77.6414, #2
                                   77.8282, #3
                                   77.7899, #4
                                   79.7502, #5
                                   79.1248, #6
                                   79.5006, #7
                                   80.2502, #8
                                   80.4164, #10
                                   78.952,  #11
                                   78.9677, #12
                                   79.0144), 398) #13
merged_tax_abund_long


merged_tax_abund_long$long <- rep(c(13.3588, #1
                                    14.4717, #2
                                    16.5766, #3
                                    15.4992, #4
                                    15.5336, #5
                                    16.0184, #6
                                    15.7149, #7
                                    22.0844, #8
                                    22.0838, #10
                                    12.0089,  #11
                                    9.4601, #12
                                    11.4401), 398) #13

merged_tax_abund_long


mergedRare7 <- subset(merged_tax_abund_long, select = c(2:7))
mergedRare7

colnames(mergedRare7)[1] <- 'species'
colnames(mergedRare7)[2] <- 'habitat'
colnames(mergedRare7)[3] <- 'NumIndiv'

mergedRare7 <- mergedRare7[, c("CodSite", "species", "lat", "long", "NumIndiv", "habitat")]
mergedRare7

mergedRare7$species <- as.factor(mergedRare7$species)
mergedRare7$lat <- as.numeric(mergedRare7$lat)
mergedRare7$long <- as.numeric(mergedRare7$long)
mergedRare7$NumIndiv <- as.integer(mergedRare7$NumIndiv)
mergedRare7$habitat <- as.factor(mergedRare7$habitat)



## try Rare 7

rareData(mergedRare7)
# Species - Sample_area - Detection_area - Abundance - Habitats

rare7_data <- rareData(mergedRare7)
# sample_area       =  3 latitudinal belts 
# detection area    =  number of latitudinal were species occur



rareForms(rare7_data)
# Species - Form

# woher die NAs? 
# warum einige 0? 

rarity_data <- rareForms(rare7_data)

rarity_data$Form
length(which(rarity_data$Form == "common"))                         # 377
length(which(rarity_data$Form == "form1"))                          # 15
length(which(rarity_data$Form == "No abundance information"))       # 6

merged_rare_rarity_m18 <- merge(rare7_data, rarity_data, by = "Species")
write.csv(merged_rare_rarity_m18, "merged_rare_rarity_m18.csv")

FORM1 <- rarity_data%>%filter(Form =="form1")
FORM1 <- left_join(FORM1, rare7_data, by = "Species")
## ? show me the rare form1 ? 


Ancyromonas_micra               2 
Aplanochytrium_blankum          2 - but not in metapr2
Blastodinium_mangini            4 
Chaetoceros_diadema_1
Chrysochromulina_rotalis
Cinetochilum_ovale
Hemiselmis_cryptochromatica
Labyrinthulomycetes_LAB7_X_sp
Labyrinthulomycetes_LAB9_X_sp
Lecudina_phyllochaetopteri
MAST-12D_sp.
Minorisa_sp
Picomonas_judraskeda
Suctoria_XX_sp
Urospora_travisiae


# drop "No abundance data"
# cause no abundance, cause they only had values in net/pump samples

rarePlotdf <- data.frame(Form = c("rare1", "common"), Abundance = c(15, 377))


plot_m18 <- ggplot(rarePlotdf, aes(x = "", y = Abundance, fill = Form)) +
  geom_col() + geom_col(color = "black") +
  geom_text(aes(label = Abundance),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#9ACD32", "#B0E0E6")) + theme_light() +
  theme(axis.title = element_blank(), plot.title=element_text(size=18, hjust = 0.5), #change font size of plot title
        #legend.text=element_text(size=12), #change font size of legend text
        #legend.title=element_text(size=14),
        legend.position = "none") +  #change font size of legend title
        labs(title = "MiSeq 18S", xlab = "")

ggsave("Pie_rare.png", plot)

## all 3 plots ###
mixed_normal <- plot_m18 + plot_n18 + plot_n16 + plot_annotation(
  title = "Occurence of the 'Rare's' on Species level", 
  theme = theme(plot.title = element_text(size = 22, hjust = 0.5)))
ggsave("try_three_rare.png", mixed_normal, height = 7, width = 15)




### --------------------------------------
### --------------------------------------
### --------------------------------------

#try rare 7 but w asv level 


# create phyloseq object

asvsm_asv <- read.csv("MiSeq_abundances_all_94.csv")
taxam_asv <- read.csv("MiSeq_taxa_all_asv_names_as_taxa.csv")                 
meta_samplesm_asv <- read.csv("meta_phyloseq_18S_miseq.csv")  

asvsm_asv <- asvsm_asv %>% tibble::column_to_rownames("otu")
taxam_asv <- taxam_asv %>% tibble::column_to_rownames("otu")
meta_samplesm_asv <- meta_samplesm_asv %>% tibble::column_to_rownames("sample")

asvsm_asv <- as.matrix(asvsm_asv)
taxam_asv <- as.matrix(taxam_asv)

OTUm_asv = otu_table(asvsm_asv, taxa_are_rows = T)
TAXm_asv = tax_table(taxam_asv)
samplesm_asv = sample_data(meta_samplesm_asv)

# make sure sample names match
sample_names(OTUm_asv)
sample_names(samplesm_asv)

species_rarem_asv <- phyloseq(OTUm_asv, TAXm_asv, samplesm_asv)
species_rarem_asv
# 2523 taxa
# 94 samples


# merge samples
species_rarem_asv_merged <- merge_samples(species_rarem_asv, "sample_merg", fun = sum)


# Extract abundance matrix from the phyloseq object
OTU1m_asv = as(otu_table(species_rarem_asv_merged), "matrix")
# transpose if necessary
if(taxa_are_rows(species_rarem_asv_merged)){OTU1m_asv <- t(OTU1m_asv)}
# Coerce to data.frame
OTUdfm_asv = as.data.frame(OTU1m_asv)


# for taxa
TAX1m_asv = as(tax_table(species_rarem_asv_merged), "matrix")
# transpose if necessary
#if(taxa_are_rows(same_merged_specsm)){TAX1m <- t(TAX1m)}
# Coerce to data.frame
TAXdfm_asv = as.data.frame(TAX1m_asv)




# try only looking at F02 + F3 cause only N for st1 and P for some stations
OTUdfm_asv <- t(OTUdfm_asv)
OTUdfm_asv <- as.data.frame(OTUdfm_asv)
str(OTUdfm_asv)
head(OTUdfm_asv)
is.data.frame(OTUdfm_asv)


OTUdfm_asv$St1 <- OTUdfm_asv$St1.F02 + OTUdfm_asv$St1.F3
OTUdfm_asv$St2 <- OTUdfm_asv$St2.F02 + OTUdfm_asv$St2.F3
OTUdfm_asv$St3 <- OTUdfm_asv$St3.F02 + OTUdfm_asv$St3.F3
OTUdfm_asv$St4 <- OTUdfm_asv$St4.F02 + OTUdfm_asv$St4.F3
OTUdfm_asv$St5 <- OTUdfm_asv$St5.F02 + OTUdfm_asv$St5.F3
OTUdfm_asv$St6 <- OTUdfm_asv$St6.F02 + OTUdfm_asv$St6.F3
OTUdfm_asv$St7 <- OTUdfm_asv$St7.F02 + OTUdfm_asv$St7.F3
OTUdfm_asv$St8 <- OTUdfm_asv$St8.F02 + OTUdfm_asv$St8.F3

OTUdfm_asv$St10 <- OTUdfm_asv$St10.F02 + OTUdfm_asv$St10.F3
OTUdfm_asv$St11 <- OTUdfm_asv$St11.F02 + OTUdfm_asv$St11.F3
OTUdfm_asv$St12 <- OTUdfm_asv$St12.F02 + OTUdfm_asv$St12.F3
OTUdfm_asv$St13 <- OTUdfm_asv$St13.F02 + OTUdfm_asv$St.13.F3

str(OTUdfm_asv)

OTUdf_reduced <- subset(OTUdfm_asv, select = -c(1:33)) # save the  merged stations
OTUdf_reduced


# merge taxa + abundance datasets 

merge(x = OTUdf_reduced, y = TAXdfm_asv, by = "row.names") # works!!!!
merged_tax_abund_masv <- merge(x = OTUdf_reduced, y = TAXdfm_asv, by = "row.names")
str(merged_tax_abund_masv)



merged_tax_abund_masv_long <- merged_tax_abund_masv %>%
  pivot_longer(cols=c("St1", "St2", "St3", "St4", "St5", "St6", "St7", "St8", "St10", "St11", "St12", "St13"),
               names_to = "Station", values_to = "Abundance")
# 30276 rows
# 2523 obs for 12 station
2523*12

merged_tax_abund_masv_long$Station <- as.factor(merged_tax_abund_masv_long$Station)
merged_tax_abund_masv_long$Abundance <- as.numeric(merged_tax_abund_masv_long$Abundance)
str(merged_tax_abund_masv_long)


## add lat + long

merged_tax_abund_masv_long$CodSite <- rep(c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13), 2523) # added stat to maybe then add lat/long? 




# try to add manually? but make sure that lat/long is according to stat
# merged_tax_abund_long$stat <- rep(c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13), 398)

merged_tax_abund_masv_long$lat <- rep(c(77.7116, #1
                                   77.6414, #2
                                   77.8282, #3
                                   77.7899, #4
                                   79.7502, #5
                                   79.1248, #6
                                   79.5006, #7
                                   80.2502, #8
                                   80.4164, #10
                                   78.952,  #11
                                   78.9677, #12
                                   79.0144), 2523) #13
merged_tax_abund_masv_long


merged_tax_abund_masv_long$long <- rep(c(13.3588, #1
                                    14.4717, #2
                                    16.5766, #3
                                    15.4992, #4
                                    15.5336, #5
                                    16.0184, #6
                                    15.7149, #7
                                    22.0844, #8
                                    22.0838, #10
                                    12.0089,  #11
                                    9.4601, #12
                                    11.4401), 2523) #13

merged_tax_abund_masv_long



mergedRare7_masv <- subset(merged_tax_abund_masv_long, select = c(2:7))
mergedRare7_masv


colnames(mergedRare7_masv)[1] <- 'species'
colnames(mergedRare7_masv)[2] <- 'habitat'
colnames(mergedRare7_masv)[3] <- 'NumIndiv'

mergedRare7_masv <- mergedRare7_masv[, c("CodSite", "species", "lat", "long", "NumIndiv", "habitat")]
mergedRare7_masv

mergedRare7_masv$species <- as.factor(mergedRare7_masv$species)
mergedRare7_masv$lat <- as.numeric(mergedRare7_masv$lat)
mergedRare7_masv$long <- as.numeric(mergedRare7_masv$long)
mergedRare7_masv$NumIndiv <- as.integer(mergedRare7_masv$NumIndiv)
mergedRare7_masv$habitat <- as.factor(mergedRare7_masv$habitat)


## try Rare 7

rareData(mergedRare7_masv)
# Species - Sample_area - Detection_area - Abundance - Habitats

rare7_data_masv <- rareData(mergedRare7_masv)
# sample_area       =  3 latitudinal belts 
# detection area    =  number of latitudinal were species occur

rare7_data_masv %>% add_row(Species='sp', Sample_area=0, Detection_area=0, Abundance=0, Habitats=0)


rareForms(rare7_data_masv)
# Species - Form

# woher die NAs? 
# warum einige 0? 

rarity_data_masv <- rareForms(rare7_data_masv)

rarity_data_masv$Form
length(which(rarity_data_masv$Form == "common"))                         # 2304
length(which(rarity_data_masv$Form == "form1"))                          # 185
length(which(rarity_data_masv$Form == "No abundance information"))       # 34

2304+185+34
colnames(rarity_data_masv)[1] <- 'species'
merged_rare_rarity_m18_asv <- merge(mergedRare7_masv, rarity_data_masv, by = "species")
write.csv(merged_rare_rarity_m18_asv, "merged_rare_rarity_m18_asvs.csv")

#FORM1_asv <- rarity_data_masv%>%filter(Form =="form1")
#FORM1_asv <- left_join(FORM1, rare7_data, by = "Species")
## ? show me the rare form1 ? 


rarePlotdf_asv <- data.frame(Form = c("rare1", "common"), Abundance = c(185, 2304))


plot_m18_asv <- ggplot(rarePlotdf_asv, aes(x = "", y = Abundance, fill = Form)) +
  geom_col() + geom_col(color = "black") +
  geom_text(aes(label = Abundance),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#9ACD32", "#B0E0E6")) + theme_light() +
  theme(axis.title = element_blank(), plot.title=element_text(size=18, hjust = 0.5), #change font size of plot title
        #legend.text=element_text(size=12), #change font size of legend text
        #legend.title=element_text(size=14),
        legend.position = "none") +  #change font size of legend title
  labs(title = "MiSeq 18S", xlab = "")

ggsave("Pie_rare_masvs.png", plot_m18)

(185/2304)*100
# 8.03 % 

## all 3 plots ###
mixed <- plot_m18_asv + plot_n18_asv + plot_n16_asv + plot_annotation(
  title = "Occurence of the 'Rare's' on ASV level", 
  theme = theme(plot.title = element_text(size = 22, hjust = 0.5)))
ggsave("try_three_rare_asv.png", mixed, height = 7, width = 15)

