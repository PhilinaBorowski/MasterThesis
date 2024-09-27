## Rare7 - MiSeq values

library(Rare7)
library(dplyr)
library(tidyr)
library(microbiome)
library(ggplot2)
library(ggpubr)
#install.packages("ggpubr")

## set wd
setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/Rare7/")

# ----------------------------------

#### PRE-work
## start filtering only ASVs with species level 

N16_data <- read.csv("NextSeq_16S_cleaned_reduced_no_1rowsum.csv", row.names = 1)

N16_spec_no_NAs <- N16_data %>% drop_na()
str(N16_spec_no_NAs)
## total of 195 ASVs left

write.csv(N16_spec_no_NAs, "NextSeq_16S_reduced_no_NAs_in_spec.csv")

# write a data frame w/ lat + long + species + numb + habitat

#'data.frame': 5 obs. of 5 variables:
#$ specie : Factor w/ 5 levels "sp.1","sp.2",..: 1 2 3 4 5
#$ lat : num -4.7 -14.7 -14.7 -24.7 -24.7
#$ long : num -42.9 -39.9 -45.9 -45.1 -45.9
#$ NumIndiv: int 2 1 3 1 32
#$ habitat : Factor w/ 2 levels "forest","savanna": 2 1 2 1 2


# create phyloseq object
# first manually create the asvs + taxa csvs

asvs <- read.csv("N16_abundances_asvs_only_spec.csv")
taxa <- read.csv("N16_taxa_only_spec.csv")                 
meta_samples <- read.csv("NextSeq_16S_metadata.csv")  

asvs <- asvs %>% tibble::column_to_rownames("otu")
taxa <- taxa %>% tibble::column_to_rownames("otu")
meta_samples <- meta_samples %>% tibble::column_to_rownames("sample")

asvs <- as.matrix(asvs)
taxa <- as.matrix(taxa)

OTU = otu_table(asvs, taxa_are_rows = T)
TAX = tax_table(taxa)
samples = sample_data(meta_samples)

# make sure sample names match
species_rare <- phyloseq(OTU, TAX, samples)
species_rare
# 195 taxa
# 36 samples
species_rare@sam_data


## ---------------------------

###################### 
COPied MISEQ SCRIPT
######################


## merge samples

same_merged_n16 <- merge_samples(species_rare, "sample_merg", fun = sum)
# In asMethod(object) : NAs introduced by coercion
# NAs in Seuqncer, DNA; depth, sampling_method, sample_merg 
# but why? cause they're characters? 

same_merged_n16 # 12 samples left w/ 195 taxa

same_merged_n16@otu_table    # st1 f02 for asv0001 : 3582
species_rare@otu_table       # st1 f02 for asv0001 : 1064 + 1526 + 992 = 3582



## merge species

same_merged_specs_n16 <- tax_glom(same_merged_n16, taxrank = "Species", NArm = F)
same_merged_specs_n16
# 12 samples w 161 taxa, same as above? 


same_merged_specs_n16@otu_table



# Extract abundance matrix from the phyloseq object
OTU1n16 = as(otu_table(same_merged_specs_n16), "matrix")
# transpose if necessary
if(taxa_are_rows(same_merged_specs_n16)){OTU1n16 <- t(OTU1n16)}
# Coerce to data.frame
OTUdfn16 = as.data.frame(OTU1n16)

write.csv((OTUdfn16), "N16_abundances_merged_spec.csv")

# for taxa
TAX1n16 = as(tax_table(same_merged_specs_n16), "matrix")
# transpose if necessary
if(taxa_are_rows(same_merged_specs_n16)){TAX1n16 <- t(TAX1n16)}
# Coerce to data.frame
TAXdfn16 = as.data.frame(TAX1n16)

write.csv((TAXdfn16), "N16_taxnomomy_merged_taxa_specs.csv")  



# delete everything but species
TAXdf_spec_n16 <- subset(TAXdfn16, select = -c(1:5))
TAXdf_spec_n16

TAXdf_spec_n16$species <- paste(TAXdf_spec_n16$Genus, TAXdf_spec_n16$Species, sep="_")
TAXdf_spec_n16 
TAXdf_spec_only_n16 <- subset(TAXdf_spec_n16, select = -c(1:2))
TAXdf_spec_only_n16


# merge taxa + abundance datasets 

OTUdfn16 <- t(OTUdfn16)

merge(x = OTUdfn16, y = TAXdf_spec_only_n16, by = "row.names") 
merged_tax_abund_n16 <- merge(x = OTUdfn16, y = TAXdf_spec_only_n16, by = "row.names")
str(merged_tax_abund_n16)


# change station names to easier ones lol
merged_tax_abund_n16 <- merged_tax_abund_n16 %>% rename("St1" = "N16-St1.F02", "St2" = "N16-St2.F02", "St3" = "N16-St3.F02", "St4" = "N16-St4.F02", "St5" = "N16-St5.F02",
                                                        "St6" = "N16-St6.F02", "St7" = "N16-St7.F02", "St8" = "N16-St8.F02", "St10" = "N16-St10.F02", 
                                                        "St11" = "N16-St11.F02", "St12" = "N16-St12.F02", "St13" = "N16-St13.F02")



## create columns for lat + long

merged_tax_abund_long_n16 <- merged_tax_abund_n16 %>%
  pivot_longer(cols=c("St1", "St2", "St3", "St4", "St5", "St6", "St7", "St8", "St10", "St11", "St12", "St13"),
               names_to = "Station", values_to = "Abundance")
# 1932 rows
# 161 obs for 12 station
161*12
# = 1932

merged_tax_abund_long_n16
merged_tax_abund_long_n16$Station <- as.factor(merged_tax_abund_long_n16$Station)
merged_tax_abund_long_n16$Abundance <- as.numeric(merged_tax_abund_long_n16$Abundance)
str(merged_tax_abund_long_n16)



# add long-/latitude 

merged_tax_abund_long_n16$CodSite <- rep(c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13), 161) 


merged_tax_abund_long_n16$lat <- rep(c(77.7116, #1
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
                                   79.0144), 161) #13
merged_tax_abund_long_n16


merged_tax_abund_long_n16$long <- rep(c(13.3588, #1
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
                                    11.4401), 161) #13

merged_tax_abund_long_n16


mergedRare7_n16 <- subset(merged_tax_abund_long_n16, select = c(2:7))
mergedRare7_n16


mergedRare7_n16 <- mergedRare7_n16 %>% rename(habitat = Station, NumIndiv = Abundance)
mergedRare7_n16 <- mergedRare7_n16[, c("CodSite", "species", "lat", "long", "NumIndiv", "habitat")]
mergedRare7_n16

mergedRare7_n16$species <- as.factor(mergedRare7_n16$species)
mergedRare7_n16$lat <- as.numeric(mergedRare7_n16$lat)
mergedRare7_n16$long <- as.numeric(mergedRare7_n16$long)
mergedRare7_n16$NumIndiv <- as.integer(mergedRare7_n16$NumIndiv)
mergedRare7_n16$habitat <- as.factor(mergedRare7_n16$habitat)



## try Rare 7

rareData(mergedRare7_n16)
# Species - Sample_area - Detection_area - Abundance - Habitats

rare7_data_n16 <- rareData(mergedRare7_n16)
# sample_area       =  3 latitudinal belts 
# detection area    =  number of latitudinal were species occur

rare7_data_n16 %>% add_row(Species='sp', Sample_area=0, Detection_area=0, Abundance=0, Habitats=0)



rareForms(rare7_data_n16)
# Species - Form
#######
# change function cause AI told me lol

##########
rarity_data_n16 <- rareForms(rare7_data_n16)

merged_rare_rarity_16 <- merge(rare7_data_n16, rarity_data_n16, by = "Species")
write.csv(merged_rare_rarity_16, "merged_rare_rarity_16.csv")

rarity_data_n16$Form
length(which(rarity_data_n16$Form == "common"))                         # 150
length(which(rarity_data_n16$Form == "form1"))                          # 11
length(which(rarity_data_n16$Form == "No abundance information"))       # 0

rarity_data_n16$Form
length(which(merged_rare_rarity_16$Form == "common"))                         # 150
length(which(merged_rare_rarity_16$Form == "form1"))                          # 11
length(which(merged_rare_rarity_16$Form == "No abundance information"))       # 0

FORM1_16 <- rarity_data_n16%>%filter(Form =="form1")
FORM1_16 <- left_join(FORM1_16, rare7_data_n16, by = "Species")
## ? show me the rare form1 ? 
Species  Form Sample_area Detection_area Abundance Habitats
1           Blastococcus_saxobsidens form1           3              3         2       12
2              Brevundimonas_bullata form1           3              3         2       12
3  Candidatus Omnitrophus_magneticus form1           3              3         2       12
4               Cellulophaga_baltica form1           3              3         2       12
5             Chryseobacterium_aahli form1           3              3         2       12
6            Gemmatimonas_aurantiaca form1           3              3         2       12
7             Mucilaginibacter_rigui form1           3              3         2       12
8                 Mycoplasma_moatsii form1           3              3         2       12
9               Rhodococcus_fascians form1           3              3         2       12
10     Sulfurospirillum_arcachonense form1           3              3         2       12
11           Temperatibacter_marinus form1           3              3         2       12




#rarePlot <- rarity_data %>% filter(!Form == "No abundance information")
# 392 cause -6 no abund

rarePlotdf_n16 <- data.frame(Form = c("rare1", "common"), Abundance = c(11, 150))
## change abundance according to rarity data


plot_n16 <- ggplot(rarePlotdf_n16, aes(x = "", y = Abundance, fill = Form)) +
  geom_col() + geom_col(color = "black") +
  geom_text(aes(label = Abundance),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#9ACD32", "#B0E0E6")) + theme_light() +
  theme(axis.title = element_blank(), plot.title=element_text(size=18, hjust = 0.5), #change font size of plot title
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=14)) +  #change font size of legend title
  labs(title = "NextSeq 16S", xlab = "")

ggsave("Pie_rare_Next16S.png", plot)



## ------------------------------
## ------------------------------
## ------------------------------



#try rare 7 but w asv level 


# create phyloseq object

asvsn16_asv <- read.csv("NextSeq16_abundances_asvs.csv")
taxan16_asv <- read.csv("NextSeq16_taxa_asvs.csv")                 
meta_samplesn16_asv <- read.csv("NextSeq_16S_metadata.csv")  

asvsn16_asv <- asvsn16_asv %>% tibble::column_to_rownames("otu")
taxan16_asv <- taxan16_asv %>% tibble::column_to_rownames("otu")
meta_samplesn16_asv <- meta_samplesn16_asv %>% tibble::column_to_rownames("sample")

asvsn16_asv <- as.matrix(asvsn16_asv)
taxan16_asv <- as.matrix(taxan16_asv)

OTUn16_asv = otu_table(asvsn16_asv, taxa_are_rows = T)
TAXn16_asv = tax_table(taxan16_asv)
samplesn16_asv = sample_data(meta_samplesn16_asv)


species_raren16_asv <- phyloseq(OTUn16_asv, TAXn16_asv, samplesn16_asv)
species_raren16_asv
# 2044 taxa
# 36 samples


# merge samples
species_raren16_asv_merged <- merge_samples(species_raren16_asv, "sample_merg", fun = sum)


# Extract abundance matrix from the phyloseq object
OTU1n16_asv = as(otu_table(species_raren16_asv_merged), "matrix")
# transpose if necessary
if(taxa_are_rows(species_raren16_asv_merged)){OTU1n16_asv <- t(OTU1n16_asv)}
# Coerce to data.frame
OTUdfn16_asv = as.data.frame(OTU1n16_asv)


# for taxa
TAXn16_asv = as(tax_table(species_raren16_asv_merged), "matrix")
# transpose if necessary
#if(taxa_are_rows(same_merged_specsm)){TAX1m <- t(TAX1m)}
# Coerce to data.frame
TAXn16_asv = as.data.frame(TAXn16_asv)




# try only looking at F02 + F3 cause only N for st1 and P for some stations
OTUdfn16_asv <- t(OTUdfn16_asv)
OTUdfn16_asv <- as.data.frame(OTUdfn16_asv)
str(OTUdfn16_asv)
head(OTUdfn16_asv)
is.data.frame(OTUdfn16_asv)



# merge taxa + abundance datasets 

merge(x = OTUdfn16_asv, y = TAXn16_asv, by = "row.names") 
merged_tax_abund_n16asv <- merge(x = OTUdfn16_asv, y = TAXn16_asv, by = "row.names")
str(merged_tax_abund_n16asv)


# change station names to easier ones lol
merged_tax_abund_n16asv <- merged_tax_abund_n16asv %>% rename("St1" = "N16-St1.F02", "St2" = "N16-St2.F02", "St3" = "N16-St3.F02", "St4" = "N16-St4.F02", "St5" = "N16-St5.F02",
                                                        "St6" = "N16-St6.F02", "St7" = "N16-St7.F02", "St8" = "N16-St8.F02", "St10" = "N16-St10.F02", 
                                                        "St11" = "N16-St11.F02", "St12" = "N16-St12.F02", "St13" = "N16-St13.F02")

merged_tax_abund_n16asv_long <- merged_tax_abund_n16asv %>%
  pivot_longer(cols=c("St1", "St2", "St3", "St4", "St5", "St6", "St7", "St8", "St10", "St11", "St12", "St13"),
               names_to = "Station", values_to = "Abundance")
# 24528 rows
# 2044 obs for 12 station
2044*12

merged_tax_abund_n16asv_long$Station <- as.factor(merged_tax_abund_n16asv_long$Station)
merged_tax_abund_n16asv_long$Abundance <- as.numeric(merged_tax_abund_n16asv_long$Abundance)
str(merged_tax_abund_n16asv_long)


## add lat + long

merged_tax_abund_n16asv_long$CodSite <- rep(c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13), 2044) # added stat to maybe then add lat/long? 




# try to add manually? but make sure that lat/long is according to stat
# merged_tax_abund_long$stat <- rep(c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13), 398)

merged_tax_abund_n16asv_long$lat <- rep(c(77.7116, #1
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
                                        79.0144), 2044) #13
merged_tax_abund_n16asv_long


merged_tax_abund_n16asv_long$long <- rep(c(13.3588, #1
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
                                         11.4401), 2044) #13

merged_tax_abund_n16asv_long



mergedRare7_n16asv <- subset(merged_tax_abund_n16asv_long, select = c(2:7))
mergedRare7_n16asv


mergedRare7_n16asv <- mergedRare7_n16asv %>% rename(habitat = Station, NumIndiv = Abundance)
mergedRare7_n16asv <- mergedRare7_n16asv[, c("CodSite", "species", "lat", "long", "NumIndiv", "habitat")]
mergedRare7_n16asv

mergedRare7_n16asv$species <- as.factor(mergedRare7_n16asv$species)
mergedRare7_n16asv$lat <- as.numeric(mergedRare7_n16asv$lat)
mergedRare7_n16asv$long <- as.numeric(mergedRare7_n16asv$long)
mergedRare7_n16asv$NumIndiv <- as.integer(mergedRare7_n16asv$NumIndiv)
mergedRare7_n16asv$habitat <- as.factor(mergedRare7_n16asv$habitat)


## try Rare 7

rareData(mergedRare7_n16asv)
# Species - Sample_area - Detection_area - Abundance - Habitats

rare7_data_n16asv <- rareData(mergedRare7_n16asv)
# sample_area       =  3 latitudinal belts 
# detection area    =  number of latitudinal were species occur

rare7_data_n16asv %>% add_row(Species='sp', Sample_area=0, Detection_area=0, Abundance=0, Habitats=0)


rareForms(rare7_data_n16asv)
# Species - Form

# woher die NAs? 
# warum einige 0? 

rarity_data_n16asv <- rareForms(rare7_data_n16asv)

rarity_data_n16asv$Form
length(which(rarity_data_n16asv$Form == "common"))                         # 1716
length(which(rarity_data_n16asv$Form == "form1"))                          # 328
length(which(rarity_data_n16asv$Form == "form7"))                          # 0
length(which(rarity_data_n16asv$Form == "No abundance information"))       # 0

1716+328
rarity_data_n16asv$species <- rarity_data_n16asv %>% rename(species = Species)
merged_rare_rarity_n16_asv <- merge(mergedRare7_n16asv, rarity_data_n16asv, by = "row.names")
write.csv(merged_rare_rarity_n16_asv, "merged_rare_rarity_n16_asvs.csv")

#FORM1_asv <- rarity_data_masv%>%filter(Form =="form1")
#FORM1_asv <- left_join(FORM1, rare7_data, by = "Species")
## ? show me the rare form1 ? 


rarePlotdf_asv_n16 <- data.frame(Form = c("rare1", "common"), Abundance = c(328, 1716))


plot_n16_asv <- ggplot(rarePlotdf_asv_n16, aes(x = "", y = Abundance, fill = Form)) +
  geom_col() + geom_col(color = "black") +
  geom_text(aes(label = Abundance),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#9ACD32", "#B0E0E6")) + theme_light() +
  theme(axis.title = element_blank(), plot.title=element_text(size=18, hjust = 0.5), #change font size of plot title
        legend.text=element_text(size=12), #change font size of legend text
        legend.title=element_text(size=14)) +  #change font size of legend title
  labs(title = "NextSeq 16S", xlab = "")

ggsave("Pie_rare_n16asvs.png", plot_n16)

(328/1716)*100
# 19.11
