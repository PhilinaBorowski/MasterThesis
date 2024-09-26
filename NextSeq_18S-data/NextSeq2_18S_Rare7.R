### Rare 7 
### NextSeq 18S

library(Rare7)
library(dplyr)
library(tidyr)
library(microbiome)
library(ggplot2)


## set wd
setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/Rare7")

#### PRE-work
## start filtering only ASVs with species level 

N18_data <- read.csv("NextSeq_18S_reduced_no_1rowsum_133.csv", row.names = 1)

N18_spec_no_NAs <- N18_data %>% drop_na()
str(N18_spec_no_NAs)
## total of 4804 ASVs left

write.csv(N18_spec_no_NAs, "NextSeq_18S_reduced_no_NAs_in_spec.csv")

# write a data frame w/ lat + long + species + numb + habitat

#'data.frame': 5 obs. of 5 variables:
#$ specie : Factor w/ 5 levels "sp.1","sp.2",..: 1 2 3 4 5
#$ lat : num -4.7 -14.7 -14.7 -24.7 -24.7
#$ long : num -42.9 -39.9 -45.9 -45.1 -45.9
#$ NumIndiv: int 2 1 3 1 32
#$ habitat : Factor w/ 2 levels "forest","savanna": 2 1 2 1 2


# create phyloseq object
# first manually create the asvs + taxa csvs


asvs18 <- read.csv("NextSeq18S_abundance.csv")
taxa18 <- read.csv("NextSeq18S_taxa.csv")                 
meta_samples18 <- read.csv("metadata_phyloseq_NextSeq_18S.csv")  

asvs18 <- asvs18 %>% tibble::column_to_rownames("otu")
taxa18 <- taxa18 %>% tibble::column_to_rownames("otu")
meta_samples18 <- meta_samples18 %>% tibble::column_to_rownames("sample")

asvs18 <- as.matrix(asvs18)
taxa18 <- as.matrix(taxa18)

OTU18 = otu_table(asvs18, taxa_are_rows = T)
TAX18 = tax_table(taxa18)
samples18 = sample_data(meta_samples18)

# make sure sample names match
species_rare18 <- phyloseq(OTU18, TAX18, samples18)
species_rare18
# 4804 taxa
# 133 samples
species_rare18@sam_data

##
## merge samples

same_merged_n18 <- merge_samples(species_rare18, "sample_merg", fun = sum)
# In asMethod(object) : NAs introduced by coercion
# NAs in Seuqncer, DNA; depth, sampling_method, sample_merg 
# but why? cause they're characters? 

same_merged_n18
# 61 samples left
# still 4804 taxa

same_merged_n18@otu_table <- t(same_merged_n18@otu_table)

same_merged_n18@otu_table    # st1 f02 for asv0001 : 7983
species_rare18@otu_table     # st1 f02 for asv0001 : 1623 + 833 + 5527 = 7983


## merge species

same_merged_specs_n18 <- tax_glom(same_merged_n18, taxrank = "Species", NArm = F)
same_merged_specs_n18
# 61 samples w 688 taxa

same_merged_specs_n18@otu_table
# more counts now for asv0001 + st1 f02
# prob due to merging


# ---
# Extract abundance matrix from the phyloseq object
OTU1n18 = as(otu_table(same_merged_specs_n18), "matrix")
# transpose if necessary
if(taxa_are_rows(same_merged_specs_n18)){OTU1n18 <- t(OTU1n18)}
# Coerce to data.frame
OTUdfn18 = as.data.frame(OTU1n18)

# for taxa
TAX1n18 = as(tax_table(same_merged_specs_n18), "matrix")
# transpose if necessary
if(taxa_are_rows(same_merged_specs_n18)){TAX1n18 <- t(TAX1n18)}
# Coerce to data.frame
TAXdfn18 = as.data.frame(TAX1n18)


# delete everything but species
TAXdfn18 <- t(TAXdfn18)
TAXdf_spec_n18 <- subset(TAXdfn18, select = -c(1:8))
TAXdf_spec_n18

# merge taxa + abundance datasets 

OTUdfn18 <- t(OTUdfn18)

merged_tax_abund_n18 <- merge(x = OTUdfn18, y = TAXdf_spec_n18, by = "row.names")
str(merged_tax_abund_n18)


# merge samples into stations

colnames(merged_tax_abund_n18) <-  gsub("N18-", "", colnames(merged_tax_abund_n18))

merged_tax_abund_n18$St1 <- merged_tax_abund_n18$St1.F02 + merged_tax_abund_n18$St1.F3# + merged_tax_abund_n18$St1.N + merged_tax_abund_n18$St1.P #merged_tax_abund_n18$St1.dDNA +
merged_tax_abund_n18$St2 <- merged_tax_abund_n18$St2.F02 + merged_tax_abund_n18$St2.F3# + merged_tax_abund_n18$St2.N + merged_tax_abund_n18$St2.P merged_tax_abund_n18$St2.dDNA + 
merged_tax_abund_n18$St3 <- merged_tax_abund_n18$St3.F02 + merged_tax_abund_n18$St3.F3# + merged_tax_abund_n18$St3.N + merged_tax_abund_n18$St3.P merged_tax_abund_n18$St3.dDNA + 
merged_tax_abund_n18$St4 <- merged_tax_abund_n18$St4.F02 + merged_tax_abund_n18$St4.F3# + merged_tax_abund_n18$St4.N + merged_tax_abund_n18$St4.P merged_tax_abund_n18$St4.dDNA +
merged_tax_abund_n18$St5 <- merged_tax_abund_n18$St5.F02 + merged_tax_abund_n18$St5.F3# + merged_tax_abund_n18$St5.N + merged_tax_abund_n18$St5.P merged_tax_abund_n18$St5.dDNA + 
merged_tax_abund_n18$St6 <- merged_tax_abund_n18$St6.F02 + merged_tax_abund_n18$St6.F3# + merged_tax_abund_n18$St6.N + merged_tax_abund_n18$St6.P merged_tax_abund_n18$St6.dDNA + 
merged_tax_abund_n18$St7 <- merged_tax_abund_n18$St7.F02 + merged_tax_abund_n18$St7.F3# + merged_tax_abund_n18$St7.N + merged_tax_abund_n18$St7.P merged_tax_abund_n18$St7.dDNA +
merged_tax_abund_n18$St8 <- merged_tax_abund_n18$St8.F02 + merged_tax_abund_n18$St8.F3# + merged_tax_abund_n18$St8.N + merged_tax_abund_n18$St8.P merged_tax_abund_n18$St8.dDNA + 
#merged_tax_abund_n18$St9 <- merged_tax_abund_n18$St9.dDNA
merged_tax_abund_n18$St10 <- merged_tax_abund_n18$St10.F02 + merged_tax_abund_n18$St10.F3# + merged_tax_abund_n18$St10.N + merged_tax_abund_n18$St10.P merged_tax_abund_n18$St10.dDNA + 
merged_tax_abund_n18$St11 <- merged_tax_abund_n18$St11.F02 + merged_tax_abund_n18$St11.F3# + merged_tax_abund_n18$St11.N + merged_tax_abund_n18$St11.P  merged_tax_abund_n18$St11.dDNA +
merged_tax_abund_n18$St12 <- merged_tax_abund_n18$St12.F02 + merged_tax_abund_n18$St12.F3# + merged_tax_abund_n18$St12.N + merged_tax_abund_n18$St12.P merged_tax_abund_n18$St12.dDNA + 
merged_tax_abund_n18$St13 <- merged_tax_abund_n18$St13.F02 + merged_tax_abund_n18$St13.F3# + merged_tax_abund_n18$St13.N + merged_tax_abund_n18$St13.P merged_tax_abund_n18$St13.dDNA + 

str(merged_tax_abund_n18)

merged_tax_abund_n18 <- subset(merged_tax_abund_n18, select = -c(1:62)) # remove the 61 sample names + row.names (ASVs)
merged_tax_abund_n18



## create columns for lat + long

merged_tax_abund_long_n18 <- merged_tax_abund_n18 %>%
  pivot_longer(cols=c("St1", "St2", "St3", "St4", "St5", "St6", "St7", "St8", "St10", "St11", "St12", "St13"),
               names_to = "Station", values_to = "Abundance")
# 8256 rows
# 13 station
8256/12
# = 688

merged_tax_abund_long_n18$Station <- as.factor(merged_tax_abund_long_n18$Station)
merged_tax_abund_long_n18$Abundance <- as.numeric(merged_tax_abund_long_n18$Abundance)
str(merged_tax_abund_long_n18)


# add long-/latitude 

merged_tax_abund_long_n18$CodSite <- rep(c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13), 688) 


merged_tax_abund_long_n18$lat <- rep(c(77.7116, #1
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
                                       79.0144), 688) #13
merged_tax_abund_long_n18


merged_tax_abund_long_n18$long <- rep(c(13.3588, #1
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
                                        11.4401), 688) #13

merged_tax_abund_long_n18



# reorder so it fits for Rare7

Rare7_18S <- merged_tax_abund_long_n18 %>% rename(habitat = Station, NumIndiv = Abundance)
Rare7_18S <- Rare7_18S[, c("CodSite", "Species", "lat", "long", "NumIndiv", "habitat")]
Rare7_18S

Rare7_18S$Species <- as.factor(Rare7_18S$Species)
Rare7_18S$lat <- as.numeric(Rare7_18S$lat)
Rare7_18S$long <- as.numeric(Rare7_18S$long)
Rare7_18S$NumIndiv <- as.integer(Rare7_18S$NumIndiv)
Rare7_18S$habitat <- as.factor(Rare7_18S$habitat)

names(Rare7_18S)[2] <- "species"

write.csv((Rare7_18S), "Rare7_data_NextSeq_18S")  

## try Rare 7

rareData(Rare7_18S)
# Species - Sample_area - Detection_area - Abundance - Habitats

rare7_data_n18 <- rareData(Rare7_18S)
# sample_area       =  3 latitudinal belts 
# detection area    =  number of latitudinal were species occur

rareForms(rare7_data_n18)
# Species - Form
# doesnt work

rare7_data_n18 %>% add_row(Species='sp', Sample_area=0, Detection_area=0, Abundance=0, Habitats=0)

#######################

rarity_data18 <- rareForms(rare7_data_n18)

merged_rare_rarity_18 <- merge(rare7_data_n18, rarity_data18, by = "Species")
write.csv(merged_rare_rarity_18, "merged_rare_rarity_18.csv")

rarity_data18$Form
length(which(rarity_data18$Form == "common"))                         # 576
length(which(rarity_data18$Form == "form1"))                          # 29
length(which(rarity_data18$Form == "No abundance information"))       # 83 - net, pump or deep dna samples

FORM1 <- rarity_data18%>%filter(Form =="form1")
FORM1 <- left_join(FORM1, rarity_data18, by = "Species")


                        Species  Form.x Form.y
1             Aegilops_tauschii  form1  form1
2            Apostomatia_XX_sp.  form1  form1
3                Archerella_sp.  form1  form1
4      Arcocellulus_cornucervis  form1  form1
5   Asterochloris_phycobiontica  form1  form1
6      Caecitellus_paraparvulus  form1  form1
7            Candida_catenulata  form1  form1
8    Chaetoceros_cf_tortissimus  form1  form1
9          Challengeron_tizardi  form1  form1
10            Cryptomonas_ovata  form1  form1
11          Diacronema_vlkianum  form1  form1
12         Diaphanoeca_undulata  form1  form1
13 Endomyxa_Novel-clade-9_X_sp.  form1  form1
14           Ephelota_gemmipara  form1  form1
15          Goussia3_balatonica  form1  form1
16          Gymnodinium_smaydae  form1  form1
17          Hypocoma_acinetarum  form1  form1
18            Ichthyodinium_sp.  form1  form1
19      Kathablepharis_japonica  form1  form1
20 Labyrinthulaceae_ANT10_3_sp.  form1  form1
21               Mesodinium_sp.  form1  form1
22           Naviculaceae_X_sp.  form1  form1
23        Orthodonellidae_X_sp.  form1  form1
24               Parvamoeba_sp.  form1  form1
25                Pirula_salina  form1  form1
26      Protoperidinium_elegans  form1  form1
27            Selenidium_fallax  form1  form1
28          Telonemia_XXXXX_sp.  form1  form1
29        Tetracystis_vinatzeri  form1  form1
30          Thalassomyces_fagei  form1  form1
31          Trochilioides_recta  form1  form1
32     Ulvales-relatives_XX_sp.



rarePlotdf18 <- data.frame(Form = c("rare1", "common"), Abundance = c(32, 576))


plot_n18 <- ggplot(rarePlotdf18, aes(x = "", y = Abundance, fill = Form)) +
  geom_col() + geom_col(color = "black") +
  geom_text(aes(label = Abundance),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#9ACD32", "#B0E0E6")) + theme_light() +
  theme(axis.title = element_blank(), plot.title=element_text(size=18, hjust = 0.5), #change font size of plot title
        #legend.text=element_text(size=12), #change font size of legend text
        #legend.title=element_text(size=14),
        legend.position = "none") +  #change font size of legend title
  labs(title = "NextSeq 18S", xlab = "")

ggsave("Pie_rare_Next18S.png", plot_n18)

(29/576)*100
# 5.04 %




## ------------------------------
## ------------------------------
## ------------------------------



#try rare 7 but w asv level 


# create phyloseq object

asvsn18_asv <- read.csv("NextSeq18_abundances_asvs.csv")
taxan18_asv <- read.csv("NextSeq18S_taxa_asvs.csv")                 
meta_samplesn18_asv <- read.csv("metadata_phyloseq_NextSeq_18S.csv")  

asvsn18_asv <- asvsn18_asv %>% tibble::column_to_rownames("otu")
taxan18_asv <- taxan18_asv %>% tibble::column_to_rownames("otu")
meta_samplesn18_asv <- meta_samplesn18_asv %>% tibble::column_to_rownames("sample")

asvsn18_asv <- as.matrix(asvsn18_asv)
taxan18_asv <- as.matrix(taxan18_asv)

OTUn18_asv = otu_table(asvsn18_asv, taxa_are_rows = T)
TAXn18_asv = tax_table(taxan18_asv)
samplesn18_asv = sample_data(meta_samplesn18_asv)


species_raren18_asv <- phyloseq(OTUn18_asv, TAXn18_asv, samplesn18_asv)
species_raren18_asv
# 10425 taxa
# 133 samples


# merge samples
species_raren18_asv_merged <- merge_samples(species_raren18_asv, "sample_merg", fun = sum)


# Extract abundance matrix from the phyloseq object
OTU1n18_asv = as(otu_table(species_raren18_asv_merged), "matrix")
# transpose if necessary
if(taxa_are_rows(species_raren18_asv_merged)){OTU1n18_asv <- t(OTU1n18_asv)}
# Coerce to data.frame
OTUdfn18_asv = as.data.frame(OTU1n18_asv)


# for taxa
TAXn18_asv = as(tax_table(species_raren18_asv_merged), "matrix")
# transpose if necessary
#if(taxa_are_rows(same_merged_specsm)){TAX1m <- t(TAX1m)}
# Coerce to data.frame
TAXn18_asv = as.data.frame(TAXn18_asv)




# try only looking at F02 + F3 cause only N for st1 and P for some stations
OTUdfn18_asv <- t(OTUdfn18_asv)
OTUdfn18_asv <- as.data.frame(OTUdfn18_asv)
str(OTUdfn18_asv)
head(OTUdfn18_asv)
is.data.frame(OTUdfn18_asv)



# merge taxa + abundance datasets 

merge(x = OTUdfn18_asv, y = TAXn18_asv, by = "row.names") 
merged_tax_abund_n18asv <- merge(x = OTUdfn18_asv, y = TAXn18_asv, by = "row.names")
str(merged_tax_abund_n18asv)

# merge samples into stations

colnames(merged_tax_abund_n18asv) <-  gsub("N18-", "", colnames(merged_tax_abund_n18asv))

is.data.frame(merged_tax_abund_n18asv)

merged_tax_abund_n18asv$St1 <- merged_tax_abund_n18asv$St1.F02 + merged_tax_abund_n18asv$St1.F3# + merged_tax_abund_n18$St1.N + merged_tax_abund_n18$St1.P #merged_tax_abund_n18$St1.dDNA +
merged_tax_abund_n18asv$St2 <- merged_tax_abund_n18asv$St2.F02 + merged_tax_abund_n18asv$St2.F3# + merged_tax_abund_n18$St2.N + merged_tax_abund_n18$St2.P merged_tax_abund_n18$St2.dDNA + 
merged_tax_abund_n18asv$St3 <- merged_tax_abund_n18asv$St3.F02 + merged_tax_abund_n18asv$St3.F3# + merged_tax_abund_n18$St3.N + merged_tax_abund_n18$St3.P merged_tax_abund_n18$St3.dDNA + 
merged_tax_abund_n18asv$St4 <- merged_tax_abund_n18asv$St4.F02 + merged_tax_abund_n18asv$St4.F3# + merged_tax_abund_n18$St4.N + merged_tax_abund_n18$St4.P merged_tax_abund_n18$St4.dDNA +
merged_tax_abund_n18asv$St5 <- merged_tax_abund_n18asv$St5.F02 + merged_tax_abund_n18asv$St5.F3# + merged_tax_abund_n18$St5.N + merged_tax_abund_n18$St5.P merged_tax_abund_n18$St5.dDNA + 
merged_tax_abund_n18asv$St6 <- merged_tax_abund_n18asv$St6.F02 + merged_tax_abund_n18asv$St6.F3# + merged_tax_abund_n18$St6.N + merged_tax_abund_n18$St6.P merged_tax_abund_n18$St6.dDNA + 
merged_tax_abund_n18asv$St7 <- merged_tax_abund_n18asv$St7.F02 + merged_tax_abund_n18asv$St7.F3# + merged_tax_abund_n18$St7.N + merged_tax_abund_n18$St7.P merged_tax_abund_n18$St7.dDNA +
merged_tax_abund_n18asv$St8 <- merged_tax_abund_n18asv$St8.F02 + merged_tax_abund_n18asv$St8.F3# + merged_tax_abund_n18$St8.N + merged_tax_abund_n18$St8.P merged_tax_abund_n18$St8.dDNA + 
#merged_tax_abund_n18$St9 <- merged_tax_abund_n18$St9.dDNA
merged_tax_abund_n18asv$St10 <- merged_tax_abund_n18asv$St10.F02 + merged_tax_abund_n18asv$St10.F3# + merged_tax_abund_n18$St10.N + merged_tax_abund_n18$St10.P merged_tax_abund_n18$St10.dDNA + 
merged_tax_abund_n18asv$St11 <- merged_tax_abund_n18asv$St11.F02 + merged_tax_abund_n18asv$St11.F3# + merged_tax_abund_n18$St11.N + merged_tax_abund_n18$St11.P  merged_tax_abund_n18$St11.dDNA +
merged_tax_abund_n18asv$St12 <- merged_tax_abund_n18asv$St12.F02 + merged_tax_abund_n18asv$St12.F3# + merged_tax_abund_n18$St12.N + merged_tax_abund_n18$St12.P merged_tax_abund_n18$St12.dDNA + 
merged_tax_abund_n18asv$St13 <- merged_tax_abund_n18asv$St13.F02 + merged_tax_abund_n18asv$St13.F3# + merged_tax_abund_n18$St13.N + merged_tax_abund_n18$St13.P merged_tax_abund_n18$St13.dDNA + 

str(merged_tax_abund_n18asv)

merged_tax_abund_n18asv <- subset(merged_tax_abund_n18asv, select = -c(1:62)) # remove the 61 sample names + row.names (ASVs)
merged_tax_abund_n18asv

merged_tax_abund_n18asv_long <- merged_tax_abund_n18asv %>%
  pivot_longer(cols=c("St1", "St2", "St3", "St4", "St5", "St6", "St7", "St8", "St10", "St11", "St12", "St13"),
               names_to = "Station", values_to = "Abundance")
# 125100 rows
# 10425 obs for 12 station
10425*12

merged_tax_abund_n18asv_long$Station <- as.factor(merged_tax_abund_n18asv_long$Station)
merged_tax_abund_n18asv_long$Abundance <- as.numeric(merged_tax_abund_n18asv_long$Abundance)
str(merged_tax_abund_n18asv_long)


## add lat + long

merged_tax_abund_n18asv_long$CodSite <- rep(c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13), 10425) # added stat to maybe then add lat/long? 




# try to add manually? but make sure that lat/long is according to stat
# merged_tax_abund_long$stat <- rep(c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13), 398)

merged_tax_abund_n18asv_long$lat <- rep(c(77.7116, #1
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
                                          79.0144), 10425) #13
merged_tax_abund_n18asv_long


merged_tax_abund_n18asv_long$long <- rep(c(13.3588, #1
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
                                           11.4401), 10425) #13

merged_tax_abund_n18asv_long



#mergedRare7_n18asv <- subset(merged_tax_abund_n18asv_long, select = c(2:7))
mergedRare7_n18asv <-merged_tax_abund_n18asv_long


mergedRare7_n18asv <- mergedRare7_n18asv %>% rename(habitat = Station, NumIndiv = Abundance)
mergedRare7_n18asv <- mergedRare7_n18asv[, c("CodSite", "species", "lat", "long", "NumIndiv", "habitat")]
mergedRare7_n18asv

mergedRare7_n18asv$species <- as.factor(mergedRare7_n18asv$species)
mergedRare7_n18asv$lat <- as.numeric(mergedRare7_n18asv$lat)
mergedRare7_n18asv$long <- as.numeric(mergedRare7_n18asv$long)
mergedRare7_n18asv$NumIndiv <- as.integer(mergedRare7_n18asv$NumIndiv)
mergedRare7_n18asv$habitat <- as.factor(mergedRare7_n18asv$habitat)


## try Rare 7

rareData(mergedRare7_n18asv)
# Species - Sample_area - Detection_area - Abundance - Habitats

rare7_data_n18asv <- rareData(mergedRare7_n18asv)
# sample_area       =  3 latitudinal belts 
# detection area    =  number of latitudinal were species occur

rare7_data_n18asv %>% add_row(Species='sp', Sample_area=0, Detection_area=0, Abundance=0, Habitats=0)


rareForms(rare7_data_n18asv)
# Species - Form

# woher die NAs? 
# warum einige 0? 

rarity_data_n18asv <- rareForms(rare7_data_n18asv)

rarity_data_n18asv$Form
length(which(rarity_data_n18asv$Form == "common"))                      # 7651
length(which(rarity_data_n18asv$Form == "form1"))                       # 1329
length(which(rarity_data_n18asv$Form == "form2"))                       # 0
length(which(rarity_data_n18asv$Form == "No abundance information"))    # 1445

7651+1329+1445
rarity_data_n18asv$species <- rarity_data_n18asv %>% rename(species = Species)
merged_rare_rarity_n18_asv <- merge(mergedRare7_n18asv, rarity_data_n18asv, by = "row.names")
write.csv(merged_rare_rarity_n18_asv, "merged_rare_rarity_n18_asvs.csv")

#FORM1_asv <- rarity_data_masv%>%filter(Form =="form1")
#FORM1_asv <- left_join(FORM1, rare7_data, by = "Species")
## ? show me the rare form1 ? 


rarePlotdf_asv_n18 <- data.frame(Form = c("rare1", "common"), Abundance = c(1329, 7651))


plot_n18_asv <- ggplot(rarePlotdf_asv_n18, aes(x = "", y = Abundance, fill = Form)) +
  geom_col() + geom_col(color = "black") +
  geom_text(aes(label = Abundance),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#9ACD32", "#B0E0E6")) + theme_light() +
  theme(axis.title = element_blank(), plot.title=element_text(size=18, hjust = 0.5), #change font size of plot title
        #legend.text=element_text(size=12), #change font size of legend text
        #legend.title=element_text(size=14),
        legend.position = "none") +  #change font size of legend title
  labs(title = "NextSeq 18S", xlab = "")

ggsave("Pie_rare_n18asvs.png", plot_n18)

(1329/7651)*100
# 17.37
