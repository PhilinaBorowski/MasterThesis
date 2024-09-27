### MiSeq - Amplicon sequences


### set wd ###
setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/MiSeq_results")

### packages

library(dplyr)

### read in files

#???
#nifH_ASVTable <- read.csv("/Users/corahoerstmann/Documents/MIO_eIMPACT/DNA/AllGenetics_Amplicons/nifH_amplicons/seqtab_all_eIMPACT_nifH.csv", row.names = 1)
#nifH_taxonomy <- read.csv2("/Users/corahoerstmann/Documents/MIO_eIMPACT/DNA/AllGenetics_Amplicons/nifH_amplicons/ASVs_Taxonomy_eIMPACT_v2_0_5.tsv", row.names = 1, sep = "\t")

#???
#X16S_ASVTable <- read.csv("/Users/corahoerstmann/Documents/MIO_eIMPACT/DNA/AllGenetics_Amplicons/16S_amplicons/seqtab_all_eIMPACT_16S_new.csv", row.names = 1)
#X16S_taxonomy <- read.csv2("/Users/corahoerstmann/Documents/MIO_eIMPACT/DNA/AllGenetics_Amplicons/16S_amplicons/taxonomy_eIMPACT_16S_new_SILVAv138.1.csv", row.names = 1)

X18S_ASVTable <- read.csv("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/MiSeq_results/MiSeq_18S_seqtab_all.csv", row.names = 1)
X18S_taxonomy <- read.csv("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/MiSeq_results/MiSeq_18S_assigned_all.csv", row.names = 1)

str(X18S_taxonomy) # 4263
str(X18S_ASVTable) # 4263

###

#Remove metazoans # only for 16S
#mitochondria, chloroplasts

#filter out mitochondria and chloroplasts
#taxonomy_18S_r <- X18S_taxonomy%>%dplyr::filter(Family == "Mitochondria") 
#X18S_taxonomy <- X18S_taxonomy%>%dplyr::filter(!rownames(X16S_taxonomy) %in% rownames(taxonomy_16S_r)) 

#taxonomy_16S_r <-X16S_taxonomy%>%dplyr::filter(Order == "Chloroplast") 
#X16S_taxonomy <- X16S_taxonomy%>%dplyr::filter(!rownames(X16S_taxonomy) %in% rownames(taxonomy_16S_r)) 

#rm(taxonomy_16S_r)

#X16S_ASVTable <- X16S_ASVTable%>%dplyr::filter(rownames(X16S_ASVTable) %in% rownames(X16S_taxonomy))



##filter out metazoans

#taxonomy_18S_o <- X18S_taxonomy%>%dplyr::filter(!Order == "Metazoa") 
#X18S_taxonomy <- X18S_taxonomy%>%dplyr::filter(rownames(X18S_taxonomy) %in% rownames(taxonomy_18S_o)) 

#taxonomy_18S_c <- X18S_taxonomy%>%dplyr::filter(!Class == "Metazoa") 
#X18S_taxonomy <- X18S_taxonomy%>%dplyr::filter(rownames(X18S_taxonomy) %in% rownames(taxonomy_18S_c)) 

#taxonomy_18S_d <- X18S_taxonomy %>% dplyr::filter(!Division == "Metazoa") 
#X18S_taxonomy <- X18S_taxonomy %>% dplyr::filter(rownames(X18S_taxonomy) %in% rownames(taxonomy_18S_d)) 

#str(taxonomy_18S_d) # 3802 



metazoa <- X18S_taxonomy %>% dplyr::filter(Division == "Metazoa") 
str(metazoa) # 152 obs. of 105

other18S <- X18S_taxonomy%>%dplyr::filter(!Division %in% metazoa$Division)

# 4263 ASVs before 
4263 - 152

write.csv(other18S, "/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/MiSeq_results/MiSeq_metabarcoding_wo_metazoa.csv", quote=FALSE )

