## upset NextSeq 18S

library(UpSetR)
library(vegan)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggpubr)



setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/upset/")
wd <- setwd("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/upset/")


## use table wo all samples
up_next18 <- read.csv("/isibhv/projects/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_52_UJOHN_006B_Eukaryotes_HE627_Run2/upset/NextSeq2_18S_merged_wo_metazoa_no_txt133.csv")
up_next18 <- as.data.frame(up_next18)

## create a binary table of presence + absence
up_18_pa <- decostand(up_next18, method = "pa")

# into dataframe
up_18_pa <- as.data.frame(up_18_pa)


# merge columns into fjords

up_18_pa$van_Mijenfjorden <- as.integer(up_18_pa$St1.dDNA|up_18_pa$St1.N|
                                          up_18_pa$St1.F02.1|up_18_pa$St1.F02.2|up_18_pa$St1.F02.3|
                                          up_18_pa$St1.F3.1|up_18_pa$St1.F3.2|up_18_pa$St1.F3.3|
                                          up_18_pa$St1.P1|up_18_pa$St1.P3|
                                          up_18_pa$St2.dDNA|up_18_pa$St2.N|
                                          up_18_pa$St2.F02.1|up_18_pa$St2.F02.2|up_18_pa$St2.F02.3|
                                          up_18_pa$St2.F3.1|up_18_pa$St2.F3.2|up_18_pa$St2.F3.3|
                                          up_18_pa$St2.P.1|up_18_pa$St2.P.2|up_18_pa$St2.P.3|
                                          up_18_pa$St3.dDNA|up_18_pa$St3.N|
                                          up_18_pa$St3.F02.1|up_18_pa$St3.F02.2|up_18_pa$St3.F02.3|
                                          up_18_pa$St3.F3.1|up_18_pa$St3.F3.2|up_18_pa$St3.F3.3|
                                          up_18_pa$St3.P.1|up_18_pa$St3.P.2|up_18_pa$St3.P.3|
                                          up_18_pa$St4.dDNA|up_18_pa$St4.N|
                                          up_18_pa$St4.F02.1|up_18_pa$St4.F02.2|up_18_pa$St4.F02.3|
                                          up_18_pa$St4.F3.1|up_18_pa$St4.F3.2|up_18_pa$St4.F3.3|
                                          up_18_pa$St4.P.1|up_18_pa$St4.P.2|up_18_pa$St4.P.3)
#43

up_18_pa$Wijdefjorden <- as.integer(up_18_pa$St5.dDNA|up_18_pa$St5.N|up_18_pa$St5.N.150|
                                      up_18_pa$St5.F02.1|up_18_pa$St5.F02.2|up_18_pa$St5.F02.3|
                                      up_18_pa$St1.F3.1|up_18_pa$St1.F3.2|up_18_pa$St5.F3.3|
                                      up_18_pa$St5.P.1|up_18_pa$St5.P.2|up_18_pa$St5.P.3|
                                      up_18_pa$St6.dDNA|up_18_pa$St6.N|
                                      up_18_pa$St6.F02.1|up_18_pa$St6.F02.2|up_18_pa$St6.F02.3|
                                      up_18_pa$St6.F3.1|up_18_pa$St6.F3.2|up_18_pa$St6.F3.3|
                                      up_18_pa$St6.P.1|up_18_pa$St6.P.2|up_18_pa$St6.P.3|
                                      up_18_pa$St7.dDNA|up_18_pa$St7.N|
                                      up_18_pa$St7.F02.1|up_18_pa$St7.F02.2|up_18_pa$St7.F02.3|
                                      up_18_pa$St7.F3.1|up_18_pa$St7.F3.2|up_18_pa$St7.F3.3|
                                      up_18_pa$St7.P.1|up_18_pa$St7.P.2|up_18_pa$St7.P.3)
# 34


up_18_pa$Rijpfjorden <- as.integer(up_18_pa$St8.dDNA|up_18_pa$St8.N|
                                     up_18_pa$St8.F02.1|up_18_pa$St8.F02.2|up_18_pa$St8.F02.3|
                                     up_18_pa$St8.F3.1|up_18_pa$St8.F3.2|up_18_pa$St8.F3.3|
                                     up_18_pa$St8.P.1|up_18_pa$St8.P.2|up_18_pa$St8.P.3|
                                     up_18_pa$St9.dDNA|
                                     up_18_pa$St10.dDNA|up_18_pa$St10.N|
                                     up_18_pa$St10.F02.1|up_18_pa$St10.F02.2|up_18_pa$St10.F02.3|
                                     up_18_pa$St10.F3.1|up_18_pa$St10.F3.2|up_18_pa$St10.F3.3|
                                     up_18_pa$St10.P.1|up_18_pa$St10.P.2|up_18_pa$St10.P.3)
# 23

up_18_pa$Kongsfjorden <- as.integer(up_18_pa$St11.dDNA|up_18_pa$St11.N|
                                      up_18_pa$St11.F02.1|up_18_pa$St11.F02.2|up_18_pa$St11.F02.3|
                                      up_18_pa$St11.F3.1|up_18_pa$St11.F3.2|up_18_pa$St11.F3.3|
                                      up_18_pa$St11.P.1|up_18_pa$St11.P.2|up_18_pa$St11.P.3|
                                      up_18_pa$St12.dDNA|up_18_pa$St12.N|
                                      up_18_pa$St12.F02.1|up_18_pa$St12.F02.2|up_18_pa$St12.F02.3|
                                      up_18_pa$St12.F3.1|up_18_pa$St12.F3.2|up_18_pa$St12.F3.3|
                                      up_18_pa$St12.P.1|up_18_pa$St12.P.2|up_18_pa$St12.P.3|
                                      up_18_pa$St13.dDNA|up_18_pa$St13.N|
                                      up_18_pa$St13.F02.1|up_18_pa$St13.F02.2|up_18_pa$St13.F02.3|
                                      up_18_pa$St13.F3.1|up_18_pa$St13.F3.2|up_18_pa$St13.F3.3|
                                      up_18_pa$St13.P.1|up_18_pa$St13.P.2|up_18_pa$St13.P.3)
# 33
33+23+34+43

head(up_18_pa)

# define the fjords and filter out the 0s
fjords_18 <- up_18_pa[,c(134:137)]
fjords_18 <- fjords_18 %>% filter_all(any_vars(. !=0))
filter_all(fjords_18, any_vars(. != 0))
is.na(fjords_18)


## upset plot
# w colours
upset(fjords_18, sets = c("Rijpfjorden", "Wijdefjorden", "Kongsfjorden", "van_Mijenfjorden"), keep.order = T, point.size = 3.8, line.size = 1.3,
      mainbar.y.label = "Fjord Intersections", sets.x.label = "Taxa per fjord", order.by = "freq",  text.scale = c(1.5, 1.5, 1.4, 1.3, 1.8, 1.5),
      sets.bar.color = c("#2fa1ff", "#8a8bff", "#bb6fff", "#e241ff"), main.bar.color = "grey30", matrix.color = "grey30", queries = list(
        list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden"), active = T, color = "#cf58ff"),
        list(query = intersects, params = list("Wijdefjorden", "Rijpfjorden"), active = T, color = "#5d96ff"))) #+
# ggplot2::ggtitle('MiSeq 18S')
ggsave("NextSeq_18S_UpSet.png", dpi = 600)



# -------------------------



# w colours 1
upset(fjords_18, sets = c("Rijpfjorden", "Wijdefjorden", "Kongsfjorden", "van_Mijenfjorden"), keep.order = T, point.size = 3.8, line.size = 1.3,
      mainbar.y.label = "Fjord Intersections", sets.x.label = "Taxa per fjord", order.by = "freq",  text.scale = c(2, 1.8, 1.4, 1.3, 1.8, 2),
      sets.bar.color = c("#0edaff", "#5e91e4", "#ab55c7", "#ff00ad"), main.bar.color = "grey30", matrix.color = "grey30", queries = list(
        list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden"), active = T, color = "#ffaaf2"),
        list(query = intersects, params = list("Wijdefjorden", "Rijpfjorden"), active = T, color = "#8abee4"),
        list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden", "Wijdefjorden"), active = T, color = "#c5b4eb"),
        list(query = intersects, params = "van_Mijenfjorden", active = T, color = "#ff00ad"),
        list(query = intersects, params = "Kongsfjorden", active = T, color = "#ab55c7"),
        list(query = intersects, params = "Wijdefjorden", active = T, color = "#5e91e4"),
        list(query = intersects, params = "Rijpfjorden", active = T, color = "#0edaff"))) #+
# ggplot2::ggtitle('MiSeq 18S')
ggsave("MiSeq_UpSet.png", dpi = 600)




# w colours 2
upset(fjords_18, sets = c("Rijpfjorden", "Wijdefjorden", "Kongsfjorden", "van_Mijenfjorden"), keep.order = T, point.size = 3.8, line.size = 1.3,
      mainbar.y.label = "Fjord Intersections", sets.x.label = "Taxa per fjord", order.by = "freq",  text.scale = c(1.5, 1.5, 1.4, 1.3, 1.8, 1.5),
      sets.bar.color = c("#2fa1ff", "#8a8bff", "#bb6fff", "#e241ff"), main.bar.color = "grey30", matrix.color = "grey30", queries = list(
        list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden"), active = T, color = "#cf58ff"),
        list(query = intersects, params = list("Wijdefjorden", "Rijpfjorden"), active = T, color = "#5d96ff"),
        list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden", "Wijdefjorden"), active = T, color = "#9677ff"),
        list(query = intersects, params = "van_Mijenfjorden", active = T, color = "#e241ff"),
        list(query = intersects, params = "Kongsfjorden", active = T, color = "#bb6fff"),
        list(query = intersects, params = "Wijdefjorden", active = T, color = "#8a8bff"),
        list(query = intersects, params = "Rijpfjorden", active = T, color = "#2fa1ff"))) 




## ----------------- 

MiSeq 


## upsetr 

#install.packages("UpSetR")
#install.packages("ComplexUpset")
library(UpSetR)
library(vegan)
library(ComplexHeatmap)
library(dplyr)
#library(ComplexUpset)
library(ggplot2)


setwd("Z:/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/UpSetR")
wd <- setwd("Z:/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/UpSetR")


## use table wo all samples
up_miseq <- read.csv("Z:/AG_John/Expeditions/HE627_Data/Sequences_MiSeq/UpSetR/MiSeq_metabarcoding_wo_metazoa_95.csv")
up_miseq <- as.data.frame(up_miseq)

#upset(fromExpression(up_miseq))

## create a binary table of presence + absence
up_pa <- decostand(up_miseq, method = "pa")

# into dataframe
up_pa <- as.data.frame(up_pa)

# merge columns into fjords

# Operators & and | perform element-wise operation producing result having length of the longer operand
# but different outocme if i use &
up_pa$van_Mijenfjorden <- as.integer(up_pa$St1.F02.1|up_pa$St1.F02.2|up_pa$St1.F02.3|
                                       up_pa$St1.F3.1|up_pa$St1.F3.2|up_pa$St1.F3.3|
                                       up_pa$St1.N|up_pa$St1.P1|                                 # deleted St1 P2
                                       up_pa$St1.P3|up_pa$St2.F02.1|up_pa$St2.F02.2|
                                       up_pa$St2.F02.3|up_pa$St2.F3.1|up_pa$St2.F3.2|
                                       up_pa$St2.F3.3|up_pa$St2.P1|up_pa$St2.P2|
                                       up_pa$St2.P3|up_pa$St3.F02.1|up_pa$St3.F02.2|
                                       up_pa$St3.F02.3|up_pa$St3.F3.1|up_pa$St3.F3.2|
                                       up_pa$St3.F3.3|up_pa$St3.P1|up_pa$St3.P2|
                                       up_pa$St3.P3|up_pa$St4.F02.1|up_pa$St4.F02.2|
                                       up_pa$St4.F02.3|up_pa$St4.F3.1|up_pa$St4.F3.2|
                                       up_pa$St4.F3.3|up_pa$St4.P1|up_pa$St4.P2|
                                       up_pa$St4.P3)

up_pa$Wijdefjorden <- as.integer(up_pa$St5.F02.1|up_pa$St5.F02.2|up_pa$St5.F02.3|
                                   up_pa$St5.F3.1|up_pa$St5.F3.2|up_pa$St5.F3.3|
                                   up_pa$St5.P1|up_pa$St5.P2|up_pa$St5.P3|
                                   up_pa$St6.F02.1|up_pa$St6.F02.2|up_pa$St6.F02.3|
                                   up_pa$St6.F3.1|up_pa$St6.F3.2|up_pa$St6.F3.3|
                                   up_pa$St6.P1|up_pa$St6.P2|up_pa$St6.P3|
                                   up_pa$St7.F02.1|up_pa$St7.F02.2|up_pa$St7.F02.3|
                                   up_pa$St7.F3.1|up_pa$St7.F3.2|up_pa$St7.F3.3|
                                   up_pa$St7.P1|up_pa$St7.P2|up_pa$St7.P3)

up_pa$Rijpfjorden <- as.integer(up_pa$St8.F02.1|up_pa$St8.F02.2|up_pa$St8.F02.3|
                                  up_pa$St8.F3.1|up_pa$St1.F3.2|up_pa$St8.F3.3|
                                  up_pa$St8.P1|up_pa$St8.P2|up_pa$St10.F02.1|
                                  up_pa$St10.F02.2|up_pa$St10.F02.3|up_pa$St10.F3.1|
                                  up_pa$St10.F3.2|up_pa$St10.F3.3)

up_pa$Kongsfjorden <- as.integer(up_pa$St11.F02.1|up_pa$St11.F02.2|up_pa$St11.F02.3|
                                   up_pa$St11.F3.1|up_pa$St11.F3.2|up_pa$St11.F3.3|
                                   up_pa$St12.F02.1|up_pa$St12.F02.2|up_pa$St12.F02.3|
                                   up_pa$St12.F3.1|up_pa$St12.F3.2|up_pa$St12.F3.3|
                                   up_pa$St13.F02.1|up_pa$St13.F02.2|up_pa$St13.F02.3|
                                   up_pa$St13.F3.1|up_pa$St13.F3.2|up_pa$St13.F3.3)

head(up_pa)




# define the fjords and filter out the 0s
fjords <- up_pa[,c(96:99)]
fjords <- fjords %>% filter_all(any_vars(. !=0))
filter_all(fjords, any_vars(. != 0))
is.na(fjords)


# w colours 1
upset(fjords, sets = c("Rijpfjorden", "Wijdefjorden", "Kongsfjorden", "van_Mijenfjorden"), keep.order = T, point.size = 3.8, line.size = 1.3,
      mainbar.y.label = "Fjord Intersections", sets.x.label = "Taxa per fjord", order.by = "freq",  text.scale = c(2, 1.8, 1.4, 1.3, 1.8, 2),
      sets.bar.color = c("#0edaff", "#5e91e4", "#ab55c7", "#ff00ad"), main.bar.color = "grey30", matrix.color = "grey30", queries = list(
        list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden"), active = T, color = "#ffaaf2"),
        list(query = intersects, params = list("Wijdefjorden", "Rijpfjorden"), active = T, color = "#8abee4"),
        list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden", "Wijdefjorden"), active = T, color = "#c5b4eb"),
        list(query = intersects, params = "van_Mijenfjorden", active = T, color = "#ff00ad"),
        list(query = intersects, params = "Kongsfjorden", active = T, color = "#ab55c7"),
        list(query = intersects, params = "Wijdefjorden", active = T, color = "#5e91e4"),
        list(query = intersects, params = "Rijpfjorden", active = T, color = "#0edaff"))) #+
# ggplot2::ggtitle('MiSeq 18S')
ggsave("MiSeq_UpSet.png", dpi = 600)



# w colours 2
upset(fjords, sets = c("Rijpfjorden", "Wijdefjorden", "Kongsfjorden", "van_Mijenfjorden"), keep.order = T, point.size = 3.8, line.size = 1.3,
      mainbar.y.label = "Fjord Intersections", sets.x.label = "Taxa per fjord", order.by = "freq",  text.scale = c(1.5, 1.5, 1.4, 1.3, 1.8, 1.5),
      sets.bar.color = c("#2fa1ff", "#8a8bff", "#bb6fff", "#e241ff"), main.bar.color = "grey30", matrix.color = "grey30", queries = list(
        list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden"), active = T, color = "#cf58ff"),
        list(query = intersects, params = list("Wijdefjorden", "Rijpfjorden"), active = T, color = "#5d96ff"),
        list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden", "Wijdefjorden"), active = T, color = "#9677ff"),
        list(query = intersects, params = "van_Mijenfjorden", active = T, color = "#e241ff"),
        list(query = intersects, params = "Kongsfjorden", active = T, color = "#bb6fff"),
        list(query = intersects, params = "Wijdefjorden", active = T, color = "#8a8bff"),
        list(query = intersects, params = "Rijpfjorden", active = T, color = "#2fa1ff"))) #+

 

## --

NextSeq 16S

## upset nextseq 16S

## 
library(UpSetR)
library(vegan)
library(ComplexHeatmap)
library(dplyr)
#library(ComplexUpset)
library(ggplot2)


setwd("Z:/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/upset/")
wd <- setwd("Z:/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/upset/")


## use table wo all samples

up_next16 <- read.csv("Z:/AG_John/Expeditions/HE627_Data/Sequences_NextSeq_Neu/Run_53_UJOHN_007_Prokaryotes_HE627_Dedmar_Jakob/HE627_16S/upset/NextSeq_16S_ampli_cleaned_just_numb.csv")
up_next16 <- as.data.frame(up_next16)


## create a binary table of presence + absence
up_nx_pa <- decostand(up_next16, method = "pa")

# into dataframe
up_nx_pa <- as.data.frame(up_nx_pa)


# merge columns into fjords

up_nx_pa$van_Mijenfjorden <- as.integer(up_nx_pa$Pro_St1.F02.1|up_nx_pa$Pro_St1.F02.2|up_nx_pa$Pro_St1.F02.3|
                                          up_nx_pa$Pro_St2.F02.1|up_nx_pa$Pro_St2.F02.2|up_nx_pa$Pro_St2.F02.3|
                                          up_nx_pa$Pro_St3.F02.1|up_nx_pa$Pro_St3.F02.2|up_nx_pa$Pro_St3.F02.3|
                                          up_nx_pa$Pro_St4.F02.1|up_nx_pa$Pro_St4.F02.2|up_nx_pa$Pro_St4.F02.3)

up_nx_pa$Wijdefjorden <- as.integer(up_nx_pa$Pro_St5.F02.1|up_nx_pa$Pro_St5.F02.2|up_nx_pa$Pro_St5.F02.3|
                                      up_nx_pa$Pro_St6.F02.1|up_nx_pa$Pro_St6.F02.2|up_nx_pa$Pro_St6.F02.3|
                                      up_nx_pa$Pro_St7.F02.1|up_nx_pa$Pro_St7.F02.2|up_nx_pa$Pro_St7.F02.3)

up_nx_pa$Rijpfjorden <- as.integer(up_nx_pa$Pro_St8.F02.1|up_nx_pa$Pro_St8.F02.2|up_nx_pa$Pro_St8.F02.3|
                                     up_nx_pa$Pro_St10.F02.1|up_nx_pa$Pro_St10.F02.2|up_nx_pa$Pro_St10.F02.3)

up_nx_pa$Kongsfjorden <- as.integer(up_nx_pa$Pro_St11.F02.1|up_nx_pa$Pro_St11.F02.2|up_nx_pa$Pro_St11.F02.3|
                                      up_nx_pa$Pro_St12.F02.1|up_nx_pa$Pro_St12.F02.2|up_nx_pa$Pro_St12.F02.3|
                                      up_nx_pa$Pro_St13.F02.1|up_nx_pa$Pro_St13.F02.2|up_nx_pa$Pro_St13.F02.3)


head(up_nx_pa)



# define the fjords and filter out the 0s
fjords_n <- up_nx_pa[,c(37:40)]
fjords_n <- fjords_n %>% filter_all(any_vars(. !=0))
filter_all(fjords_n, any_vars(. != 0))
is.na(fjords_n)


## upset plot
# w colours
upset(fjords_n, sets = c("Rijpfjorden", "Wijdefjorden", "Kongsfjorden", "van_Mijenfjorden"), keep.order = T, point.size = 3.8, line.size = 1.3,
      mainbar.y.label = "Fjord Intersections", sets.x.label = "Taxa per fjord", order.by = "freq",  text.scale = c(1.5, 1.5, 1.4, 1.3, 1.8, 1.5),
      sets.bar.color = c("#2fa1ff", "#8a8bff", "#bb6fff", "#e241ff"), main.bar.color = "grey30", matrix.color = "grey30", queries = list(
        list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden"), active = T, color = "#cf58ff"),
        list(query = intersects, params = list("Wijdefjorden", "Rijpfjorden"), active = T, color = "#5d96ff"))) #+
# ggplot2::ggtitle('MiSeq 18S')
ggsave("NextSeq_16S_UpSet.png", dpi = 600)


#------------------

# w colours 1
upset(fjords_n, sets = c("Rijpfjorden", "Wijdefjorden", "Kongsfjorden", "van_Mijenfjorden"), keep.order = T, point.size = 3.8, line.size = 1.3,
      mainbar.y.label = "Fjord Intersections", sets.x.label = "Taxa per fjord", order.by = "freq",  text.scale = c(2, 1.8, 1.4, 1.3, 1.8, 2),
      sets.bar.color = c("#0edaff", "#5e91e4", "#ab55c7", "#ff00ad"), main.bar.color = "grey30", matrix.color = "grey30", queries = list(
        list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden"), active = T, color = "#ffaaf2"),
        list(query = intersects, params = list("Wijdefjorden", "Rijpfjorden"), active = T, color = "#8abee4"),
        list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden", "Wijdefjorden"), active = T, color = "#c5b4eb"),
        list(query = intersects, params = "van_Mijenfjorden", active = T, color = "#ff00ad"),
        list(query = intersects, params = "Kongsfjorden", active = T, color = "#ab55c7"),
        list(query = intersects, params = "Wijdefjorden", active = T, color = "#5e91e4"),
        list(query = intersects, params = "Rijpfjorden", active = T, color = "#0edaff"))) 


# w colours 2
set.size(5,8)
upset(fjords_n, sets = c("Rijpfjorden", "Wijdefjorden", "Kongsfjorden", "van_Mijenfjorden"), keep.order = T, point.size = 3.8, line.size = 1.3,
      mainbar.y.label = "Fjord Intersections", sets.x.label = "Taxa per fjord", order.by = "freq",  text.scale = c(1.5, 1.5, 1.4, 1.3, 1.8, 1.5),
      sets.bar.color = c("#2fa1ff", "#8a8bff", "#bb6fff", "#e241ff"), main.bar.color = "grey30", matrix.color = "grey30", queries = list(
        list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden"), active = T, color = "#cf58ff"),
        list(query = intersects, params = list("Wijdefjorden", "Rijpfjorden"), active = T, color = "#5d96ff"),
        list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden", "Wijdefjorden"), active = T, color = "#9677ff"),
        list(query = intersects, params = "van_Mijenfjorden", active = T, color = "#e241ff"),
        list(query = intersects, params = "Kongsfjorden", active = T, color = "#bb6fff"),
        list(query = intersects, params = "Wijdefjorden", active = T, color = "#8a8bff"),
        list(query = intersects, params = "Rijpfjorden", active = T, color = "#2fa1ff")))





figure <- ggarrange(plot1, plot2, plot3,
                    labels = c("A", "B", "C"),
                    ncol = 2, nrow = 2)
figure

(plot1 | plot2) + plot_layout(guides='collect')



set.size(8, 3)
(
  upset(fjords_n, sets = c("Rijpfjorden", "Wijdefjorden", "Kongsfjorden", "van_Mijenfjorden"), keep.order = T, point.size = 3.8, line.size = 1.3,
      mainbar.y.label = "Fjord Intersections", sets.x.label = "Taxa per fjord", order.by = "freq",  text.scale = c(1.5, 1.5, 1.4, 1.3, 1.8, 1.5),
      sets.bar.color = c("#0edaff", "#5e91e4", "#ab55c7", "#ff00ad"), main.bar.color = "grey30", matrix.color = "grey30", queries = list(
        list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden"), active = T, color = "#ffaaf2"),
        list(query = intersects, params = list("Wijdefjorden", "Rijpfjorden"), active = T, color = "#8abee4"),
        list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden", "Wijdefjorden"), active = T, color = "#c5b4eb"),
        list(query = intersects, params = "van_Mijenfjorden", active = T, color = "#ff00ad"),
        list(query = intersects, params = "Kongsfjorden", active = T, color = "#ab55c7"),
        list(query = intersects, params = "Wijdefjorden", active = T, color = "#5e91e4"),
        list(query = intersects, params = "Rijpfjorden", active = T, color = "#0edaff"))) 
  + upset(fjords_n, sets = c("Rijpfjorden", "Wijdefjorden", "Kongsfjorden", "van_Mijenfjorden"), keep.order = T, point.size = 3.8, line.size = 1.3,
        mainbar.y.label = "Fjord Intersections", sets.x.label = "Taxa per fjord", order.by = "freq",  text.scale = c(1.5, 1.5, 1.4, 1.3, 1.8, 1.5),
        sets.bar.color = c("#2fa1ff", "#8a8bff", "#bb6fff", "#e241ff"), main.bar.color = "grey30", matrix.color = "grey30", queries = list(
          list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden"), active = T, color = "#cf58ff"),
          list(query = intersects, params = list("Wijdefjorden", "Rijpfjorden"), active = T, color = "#5d96ff"),
          list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden", "Wijdefjorden"), active = T, color = "#9677ff"),
          list(query = intersects, params = "van_Mijenfjorden", active = T, color = "#e241ff"),
          list(query = intersects, params = "Kongsfjorden", active = T, color = "#bb6fff"),
          list(query = intersects, params = "Wijdefjorden", active = T, color = "#8a8bff"),
          list(query = intersects, params = "Rijpfjorden", active = T, color = "#2fa1ff"))) 
)

