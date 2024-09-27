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
      mainbar.y.label = "Fjord Intersections", sets.x.label = "Taxa per fjord", order.by = "freq",  text.scale = c(1.5, 1.5, 1.4, 1.3, 1.8, 1.5),
      sets.bar.color = c("#0edaff", "#5e91e4", "#ab55c7", "#ff00ad"), main.bar.color = "grey30", matrix.color = "grey30", queries = list(
        list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden"), active = T, color = "#ffaaf2"),
        list(query = intersects, params = list("Wijdefjorden", "Rijpfjorden"), active = T, color = "#8abee4"),
        list(query = intersects, params = list("van_Mijenfjorden", "Kongsfjorden", "Wijdefjorden"), active = T, color = "#c5b4eb"),
        list(query = intersects, params = "van_Mijenfjorden", active = T, color = "#ff00ad"),
        list(query = intersects, params = "Kongsfjorden", active = T, color = "#ab55c7"),
        list(query = intersects, params = "Wijdefjorden", active = T, color = "#5e91e4"),
        list(query = intersects, params = "Rijpfjorden", active = T, color = "#0edaff"))) #+
# ggplot2::ggtitle('MiSeq 18S')



# w colours 2
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
