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


upset(fjords, sets = c("van_Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden"), mb.ratio = c(45, 55), order.by = "freq")
# works


upset(fjords, sets = c("van_Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden"), mb.ratio = c(45, 55), order.by = "degree")


upset(fjords)

upset(fjords, sets = c("van_Mijenfjorden", "Kongsfjorden", "Wijdefjorden", "Rijpfjorden"), keep.order = T)
upset(fjords, sets = c("Rijpfjorden", "Wijdefjorden", "Kongsfjorden", "van_Mijenfjorden"), keep.order = T)


upset(fjords, sets = c("Rijpfjorden", "Wijdefjorden", "Kongsfjorden", "van_Mijenfjorden"), keep.order = T, point.size = 4, line.size = 1.5,
      mainbar.y.label = "Fjord Intersections", sets.x.label = "Taxa per fjord")


# group by sets
upset(fjords, sets = c("Rijpfjorden", "Wijdefjorden", "Kongsfjorden", "van_Mijenfjorden"), keep.order = T, point.size = 4, line.size = 1.5,
      mainbar.y.label = "Fjord Intersections", sets.x.label = "Taxa per fjord", group.by = "sets") #title = "MiSeq 18S")




# w colours 1
upset(fjords, sets = c("Rijpfjorden", "Wijdefjorden", "Kongsfjorden", "van_Mijenfjorden"), keep.order = T, point.size = 3.8, line.size = 1.3,
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

